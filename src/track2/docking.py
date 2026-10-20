"""AutoDock Vina docking with a mandatory validation per target (rules: `protocol/addendum_docking.md`).

    python -m track2.docking --stage validate     # redocking and enrichment for every target (decides which targets may be used)
    python -m track2.docking --stage screen       # dock every approved small molecule on the validated targets (hit threshold from the validation)
    python -m track2.docking --stage top          # dock the top candidates on every validated target (context)

Needs `tools/vina_1.2.7_win.exe` (or `vina` on the PATH), Open Babel (`pip install openbabel-wheel`) and RDKit. Docking scores are hypotheses about binding, not
measurements; a target whose validation fails is not used.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
DOCK = Path("data/docking")
VINA = Path("tools/vina_1.2.7_win.exe")
EXHAUSTIVENESS, NUM_MODES, N_DECOYS, MAX_HEAVY, MAX_ROT = 8, 9, 200, 50, 15
REDOCK_RMSD_MAX, AUC_MIN, MIN_ACTIVES, PRECISION_MIN = 2.0, 0.70, 4, 0.5
SEED = 20261007


def vina_exe() -> str:
    if VINA.exists():
        return str(VINA)
    found = shutil.which("vina")
    if not found:
        sys.exit("AutoDock Vina not found: put vina_1.2.7_win.exe in tools/ (see tools/vina_download_record.json) or install it on the PATH.")
    return found


def _ob():
    from openbabel import openbabel as ob  # imported lazily: optional dependency
    return ob


def pdb_to_pdbqt(src: Path, dst: Path, receptor: bool) -> None:
    ob = _ob()
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdbqt")
    if receptor:
        conv.AddOption("r", ob.OBConversion.OUTOPTIONS)  # rigid: no torsion tree
    mol = ob.OBMol()
    conv.ReadFile(mol, str(src))
    mol.DeleteHydrogens()
    if receptor:
        mol.CorrectForPH(7.4)
    mol.AddHydrogens()
    conv.WriteFile(mol, str(dst))


def smiles_to_pdbqt(smiles: str, dst: Path, seed: int = SEED) -> bool:
    """RDKit 3D embedding (fixed seed, MMFF) then Open Babel PDBQT. False if the molecule cannot be embedded."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:  # noqa: BLE001 - optimisation is optional
        pass
    ob = _ob()
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("mol", "pdbqt")
    obmol = ob.OBMol()
    conv.ReadString(obmol, Chem.MolToMolBlock(mol))
    conv.WriteFile(obmol, str(dst))
    return True


def heavy_ok(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return bool(mol) and mol.GetNumHeavyAtoms() <= MAX_HEAVY and rdMolDescriptors.CalcNumRotatableBonds(mol) <= MAX_ROT


def pdbqt_atoms(path: Path, model: int = 1) -> np.ndarray:
    coords, current = [], 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL"):
            current = int(line.split()[1])
        if (line.startswith("ATOM") or line.startswith("HETATM")) and (current in (0, model)) and line[77:79].strip() not in ("H", "HD"):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)


def run_vina(receptor: Path, ligand: Path, centre: list[float], box: float, out: Path, seed: int = 1) -> float | None:
    cmd = [vina_exe(), "--receptor", str(receptor), "--ligand", str(ligand), "--center_x", f"{centre[0]}", "--center_y", f"{centre[1]}", "--center_z", f"{centre[2]}",
           "--size_x", f"{box}", "--size_y", f"{box}", "--size_z", f"{box}", "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES), "--seed", str(seed),
           "--cpu", "1", "--out", str(out)]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    m = re.search(r"^REMARK VINA RESULT:\s+(-?\d+\.\d+)", out.read_text(encoding="utf-8") if out.exists() else "", re.M)
    return float(m.group(1)) if m else None


def prepare_target(name: str, rec: dict[str, Any]) -> dict[str, Any]:
    d = DOCK / name
    pdb_to_pdbqt(d / "receptor.pdb", d / "receptor.pdbqt", receptor=True)
    pdb_to_pdbqt(d / "crystal_ligand.pdb", d / "crystal_ligand.pdbqt", receptor=False)
    return {"dir": d, "centre": rec["selected"]["centre"], "box": rec["selected"]["box"]}


def redock(name: str, info: dict[str, Any]) -> dict[str, Any]:
    d = info["dir"]
    out = d / "redock_out.pdbqt"
    score = run_vina(d / "receptor.pdbqt", d / "crystal_ligand.pdbqt", info["centre"], info["box"], out)
    crystal, pose = pdbqt_atoms(d / "crystal_ligand.pdbqt"), pdbqt_atoms(out, 1)
    rmsd = float(math.sqrt(((crystal - pose) ** 2).sum(axis=1).mean())) if len(crystal) == len(pose) and len(pose) else None
    return {"best_affinity": score, "rmsd_best_pose": rmsd, "n_heavy_atoms": int(len(crystal)), "passes": rmsd is not None and rmsd <= REDOCK_RMSD_MAX}


def dock_set(name: str, info: dict[str, Any], ligands: dict[str, str], workers: int = 8) -> dict[str, float | None]:
    """Best Vina affinity per ligand id; ligand preparation is per call and temporary."""
    d = info["dir"]

    def one(item: tuple[str, str]) -> tuple[str, float | None]:
        lid, smiles = item
        with tempfile.TemporaryDirectory() as tmp:
            lig, out = Path(tmp) / "l.pdbqt", Path(tmp) / "o.pdbqt"
            if not smiles_to_pdbqt(smiles, lig):
                return lid, None
            return lid, run_vina(d / "receptor.pdbqt", lig, info["centre"], info["box"], out)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        return dict(ex.map(one, ligands.items()))


def auc(scores_active: list[float], scores_decoy: list[float]) -> float:
    """ROC AUC of the Vina affinity (more negative = better): probability that an active scores better than a decoy (ties count half)."""
    wins = sum((a < d) + 0.5 * (a == d) for a in scores_active for d in scores_decoy)
    return wins / (len(scores_active) * len(scores_decoy))


def threshold(scores_active: list[float], scores_decoy: list[float]) -> float | None:
    """Least negative affinity cut-off whose precision on the validation set is >= PRECISION_MIN (needs at least 3 actives under the cut-off)."""
    cuts = sorted(set(scores_active + scores_decoy))
    best = None
    for c in cuts:
        tp, fp = sum(a <= c for a in scores_active), sum(d <= c for d in scores_decoy)
        if tp >= 3 and tp / (tp + fp) >= PRECISION_MIN:
            best = c
    return best


def validation_sets(universe: dict[str, Any], e1: dict[str, Any], uniprot: str) -> tuple[dict[str, str], dict[str, str]]:
    """Actives (curated or measured link to the protein) and a fixed random sample of decoys (other approved small molecules)."""
    ok = {d["drug_id"]: d for d in universe["drugs"] if d.get("smiles") and d.get("molecule_type") == "Small molecule" and heavy_ok(d["smiles"])}
    actives = {i: d["smiles"] for i, d in ok.items() if uniprot in d["human_targets"] or uniprot in e1.get(i, {})}
    rest = sorted(set(ok) - set(actives))
    sample = random.Random(SEED).sample(rest, min(N_DECOYS, len(rest)))
    return actives, {i: ok[i]["smiles"] for i in sample}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("validate", "screen", "top"), required=True)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--targets", nargs="*")
    args = parser.parse_args()
    structures = json.loads((DOCK / "structures.json").read_text(encoding="utf-8"))
    universe = json.loads(Path("data/universe.json").read_text(encoding="utf-8"))
    e1 = json.loads(Path("data/dti_e1_all.json").read_text(encoding="utf-8"))["edges"]
    names = [t for t in structures if structures[t]["selected"] and (not args.targets or t in args.targets)]
    if args.stage == "validate":
        result: dict[str, Any] = {"rules": {"redock_rmsd_max": REDOCK_RMSD_MAX, "auc_min": AUC_MIN, "min_actives": MIN_ACTIVES, "precision_min": PRECISION_MIN,
                                            "exhaustiveness": EXHAUSTIVENESS, "n_decoys": N_DECOYS, "seed": SEED}, "targets": {}}
        for name in names:
            info = prepare_target(name, structures[name])
            rd = redock(name, info)
            actives, decoys = validation_sets(universe, e1, structures[name]["uniprot"])
            entry: dict[str, Any] = {"structure": structures[name]["selected"]["pdb_id"], "redock": rd, "n_actives": len(actives), "n_decoys": len(decoys)}
            print(name, entry["structure"], "redock rmsd", rd["rmsd_best_pose"], "actives", len(actives), flush=True)
            if len(actives) >= MIN_ACTIVES:
                a, dcy = dock_set(name, info, actives, args.workers), dock_set(name, info, decoys, args.workers)
                sa, sd = [v for v in a.values() if v is not None], [v for v in dcy.values() if v is not None]
                entry.update({"auc": auc(sa, sd) if sa and sd else None, "n_actives_docked": len(sa), "n_decoys_docked": len(sd), "active_scores": a, "decoy_scores": dcy,
                              "threshold_affinity": threshold(sa, sd) if sa and sd else None})
            entry["passes_redock"] = rd["passes"]
            entry["passes_enrichment"] = bool(entry.get("auc") is not None and entry["auc"] >= AUC_MIN and entry.get("n_actives_docked", 0) >= MIN_ACTIVES)
            entry["validated"] = bool(entry["passes_redock"] and entry["passes_enrichment"] and entry.get("threshold_affinity") is not None)
            result["targets"][name] = entry
            print("   ", {k: entry.get(k) for k in ("auc", "passes_redock", "passes_enrichment", "threshold_affinity", "validated")}, flush=True)
        (args.out / "docking_validation.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        return
    validation = json.loads((args.out / "docking_validation.json").read_text(encoding="utf-8"))["targets"]
    usable = [t for t in names if validation.get(t, {}).get("validated")]
    if not usable:
        print("No target passed the docking validation: nothing to screen (this is a result, see out/docking_validation.json).")
        return
    if args.stage == "screen":
        ok = {d["drug_id"]: d["smiles"] for d in universe["drugs"] if d.get("smiles") and d.get("molecule_type") == "Small molecule" and heavy_ok(d["smiles"])}
        result = {}
        for name in usable:
            info = prepare_target(name, structures[name])
            result[name] = {"threshold_affinity": validation[name]["threshold_affinity"], "scores": dock_set(name, info, ok, args.workers)}
            print(name, "screened", len(ok), flush=True)
        (args.out / "docking_screen.json").write_text(json.dumps(result, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    else:
        pool = json.loads((args.out / "proposals.json").read_text(encoding="utf-8"))["proposals"]
        smiles = {d["drug_id"]: d["smiles"] for d in universe["drugs"] if d.get("smiles")}
        ligands = {r["drug_id"]: smiles[r["drug_id"]] for r in pool if r["drug_id"] in smiles and heavy_ok(smiles[r["drug_id"]])}
        result = {name: {"threshold_affinity": validation[name]["threshold_affinity"], "scores": dock_set(name, prepare_target(name, structures[name]), ligands, args.workers)} for name in usable}
        (args.out / "docking_top.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        print("top candidates docked:", len(ligands), "on", usable)


if __name__ == "__main__":
    main()
