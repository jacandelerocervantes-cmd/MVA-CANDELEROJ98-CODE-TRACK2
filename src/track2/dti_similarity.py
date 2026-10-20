"""E2 links: drug-target links predicted by ligand similarity, with a retrospective validation (`protocol/addendum_dti.md`).

    python -m track2.dti_similarity --out out

For each protein of the modules with at least 10 reference ligands (ChEMBL, pChEMBL >= 6, the 300 strongest) the score of a drug
is its maximum Tanimoto similarity (Morgan radius 2, 2048 bits) to a reference ligand, excluding the drug itself. tau is the
smallest threshold whose precision is >= 0.5 on known links (E0 + E1) with the drug removed from the references. If no tau reaches
0.5 the layer is not used. Predictions are hypotheses about binding, not measurements.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

from track2.http_cache import get_json
from track2.dti_edges import BASE, PCHEMBL_MIN, WORKERS

RDLogger.DisableLog("rdApp.*")
MIN_LIGANDS, MAX_LIGANDS, PRECISION_MIN = 10, 300, 0.5
TAUS = [round(0.30 + 0.025 * i, 3) for i in range(0, 25)]  # 0.30 ... 0.90


def fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) if mol else None


def uniprot_to_target(accessions: list[str], get=get_json) -> dict[str, str]:
    """UniProt -> ChEMBL single-protein target id (batches of 50)."""
    out: dict[str, str] = {}
    for i in range(0, len(accessions), 50):
        chunk = accessions[i:i + 50]
        data = get(f"{BASE}/target.json?target_components__accession__in={','.join(chunk)}&target_type=SINGLE%20PROTEIN&limit=1000")
        for t in data["targets"]:
            comps = t.get("target_components") or []
            if len(comps) == 1 and comps[0].get("accession") in chunk:
                out[comps[0]["accession"]] = t["target_chembl_id"]
    return out


def reference_ligands(target_id: str, get=get_json) -> list[tuple[str, str]]:
    """(molecule_chembl_id, smiles) of the strongest ligands (pChEMBL >= 6, exact relation, binding assay)."""
    url = (f"{BASE}/activity.json?target_chembl_id={target_id}&pchembl_value__gte={PCHEMBL_MIN}&assay_type=B&standard_relation=%3D"
           f"&order_by=-pchembl_value&limit={MAX_LIGANDS}")
    seen: dict[str, str] = {}
    for a in get(url)["activities"]:
        if a.get("canonical_smiles") and a["molecule_chembl_id"] not in seen:
            seen[a["molecule_chembl_id"]] = a["canonical_smiles"]
    return list(seen.items())


def score_matrix(drug_fps: list[Any], drug_ids: list[str], ligands: dict[str, list[tuple[str, Any]]]) -> dict[str, np.ndarray]:
    """Per target: max Tanimoto of each drug to the target's ligands (the drug's own entry excluded)."""
    out = {}
    for target, ligs in ligands.items():
        ids = [m for m, _ in ligs]
        fps = [f for _, f in ligs]
        scores = np.zeros(len(drug_fps))
        for i, fp in enumerate(drug_fps):
            if fp is None:
                scores[i] = np.nan
                continue
            sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, fps))
            for j, m in enumerate(ids):
                if m == drug_ids[i]:
                    sims[j] = -1.0
            scores[i] = sims.max()
        out[target] = scores
    return out


def validate(scores: dict[str, np.ndarray], positives: set[tuple[int, str]], n_drugs: int) -> dict[str, Any]:
    """Pooled precision/recall per threshold over all (drug, target) pairs of the targets tested."""
    y, s = [], []
    for target, sc in scores.items():
        for i in range(n_drugs):
            if np.isnan(sc[i]):
                continue
            y.append((i, target) in positives)
            s.append(sc[i])
    y_arr, s_arr = np.array(y), np.array(s)
    rows = []
    for tau in TAUS:
        sel = s_arr >= tau
        tp, fp = int((y_arr & sel).sum()), int((~y_arr & sel).sum())
        rows.append({"tau": tau, "predicted": tp + fp, "true_positive": tp, "precision": tp / (tp + fp) if tp + fp else None,
                     "recall": tp / int(y_arr.sum()) if y_arr.sum() else None})
    ok = [r for r in rows if r["precision"] is not None and r["precision"] >= PRECISION_MIN and r["true_positive"] >= 5]
    return {"n_pairs": int(len(y_arr)), "n_positive": int(y_arr.sum()), "base_rate": float(y_arr.mean()), "by_tau": rows,
            "tau": ok[0]["tau"] if ok else None, "valid": bool(ok)}


def heldout(scores: dict[str, np.ndarray], positives: set[tuple[int, str]], n_drugs: int, repeats: int = 20, seed: int = 7) -> dict[str, Any]:
    """Out-of-sample check: split the targets in two halves at random, choose tau on one half by the same rule and measure the precision at that tau on the other half.

    The in-sample tau of `validate` is the first threshold with precision >= 0.5 on the same pairs, so it is optimistic; this gives the precision one would get on targets
    that played no part in choosing it.
    """
    rng = np.random.default_rng(seed)
    names = sorted(scores)
    precisions, taus, tp_total, n_pred_total = [], [], 0, 0
    for _ in range(repeats):
        order = rng.permutation(len(names))
        half = len(names) // 2
        a = {names[i]: scores[names[i]] for i in order[:half]}
        b = {names[i]: scores[names[i]] for i in order[half:]}
        pa = {(i, t) for (i, t) in positives if t in a}
        pb = {(i, t) for (i, t) in positives if t in b}
        va = validate(a, pa, n_drugs)
        if not va["valid"]:
            continue
        tau = va["tau"]
        taus.append(tau)
        tp = fp = 0
        for t, sc in b.items():
            for i in range(n_drugs):
                if not np.isnan(sc[i]) and sc[i] >= tau:
                    tp += (i, t) in pb
                    fp += (i, t) not in pb
        if tp + fp:
            precisions.append(tp / (tp + fp))
            tp_total += tp
            n_pred_total += tp + fp
    return {"repeats": repeats, "n_valid_splits": len(taus), "tau_mean": float(np.mean(taus)) if taus else None,
            "heldout_precision_mean": float(np.mean(precisions)) if precisions else None,
            "heldout_precision_p10_p90": [float(np.percentile(precisions, 10)), float(np.percentile(precisions, 90))] if precisions else None,
            "heldout_pooled_precision": tp_total / n_pred_total if n_pred_total else None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--e1", type=Path, default=Path("data/dti_e1_all.json"))
    parser.add_argument("--module-edges", type=Path, default=Path("out/dti_e1_module_edges.json"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    e1 = json.loads(args.e1.read_text(encoding="utf-8"))["edges"]
    from track2.dti_edges import load_interactome, module_uniprots
    interactome = load_interactome(700)
    module_uniprot = module_uniprots(Path("config/mechanisms_from_submission"), interactome)
    mapping = uniprot_to_target(sorted(module_uniprot))
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        ligs_raw = dict(zip(mapping, ex.map(lambda u: reference_ligands(mapping[u]), mapping)))
    ligands = {u: [(m, fingerprint(smi)) for m, smi in v] for u, v in ligs_raw.items() if len(v) >= MIN_LIGANDS}
    ligands = {u: [(m, f) for m, f in v if f is not None] for u, v in ligands.items()}
    ligands = {u: v for u, v in ligands.items() if len(v) >= MIN_LIGANDS}
    drugs = [d for d in universe["drugs"] if d.get("smiles") and d.get("molecule_type") == "Small molecule"]
    drug_ids = [d["drug_id"] for d in drugs]
    fps = [fingerprint(d["smiles"]) for d in drugs]
    scores = score_matrix(fps, drug_ids, ligands)
    index = {d: i for i, d in enumerate(drug_ids)}
    positives: set[tuple[int, str]] = set()
    for d in drugs:
        for u in set(d["human_targets"]) | set(e1.get(d["drug_id"], {})):
            if u in ligands:
                positives.add((index[d["drug_id"]], u))
    validation = validate(scores, positives, len(drugs))
    validation["heldout_by_target"] = heldout(scores, positives, len(drugs))
    validation.update({"n_targets_with_ligands": len(ligands), "n_module_proteins": len(module_uniprot), "precision_min": PRECISION_MIN})
    edges: dict[str, dict[str, float]] = {}
    if validation["valid"]:
        tau = validation["tau"]
        for u, sc in scores.items():
            for i in np.where(sc >= tau)[0]:
                known = u in set(drugs[i]["human_targets"]) | set(e1.get(drug_ids[i], {}))
                if not known:  # only links that are not already curated or measured
                    edges.setdefault(drug_ids[i], {})[u] = round(float(sc[i]), 3)
    (args.out / "dti_e2_validation.json").write_text(json.dumps(validation, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    (args.out / "dti_e2_edges.json").write_text(json.dumps({"tau": validation["tau"], "edges": edges}, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print({k: v for k, v in validation.items() if k != "by_tau"}, "| drugs with new predicted links:", len(edges))


if __name__ == "__main__":
    main()
