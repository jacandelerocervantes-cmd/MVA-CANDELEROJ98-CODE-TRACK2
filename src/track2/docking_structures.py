"""Choose and prepare one PDB structure per docking target by the rule of `protocol/addendum_docking.md`.

    python -m track2.docking_structures

For each target: among the human X-ray entries of the RCSB PDB for the UniProt accession with resolution <= 2.8 angstrom, take the best-resolution one that has a
drug-like ligand (15 to 60 heavy atoms, not a metal, solvent or buffer) bound in the protein chain. Writes the receptor (chain without waters, ligands or alternative
locations) and the co-crystal ligand to `data/docking/<target>/`, and the search box (centre of the ligand, 22 angstrom). Source: RCSB PDB (CC0).
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import numpy as np

TARGETS = {"CDK1": "P06493", "AURKA": "O14965", "AURKB": "Q96GD4", "PLK1": "P53350", "PLK4": "O00444", "TTK": "P33981", "NEK7": "Q8TDX7"}
OUT = Path("data/docking")
MAX_RES = 2.8
BOX = 22.0
IGNORED = {"HOH", "DOD", "GOL", "EDO", "PEG", "PGE", "PG4", "P6G", "1PE", "DMS", "SO4", "PO4", "ACT", "FMT", "CL", "NA", "MG", "ZN", "CA", "K", "MN", "IOD", "BME", "TRS", "MES", "EPE", "CIT", "TLA", "MPD", "IMD", "NO3", "FLC", "ACY", "EOH", "NH4", "BR", "SCN", "CD", "NI", "CO", "CU", "FE", "SEP", "TPO", "PTR"}
UA = {"User-Agent": "hselith-track2/0.1 (research)"}


def _post(url: str, payload: dict[str, Any]) -> Any:
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json", **UA})
    return json.loads(urllib.request.urlopen(req, timeout=120).read())


def candidate_entries(uniprot: str) -> list[str]:
    """PDB ids with the protein, human, X-ray, resolution <= MAX_RES, best resolution first."""
    query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator": "exact_match", "value": uniprot}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": MAX_RES}}]},
        "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 40}, "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]}}
    return [r["identifier"] for r in _post("https://search.rcsb.org/rcsbsearch/v2/query", query).get("result_set", [])]


def parse_pdb(text: str) -> tuple[list[str], list[str]]:
    atoms = [l for l in text.splitlines() if l.startswith("ATOM")]
    hets = [l for l in text.splitlines() if l.startswith("HETATM")]
    return atoms, hets


def druglike_ligand(hets: list[str], chain: str) -> tuple[str, list[str]] | None:
    """The largest drug-like HETATM group (by heavy atoms) that sits in contact with the chosen chain; chain contact is approximated by the same chain id."""
    groups: dict[tuple[str, str, str], list[str]] = {}
    for l in hets:
        name, ch, res = l[17:20].strip(), l[21], l[22:27].strip()
        if name in IGNORED or ch != chain or l[16] not in (" ", "A"):
            continue
        groups.setdefault((name, ch, res), []).append(l)
    best = None
    for (name, _, _), lines in groups.items():
        heavy = [l for l in lines if l[76:78].strip() not in ("H", "D")]
        if 15 <= len(heavy) <= 60 and (best is None or len(heavy) > len(best[1])):
            best = (name, heavy)
    return best


def fetch_pdb(pdb_id: str, retries: int = 5) -> str | None:
    """PDB-format text of an entry, cached in `data/docking/_pdb/`; None if the entry has no PDB-format file (404)."""
    cache = OUT / "_pdb" / f"{pdb_id}.pdb"
    if cache.exists():
        return cache.read_text(encoding="utf-8", errors="replace")
    for attempt in range(retries):
        try:
            text = urllib.request.urlopen(urllib.request.Request(f"https://files.rcsb.org/download/{pdb_id}.pdb", headers=UA), timeout=120).read().decode("utf-8", "replace")
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(text, encoding="utf-8", newline=chr(10))
            return text
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return None
            time.sleep(3 + 3 * attempt)
        except urllib.error.URLError:
            time.sleep(3 + 3 * attempt)
    raise RuntimeError(f"Could not download {pdb_id} after {retries} attempts")


def select(uniprot: str) -> dict[str, Any] | None:
    for pdb_id in candidate_entries(uniprot):
        text = fetch_pdb(pdb_id)
        if text is None:  # very large entries exist only as mmCIF: skipped (the PDB-format parser cannot read them)
            continue
        atoms, hets = parse_pdb(text)
        dbref = [l for l in text.splitlines() if l.startswith("DBREF") and uniprot in l]
        chains = sorted({l[12] for l in dbref}) or sorted({l[21] for l in atoms})
        for chain in chains:
            lig = druglike_ligand(hets, chain)
            if lig:
                receptor = [l for l in atoms if l[21] == chain and l[16] in (" ", "A")]
                coords = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lig[1]])
                res_line = next((l for l in text.splitlines() if l.startswith("REMARK   2 RESOLUTION")), "")
                return {"pdb_id": pdb_id, "chain": chain, "ligand": lig[0], "n_ligand_heavy_atoms": len(lig[1]), "centre": [round(float(x), 3) for x in coords.mean(axis=0)],
                        "box": BOX, "resolution_line": res_line.strip(), "receptor_lines": receptor, "ligand_lines": lig[1]}
    return None


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    record = {}
    for target, uniprot in TARGETS.items():
        sel = select(uniprot)
        if not sel:
            record[target] = {"uniprot": uniprot, "selected": None, "reason": "no X-ray entry <= 2.8 angstrom with a drug-like ligand in the chain"}
            print(target, "no structure")
            continue
        d = OUT / target
        d.mkdir(exist_ok=True)
        (d / "receptor.pdb").write_text("\n".join(sel.pop("receptor_lines")) + "\nEND\n", encoding="utf-8", newline="\n")
        (d / "crystal_ligand.pdb").write_text("\n".join(sel.pop("ligand_lines")) + "\nEND\n", encoding="utf-8", newline="\n")
        record[target] = {"uniprot": uniprot, "selected": sel}
        print(target, sel["pdb_id"], sel["chain"], sel["ligand"], sel["n_ligand_heavy_atoms"], sel["centre"], sel["resolution_line"][-25:])
    (OUT / "structures.json").write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
