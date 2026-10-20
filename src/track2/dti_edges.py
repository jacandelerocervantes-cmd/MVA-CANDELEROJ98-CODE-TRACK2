"""E1 links: measured ChEMBL bioactivity of approved small-molecule drugs on human single-protein targets (`protocol/addendum_dti.md`).

    python -m track2.dti_edges --universe data/universe.json --out out

For every approved drug with a structure, activities with pChEMBL >= 6.0 (<= 1 uM), exact relation, binding assay, human target,
are collected and mapped to UniProt. All targets are cached in `data/dti_e1_all.json`; the public output keeps only the proteins
of the submission-derived modules (`out/dti_e1_module_edges.json`). Nothing is predicted here.
"""

from __future__ import annotations

import argparse
import json
import statistics
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from track2.http_cache import get_json
from track2.modules import build_module, load_mechanism
from track2.string_net import load_interactome

BASE = "https://www.ebi.ac.uk/chembl/api/data"
PCHEMBL_MIN = 6.0
TYPES = "IC50,Ki,Kd,EC50,Potency"
WORKERS = 4


def drug_activities(drug_id: str, get=get_json) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    offset = 0
    while True:
        url = (f"{BASE}/activity.json?molecule_chembl_id={drug_id}&pchembl_value__gte={PCHEMBL_MIN}&assay_type=B"
               f"&target_organism=Homo%20sapiens&standard_relation=%3D&standard_type__in={TYPES}&limit=1000&offset={offset}")
        data = get(url)
        out.extend(data["activities"])
        if not data["page_meta"]["next"]:
            return out
        offset += 1000


def summarise(activities: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    per: dict[str, list[float]] = {}
    for a in activities:
        try:
            value = float(a["pchembl_value"])
            target = a["target_chembl_id"]
        except (KeyError, TypeError, ValueError):
            continue
        per.setdefault(target, []).append(value)
    return {t: {"n": len(v), "median_pchembl": round(statistics.median(v), 2)} for t, v in per.items()}


def target_to_uniprot(target_ids: list[str], get=get_json) -> dict[str, str]:
    """ChEMBL single-protein targets -> UniProt accession (batches of 100)."""
    out: dict[str, str] = {}
    ids = sorted(set(target_ids))
    for i in range(0, len(ids), 100):
        data = get(f"{BASE}/target.json?target_chembl_id__in={','.join(ids[i:i + 100])}&limit=1000")
        for t in data["targets"]:
            comps = t.get("target_components") or []
            if t.get("target_type") == "SINGLE PROTEIN" and len(comps) == 1 and comps[0].get("accession"):
                out[t["target_chembl_id"]] = comps[0]["accession"]
    return out


def module_nodes(mechanism_dir: Path, interactome) -> set[int]:
    nodes: set[int] = set()
    for p in sorted(mechanism_dir.glob("*.json")):
        mech = load_mechanism(p)
        if mech.get("status") == "mechanism_candidate":
            nodes |= {int(n) for n in build_module(mech, interactome)["nodes"]}
    return nodes


def module_uniprots(mechanism_dir: Path, interactome) -> set[str]:
    """Every UniProt accession that maps to a node of a module (several accessions can share one STRING node)."""
    by_node: dict[int, set[str]] = {}
    for u, n in interactome.uniprot_to_node.items():
        by_node.setdefault(n, set()).add(u)
    out: set[str] = set()
    for p in sorted(mechanism_dir.glob("*.json")):
        mech = load_mechanism(p)
        if mech.get("status") == "mechanism_candidate":
            for n in build_module(mech, interactome)["nodes"]:
                out |= by_node.get(int(n), set())
    return out


def collect(universe: dict[str, Any], get=get_json, workers: int = WORKERS) -> dict[str, Any]:
    drugs = [d for d in universe["drugs"] if d.get("smiles") and d.get("molecule_type") == "Small molecule"]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        acts = dict(zip((d["drug_id"] for d in drugs), ex.map(lambda d: summarise(drug_activities(d["drug_id"], get)), drugs)))
    mapping = target_to_uniprot([t for v in acts.values() for t in v], get)
    edges = {d: {mapping[t]: v for t, v in tv.items() if t in mapping} for d, tv in acts.items()}
    return {"n_drugs_queried": len(drugs), "chembl_version": universe["chembl_version"], "pchembl_min": PCHEMBL_MIN, "edges": edges}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    result = collect(universe)
    Path("data").mkdir(exist_ok=True)
    Path("data/dti_e1_all.json").write_text(json.dumps(result, sort_keys=True), encoding="utf-8", newline="\n")
    interactome = load_interactome(700)
    module_uniprot = module_uniprots(args.mechanism_dir, interactome)
    kept = {d: {u: v for u, v in e.items() if u in module_uniprot} for d, e in result["edges"].items()}
    kept = {d: e for d, e in kept.items() if e}
    summary = {"n_drugs_queried": result["n_drugs_queried"], "n_module_nodes": len(module_nodes(args.mechanism_dir, interactome)), "n_module_accessions": len(module_uniprot), "n_drugs_with_module_link": len(kept),
               "n_links": sum(len(e) for e in kept.values()), "n_drugs_any_target": sum(1 for e in result["edges"].values() if e)}
    (args.out / "dti_e1_module_edges.json").write_text(json.dumps({"summary": summary, "edges": kept}, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(summary)


if __name__ == "__main__":
    main()
