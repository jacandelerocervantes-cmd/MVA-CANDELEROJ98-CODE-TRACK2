"""New repositioning hypotheses from the expanded drug-target links (rules: `protocol/addendum_dti.md`).

    python -m track2.novel_hypotheses --variant e0_e1 --out out

A drug is a **new** hypothesis if, with the expanded links, z <= -2 in at least one module at STRING 700 and it either was not scorable
before or its curated-only z was > -2. Robustness uses the nine configurations of `addendum_robustness.md` on the expanded ranking.
The pool also keeps non-new drugs that are in the top 20 of at least 5 configurations (flagged `is_new = false`). Each row lists the module proteins it is linked to and where each link comes from (curated, measured, predicted).
Not a recommendation: whether any drug is appropriate for the child is a clinical decision.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path
from typing import Any

from track2.consensus import SCHEMES, THRESHOLDS, TOP_N, ROBUST_MIN, FRAGILE_MAX, combined, top, weights
from track2.labels import label_info
from track2.literature_context import classify
from track2.modules import build_module, load_mechanism
from track2.pediatric_review import load_review, reviewed_profile
from track2.proximity import min_distance_to_set
from track2.string_net import load_interactome

PROMISCUOUS_MODULE_TARGETS = 10
BROAD_PROFILE_TARGETS = 20  # a drug with more measured or curated targets than this is flagged: its closeness may reflect breadth, not specificity
EXCLUDED_ATC = ("N01", "N05C")
TUMOUR_ATC = ("L01", "L02")


def drug_class(atc: list[str]) -> str:
    if any(a.startswith(EXCLUDED_ATC) for a in atc):
        return "excluded_acute_use_anaesthetic_sedative"
    if any(a.startswith(TUMOUR_ATC) for a in atc):
        return "antineoplastic_tumour_context"
    return "constitutional_hypothesis"


def module_proteins(mechanisms: list[dict[str, Any]], interactome) -> dict[str, set[str]]:
    by_node: dict[int, set[str]] = {}
    for u, n in interactome.uniprot_to_node.items():
        by_node.setdefault(n, set()).add(u)  # several accessions can share one STRING node
    return {m["id"]: {u for n in build_module(m, interactome)["nodes"] for u in by_node.get(int(n), set())} for m in mechanisms}


def links(drug: dict[str, Any], e1: dict[str, Any], e2: dict[str, float], module_uni: set[str], symbols: dict[str, str],
          actions: dict[str, str]) -> list[dict[str, Any]]:
    out = []
    for u in sorted(module_uni & (set(drug["human_targets"]) | set(e1) | set(e2))):
        source = "curated" if u in drug["human_targets"] else "measured" if u in e1 else "predicted"
        out.append({"uniprot": u, "gene": symbols.get(u, u), "source": source, "action": actions.get(u),
                    "median_pchembl": (e1.get(u) or {}).get("median_pchembl"), "n_measurements": (e1.get(u) or {}).get("n"), "similarity": e2.get(u)})
    return out


def curated_actions(drug: dict[str, Any], targets: dict[str, Any]) -> dict[str, str]:
    """UniProt -> ChEMBL action type (INHIBITOR, AGONIST, ...) from the curated mechanism table."""
    out: dict[str, str] = {}
    for m in drug["mechanisms"]:
        for u in (targets.get(m.get("target_id")) or {}).get("uniprot", []):
            if m.get("action"):
                out[u] = m["action"]
    return out


def proximal_targets(drug: dict[str, Any], e1: dict[str, Any], interactome, distance: dict[str, Any], symbols: dict[str, str], keep: int = 5) -> list[dict[str, Any]]:
    """The drug's targets closest to any module (STRING distance in steps), with their source and measured pChEMBL."""
    rows = []
    for u in set(drug["human_targets"]) | set(e1):
        node = interactome.uniprot_to_node.get(u)
        if node is None:
            continue
        mid, d = min(((m, int(dist[node])) for m, dist in distance.items()), key=lambda x: x[1])
        rows.append({"uniprot": u, "gene": symbols.get(u, u), "distance": d, "module": mid.removeprefix("sub_"),
                     "source": "curated" if u in drug["human_targets"] else "measured", "median_pchembl": (e1.get(u) or {}).get("median_pchembl"), "n_measurements": (e1.get(u) or {}).get("n")})
    return sorted(rows, key=lambda r: (r["distance"], r["gene"]))[:keep]


def run(out: Path, mechanism_dir: Path, universe: dict[str, Any], e1_all: dict[str, Any], e2: dict[str, Any], variant: str,
        config: dict[str, Any], review: dict[str, Any], get_label=label_info, top_n: int = TOP_N, pool_min: int = 5) -> dict[str, Any]:
    mechanisms = [m for m in (load_mechanism(p) for p in sorted(mechanism_dir.glob("*.json"))) if m.get("status") == "mechanism_candidate"]
    base = {m["id"]: json.loads((out / f"threshold_scores_{m['id']}.json").read_text(encoding="utf-8")) for m in mechanisms}
    aug = {m["id"]: json.loads((out / f"dti_threshold_scores_{variant}_{m['id']}.json").read_text(encoding="utf-8")) for m in mechanisms}
    wts = weights(mechanisms)
    table = {(t, s): combined({mid: aug[mid][t] for mid in aug}, wts[s]) for t, s in itertools.product(THRESHOLDS, SCHEMES)}
    tops = {k: top(v, top_n) for k, v in table.items()}
    counts: dict[str, int] = {}
    for lst in tops.values():
        for d in lst:
            counts[d] = counts.get(d, 0) + 1
    interactome = load_interactome(700)
    mod_uni = module_proteins(mechanisms, interactome)
    node_symbol: dict[int, str] = {}
    for sym, n in interactome.symbol_to_node.items():
        node_symbol.setdefault(n, sym)
    symbols = {u: node_symbol[n] for u, n in interactome.uniprot_to_node.items() if n in node_symbol}
    drugs = {d["drug_id"]: d for d in universe["drugs"]}
    edges1 = e1_all["edges"]
    mod_nodes = {m["id"]: build_module(m, interactome)["nodes"] for m in mechanisms}
    distance = {mid: min_distance_to_set(interactome.adjacency, nodes) for mid, nodes in mod_nodes.items()}
    rows = []
    for d_id, drug in drugs.items():
        z_aug = {mid: aug[mid]["700"].get(d_id) for mid in aug}
        z_base = {mid: base[mid]["700"].get(d_id) for mid in base}
        za = [v for v in z_aug.values() if v is not None]
        zb = [v for v in z_base.values() if v is not None]
        if not za:
            continue
        is_new = min(za) <= -2.0 and not (zb and min(zb) <= -2.0)  # close only thanks to the expanded links (or not scorable before)
        if not is_new and counts.get(d_id, 0) < pool_min:
            continue  # neither new nor in the top-N of at least `pool_min` configurations
        actions = curated_actions(drug, universe.get("targets", {}))
        all_links = {mid: links(drug, edges1.get(d_id, {}), (e2.get("edges") or {}).get(d_id, {}), mod_uni[mid], symbols, actions) for mid in mod_uni}
        n_linked = len({x["uniprot"] for v in all_links.values() for x in v})
        proximal = proximal_targets(drug, edges1.get(d_id, {}), interactome, distance, symbols)
        n = counts.get(d_id, 0)
        rows.append({
            "drug_id": d_id, "name": drug["name"], "atc": drug["atc"], "class": drug_class(drug["atc"]), "is_new": is_new,
            "black_box_warning": bool(drug.get("black_box_warning")),
            "scorable_before": bool(zb), "curated_best_z": round(min(zb), 2) if zb else None, "expanded_best_z": round(min(za), 2),
            "z_by_module_expanded": {k: (round(v, 2) if v is not None else None) for k, v in z_aug.items()},
            "robustness_count": n, "of": len(tops), "robustness": "robust" if n >= ROBUST_MIN else "fragile" if n <= FRAGILE_MAX else "intermediate",
            "n_module_proteins_linked": n_linked, "promiscuous": n_linked > PROMISCUOUS_MODULE_TARGETS, "n_targets_expanded": len(set(drug["human_targets"]) | set(edges1.get(d_id, {}))),
            "broad_target_profile": len(set(drug["human_targets"]) | set(edges1.get(d_id, {}))) > BROAD_PROFILE_TARGETS, "proximal_targets": proximal,
            "links": all_links, "oral": drug.get("oral"), "parenteral": drug.get("parenteral"), "topical": drug.get("topical"),
            "mechanism_of_action": sorted({m["mechanism"] for m in drug["mechanisms"] if m["mechanism"]}),
        })
    rows.sort(key=lambda r: (-r["robustness_count"], r["expanded_best_z"], r["drug_id"]))
    # labels (network, cached) only for the ones that could be proposed
    for r in rows:
        if r["class"] == "excluded_acute_use_anaesthetic_sedative" or r["robustness_count"] < 5:
            r["label"] = None
            continue
        lab = get_label(r["name"] or "")
        entry = review.get(r["drug_id"])
        profile, source = reviewed_profile(lab["age_profile"], entry)
        min_age = entry["min_age_years"] if entry and entry["verdict"] == "established" else lab["pediatric"]["min_age_years"]
        r["label"] = {"found": lab["found"], "pediatric_status": lab["pediatric"]["status"], "min_age_years": min_age, "verdict": (entry or {}).get("verdict"), "caution": (entry or {}).get("caution"),
                      "bands": profile, "source": source, "excerpt": lab["pediatric"]["excerpt"][:300],
                      "excipient_alerts": [a["excipient"] for a in lab.get("excipient_alerts", [])]}
        r["literature_context"] = [{"class": e["class"], "direction_hint": e["direction_hint"], "pmids": e["pmids"]}
                                   for e in classify({"mechanism_of_action": r["mechanism_of_action"], "risk_flags": ["ANTINEOPLASTIC_AGENT"] if r["class"].startswith("antineoplastic") else []}, config)]
    summary = {"variant": variant, "n_pool": len(rows), "n_new_hypotheses": sum(1 for r in rows if r["is_new"]),
               "n_robust_new": sum(1 for r in rows if r["is_new"] and r["robustness"] == "robust"),
               "n_new_not_scorable_before": sum(1 for r in rows if r["is_new"] and not r["scorable_before"]),
               "by_class_new": {c: sum(1 for r in rows if r["is_new"] and r["class"] == c) for c in ("constitutional_hypothesis", "antineoplastic_tumour_context", "excluded_acute_use_anaesthetic_sedative")},
               "n_promiscuous": sum(1 for r in rows if r["promiscuous"])}
    return {"summary": summary, "rule": {"configurations": len(tops), "top_n": top_n, "pool_min": pool_min, "robust_min": ROBUST_MIN}, "drugs": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=("e0_e1", "e0_e1_e2"), required=True)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--top-n", type=int, default=TOP_N, help="rank cut-off of the proximity lines (20 strict; 50 for the wider ring of addendum_wide.md)")
    parser.add_argument("--pool-min", type=int, default=5)
    parser.add_argument("--tag", default="", help="suffix of the output files, for example _wide")
    args = parser.parse_args()
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    e1_all = json.loads(Path("data/dti_e1_all.json").read_text(encoding="utf-8"))
    e2_path = args.out / "dti_e2_edges.json"
    e2 = json.loads(e2_path.read_text(encoding="utf-8")) if args.variant == "e0_e1_e2" and e2_path.exists() else {"edges": {}}
    config = json.loads(Path("config/literature_context.json").read_text(encoding="utf-8"))
    result = run(args.out, args.mechanism_dir, universe, e1_all, e2, args.variant, config, load_review(), top_n=args.top_n, pool_min=args.pool_min)
    (args.out / f"novel_hypotheses_{args.variant}{args.tag}.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(args.out / f"novel_hypotheses_{args.variant}{args.tag}.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["drug_id", "name", "is_new", "class", "robustness", "robustness_count", "scorable_before", "curated_best_z", "expanded_best_z",
                         "n_module_proteins_linked", "promiscuous", "oral", "parenteral", "topical", "pediatric_status", "min_age_years"])
        for r in result["drugs"]:
            lab = r.get("label") or {}
            writer.writerow([r["drug_id"], r["name"], r["is_new"], r["class"], r["robustness"], r["robustness_count"], r["scorable_before"], r["curated_best_z"],
                             r["expanded_best_z"], r["n_module_proteins_linked"], r["promiscuous"], r["oral"], r["parenteral"], r["topical"],
                             lab.get("pediatric_status", ""), lab.get("min_age_years") if lab.get("min_age_years") is not None else ""])
    print(result["summary"])


if __name__ == "__main__":
    main()
