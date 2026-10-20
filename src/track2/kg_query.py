"""Line A applied to the patient (only valid if `kg_validation` passed): personalised PageRank from the disease node and the patient's HPO phenotype nodes.

    python -m track2.kg_query --out out [--disease "mosaic variegated aneuploidy syndrome"]

Three restart sets are ranked so the source of the signal is visible: the disease node alone, the phenotype nodes alone, and both (the one that decides). A drug is
*supported by A* if it is a drug of the register universe (matched by name) among the `TOP` best of those under the combined restart. Writes `out/kg_query.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from track2.fda_approved import candidate_names
from track2.knowledge_graph import Graph, load_graph, personalised_pagerank, restart_vector
from track2.phenotype_module import CORE_TERMS, NONSPECIFIC_TERMS, query_terms

TOP = 50
DEFAULT_DISEASE = "mosaic variegated aneuploidy syndrome"


def hpo_node(graph: Graph, hpo_id: str) -> int | None:
    """PrimeKG stores the number of the HPO identifier (HP:0002859 -> '2859') as the id of an effect/phenotype node."""
    number = str(int(hpo_id.split(":")[1]))
    hit = np.flatnonzero((graph.node_type == "effect/phenotype") & (graph.node_id == number))
    return int(hit[0]) if len(hit) else None


def universe_match(graph: Graph, universe: dict[str, Any]) -> dict[int, dict[str, str]]:
    """Drug node -> universe drug, by upper-case name or name variant."""
    by_name: dict[str, dict[str, str]] = {}
    for d in universe["drugs"]:
        for variant in candidate_names(d["name"]):
            by_name.setdefault(variant, {"drug_id": d["drug_id"], "name": d["name"]})
    out = {}
    for node in np.flatnonzero(graph.node_type == "drug"):
        for variant in candidate_names(str(graph.node_name[node]).upper()):
            if variant in by_name:
                out[int(node)] = by_name[variant]
                break
    return out


def rank(graph: Graph, restart_nodes: dict[int, float], matched: dict[int, dict[str, str]]) -> list[dict[str, Any]]:
    p = personalised_pagerank(graph.adj, restart_vector(graph.adj.shape[0], restart_nodes)[:, None])[:, 0]
    nodes = sorted(matched, key=lambda n: (-p[n], n))
    return [{"rank": i + 1, "drug_id": matched[n]["drug_id"], "name": matched[n]["name"], "score": float(p[n])} for i, n in enumerate(nodes)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--sources", type=Path, default=Path("data/sources/primekg"))
    parser.add_argument("--patient", type=Path, default=Path("config/patient.json"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--disease", default=DEFAULT_DISEASE)
    parser.add_argument("--e5", action="store_true", help="addendum 11, E5: leave out the parental term; output goes to kg_query_e5*.json")
    parser.add_argument("--leave-out", nargs="*", default=[], help="HPO ids left out (ablation, with --e5)")
    parser.add_argument("--only", choices=["core", "nonspecific"], default=None, help="restrict the phenotype nodes (with --e5)")
    args = parser.parse_args()
    graph = load_graph(args.sources / "nodes.tab", args.sources / "edges.csv")
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    patient = json.loads(args.patient.read_text(encoding="utf-8"))
    terms = query_terms(patient, set(args.leave_out), {"core": CORE_TERMS, "nonspecific": NONSPECIFIC_TERMS}.get(args.only)) if args.e5 else [t["id"] for t in patient["hpo_terms"]]
    disease = np.flatnonzero((graph.node_type == "disease") & (np.char.lower(graph.node_name) == args.disease.lower()))
    if not len(disease):
        raise SystemExit(f"disease node '{args.disease}' not found in the graph")
    d_node = int(disease[0])
    pheno = {t: hpo_node(graph, t) for t in terms}
    p_nodes = {n: 1.0 for n in pheno.values() if n is not None}
    matched = universe_match(graph, universe)
    rankings = {"disease_only": rank(graph, {d_node: 1.0}, matched), "phenotypes_only": rank(graph, p_nodes, matched),
                "combined": rank(graph, {d_node: 1.0, **p_nodes}, matched)}
    top = {k: {r["drug_id"] for r in v[:TOP]} for k, v in rankings.items()}
    result = {"disease_node": {"index": d_node, "name": str(graph.node_name[d_node]), "degree_without_drug_disease_edges": int(graph.adj[d_node].nnz)},
              "phenotype_nodes": {t: ({"index": n, "name": str(graph.node_name[n]), "degree": int(graph.adj[n].nnz)} if n is not None else None) for t, n in pheno.items()},
              "n_drugs_matched_to_universe": len(matched), "top": TOP,
              "overlap_top": {"disease_only_and_phenotypes_only": len(top["disease_only"] & top["phenotypes_only"]),
                              "combined_and_disease_only": len(top["combined"] & top["disease_only"]), "combined_and_phenotypes_only": len(top["combined"] & top["phenotypes_only"])},
              "rankings": {k: v[:TOP * 2] for k, v in rankings.items()}, "supported_by_A": rankings["combined"][:TOP],
              "rank_combined": {r["drug_id"]: r["rank"] for r in rankings["combined"]}}
    name = "kg_query.json"
    if args.e5:
        name = "kg_query_e5" + (f"_{args.only}" if args.only else "") + ("_minus_" + "_".join(x.replace("HP:", "") for x in args.leave_out) if args.leave_out else "") + ".json"
    (args.out / name).write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("disease node degree:", result["disease_node"]["degree_without_drug_disease_edges"], "| phenotype nodes found:", sum(1 for v in pheno.values() if v is not None), "of", len(terms),
          "| matched drugs:", len(matched), "| overlaps:", result["overlap_top"])
    print("top 25 combined:", [r["name"] for r in rankings["combined"][:25]])


if __name__ == "__main__":
    main()
