"""Specificity of the patient query in line A (protocol/addendum_four_lines.md, "Three further analyses", item 2).

    python -m track2.kg_nullquery --out out [--queries 500]

Random queries with the same shape as the patient's: one disease node whose degree is within 30 % of the real disease node's and one phenotype node per real phenotype node with degree
within 30 % of it. Same graph, restart weights and ranking among the drugs matched to the register universe. Reports how often fostamatinib is in the top 5 and top 50 of a random query,
and, for each drug of the real top 50, the fraction of random queries that also put it in their top 50 (generic drugs have a high fraction). Writes `out/kg_nullquery.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from track2.kg_query import DEFAULT_DISEASE, TOP, hpo_node, universe_match
from track2.knowledge_graph import Graph, load_graph, personalised_pagerank, restart_vector
from track2.phenotype_module import query_terms

SEED, TOLERANCE, CHUNK = 20261007, 0.30, 25


def matched_pool(graph: Graph, node_type: str, degree: float, exclude: set[int]) -> np.ndarray:
    """Nodes of one type whose degree is within TOLERANCE of `degree`, without `exclude`."""
    deg = np.asarray(graph.adj.sum(axis=1)).ravel()
    ok = (graph.node_type == node_type) & (deg >= degree * (1 - TOLERANCE)) & (deg <= degree * (1 + TOLERANCE))
    return np.array([n for n in np.flatnonzero(ok) if int(n) not in exclude])


def draw_queries(graph: Graph, disease: int, phenos: list[int], n: int, seed: int = SEED) -> list[dict[int, float]]:
    rng = np.random.default_rng(seed)
    deg = np.asarray(graph.adj.sum(axis=1)).ravel()
    d_pool = matched_pool(graph, "disease", float(deg[disease]), {disease})
    p_pools = [matched_pool(graph, "effect/phenotype", float(deg[p]), set(phenos)) for p in phenos]
    if len(d_pool) == 0 or any(len(p) == 0 for p in p_pools):
        raise SystemExit("no degree-matched nodes for the null queries")
    queries = []
    for _ in range(n):
        q = {int(rng.choice(d_pool)): 1.0}
        for pool in p_pools:  # one phenotype node per real one, all distinct, so the query has the same number of nodes as the patient's
            free = [int(x) for x in pool if int(x) not in q]
            q[int(rng.choice(free))] = 1.0
        queries.append(q)
    return queries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--sources", type=Path, default=Path("data/sources/primekg"))
    parser.add_argument("--patient", type=Path, default=Path("config/patient.json"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--queries", type=int, default=500)
    parser.add_argument("--e5", action="store_true", help="addendum 11, E5: leave out the parental term; output goes to kg_nullquery_e5.json")
    args = parser.parse_args()
    graph = load_graph(args.sources / "nodes.tab", args.sources / "edges.csv")
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    matched = universe_match(graph, universe)
    m_nodes = np.array(sorted(matched))
    disease = int(np.flatnonzero((graph.node_type == "disease") & (np.char.lower(graph.node_name) == DEFAULT_DISEASE))[0])
    patient = json.loads(args.patient.read_text(encoding="utf-8"))
    terms = query_terms(patient) if args.e5 else [t["id"] for t in patient["hpo_terms"]]
    phenos = [p for p in (hpo_node(graph, t) for t in terms) if p is not None]
    n = graph.adj.shape[0]
    real = personalised_pagerank(graph.adj, restart_vector(n, {disease: 1.0, **{p: 1.0 for p in phenos}})[:, None])[:, 0]
    real_top = [matched[int(m_nodes[i])]["drug_id"] for i in np.argsort(-real[m_nodes], kind="stable")[:TOP]]
    fosta_nodes = [i for i, node in enumerate(m_nodes) if matched[int(node)]["name"] == "FOSTAMATINIB"]
    queries = draw_queries(graph, disease, phenos, args.queries)
    top50_count: dict[str, int] = {}
    fosta_ranks = []
    for start in range(0, len(queries), CHUNK):
        part = queries[start:start + CHUNK]
        scores = personalised_pagerank(graph.adj, np.column_stack([restart_vector(n, q) for q in part]))
        for k in range(len(part)):
            order = np.argsort(-scores[m_nodes, k], kind="stable")
            for i in order[:TOP]:
                d = matched[int(m_nodes[i])]["drug_id"]
                top50_count[d] = top50_count.get(d, 0) + 1
            pos = np.empty(len(order), dtype=int)
            pos[order] = np.arange(1, len(order) + 1)
            fosta_ranks.append(int(min(pos[i] for i in fosta_nodes)))
    ranks = np.array(fosta_ranks)
    names = {d["drug_id"]: d["name"] for d in universe["drugs"]}
    frac = {d: top50_count.get(d, 0) / len(queries) for d in dict.fromkeys(real_top)}
    result: dict[str, Any] = {"queries": len(queries), "tolerance": TOLERANCE, "n_phenotype_nodes": len(phenos), "real_top50_unique_drugs": len(frac),
                              "fostamatinib": {"real_rank_note": "see kg_query.json", "random_top5": float((ranks <= 5).mean()), "random_top50": float((ranks <= 50).mean()),
                                               "median_rank_random": float(np.median(ranks)), "empirical_p_rank_le_2": float((1 + (ranks <= 2).sum()) / (1 + len(ranks)))},
                              "fraction_of_random_queries_with_drug_in_top50": {names[d]: round(f, 3) for d, f in sorted(frac.items(), key=lambda kv: -kv[1])},
                              "n_real_top50_drugs_generic_at_least_5_percent": int(sum(1 for f in frac.values() if f >= 0.05)),
                              "n_real_top50_drugs_specific_below_5_percent": int(sum(1 for f in frac.values() if f < 0.05))}
    (args.out / ("kg_nullquery_e5.json" if args.e5 else "kg_nullquery.json")).write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(result["fostamatinib"], "| generic drugs in the real top 50:", result["n_real_top50_drugs_generic_at_least_5_percent"], "of", len(frac))


if __name__ == "__main__":
    main()
