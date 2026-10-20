"""Line A, what is actually used (protocol/addendum_four_lines.md, "Sensitivity analyses after the independent review", item 3).

    python -m track2.kg_extra --out out

On the same 200 diseases of `kg_validation`: (a) recall of the indicated drugs in the top 50 of the personalised PageRank (all drug nodes) and in the top 50 of the drugs matched to
the register universe; (b) a degree-only baseline (score = degree of the drug) and the per-disease excess of the AUC over it, with a bootstrap interval over diseases; (c) the same
with the restart on the disease node plus its own phenotype nodes (`disease_phenotype_positive`). Reported as limitations; the pass rule of line A is not changed. Writes `out/kg_extra.json`.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from track2.benchmark_stats import auc_from_scores, bootstrap_mean_ci
from track2.kg_query import universe_match
from track2.kg_validation import CHUNK, MIN_DRUGS, sample_diseases
from track2.knowledge_graph import Graph, load_graph, personalised_pagerank, restart_vector

TOP = 50


def phenotype_neighbours(edges_csv: Path, graph: Graph, diseases: set[int]) -> dict[int, list[int]]:
    """Disease node -> its phenotype nodes (relation `disease_phenotype_positive`), read from the raw edge file."""
    out: dict[int, list[int]] = {d: [] for d in diseases}
    with open(edges_csv, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        h = next(reader)
        ir, ix, iy = h.index("relation"), h.index("x_index"), h.index("y_index")
        for row in reader:
            if row[ir] == "disease_phenotype_positive":
                x, y = int(row[ix]), int(row[iy])
                if x in out and graph.node_type[y] == "effect/phenotype":
                    out[x].append(y)
                elif y in out and graph.node_type[x] == "effect/phenotype":
                    out[y].append(x)
    return out


def recall_at(scores: np.ndarray, positive: np.ndarray, top: int = TOP) -> float:
    order = np.argsort(-scores, kind="stable")[:top]
    return float(positive[order].sum() / positive.sum())


def evaluate(graph: Graph, diseases: list[int], phenos: dict[int, list[int]], matched: set[int]) -> dict[str, Any]:
    drug_nodes = np.flatnonzero(graph.node_type == "drug")
    in_universe = np.array([int(n) in matched for n in drug_nodes])
    degree = np.asarray(graph.adj.sum(axis=1)).ravel()[drug_nodes].astype(float)
    n = graph.adj.shape[0]
    rows = []
    for variant in ("disease_only", "disease_plus_own_phenotypes"):
        for start in range(0, len(diseases), CHUNK):
            chunk = diseases[start:start + CHUNK]
            seeds = [{d: 1.0, **({p: 1.0 for p in phenos.get(d, [])} if variant != "disease_only" else {})} for d in chunk]
            scores = personalised_pagerank(graph.adj, np.column_stack([restart_vector(n, s) for s in seeds]))
            for k, d in enumerate(chunk):
                positive = np.isin(drug_nodes, list(graph.indications[d]))
                if positive.sum() < MIN_DRUGS:
                    continue
                s = scores[drug_nodes, k]
                row = {"variant": variant, "disease": int(d), "auc": auc_from_scores(s, positive), "auc_degree_baseline": auc_from_scores(degree, positive),
                       "recall50_all": recall_at(s, positive), "recall50_degree_baseline": recall_at(degree, positive)}
                if (positive & in_universe).sum() >= 1:
                    row["recall50_universe"] = recall_at(np.where(in_universe, s, -1.0), positive & in_universe)
                    row["recall50_universe_degree_baseline"] = recall_at(np.where(in_universe, degree, -1.0), positive & in_universe)
                rows.append(row)
    out: dict[str, Any] = {}
    for variant in ("disease_only", "disease_plus_own_phenotypes"):
        r = [x for x in rows if x["variant"] == variant]
        excess = np.array([x["auc"] - x["auc_degree_baseline"] for x in r])
        lo, hi = bootstrap_mean_ci(excess)
        u = [x for x in r if "recall50_universe" in x]
        out[variant] = {"n_diseases": len(r), "mean_auc": float(np.mean([x["auc"] for x in r])), "mean_auc_degree_baseline": float(np.mean([x["auc_degree_baseline"] for x in r])),
                        "mean_excess_over_degree": float(excess.mean()), "excess_ci95": [lo, hi], "diseases_where_pagerank_beats_degree": int((excess > 0).sum()),
                        "mean_recall50_all": float(np.mean([x["recall50_all"] for x in r])), "mean_recall50_degree_baseline": float(np.mean([x["recall50_degree_baseline"] for x in r])),
                        "n_with_universe_drug": len(u), "mean_recall50_universe": float(np.mean([x["recall50_universe"] for x in u])),
                        "mean_recall50_universe_degree_baseline": float(np.mean([x["recall50_universe_degree_baseline"] for x in u])),
                        "mean_n_phenotype_nodes": float(np.mean([len(phenos.get(x["disease"], [])) for x in r]))}
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--sources", type=Path, default=Path("data/sources/primekg"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    args = parser.parse_args()
    graph = load_graph(args.sources / "nodes.tab", args.sources / "edges.csv")
    diseases = sample_diseases(graph)
    phenos = phenotype_neighbours(args.sources / "edges.csv", graph, set(diseases))
    matched = set(universe_match(graph, json.loads(args.universe.read_text(encoding="utf-8"))))
    result = {"top": TOP, "n_diseases_sampled": len(diseases), "results": evaluate(graph, diseases, phenos, matched)}
    (args.out / "kg_extra.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    for k, v in result["results"].items():
        print(k, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in v.items() if a != "excess_ci95"}, "excess CI", [round(x, 3) for x in v["excess_ci95"]])


if __name__ == "__main__":
    main()
