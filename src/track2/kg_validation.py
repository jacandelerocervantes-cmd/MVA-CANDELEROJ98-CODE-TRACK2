"""Validation of line A (rule: `protocol/addendum_four_lines.md`): hidden known indications.

    python -m track2.kg_validation --out out

200 diseases (random sample, seed 20261007) with at least three drugs linked by an `indication` edge. Because the graph used for ranking has no drug-disease edge at all, every
indication is hidden. For each disease the personalised PageRank from the disease node ranks all drug nodes and the ROC AUC of its indicated drugs against the other drug
nodes is computed. The line is used only if the mean AUC >= 0.70 and the lower bound of the 95 % bootstrap interval over diseases > 0.55. Writes `out/kg_validation.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from track2.benchmark_stats import DiseaseData, auc_from_scores, drug_bins, summarize
from track2.knowledge_graph import Graph, load_graph, personalised_pagerank, restart_vector

N_DISEASES, MIN_DRUGS, SEED = 200, 3, 20261007
AUC_MIN, CI_LOWER_MIN = 0.70, 0.55
CHUNK = 40


def sample_diseases(graph: Graph, n: int = N_DISEASES, seed: int = SEED) -> list[int]:
    eligible = sorted(d for d, drugs in graph.indications.items() if len(drugs) >= MIN_DRUGS and graph.node_type[d] == "disease")
    rng = np.random.default_rng(seed)
    return sorted(rng.choice(eligible, size=min(n, len(eligible)), replace=False).tolist())


def stratify(per: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Post hoc, not a decision rule: mean AUC by tercile of the disease node's degree (an ultrarare disease has few edges)."""
    rows = sorted(per.values(), key=lambda v: v["degree"])
    out = {}
    for name, part in zip(("lowest_degree_third", "middle_third", "highest_degree_third"), np.array_split(np.arange(len(rows)), 3)):
        if len(part):
            sel = [rows[i] for i in part]
            out[name] = {"n": len(sel), "degree_range": [sel[0]["degree"], sel[-1]["degree"]], "mean_auc": round(float(np.mean([v["auc"] for v in sel])), 4)}
    return out


def validate(graph: Graph, diseases: list[int]) -> dict[str, Any]:
    drug_nodes = np.flatnonzero(graph.node_type == "drug")
    degree = np.asarray(graph.adj.sum(axis=1)).ravel()[drug_nodes]
    bins = drug_bins(degree)
    items, per = [], {}
    n = graph.adj.shape[0]
    for start in range(0, len(diseases), CHUNK):
        chunk = diseases[start:start + CHUNK]
        scores = personalised_pagerank(graph.adj, np.column_stack([restart_vector(n, {d: 1.0}) for d in chunk]))
        for k, d in enumerate(chunk):
            positive = np.isin(drug_nodes, list(graph.indications[d]))
            if positive.sum() < MIN_DRUGS:
                continue
            s = scores[drug_nodes, k]
            auc = auc_from_scores(s, positive)
            per[str(d)] = {"name": str(graph.node_name[d]), "n_known": int(positive.sum()), "auc": round(auc, 4), "degree": int(graph.adj[d].nnz)}
            items.append(DiseaseData(s, positive, bins))
    summary = summarize([v["auc"] for v in per.values()], items)
    validated = bool(summary["mean_auc"] >= AUC_MIN and summary["ci95_low"] > CI_LOWER_MIN)
    by_degree = stratify(per)
    return {"rule": {"n_diseases": N_DISEASES, "min_drugs": MIN_DRUGS, "auc_min": AUC_MIN, "ci_lower_min": CI_LOWER_MIN, "seed": SEED},
            "n_evaluated": len(per), "summary": summary, "validated": validated,
            "post_hoc_by_disease_degree": by_degree, "per_disease": per}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--sources", type=Path, default=Path("data/sources/primekg"))
    args = parser.parse_args()
    graph = load_graph(args.sources / "nodes.tab", args.sources / "edges.csv")
    print("graph:", graph.adj.shape[0], "nodes,", graph.adj.nnz // 2, "undirected edges without drug-disease links;", len(graph.indications), "diseases with indications", flush=True)
    result = validate(graph, sample_diseases(graph))
    (args.out / "kg_validation.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    s = result["summary"]
    print("post hoc by disease degree:", result["post_hoc_by_disease_degree"])
    print("evaluated", result["n_evaluated"], "| mean AUC", round(s["mean_auc"], 3), "| CI", round(s["ci95_low"], 3), round(s["ci95_high"], 3),
          "| permutation p", s["permutation"]["p_one_sided"], "| line A validated:", result["validated"])


if __name__ == "__main__":
    main()
