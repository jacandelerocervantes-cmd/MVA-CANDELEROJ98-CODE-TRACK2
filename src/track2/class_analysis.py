"""Class-level analysis (addendum 11, E10): are drugs of one mechanism class ranked higher than a degree-matched random class?

    python -m track2.class_analysis --out out [--permutations 10000]

Unit = class of approved drugs: (1) ChEMBL mechanism-of-action text of the drug, (2) ATC level 4 (first five characters). A drug can belong to several classes of definition (1).
For every class with >= 3 drugs scored by a line, the statistic is the mean goodness percentile of its members (1 = best of all drugs scored by the line). The null replaces each
member by a random drug of the same degree decile (without replacement), 10,000 times, so a class made of well-connected hubs is compared with other hubs. One-sided p, then
Benjamini-Hochberg FDR across the classes of that definition within the line; passes at FDR < 0.10.
Lines: A (PageRank rank of the patient query without the parental term, PrimeKG degree), B2 (BUB1B-only module, links with >= 2 records, mean of -z over the three STRING thresholds, STRING
mean target degree), B1 (reported only; it did not validate) and C (not evaluated: no usable hits). A class is *supported* if it passes in >= 2 usable lines. Evidence level: 0 = passes in
one line or graph score only; 1 = passes in two independent usable lines. Levels 2 and 3 need a measured effect or a verified sentence and are never assigned here.
Writes `out/class_analysis.json`.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import rankdata

SEED, MIN_CLASS, FDR, N_BINS = 20261008, 3, 0.10, 10
SKIP_MOA = {"unknown", "diagnostic agent", "radioactive agent", "imaging agent"}


def classes_of(universe: dict[str, Any]) -> dict[str, dict[str, set[str]]]:
    moa: dict[str, set[str]] = defaultdict(set)
    atc4: dict[str, set[str]] = defaultdict(set)
    for d in universe["drugs"]:
        for m in d.get("mechanisms") or []:
            text = (m.get("mechanism") or "").strip()
            if text and text.lower() not in SKIP_MOA:
                moa[text].add(d["drug_id"])
        for a in d.get("atc") or []:
            if len(a) >= 5:
                atc4[a[:5]].add(d["drug_id"])
    return {"moa": dict(moa), "atc4": dict(atc4)}


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    adj = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


def evaluate_classes(goodness: dict[str, float], degree: dict[str, float], classes: dict[str, set[str]], permutations: int, seed: int = SEED) -> dict[str, Any]:
    """goodness: drug -> higher is better; degree: drug -> degree. Returns per-class statistics and FDR."""
    ids = sorted(set(goodness) & set(degree))
    pct = rankdata([goodness[i] for i in ids]) / len(ids)
    pos = {d: k for k, d in enumerate(ids)}
    order = np.argsort([degree[i] for i in ids], kind="stable")
    bin_of = np.empty(len(ids), dtype=int)
    for b, members in enumerate(np.array_split(order, N_BINS)):
        bin_of[members] = b
    pools = [np.flatnonzero(bin_of == b) for b in range(N_BINS)]
    rng = np.random.default_rng(seed)
    rows = []
    for name, drugs in sorted(classes.items()):
        members = np.array(sorted(pos[d] for d in drugs if d in pos), dtype=int)
        if len(members) < MIN_CLASS:
            continue
        observed = float(pct[members].mean())
        null_sum = np.zeros(permutations)
        for b in range(N_BINS):
            k = int((bin_of[members] == b).sum())
            if k:
                pool = pools[b]
                pick = rng.random((permutations, len(pool))).argpartition(k - 1, axis=1)[:, :k] if k < len(pool) else np.tile(np.arange(len(pool)), (permutations, 1))
                null_sum += pct[pool][pick].sum(axis=1)
        null = null_sum / len(members)
        rows.append({"class": name, "n_members": int(len(members)), "mean_percentile": observed, "null_mean": float(null.mean()), "p": float((1 + (null >= observed).sum()) / (1 + permutations)),
                     "members": sorted(ids[m] for m in members)})
    if rows:
        adj = bh(np.array([r["p"] for r in rows]))
        for r, q in zip(rows, adj):
            r["fdr"] = float(q)
            r["passes"] = bool(q < FDR)
    return {"n_drugs_scored": len(ids), "n_classes_tested": len(rows), "n_passing": sum(1 for r in rows if r["passes"]), "classes": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--sources", type=Path, default=Path("data/sources/primekg"))
    args = parser.parse_args()
    from track2.benchmark import drug_mean_degree
    from track2.kg_query import universe_match
    from track2.knowledge_graph import load_graph
    from track2.string_net import load_interactome

    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    names = {d["drug_id"]: d["name"] for d in universe["drugs"]}
    classes = classes_of(universe)
    # Line A: goodness = -(best rank); degree = PrimeKG degree (largest over the nodes that map to the drug)
    rank = json.loads((args.out / "kg_query_e5.json").read_text(encoding="utf-8"))["rank_combined"]
    best: dict[str, int] = {}
    for d, r in rank.items():
        best[d] = min(r, best.get(d, r))
    graph = load_graph(args.sources / "nodes.tab", args.sources / "edges.csv")
    deg_a: dict[str, float] = {}
    for node, m in universe_match(graph, universe).items():
        deg_a[m["drug_id"]] = max(deg_a.get(m["drug_id"], 0), float(graph.adj[node].nnz))
    # Lines B: goodness = mean of -z over the three STRING thresholds; degree = mean STRING degree of the drug's targets
    net = load_interactome(700)
    deg_b = drug_mean_degree(universe, net)

    def mean_minus_z(table: dict[str, dict[str, float]]) -> dict[str, float]:
        drugs = set.intersection(*(set(z) for z in table.values()))
        return {d: -float(np.mean([table[t][d] for t in table])) for d in drugs}

    b2 = mean_minus_z(json.loads((args.out / "core_n2" / "core_threshold_scores_sub_BUB1B.json").read_text(encoding="utf-8")))
    b1 = mean_minus_z({t: z for t, z in json.loads((args.out / "phenotype_threshold_scores_e5.json").read_text(encoding="utf-8"))["min2_terms"].items()})
    lines = {"A": ({d: -r for d, r in best.items()}, deg_a, True), "B2": (b2, deg_b, True), "B1": (b1, deg_b, False)}
    result: dict[str, Any] = {"permutations": args.permutations, "fdr": FDR, "min_class": MIN_CLASS, "lines": {}, "usable_lines": ["A", "B2"], "C": "not evaluated: no usable hits"}
    for definition, cls in classes.items():
        for line, (goodness, degree, _usable) in lines.items():
            result["lines"][f"{definition}/{line}"] = evaluate_classes(goodness, degree, cls, args.permutations)
            r = result["lines"][f"{definition}/{line}"]
            print(definition, line, "drugs", r["n_drugs_scored"], "classes", r["n_classes_tested"], "passing", r["n_passing"], flush=True)
    supported = []
    for definition in classes:
        passing = {line: {r["class"] for r in result["lines"][f"{definition}/{line}"]["classes"] if r["passes"]} for line in lines}
        usable = [passing[ln] for ln in ("A", "B2")]
        for c in sorted(usable[0] & usable[1]):
            members = next(r["members"] for r in result["lines"][f"{definition}/A"]["classes"] if r["class"] == c)
            supported.append({"definition": definition, "class": c, "evidence_level": 1, "lines": ["A", "B2"], "members": [names[m] for m in members]})
        result[f"single_line_{definition}"] = {ln: sorted(passing[ln] - (usable[0] & usable[1] if ln in ("A", "B2") else set())) for ln in lines}
    result["supported_in_two_usable_lines"] = supported
    (args.out / "class_analysis.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("classes supported in A and B2:", len(supported), [s["class"] for s in supported][:15])


if __name__ == "__main__":
    main()
