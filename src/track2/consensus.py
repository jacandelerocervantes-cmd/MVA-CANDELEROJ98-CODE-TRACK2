"""Consensus over STRING thresholds and weighting schemes (rules: `protocol/addendum_robustness.md`, analyses B and C).

    python -m track2.consensus --out out

Reads `out/threshold_scores_<id>.json` (z of every scorable drug at STRING 400/700/900, one file per gene module) and the
module configurations. Nine configurations = 3 thresholds x 3 weighting schemes (W1 epcr weights, W2 equal weights, W3 only
the genes with a ClinVar pathogenic variant). A drug is robust if it is in the top 20 of the combined ranking in at least 7 of
the 9 configurations, fragile if in 3 or fewer. Nothing here changes the scores: it only reports how stable they are.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path
from typing import Any

from scipy.stats import spearmanr

THRESHOLDS = ("700", "400", "900")
SCHEMES = ("W1_epcr", "W2_equal", "W3_pathogenic_gene_only")
TOP_N = 20
ROBUST_MIN, FRAGILE_MAX = 7, 3


def _is_pathogenic(mechanism: dict[str, Any]) -> bool:
    return any("pathogenic" in (c.get("classification") or "").lower() for c in mechanism.get("evidence_from_submission", {}).get("clinvar", []))


def weights(mechanisms: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    out = {"W1_epcr": {m["id"]: float(m["weight"]) for m in mechanisms}, "W2_equal": {m["id"]: 1.0 for m in mechanisms}}
    path = {m["id"]: 1.0 for m in mechanisms if _is_pathogenic(m)}
    out["W3_pathogenic_gene_only"] = path or out["W1_epcr"]
    return out


def combined(z_by_module: dict[str, dict[str, float]], w: dict[str, float]) -> dict[str, float]:
    scores: dict[str, float] = {}
    for module, weight in w.items():
        for drug, z in z_by_module[module].items():
            scores[drug] = scores.get(drug, 0.0) + weight * max(0.0, -z)
    return scores


def top(scores: dict[str, float], n: int = TOP_N) -> list[str]:
    ranked = sorted((d for d, v in scores.items() if v > 0), key=lambda d: (-scores[d], d))
    return ranked[:n]


def _jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if sa | sb else 0.0


def run(out: Path, mechanism_dir: Path, universe: dict[str, Any]) -> dict[str, Any]:
    mechanisms = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(mechanism_dir.glob("*.json"))]
    mechanisms = [m for m in mechanisms if m.get("status") == "mechanism_candidate" and (out / f"threshold_scores_{m['id']}.json").exists()]
    scores_by_threshold = {m["id"]: json.loads((out / f"threshold_scores_{m['id']}.json").read_text(encoding="utf-8")) for m in mechanisms}
    names = {d["drug_id"]: d["name"] for d in universe["drugs"]}
    wts = weights(mechanisms)
    table: dict[tuple[str, str], dict[str, float]] = {}
    for thr, scheme in itertools.product(THRESHOLDS, SCHEMES):
        z = {mid: scores_by_threshold[mid][thr] for mid in scores_by_threshold}
        table[(thr, scheme)] = combined(z, wts[scheme])
    tops = {k: top(v) for k, v in table.items()}
    counts: dict[str, int] = {}
    for lst in tops.values():
        for d in lst:
            counts[d] = counts.get(d, 0) + 1
    drugs = []
    for d, n in sorted(counts.items(), key=lambda kv: (-kv[1], names.get(kv[0]) or "")):
        cls = "robust" if n >= ROBUST_MIN else "fragile" if n <= FRAGILE_MAX else "intermediate"
        drugs.append({"drug_id": d, "name": names.get(d), "robustness_count": n, "of": len(tops), "class": cls,
                      "combined_700_W1": round(table[("700", "W1_epcr")].get(d, 0.0), 4),
                      "in_top20_by_config": {f"{t}|{s}": d in tops[(t, s)] for (t, s) in tops}})
    agreement: dict[str, Any] = {"schemes_at_700": {}, "thresholds_under_W1": {}}
    for a, b in itertools.combinations(SCHEMES, 2):
        common = sorted(set(table[("700", a)]) & set(table[("700", b)]))
        agreement["schemes_at_700"][f"{a} vs {b}"] = {
            "jaccard_top20": _jaccard(tops[("700", a)], tops[("700", b)]),
            "spearman": float(spearmanr([table[("700", a)][d] for d in common], [table[("700", b)][d] for d in common]).statistic)}
    for t in ("400", "900"):
        common = sorted(set(table[("700", "W1_epcr")]) & set(table[(t, "W1_epcr")]))
        agreement["thresholds_under_W1"][f"700 vs {t}"] = {
            "jaccard_top20": _jaccard(tops[("700", "W1_epcr")], tops[(t, "W1_epcr")]),
            "spearman": float(spearmanr([table[("700", "W1_epcr")][d] for d in common], [table[(t, "W1_epcr")][d] for d in common]).statistic)}
    summary = {c: sum(1 for x in drugs if x["class"] == c) for c in ("robust", "intermediate", "fragile")}
    return {"rule": {"configurations": len(tops), "top_n": TOP_N, "robust_min": ROBUST_MIN, "fragile_max": FRAGILE_MAX},
            "modules": [m["id"] for m in mechanisms], "weights": wts, "summary": summary, "agreement": agreement, "drugs": drugs}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    args = parser.parse_args()
    result = run(args.out, args.mechanism_dir, json.loads(args.universe.read_text(encoding="utf-8")))
    (args.out / "consensus_ranking.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(args.out / "consensus_ranking.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["drug_id", "name", "robustness_count", "of", "class", "combined_700_W1"])
        for d in result["drugs"]:
            writer.writerow([d["drug_id"], d["name"], d["robustness_count"], d["of"], d["class"], d["combined_700_W1"]])
    print(result["summary"], json.dumps(result["agreement"], indent=1))


if __name__ == "__main__":
    main()
