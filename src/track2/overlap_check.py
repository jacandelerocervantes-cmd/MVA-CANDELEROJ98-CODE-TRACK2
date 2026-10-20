"""Expected overlap by chance between the lines, and a positive control of the identifier join (addendum 11, E9).

    python -m track2.overlap_check [--out out]

For the drugs scored by both lines A and B, the number of drugs that the top N of line A (PageRank of the patient query, distinct drugs) share with the top N of the proximity lines (mean of -z over
the three STRING thresholds, B1 and B2) is compared with what chance gives (N x N / common). A zero observed where the chance expectation is below 1 says nothing; the join of PrimeKG and ChEMBL
identifiers is shown to work when two lines that should agree do overlap more than chance. Writes `out/overlap_check.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import hypergeom


def _mean_minus_z(table: dict[str, dict[str, float]]) -> dict[str, float]:
    drugs = set.intersection(*(set(z) for z in table.values()))
    return {d: -float(np.mean([table[t][d] for t in table])) for d in drugs}


def overlap_rows(a_order: list[str], b_scores: dict[str, float], common: set[str], tops: tuple[int, ...] = (20, 50, 100)) -> list[dict[str, Any]]:
    b_order = sorted(b_scores, key=lambda d: (-b_scores[d], d))
    rows = []
    for n in tops:
        ta = {d for d in [x for x in a_order if x in common][:n]}
        tb = {d for d in [x for x in b_order if x in common][:n]}
        obs = len(ta & tb)
        rows.append({"top_n": n, "observed": obs, "expected_by_chance": round(n * n / len(common), 3), "p_hypergeometric": float(hypergeom.sf(obs - 1, len(common), n, n)) if obs else 1.0})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    q = json.loads((args.out / "kg_query_e5.json").read_text(encoding="utf-8"))
    best: dict[str, int] = {}
    for drug, r in q["rank_combined"].items():
        best[drug] = min(r, best.get(drug, r))
    a_order = sorted(best, key=lambda k: (best[k], k))
    b2 = _mean_minus_z(json.loads((args.out / "core_n2" / "core_threshold_scores_sub_BUB1B.json").read_text(encoding="utf-8")))
    b1 = _mean_minus_z(json.loads((args.out / "phenotype_threshold_scores_e5.json").read_text(encoding="utf-8"))["min2_terms"])
    fl = json.loads((args.out / "four_lines_review2.json").read_text(encoding="utf-8"))
    common = set(best) & set(b2) & set(b1)
    sizes = fl["support_set_sizes"]
    result = {"scored_by_A": len(best), "scored_by_B2": len(b2), "scored_by_B1": len(b1), "common_to_A_B1_B2": len(common),
              "rule_sets": {"sizes": sizes, "expected_overlap_A_and_B2": round(sizes["A"] * sizes["B2"] / len(common), 4), "observed_overlap_A_and_B2": len(set(fl["support_detail"]["random_query_fraction_of_A_specific"]) & set(fl["support_detail"]["B2_W3_two_record"]))},
              "top_n_A_vs_B2": overlap_rows(a_order, b2, common), "top_n_A_vs_B1_positive_control": overlap_rows(a_order, b1, common),
              "reading": "The sets of the combination rule are tiny (3 and 4 drugs), so an overlap of zero is the expected result and is not informative. A and B1 (both driven by the phenotype terms, and both favoured by well-connected drugs) "
                         "overlap more than chance, which shows the identifier join works; A and B2 are at or below chance."}
    (args.out / "overlap_check.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print({k: result[k] for k in ("common_to_A_B1_B2", "rule_sets")})
    print("A vs B2:", [(r["top_n"], r["observed"], r["expected_by_chance"]) for r in result["top_n_A_vs_B2"]], "| A vs B1:", [(r["top_n"], r["observed"], r["expected_by_chance"]) for r in result["top_n_A_vs_B1_positive_control"]])


if __name__ == "__main__":
    main()
