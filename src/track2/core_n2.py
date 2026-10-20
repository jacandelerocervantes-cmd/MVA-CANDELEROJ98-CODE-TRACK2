"""Sensitivity analysis of line B (protocol/addendum_four_lines.md, "Sensitivity analyses after the independent review", item 1).

    python -m track2.core_n2 --out out/core_n2

Runs `core_sensitivity` again (mitosis modules without generic complexes, nine configurations, same seed) but with the measured drug-target links limited to those with at least
two activity records (`MIN_RECORDS`, the rule already used in the plausibility tiers) and without the similarity-predicted links. Compares the robust drugs with the original run
(`out/core_sensitivity.json`). Writes `<out>/core_sensitivity_n2.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from track2.core_sensitivity import run

MIN_RECORDS = 2


def restrict_links(e1: dict[str, dict[str, dict[str, Any]]], min_records: int = MIN_RECORDS) -> dict[str, dict[str, dict[str, Any]]]:
    """Keep only the measured links backed by at least `min_records` activity records; drugs left without links are dropped."""
    out = {}
    for drug, targets in e1.items():
        kept = {t: v for t, v in targets.items() if (v.get("n") or 0) >= min_records}
        if kept:
            out[drug] = kept
    return out


def compare(original: dict[str, Any], strict: dict[str, Any]) -> dict[str, Any]:
    """Which drugs robust in the original run stay robust, with their counts in both runs."""
    new = {d["drug_id"]: d["robustness_count_core"] for d in strict["drugs"]}
    rows = [{"drug_id": d["drug_id"], "name": d["name"], "count_original": d["robustness_count_core"], "count_n2_no_E2": new.get(d["drug_id"], 0),
             "robust_original": d["L1_core"], "robust_n2_no_E2": new.get(d["drug_id"], 0) >= strict["rule"]["robust_min"]} for d in original["drugs"] if d["L1_core"]]
    return {"originally_robust": rows, "n_robust_original": original["n_robust_core"], "n_robust_n2_no_E2": strict["n_robust_core"],
            "robust_in_strict_run": [d["name"] for d in strict["drugs"] if d["L1_core"]]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out/core_n2"))
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--original", type=Path, default=Path("out/core_sensitivity.json"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    universe = json.loads(Path("data/universe.json").read_text(encoding="utf-8"))
    e1 = restrict_links(json.loads(Path("data/dti_e1_all.json").read_text(encoding="utf-8"))["edges"])
    strict = run(args.out, args.mechanism_dir, universe, e1, {})
    result = {"min_records": MIN_RECORDS, "n_drugs_with_links": len(e1), "strict": strict, "comparison": compare(json.loads(args.original.read_text(encoding="utf-8")), strict)}
    (args.out / "core_sensitivity_n2.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    for r in result["comparison"]["originally_robust"]:
        print(r["name"], "original", r["count_original"], "| n>=2 without E2:", r["count_n2_no_E2"], "robust" if r["robust_n2_no_E2"] else "NOT robust")
    print("robust in the strict run:", result["comparison"]["robust_in_strict_run"])


if __name__ == "__main__":
    main()
