"""Validation tests P1, P3 and P4 on a submission-derived module (rules: `protocol/addendum_robustness.md`, analysis A).

    python -m track2.validate_submission_modules --mechanism config/mechanisms_from_submission/sub_BUB1B.json --out out

Same code and decision rules as `track2.validation` (frozen protocol); only the module changes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from track2.modules import build_module, load_mechanism
from track2.run_pipeline import load_universe
from track2.string_net import load_interactome
from track2.validation import p1_leave_one_out, p3_robustness

P1_MAX_SEEDS = 120


def validate(path: Path, out: Path, universe: dict) -> None:
    mechanism = load_mechanism(path)
    interactome = load_interactome(700)
    module = build_module(mechanism, interactome)
    results = {"protocol": "protocol/frozen_protocol.md + protocol/addendum_robustness.md (analysis A)", "mechanism": mechanism["id"],
               "P1_leave_one_out": p1_leave_one_out(interactome, module["nodes"], P1_MAX_SEEDS),
               "P3_P4_robustness": p3_robustness(mechanism, universe)}
    (out / f"validation_{mechanism['id']}.json").write_text(json.dumps(results, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(mechanism["id"], results["P1_leave_one_out"]["auc"], results["P3_P4_robustness"]["verdict"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanism", type=Path, help="one module file")
    parser.add_argument("--all", action="store_true", help="every module file of --mechanism-dir")
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    if not args.all and not args.mechanism:
        parser.error("give --mechanism FILE or --all")
    paths = sorted(args.mechanism_dir.glob("*.json")) if args.all else [args.mechanism]
    universe = load_universe()
    for path in paths:
        validate(path, args.out, universe)


if __name__ == "__main__":
    main()
