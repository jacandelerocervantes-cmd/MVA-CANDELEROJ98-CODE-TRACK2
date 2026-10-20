"""Proximity of every approved drug to each submission-derived module with expanded drug-target links (`protocol/addendum_dti.md`).

    python -m track2.augmented_proximity --variant e0_e1 --out out
    python -m track2.augmented_proximity --variant e0_e1_e2 --out out

Variants: `e0_e1` = curated + measured links; `e0_e1_e2` = also similarity-predicted links (only if the retrospective validation passed).
Writes `out/dti_threshold_scores_<variant>_<module>.json` = {STRING score: {drug_id: z}} for STRING 700, 400 and 900. Same null model,
iterations and seed as the main analysis; only the target sets of the drugs change.
"""

from __future__ import annotations

import argparse
import copy
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from track2.modules import load_mechanism
from track2.run_pipeline import STRING_SCORES, run_mechanism

VARIANTS = ("e0_e1", "e0_e1_e2")


def augmented_universe(universe: dict[str, Any], e1: dict[str, dict[str, Any]], e2: dict[str, dict[str, float]] | None) -> dict[str, Any]:
    out = copy.deepcopy(universe)
    for drug in out["drugs"]:
        extra = set(e1.get(drug["drug_id"], {}))
        if e2:
            extra |= set(e2.get(drug["drug_id"], {}))
        drug["human_targets"] = sorted(set(drug["human_targets"]) | extra)
    return out


def _run(args: tuple[str, str, str, dict[str, Any], int]) -> tuple[str, dict[str, dict[str, float]]]:
    variant, mech_path, out_dir, universe, iterations = args
    mech = load_mechanism(Path(mech_path))
    result = {}
    for score in STRING_SCORES:
        rows = run_mechanism(mech, universe, score, iterations)["rows"]
        result[str(score)] = {r["drug_id"]: round(r["proximity"]["z"], 6) for r in rows}
    Path(out_dir, f"dti_threshold_scores_{variant}_{mech['id']}.json").write_text(json.dumps(result, sort_keys=True), encoding="utf-8", newline="\n")
    return mech["id"], {k: len(v) for k, v in result.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=VARIANTS, required=True)
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--e1", type=Path, default=Path("data/dti_e1_all.json"))
    parser.add_argument("--e2", type=Path, default=Path("out/dti_e2_edges.json"))
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    e1 = json.loads(args.e1.read_text(encoding="utf-8"))["edges"]
    e2 = None
    if args.variant == "e0_e1_e2":
        e2_file = json.loads(args.e2.read_text(encoding="utf-8"))
        if e2_file["tau"] is None:
            raise SystemExit("E2 did not pass its retrospective validation (no tau): variant e0_e1_e2 is not run.")
        e2 = e2_file["edges"]
    aug = augmented_universe(universe, e1, e2)
    paths = [p for p in sorted(args.mechanism_dir.glob("*.json")) if load_mechanism(p).get("status") == "mechanism_candidate"]
    jobs = [(args.variant, str(p), str(args.out), aug, args.iterations) for p in paths]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for mid, counts in ex.map(_run, jobs):
            print(mid, counts, flush=True)


if __name__ == "__main__":
    main()
