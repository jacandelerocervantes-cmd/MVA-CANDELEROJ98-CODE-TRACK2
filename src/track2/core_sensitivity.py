"""Sensitivity analysis without generic complexes (rule: `protocol/addendum_core.md`).

    python -m track2.core_sensitivity --out out

Rebuilds each submission-derived module without proteasome subunits and tubulins, recomputes the proximity of all approved drugs (expanded links, STRING 400/700/900) and
reports for every drug the nine-configuration robustness count with the core modules (`L1_core`). Writes `out/core_threshold_scores_<module>.json` and
`out/core_sensitivity.json`.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

from track2.augmented_proximity import augmented_universe
from track2.consensus import SCHEMES, THRESHOLDS, ROBUST_MIN, combined, top, weights
from track2.modules import load_mechanism
from track2.run_pipeline import STRING_SCORES, run_mechanism

EXCLUDE = r"^(PSM[ABCDG][0-9]+|ADRM1|POMP|TUBA\w+|TUBB\w*|TUBG[12])$"


def core_mechanism(mech: dict[str, Any]) -> dict[str, Any]:
    out = json.loads(json.dumps(mech))
    out["id"] = "core_" + mech["id"]
    out["seeds"]["exclude_symbol_regex"] = EXCLUDE
    return out


def run(out_dir: Path, mech_dir: Path, universe: dict[str, Any], e1: dict[str, Any], e2: dict[str, Any], iterations: int = 1000) -> dict[str, Any]:
    aug = augmented_universe(universe, e1, e2)
    mechs = [m for m in (load_mechanism(p) for p in sorted(mech_dir.glob("*.json"))) if m.get("status") == "mechanism_candidate"]
    cores = [core_mechanism(m) for m in mechs]
    by = {}
    for m in cores:
        res = {}
        for score in STRING_SCORES:
            rows = run_mechanism(m, aug, score, iterations)["rows"]
            res[str(score)] = {r["drug_id"]: round(r["proximity"]["z"], 6) for r in rows}
        (out_dir / f"core_threshold_scores_{m['id'].removeprefix('core_')}.json").write_text(json.dumps(res, sort_keys=True), encoding="utf-8", newline="\n")
        by[m["id"]] = res
    wts = weights(mechs)  # the same weights, keyed by the original module ids
    key = {c["id"]: c["id"].removeprefix("core_") for c in cores}
    table = {(t, s): combined({key[c]: by[c][t] for c in by}, wts[s]) for t, s in itertools.product(THRESHOLDS, SCHEMES)}
    counts: dict[str, int] = {}
    for sc in table.values():
        for d in top(sc):
            counts[d] = counts.get(d, 0) + 1
    names = {d["drug_id"]: d["name"] for d in universe["drugs"]}
    drugs = [{"drug_id": d, "name": names.get(d), "robustness_count_core": n, "L1_core": n >= ROBUST_MIN} for d, n in sorted(counts.items(), key=lambda kv: (-kv[1], names.get(kv[0]) or ""))]
    return {"rule": {"exclude_regex": EXCLUDE, "configurations": len(table), "robust_min": ROBUST_MIN}, "n_robust_core": sum(1 for x in drugs if x["L1_core"]), "drugs": drugs}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    args = parser.parse_args()
    universe = json.loads(Path("data/universe.json").read_text(encoding="utf-8"))
    e1 = json.loads(Path("data/dti_e1_all.json").read_text(encoding="utf-8"))["edges"]
    e2 = json.loads((args.out / "dti_e2_edges.json").read_text(encoding="utf-8")).get("edges") or {}
    result = run(args.out, args.mechanism_dir, universe, e1, e2)
    (args.out / "core_sensitivity.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(result["n_robust_core"], "drugs robust with the core modules;", [d["name"] for d in result["drugs"] if d["L1_core"]][:25])


if __name__ == "__main__":
    main()
