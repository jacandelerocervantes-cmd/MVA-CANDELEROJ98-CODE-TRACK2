"""Empirical robustness of the drug-selection rules (R1-R4, addendum 11 / E16) to errors in the AI-assigned evidence labels.

    PYTHONPATH=src python -m track2.empirical_label_robustness [--out out] [--universe data/universe.json] [--reps 1000] [--seed 20261014]

The labels (evidence_level 0-3, direction) of the verified literature sentences were assigned by AI search tasks and were not reviewed by a specialist.
This module asks how stable the final selection is when a fraction of those labels is wrong.

Method (deterministic: fixed seed, one independent random stream per error rate):
  * Entries: all verified entries of `literature_summary.collect` (drug names already mapped to the names of the universe).
  * For each error rate p in {5, 10, 20, 30} percent and each of `reps` (>= 1000) replicates, exactly round(p * N) distinct entries are
    corrupted (N = number of verified entries, drawn without replacement). A corrupted entry receives one of three error types with equal
    probability: (a) evidence_level +-1 (random sign, kept within 0..3: at 0 only +1, at 3 only -1); (b) direction flip may_help <-> may_harm;
    (c) both. Entries whose direction is not may_help/may_harm (unclear, no_effect, selectively_kills_aneuploid_cells) cannot be flipped and
    always receive a level change instead. The `use` flag, section, and consequence tag are left unchanged (the flag `use=caution` therefore
    keeps producing cautions, as in the original rules).
  * The selection is recomputed with the rules of `literature_summary.build` (classify() is imported and called; the grouping and caution
    logic of build() is replicated in `decide` because build() also needs signals/label files that do not depend on the labels).
  * Per drug: fraction of replicates in which it stays proposed in the same group as in the baseline. Expected number of the baseline
    proposed drugs that survive, with a 95 percent percentile interval over the replicates. Caution set: per baseline-caution drug, the
    fraction of replicates in which it stays in the cautions list, plus the mean number of new cautions.
Writes out/empirical/label_robustness.json.
"""

from __future__ import annotations

import argparse
import copy
import json
import random
from pathlib import Path
from typing import Any

from track2 import literature_summary as LS

RATES = (0.05, 0.10, 0.20, 0.30)
FLIPPABLE = {"may_help": "may_harm", "may_harm": "may_help"}


def decide(data: dict[str, dict[str, list[dict[str, Any]]]], approved: set[str]) -> tuple[dict[str, str], set[str]]:
    """Return ({drug: group} for proposed drugs, set of caution drugs); same logic as literature_summary.build."""
    proposed: dict[str, str] = {}
    caution: set[str] = set()
    for name in sorted(n for n in data if n in approved):
        c = LS.classify(data[name]["entries"])
        group = "constitutional" if c["constitutional"] else ("tumour_context" if c["tumour"] else None)
        if name in LS.EDITORIAL:
            group = None
        if c["caution"] and name not in ("BORTEZOMIB", "CHLOROQUINE"):
            group = None
        if name in LS.SIBLINGS:
            group = "class_sibling"
        if group:
            proposed[name] = group
        if c["caution"] or (c["harms"] and name not in proposed):
            caution.add(name)
    return proposed, caution


def corrupt(data: dict[str, dict[str, list[dict[str, Any]]]], rate: float, rng: random.Random) -> dict[str, dict[str, list[dict[str, Any]]]]:
    """Return a copy of `data` with round(rate * N) entries corrupted as documented in the module docstring."""
    refs = [(n, i) for n in sorted(data) for i in range(len(data[n]["entries"]))]
    k = round(rate * len(refs))
    new = {n: {"entries": [dict(r) for r in v["entries"]]} for n, v in data.items()}
    for n, i in rng.sample(refs, k):
        r = new[n]["entries"][i]
        kind = rng.choice(("level", "flip", "both"))
        if r["direction"] not in FLIPPABLE:
            kind = "level"
        if kind in ("flip", "both"):
            r["direction"] = FLIPPABLE[r["direction"]]
        if kind in ("level", "both"):
            steps = [s for s in (-1, 1) if 0 <= r["level"] + s <= 3]
            r["level"] += rng.choice(steps)
    return new


def interval(values: list[int]) -> list[float]:
    s = sorted(values)
    n = len(s)
    return [float(s[int(0.025 * (n - 1))]), float(s[int(round(0.975 * (n - 1)))])]


def simulate(data: dict[str, dict[str, list[dict[str, Any]]]], approved: set[str], reps: int = 1000, seed: int = 20261014, rates: tuple[float, ...] = RATES) -> dict[str, Any]:
    base_prop, base_caut = decide(data, approved)
    n_entries = sum(len(v["entries"]) for v in data.values())
    result: dict[str, Any] = {"baseline_proposed": dict(sorted(base_prop.items())), "baseline_caution": sorted(base_caut), "n_entries": n_entries,
                              "replicates": reps, "seed": seed, "rates": {}}
    for rate in rates:
        rng = random.Random(f"{seed}-{rate}")
        stay = {d: 0 for d in base_prop}
        cstay = {d: 0 for d in base_caut}
        survivors: list[int] = []
        new_caut = 0
        for _ in range(reps):
            prop, caut = decide(corrupt(data, rate, rng), approved)
            ok = [d for d, g in base_prop.items() if prop.get(d) == g]
            for d in ok:
                stay[d] += 1
            for d in base_caut & caut:
                cstay[d] += 1
            survivors.append(len(ok))
            new_caut += len(caut - base_caut)
        result["rates"][f"{round(rate * 100)}"] = {
            "n_corrupted_entries": round(rate * n_entries),
            "survival_fraction_by_drug": {d: stay[d] / reps for d in sorted(stay)},
            "expected_survivors": sum(survivors) / reps, "survivors_95_interval": interval(survivors), "n_baseline_proposed": len(base_prop),
            "caution_retained_fraction_by_drug": {d: cstay[d] / reps for d in sorted(cstay)}, "mean_new_cautions": new_caut / reps}
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--reps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261014)
    args = parser.parse_args()
    approved = {d["name"].upper() for d in json.loads(args.universe.read_text(encoding="utf-8"))["drugs"]}
    data = LS.collect(args.out / "literature_hypotheses")
    result = simulate(data, approved, max(args.reps, 1000), args.seed)
    result["method"] = __doc__
    dest = args.out / "empirical" / "label_robustness.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print("baseline proposed:", result["baseline_proposed"])
    for r, v in result["rates"].items():
        print(f"{r}%: expected survivors {v['expected_survivors']:.2f} of {v['n_baseline_proposed']}, 95% interval {v['survivors_95_interval']}")
        print("   ", {d: round(f, 3) for d, f in v["survival_fraction_by_drug"].items()})


if __name__ == "__main__":
    main()
