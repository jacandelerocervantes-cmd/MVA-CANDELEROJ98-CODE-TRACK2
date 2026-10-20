"""How much does the mitosis line depend on the weight given to the secondary Track 1 genes? (protocol/addendum_four_lines.md, "Three further analyses", item 1.)

    python -m track2.track1_scenarios --out out

For each weighting scheme separately (W1 weights from Track 1, W2 equal weights, W3 the pathogenic gene only = the single-hit-like scenario) a drug is *scheme-robust* if it is among the
20 best at all three STRING thresholds. Reports the scheme-robust set of each scheme, for the original links and for the two-record variant (`out/core_n2/`), and the drugs that are
scheme-robust in all three schemes. Writes `out/track1_scenarios.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from track2.consensus import SCHEMES, THRESHOLDS, combined, top, weights


def scheme_robust(scores: dict[str, dict[str, dict[str, float]]], wts: dict[str, dict[str, float]]) -> dict[str, set[str]]:
    """scores[module][threshold][drug] = z. Returns scheme -> drugs in the top 20 at every threshold."""
    out = {}
    for scheme in SCHEMES:
        lists = [set(top(combined({m: scores[m][t] for m in scores}, wts[scheme]))) for t in THRESHOLDS]
        out[scheme] = set.intersection(*lists)
    return out


def load_scores(folder: Path, prefix: str, module_ids: list[str]) -> dict[str, dict[str, dict[str, float]]]:
    return {m: json.loads((folder / f"{prefix}{m}.json").read_text(encoding="utf-8")) for m in module_ids}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    args = parser.parse_args()
    mechs = [m for m in (json.loads(p.read_text(encoding="utf-8")) for p in sorted(args.mechanism_dir.glob("*.json"))) if m.get("status") == "mechanism_candidate"]
    ids = [m["id"] for m in mechs]
    wts = weights(mechs)
    names = {d["drug_id"]: d["name"] for d in json.loads(args.universe.read_text(encoding="utf-8"))["drugs"]}
    result: dict[str, Any] = {"weights": wts, "variants": {}}
    for label, folder in (("original_links", args.out), ("two_record_links_no_predicted", args.out / "core_n2")):
        robust = scheme_robust(load_scores(folder, "core_threshold_scores_", ids), wts)
        common = set.intersection(*robust.values())
        result["variants"][label] = {"scheme_robust": {s: sorted(names.get(d, d) for d in v) for s, v in robust.items()},
                                     "robust_in_all_schemes": sorted(names.get(d, d) for d in common)}
    (args.out / "track1_scenarios.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    for label, v in result["variants"].items():
        print(label)
        for s, d in v["scheme_robust"].items():
            print("  ", s, d)
        print("   in all schemes:", v["robust_in_all_schemes"])


if __name__ == "__main__":
    main()
