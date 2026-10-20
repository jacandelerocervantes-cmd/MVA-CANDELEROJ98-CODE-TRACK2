"""Paediatric-labelled arm (rule: `protocol/addendum_pediatric.md`).

    python -m track2.pediatric_arm --out out

Pool = the 300 best drugs by combined score (STRING 700, epcr weights); keep those whose US label states established paediatric
use (`ESTABLISHED`) or mixed use to verify (`MIXED`). Each is reported with its rank in the nine configurations of
`addendum_robustness.md`; stable = rank within the best 100 in at least 7 of 9. These drugs were selected by label, not by
distance, so their proximity is weaker than for the main candidates and the report shows it. Not a recommendation or a dose.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import statistics
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from track2.consensus import SCHEMES, THRESHOLDS, combined, weights
from track2.labels import label_info
from track2.literature_context import classify
from track2.pediatric_review import load_review

POOL_N = 300
STABLE_RANK, STABLE_MIN = 100, 7


def ranks(scores: dict[str, float]) -> dict[str, int]:
    order = sorted(scores, key=lambda d: (-scores[d], d))
    return {d: i + 1 for i, d in enumerate(order)}


def run(out: Path, mechanism_dir: Path, universe: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    mechanisms = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(mechanism_dir.glob("*.json"))]
    mechanisms = [m for m in mechanisms if m.get("status") == "mechanism_candidate" and (out / f"threshold_scores_{m['id']}.json").exists()]
    by_thr = {m["id"]: json.loads((out / f"threshold_scores_{m['id']}.json").read_text(encoding="utf-8")) for m in mechanisms}
    wts = weights(mechanisms)
    table = {(t, s): combined({mid: by_thr[mid][t] for mid in by_thr}, wts[s]) for t, s in itertools.product(THRESHOLDS, SCHEMES)}
    rank_tables = {k: ranks({d: v for d, v in sc.items()}) for k, sc in table.items()}
    main = table[("700", "W1_epcr")]
    pool = [d for d in sorted(main, key=lambda d: (-main[d], d)) if main[d] > 0][:POOL_N]
    drugs = {d["drug_id"]: d for d in universe["drugs"]}
    review = load_review()
    with ThreadPoolExecutor(max_workers=4) as pool_ex:
        labels = dict(zip(pool, pool_ex.map(lambda d: label_info(drugs[d]["name"] or ""), pool)))
    flags: dict[str, list[str]] = {}
    for m in mechanisms:
        for r in json.loads((out / f"scores_{m['id']}.json").read_text(encoding="utf-8"))["rows"]:
            flags.setdefault(r["drug_id"], [])
            flags[r["drug_id"]] = sorted(set(flags[r["drug_id"]]) | set(r["risk_flags"]))
    rows = []
    for d in pool:
        lab = labels[d]
        ped = lab["pediatric"]
        if ped["status"] not in ("ESTABLISHED", "MIXED"):
            continue
        per_cfg = {f"{t}|{s}": rank_tables[(t, s)].get(d) for t, s in rank_tables}
        valid = [r for r in per_cfg.values() if r is not None]
        n_good = sum(1 for r in valid if r <= STABLE_RANK)
        moa = sorted({m["mechanism"] for m in drugs[d]["mechanisms"] if m["mechanism"]})
        z700 = {mid: by_thr[mid]["700"].get(d) for mid in by_thr}
        zs = [v for v in z700.values() if v is not None]
        entry = review.get(d)
        final = ("established_reviewed" if entry and entry["verdict"] == "established" else
                 "not_established_reviewed" if entry and entry["verdict"] == "not_established" else "heuristic_unreviewed")
        rows.append({
            "review": entry, "final_label_status": final, "min_z_700": round(min(zs), 2) if zs else None,
            "meaningful_proximity": bool(zs) and min(zs) <= -2.0,
            "headline": final == "established_reviewed" and bool(zs) and min(zs) <= -2.0,
            "drug_id": d, "name": drugs[d]["name"], "combined_700_W1": round(main[d], 4), "rank_700_W1": rank_tables[("700", "W1_epcr")][d],
            "z_700_by_module": {k: (round(v, 2) if v is not None else None) for k, v in z700.items()},
            "ranks_by_configuration": per_cfg, "median_rank": statistics.median(valid), "n_configs_rank_le_100": n_good,
            "stability": "stable" if n_good >= STABLE_MIN else "unstable",
            "pediatric_label_status": ped["status"], "min_age_years": ped["min_age_years"], "age_bands": lab["age_profile"],
            "pediatric_excerpt": ped["excerpt"], "mechanism_of_action": moa, "risk_flags": flags.get(d, []),
            "excipient_alerts": [a["excipient"] for a in lab.get("excipient_alerts", [])],
            "literature_context": [{"class": e["class"], "direction_hint": e["direction_hint"], "pmids": e["pmids"]}
                                   for e in classify({"mechanism_of_action": moa, "risk_flags": flags.get(d, [])}, config)],
        })
    rows.sort(key=lambda r: (r["median_rank"], r["drug_id"]))
    return {"rule": {"pool": POOL_N, "stable_rank": STABLE_RANK, "stable_min_configs": STABLE_MIN}, "pool_size": len(pool),
            "n_with_label": sum(1 for d in pool if labels[d]["found"]), "n_paediatric": len(rows),
            "n_stable": sum(1 for r in rows if r["stability"] == "stable"),
            "n_headline": sum(1 for r in rows if r["headline"]), "headline": [r["name"] for r in rows if r["headline"]], "drugs": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--config", type=Path, default=Path("config/literature_context.json"))
    args = parser.parse_args()
    result = run(args.out, args.mechanism_dir, json.loads(args.universe.read_text(encoding="utf-8")), json.loads(args.config.read_text(encoding="utf-8")))
    (args.out / "pediatric_arm.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(args.out / "pediatric_arm.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["drug_id", "name", "headline", "final_label_status", "min_z_700", "stability", "median_rank", "n_configs_rank_le_100", "combined_700_W1", "pediatric_label_status", "min_age_years", "risk_flags"])
        for r in result["drugs"]:
            writer.writerow([r["drug_id"], r["name"], r["headline"], r["final_label_status"], r["min_z_700"], r["stability"], r["median_rank"], r["n_configs_rank_le_100"], r["combined_700_W1"],
                             r["pediatric_label_status"], r["min_age_years"] if r["min_age_years"] is not None else "", ";".join(r["risk_flags"])])
    print({k: v for k, v in result.items() if k != "drugs"})


if __name__ == "__main__":
    main()
