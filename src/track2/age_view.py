"""Candidates by age band: the patient's age is unknown, so every candidate is shown for every band.

    python -m track2.age_view --out out                      # all bands (the case today: age not documented)
    python -m track2.age_view --out out --age-years 7        # the same tables plus a view for the child band
    python -m track2.age_view --out out --age-band infant

Reads the main candidates (`candidates_enriched.json`, with robustness) and the extrapolation arm
(`arm_extrapolation/candidates.json`). The label status per band comes from the heuristic DailyMed extraction and must be
verified by hand. This is not dosing advice: it only says whether a US label supports use in an age band.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from track2.labels import BANDS, band_for_age
from track2.pediatric_review import load_review, reviewed_profile

ALL_BANDS = list(BANDS) + ["adult"]
USE = {"LABELLED", "LABELLED_USE", "AGE_NOT_STATED"}
NO_USE = {"NO_LABELLED_USE"}


def band_group(status: str) -> str:
    return "labelled_use" if status in USE else "no_labelled_use" if status in NO_USE else "unknown_or_verify"


def view(candidates: list[dict[str, Any]], arm: str, review: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        robustness = c.get("layers", {}).get("robustness", {}).get("class", "")
        age, source = reviewed_profile(c["age_profile"], (review or {}).get(c["drug_id"]))
        rows.append({"arm": arm, "drug_id": c["drug_id"], "name": c["name"], "robustness": robustness, "label_source": source,
                     "bands": {b: age.get(b, "UNKNOWN") for b in ALL_BANDS},
                     "excipient_alerts": [a["excipient"] for a in c["excipients"]["alerts"]],
                     "risk_flags": c["risk_flags"]})
    return rows


def by_band(rows: list[dict[str, Any]]) -> dict[str, dict[str, list[str]]]:
    out: dict[str, dict[str, list[str]]] = {}
    for band in ALL_BANDS:
        groups: dict[str, list[str]] = {"labelled_use": [], "no_labelled_use": [], "unknown_or_verify": []}
        for r in rows:
            groups[band_group(r["bands"][band])].append(r["name"])
        out[band] = {k: sorted(v) for k, v in groups.items()}
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--age-years", type=float)
    group.add_argument("--age-band", choices=ALL_BANDS)
    args = parser.parse_args()
    reg = json.loads((args.out / "candidate_register.json").read_text(encoding="utf-8"))["drugs"]
    rows = [{"arm": ";".join(r["sources"]), "drug_id": r["drug_id"], "name": r["name"], "bands": r["age_bands"], "excipient_alerts": r.get("excipient_alerts") or [],
             "boxed_warning": bool(r.get("boxed_warning"))} for r in reg]
    selected = args.age_band or (band_for_age(args.age_years) if args.age_years is not None else None)
    result: dict[str, Any] = {
        "note": "Label support per age band, from the openFDA / DailyMed label (heuristic extraction, replaced by the hand-read verdicts of config/pediatric_label_review.json where present). Not dosing advice. "
                "The patient's age is not documented, so all bands are shown.",
        "selected_band": selected, "register": by_band(rows),
    }
    (args.out / "age_views.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(args.out / "age_views.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["sources", "drug_id", "name", *ALL_BANDS, "excipient_alerts", "boxed_warning"])
        for r in rows:
            writer.writerow([r["arm"], r["drug_id"], r["name"], *[r["bands"].get(b, "UNKNOWN") for b in ALL_BANDS], ";".join(r["excipient_alerts"]), r["boxed_warning"]])
    if selected:
        chosen = [r for r in rows if band_group(r["bands"].get(selected, "UNKNOWN")) != "no_labelled_use"]
        with open(args.out / f"age_view_{selected}.csv", "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["sources", "drug_id", "name", "label_status", "excipient_alerts"])
            for r in chosen:
                writer.writerow([r["arm"], r["drug_id"], r["name"], r["bands"].get(selected, "UNKNOWN"), ";".join(r["excipient_alerts"])])
        print(f"{selected}: {len(chosen)} candidates are not labelled-as-unsupported for this band")
    print({b: {k: len(v) for k, v in g.items()} for b, g in result["register"].items()})


if __name__ == "__main__":
    main()
