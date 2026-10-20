"""Manual review of paediatric label status (`config/pediatric_label_review.json`) applied over the heuristic extraction.

The heuristic in `labels.py` marks a label MIXED whenever it contains both "established" and "not established" phrases, so
it cannot tell "established for one indication" from "not established". For the drugs that matter, a reading of the label excerpt
is recorded in the review file; the review overrides the heuristic band profile, and drugs without a review keep the heuristic.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from track2.labels import BANDS

REVIEW_PATH = Path("config/pediatric_label_review.json")


def load_review(path: Path = REVIEW_PATH) -> dict[str, dict[str, Any]]:
    return json.loads(path.read_text(encoding="utf-8"))["entries"] if path.exists() else {}


def reviewed_profile(profile: dict[str, str], entry: dict[str, Any] | None) -> tuple[dict[str, str], str]:
    """Band profile after the review and its source: `reviewed` or `heuristic`."""
    if not entry or entry["verdict"] not in ("established", "not_established"):
        return profile, "heuristic"
    out = dict(profile)
    for band, (lo, hi) in BANDS.items():
        if entry["verdict"] == "established" and entry.get("min_age_years") is None:
            out[band] = "AGE_NOT_STATED"  # the label names children without an age: use is labelled, the age range is not
        else:
            out[band] = "LABELLED_USE" if entry["verdict"] == "established" and hi > entry["min_age_years"] else "NO_LABELLED_USE"
    return out, "reviewed"
