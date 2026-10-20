"""Age groups for the submission (0 to <3, 3 to <12, 12 and older) derived from the five label bands used in the methods.

The methods keep the five bands of `labels.BANDS` (neonate, infant, child, adolescent, adult). For the report the bands are grouped as
0-<3 years, 3-<12 years and 12 years and older (adolescents and adults). A group status is computed from the youngest age stated by the label
when it is known, and from the band profile otherwise:

- `labelled`   the label supports use across the whole group (youngest labelled age <= start of the group);
- `partly`     the label supports use from an age inside the group (the age is reported);
- `not_labelled` the label says paediatric use is not established (or the adult label only, for the 12+ group, is `adult_only`);
- `verify`     the label is mixed or could not be read, or its paediatric section only refers to paediatric information in another section without stating that paediatric use is established; verify by hand.

This is label support, not a recommendation and not dosing advice.
"""

from __future__ import annotations

from typing import Any

GROUPS = {"0_to_under_3": (0.0, 3.0), "3_to_under_12": (3.0, 12.0), "12_and_older": (12.0, 120.0)}
USE = {"LABELLED", "LABELLED_USE"}


def _by_min_age(min_age: float, lo: float, hi: float) -> str:
    if min_age <= lo:
        return "labelled"
    if min_age < hi:
        return f"partly (from {min_age:g} y)"
    return "not_labelled"


def group_status(profile: dict[str, str], min_age: float | None, verdict: str | None = None) -> dict[str, str]:
    """Status per group. `verdict` is the reviewed verdict (`established`, `not_established`) when the label was read by hand."""
    out: dict[str, str] = {}
    if any(profile.get(b) == "AGE_NOT_STATED" for b in ("neonate", "infant", "child", "adolescent")):
        return {name: "age_not_stated" for name in GROUPS}  # labelled for children, minimum age not stated
    for name, (lo, hi) in GROUPS.items():
        if verdict == "not_established":
            out[name] = "adult_only" if name == "12_and_older" and profile.get("adult") in USE else "not_labelled"
        elif min_age is not None and any(profile.get(b) in USE for b in ("neonate", "infant", "child", "adolescent")):
            out[name] = _by_min_age(min_age, lo, hi)
        else:
            bands = {"0_to_under_3": ("neonate", "infant"), "3_to_under_12": ("child",), "12_and_older": ("adolescent",)}[name]
            vals = [profile.get(b, "UNKNOWN") for b in bands]
            if all(v in USE for v in vals):
                out[name] = "labelled"
            elif any(v in USE for v in vals):
                out[name] = "partly"
            elif all(v == "NO_LABELLED_USE" for v in vals):
                out[name] = "adult_only" if name == "12_and_older" and profile.get("adult") in USE else "not_labelled"
            else:
                out[name] = "verify"
    return out


def summarise(rows: list[dict[str, Any]]) -> dict[str, dict[str, list[str]]]:
    """For each group, drug names by status class (labelled / partly / not_labelled_or_adult_only / verify)."""
    result: dict[str, dict[str, list[str]]] = {}
    for g in GROUPS:
        buckets: dict[str, list[str]] = {"labelled": [], "partly": [], "not_labelled_or_adult_only": [], "verify": []}
        for r in rows:
            st = r["age_groups"][g]
            key = "labelled" if st in ("labelled", "age_not_stated") else "partly" if st.startswith("partly") else "verify" if st == "verify" else "not_labelled_or_adult_only"
            buckets[key].append(r["name"])
        result[g] = {k: sorted(v) for k, v in buckets.items()}
    return result
