"""Descriptive toxicity profile from the FDA Adverse Event Reporting System (openFDA, CC0).

For each drug: total number of spontaneous reports and the share of its 100 most reported reaction terms
(MedDRA preferred terms) that fall in organ-system groups relevant to the proband (renal, hepatic,
haematological, cardiac, neurological, hypersensitivity). This is **descriptive only**:
spontaneous reports have no denominator, are confounded by the indication (for example the disease treated
appears among the reactions) and cannot give an incidence or compare drugs. It complements, and never replaces,
the label text and a clinical review.
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any

from track2.http_cache import get_json

OPENFDA = "https://api.fda.gov/drug/event.json"
TOP_TERMS = 100
GROUPS = {
    "renal": r"(renal|kidney|nephro|creatinine|proteinuria|hypercalc|calcul)",
    "hepatic": r"(hepat|liver|transaminase|bilirubin|jaundice)",
    "haematological": r"(neutropeni|thrombocytopeni|anaemi|anemi|leukopeni|pancytopeni|myelosuppress|cytopeni)",
    "cardiac": r"(cardiac|cardio|arrhythm|qt |heart|tachycard|bradycard)",
    "neurological": r"(neuropath|seizure|encephalopath|neurotox)",
    "hypersensitivity": r"(hypersensitiv|anaphyla|angioedema|rash|urticaria|stevens|epidermal)",
}


def _search(name: str) -> str:
    return f'patient.drug.openfda.generic_name.exact:"{name.upper()}"'


_SAFE = ':+"'


def _url(name: str, extra: str) -> str:
    query = urllib.parse.quote(_search(name), safe=_SAFE)
    return f"{OPENFDA}?search={query}{extra}"


def classify_terms(terms: list[dict[str, Any]]) -> dict[str, float]:
    """Share of the reports of the listed terms that belong to each organ-system group."""
    total = sum(t["count"] for t in terms) or 1
    shares = {}
    for group, pattern in GROUPS.items():
        hit = sum(t["count"] for t in terms if re.search(pattern, t["term"].lower()))
        shares[group] = round(hit / total, 4)
    return shares


def toxicity_profile(drug_name: str) -> dict[str, Any]:
    """openFDA FAERS summary for a drug name; `available = False` when openFDA has no reports under that name."""
    name = re.sub(r"\b(ANHYDROUS|HYDROCHLORIDE|CITRATE|MESYLATE|SULFATE|PHOSPHATE)\b", "", drug_name.upper()).split(",")[0].strip()
    try:
        total = get_json(_url(name, "&limit=1"))["meta"]["results"]["total"]
        terms = get_json(_url(name, f"&count=patient.reaction.reactionmeddrapt.exact&limit={TOP_TERMS}"))["results"]
    except LookupError:
        return {"available": False, "source": "openFDA FAERS", "reason": "openFDA has no reports under this name"}
    return {
        "available": True, "source": "openFDA FAERS (CC0)", "name_queried": name, "total_reports": total,
        "top_terms": [{"term": t["term"], "count": t["count"]} for t in terms[:10]],
        "group_shares_of_top_terms": classify_terms(terms),
        "caveat": "Spontaneous reports: no denominator, no incidence, confounded by indication; descriptive only.",
    }
