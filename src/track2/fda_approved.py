"""FDA eligibility of a drug (rule: `protocol/addendum_plausibility.md`): an openFDA Drugs@FDA record whose active ingredient matches and at
least one product with marketing status `Prescription` (currently marketed). Discontinued-only or unmatched drugs are not eligible.

A service error raises; it is never read as "not eligible" or "eligible".
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any, Callable

from track2.http_cache import get_json

URL = "https://api.fda.gov/drug/drugsfda.json?search={q}&limit=100"
SALT_WORDS = {"ANHYDROUS", "HYDROCHLORIDE", "SODIUM", "CITRATE", "MESYLATE", "SULFATE", "MALEATE", "TARTRATE", "ACETATE", "PHOSPHATE", "BROMIDE", "POTASSIUM",
              "CALCIUM", "FUMARATE", "SUCCINATE", "HYDROBROMIDE", "MONOHYDRATE", "DIHYDRATE", "TOSYLATE", "BESYLATE", "NITRATE", "CHLORIDE", "DISODIUM"}


def candidate_names(name: str) -> list[str]:
    """The name as given, then without salt or hydrate words, then the first word (ordered, no duplicates)."""
    upper = re.sub(r"\s+", " ", name.upper().strip())
    stripped = " ".join(w for w in upper.split() if w not in SALT_WORDS)
    out = [upper]
    for n in (stripped, upper.split()[0] if upper else ""):
        if n and n not in out:
            out.append(n)
    return out


def lookup(name: str, get: Callable[[str], Any] = get_json) -> dict[str, Any]:
    for candidate in candidate_names(name):
        query = urllib.parse.quote(f'products.active_ingredients.name:"{candidate}"')
        try:
            data = get(URL.format(q=query))
        except LookupError:  # openFDA answers 404 when nothing matches
            continue
        statuses = sorted({p.get("marketing_status") or "None" for r in data["results"] for p in r["products"]})
        if "Prescription" in statuses:
            apps = sorted({r["application_number"] for r in data["results"]})
            return {"eligible": True, "matched_name": candidate, "marketing_statuses": statuses, "applications": apps[:5], "n_applications": len(apps)}
        return {"eligible": False, "matched_name": candidate, "marketing_statuses": statuses, "applications": [], "n_applications": 0,
                "reason": "no product with marketing status Prescription"}
    return {"eligible": False, "matched_name": None, "marketing_statuses": [], "applications": [], "n_applications": 0, "reason": "no Drugs@FDA record matches the active ingredient"}
