"""Drug labels (DailyMed/FDA): paediatric use, minimum age and safety flags by keyword.

The extraction is **heuristic**: it returns the label excerpt and a conservative classification that
must be verified by hand for the final candidates. If no label is found, the status is UNKNOWN (neither approved
nor rejected).
"""

from __future__ import annotations

import re
import urllib.parse
import xml.etree.ElementTree as ET
from typing import Any

from track2.http_cache import CACHE_DIR, get_bytes, get_json

DAILYMED = "https://dailymed.nlm.nih.gov/dailymed/services/v2"
H = "{urn:hl7-org:v3}"
SECTIONS = {"pediatric": "34081-0", "contraindications": "34070-3", "warnings": "43685-7", "boxed": "34066-1", "indications": "34067-9"}
BANDS = {"neonate": (0.0, 1 / 12), "infant": (1 / 12, 2.0), "child": (2.0, 12.0), "adolescent": (12.0, 18.0)}
KEYWORDS = {
    "RENAL": r"(renal|kidney|nephro)[^.]{0,100}(toxicit|failure|impairment|injur|insufficien)",
    "HEPATIC": r"(hepat|liver)[^.]{0,100}(toxicit|failure|impairment|injur|enzyme)",
    "HEMATOLOGIC": r"(myelosuppress|neutropeni|thrombocytopeni|cytopeni|anemi)",
    "EMBRYO_FETAL": r"(embryo|fetal)[^.]{0,60}(toxicit|harm)",
    "NEPHROCALCINOSIS_OR_STONES": r"(nephrocalcinosis|nephrolithiasis|kidney stone|renal calcul|hypercalciuria|hypercalcemia)",
    "HYPERSENSITIVITY": r"(hypersensitivity|anaphyla|angioedema|stevens-johnson|toxic epidermal)",
    "CARDIAC": r"(qt prolongation|cardiotox|cardiomyopath|arrhythmia|heart failure)",
    "NEUROLOGIC": r"(neurotox|peripheral neuropathy|seizure)",
    "PULMONARY": r"(pulmonary toxicity|interstitial lung|pneumonitis)",
}
# Author-curated list of excipients with well-known label-documented concerns; it is a prompt for pharmacist
# review, not an exhaustive allergen database.
EXCIPIENT_ALERTS = {
    "BENZYL ALCOHOL": "gasping syndrome in neonates and premature infants",
    "PROPYLENE GLYCOL": "accumulation in neonates and in renal impairment",
    "POLYSORBATE 80": "hypersensitivity reactions reported",
    "POLYOXYL": "polyoxyl castor oil vehicles are associated with anaphylactoid reactions",
    "CASTOR OIL": "castor oil vehicles are associated with anaphylactoid reactions",
    "SULFITE": "sulfite sensitivity",
    "METABISULFITE": "sulfite sensitivity",
    "LACTOSE": "lactose intolerance",
    "WHEAT": "gluten-containing starch",
    "PEANUT": "peanut allergy",
    "SOYBEAN": "soy allergy",
    "TARTRAZINE": "dye hypersensitivity",
    "FD&C YELLOW NO. 5": "dye hypersensitivity",
    "CYCLODEXTRIN": "renal accumulation of the vehicle in renal impairment",
    "BETADEX": "sulfobutylether-beta-cyclodextrin vehicle: renal accumulation in renal impairment",
    "SULFOBUTYL": "sulfobutylether-beta-cyclodextrin vehicle: renal accumulation in renal impairment",
    "GELATIN": "animal-derived excipient",
    "LATEX": "latex allergy",
}


def _clean_name(name: str) -> str:
    return re.sub(r"\b(ANHYDROUS|HYDROCHLORIDE|SODIUM|CITRATE|MESYLATE|SULFATE|PHOSPHATE)\b", "", name.upper()).split(",")[0].strip()


def find_setid(name: str) -> tuple[str, str] | None:
    base = _clean_name(name).lower()
    if not base:
        return None
    data = get_json(f"{DAILYMED}/spls.json?drug_name={urllib.parse.quote(base)}&pagesize=25")
    for item in data.get("data", []):
        if base.split()[0] in item["title"].lower():
            return item["setid"], item["title"]
    return None


def _spl_xml(setid: str) -> ET.Element:
    path = CACHE_DIR / f"spl_{setid}.xml"
    if not path.exists():
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        path.write_bytes(get_bytes(f"{DAILYMED}/spls/{setid}.xml", accept="*/*"))
    return ET.fromstring(path.read_bytes())


def _sections(root: ET.Element) -> dict[str, str]:
    out: dict[str, list[str]] = {k: [] for k in SECTIONS}
    inverse = {v: k for k, v in SECTIONS.items()}
    for section in root.iter(H + "section"):
        code = section.find(H + "code")
        if code is not None and code.get("code") in inverse:
            out[inverse[code.get("code")]].append(re.sub(r"\s+", " ", "".join(section.itertext())).strip())
    return {k: " ".join(v) for k, v in out.items()}


def classify_pediatric(text: str) -> dict[str, Any]:
    if not text:
        return {"status": "UNKNOWN", "min_age_years": None, "excerpt": ""}
    t = text.lower()
    established = bool(re.search(r"(safety and (effectiveness|efficacy)[^.]*?(have|has) been established)|(approved for[^.]*pediatric)|(indicated[^.]*pediatric)|(pediatric patients[^.]*(\d+)\s*(years|months))", t))
    not_established = bool(re.search(r"(not been established|not established|not recommended|contraindicated in (pediatric|children))", t))
    status = "MIXED" if established and not_established else "ESTABLISHED" if established else "NOT_ESTABLISHED" if not_established else "UNKNOWN"
    ages = []
    for number, unit in re.findall(r"(\d+(?:\.\d+)?)\s*(years?|months?|weeks?)\s*(?:of age|and older|or older|to|and)", t):
        value = float(number)
        ages.append(value / 12 if unit.startswith("month") else value / 52 if unit.startswith("week") else value)
    return {"status": status, "min_age_years": min(ages) if ages else None, "excerpt": text[:400]}


def band_for_age(years: float) -> str:
    """Age band of a patient (the same bands the label profile uses; 18 years and above is `adult`)."""
    for band, (lo, hi) in BANDS.items():
        if lo <= years < hi:
            return band
    return "adult"


def band_profile(pediatric: dict[str, Any], has_label: bool) -> dict[str, str]:
    """Label status per age band (says nothing about this patient)."""
    profile = {"adult": "LABELLED" if has_label else "UNKNOWN"}
    for band, (lo, hi) in BANDS.items():
        if pediatric["status"] == "MIXED":
            profile[band] = "VERIFY_MANUALLY"  # the label combines established and non-established uses
        elif pediatric["status"] == "ESTABLISHED" and pediatric["min_age_years"] is not None:
            profile[band] = "LABELLED_USE" if hi > pediatric["min_age_years"] else "NO_LABELLED_USE"
        elif pediatric["status"] == "ESTABLISHED":
            profile[band] = "LABELLED_AGE_NOT_PARSED"
        elif pediatric["status"] == "NOT_ESTABLISHED":
            profile[band] = "NO_LABELLED_USE"
        else:
            profile[band] = "UNKNOWN"
    return profile


def inactive_ingredients(root: ET.Element) -> list[str]:
    names = []
    for ingredient in root.iter(H + "ingredient"):
        if ingredient.get("classCode") == "IACT":
            node = ingredient.find(H + "ingredientSubstance/" + H + "name")
            if node is not None and node.text:
                names.append(node.text.strip().upper())
    return sorted(set(names))


EXCIPIENT_EXACT = {  # whole-name matches only: "polyvinyl alcohol" or "isopropyl alcohol" must not trigger the ethanol alert
    "ALCOHOL": "ethanol content, relevant in children", "DEHYDRATED ALCOHOL": "ethanol content, relevant in children",
    "ALCOHOL, DEHYDRATED": "ethanol content, relevant in children", "ETHANOL": "ethanol content, relevant in children",
}


def excipient_alerts(names: list[str]) -> list[dict[str, str]]:
    alerts = []
    for name in names:
        if name in EXCIPIENT_EXACT:
            alerts.append({"excipient": name, "concern": EXCIPIENT_EXACT[name]})
            continue
        for key, reason in EXCIPIENT_ALERTS.items():
            if key in name:
                alerts.append({"excipient": name, "concern": reason})
                break
    return alerts


def keyword_excerpts(text: str, flags: list[str]) -> dict[str, str]:
    """First label sentence that triggered each keyword flag (for the evidence trail)."""
    out = {}
    for flag in flags:
        match = re.search(KEYWORDS[flag], text, re.I)
        if match:
            start = max(0, match.start() - 90)
            out[flag] = text[start:match.end() + 90].strip()
    return out


def label_info(name: str) -> dict[str, Any]:
    empty = classify_pediatric("")
    found = find_setid(name)
    if not found:
        return {"found": False, "pediatric": empty, "age_profile": band_profile(empty, False), "keyword_flags": [],
                "keyword_excerpts": {}, "contraindicated_hypersensitivity": False, "inactive_ingredients": [], "excipient_alerts": []}
    setid, title = found
    root = _spl_xml(setid)
    sections = _sections(root)
    risk_text = " ".join([sections["boxed"], sections["contraindications"], sections["warnings"]])
    flags = sorted(k for k, pattern in KEYWORDS.items() if re.search(pattern, risk_text.lower()))
    pediatric = classify_pediatric(sections["pediatric"])
    inactive = inactive_ingredients(root)
    return {
        "found": True, "setid": setid, "title": title, "pediatric": pediatric, "keyword_flags": flags,
        "keyword_excerpts": keyword_excerpts(risk_text, flags),
        "contraindicated_hypersensitivity": bool(re.search(KEYWORDS["HYPERSENSITIVITY"], sections["contraindications"].lower())),
        "inactive_ingredients": inactive, "excipient_alerts": excipient_alerts(inactive),
        "has_boxed_warning": bool(sections["boxed"]), "age_profile": band_profile(pediatric, True),
        "heuristic": True,
    }
