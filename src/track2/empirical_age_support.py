"""Empirical age support of the US label under an unknown patient age (0-18 years).

LABEL SUPPORT, NOT A RECOMMENDATION. The age of the patient is unknown and is never assumed. For each drug of the paediatric label table:

- L_d(a), a in [0, 18]: 1 where the label states use (`established`: from `min_age_years`, or the whole range when no age is stated),
  0 where the label says paediatric use is not established or the label is adult only, and unknown for the verdict `verify`.
- H_d(a): 1 at ages covered by an age-specific caution found by regular expressions in the boxed warning, warnings, contraindications,
  paediatric use and dosage sections ('under 6 years', 'less than 13 years', 'neonates', ...). Each kept sentence is a verbatim substring of
  the label text (whitespace normalised). Weight cut-offs ('up to 20 kg') are kept as flagged text and are NOT converted to ages.
- C_d = measure{a in [0, 18]: L_d(a) = 1 and H_d(a) = 0} / 18, as a lower and an upper bound (unknown L: 0 for the lower, 1 for the upper bound),
  and the same fraction inside each age group of `age_groups.GROUPS` (clipped at 18 years).

No dose is computed or recommended. The regular expressions are heuristic: a caution about one indication is applied to the whole range.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

from track2 import age_groups, labels

MAX_AGE = 18.0
HEADER = {
    "purpose": "Label support under an unknown patient age: the fraction of the paediatric age range 0-18 years over which the US label supports use and no age-specific label warning contradicts it.",
    "not_a_recommendation": "This is NOT a recommendation of treatment and NOT dosing advice. No dose is computed. The age of the patient is unknown and is never assumed.",
    "definition": "C = measure{a in [0,18] : L(a)=1 and H(a)=0} / 18; L from the reviewed label verdict (verify = unknown: lower bound 0, upper bound 1); H from age-specific cautions extracted by regular expressions (verbatim sentences). Weight cut-offs are flagged, not converted.",
    "heuristic": "Extraction is by regular expressions over label sections and must be verified by hand; a caution about one indication is applied to the whole age range.",
}

SECTION_CODES = {
    "boxed_warning": ("34066-1",),
    "warnings": ("43685-7", "34071-1", "42232-9", "43684-0"),
    "contraindications": ("34070-3",),
    "pediatric_use": ("34081-0",),
    "dosage": ("34068-7",),
}

# ------------------------------------------------------------------ L_d --------------------------------------------------------


def label_support(verdict: str, min_age: float | None) -> dict[str, Any]:
    """Step function L(a) on [0, 18]: {'state': 'known' | 'unknown', 'intervals': [[lo, hi], ...] where L = 1}."""
    if verdict == "verify":
        return {"state": "unknown", "intervals": []}
    if verdict == "established":
        start = 0.0 if min_age is None else max(0.0, float(min_age))
        return {"state": "known", "intervals": [[start, MAX_AGE]] if start < MAX_AGE else []}
    return {"state": "known", "intervals": []}  # not_established, no_us_label or anything else: no stated use


# ------------------------------------------------------------------ intervals --------------------------------------------------


def merge(intervals: list[list[float]]) -> list[list[float]]:
    out: list[list[float]] = []
    for lo, hi in sorted([max(0.0, a), min(MAX_AGE, b)] for a, b in intervals if min(MAX_AGE, b) > max(0.0, a)):
        if out and lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


def measure(intervals: list[list[float]], lo: float = 0.0, hi: float = MAX_AGE) -> float:
    return sum(max(0.0, min(b, hi) - max(a, lo)) for a, b in merge(intervals))


def _intersection_measure(a: list[list[float]], b: list[list[float]], lo: float, hi: float) -> float:
    total = 0.0
    for x0, x1 in merge(a):
        for y0, y1 in merge(b):
            total += max(0.0, min(x1, y1, hi) - max(x0, y0, lo))
    return total


def supported_measure(L_intervals: list[list[float]], H_intervals: list[list[float]], lo: float = 0.0, hi: float = MAX_AGE) -> float:
    """measure{a in [lo, hi]: L = 1 and H = 0}."""
    return measure(L_intervals, lo, hi) - _intersection_measure(L_intervals, H_intervals, lo, hi)


def support_bounds(L: dict[str, Any], H_intervals: list[list[float]], lo: float = 0.0, hi: float = MAX_AGE) -> dict[str, float]:
    """Lower/upper fraction of [lo, hi] supported (unknown L: lower 0, upper 1 before removing H)."""
    hi = min(hi, MAX_AGE)
    width = hi - lo
    if L["state"] == "unknown":
        lower, upper = 0.0, supported_measure([[0.0, MAX_AGE]], H_intervals, lo, hi)
    else:
        lower = upper = supported_measure(L["intervals"], H_intervals, lo, hi)
    return {"lower": round(lower / width, 6), "upper": round(upper / width, 6), "width_years": width}


# ------------------------------------------------------------------ cautions ---------------------------------------------------

NUM = r"(\d+(?:\.\d+)?)"
UNIT = r"(years?|yrs?|months?|weeks?|days?)\b"
AGE_CONTEXT = re.compile(r"\b(age|aged|old|child|children|pediatric|paediatric|infant|infants|neonat\w*|newborn\w*|adolescent\w*|younger|juvenile|patients?|premature|preterm)\b", re.I)
NOT_AGE_AFTER = r"(?!\s*(?:of\s+(?:treatment|therapy|dosing|exposure|use|prior|gestation|gestational|continuous)|gestation|gestational|after|post|since|ago|prior|duration|of\s+the\s+treatment))"
LT = r"(?:under|below|less\s+than|younger\s+than|lower\s+than|fewer\s+than|before|<|&lt;)"
CUE = re.compile(
    r"not\s+recommended|contraindicat|should\s+not|must\s+not|do\s+not|avoid|not\s+indicated|not\s+(?:been\s+)?established|not\s+approved|not\s+been\s+studied|not\s+studied"
    r"|have\s+not\s+been|has\s+not\s+been|risk\s+of|fatal|death|gasping|toxicit|more\s+frequent|higher\s+(?:incidence|rate)|greater\s+(?:risk|incidence|sensitivity)"
    r"|increased\s+(?:risk|incidence)|caution|warning|not\s+to\s+be\s+used|safety\s+and\s+(?:effectiveness|efficacy)[^.]{0,40}not", re.I)
RANGE_BEFORE = re.compile(r"\d(?:\.\d+)?\s*(?:years?|months?|weeks?|days?)?\s*(?:to|-|–|and)\s*$", re.I)
ANIMAL = re.compile(r"\b(?:rats?|dogs?|mice|mouse|monkeys?|animals?|rabbits?|preclinical|juvenile\s+(?:toxicity|animal))\b", re.I)
MATERNAL = re.compile(r"born\s+to|pregnan|maternal|\bwomen\b|breast|lactat|in\s+utero|fetal|fetus", re.I)
MAX_SENTENCE = 700
DOSE_TABLE = re.compile(r"mg/m|m2|\bBSA\b|Recommended Daily Dosage", re.I)
TIME_NOT_AGE = re.compile(r"\d+\s*(?:years?|months?|weeks?)\s*(?:previously|ago|post)|\d+\s*(?:years?|months?|weeks?)\s+previously|post-", re.I)
STRONG_RANGE_CUE = re.compile(r"not\s+recommended|contraindicat|should\s+not(?!\s+be\s+less)|must\s+not|avoid|not\s+indicated|fatal|death|gasping|greater\s+risk|increased\s+risk", re.I)
POSITIVE_ONLY = re.compile(r"(?:have|has|been)\s+(?:been\s+)?established", re.I)
NEGATION = re.compile(r"\bnot\b[^.]{0,30}established|not\s+recommended|contraindicat|should\s+not|avoid|risk|death|fatal|toxicit|caution|warning|gasping|more\s+frequent", re.I)


def _years(value: str, unit: str) -> float:
    v = float(value)
    u = unit.lower()
    if u.startswith("month"):
        return v / 12.0
    if u.startswith("week"):
        return v / 52.1429
    if u.startswith("day"):
        return v / 365.25
    return v


def _age_intervals(sentence: str) -> list[tuple[list[float], str]]:
    """(interval, matched text) for each age cut-off found in the sentence."""
    found: list[tuple[list[float], str]] = []
    s = sentence
    for m in re.finditer(rf"{LT}\s*(?:the\s+age\s+of\s+)?{NUM}\s*{UNIT}{NOT_AGE_AFTER}", s, re.I):
        if RANGE_BEFORE.search(s[: m.start()]):
            continue  # 'aged 5 to <17 years' describes a studied range, not a cut-off
        found.append(([0.0, _years(m.group(1), m.group(2))], m.group(0)))
    for m in re.finditer(rf"(?<![\w.]){NUM}\s*(?:-\s*)?{UNIT}\s+(?:of\s+age\s+)?(?:and|or)\s+(?:younger|under|below|less)", s, re.I):
        found.append(([0.0, _years(m.group(1), m.group(2))], m.group(0)))
    for m in re.finditer(rf"up\s+to\s+(?:the\s+age\s+of\s+)?{NUM}\s*{UNIT}\s+(?:of\s+age|old)", s, re.I):
        found.append(([0.0, _years(m.group(1), m.group(2))], m.group(0)))
    for m in re.finditer(r"≤\s*" + NUM + r"\s*" + UNIT, s, re.I):
        found.append(([0.0, _years(m.group(1), m.group(2))], m.group(0)))
    if STRONG_RANGE_CUE.search(s):
        for m in re.finditer(rf"(?:aged?\s+|patients\s+|children\s+|infants\s+)?{NUM}\s*(?:{UNIT})?\s*(?:to|-|–|and)\s*{NUM}\s*{UNIT}(?=\s+(?:of\s+age|old)|[\s,.;)])", s, re.I):
            unit1 = m.group(2) or m.group(4)
            lo, hi = _years(m.group(1), unit1), _years(m.group(3), m.group(4))
            if lo < hi and re.search(r"age|aged|old|child|patient|infant", s[max(0, m.start() - 40): m.end() + 40], re.I):
                found.append(([lo, hi], m.group(0)))
    return found


def _neonate_interval(sentence: str) -> list[tuple[list[float], str]]:
    if MATERNAL.search(sentence):
        return []
    return [([0.0, 28 / 365.25], m.group(0)) for m in re.finditer(r"\b(neonat\w*|newborn\w*|premature\s+(?:infant|neonat)\w*|preterm\s+\w+)", sentence, re.I)][:1]


def _weight_cutoffs(sentence: str) -> list[dict[str, Any]]:
    out = []
    pats = [rf"weight\s+of\s+{NUM}\s*kg\b", rf"(?:up\s+to|less\s+than|under|below|not\s+exceeding|<|≤)\s*{NUM}\s*kg\b", rf"{NUM}\s*kg\s*(?:or\s+less|and\s+under|or\s+below|or\s+lighter)"]
    for p in pats:
        for m in re.finditer(r"(?<![/\w])" + p, sentence, re.I):
            if float(m.group(1)) <= 60:
                out.append({"cutoff_kg": float(m.group(1)), "matched_text": m.group(0)})
    return out


def sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?])\s+(?=[A-Z0-9\[(“\"])", text) if s]


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _kind(sentence: str) -> str:
    if re.search(r"contraindicat", sentence, re.I):
        return "contraindicated"
    if re.search(r"not\s+recommended|should\s+not|must\s+not|do\s+not|avoid|not\s+indicated|not\s+to\s+be\s+used", sentence, re.I):
        return "not_recommended"
    if re.search(r"not\s+(?:been\s+)?established|not\s+approved|not\s+(?:been\s+)?studied|have\s+not\s+been|has\s+not\s+been", sentence, re.I):
        return "not_established"
    return "risk_or_warning"


def extract_cautions(sections: dict[str, str]) -> list[dict[str, Any]]:
    """Age- and weight-specific cautions with their verbatim sentence. Every kept sentence is asserted to be a substring of its section."""
    out: list[dict[str, Any]] = []
    seen: set[str] = set()
    for name in sorted(sections):
        text = _norm(sections[name] or "")
        for sent in sentences(text):
            if sent in seen or len(sent) > MAX_SENTENCE or ANIMAL.search(sent) or not CUE.search(sent) or not AGE_CONTEXT.search(sent) and not _weight_cutoffs(sent):
                continue
            if POSITIVE_ONLY.search(sent) and not NEGATION.search(sent):
                continue
            ages = [(iv, t) for iv, t in (_age_intervals(sent) if AGE_CONTEXT.search(sent) else []) if iv[1] <= MAX_AGE]  # >18 y: adult comparison, not paediatric
            if CUE.search(sent) and AGE_CONTEXT.search(sent) and not ages:
                ages = _neonate_interval(sent)
            weights = _weight_cutoffs(sent)
            if not ages and not weights:
                continue
            assert sent in text, f"sentence not verbatim in section {name}"
            seen.add(sent)
            out.append({
                "section": name, "sentence": sent, "kind": _kind(sent),
                "age_intervals_years": [iv for iv, _ in ages], "age_matches": [t for _, t in ages],
                "weight_cutoffs": weights,
                "weight_flag_not_converted": bool(weights),
            })
    return out


def h_intervals(cautions: list[dict[str, Any]]) -> list[list[float]]:
    return merge([iv for c in cautions for iv in c["age_intervals_years"]])


# ------------------------------------------------------------------ build ------------------------------------------------------


def split_cautions(cautions: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Hazards (reduce support), weight flags (never converted to ages) and dropped sentences with the reason."""
    hazards: list[dict[str, Any]] = []
    weights: list[dict[str, Any]] = []
    dropped: list[dict[str, Any]] = []
    for c in cautions:
        if c["weight_cutoffs"]:
            weights.append({k: c[k] for k in ("section", "sentence", "weight_cutoffs")})
        if not c["age_intervals_years"]:
            continue
        reason = None
        if c["kind"] == "not_established":
            reason = "restates L_d (paediatric use not established below an age) or refers to another indication; not a hazard"
        elif c["section"] == "dosage" and DOSE_TABLE.search(c["sentence"]):
            reason = "dosing table / dose-by-body-surface text, not an age-specific hazard statement"
        elif TIME_NOT_AGE.search(c["sentence"]):
            reason = "number is a time since an event or in an adult trial (previously, post), not a patient age"
        if reason:
            dropped.append({"section": c["section"], "sentence": c["sentence"], "age_intervals_years": c["age_intervals_years"], "reason": reason})
        else:
            hazards.append({k: c[k] for k in ("section", "sentence", "kind", "age_intervals_years", "weight_cutoffs")})
    return {"hazards": hazards, "weight_flags": weights, "dropped": dropped}


def drug_record(name: str, drug_id: str, review: dict[str, Any], sections: dict[str, str], setid: str | None) -> dict[str, Any]:
    verdict = review.get("verdict", "verify")
    min_age = review.get("min_age_years")
    L = label_support(verdict, min_age)  # from the reviewed config only (indication-specific)
    parts = split_cautions(extract_cautions(sections))
    H = h_intervals(parts["hazards"])
    H_old = merge(H + [iv for d in parts["dropped"] if d["reason"].startswith("restates") for iv in d["age_intervals_years"]])
    groups = {g: support_bounds(L, H, lo, min(hi, MAX_AGE)) for g, (lo, hi) in age_groups.GROUPS.items()}
    return {
        "drug_id": drug_id, "name": name, "setid": setid, "verdict": verdict, "min_age_years": min_age, "scope": review.get("scope"),
        "L": L, "H_intervals_years": H, **parts,
        "n_hazards": len(parts["hazards"]), "n_dropped": len(parts["dropped"]), "n_weight_flags": len(parts["weight_flags"]),
        "C": support_bounds(L, H), "C_by_group": groups,
        "sensitivity_C_if_not_established_sentences_counted_as_hazards": support_bounds(L, H_old),
    }


def read_sections(setid: str) -> dict[str, str]:
    root = labels._spl_xml(setid)
    inverse = {c: k for k, codes in SECTION_CODES.items() for c in codes}
    out: dict[str, list[str]] = {k: [] for k in SECTION_CODES}
    for section in root.iter(labels.H + "section"):
        code = section.find(labels.H + "code")
        if code is not None and code.get("code") in inverse:
            t = _norm("".join(section.itertext()))
            if t and t not in out[inverse[code.get("code")]]:
                out[inverse[code.get("code")]].append(t)
    return {k: " ".join(v) for k, v in out.items()}


def build(table: dict[str, Any], config: dict[str, Any], loader=read_sections) -> dict[str, Any]:
    entries = config.get("entries", config)
    by_name = {e["name"].upper(): (k, e) for k, e in entries.items()}
    drugs = []
    for row in table["farmacos"]:
        key, review = by_name.get(row["name"].upper(), (row["drug_id"], row.get("revision_manual") or {}))
        setid = (row.get("etiqueta") or {}).get("setid")
        sections = loader(setid) if setid else {}
        drugs.append(drug_record(row["name"], row.get("drug_id", key), review, sections, setid))
    drugs.sort(key=lambda d: d["name"])
    return {"header": HEADER, "age_range_years": [0, MAX_AGE], "groups": {g: list(v) for g, v in age_groups.GROUPS.items()},
            "n_drugs": len(drugs), "n_hazards_kept": sum(d["n_hazards"] for d in drugs), "n_sentences_dropped": sum(d["n_dropped"] for d in drugs),
            "n_weight_flags": sum(d["n_weight_flags"] for d in drugs), "drugs": drugs}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", type=Path, default=Path("out/age_label_table.json"))
    p.add_argument("--config", type=Path, default=Path("config/pediatric_label_review.json"))
    p.add_argument("--out", type=Path, default=Path("out/empirical/age_support.json"))
    a = p.parse_args(argv)
    result = build(json.loads(a.table.read_text(encoding="utf-8")), json.loads(a.config.read_text(encoding="utf-8")))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=1, ensure_ascii=False, sort_keys=True), encoding="utf-8", newline="\n")
    print("drugs:", result["n_drugs"], "hazards:", result["n_hazards_kept"], "dropped:", result["n_sentences_dropped"], "weight flags:", result["n_weight_flags"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
