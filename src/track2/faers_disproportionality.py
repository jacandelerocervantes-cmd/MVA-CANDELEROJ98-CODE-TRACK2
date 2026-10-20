"""Pharmacovigilance disproportionality (PRR / ROR) from openFDA FAERS counts, per candidate drug and organ-system group.

For every candidate drug and every outcome group of `track2.toxicity.GROUPS` the 2x2 table is built from
report-level counts (one report counts once, however many terms of the group it carries):

    a = reports with the drug AND the outcome group      b = reports with the drug, without the group
    c = reports with the group, without the drug         d = reports with neither

obtained from four count queries: N_drug, a, G (group total in all of FAERS) and N (all FAERS reports):
b = N_drug - a, c = G - a, d = N - N_drug - c.

How the openFDA group definition is built (openFDA cannot evaluate a regular expression): the MedDRA preferred
terms (PT) are listed with `count=patient.reaction.reactionmeddrapt.exact` (at most 999 terms per query without an
API key) for FAERS as a whole and for each drug; the group is the union of the PTs that match the regular
expression of `toxicity.GROUPS` in ANY of those lists. The same PT list is then used for every drug, so a, c and G
share one definition. PTs outside the top 999 of both FAERS and every drug are not in the definition (the
definition is therefore conservative: it can only miss rare PTs). PTs with characters other than letters, digits,
space and hyphen (a handful, for example an apostrophe-like '^') are left out because openFDA rejects them.

Statistics: PRR = (a/(a+b)) / (c/(c+d)), ROR = (a*d)/(b*c), both with 95 % CI on the log scale; chi-square with
Yates correction. Conventional signal rule (Evans et al. 2001): PRR >= 2 AND a >= 3 AND chi2 >= 4. A second flag,
`ror_signal`, is lower bound of the ROR 95 % CI > 1 with a >= 3. When b or c is 0 the Haldane-Anscombe correction
(+0.5 in every cell) is applied and flagged; when a = 0 no estimate is given.

This is a hypothesis-generating screening of SPONTANEOUS reports: no incidence, confounding by indication,
reporting bias, duplicates, no causality. Not for medical decisions.

    python -m track2.faers_disproportionality --candidates out/candidates.json --out out
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
import urllib.parse
from collections.abc import Callable
from pathlib import Path
from typing import Any

from track2 import http_cache
from track2.toxicity import GROUPS

OPENFDA = "https://api.fda.gov/drug/event.json"
PT_FIELD = "patient.reaction.reactionmeddrapt.exact"
TOP_TERMS = 999  # maximum `limit` of a count query without an API key
PAUSE_SECONDS = 0.3  # < 240 requests per minute without a key
SIGNAL_RULE = {"prr_min": 2.0, "a_min": 3, "chi2_min": 4.0, "chi2": "Yates-corrected"}
Z95 = 1.959964
CAVEATS = [
    "Spontaneous reports: no denominator, no incidence; reporting is influenced by media, drug age and indication.",
    "Confounding by indication: the disease treated (or its complications) is reported with the drug.",
    "Duplicate and incomplete reports; a drug is matched by its generic name and may be suspect or concomitant.",
    "A disproportionality signal is a hypothesis for expert pharmacovigilance review, not evidence of causality.",
    "The outcome groups are regular-expression groups of MedDRA preferred terms (see module doc); not for medical decisions.",
]
_SALTS = r"\b(ANHYDROUS|HYDROCHLORIDE|CITRATE|MESYLATE|SULFATE|PHOSPHATE)\b"
_SAFE = ':+"()'
_QUERY_SAFE = re.compile(r"[A-Za-z0-9 \-]+")  # other characters (e.g. '^', quotes) make openFDA answer 400: such PTs are left out
Getter = Callable[[str], Any]


def clean_name(drug_name: str) -> str:
    """Same name cleaning as `toxicity.toxicity_profile`."""
    return re.sub(_SALTS, "", drug_name.upper()).split(",")[0].strip()


def drug_search(name: str) -> str:
    return f'patient.drug.openfda.generic_name.exact:"{name}"'


def _url(search: str | None, extra: str) -> str:
    prefix = f"?search={urllib.parse.quote(search, safe=_SAFE)}" if search else "?"
    return f"{OPENFDA}{prefix}{'&' if search else ''}{extra}"


def group_clause(terms: list[str]) -> str:
    return f"{PT_FIELD}:(" + "+".join(f'"{t}"' for t in terms) + ")"


def throttled_get(url: str) -> Any:
    """`http_cache.get_json` plus a pause after every request that really goes to the network."""
    cached = http_cache._cache_path(url, ".json").exists()  # noqa: SLF001 - same cache layout
    result = http_cache.get_json(url)
    if not cached:
        time.sleep(PAUSE_SECONDS)
    return result


def count_reports(search: str | None, get: Getter) -> int:
    """Number of reports matching a search (0 when openFDA answers 404 'no matches')."""
    try:
        return int(get(_url(search, "limit=1"))["meta"]["results"]["total"])
    except LookupError:
        return 0


def top_terms(search: str | None, get: Getter) -> list[str]:
    try:
        rows = get(_url(search, f"count={PT_FIELD}&limit={TOP_TERMS}"))["results"]
    except LookupError:
        return []
    return [r["term"] for r in rows]


def build_group_terms(term_lists: list[list[str]]) -> dict[str, list[str]]:
    """Union over all lists of the PTs matching each group's regular expression (sorted)."""
    seen = {t for terms in term_lists for t in terms}
    return {g: sorted(t for t in seen if re.search(p, t.lower()) and _QUERY_SAFE.fullmatch(t)) for g, p in GROUPS.items()}


def disproportionality(a: int, b: int, c: int, d: int) -> dict[str, Any]:
    """PRR, ROR (95 % CI), chi-square and the signal flags of one 2x2 table."""
    table = {"a": a, "b": b, "c": c, "d": d}
    if min(table.values()) < 0:
        raise ValueError(f"negative cell in 2x2 table: {table}")
    out: dict[str, Any] = {"table": table, "signal": False, "ror_signal": False}
    if a == 0:
        return {**out, "estimable": False, "note": "no report of the drug with this outcome group; no estimate"}
    corrected = 0 in (b, c, d)
    fa, fb, fc, fd = (x + 0.5 if corrected else float(x) for x in (a, b, c, d))
    n = fa + fb + fc + fd
    prr = (fa / (fa + fb)) / (fc / (fc + fd))
    se_prr = math.sqrt(1 / fa - 1 / (fa + fb) + 1 / fc - 1 / (fc + fd))
    ror = (fa * fd) / (fb * fc)
    se_ror = math.sqrt(1 / fa + 1 / fb + 1 / fc + 1 / fd)
    chi2 = n * (max(abs(fa * fd - fb * fc) - n / 2, 0.0)) ** 2 / ((fa + fb) * (fc + fd) * (fa + fc) * (fb + fd))
    out.update({
        "estimable": True, "haldane_correction_applied": corrected,
        "prr": round(prr, 3), "prr_ci95": [round(math.exp(math.log(prr) + s * Z95 * se_prr), 3) for s in (-1, 1)],
        "ror": round(ror, 3), "ror_ci95": [round(math.exp(math.log(ror) + s * Z95 * se_ror), 3) for s in (-1, 1)],
        "chi2_yates": round(chi2, 3),
    })
    out["signal"] = bool(prr >= SIGNAL_RULE["prr_min"] and a >= SIGNAL_RULE["a_min"] and chi2 >= SIGNAL_RULE["chi2_min"])
    out["ror_signal"] = bool(out["ror_ci95"][0] > 1 and a >= SIGNAL_RULE["a_min"])
    return out


def drug_table_inputs(name: str, group_terms: dict[str, list[str]], get: Getter) -> tuple[int, dict[str, int]]:
    """N_drug and, per group, `a` (reports with the drug and any PT of the group)."""
    search = drug_search(name)
    n_drug = count_reports(search, get)
    a = {g: count_reports(f"{search} AND {group_clause(terms)}", get) if terms else 0 for g, terms in group_terms.items()}
    return n_drug, a


def assemble(name: str, n_drug: int, a: dict[str, int], group_totals: dict[str, int], n_total: int) -> dict[str, Any]:
    groups = {}
    for g, a_g in a.items():
        b, c = n_drug - a_g, group_totals[g] - a_g
        d = n_total - n_drug - c
        groups[g] = disproportionality(a_g, b, c, d)
    return {"available": True, "name_queried": name, "n_drug_reports": n_drug, "groups": groups}


def run(candidates: list[dict[str, Any]], get: Getter = throttled_get, log: Callable[[str], None] = print) -> dict[str, Any]:
    names = {c["drug_id"]: clean_name(c["name"] or "") for c in candidates}
    n_total = count_reports(None, get)
    log(f"FAERS total reports: {n_total}")
    # pass 1: drug existence and PT lists (group definition needs every list)
    totals: dict[str, int] = {}
    lists = [top_terms(None, get)]
    for drug_id, name in names.items():
        totals[drug_id] = count_reports(drug_search(name), get) if name else 0
        if totals[drug_id]:
            lists.append(top_terms(drug_search(name), get))
        log(f"pass1 {drug_id} {name}: {totals[drug_id]} reports")
    group_terms = build_group_terms(lists)
    group_totals = {g: count_reports(group_clause(t), get) if t else 0 for g, t in group_terms.items()}
    log("group sizes: " + json.dumps({g: [len(group_terms[g]), group_totals[g]] for g in group_terms}))
    results: dict[str, Any] = {}
    for drug_id, name in names.items():
        if not totals[drug_id]:
            results[drug_id] = {"available": False, "name_queried": name, "reason": "openFDA has no reports under this generic name"}
            continue
        n_drug, a = drug_table_inputs(name, group_terms, get)
        results[drug_id] = assemble(name, n_drug, a, group_totals, n_total)
        log(f"pass2 {drug_id} {name}: " + ",".join(g for g, r in results[drug_id]["groups"].items() if r["signal"]))
    return {
        "source": "openFDA drug/event (FAERS, CC0)", "faers_total_reports": n_total, "signal_rule": SIGNAL_RULE,
        "group_definition": {"method": "union of MedDRA PTs matching the toxicity.GROUPS regex among the top 999 PTs of FAERS and of each drug",
                             "pt_count": {g: len(t) for g, t in group_terms.items()}, "group_total_reports": group_totals},
        "caveats": CAVEATS, "drugs": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=Path("out/candidates.json"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    candidates = json.loads(args.candidates.read_text(encoding="utf-8"))
    result = run(candidates, log=lambda m: print(m, flush=True))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "faers_disproportionality.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    available = [d for d in result["drugs"].values() if d["available"]]
    sys.stdout.write(f"available {len(available)}/{len(result['drugs'])}\n")


if __name__ == "__main__":
    main()
