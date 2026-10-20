"""Planning power and sample-size table for the bench experiments proposed in report section S4.

The report proposes, for each drug hypothesis, an experiment that would refute it: aneuploidy by chromosome
spreads or single-cell sequencing (plus senescence, phospho-S6 and myotube readouts) in BUB1B-hypomorphic human
fibroblasts or myoblasts, drug against vehicle. This module sizes the aneuploidy readout only, because it is the
only readout for which a published baseline in MVA patient cells could be sourced and checked.

THESE ARE PLANNING NUMBERS, NOT RESULTS. No experiment was run and no patient data are used.

Method
------
* Outcome: fraction of aneuploid metaphases (or cells) per arm. Null: equal fractions. Alternative: the drug arm
  has the baseline fraction p0 reduced by a relative 25, 50 or 75 percent (p1 = p0 * (1 - r)).
* Test: two-sided, alpha 0.05, two independent proportions. Primary size: arcsine (Cohen's h) formula,
  n = 2 * ((z_{1-alpha/2} + z_{power}) / h)^2 per arm, with h = |2 asin(sqrt(p0)) - 2 asin(sqrt(p1))|.
  Check: the exact conditional (Fisher) test power is computed by enumeration of the two binomials, at the
  arcsine size and, when the size is small enough, as the smallest n whose exact power reaches the target.
* Clustering: cells within one culture are correlated. Design effect DE = 1 + (m - 1) * ICC for m cells scored
  per independent culture; independent cultures per arm k = max(3, ceil(n * DE / m)); cells per arm = k * m.
  ICC is not sourced (no between-culture variance was found); it is a stated assumption reported as the range
  0.05 to 0.3.

Command: ``python -m track2.empirical_power`` (needs network once: re-fetches and verifies every quote, and fails
loudly if any quote is not found in its source text).
"""

from __future__ import annotations

import html
import json
import math
import re
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import binom, hypergeom, norm

OUT_PATH = Path("out/empirical/power_analysis.json")
EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"

ALPHA = 0.05
POWERS = (0.8, 0.9)
RELATIVE_REDUCTIONS = (0.25, 0.50, 0.75)
ICC_VALUES = (0.05, 0.1, 0.2, 0.3)
CELLS_PER_CULTURE = (50, 100, 200)
HEADLINE_CELLS_PER_CULTURE = 50
MIN_CULTURES_PER_ARM = 3
EXACT_SCAN_MAX_N = 400  # smallest exact n is searched only when the arcsine size is at most this
EXACT_POWER_MAX_N = 2500  # exact power is evaluated at the arcsine size only up to this n

# Verbatim quotes. "pmcid" None means the quote comes from the abstract record of the PMID.
SOURCES: list[dict[str, Any]] = [
    {"id": "BUB1B_MVA1_patient", "pmid": "35804254", "pmcid": "PMC9616775",
     "cell_type": "cultured patient cells (karyotype of one BUB1B MVA1 patient)",
     "quote": "Cytogenetic analysis revealed a normal karyotype in 88 cells (46, XX) and various aneuploidies in 12 cells.",
     "value": {"aneuploid": 12, "total": 100}},
    {"id": "MAD1L1_PBMC_metaphases", "pmid": "36322655", "pmcid": "PMC9629740",
     "cell_type": "peripheral blood mononuclear cells, metaphase spreads (biallelic MAD1L1, not BUB1B)",
     "quote": "identified 39% of aneuploid metaphases in the proband",
     "value": {"fraction": 0.39}},
    {"id": "MAD1L1_single_cell_seq", "pmid": "36322655", "pmcid": "PMC9629740",
     "cell_type": "peripheral blood cells, shallow single-cell sequencing (biallelic MAD1L1, not BUB1B)",
     "quote": "identified whole chromosome aneuploidies in 10 of 23 (43.5%) proband cells but not in controls",
     "value": {"aneuploid": 10, "total": 23}},
    {"id": "SMC5_MVA_lymphocytes", "pmid": "41374403", "pmcid": "PMC12691382",
     "cell_type": "lymphocyte karyotype (biallelic SMC5, not BUB1B)",
     "quote": "Karyotyping of 71 metaphase cells revealed abnormal ploidy in 36 cells (50.7%), specifically, 9 cells (12.7%) "
              "showed near-tetraploidy and 27 cells (38.0%) exhibited MVA",
     "value": {"aneuploid": 27, "total": 71}},
    {"id": "CEP192_MVA_lymphocytes", "pmid": "37981762", "pmcid": "PMC10716027",
     "cell_type": "lymphocytes, GTG-banding, two siblings (biallelic CEP192, not BUB1B)",
     "quote": "GTG-banding in III-1 and III-2 revealed that 15.52%–24.14% of lymphocytes were tetraploidy cells, "
              "26.44%–27.27% were MVA cells",
     "value": {"fraction_range": [0.2644, 0.2727]}},
    {"id": "MAD2L1BP_low_level", "pmid": "37796616", "pmcid": "PMC10721328",
     "cell_type": "metaphases from archived slides (biallelic MAD2L1BP, not BUB1B); low-level mosaic example",
     "quote": "Reanalysis of 200 metaphases from 3 preserved slides archived 8 years ago revealed moderate mosaic aneuploidy, "
              "with 10% of cells showing chromosome (mostly) losses and gains",
     "value": {"aneuploid": 20, "total": 200}},
    {"id": "MVA_definition_threshold", "pmid": "37796616", "pmcid": "PMC10721328",
     "cell_type": "definition of MVA",
     "quote": "constitutional mosaicism for aneuploidy (chromosome gains and losses of variable proportion, but usually "
              "affecting ≥25% of cells)",
     "value": {"fraction_lower_bound": 0.25}},
    {"id": "MVA_lymphocytes_about_half", "pmid": "18548531", "pmcid": None,
     "cell_type": "cultured lymphocytes of two MVA sibs without BUB1B mutation",
     "quote": "Around one half of the cultured lymphocytes from our patients had aneuploidies ranging from nullisomies to heptasomies",
     "value": {"fraction": 0.5}},
    {"id": "tetraploid_background_children", "pmid": "41374403", "pmcid": "PMC12691382",
     "cell_type": "control children, tetraploid cells only (not aneuploidy); background example",
     "quote": "we analyzed 100 metaphase cells from each of 10 children with congenital structural abnormalities. Tetraploid cells "
              "were rare—observed in only one cell from a single individual (0.1% per sample, on average)",
     "value": {"tetraploid_fraction": 0.001}},
]

# Baseline scenarios used in the table: (key, p0, source ids, role).
SCENARIOS: list[dict[str, Any]] = [
    {"key": "S_low_BUB1B", "p0": 0.12, "sources": ["BUB1B_MVA1_patient"],
     "role": "only BUB1B-specific patient value found (12 of 100 cells); conservative (small effect sizes need more cells)"},
    {"key": "S_threshold", "p0": 0.25, "sources": ["MVA_definition_threshold"],
     "role": "lower bound of the usual MVA range (>=25 percent of cells); non-BUB1B and BUB1B patients combined"},
    {"key": "S_high", "p0": 0.39, "sources": ["MAD1L1_PBMC_metaphases", "SMC5_MVA_lymphocytes"],
     "role": "upper end of published MVA lymphocyte values (39 percent, 38 percent); non-BUB1B genes"},
]

NOT_SOURCED = [
    "Aneuploid fraction in euploid control human fibroblasts or lymphocytes (background rate): no numeric value was "
    "found in the sources searched; the MAD1L1 paper only states aneuploidy was 'not in controls'. The drug-vs-vehicle "
    "design does not need it, but a patient-vs-control comparison would.",
    "Between-culture variability and intraclass correlation of aneuploid fractions in MVA/BUB1B cells: not found. ICC "
    "0.05-0.3 is an assumption.",
    "Aneuploid fraction in BUB1B-hypomorphic fibroblasts or myoblasts specifically: Hanks 2004 (PMID 15475955), Matsuura 2006 "
    "(PMID 16411201) and Suijkerbuijk 2010 (PMID 20516114) abstracts, or full texts when available, gave no usable fraction.",
    "Baselines and variances for senescence (SA-beta-gal, p21), phospho-S6 and myotube differentiation readouts: not "
    "sourced, not sized here. Size them from pilot data (continuous outcomes need a standardised effect size and a "
    "between-culture variance).",
]


# ----------------------------------------------------------------------------------------------- quotes

def normalise(text: str) -> str:
    """HTML-unescape, collapse whitespace; used on both the quote and the source text."""
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def verify_quote(quote: str, source_text: str) -> bool:
    return normalise(quote) in normalise(source_text)


def assert_quote(quote: str, source_text: str, label: str) -> None:
    if not verify_quote(quote, source_text):
        raise AssertionError(f"Quote not found in source text ({label}): {quote!r}")


def _strip_tags(xml: str) -> str:
    return re.sub(r"<[^>]+>", " ", xml)


def fetch_source_text(source: dict[str, Any]) -> str:
    """Abstract (PMID record) or open-access full text (PMC) from Europe PMC; errors raise, never pass silently."""
    from track2.http_cache import get_bytes  # lazy: keeps the module importable without network code paths

    time.sleep(0.25)  # <= 4 requests/s
    if source["pmcid"]:
        url = f"{EUROPE_PMC}/{source['pmcid']}/fullTextXML"
        return _strip_tags(get_bytes(url, accept="application/xml").decode("utf-8"))
    url = f"{EUROPE_PMC}/search?query=EXT_ID:{source['pmid']}%20AND%20SRC:MED&format=json&resultType=core"
    result = json.loads(get_bytes(url))["resultList"]["result"]
    if not result:
        raise LookupError(f"PMID {source['pmid']} not found in Europe PMC")
    return _strip_tags((result[0].get("title") or "") + " " + (result[0].get("abstractText") or ""))


def verify_all_sources(texts: dict[str, str] | None = None) -> list[dict[str, Any]]:
    """Check every quote; `texts` (keyed by source id) allows offline use. Raises on the first miss."""
    out = []
    for src in SOURCES:
        text = texts[src["id"]] if texts is not None else fetch_source_text(src)
        assert_quote(src["quote"], text, f"{src['id']} PMID {src['pmid']}")
        out.append({**src, "verified": True, "verified_against": "PMC full text" if src["pmcid"] else "PubMed abstract (Europe PMC)"})
    return out


# ----------------------------------------------------------------------------------------------- statistics

def cohens_h(p1: float, p2: float) -> float:
    return abs(2 * math.asin(math.sqrt(p1)) - 2 * math.asin(math.sqrt(p2)))


def arcsine_n_per_arm(p1: float, p2: float, alpha: float = ALPHA, power: float = 0.8) -> int:
    """Two-sided arcsine (Cohen's h) sample size per arm: n = 2 * ((z_a + z_b) / h)^2, rounded up."""
    if not (0 < p1 < 1 and 0 < p2 < 1) or p1 == p2:
        raise ValueError("proportions must differ and lie strictly inside (0, 1)")
    h = cohens_h(p1, p2)
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    return math.ceil(2 * (z / h) ** 2)


def arcsine_power(p1: float, p2: float, n: int, alpha: float = ALPHA) -> float:
    h = cohens_h(p1, p2)
    return float(norm.cdf(h * math.sqrt(n / 2) - norm.ppf(1 - alpha / 2)))


def _fisher_reject_matrix(n: int, alpha: float) -> np.ndarray:
    """reject[x1, x2] for the two-sided Fisher exact test with n cells per arm (minimum-likelihood p-value)."""
    reject = np.zeros((n + 1, n + 1), dtype=bool)
    for s in range(2 * n + 1):
        lo, hi = max(0, s - n), min(n, s)
        xs = np.arange(lo, hi + 1)
        pmf = hypergeom.pmf(xs, 2 * n, n, s)
        # two-sided p of x = sum of pmf(y) <= pmf(x) (relative tolerance as in scipy.stats.fisher_exact)
        order = np.argsort(pmf)
        cum = np.cumsum(pmf[order])
        pvals = np.empty_like(pmf)
        sorted_pmf = pmf[order]
        # ties: every x with equal pmf gets the cumulative sum up to the last tie
        last = np.searchsorted(sorted_pmf, sorted_pmf * (1 + 1e-7), side="right") - 1
        pvals[order] = cum[last]
        reject[xs, s - xs] = pvals <= alpha
    return reject


def exact_power(p1: float, p2: float, n: int, alpha: float = ALPHA) -> float:
    """Exact power of the two-sided Fisher exact test for n cells per arm, by full enumeration."""
    xs = np.arange(n + 1)
    joint = np.outer(binom.pmf(xs, n, p1), binom.pmf(xs, n, p2))
    return float(joint[_fisher_reject_matrix(n, alpha)].sum())


def exact_min_n(p1: float, p2: float, target: float, alpha: float = ALPHA, n_start: int = 2, n_stop: int | None = None) -> int | None:
    """Smallest n with exact power >= target (first crossing; exact power is saw-toothed), or None if not reached by n_stop."""
    n_stop = n_stop or 10_000
    for n in range(n_start, n_stop + 1):
        if exact_power(p1, p2, n, alpha) >= target:
            return n
    return None


def design_effect(cells_per_culture: int, icc: float) -> float:
    return 1 + (cells_per_culture - 1) * icc


def cultures_and_cells(n_independent: int, cells_per_culture: int, icc: float, min_cultures: int = MIN_CULTURES_PER_ARM) -> dict[str, Any]:
    """Cultures k and cells per arm needed when cells within a culture share an intraclass correlation."""
    de = design_effect(cells_per_culture, icc)
    k = max(min_cultures, math.ceil(n_independent * de / cells_per_culture))
    return {"design_effect": round(de, 4), "cultures_per_arm": k, "cells_per_culture": cells_per_culture,
            "cells_per_arm": k * cells_per_culture}


# ----------------------------------------------------------------------------------------------- table

def build_table(exact: bool = True) -> list[dict[str, Any]]:
    rows = []
    for sc in SCENARIOS:
        for r in RELATIVE_REDUCTIONS:
            p0 = sc["p0"]
            p1 = round(p0 * (1 - r), 6)
            for power in POWERS:
                n = arcsine_n_per_arm(p0, p1, ALPHA, power)
                row: dict[str, Any] = {
                    "scenario": sc["key"], "baseline_fraction": p0, "relative_reduction": r, "drug_arm_fraction": p1,
                    "alpha_two_sided": ALPHA, "power_target": power, "cohens_h": round(cohens_h(p0, p1), 5),
                    "n_cells_per_arm_independent": n,
                    "arcsine_power_at_n": round(arcsine_power(p0, p1, n), 4),
                }
                if exact:
                    row["exact_fisher_power_at_arcsine_n"] = round(exact_power(p0, p1, n), 4) if n <= EXACT_POWER_MAX_N else None
                    row["exact_fisher_min_n"] = (exact_min_n(p0, p1, power, n_start=max(2, int(n * 0.7)), n_stop=int(n * 1.5))
                                                 if n <= EXACT_SCAN_MAX_N else None)
                clustered = []
                for m in CELLS_PER_CULTURE:
                    for icc in ICC_VALUES:
                        clustered.append({"icc": icc, **cultures_and_cells(n, m, icc)})
                row["clustered"] = clustered
                rows.append(row)
    return rows


def markdown_summary(rows: list[dict[str, Any]]) -> str:
    m = HEADLINE_CELLS_PER_CULTURE
    lo, hi = min(ICC_VALUES), max(ICC_VALUES)
    head = (f"Planning numbers, not results. Two-sided alpha {ALPHA}, arcsine sample size per arm; clustering with {m} cells per culture, "
            f"ICC {lo}-{hi} (assumed). Independent cells = no clustering. Cultures = independent cultures per arm (minimum {MIN_CULTURES_PER_ARM}).")
    lines = [head, "",
             "| Baseline aneuploid fraction | Relative reduction (drug arm) | Power | Independent cells per arm | Cultures per arm (ICC 0.05 / 0.3) | Cells per arm (ICC 0.05 / 0.3) |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        c = {x["icc"]: x for x in r["clustered"] if x["cells_per_culture"] == m}
        lines.append(
            f"| {r['baseline_fraction']:.2f} | {int(r['relative_reduction'] * 100)}% ({r['drug_arm_fraction']:.4f}) | {r['power_target']} | "
            f"{r['n_cells_per_arm_independent']} | {c[lo]['cultures_per_arm']} / {c[hi]['cultures_per_arm']} | "
            f"{c[lo]['cells_per_arm']} / {c[hi]['cells_per_arm']} |")
    return "\n".join(lines)


def build_payload(verified_sources: list[dict[str, Any]], exact: bool = True) -> dict[str, Any]:
    rows = build_table(exact=exact)
    return {
        "status": "PLANNING NUMBERS, NOT RESULTS. No experiment was run; these figures size a future refutation experiment from "
                  "published patient-cell baselines and an assumed intraclass correlation.",
        "scope": "Aneuploidy readout (chromosome spreads or single-cell sequencing) only; drug vs vehicle in BUB1B-hypomorphic "
                 "human fibroblasts or myoblasts. Report section S4 'Experiment that would refute it'.",
        "methods": {
            "test": "two independent proportions, two-sided, alpha 0.05",
            "primary_size": "arcsine (Cohen's h): n = 2 * ((z_{1-alpha/2} + z_power) / h)^2 per arm",
            "exact_check": "two-sided Fisher exact test power by full enumeration of the two binomials (exact_fisher_power_at_arcsine_n; "
                           "exact_fisher_min_n = first n reaching the target, only computed when the arcsine n <= %d)" % EXACT_SCAN_MAX_N,
            "effect": "relative reduction r of the baseline fraction p0: p1 = p0 * (1 - r), r in 25, 50, 75 percent",
            "clustering": "design effect DE = 1 + (m - 1) * ICC; cultures per arm k = max(%d, ceil(n * DE / m)); cells per arm = k * m"
                          % MIN_CULTURES_PER_ARM,
            "icc_range_assumed": [min(ICC_VALUES), max(ICC_VALUES)],
            "icc_values": list(ICC_VALUES), "cells_per_culture_values": list(CELLS_PER_CULTURE),
            "powers": list(POWERS),
        },
        "inputs": {
            "alpha": ALPHA, "powers": list(POWERS), "relative_reductions": list(RELATIVE_REDUCTIONS),
            "baseline_scenarios": SCENARIOS,
            "assumptions": ["ICC 0.05-0.3 (not sourced)", "the same baseline fraction applies to BUB1B-hypomorphic fibroblasts/myoblasts as to "
                            "the lymphocyte/PBMC patient values quoted (not established)",
                            "no loss of cells or cultures; no multiplicity correction for several readouts or several drugs"],
        },
        "sources": verified_sources,
        "not_sourced": NOT_SOURCED,
        "caveats": [
            "The arcsine size is approximate and optimistic for these small proportions: the exact Fisher test reaches only "
            "about 70-88 percent power at the arcsine n for a 0.8-0.9 target (see exact_fisher_power_at_arcsine_n); where "
            "exact_fisher_min_n is given it is about 5-20 percent larger. Use the larger of the two when planning.",
            "Baselines come from patient lymphocytes/PBMCs (mostly non-BUB1B MVA genes) and single cases; fibroblast or myoblast "
            "cultures may differ, and cultured cells may lose aneuploid clones.",
            "Each scenario is a published point estimate from one patient or family, not a population mean; variability between "
            "patients is not estimable from these sources.",
            "Relative reductions are the effects to detect; whether a drug could plausibly lower aneuploidy that much is not "
            "assessed. For the report's refutation logic the failure case is a rise in aneuploidy, which this two-sided test also detects.",
        ],
        "table": rows,
        "markdown_summary": markdown_summary(rows),
    }


def write_payload(payload: dict[str, Any], path: Path = OUT_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


def main() -> None:
    verified = verify_all_sources()  # raises if any quote is missing
    write_payload(build_payload(verified))
    print(f"Wrote {OUT_PATH} ({len(verified)} quotes verified)")


if __name__ == "__main__":
    main()
