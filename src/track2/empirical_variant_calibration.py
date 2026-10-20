"""Map computational missense-predictor scores to ACMG/AMP PP3/BP4 evidence strength
using PUBLISHED calibrated thresholds (no thresholds invented here).

Sources (both read in full text on PMC, quotes verified by code, see ``verify_quotes``):
  * Pejaver V et al. Am J Hum Genet 2022;109:2163-2177. PMID 36413997, PMC9748256.
    Table 2 (CADD v1.6 PHRED, PolyPhen-2, REVEL, SIFT; 4 pathogenic + 4 benign intervals).
  * Bergquist T et al. Genet Med 2025. PMID 40084623, PMC12208618.
    Table 1 (AlphaMissense; point-based intervals +1..+4 / -1..-4).

Table-row reading (documented, because the flattened table loses column headers):
  Pejaver Table 2 columns: Benign VeryStrong, Strong, Moderate, Supporting | Pathogenic
  Supporting, Moderate, Strong, VeryStrong ("-" = tool did not reach that strength).
  Bergquist Table 1 columns: Strong(-4), (-3), Moderate(-2), Supporting(-1), Indeterminate,
  Supporting(+1), Moderate(+2), (+3), Strong(+4).
  Bergquist intervals are inclusive and rounded to 3 decimals, leaving tiny gaps
  (e.g. 0.099 -> 0.100); gaps are bridged here (the next interval starts at its lower bound).

No patient-level data are used; inputs are protein-level predictor scores only.
Run ``python -m track2.empirical_variant_calibration`` (PYTHONPATH=src, needs network on
first run) to fetch/verify the sources and write out/empirical/variant_calibration.json.
"""
from __future__ import annotations

import hashlib
import html
import json
import re
import sys
import urllib.request
from pathlib import Path

INF = float("inf")
ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "out" / "empirical"
CACHE_DIR = OUT_DIR / "source_cache"
OUT_JSON = OUT_DIR / "variant_calibration.json"

SOURCES = {
    "pejaver2022": {
        "citation": "Pejaver V et al. Am J Hum Genet 2022;109(12):2163-2177",
        "pmid": "36413997", "pmcid": "PMC9748256",
        "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9748256/fullTextXML",
    },
    "bergquist2025": {
        "citation": "Bergquist T et al. Genet Med 2025 (ClinGen calibration of AlphaMissense, ESM1b, VARITY)",
        "pmid": "40084623", "pmcid": "PMC12208618",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12208618/",
    },
}

# Verbatim quotes (compared after whitespace normalisation, see normalise()).
QUOTES = {
    "pejaver_method": ("pejaver2022",
        "Scores that satisfied the posterior probability thresholds presented in Table 1 "
        "were then deemed to provide the corresponding strengths of evidence ( Table 2 )."),
    "pejaver_cadd": ("pejaver2022",
        "CADD – ≤0.15 (0.15, 17.3] (17.3, 22.7] [25.3, 28.1) ≥28.1 – –"),
    "pejaver_revel": ("pejaver2022",
        "REVEL ≤0.003 (0.003, 0.016] (0.016, 0.183] (0.183, 0.290] [0.644, 0.773) "
        "[0.773, 0.932) ≥0.932 –"),
    "pejaver_polyphen": ("pejaver2022",
        "PolyPhen2 – – ≤0.009 (0.009, 0.113] [0.978, 0.999) ≥0.999 – –"),
    "pejaver_sift": ("pejaver2022",
        "SIFT – – ≥0.327 [0.080, 0.327) (0, 0.001] 0 – –"),
    "bergquist_alphamissense": ("bergquist2025",
        "AlphaMissense - ≤ 0.070 [0.071, 0.099] [0.100, 0.169] [0.170, 0.791] "
        "[0.792, 0.905] [0.906, 0.971] [0.972, 0.989] ≥ 0.990"),
    "bergquist_default_threshold": ("bergquist2025",
        "In fact, the recommended thresholds for AlphaMissense (0.564) and ESM1b (−7.5) "
        "do not meet the Supporting level of evidence for pathogenicity or benignity, "
        "based on our calibration."),
}


def _iv(side, strength, points, lo, hi, lo_inc, hi_inc):
    return {"side": side, "strength": strength, "points": points,
            "lo": lo, "hi": hi, "lo_inclusive": lo_inc, "hi_inclusive": hi_inc}


def _P(strength, pts, lo, hi, li=True, hi_inc=False):
    return _iv("PP3", strength, pts, lo, hi, li, hi_inc)


def _B(strength, pts, lo, hi, li=False, hi_inc=True):
    return _iv("BP4", strength, pts, lo, hi, li, hi_inc)


# Pejaver strengths: VeryStrong=8 points, Strong=4, Moderate=2, Supporting=1.
THRESHOLDS = {
    "CADD": {"source": "pejaver2022", "quote_id": "pejaver_cadd", "version_note": "CADD v1.6 PHRED-scaled",
             "intervals": [
                 _B("Strong", 4, -INF, 0.15, True, True),
                 _B("Moderate", 2, 0.15, 17.3), _B("Supporting", 1, 17.3, 22.7),
                 _P("Supporting", 1, 25.3, 28.1), _P("Moderate", 2, 28.1, INF)]},
    "REVEL": {"source": "pejaver2022", "quote_id": "pejaver_revel", "version_note": "",
              "intervals": [
                  _B("VeryStrong", 8, -INF, 0.003, True, True),
                  _B("Strong", 4, 0.003, 0.016), _B("Moderate", 2, 0.016, 0.183),
                  _B("Supporting", 1, 0.183, 0.290),
                  _P("Supporting", 1, 0.644, 0.773), _P("Moderate", 2, 0.773, 0.932),
                  _P("Strong", 4, 0.932, INF)]},
    "PolyPhen2": {"source": "pejaver2022", "quote_id": "pejaver_polyphen", "version_note": "",
                  "intervals": [
                      _B("Moderate", 2, -INF, 0.009, True, True), _B("Supporting", 1, 0.009, 0.113),
                      _P("Supporting", 1, 0.978, 0.999), _P("Moderate", 2, 0.999, INF)]},
    "SIFT": {"source": "pejaver2022", "quote_id": "pejaver_sift", "version_note": "",
             "intervals": [
                 _B("Moderate", 2, 0.327, INF, True, True), _B("Supporting", 1, 0.080, 0.327, True, False),
                 _P("Supporting", 1, 0.0, 0.001, False, True), _P("Moderate", 2, 0.0, 0.0, True, True)]},
    "AlphaMissense": {"source": "bergquist2025", "quote_id": "bergquist_alphamissense",
                      "version_note": "gaps between rounded inclusive intervals bridged",
                      "intervals": [
                          _B("Moderate(3pt)", 3, -INF, 0.070, True, True),
                          _B("Moderate", 2, 0.070, 0.100, False, False),
                          _B("Supporting", 1, 0.100, 0.170, True, False),
                          _P("Supporting", 1, 0.792, 0.906), _P("Moderate", 2, 0.906, 0.972),
                          _P("Moderate(3pt)", 3, 0.972, 0.990), _P("Strong", 4, 0.990, INF, True, True)]},
}
for _t in THRESHOLDS.values():
    _s = SOURCES[_t["source"]]
    _t["pmid"], _t["citation"] = _s["pmid"], _s["citation"]
    _t["quote"] = QUOTES[_t["quote_id"]][1]


def _inside(iv, x):
    lo_ok = x >= iv["lo"] if iv["lo_inclusive"] else x > iv["lo"]
    hi_ok = x <= iv["hi"] if iv["hi_inclusive"] else x < iv["hi"]
    return lo_ok and hi_ok


def strength(predictor: str, score: float) -> dict:
    """Return the calibrated PP3/BP4 strength for one score (or 'none' if indeterminate)."""
    if predictor not in THRESHOLDS:
        raise KeyError(f"no published calibration for {predictor!r}")
    t = THRESHOLDS[predictor]
    for iv in t["intervals"]:
        if _inside(iv, score):
            rng = ("[" if iv["lo_inclusive"] else "(") + f"{iv['lo']}, {iv['hi']}" + ("]" if iv["hi_inclusive"] else ")")
            return {"predictor": predictor, "score": score, "criterion": iv["side"],
                    "strength": f"{iv['side']}_{iv['strength']}", "points": iv["points"],
                    "threshold_crossed": rng, "pmid": t["pmid"]}
    return {"predictor": predictor, "score": score, "criterion": None, "strength": "none",
            "points": 0, "threshold_crossed": None, "pmid": t["pmid"]}


# ---------------------------------------------------------------- source verification
class QuoteNotFound(RuntimeError):
    pass


def normalise(text: str) -> str:
    text = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", text, flags=re.S)
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def fetch_source(key: str) -> str:
    req = urllib.request.Request(SOURCES[key]["url"], headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return normalise(r.read().decode("utf-8"))


def verify_quotes(texts: dict) -> dict:
    """Raise QuoteNotFound unless every quote is a substring of its normalised source."""
    missing = [qid for qid, (src, q) in QUOTES.items() if normalise(q) not in texts[src]]
    if missing:
        raise QuoteNotFound(f"quotes not found in source text: {missing}")
    return {src: hashlib.sha256(t.encode("utf-8")).hexdigest() for src, t in texts.items()}


def load_cached_sources() -> dict:
    return {k: (CACHE_DIR / f"{k}.txt").read_text(encoding="utf-8") for k in SOURCES}


# ---------------------------------------------------------------- application
FACTS = {  # protein-level scores from the Track 1 report; the only inputs
    "variant": "BUB1B p.Asn1002Lys (MANE NM_001211.6)",
    "scores": {"CADD": 24.5, "SIFT": 0.01, "PolyPhen2": 0.997, "AlphaMissense": 0.923, "REVEL": 0.472},
    "gnomad_v4_exome_af": 6.8e-07,
}


def build_report(source_hashes: dict) -> dict:
    per = {}
    for p, s in FACTS["scores"].items():
        r = strength(p, s)
        r["threshold_source_quote"] = THRESHOLDS[p]["quote"]
        r["source"] = THRESHOLDS[p]["citation"]
        per[p] = r
    return {
        "variant": FACTS["variant"],
        "inputs": {"scores": FACTS["scores"], "gnomad_v4_exome_af": FACTS["gnomad_v4_exome_af"]},
        "per_predictor": per,
        "overall_statement": (
            "Calibrated strengths disagree across predictors: AlphaMissense 0.923 reaches PP3_Moderate "
            "(Bergquist 2025) and PolyPhen-2 0.997 reaches PP3_Supporting (Pejaver 2022), whereas REVEL 0.472, "
            "CADD-PHRED 24.5 and SIFT 0.01 fall in indeterminate gaps and give no PP3 strength. "
            "REVEL, the ClinGen-recommended tool of Pejaver 2022, is indeterminate. Choosing the highest-scoring "
            "tool after seeing all scores is post hoc; the defensible statement is 'PP3 at most Moderate "
            "(AlphaMissense only), none by REVEL', and the tools are correlated, so strengths must not be summed. "
            "The default AlphaMissense cut-off (0.564) is not itself calibrated evidence (Bergquist 2025)."),
        "acmg_combination_caveat": (
            "Even the most favourable reading (PM2 for rarity at Supporting, as downgraded by ClinGen SVI, plus "
            "PP3_Moderate) does not classify the variant: it remains a variant of uncertain significance "
            "(Likely pathogenic needs more points than this; point system from Tavtigian 2020, not re-verified "
            "here). PM2 and PP3 are weak lines of evidence. Phase (cis/trans relative to a second BUB1B variant) "
            "is unresolved, so no inference on biallelic status or disease can be drawn."),
        "acmg_classification": "Variant of uncertain significance (unchanged); phase unresolved",
        "sources": {k: {**{x: v[x] for x in ("citation", "pmid", "pmcid", "url")},
                        "normalised_text_sha256": source_hashes[k]} for k, v in SOURCES.items()},
        "quotes_verified_by_code": True,
        "quotes": {q: {"source": s, "text": t} for q, (s, t) in QUOTES.items()},
    }


def main() -> int:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    texts = {k: fetch_source(k) for k in SOURCES}
    hashes = verify_quotes(texts)  # raises QuoteNotFound -> non-zero exit, nothing written
    for k, t in texts.items():
        (CACHE_DIR / f"{k}.txt").write_text(t, encoding="utf-8")
    OUT_JSON.write_text(json.dumps(build_report(hashes), indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
