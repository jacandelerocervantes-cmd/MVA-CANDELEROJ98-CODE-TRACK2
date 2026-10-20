"""Offline tests for track2.empirical_variant_calibration (uses the cached source text
written by the module's main(); no network)."""
import json
from pathlib import Path

import pytest

from track2 import empirical_variant_calibration as m

pytestmark = pytest.mark.skipif(not Path("out/empirical/source_cache").exists(), reason="the cached source texts are not part of the repository (copyright); run the module to rebuild them")


def test_published_quotes_found_in_cached_sources():
    texts = m.load_cached_sources()  # fails loudly if the cache is missing
    m.verify_quotes(texts)  # raises QuoteNotFound otherwise


def test_missing_quote_fails_loudly():
    texts = {k: "unrelated text" for k in m.SOURCES}
    with pytest.raises(m.QuoteNotFound):
        m.verify_quotes(texts)


@pytest.mark.parametrize("pred,score,expected", [
    ("REVEL", 0.472, "none"), ("REVEL", 0.644, "PP3_Supporting"), ("REVEL", 0.773, "PP3_Moderate"),
    ("REVEL", 0.932, "PP3_Strong"), ("REVEL", 0.003, "BP4_VeryStrong"), ("REVEL", 0.2, "BP4_Supporting"),
    ("CADD", 24.5, "none"), ("CADD", 25.3, "PP3_Supporting"), ("CADD", 28.1, "PP3_Moderate"),
    ("CADD", 20, "BP4_Supporting"),
    ("SIFT", 0.01, "none"), ("SIFT", 0.0, "PP3_Moderate"), ("SIFT", 0.001, "PP3_Supporting"),
    ("PolyPhen2", 0.997, "PP3_Supporting"), ("PolyPhen2", 0.999, "PP3_Moderate"),
    ("AlphaMissense", 0.923, "PP3_Moderate"), ("AlphaMissense", 0.564, "none"),
    ("AlphaMissense", 0.792, "PP3_Supporting"), ("AlphaMissense", 0.99, "PP3_Strong"),
    ("AlphaMissense", 0.05, "BP4_Moderate(3pt)"),
])
def test_strength_mapping(pred, score, expected):
    assert m.strength(pred, score)["strength"] == expected


def test_unknown_predictor_raises():
    with pytest.raises(KeyError):
        m.strength("Foo", 1.0)


def test_report_for_variant_facts():
    rep = m.build_report({k: "x" for k in m.SOURCES})
    s = {p: r["strength"] for p, r in rep["per_predictor"].items()}
    assert s == {"CADD": "none", "SIFT": "none", "PolyPhen2": "PP3_Supporting",
                 "AlphaMissense": "PP3_Moderate", "REVEL": "none"}
    assert "unresolved" in rep["acmg_combination_caveat"]
    assert "uncertain significance" in rep["acmg_classification"]


def test_written_json_consistent():
    data = json.loads(m.OUT_JSON.read_text(encoding="utf-8"))
    assert data["quotes_verified_by_code"] is True
    assert data["per_predictor"]["AlphaMissense"]["strength"] == "PP3_Moderate"
