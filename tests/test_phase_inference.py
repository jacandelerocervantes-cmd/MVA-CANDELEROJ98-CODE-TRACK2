"""Offline tests for phase_inference (no network, no patient data)."""

import json
import math

import pytest

from track2 import phase_inference as p


def test_prior_is_one_half_for_any_de_novo_mix():
    for pi in (0, 0.1, 1):
        for q in (0, 0.5, 0.8, 1):
            assert p.prior_trans(pi, q) == pytest.approx(0.5)


def test_posterior_closed_form_and_monotonic():
    assert p.posterior_trans(0.9, 1e-6) == pytest.approx(0.9 / (0.9 + 1e-6))
    assert p.posterior_trans(1.0, 1e-3) > p.posterior_trans(0.5, 1e-3)
    assert p.posterior_trans(0.9, 1e-6) > p.posterior_trans(0.9, 1e-3)
    assert p.posterior_trans(0.9, 1e-6, c=0.0) == pytest.approx(0.5)  # BUB1B not causal: phase uninformative


def test_posterior_falls_with_p_causal():
    vals = [r["posterior_trans"] for r in p.causal_table(0.9, 1e-6)]
    assert vals == sorted(vals, reverse=True) or vals == sorted(vals)
    assert vals[0] > 0.999 and vals[-1] == pytest.approx(0.5)
    t = p.causal_threshold(0.9, 0.9, 1e-6)
    assert 0 < t < 1
    assert p.posterior_trans(0.9, 1e-6, t) == pytest.approx(0.9, abs=1e-6)
    assert p.causal_threshold(0.99999999, 0.5, 1e-3) is None


def test_sole_causal_check():
    ok = {"genes": {"BUB1B": {"verdict": "causal candidate"}, "X": {"verdict": "contextual"}}}
    bad = {"genes": {"BUB1B": {"verdict": "causal candidate"}, "X": {"verdict": "causal candidate"}}}
    assert p.bub1b_sole_causal(ok) and not p.bub1b_sole_causal(bad)


def test_pattern_likelihoods_sum_to_one_and_posteriors():
    for pi in (0, 0.05):
        pl = p.pattern_likelihoods(pi, 0.5, 0.01)
        for h in ("trans", "cis"):
            assert sum(pl[h].values()) == pytest.approx(1.0)
    assert p.posterior_after_parental_test("T") == 1.0
    assert p.posterior_after_parental_test("C") == 0.0
    assert p.posterior_after_parental_test("DN_M", 0.5, 0.05, 0.5) == pytest.approx(0.5)
    assert p.posterior_after_parental_test("DN_M", 0.5, 0.05, 1.0) == 1.0
    assert p.posterior_after_parental_test("T", 0.5, 0.0, 0.5, 0.01) == pytest.approx(0.99)


def test_long_read_formulas():
    assert p.p_span_given_cover(10_000, 2_500) == pytest.approx(0.75)
    assert p.p_span_given_cover(10_000, 10_000) == 0.0
    assert p.p_span_given_cover(10_000, 60_000) == 0.0
    lam = 20 * 0.75
    assert p.p_phase_decisive(20, 10_000, 2_500, 1) == pytest.approx(1 - math.exp(-lam))
    assert p.p_phase_decisive(20, 10_000, 2_500, 3) == pytest.approx(1 - math.exp(-lam) * (1 + lam + lam**2 / 2))
    d = p.min_depth(10_000, 2_500, 1, 0.95)
    assert d == pytest.approx(-math.log(0.05) / 0.75, rel=1e-6)
    assert p.min_depth(25_000, p.BUB1B_GENE_LENGTH_BP) is None
    assert p.min_depth(25_000, 5_000, 3) > p.min_depth(25_000, 5_000, 1)
    assert p.min_depth(10_000, 8_000) > p.min_depth(25_000, 8_000)


def test_build_is_json_serialisable_and_has_no_patient_fields():
    diff = {"genes": {"BUB1B": {"verdict": "causal candidate"}}}
    doc = p.build(diff)
    s = json.dumps(doc, sort_keys=True)
    lo, hi = doc["headline"]["range_over_penetrance_and_alt_cause_grid"]
    assert 0.99 < lo <= doc["headline"]["posterior_trans"] <= hi < 1
    assert doc["long_read_phasing"]["gene_body_length_bp"] == 60154
    for forbidden in ("chr15", "GT", "T>G", "QUAL"):
        assert forbidden not in s
