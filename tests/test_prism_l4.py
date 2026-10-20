from __future__ import annotations

import numpy as np

from track2 import prism_l4


def test_bh_is_monotone_and_bounded() -> None:
    adj = prism_l4.bh(np.array([0.001, 0.04, 0.03, 0.5]))
    assert adj.max() <= 1.0 and adj[0] == min(adj) and np.all(adj >= np.array([0.001, 0.04, 0.03, 0.5]) - 1e-12)


def test_match_rows_uses_synonyms_and_salt_free_names() -> None:
    keys = {"BRD:1": ("bortezomib", "BORTEZOMIB|VELCADE"), "BRD:2": ("ixazomib", "IXAZOMIB")}
    out = prism_l4.match_rows({"D1": "BORTEZOMIB", "D2": "IXAZOMIB CITRATE", "D3": "NOPE"}, keys)
    assert out == {"D1": ["BRD:1"], "D2": ["BRD:2"]}


def _data(effect: float, n: int = 200, seed: int = 1):
    rng = np.random.default_rng(seed)
    score = rng.normal(size=n)
    viability = effect * score + rng.normal(scale=0.5, size=n)  # effect < 0 means weak-score lines are less viable? see test below
    lo, hi = np.quantile(score, 0.25), np.quantile(score, 0.75)
    return viability, score, score <= lo, score >= hi, np.array(["A", "B"] * (n // 2))


def test_test_row_direction_and_adjustment() -> None:
    v, s, weak, strong, lin = _data(0.8)  # viability rises with the SAC score: weak-SAC lines are more sensitive (lower viability)
    r = prism_l4.test_row(v, s, weak, strong, lin)
    assert r["tested"] and r["diff_weak_minus_strong"] < 0 and r["p"] < 0.001 and r["p_adjusted"] < 0.001 and r["adjusted_coef_per_sd"] > 0
    flat = prism_l4.test_row(*_data(0.0, seed=3)[:1], *_data(0.0, seed=3)[1:])
    assert flat["p"] > 0.01


def test_too_few_lines_is_not_tested() -> None:
    v, s, weak, strong, lin = _data(0.5, n=20)
    assert prism_l4.test_row(v, s, weak, strong, lin) == {"tested": False}


def test_call_drug_cases() -> None:
    sens = {"fdr": 0.01, "fdr_lineage_prolif": 0.02, "diff_weak_minus_strong": -0.5, "coef_lineage_prolif_per_sd": 0.4}  # weak-SAC lines die more
    res = {"fdr": 0.01, "fdr_lineage_prolif": 0.02, "diff_weak_minus_strong": 0.5, "coef_lineage_prolif_per_sd": -0.4}
    none = {"fdr": 0.9, "fdr_lineage_prolif": 0.9, "diff_weak_minus_strong": 0.1, "coef_lineage_prolif_per_sd": -0.1}
    assert prism_l4.call_drug([sens, none]) == "supports_sensitivity_when_SAC_weak"
    assert prism_l4.call_drug([res]) == "supports_resistance_when_SAC_weak"
    assert prism_l4.call_drug([sens, res]) == "mixed" and prism_l4.call_drug([none]) == "no_association"
    # significant unadjusted but not after adjusting for lineage and proliferation
    assert prism_l4.call_drug([{**res, "fdr_lineage_prolif": 0.5}]) == "confounded_by_proliferation_or_lineage"
    # significant in both but the adjusted coefficient has the opposite sign: not confirmed
    assert prism_l4.call_drug([{**res, "coef_lineage_prolif_per_sd": 0.3}]) == "confounded_by_proliferation_or_lineage"


def test_proliferation_confounding_is_detected_by_the_adjusted_model() -> None:
    rng = np.random.default_rng(5)
    n = 400
    prolif = rng.normal(size=n)
    score = 0.9 * prolif + 0.3 * rng.normal(size=n)  # the checkpoint score is mostly proliferation
    viability = 0.8 * prolif + rng.normal(scale=0.5, size=n)  # the drug depends on proliferation only
    lo, hi = np.quantile(score, 0.25), np.quantile(score, 0.75)
    r = prism_l4.test_row(viability, score, score <= lo, score >= hi, np.array(["A", "B"] * (n // 2)), prolif)
    assert r["p"] < 1e-6 and abs(r["coef_lineage_prolif_per_sd"]) < abs(r["adjusted_coef_per_sd"])  # the effect shrinks once proliferation is in the model
    assert r["p_lineage_prolif"] > 0.01
