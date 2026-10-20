from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from track2 import empirical_power as ep


def test_arcsine_matches_textbook_value() -> None:
    # 0.5 vs 0.4, alpha 0.05, power 0.8: about 387 per arm
    assert abs(ep.arcsine_n_per_arm(0.5, 0.4, 0.05, 0.8) - 387) <= 2
    # symmetric in the two proportions
    assert ep.arcsine_n_per_arm(0.5, 0.4) == ep.arcsine_n_per_arm(0.4, 0.5)


def test_arcsine_monotone_in_power_and_effect() -> None:
    assert ep.arcsine_n_per_arm(0.12, 0.06, power=0.9) > ep.arcsine_n_per_arm(0.12, 0.06, power=0.8)
    assert ep.arcsine_n_per_arm(0.12, 0.09) > ep.arcsine_n_per_arm(0.12, 0.06) > ep.arcsine_n_per_arm(0.12, 0.03)
    n = ep.arcsine_n_per_arm(0.25, 0.125, power=0.8)
    assert ep.arcsine_power(0.25, 0.125, n) >= 0.8 > ep.arcsine_power(0.25, 0.125, n - 1)


def test_invalid_proportions_fail_loudly() -> None:
    with pytest.raises(ValueError):
        ep.arcsine_n_per_arm(0.3, 0.3)
    with pytest.raises(ValueError):
        ep.arcsine_n_per_arm(0.0, 0.2)


def test_exact_power_is_valid_and_below_arcsine_for_small_proportions() -> None:
    p1, p2, n = 0.4, 0.1, 30
    power = ep.exact_power(p1, p2, n)
    assert 0.5 < power < 1.0
    # type I error of the exact test never exceeds alpha
    assert ep.exact_power(0.2, 0.2, n) <= ep.ALPHA + 1e-9
    # exact (conservative) power does not exceed the arcsine approximation by much
    assert power <= ep.arcsine_power(p1, p2, n) + 0.05
    m = ep.exact_min_n(p1, p2, 0.8, n_start=5, n_stop=80)
    assert m is not None and ep.exact_power(p1, p2, m) >= 0.8


def test_exact_fisher_agrees_with_scipy_on_one_table() -> None:
    from scipy.stats import fisher_exact

    n = 12
    reject = ep._fisher_reject_matrix(n, 0.05)
    for x1, x2 in [(0, 8), (1, 7), (3, 3), (2, 9), (10, 2)]:
        p = fisher_exact([[x1, n - x1], [x2, n - x2]])[1]
        assert bool(reject[x1, x2]) == (p <= 0.05)


def test_design_effect_and_cultures() -> None:
    assert ep.design_effect(50, 0.0) == 1.0
    assert math.isclose(ep.design_effect(50, 0.1), 5.9)
    out = ep.cultures_and_cells(348, 50, 0.05)
    assert out["design_effect"] == 3.45
    assert out["cultures_per_arm"] == math.ceil(348 * 3.45 / 50)
    assert out["cells_per_arm"] == out["cultures_per_arm"] * 50
    # more clustering never needs fewer cells; the floor of 3 cultures applies
    assert ep.cultures_and_cells(348, 50, 0.3)["cells_per_arm"] > out["cells_per_arm"]
    assert ep.cultures_and_cells(5, 50, 0.0)["cultures_per_arm"] == ep.MIN_CULTURES_PER_ARM


def test_quote_verification_normalises_whitespace_and_fails_loudly() -> None:
    text = "Cytogenetic analysis  revealed a normal\nkaryotype in 88 cells (46, XX) and various aneuploidies in 12 cells."
    assert ep.verify_quote("Cytogenetic analysis revealed a normal karyotype in 88 cells", text)
    assert not ep.verify_quote("in 13 cells", text)
    with pytest.raises(AssertionError):
        ep.assert_quote("in 13 cells", text, "unit")


def test_sources_are_complete_and_offline_verification_works() -> None:
    ids = [s["id"] for s in ep.SOURCES]
    assert len(ids) == len(set(ids))
    assert all(s["pmid"].isdigit() and s["quote"] for s in ep.SOURCES)
    verified = ep.verify_all_sources({s["id"]: " " + s["quote"] + " " for s in ep.SOURCES})
    assert all(v["verified"] for v in verified)
    with pytest.raises(AssertionError):
        ep.verify_all_sources({s["id"]: "unrelated text" for s in ep.SOURCES})
    # every scenario baseline points to a declared source
    for sc in ep.SCENARIOS:
        assert set(sc["sources"]) <= set(ids)
    # the BUB1B baseline equals the quoted 12 aneuploid of 100 cells
    bub = next(s for s in ep.SOURCES if s["id"] == "BUB1B_MVA1_patient")
    assert bub["value"] == {"aneuploid": 12, "total": 100}
    assert ep.SCENARIOS[0]["p0"] == 0.12


def test_table_without_exact_and_markdown() -> None:
    rows = ep.build_table(exact=False)
    assert len(rows) == len(ep.SCENARIOS) * len(ep.RELATIVE_REDUCTIONS) * len(ep.POWERS)
    r = next(x for x in rows if x["baseline_fraction"] == 0.12 and x["relative_reduction"] == 0.5 and x["power_target"] == 0.8)
    assert r["drug_arm_fraction"] == 0.06
    assert r["n_cells_per_arm_independent"] == ep.arcsine_n_per_arm(0.12, 0.06, power=0.8)
    assert len(r["clustered"]) == len(ep.CELLS_PER_CULTURE) * len(ep.ICC_VALUES)
    md = ep.markdown_summary(rows)
    assert md.startswith("Planning numbers, not results") and md.count("\n| 0.") == len(rows)


def test_written_json_if_present_is_consistent() -> None:
    path = Path("out/empirical/power_analysis.json")
    if not path.exists():
        pytest.skip("power_analysis.json not generated yet")
    raw = path.read_text(encoding="utf-8")
    data = json.loads(raw)
    assert raw == json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    assert "NOT RESULTS" in data["status"]
    assert all(s["verified"] for s in data["sources"])
    assert data["markdown_summary"].startswith("Planning numbers, not results")
    first = data["table"][0]
    assert first["n_cells_per_arm_independent"] == ep.arcsine_n_per_arm(first["baseline_fraction"], first["drug_arm_fraction"], power=first["power_target"])
