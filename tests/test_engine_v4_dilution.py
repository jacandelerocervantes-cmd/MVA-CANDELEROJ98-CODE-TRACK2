"""Offline synthetic tests of engine_v4_dilution (no network, no cache files, no model files)."""
import math

import numpy as np
import pytest

from track2 import engine_v4_dilution as dl


def test_subsample_size_floor_and_cap():
    assert dl.n_sub(300, 1.0) == 300
    assert dl.n_sub(300, 0.3) == 90
    assert dl.n_sub(300, 0.1) == 30
    assert dl.n_sub(300, 0.03) == 9
    assert dl.n_sub(300, 0.01) == 5            # floor binds (3 -> 5)
    assert dl.n_sub(11, 0.1) == 5
    assert dl.n_sub(3, 0.3) == 3               # never more than n
    assert dl.n_sub(0, 0.5) == 0


def test_subsample_deterministic_sorted_and_replicates_differ():
    a = dl.subsample_indices(300, 0.1, 7, 2, 0)
    assert np.array_equal(a, dl.subsample_indices(300, 0.1, 7, 2, 0))
    assert len(a) == 30 and len(set(a.tolist())) == 30 and np.all(np.diff(a) > 0)
    assert not np.array_equal(a, dl.subsample_indices(300, 0.1, 7, 2, 1))
    assert not np.array_equal(a, dl.subsample_indices(300, 0.1, 8, 2, 0))
    assert np.array_equal(dl.subsample_indices(40, 1.0, 1, 0, 0), np.arange(40))


def test_sample_keys_one_replicate_at_full_and_five_otherwise():
    ks = dl.sample_keys()
    assert ks[0] == (0, 0) and len(ks) == 1 + 4 * dl.N_REPS
    assert dl.key_of(0, 3) == (0, 0) and dl.key_of(2, 3) == (2, 3)


def test_design_constants_fixed_by_protocol():
    assert dl.FRACTIONS == [1.0, 0.3, 0.1, 0.03, 0.01] and dl.FLOOR == 5 and dl.N_REPS == 5 and dl.SEED == 20261015


def test_unvalidated_registry_complete():
    for k, p in dl.UNVALIDATED.items():
        assert p["status"] in ("VALIDATED", "UNVALIDATED"), k
        assert "conservative_default" in p, k


def test_cache_only_getter_refuses_network():
    with pytest.raises(RuntimeError):
        dl.cache_only_getter("u", {})


def test_break_even_definition():
    g = np.logspace(math.log10(5), math.log10(300), 200)
    be = dl.break_even(g, np.where(g >= 40, 0.05, -0.05), np.where(g >= 20, 0.1, -0.1))
    assert 40 <= be["value"] < 42 and 20 <= be["point_estimate_crossing"] < 21
    # a dip after the first positive stretch: the smallest N with a positive lower bound ALL the way up
    low = np.where((g >= 10) & (g < 50), 0.05, np.where(g >= 100, 0.05, -0.05))
    assert dl.break_even(g, low, low)["value"] >= 100
    assert dl.break_even(g, np.full(200, 0.1), np.full(200, 0.1))["value"] is None
    assert "below the smallest" in dl.break_even(g, np.full(200, 0.1), np.full(200, 0.1))["statement"]
    assert "not reached" in dl.break_even(g, np.full(200, -0.1), np.full(200, -0.1))["statement"]


def test_placement_relative_to_break_even():
    be = {"value": 40.0, "statement": "lower bound above 0 from this volume upward"}
    assert dl.placement(300, be).startswith("above") and "retrieval cap" in dl.placement(300, be)
    assert dl.placement(11, be).startswith("below")
    assert "not reached" in dl.placement(300, {"value": None, "statement": "not reached: x"})


def _synthetic_results(n_dis=12, slope=0.1, noise=0.0):
    """Frozen-style results where the model AUC rises with log10(n): gap = slope * (log10 n - 1) (+ noise)."""
    rng = np.random.default_rng(3)
    res, members = {}, {}
    for i in range(n_dis):
        d = f"d{i}"
        members[f"g{i}"] = [d]
        ns = {fi: [dl.n_sub(300, f)] * (1 if fi == 0 else dl.N_REPS) for fi, f in enumerate(dl.FRACTIONS)}
        met = {}
        for m in dl.MODELS:
            met[m] = {fi: [(0.5 + slope * (math.log10(n) - 1) + noise * rng.normal(), 0.5, 0.5) for n in v] for fi, v in ns.items()}
        res[d] = {"n_sub": ns, "metrics": met, "popularity": (0.5, 0.1, 0.1)}
    return res, members


def test_gap_curve_recovers_slope_and_break_even():
    res, members = _synthetic_results(slope=0.1, noise=0.0)
    fit = dl.fit_gap_curve(res, members, "B0", n_boot=200)
    assert fit["slope_per_log10_papers"] == pytest.approx(0.1, abs=1e-9)
    assert fit["intercept"] == pytest.approx(-0.1, abs=1e-9)
    # gap = 0.1 * (log10 n - 1) crosses 0 at n = 10; zero noise => the lower bound equals the fit
    be = fit["break_even"]
    assert 10 <= be["value"] < 10.6
    assert fit["n_observations"] == 12 * (1 + 4 * dl.N_REPS)


def test_gap_curve_with_no_signal_has_no_break_even():
    res, members = _synthetic_results(slope=0.0, noise=0.05)
    fit = dl.fit_gap_curve(res, members, "E", n_boot=300)
    assert fit["break_even"]["value"] is None and "not reached" in fit["break_even"]["statement"]


def test_table_for_set_columns_and_gap_sign():
    n = 8
    res, members = _synthetic_results(n_dis=n, slope=0.1, noise=0.0)
    rng = np.random.default_rng(0)
    y = {}
    for i in range(n):
        v = np.zeros(60)
        v[rng.choice(60, 3, replace=False)] = 1
        y[f"d{i}"] = v
    t = dl.table_for_set(res, members, y)
    assert t["n_groups"] == n and set(t["fractions"]) == {"1.0", "0.3", "0.1", "0.03", "0.01"}
    row = t["fractions"]["1.0"]["B0"]
    assert row["auc_gap_vs_popularity"]["estimate"] == pytest.approx(0.1 * (math.log10(300) - 1), abs=1e-5)
    low = t["fractions"]["0.01"]["B0"]["auc_gap_vs_popularity"]["estimate"]
    assert low == pytest.approx(0.1 * (math.log10(5) - 1), abs=1e-5) and low < row["auc_gap_vs_popularity"]["estimate"]
    assert t["fractions"]["0.01"]["n_papers_mean"] == 5 and t["fractions"]["0.01"]["B0"]["n_replicates"] == dl.N_REPS
    assert t["fractions"]["1.0"]["B0"]["n_replicates"] == 1


def test_thinned_features_has_no_gold_argument_and_f1_is_full():
    import inspect
    assert "gold" not in " ".join(inspect.signature(dl.thinned_features).parameters)
    assert "y" not in inspect.signature(dl.thinned_features).parameters
    assert list(inspect.signature(dl.score_diseases).parameters)[:4] == ["diseases", "feats", "heads", "y"]


def test_heldout_lock_refuses_second_scoring(tmp_path, monkeypatch):
    lock = tmp_path / "x.lock"
    lock.write_text("{}")
    monkeypatch.setattr(dl, "LOCK", lock)
    with pytest.raises(RuntimeError, match="lock"):
        dl.run_heldout()
