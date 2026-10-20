"""Offline synthetic tests of engine_v4_ablation (no network, no cache files, no model files)."""
import numpy as np
import pytest

from track2 import engine_v4 as v4
from track2 import engine_v4_ablation as ab
from track2 import engine_v4_dilution as dl


def test_column_groups_partition_b0_features():
    assert sorted(ab.LIT_COLS + ab.PRIOR_COLS) == list(range(8)) and not set(ab.LIT_COLS) & set(ab.PRIOR_COLS)
    assert ab.PRIOR_COLS == [4, 6, 7] and ab.LIT_COLS == [0, 1, 2, 3, 5]


def test_unvalidated_registry_complete():
    for k, p in ab.UNVALIDATED.items():
        assert p["status"] in ("VALIDATED", "UNVALIDATED"), k
        assert "conservative_default" in p, k


def _head(rng, X, w):
    y = (X @ w + rng.normal(size=len(X)) > 0).astype(float)
    h = v4.head_params(v4.fit_head("logit", X, y), list(range(8)))
    return h


def test_mean_replacement_and_frozen_scorers():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 8)) + 2.0
    h = _head(rng, X, rng.normal(size=8))
    mean = X.mean(0)
    assert np.abs(mean - np.asarray(h["mu"])).max() < 1e-12            # DEV mean == frozen standardisation mean
    po, lo, b0 = (ab.frozen_scorer(m, h, mean) for m in ("PRIOR_ONLY", "LIT_ONLY", "B0"))
    Xn = rng.normal(size=(50, 8)) + 5.0
    assert np.allclose(b0("d", Xn), v4.head_decision(h, Xn))
    # PRIOR_ONLY ignores literature columns completely; LIT_ONLY ignores prior columns completely
    Xl = Xn.copy()
    Xl[:, ab.LIT_COLS] = rng.normal(size=(50, 5)) * 9
    assert np.allclose(po("d", Xn), po("d", Xl))
    Xp = Xn.copy()
    Xp[:, ab.PRIOR_COLS] = rng.normal(size=(50, 3)) * 9
    assert np.allclose(lo("d", Xn), lo("d", Xp))
    # replaced columns contribute exactly the intercept (standardised value 0)
    assert np.allclose(po("d", np.tile(mean, (3, 1))), h["intercept"])
    # input not mutated
    X0 = Xn.copy()
    po("d", Xn)
    assert np.array_equal(X0, Xn)


def test_refit_heads_group_folds_and_oof(monkeypatch):
    rng = np.random.default_rng(1)
    dev = [f"d{i}" for i in range(20)]
    members = {f"g{i}": [d] for i, d in enumerate(dev)}
    n = 80
    Z = {d: rng.normal(size=(n, 8)) for d in dev}
    y = {d: (Z[d][:, 4] + rng.normal(size=n) > 1).astype(float) for d in dev}
    classes = np.array([f"c{i}" for i in range(n)])
    oof, full = ab.fit_refit_heads(dev, members, Z, y, classes, ab.PRIOR_COLS)
    assert set(oof) == set(dev) and len(full["coef"]) == 3
    # a disease never scored by a head trained on its own fold: heads differ across folds
    assert len({id(h) for h in oof.values()}) == 5
    # our f=1 OOF scores equal engine_v4.oof_head (nested Platt does not change raw scores)
    assert ab.oof_f1_check(dev, members, Z, y, classes, ab.PRIOR_COLS, oof) < 1e-9


def _synthetic_res(prior_gap, lit_gap, b0_gap, n_dis=10):
    res, members = {}, {}
    for i in range(n_dis):
        d = f"d{i}"
        members[f"g{i}"] = [d]
        ns = {fi: [dl.n_sub(300, f)] * (1 if fi == 0 else dl.N_REPS) for fi, f in enumerate(dl.FRACTIONS)}
        gaps = {"B0": b0_gap, "PRIOR_ONLY": prior_gap, "LIT_ONLY": lit_gap}
        res[d] = {"n_sub": ns, "popularity": (0.5, 0.1, 0.1),
                  "metrics": {m: {fi: [(0.5 + gaps[m] + 0.001 * i, 0.2, 0.2)] * len(v) for fi, v in ns.items()} for m in gaps}}
    return res, members


def _tab(prior_gap, lit_gap, b0_gap):
    res, members = _synthetic_res(prior_gap, lit_gap, b0_gap)
    return ab.table(res, members, {}, ("B0", "PRIOR_ONLY", "LIT_ONLY"), [("B0", "PRIOR_ONLY")])


def test_table_paired_difference():
    t = _tab(0.05, 0.02, 0.08)
    row = t["fractions"]["0.01"]
    assert row["paired"]["B0_minus_PRIOR_ONLY"]["auc"]["estimate"] == pytest.approx(0.03, abs=1e-5)
    assert row["PRIOR_ONLY"]["auc_gap_vs_popularity"]["estimate"] == pytest.approx(0.05 + 0.0045, abs=1e-5)
    assert set(t["fractions"]) == {"1.0", "0.3", "0.1", "0.03", "0.01"}


def test_decision_rule_four_outcomes():
    # prior reproduces B0 (identical), literature-only gap lower bound > 0 -> both
    d = ab.decide(_tab(0.08, 0.04, 0.08))
    assert d["H_PRIOR_gap_explained_by_prior"]["holds"] and d["H_LIT_literature_contributes_at_low_volume"]["holds"]
    assert d["outcome"].startswith("both")
    # prior alone far below B0 (B0 CI tight, excluded), literature positive -> H-LIT only
    d = ab.decide(_tab(0.0, 0.05, 0.08))
    assert not d["H_PRIOR_gap_explained_by_prior"]["holds"] and d["outcome"].startswith("H-LIT only")
    # literature-only gap negative, prior reproduces -> H-PRIOR only
    d = ab.decide(_tab(0.08, -0.05, 0.08))
    assert d["outcome"].startswith("H-PRIOR only") and d["H_LIT_literature_contributes_at_low_volume"]["status"] == "does not hold"
    d = ab.decide(_tab(0.0, -0.05, 0.08))
    assert d["outcome"].startswith("neither")


def test_h_lit_needs_both_fractions():
    t = _tab(0.08, 0.04, 0.08)
    t["fractions"]["0.03"]["LIT_ONLY"]["auc_gap_vs_popularity"]["ci95"] = [-0.01, 0.1]
    d = ab.decide(t)["H_LIT_literature_contributes_at_low_volume"]
    assert not d["holds"] and d["status"] == "partial"


def test_heldout_lock_refuses_second_scoring(tmp_path, monkeypatch):
    lock = tmp_path / "x.lock"
    lock.write_text("{}")
    monkeypatch.setattr(ab, "LOCK", lock)
    with pytest.raises(RuntimeError, match="lock"):
        ab.run_heldout()
