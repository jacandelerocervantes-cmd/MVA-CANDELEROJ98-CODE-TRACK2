from __future__ import annotations

import numpy as np

from track2 import line_c_adjusted as L


def _data(seed: int = 1, n: int = 240):
    rng = np.random.default_rng(seed)
    lineage = np.array(["A"] * 80 + ["B"] * 80 + ["C"] * 80, dtype=object)
    offset = np.array([0.0, 1.0, 2.0])[np.repeat([0, 1, 2], 80)]
    score = rng.normal(size=n) + offset  # score differs between lineages
    prolif = rng.normal(size=n)
    return rng, lineage, offset, score, prolif


def test_within_perms_stay_in_group_and_are_permutations() -> None:
    g = np.array(["a"] * 5 + ["b"] * 7, dtype=object)
    p = L.within_perms(g, 20, np.random.default_rng(0))
    assert all(sorted(r) == list(range(12)) for r in p)
    assert (p[:, :5] < 5).all() and (p[:, 5:] >= 5).all()


def test_pure_lineage_confounding_is_removed() -> None:
    rng, lineage, offset, score, prolif = _data()
    y = offset + rng.normal(size=len(score)) * 0.5  # y depends on lineage only
    r = L.partial_perm(y, score, prolif, lineage, n_perm=200)
    assert abs(r["rho_adj"]) < 0.25
    assert stats_rho(y, score) > 0.3  # the unadjusted association is spurious but large


def stats_rho(a, b) -> float:
    from scipy import stats
    return float(stats.spearmanr(a, b).statistic)


def test_true_within_lineage_effect_is_kept_and_sign_matches() -> None:
    rng, lineage, offset, score, prolif = _data(2)
    y = 0.8 * score + offset + rng.normal(size=len(score)) * 0.5
    r = L.partial_perm(y, score, prolif, lineage, n_perm=200)
    assert r["rho_adj"] > 0.4 and r["p_perm"] < 0.02
    w = L.within_lineage(y, score, prolif, lineage, min_n=30)
    assert w["n_lineages"] == 3 and w["pooled_rho"] > 0.4 and w["p_random"] < 0.01


def test_proliferation_confounding_is_removed() -> None:
    rng, lineage, offset, score, prolif = _data(3)
    score = score + 2 * prolif  # score tracks proliferation
    y = prolif + rng.normal(size=len(score)) * 0.5  # y depends on proliferation only
    r = L.partial_perm(y, score, prolif, lineage, n_perm=200)
    assert abs(r["rho_adj"]) < 0.2
    assert stats_rho(y, score) > 0.5


def test_analyse_bh_survival_and_not_estimable() -> None:
    rng, lineage, offset, score, prolif = _data(4)
    Y = {"TRUE": 0.8 * score + offset + rng.normal(size=240) * 0.5, "NULL": rng.normal(size=240), "FEW": np.where(np.arange(240) < 50, 0.1, np.nan)}
    res = L.analyse(Y, score, prolif, lineage, n_perm=200, min_lines=100)
    assert not res["FEW"]["estimable"] and L.survives(res["FEW"], 0.3) == "not estimable"
    assert L.survives(res["TRUE"], 0.3) == "survives"
    assert L.survives(res["TRUE"], -0.3) == "does not survive"  # opposite original sign
    assert L.survives(res["NULL"], 0.3) == "does not survive"
    res2 = L.analyse(Y, score, prolif, lineage, exclude={"C"}, n_perm=100, min_lines=100)
    assert res2["TRUE"]["n_lines"] == 160


def test_lolo_reports_sign_retention() -> None:
    rng, lineage, offset, score, prolif = _data(5)
    y = 0.8 * score + offset + rng.normal(size=240) * 0.5
    o = L.lolo(y, score, prolif, lineage, 1.0, min_n=30, n_perm=100)
    assert o["n_omissions"] == 3 and o["n_sign_kept"] == 3
