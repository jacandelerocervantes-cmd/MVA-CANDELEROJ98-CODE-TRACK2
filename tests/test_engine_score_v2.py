"""Offline synthetic tests of track2.engine_score_v2 (no network, no real data)."""
import json

import numpy as np
import pytest

from track2 import engine_score_v2 as v2


def test_gold_vectors_union_and_case():
    index = {"aspirin": 0, "ibuprofen": 1, "caffeine": 2, "zinc": 3}
    dis = {"gold1": [{"name": "Aspirin"}], "gold2": [{"name": "zinc"}], "gold3": [{"name": "caffeine"}, {"name": "Aspirin"}, {"name": "unknown"}]}
    g = v2.gold_vectors(dis, index, 4)
    assert g["gold1"].tolist() == [1, 0, 0, 0]
    assert g["gold3"].tolist() == [1, 0, 1, 0]
    assert g["union"].tolist() == [1, 0, 1, 0]
    assert g["gold2"].tolist() == [0, 0, 0, 1]


def test_mechanism_majority_rule():
    genes = [{"mechanism_class_general": c} for c in ["a", "a", "a", "b", "other"]]
    r = v2.disease_mech_class(genes)
    assert r["majority"] == "a" and r["plurality"] == "a"
    genes = [{"mechanism_class_general": c} for c in ["a", "a", "b", "b"]]  # 50 % is not a majority
    r = v2.disease_mech_class(genes)
    assert r["majority"] is None and r["plurality"] == "a"
    assert v2.disease_mech_class([]) == {"majority": None, "plurality": None}


def test_common_flags_rule():
    ds = [{"name": "Alzheimer disease", "inheritance": ["Autosomal dominant"]},
          {"name": "toy complex", "inheritance": ["Multigenic/multifactorial"]},
          {"name": "toy sporadic", "inheritance": ["Sporadic"]},
          {"name": "toy rare", "inheritance": ["Autosomal recessive"]}]
    f = v2.common_flags(ds)
    assert f["Alzheimer disease"] == {"S1": True, "S2": True}
    assert f["toy complex"] == {"S1": False, "S2": True}
    assert f["toy sporadic"]["S2"] is False and f["toy rare"]["S2"] is False


def test_request_guard_counts_and_aborts(tmp_path):
    calls = []
    g = v2.RequestGuard(lambda u, p: calls.append(1) or {"ok": 1}, tmp_path / "p.json", limit=3, every=1)
    for _ in range(3):
        assert g("u", {}) == {"ok": 1}
    assert json.loads((tmp_path / "p.json").read_text())["requests"] == 3
    with pytest.raises(RuntimeError):
        g("u", {})
    assert len(calls) == 3


def test_class_auc_table_groups_and_skips_nan():
    auc = {"d1": 0.9, "d2": 0.7, "d3": 0.5, "d4": float("nan"), "d5": 0.6}
    pop = {k: 0.5 for k in auc}
    mech = {"d1": {"majority": "x"}, "d2": {"majority": "x"}, "d3": {"majority": "y"}, "d4": {"majority": "y"}, "d5": {"majority": None}}
    t = v2.class_auc_table(auc, pop, mech, "majority", n_boot=50)
    assert set(t) == {"x", "y"} and t["x"]["n_diseases"] == 2 and t["y"]["n_diseases"] == 1
    assert abs(t["x"]["engine_auc"]["estimate"] - 0.8) < 1e-9
    assert abs(t["x"]["engine_minus_popularity"]["estimate"] - 0.3) < 1e-9


def _toy_per(n_dis=8, n_drugs=40, seed=0):
    rng = np.random.default_rng(seed)
    per = {}
    for i in range(n_dis):
        d = f"d{i}"
        g1 = np.zeros(n_drugs); g1[rng.choice(n_drugs, 3, replace=False)] = 1
        g3 = np.zeros(n_drugs) if i % 3 == 0 else g1.copy()
        X = rng.normal(size=(n_drugs, 8)); X[:, 1] += 2 * g1
        per[d] = {"gold1": g1, "gold2": g1.copy(), "gold3": g3, "union": np.maximum(g1, g3),
                  "primary": {"X": X, "popularity": rng.random(n_drugs), "proxy": rng.random(n_drugs)},
                  "temporal": {"X": X.copy(), "popularity": rng.random(n_drugs), "proxy": rng.random(n_drugs)}}
    classes = np.array([f"C{i // 2}" for i in range(n_drugs)])
    return per, classes


def test_analyse_set_and_criteria_structure():
    per, classes = _toy_per()
    r = v2.analyse_set(per, classes, sorted(per), "gold1")
    assert set(r["arms"]) == {"A", "B", "C", "D"} and r["label"] == v2.LABEL
    dec = r["decision_arm_A"]
    assert set(k for k in dec if k.startswith("criterion")) == {
        "criterion_1_auc_minus_popularity_ge_0.05_and_ci_lower_gt_0", "criterion_2_calibration_slope_in_0.8_1.2",
        "criterion_3_precision_at_10_above_best_baseline"}
    assert r["arms"]["A"]["auc"]["estimate"] > 0.6
    # gold 3 empty for d0, d3, d6: they are dropped
    r3 = v2.analyse_set(per, classes, sorted(per), "gold3")
    assert r3["diseases_excluded_empty_gold"] == ["d0", "d3", "d6"] and len(r3["diseases"]) == 5


def test_temporal_only_no_prior_arms():
    per, classes = _toy_per()
    t = v2.analyse_temporal(per, classes, sorted(per))
    assert set(t["all_gold2"]["arms"]) == {"B", "D"}
    assert t["excluding_gold1_drugs"] is None or set(t["excluding_gold1_drugs"]["arms"]) == {"B", "D"}
