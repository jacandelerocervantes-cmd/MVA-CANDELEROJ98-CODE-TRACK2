"""Offline synthetic tests of engine_v3 (no network, no cache files)."""
import json

import numpy as np
import pytest

from track2 import engine_score as es
from track2 import engine_v3 as v3


def _diseases(n_new=27):
    ds = [{"name": f"v1_{i}", "origin": "v1"} for i in range(17)] + [{"name": v3.WILMS, "origin": "v2_new"}]
    ds += [{"name": f"new_{i:02d}", "origin": "v2_new"} for i in range(n_new)]
    return ds


def test_split_rule_sizes_and_determinism():
    s = v3.make_split(_diseases())
    assert len(s["dev"]) == 17 + 1 + 13 and len(s["heldout"]) == 14
    assert v3.WILMS in s["dev"] and not set(s["dev"]) & set(s["heldout"])
    assert s == v3.make_split(list(reversed(_diseases())))


def test_split_matches_protocol_text():
    from pathlib import Path
    txt = (Path(__file__).resolve().parents[1] / "protocol" / "addendum_engine_v3.md").read_text(encoding="utf-8")
    b = json.loads((Path(__file__).resolve().parents[1] / "out" / "engine" / "benchmark_v2.json").read_text(encoding="utf-8"))
    s = v3.make_split(b["diseases"])
    assert len(s["heldout"]) == 14 and all(h in txt for h in s["heldout"])


def test_f1_idf_penalises_popular_drug():
    c = np.array([5.0, 5.0]); nb = np.array([2.0, 2.0]); nh = np.array([1.0, 1.0])
    pop = np.array([100.0, 1_000_000.0])
    f = v3.f1_block(c, nb, nh, pop, 50, 40_000_000.0)
    assert f[0, 0] > f[1, 0] and f[0, 3] > f[1, 3]


def test_class_genes_and_keywords():
    genes = [{"symbol": "A", "mechanism_class_general": "mitochondrial"}, {"symbol": "B", "mechanism_class_general": "mitochondrial"},
             {"symbol": "C", "mechanism_class_general": "other"}, {"symbol": "D", "mechanism_class_general": "metabolic_enzyme"}]
    cls, syms = v3.disease_class_genes(genes)
    assert cls == "mitochondrial" and syms == ["A", "B"]
    assert v3.disease_class_genes([{"symbol": "X", "mechanism_class_general": "other"}]) == (None, [])
    assert v3.keyword_regex("other") is None and v3.keyword_regex("mitochondrial").search("Mitochondrial dysfunction")


def test_class_counts_synthetic():
    from track2 import empirical_literature_control as elc
    names = ["aspirin", "zinc"]
    dre = elc.drug_regex(names)
    papers = [{"pid": "p1", "text": "GENE1 and GENE2 in mitochondrial disease. Aspirin improves mitochondrial function."},
              {"pid": "p2", "text": "GENE1 only. Zinc was tested in patients."}]
    cc = v3.class_counts(papers, ["GENE1", "GENE2"], ["GENE2"], v3.keyword_regex("mitochondrial"), dre, es.name_gate(names))
    assert cc["cg"].get("aspirin") == 1 and "zinc" not in cc["cg"] and cc["kw"].get("aspirin") == 1 and cc["n_cg"] == 1


def test_f3_leave_one_out_does_not_use_own_label():
    k = 3
    cls_id = np.array([0, 0, 1, 1, 2, 2])
    y = {"a": np.array([1, 0, 0, 0, 0, 0.]), "b": np.array([0, 0, 0, 0, 1, 0.]), "c": np.array([1, 0, 0, 0, 0, 0.])}
    ref = v3.F3Ref(["a", "b", "c"], y, cls_id, k)
    full, loo_a = ref.table(None), ref.table("a")
    assert full[0] > 0
    only_bc = v3.F3Ref(["b", "c"], y, cls_id, k).table(None)
    assert np.allclose(loo_a, only_bc)
    assert np.allclose(v3.F3Ref([], {}, cls_id, k).table(None), 0)


def test_pick_C_one_se_rule():
    grid = (0.01, 0.1, 1.0)
    aucs = {0.01: {"a": 0.70, "b": 0.74, "c": 0.72}, 0.1: {"a": 0.72, "b": 0.76, "c": 0.74}, 1.0: {"a": 0.72, "b": 0.76, "c": 0.75}}
    C, diag = v3.pick_C(grid, aucs)
    assert C == 0.1 and diag["best_C"] == 1.0


def test_select_config_prefers_simpler_within_se():
    stats = {"base": {"mean_auc": 0.70, "se": 0.03, "complexity": (0, 8)},
             "F1": {"mean_auc": 0.71, "se": 0.03, "complexity": (1, 12)},
             "F1+F2": {"mean_auc": 0.73, "se": 0.03, "complexity": (2, 15)},
             "F3": {"mean_auc": 0.60, "se": 0.03, "complexity": (1, 9)}}
    s = v3.select_config(stats)
    assert s["best_by_mean"] == "F1+F2" and s["selected"] == "base"
    stats["base"]["mean_auc"] = 0.60
    assert v3.select_config(stats)["selected"] == "F1"


def _synthetic(n_dis=12, n_drug=200, seed=0):
    rng = np.random.default_rng(seed)
    classes = np.array([f"c{i % 20}" for i in range(n_drug)])
    cls_names, cls_id = np.unique(classes, return_inverse=True)
    Z, y = {}, {}
    for i in range(n_dis):
        yy = np.zeros(n_drug); yy[rng.choice(n_drug, 8, replace=False)] = 1
        X = rng.normal(size=(n_drug, 15))
        X[:, 1] += 1.5 * yy
        X[:, 9] += 0.8 * yy
        Z[f"d{i:02d}"], y[f"d{i:02d}"] = X, yy
    return Z, y, classes, cls_id, len(cls_names)


def test_fit_cv_freeze_roundtrip_and_signal():
    Z, y, classes, cls_id, k = _synthetic()
    dev = sorted(Z)
    res, chosen, folds = v3.outer_cv(("F1", "F3"), dev, Z, y, classes, cls_id, k, (0.1, 1.0))
    assert set(res) == set(dev) and len(chosen) == len(folds)
    auc = np.mean([es.auc_mw(res[d]["score"], y[d]) for d in dev])
    assert auc > 0.75
    m = v3.fit_model(("F1", "F3"), 1.0, dev, Z, y, cls_id, k)
    fr = v3.freeze(m, (0.0, 1.0), [str(c) for c in np.unique(classes)], {"dev": dev, "heldout": []})
    t = v3.thaw(json.loads(json.dumps(fr)))
    s0 = v3.score_model(m, Z["d00"], cls_id)
    assert np.allclose(v3.predict_frozen(t, Z["d00"], cls_id)["score"], s0)
    fr["model"]["C"] = 99.0
    with pytest.raises(RuntimeError):
        v3.thaw(fr)


def test_base_config_matches_v1_scorer():
    Z, y, classes, cls_id, k = _synthetic(6)
    dev = sorted(Z)
    m = v3.fit_model((), 1.0, dev, Z, y, cls_id, k)
    score, _ = es._fit_scorer([Z[d][:, :8] for d in dev], [y[d] for d in dev], [np.ones(len(classes), bool)] * len(dev))
    assert np.allclose(v3.score_model(m, Z["d01"], cls_id), score(Z["d01"][:, :8]), atol=1e-6)


def test_rank_table_and_evaluate():
    names = ["a", "b", "c"]
    assert v3.rank_table(np.array([0.1, 0.9, 0.5]), names, ["b", "a", "zzz"]) == {"b": {"rank": 1, "score": 0.9}, "a": {"rank": 3, "score": 0.1}, "zzz": None}
    Z, y, classes, cls_id, k = _synthetic(6)
    ds = sorted(Z)
    m = v3.fit_model((), 1.0, ds, Z, y, cls_id, k)
    pr = {d: {"score": v3.score_model(m, Z[d], cls_id), "platt": 1 / (1 + np.exp(-v3.score_model(m, Z[d], cls_id))), "prev": np.full(len(cls_id), 0.04)} for d in ds}
    pop = {d: np.random.default_rng(1).random(len(cls_id)) for d in ds}
    ev = v3.evaluate(ds, {"m": pr, "m2": pr}, y, pop, pop, reference="m")
    assert "decision" in ev["models"]["m"] and ev["models"]["m2"]["auc_minus_m"]["estimate"] == 0.0
