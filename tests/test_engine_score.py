"""Offline synthetic tests of track2.engine_score (no network, no real data)."""
import numpy as np

from track2 import engine_score as es
from track2 import empirical_literature_control as elc


def test_auc_precision_recall_ties():
    y = np.array([1, 0, 0, 1, 0])
    assert es.auc_mw(np.array([5.0, 1, 2, 4, 3]), y) == 1.0
    assert es.auc_mw(np.zeros(5), y) == 0.5
    assert np.isnan(es.auc_mw(np.zeros(3), np.zeros(3)))
    assert es.precision_at_k(np.array([5.0, 1, 2, 4, 3]), y, 2) == 1.0
    # all tied: expected hits at k = k * prevalence
    assert abs(es.precision_at_k(np.zeros(5), y, 2) - 0.4) < 1e-12
    assert abs(es.recall_at_k(np.zeros(5), y, 2) - 0.4) < 1e-12
    # agrees with the earlier proxy's exact tie-averaged precision
    sc = {"a": 1.0, "b": 1.0, "c": 0.0, "d": 1.0}
    arr = np.array([sc[k] for k in "abcd"])
    assert abs(es.precision_at_k(arr, np.array([1, 0, 0, 1]), 2) - elc.precision_at_k(sc, {"a", "d"}, list("abcd"), 2)) < 1e-12
    assert es.random_metrics(y) == (0.5, 0.4, 1.0)


def test_logit_fit_recovers_slope_and_weights():
    rng = np.random.default_rng(0)
    x = rng.normal(size=20000)
    y = (rng.random(20000) < 1 / (1 + np.exp(-(0.5 + 2.0 * x)))).astype(float)
    a, b = es.logit_fit(x, y)
    assert abs(a - 0.5) < 0.1 and abs(b - 2.0) < 0.15
    a2, b2 = es.logit_fit(x, y, np.full(len(x), 3.0))
    assert abs(b2 - b) < 1e-3


def test_folds_are_by_disease_and_deterministic():
    ds = [f"d{i}" for i in range(17)]
    f1, f2 = es.make_folds(ds), es.make_folds(ds)
    assert f1 == f2 and len(f1) == 5
    flat = [d for f in f1 for d in f]
    assert sorted(flat) == sorted(ds) and len(set(flat)) == 17
    assert len(es.make_folds([f"d{i}" for i in range(6)])) < 5


def test_class_never_split_and_oof_signal():
    rng = np.random.default_rng(1)
    n_drugs, ds = 60, [f"d{i}" for i in range(10)]
    classes = np.array([f"C{i // 3}" for i in range(n_drugs)])
    X, y = {}, {}
    for d in ds:
        pos = rng.choice(n_drugs, 3, replace=False)
        lab = np.zeros(n_drugs); lab[pos] = 1
        f = rng.normal(size=(n_drugs, 3)); f[:, 0] += 3 * lab
        X[d], y[d] = f, lab
    folds = es.make_folds(ds)
    # exclusion mask used in training never keeps a class that holds a test-fold gold drug
    for fold in folds:
        bad = es._gold_classes(fold, y, classes)
        keep = es._class_mask(classes, bad)
        assert not (set(classes[keep]) & bad)
    pred = es.oof_predictions(X, y, classes, folds)
    assert set(pred) == set(ds)
    aucs = [es.auc_mw(pred[d]["score"], y[d]) for d in ds]
    assert np.mean(aucs) > 0.8
    assert all((0 <= pred[d]["platt"]).all() and (pred[d]["platt"] <= 1).all() for d in ds)


def test_literature_counts_and_features():
    dre = elc.drug_regex(["rapamycin", "aspirin"])
    papers = [
        {"pid": "MED:1", "year": 2010, "text": "CFTR study. Rapamycin improved lung function in patients with CFTR mutations."},
        {"pid": "MED:2", "year": 2011, "text": "CFTR in mice. Aspirin is toxic and caused adverse effects in mice."},
        {"pid": "MED:3", "year": 2012, "text": "No gene here. Rapamycin improved outcomes in patients."},
        {"pid": "MED:4", "year": 2013, "text": "cftr lower case only. Rapamycin was used."},
    ]
    lit = es.literature_counts(papers, ["CFTR"], dre)
    c = lit["counts"]
    assert lit["n_papers_with_gene"] == 2
    assert c["rapamycin"]["n_papers"] == 1 and c["rapamycin"]["n_benefit"] == 1 and c["rapamycin"]["n_human"] == 1
    assert c["aspirin"]["n_papers"] == 1 and c["aspirin"]["n_harm"] >= 1 and c["aspirin"]["n_benefit"] == 0
    names = ["aspirin", "rapamycin", "zzz"]
    X = es.disease_matrix(lit, names, {"aspirin": 1000, "rapamycin": 10, "zzz": 0}, 10_000, None)
    # same co-mention, lower popularity -> larger excess
    assert X[1, 5] > X[0, 5] > X[2, 5]
    assert es.arm_columns(False, False) == [0, 1, 2, 3] and es.arm_columns(True, True) == [0, 1, 2, 3, 4, 5, 6, 7]


def test_query_has_no_drug_and_fetch_uses_cache(tmp_path):
    q = es.build_gene_query(["GBA", "CFTR"], 2018)
    assert 'TITLE_ABS:"CFTR"' in q and "PUB_YEAR:[1800 TO 2018]" in q
    calls = []

    def getter(url, params):
        page = 0 if params["cursorMark"] == "*" else int(params["cursorMark"])
        calls.append(page)
        items = [{"source": "MED", "id": str(page * 100 + i), "pubYear": "2000", "title": "t", "abstractText": "a"} for i in range(100)]
        return {"resultList": {"result": items}, "nextCursorMark": str(page + 1)}

    p = es.fetch_papers(["GBA"], 2018, tmp_path, getter)
    assert len(p) == 300 and len({x["pid"] for x in p}) == 300 and calls == [0, 1, 2]
    stuck = es.fetch_papers(["GBB"], 2018, tmp_path, lambda u, q: {"resultList": {"result": [{"source": "MED", "id": "1", "abstractText": "a"}] * 100},
                                                                      "nextCursorMark": "*"})
    assert len(stuck) == 1  # a cursor that does not advance cannot loop forever
    n = len(calls)
    es.fetch_papers(["GBA"], 2018, tmp_path, getter)
    assert len(calls) == n


def test_prior_features():
    targets = {"d1": {"A", "B"}, "d2": {"Z"}}
    ppi = {"A": {"X"}, "C": {"Z"}}
    f = es.prior_features(["A", "C"], ["d1", "d2", "d3"], targets, ppi)
    assert f[0, 0] > 0 and f[2].sum() == 0 and f[1, 0] == 0 and f[1, 1] > 0


def test_analyse_decide_and_bootstrap_determinism():
    rng = np.random.default_rng(2)
    n_drugs, ds = 80, [f"d{i}" for i in range(12)]
    classes = np.array([f"C{i}" for i in range(n_drugs)])
    X, y, pop, prx = {}, {}, {}, {}
    for d in ds:
        lab = np.zeros(n_drugs); lab[rng.choice(n_drugs, 4, replace=False)] = 1
        f = rng.normal(size=(n_drugs, 8)); f[:, 1] += 2 * lab
        X[d], y[d] = f, lab
        pop[d] = rng.normal(size=n_drugs)
        prx[d] = rng.integers(0, 3, n_drugs).astype(float)
    arms = {"A": (True, True), "D": (False, False)}
    r1 = es.analyse(ds, X, y, pop, prx, classes, arms)
    r2 = es.analyse(ds, X, y, pop, prx, classes, arms)
    assert r1 == r2
    a = r1["arms"]["A"]
    assert a["auc"]["estimate"] > 0.7 and a["auc"]["ci95"][0] <= a["auc"]["estimate"] <= a["auc"]["ci95"][1]
    assert a["auc_minus_popularity"]["estimate"] > 0.05
    dec = es.decide(a)
    assert set(dec) == {"criterion_1_auc_minus_popularity_ge_0.05_and_ci_lower_gt_0", "criterion_2_calibration_slope_in_0.8_1.2",
                        "criterion_3_precision_at_10_above_best_baseline", "verdict"}
    assert isinstance(dec["criterion_1_auc_minus_popularity_ge_0.05_and_ci_lower_gt_0"]["met"], bool)
    assert 0 < a["calibration_platt"]["brier"]["estimate"] < 0.2


def test_decide_reports_failure_whatever_the_result():
    arm = {"auc_minus_popularity": {"estimate": 0.01, "ci95": [-0.1, 0.1]},
           "calibration_platt": {"calibration_slope": {"estimate": 0.5}},
           "precision_at_10_minus_best_baseline": {"estimate": 0.0}}
    d = es.decide(arm)
    assert d["verdict"].startswith("not supported")
    assert not any(d[k]["met"] for k in d if k != "verdict")
