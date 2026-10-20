import random

from track2 import empirical_literature_control as m

UNI = ["aspirin", "sirolimus", "metformin", "ivacaftor", "caffeine", "warfarin", "digoxin",
       "atenolol", "baclofen", "lithium", "tamoxifen", "heparin"]


def test_sentence_scoring():
    dre = m.drug_regex(UNI)
    abs_ = [("1", "Sirolimus improved mTOR signalling in cells. Aspirin was given. "
                  "Metformin-like effects were absent. IVACAFTOR restored CFTR function and sirolimus rescued CFTR.")]
    c = m.score_abstracts(abs_, ["mTOR", "CFTR"], dre)
    assert c == {"sirolimus": 2, "ivacaftor": 1}


def test_whole_word_and_requires_all_three():
    dre = m.drug_regex(UNI)
    c = m.score_abstracts([("1", "Prometformin improved CFTR. Aspirin improved outcome. Lithium binds CFTR.")],
                          ["CFTR"], dre)
    assert c == {}


def test_auc_and_precision():
    sc = {"aspirin": 3, "sirolimus": 1}
    assert m.auc_ties(sc, {"aspirin", "sirolimus"}, UNI) == 1.0
    assert m.auc_ties({}, {"aspirin"}, UNI) == 0.5
    assert m.precision_at_k(sc, {"aspirin"}, UNI, 10) == 0.1
    assert abs(m.precision_at_k({}, {"aspirin", "warfarin"}, UNI, 10) - 2 / 12) < 1e-9


def test_wilson_and_boot():
    p, lo, hi = m.wilson(5, 10)
    assert lo < p < hi
    b = m.boot_ci([0.2, 0.4, 0.6], random.Random(1), 200)
    assert b["lo"] <= b["mean"] <= b["hi"]


def test_evaluate_and_aggregate_deterministic():
    sc = {"aspirin": 4, "sirolimus": 2, "caffeine": 1}
    r1 = m.evaluate_disease(sc, {"aspirin", "sirolimus"}, UNI, random.Random(3), 200)
    r2 = m.evaluate_disease(sc, {"aspirin", "sirolimus"}, UNI, random.Random(3), 200)
    assert r1 == r2 and r1["gold_hits"] == 2 and r1["auc"] == 1.0
    assert r1["enrichment_hits"] > 1
    assert r1["nongold_pass"] == 1
    a = m.aggregate({"d1": r1, "d2": r2}, random.Random(1))
    assert a["n_diseases"] == 2 and a["gold_hit_rate_pooled"]["rate"] == 1.0


def test_select_diseases_seeded():
    avail = {c[0]: {"x"} for c in m.CANDIDATES}
    a = m.select_diseases(avail, 1, 5)
    assert a == m.select_diseases(avail, 1, 5) and len(a) == 5
    assert m.select_diseases({}, 1, 5) == []


def test_epmc_cache_offline(tmp_path):
    calls = []

    def getter(url, params):
        calls.append(params["page"])
        return {"resultList": {"result": [{"id": "9", "abstractText": "x."}]}}

    r = m.epmc_fetch('"A" OR "B"', tmp_path, getter=getter)
    assert r == [("9", "x.")] and m.epmc_fetch('"A" OR "B"', tmp_path, getter=getter) == r
    assert calls == [1]
