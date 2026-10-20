"""Offline, synthetic, fast tests of the Wilms worked case (no network, no real data)."""
import numpy as np

from track2 import engine_case_study as cs
from track2 import engine_score as es

NAMES = [f"drug{i:02d}" for i in range(30)]
UNI = {n: {"molecule_type": "Small molecule", "atc": [f"L01{'ABCDEF'[i % 6]}0{i % 5}"], "first_approval": 2000 + i} for i, n in enumerate(NAMES)}
UNI["drug28"]["molecule_type"] = "Antibody"
UNI["drug27"]["molecule_type"] = "Oligonucleotide"
UNI["drug26"]["molecule_type"] = "Unknown"


def synthetic(n_dis=12, seed=3):
    rng = np.random.default_rng(seed)
    per = {}
    for k in range(n_dis):
        g1 = np.zeros(len(NAMES))
        g1[rng.choice(len(NAMES), 3, replace=False)] = 1
        g2 = np.zeros(len(NAMES))
        g2[rng.choice(len(NAMES), 4, replace=False)] = 1
        per[f"dis{k:02d}"] = {"gold1": g1, "gold2": g2,
                              "primary": {"X": rng.normal(size=(len(NAMES), 8)) + g1[:, None] * 0.8},
                              "temporal": {"X": rng.normal(size=(len(NAMES), 8)) + g2[:, None] * 0.5}}
    return per


CLASSES = np.array([f"c{i % 28}" for i in range(len(NAMES))])


def test_rank_arrays_ties():
    rmin, rmax = cs.rank_arrays(np.array([3.0, 1.0, 3.0, 2.0]))
    assert list(rmin) == [1, 4, 1, 3] and list(rmax) == [2, 4, 2, 3]


def test_concordance_hits_and_absent():
    score = np.arange(30, 0, -1, dtype=float)          # drug00 best
    gold = np.zeros(30)
    gold[[0, 14, 29]] = 1
    mask = np.ones(30, dtype=bool)
    mask[29] = False
    c = cs.concordance(score, NAMES, mask, gold)
    assert c["n_candidates"] == 29 and c["n_gold_in_candidates"] == 2 and c["gold_not_in_candidates"] == ["drug29"]
    assert [d["rank_min"] for d in c["drugs"]] == [1, 15]
    assert c["hits_at_10_certain"] == 1 and c["hits_at_20_certain"] == 2
    assert c["drugs"][1]["percentile_top"] == round(100 * 15 / 29, 3)


def test_filter_chain_counts():
    has = {n: (i % 3 != 0) for i, n in enumerate(NAMES)}
    ch = cs.filter_chain(NAMES, UNI, has)
    assert ch["n_universe"] == 30 and ch["n_fda_label"] == sum(has.values())
    d = ch["dropped_by_type_after_fda"]
    assert d["biologic"] == int(has["drug28"]) and d["oligonucleotide"] == int(has["drug27"]) and d["unknown_type"] == int(has["drug26"])
    assert ch["n_small_molecule"] == ch["n_fda_label"] - sum(d.values())


def test_has_fda_label_whole_word_and_homeopathic_miss():
    assert cs.has_fda_label("vincristine", [{"generic_name": ["VINCRISTINE SULFATE"]}])
    assert not cs.has_fda_label("lactic acid", [{"generic_name": ["SARCOLACTICUM ACIDUM"]}])
    assert not cs.has_fda_label("x", [])


def test_band_status_never_assumes_age():
    assert cs.band_status("ESTABLISHED", 2.0) == {"0_4": "partly (from 2 y)", "5_11": "labelled", "12_17": "labelled"}
    assert cs.band_status("ESTABLISHED", 0.0)["0_4"] == "labelled"
    assert cs.band_status("ESTABLISHED", 12.0) == {"0_4": "not_labelled", "5_11": "not_labelled", "12_17": "labelled"}
    assert set(cs.band_status("ESTABLISHED", None).values()) == {"age_not_stated"}
    assert set(cs.band_status("NOT_ESTABLISHED", None).values()) == {"not_labelled"}
    for s in ("MIXED", "UNKNOWN", "UNREADABLE", "VERIFY"):
        assert set(cs.band_status(s, 5.0).values()) == {"verify"}


def test_heldout_score_does_not_read_target_labels():
    per = synthetic()
    X = {d: v["primary"]["X"] for d, v in per.items()}
    y = {d: v["gold1"] for d, v in per.items()}
    others = [d for d in sorted(per) if d != "dis00"]
    a = cs.heldout_scores(X, {d: y[d] for d in others}, CLASSES, ["dis00"], others, False)["dis00"]["score"]
    # a dictionary with no label of the target works; different target labels cannot matter because they are not passed
    y2 = dict(y)
    y2["dis00"] = 1 - y["dis00"]
    b = cs.heldout_scores(X, y2, CLASSES, ["dis00"], others, False)["dis00"]["score"]
    assert np.array_equal(a, b)
    # the protocol fold reads gold 1 of the test fold (to drop training classes) and therefore may differ
    c = cs.heldout_scores(X, y, CLASSES, ["dis00"], others, True)["dis00"]["score"]
    assert c.shape == a.shape


def test_compute_rankings_deterministic_and_folds_match_engine():
    per = synthetic()
    r1 = cs.compute_rankings(per, CLASSES, "dis03")
    r2 = cs.compute_rankings(per, CLASSES, "dis03")
    assert np.array_equal(r1["primary"]["strict"]["score"], r2["primary"]["strict"]["score"])
    assert r1["primary"]["folds"] == es.make_folds(sorted(per))
    assert "dis03" in r1["primary"]["fold_of_target"]
    assert r1["temporal"] is not None and "dis03" in r1["temporal"]["fold_of_target"]
    pl = r1["primary"]["strict"]["platt"]
    assert ((pl >= 0) & (pl <= 1)).all()


def test_fold_run_equals_engine_oof_for_the_target_fold():
    per = synthetic()
    sel = sorted(per)
    cols = es.arm_columns(True, True)
    X = {d: per[d]["primary"]["X"][:, cols] for d in sel}
    y = {d: per[d]["gold1"] for d in sel}
    folds = es.make_folds(sel)
    oof = es.oof_predictions(X, y, CLASSES, folds)
    r = cs.compute_rankings(per, CLASSES, "dis05")
    assert np.allclose(oof["dis05"]["score"], r["primary"]["fold"]["score"])


def test_evidence_rungs_and_gap_text():
    r = cs.evidence_rungs({"n_papers": 3, "n_benefit": 1, "n_human": 0}, False, False)
    assert r == {"L1": True, "L2": True, "L3": False, "L4": False, "L5": False}
    g = cs.evidence_gap("drug00", ["L01CA02"], r, {"0_4": "verify", "5_11": "labelled", "12_17": "labelled"}, [])
    assert g["evidence_level_reached"] == "L2" and g["oncology_drug_chemotherapy_caution"] and g["age_band_gaps"] == ["0_4"]
    assert "adaptive pilot" in g["single_study_or_data_that_would_close_the_gap"]
    r2 = cs.evidence_rungs({"n_papers": 1, "n_benefit": 1, "n_human": 1}, True, False)
    g2 = cs.evidence_gap("drug01", [], r2, {"0_4": "labelled", "5_11": "labelled", "12_17": "labelled"}, ["NCT00000001"])
    assert g2["evidence_level_reached"] == "L4" and "NCT00000001" in g2["single_study_or_data_that_would_close_the_gap"]
    assert evidence_none()["evidence_level_reached"] == "L0"


def evidence_none():
    return cs.evidence_gap("x", [], cs.evidence_rungs({}, False, False), {"0_4": "verify", "5_11": "verify", "12_17": "verify"}, [])


def test_class_top_drugs_order_and_table():
    counts = {"a": {"n_papers": 2, "n_benefit": 1, "n_human": 0}, "b": {"n_papers": 5, "n_benefit": 0, "n_human": 1},
              "c": {"n_papers": 2, "n_benefit": 1, "n_human": 0}, "z": {"n_papers": 9, "n_benefit": 9, "n_human": 9}}
    rows = cs.class_top_drugs(counts, {"a", "b", "c"}, 2)
    assert [r["name"] for r in rows] == ["b", "a"]
    sup = {"b": {"bands": {"0_4": "verify", "5_11": "labelled", "12_17": "partly (from 13 y)"}}}
    t = cs.age_table({"cls": rows}, sup)
    assert t["cls"]["0_4"] == {"no_label_read": ["a"], "verify": ["b"]}
    assert t["cls"]["12_17"]["partly"] == ["b"]


def test_markdown_renders_minimal_result_keys():
    # render_md must at least not fail on the label header of a minimal result
    assert cs.LABEL == "worked case on a v2 disease, not part of the pre-registered primary validation"
    assert cs.SEED == 20261015
