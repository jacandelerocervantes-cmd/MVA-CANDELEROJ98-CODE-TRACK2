"""Offline synthetic tests of engine_v4_llm (no network, no cache files, no LLM, no model files)."""
import inspect
import json

import numpy as np
import pytest

from track2 import engine_v4 as v4
from track2 import engine_v4_dilution as dl
from track2 import engine_v4_llm as m


# ----------------------------------------------------------------------------- registry and constants
def test_unvalidated_registry_complete():
    assert v4.check_registry(m.UNVALIDATED) == []
    for k, p in m.UNVALIDATED.items():
        assert p["status"] == "UNVALIDATED" and "conservative_default" in p, k


def test_constants_follow_the_protocol():
    assert (m.TOP_K, m.MAX_PAPERS_PER_PAIR, m.MAX_SENT_PER_PAPER, m.SHARD_SIZE, m.KAPPA_GATE) == (30, 3, 2, 100, 0.6)
    assert m.PRIOR_COLS == [4, 6, 7] and m.MODEL_COLS["PL"] == [4, 6, 7, 8, 9, 10] and m.MODEL_COLS["PO"] == [4, 6, 7]
    assert m.ITEM_FIELDS == ("id", "drug", "disease", "sentence", "pmid")


def test_feature_construction_has_no_gold_argument():
    for fn in (m.llm_features, m.paper_drug_sentences, m.pair_papers, m.union_pair_papers, m.build_pair_index, m.make_items, m.consensus, m.label_table):
        bad = [p for p in inspect.signature(fn).parameters if "gold" in p.lower() or p.lower() in ("y", "labels_gold", "score")]
        assert not bad, (fn.__name__, bad)


# ----------------------------------------------------------------------------- candidates
def test_candidate_set_is_top_k_union_gold_with_provenance():
    score = np.arange(100, dtype=float)[::-1]          # index 0 has the highest score
    cand, prov = m.candidate_indices(score, [50, 99, 3], top_k=5)
    assert cand == sorted({0, 1, 2, 3, 4, 50, 99})
    assert prov[3]["top"] is True and prov[50]["top"] is False


# ----------------------------------------------------------------------------- thinning reuse
def _per_paper(n=40, seed=1, drugs=("aspirin", "sirolimus", "zinc")):
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        d = {}
        for g in drugs:
            if rng.random() < 0.3:
                d[g] = [f"S{i}-{g}-a.", f"S{i}-{g}-b."]
        out.append(d)
    return out


def test_pair_papers_caps_at_three_in_retrieval_order():
    pp = [{"a": ["x"]} for _ in range(10)]
    assert m.pair_papers(pp, range(10)) == {"a": [0, 1, 2]}
    assert m.pair_papers(pp, [7, 4, 9, 2, 5]) == {"a": [2, 4, 5]}       # order of the kept set, not of the argument


def test_union_serves_every_sample_exactly():
    """The annotated union must contain, for every sample, the first-3 kept papers of every drug, and picking the first 3 among the
    kept papers of the union must give the same papers as picking them among all papers (the reuse argument of protocol section 2)."""
    pp = _per_paper(60)
    papers = [{"pid": f"MED:{1000 + i}"} for i in range(60)]
    rank = 7
    index = m.build_pair_index(pp, papers, rank)
    for fi, rep in dl.sample_keys():
        kept = dl.subsample_indices(len(pp), dl.FRACTIONS[fi], rank, fi, rep)
        full = m.pair_papers(pp, kept)
        for drug, lst in full.items():
            from_union = [e["i"] for e in index[drug] if e["i"] in set(kept.tolist())][:m.MAX_PAPERS_PER_PAIR]
            assert from_union == lst, (fi, rep, drug)
            assert set(lst) <= {e["i"] for e in index[drug]}


def test_features_reuse_labels_of_kept_papers_only():
    pp = [{"aspirin": ["S0."]}, {"aspirin": ["S1."]}, {"aspirin": ["S2."]}, {"aspirin": ["S3."]}]
    papers = [{"pid": f"MED:{i}"} for i in range(4)]
    index = {"aspirin": [{"i": i, "pid": papers[i]["pid"], "sents": pp[i]["aspirin"]} for i in range(4)]}
    lab = {("D", "aspirin", str(i), f"S{i}."): {"relevant": True, "direction": "benefit", "reliability": "robust"} for i in range(4)}
    L = m.llm_features("D", index, {0, 1, 2, 3}, lab, {"aspirin": 0}, 2)
    assert np.allclose(L[0], [np.log1p(3), 0, np.log1p(3)]) and (L[1] == 0).all()     # cap of 3 papers
    L2 = m.llm_features("D", index, {2, 3}, lab, {"aspirin": 0}, 2)
    assert np.allclose(L2[0], [np.log1p(2), 0, np.log1p(2)])                         # only the kept papers
    lab[("D", "aspirin", "2", "S2.")] = {"relevant": False, "direction": "benefit", "reliability": "robust"}
    lab[("D", "aspirin", "3", "S3.")] = {"relevant": True, "direction": None, "reliability": "unresolved"}
    assert (m.llm_features("D", index, {2, 3}, lab, {"aspirin": 0}, 2) == 0).all()   # irrelevant and unresolved are not counted


# ----------------------------------------------------------------------------- items, shards
def _raw_items(n):
    return [{"disease": f"Dis{k % 7}", "drug": f"drug{k % 11}", "pmid": str(1000 + k), "sentence": f"Sentence number {k}."} for k in range(n)]


def test_shards_have_only_item_fields_at_most_100_and_are_deterministic():
    items = m.shuffle_and_id(_raw_items(257), "dev")
    m.assert_item_fields(items)
    assert all(not ({"gold", "score", "is_gold", "in_gold1"} & set(it)) for it in items)
    sh = m.shard(items)
    assert [len(s) for s in sh] == [100, 100, 57]
    assert items == m.shuffle_and_id(list(reversed(_raw_items(257))), "dev")            # independent of input order
    assert len({it["id"] for it in items}) == 257 and items[0]["id"].startswith("dev-")
    assert m.shuffle_and_id(_raw_items(5), "heldout")[0]["id"].startswith("ho-")
    with pytest.raises(AssertionError):
        m.assert_item_fields([{**items[0], "gold": 1}])


def test_dedupe_and_pid_roundtrip():
    its = _raw_items(3) + _raw_items(3)
    assert len(m.dedupe(its)) == 3
    assert m.pmid_of("MED:123") == "123" and m.pid_of("123") == "MED:123" and m.pmid_of("PPR:55") == "PPR:55" and m.pid_of("PPR:55") == "PPR:55"


def test_token_estimate_is_computed_from_the_text():
    items = m.shuffle_and_id(_raw_items(10), "dev")
    t = m.estimate_tokens(items)
    assert t["input_tokens_items"] == -(-len(json.dumps(items, ensure_ascii=False)) // 4) and t["output_tokens_two_tasks"] == 480


# ----------------------------------------------------------------------------- verbatim verification
def test_verify_sentences_statuses_and_squash_rule():
    texts = {"MED:1": "Title one. <i>Sirolimus</i> rescued  growth &amp; survival in mice.", "MED:2": None}
    its = [{"pmid": "1", "sentence": "Sirolimus rescued growth & survival in mice."},
           {"pmid": "1", "sentence": "Sirolimus cured everything."},
           {"pmid": "2", "sentence": "Anything."},
           {"pmid": "9", "sentence": "Missing paper."}]
    assert m.verify_sentences(its, texts) == ["verified", "not_found", "not_checkable", "not_checkable"]


def test_fetch_source_texts_batches_by_source_and_marks_missing():
    calls = []

    def getter(url, params):
        calls.append(params["query"])
        ids = [t.split(":")[1].split(")")[0].split(" ")[0] for t in params["query"].split(" OR ")]
        res = [{"source": "MED", "id": i, "title": f"T{i}", "abstractText": f"A{i}"} for i in ids if i != "3"]
        return {"resultList": {"result": res}}
    cache = {}
    out = m.fetch_source_texts([f"MED:{i}" for i in range(1, 121)], getter, cache)
    assert len(calls) == 3 and all("SRC:MED" in c for c in calls)                          # 120 ids / 50 per batch
    assert out["MED:1"] == "T1 A1" and out["MED:3"] is None
    m.fetch_source_texts(["MED:1", "MED:2"], getter, cache)
    assert len(calls) == 3                                                                 # served from the cache


# ----------------------------------------------------------------------------- labels and consensus
def test_consensus_relevance_is_a_conjunction_and_direction_a_two_of_three_vote():
    c = m.consensus("yes", "yes", "benefit", "benefit", "Drug improved survival.", "drug")
    assert c["relevant"] and c["direction"] == "benefit" and c["reliability"] in ("robust", "majority")
    assert not m.consensus("yes", "unclear", "benefit", "benefit", "Drug improved survival.", "drug")["relevant"]
    c3 = m.consensus("yes", "yes", "benefit", "harm", "The drug was discussed.", "drug")      # rule says unclear -> all differ
    assert c3["reliability"] == "unresolved" and c3["direction"] is None


def test_load_annotations_validates_ids_and_vocabulary(tmp_path):
    items = m.shuffle_and_id(_raw_items(3), "dev")
    rows = [{"id": it["id"], "direction": "benefit", "relevance": "yes"} for it in items]
    p = tmp_path / "a.json"
    p.write_text(json.dumps({"labels": rows}), encoding="utf-8")
    assert set(m.load_annotations(p, items)) == {it["id"] for it in items}
    p.write_text(json.dumps({"labels": rows[:2]}), encoding="utf-8")
    with pytest.raises(ValueError):
        m.load_annotations(p, items)
    rows[0]["relevance"] = "maybe"
    p.write_text(json.dumps({"labels": rows}), encoding="utf-8")
    with pytest.raises(ValueError):
        m.load_annotations(p, items)


# ----------------------------------------------------------------------------- agreement gate
def test_kappa_gate_and_generic_agreement():
    a = ["yes", "no", "unclear", "yes", "no", "yes"] * 5
    groups = [f"p{i}" for i in range(len(a))]
    r = m.agreement_generic(a, a, groups, n_boot=50)
    assert r["kappa"]["estimate"] == pytest.approx(1.0) and r["accuracy"]["estimate"] == 1.0
    b = ["no"] * len(a)
    r2 = m.agreement_generic(a, b, groups, n_boot=50)
    assert r2["kappa"]["estimate"] < 0.6 and not m.relevance_gate(r2["kappa"]["estimate"])
    assert m.relevance_gate(0.6) and m.relevance_gate(0.9) and not m.relevance_gate(0.599) and not m.relevance_gate(None) and not m.relevance_gate(float("nan"))


# ----------------------------------------------------------------------------- folds, nested CV, scoring
def _synthetic(n_dis=12, n_c=24, seed=0, signal=1.5):
    rng = np.random.default_rng(seed)
    Z, y, cls, members = {}, {}, {}, {}
    for k in range(n_dis):
        yy = (rng.random(n_c) < 0.25).astype(float)
        yy[0], yy[1] = 1.0, 0.0
        X = rng.normal(size=(n_c, 11))
        X[:, 8] += signal * yy                                  # an informative LLM column
        Z[f"d{k}"], y[f"d{k}"], cls[f"d{k}"] = X, yy, np.array([f"c{(i * 7 + k) % 60}" for i in range(n_c)])
        members[f"g{k}"] = [f"d{k}"]
    return Z, y, cls, members


def test_k_group_folds_never_split_a_group():
    members = {"g1": ["a", "b"], "g2": ["c"], "g3": ["d", "e"], "g4": ["f"], "g5": ["g"], "g6": ["h"]}
    folds = m.k_group_folds(members, 3, 5)
    assert sorted(d for f in folds for d in f) == sorted(d for v in members.values() for d in v)
    for g, ds in members.items():
        assert sum(1 for f in folds if set(ds) <= set(f)) == 1


def test_nested_cv_picks_a_grid_C_and_heads_use_training_folds_only():
    Z, y, cls, members = _synthetic()
    dev = sorted(Z)
    heads, chosen = m.nested_oof_heads(dev, members, Z, y, cls, m.MODEL_COLS["PL"])
    assert set(heads) == set(dev) and all(c in m.C_GRID for c in chosen) and len(chosen) == 5
    C, scores = m.select_C(dev, members, Z, y, cls, m.MODEL_COLS["PL"])
    assert C in m.C_GRID and set(scores) == {str(c) for c in m.C_GRID}
    fin = m.fit_final(dev, members, Z, y, cls, m.MODEL_COLS["PL"])
    assert fin["cols"] == m.MODEL_COLS["PL"] and len(fin["coef"]) == 6
    # the informative LLM column (index 8 -> position 3 in the PL columns) gets a positive weight
    assert fin["coef"][3] > 0


def test_pl_beats_po_when_the_llm_column_is_informative_and_not_otherwise():
    for signal, expect_positive in ((2.0, True), (0.0, False)):
        Z, y, cls, members = _synthetic(n_dis=20, signal=signal, seed=3)
        dev = sorted(Z)
        hp, _ = m.nested_oof_heads(dev, members, Z, y, cls, m.MODEL_COLS["PL"])
        ho, _ = m.nested_oof_heads(dev, members, Z, y, cls, m.MODEL_COLS["PO"])
        a = np.mean([m.es.auc_mw(v4.head_decision(hp[d], Z[d][:, m.MODEL_COLS["PL"]]), y[d]) for d in dev])
        b = np.mean([m.es.auc_mw(v4.head_decision(ho[d], Z[d][:, m.MODEL_COLS["PO"]]), y[d]) for d in dev])
        assert (a - b > 0.05) if expect_positive else (a - b < 0.05)


def test_score_models_pd_and_decision_rule_on_synthetic_heldout():
    Z, y, cls, members = _synthetic(n_dis=14, signal=2.5, seed=11)
    dis = sorted(Z)
    heads = {}
    for mod in ("PL", "PO", "B0L"):
        h, _ = m.nested_oof_heads(dis, members, Z, y, cls, m.MODEL_COLS[mod])
        heads[mod] = h
    samples = {d: {k: Z[d] for k in dl.sample_keys()} for d in dis}
    n_sub = {d: {k: 5 for k in dl.sample_keys()} for d in dis}
    b0 = v4.head_params(v4.fit_head("logit", np.vstack([Z[d][:, :8] for d in dis]), np.concatenate([y[d] for d in dis])), list(range(8)))
    mean = np.vstack([Z[d][:, :8] for d in dis]).mean(0)
    sc = m._scorers(heads, b0, mean)
    pop = {d: np.zeros(len(y[d])) for d in dis}
    res = m.score_models_pd(dis, samples, n_sub, sc, y, pop)
    tab = m.ab.table(res, members, y, m.MODELS, m.PAIRS)
    row = tab["fractions"]["0.01"]["paired"]["PL_minus_PO"]["auc"]
    assert row["estimate"] > 0 and len(row["ci95"]) == 2
    dec = m.decide(tab, True)
    assert set(dec["per_fraction"]) == {"0.03", "0.01"} and dec["hypothesis"] in ("supported", "partial, not supported", "not supported")
    assert dec["supported"] == (dec["per_fraction"]["0.03"]["lower_above_0"] and dec["per_fraction"]["0.01"]["lower_above_0"])
    assert m.decide(tab, False)["supported"] is False


def test_decision_rule_requires_both_fractions():
    def tab(lo03, lo01):
        mk = lambda lo: {"auc": {"estimate": 0.02, "ci95": [lo, 0.05]}}
        return {"fractions": {"0.03": {"paired": {"PL_minus_PO": mk(lo03)}}, "0.01": {"paired": {"PL_minus_PO": mk(lo01)}}}}
    assert m.decide(tab(0.001, 0.002), True)["hypothesis"] == "supported"
    assert m.decide(tab(0.001, -0.002), True)["hypothesis"] == "partial, not supported"
    assert m.decide(tab(-0.001, 0.0), True)["hypothesis"] == "not supported"


def test_coverage_counts_pairs_with_a_counted_sentence():
    L = np.zeros((6, 3))
    L[0, 0], L[4, 1] = 1.0, 1.0
    y = np.array([1, 1, 0, 0, 0, 0.0])
    c = m.coverage(L, y)
    assert c["all"] == pytest.approx(2 / 6) and c["positive"] == pytest.approx(0.5) and c["negative"] == pytest.approx(0.25)


def test_no_disease_specific_constant_in_the_module():
    src = inspect.getsource(m)
    for word in ("Wilms", "Fanconi", "BUB1B", "microcephaly"):
        assert word not in src


def test_join_period_repair_is_narrow_and_idempotent():
    src_with_period = "Sirolimus rescues growth in mice. Abstract text follows here."
    assert m.repair_join_period("Sirolimus rescues growth in mice..", src_with_period) == "Sirolimus rescues growth in mice."
    src_no_period = "Sirolimus rescues growth in mice <h4>Background</h4> Abstract text."
    assert m.repair_join_period("Sirolimus rescues growth in mice.", src_no_period) == "Sirolimus rescues growth in mice"
    assert m.repair_join_period("Sirolimus rescues growth in mice.", src_with_period) == "Sirolimus rescues growth in mice."      # already verbatim
    assert m.repair_join_period("Sirolimus cures all mice.", src_with_period) == "Sirolimus cures all mice."                       # not rescued: stays unverified
    assert m.repair_join_period("Anything.", None) == "Anything."
    idx = {"D": {"pairs": {"sirolimus": [{"i": 0, "pid": "MED:1", "sents": ["Sirolimus rescues growth in mice.."]}]}}}
    assert m.repair_pair_index(idx, {"MED:1": src_with_period}) == 1 and m.repair_pair_index(idx, {"MED:1": src_with_period}) == 0
    items = m.make_items("D", idx["D"]["pairs"])
    assert m.verify_sentences(items, {"MED:1": src_with_period}) == ["verified"]
