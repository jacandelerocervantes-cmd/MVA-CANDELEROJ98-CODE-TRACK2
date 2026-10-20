"""Offline synthetic tests of engine_v4_llm_audit (no network, no cache files, no labels, no model files)."""
import inspect
import json

from track2 import empirical_literature_control as elc
from track2 import engine_v4 as v4
from track2 import engine_v4_llm_audit as a


def test_registry_complete():
    assert v4.check_registry(a.UNVALIDATED) == []
    for k, p in a.UNVALIDATED.items():
        assert p["status"] == "UNVALIDATED" and "conservative_default" in p, k


def test_no_heldout_access_in_source():
    src = inspect.getsource(a)
    for bad in ("run_heldout", "heldout_v4_llm_scored", "opus_heldout", "sonnet_heldout", "ho_shard"):
        assert bad not in src, bad
    assert 'v["split"] == "dev"' in src and '_load_split_items("dev")' in src


def test_salt_and_plain_inn():
    assert a.salt_stripped("CLEMIZOLE HYDROCHLORIDE") == "clemizole"
    assert a.salt_stripped("VALPROIC ACID") is None
    assert a.salt_stripped("ASPIRIN") is None
    assert a.is_plain_inn("ASPIRIN") and not a.is_plain_inn("FLUDARABINE PHOSPHATE") and not a.is_plain_inn("X (Y)")
    assert a.name_variants("5-fluorouracil") == {"5 fluorouracil"}


def test_synonyms_sources(tmp_path):
    ct = tmp_path / "c.json"
    ct.write_text(json.dumps({"studies": [{"protocolSection": {"armsInterventionsModule": {"interventions": [
        {"name": "Aspirin", "otherNames": ["ASA tablets", "ab"]}]}}}]}), encoding="utf-8")
    ch = tmp_path / "m.json"
    ch.write_text(json.dumps({"pref_name": "ASPIRIN", "molecule_synonyms": [{"molecule_synonym": "Bayer Aspirin"}]}), encoding="utf-8")
    syn = a.build_synonyms(["ASPIRIN", "CLEMIZOLE HYDROCHLORIDE"], [ct], [ch])
    assert set(syn["aspirin"]) == {"asa tablets", "bayer aspirin"}
    assert syn["clemizole hydrochloride"]["clemizole"] == {"salt_stripped"}


def test_crosstab_and_kappa():
    r = a.crosstab(["yes", "no", "no", "yes"], ["yes", "no", "yes", "yes"], a.VOC_REL)
    assert r["n"] == 4 and r["agreement"] == 0.75 and r["rows_first_columns_second"]["no"]["yes"] == 1
    assert abs(r["kappa"] - 0.5) < 1e-9


def test_sample_deterministic():
    pairs = [("d", f"x{i}") for i in range(40)]
    s1, s2 = a.sample_pairs(pairs, 10), a.sample_pairs(list(reversed(pairs)), 10)
    assert s1 == s2 and len(s1) == 10 and a.sample_pairs(pairs, 100) == sorted(pairs)


def test_queries():
    qd = a.query_disease_terms("Alzheimer disease", ["Alzheimer's disease", "Alzheimer disease 5", "EOFAD", "late onset, x"])
    assert qd[0] == "Alzheimer disease" and "Alzheimer disease 5" not in qd and "late onset, x" not in qd
    q = a.direct_query(["aspirin"], qd)
    assert 'TITLE_ABS:"aspirin"' in q and "PUB_YEAR:[1800 TO 2024]" in q and "HAS_ABSTRACT:y" in q


def _dt():
    papers = [{"pid": "MED:1", "text": "Aspirin helped patients with disease Y. Aspirin was well tolerated. Aspirin dose was low. BRCA1 was mutated."},
              {"pid": "MED:2", "text": "Acetylsalicylic acid improved disease Y."},
              {"pid": "MED:3", "text": "Aspirin treated disease Y. BRCA1 noted. Aspirin again."},
              {"pid": "MED:4", "text": "Aspirin in disease Y. BRCA1 noted."}]
    dre = elc.drug_regex(["aspirin", "ibuprofen"])
    return a.DiseaseText(papers, ["BRCA1"], dre, None, elc.term_regex(["disease Y"]), elc.term_regex(["disease Y", "Y syndrome"]))


def test_extract_caps_and_gate():
    dt = _dt()
    m = a.match_prod(dt, "aspirin")
    base = a.extract_pair(dt, range(4), m, 3, 2, True)
    assert [i for i, _ in base] == [0, 2, 3] and [len(j) for _, j in base] == [2, 2, 1]
    assert [len(j) for _, j in a.extract_pair(dt, range(4), m, 3, None, True)] == [3, 2, 1]
    assert [i for i, _ in a.extract_pair(dt, range(4), m, None, None, False)] == [0, 2, 3]
    assert [i for i, _ in a.extract_pair(dt, [1, 2, 3], m, 1, 2, True)] == [2]
    ms = a.match_syn(dt, "aspirin", ["acetylsalicylic acid"])
    assert [i for i, _ in a.extract_pair(dt, range(4), ms, None, None, False)] == [0, 1, 2, 3]
    assert [i for i, _ in a.extract_pair(dt, range(4), ms, None, None, True)] == [0, 2, 3]   # paper 1 has no gene symbol


def test_records_and_metrics():
    dt = _dt()
    ext = a.extract_pair(dt, range(4), a.match_prod(dt, "aspirin"), None, None, True)
    recs = a.records_of(dt, ext)
    m = a.pair_metrics(recs + recs, "aspirin")          # duplicates are removed
    assert m["sent_drug"] == len(recs) == 6 and m["sent_drug_disease_name"] == 3 and m["papers_drug"] == 3 and m["papers_drug_and_disease"] == 3
    s = a.summarise_variant([m, {k: 0 for k in m}])
    assert s["sent_drug"]["total"] == 6 and s["sent_drug"]["pairs_with_ge1"] == 1 and s["n_pairs"] == 2


def test_stage_flags_chain_and_thinning():
    entries = [{"i": 0, "pid": "MED:1", "sents": ["s1", "s2"]}, {"i": 2, "pid": "MED:3", "sents": ["s3"]}]
    table = {("D", "x", "1", "s1"): a.label_flags("yes", "no", "benefit"), ("D", "x", "1", "s2"): a.label_flags("yes", "yes", "benefit"),
             ("D", "x", "3", "s3"): a.label_flags("no", "no", None)}
    full = a.pair_stage_flags(entries, {0, 1, 2}, True, True, "D", "x", table)
    assert all(full[s] for s in a.STAGES) and full["S6_yes_both_and_benefit"]
    thin = a.pair_stage_flags(entries, {2}, True, True, "D", "x", table)
    assert thin["S2_sentence_extracted"] and thin["S3_verified_sentence"] and not thin["S4_relevance_yes_either"]
    none = a.pair_stage_flags(entries, {1}, False, False, "D", "x", table)
    assert not none["S1a_drug_named_in_a_kept_paper"] and not none["S2_sentence_extracted"]
    only_either = a.label_flags("yes", "unclear", "benefit")
    assert only_either["S4_relevance_yes_either"] and not only_either["S5_relevance_yes_both"] and only_either["S6b_yes_either_and_benefit"]
    assert not a.label_flags("yes", "yes", "harm")["S6_yes_both_and_benefit"]


def test_largest_losses():
    c = dict(zip(a.STAGES, [100, 60, 50, 50, 49, 10, 4, 3]))
    out = a.largest_losses(c)
    assert out["largest_absolute_loss"]["to"] == a.STAGES[1] and out["lowest_retention"]["to"] == a.STAGES[5]


def test_recall_describe_language():
    assert a.recall(["a", "b", "c", "d"], {"a", "c"}) == 0.5 and a.recall([], {"a"}) is None
    assert a.describe([1, 2, 3])["median"] == 2 and a.describe([])["n"] == 0
    ls = a.language_summary(["eng", "eng", "fre", None])
    assert ls["non_english"] == 2 and ls["counts"]["unknown"] == 1


def test_epmc_cache_and_gene_paging(tmp_path):
    calls = []

    def getter(url, params):
        calls.append(params["cursorMark"])
        first = params["cursorMark"] == "*"
        n, base = (100, 0) if first else (20, 100)
        res = [{"source": "MED", "id": str(base + i), "pubYear": "2000", "language": "eng" if i % 4 else "fre", "title": "t", "abstractText": "abs"} for i in range(n)]
        return {"hitCount": 120, "nextCursorMark": "next" if first else params["cursorMark"], "resultList": {"result": res}}
    g = a.fetch_gene_records(["BRCA1"], 300, getter, tmp_path)
    assert len(g["records"]) == 120 and g["hitCount"] == 120 and calls == ["*", "next"]
    n_calls = len(calls)
    a.fetch_gene_records(["BRCA1"], 300, getter, tmp_path)           # fully cached
    assert len(calls) == n_calls
    assert a.language_summary([r["language"] for r in g["records"]])["non_english"] == 30


def test_q3_pool_and_sample():
    rec = {}
    gold = {("D", "x"): {"g1": True}, ("D", "y"): {"g1": True}, ("D", "z"): {"g1": False}}
    for i, (drug, o, s) in enumerate([("x", "no", "unclear"), ("x", "unclear", "no"), ("y", "yes", "no"), ("z", "no", "no")]):
        rec[f"i{i}"] = {"disease": "D", "drug": drug, "sentence": "s", "pmid": "1", "o_rel": o, "s_rel": s, "o_dir": "unclear", "s_dir": "unclear", "cons_dir": None}
    pool, by_pair = a.q3a_pool(rec, gold)
    assert pool == [("D", "x")] and len(by_pair) == 2
    s = a.q3a_sample(rec, gold, 5)
    assert s["n_pairs_no_yes_from_either_annotator"] == 1 and s["ids"][0] in ("i0", "i1")
    b = a.q3b(rec, gold)
    assert b["gold_positive_pairs"]["n_items"] == 3 and b["other_pairs"]["n_items"] == 1


def test_pilot_disagreement_counts():
    ids = [f"p{i}" for i in range(4)]
    rec = {i: {"o_rel": o, "s_rel": s, "o_dir": "unclear", "s_dir": "unclear"} for i, o, s in zip(ids, ["yes", "no", "no", "yes"], ["no", "yes", "no", "yes"])}
    pilot = {"opus": {i: {"relevance": r, "direction": "unclear"} for i, r in zip(ids, ["yes", "yes", "no", "no"])},
             "sonnet": {i: {"relevance": r, "direction": "unclear"} for i, r in zip(ids, ["yes", "yes", "no", "yes"])}}
    out = a.q3d(rec, pilot)
    d = out["relevance_disagreement"]
    assert d["v1_items"] == 2 and d["removed_by_v2"] == 2 and d["newly_created_by_v2"] == 1 and d["kept_in_v2"] == 0
    assert out["yes_counts"]["v2"]["yes_both"] == 2


def test_feature_functions_have_no_gold_argument():
    for fn in (a.pair_stage_flags, a.extract_pair, a.label_flags):
        assert not any("gold" in p for p in inspect.signature(fn).parameters)
