"""Offline, synthetic tests of track2.engine_v5 (no network, no project data)."""
import inspect
import json

import pytest

from track2 import engine_v5 as e5


def test_unvalidated_registry_complete():
    assert "k_records_per_pair" in e5.UNVALIDATED and e5.UNVALIDATED["k_records_per_pair"]["value"] == 15
    for k, v in e5.UNVALIDATED.items():
        assert v["status"] == "UNVALIDATED", k
        assert "conservative_default" in v, k


def test_guideline_hashes_frozen():
    assert e5.v4.sha256_file(e5.REL_GUIDELINE) == e5.GUIDELINE_SHA256["relevance_v2"]
    assert e5.v4.sha256_file(e5.DIR_GUIDELINE) == e5.GUIDELINE_SHA256["direction"]
    text = e5.PROTOCOL.read_text(encoding="utf-8")
    assert all(h in text for h in e5.GUIDELINE_SHA256.values())


def test_ordered_drug_terms_priority_and_filters():
    syn = {"foo": {"chembl_synonyms"}, "abc-123": {"ctgov_other_names"}, "brandx": {"ctgov_other_names"}, "ab": {"chembl_synonyms"},
           "drugname hcl": {"salt_stripped", "chembl_synonyms"}, 'bad"quote': {"chembl_synonyms"}, "zzzz": {"chembl_synonyms"}}
    t = e5.ordered_drug_terms("Drugname Sodium", syn)
    assert t[0] == "drugname sodium"
    assert t.index("drugname hcl") < t.index("zzzz") < t.index("brandx") < t.index("abc-123")   # source priority, then digit-free first
    assert "ab" not in t and 'bad"quote' not in t and "foo" not in t[:1] and len(t) == len(set(t))
    assert e5.ordered_drug_terms("abc", {}) == []                                              # name under 4 characters


def test_build_query_has_no_year_and_requires_abstract():
    q = e5.build_query(["aspirin", "acetylsalicylic acid"], ["Kawasaki disease", "mucocutaneous lymph node syndrome"])
    assert q.startswith('(TITLE_ABS:"aspirin" OR TITLE_ABS:"acetylsalicylic acid") AND (TITLE_ABS:"Kawasaki disease"')
    assert q.endswith("AND HAS_ABSTRACT:y") and "PUB_YEAR" not in q


def test_search_pair_cached_and_page_size(tmp_path):
    calls = []

    def getter(url, params):
        calls.append(params)
        return {"hitCount": 40, "nextCursorMark": "AoE=", "resultList": {"result": [
            {"source": "MED", "id": str(i), "title": f"T{i}", "abstractText": f"A{i} aspirin helps.", "pubYear": "2001"} for i in range(15)]}}
    a = e5.search_pair("q1", getter, tmp_path, k=15)
    b = e5.search_pair("q1", getter, tmp_path, k=15)
    assert a == b and len(calls) == 1 and calls[0]["pageSize"] == 15 and calls[0]["cursorMark"] == "*"
    assert a["hitCount"] == 40 and len(a["records"]) == 15


def test_record_sentences_cap_and_synonyms():
    rx = e5.elc.term_regex(["aspirin", "acetylsalicylic acid"])
    text = "Study of Aspirin. Aspirin lowered fever. Other thing. Acetylsalicylic acid was safe. Aspirin again. Aspirin last."
    s = e5.record_sentences(text, rx)
    assert len(s) == 3 and s[0] == "Study of Aspirin." and "Acetylsalicylic acid was safe." in s
    assert e5.record_sentences("Nothing here.", rx) == []


def test_thin_records_deterministic_and_fractions():
    assert e5.thin_records(15, 1.0, 3, 4, 0, 0) == list(range(15))
    a = e5.thin_records(15, 0.3, 3, 4, 1, 2)
    assert a == e5.thin_records(15, 0.3, 3, 4, 1, 2) and a == sorted(a)
    assert a != e5.thin_records(15, 0.3, 3, 4, 1, 3) or a != e5.thin_records(15, 0.3, 3, 5, 1, 2)
    kept = {f: sum(len(e5.thin_records(15, f, d, g, 2, r)) for d in range(40) for g in range(10) for r in range(2)) / 800 for f in (0.3, 0.03, 0.01)}
    assert 3.8 < kept[0.3] < 5.2 and 0.25 < kept[0.03] < 0.65 and kept[0.01] < 0.3 and kept[0.03] > kept[0.01]


def _index():
    return {"D1": {"drugA": {"records_all": ["MED:1", "MED:2", "MED:3"], "records": [{"pid": "MED:1", "sents": ["Drug A works."]}, {"pid": "MED:2", "sents": ["A again.", "A more."]}]},
                   "drugB": {"records_all": [], "records": []}},
            "D2": {"drugA": {"records_all": ["MED:9"], "records": [{"pid": "MED:9", "sents": ["Drug A in D2."]}]}}}


def test_items_ids_shards_and_fields():
    idx = _index()
    its = e5.make_items("D1", idx["D1"]) + e5.make_items("D2", idx["D2"])
    assert len(its) == 4
    items = e5.shuffle_and_id(its, "dev")
    e5.assert_item_fields(items)
    assert all(i["id"].startswith("dev5-") for i in items) and sorted(i["id"] for i in items) == [f"dev5-0000{n}" for n in range(4)]
    assert all(i["id"].startswith("ho5-") for i in e5.shuffle_and_id(its, "heldout"))
    assert items == e5.shuffle_and_id(its, "dev")                                             # deterministic
    assert not any(k in i for i in items for k in ("gold", "score", "g1", "top", "label"))
    big = [{"id": f"dev5-{n}", "drug": "x", "disease": "y", "sentence": "s", "pmid": "1"} for n in range(250)]
    sh = e5.lm.shard(big, e5.SHARD_SIZE)
    assert [len(s) for s in sh] == [100, 100, 50]
    assert e5.SHARD_SIZE == 100


def test_pair_stats_by_gold():
    idx = _index()
    gold = {("D1", "drugA"): {"g1": True, "g2": False, "g3": False}, ("D1", "drugB"): {"g1": False, "g2": False, "g3": False},
            ("D2", "drugA"): {"g1": False, "g2": False, "g3": True}}
    st = e5.pair_stats(idx, ["D1", "D2"], gold)
    assert st["all"]["n_pairs"] == 3 and st["all"]["n_pairs_with_ge1_direct_record"] == 2 and st["all"]["n_direct_records"] == 4
    assert st["gold1"]["n_pairs"] == 1 and st["gold1"]["n_verified_sentences"] == 3 and st["not_gold1"]["n_pairs"] == 2
    assert st["no_gold"]["n_pairs_with_ge1_direct_record"] == 0 and st["any_gold"]["n_pairs"] == 2
    assert e5.pair_stats(idx, ["D1"])["all"]["n_pairs_scarce_0_to_2_records"] == 1
    assert list(e5.pair_stats(idx, ["D1"])) == ["all"]                                         # gold keys only when asked


def test_verification_and_join_repair_reuse():
    texts = {"MED:1": "A title. Drug A works."}
    src = texts["MED:1"]
    assert e5.lm.verify_sentences([{"pmid": "1", "sentence": "Drug A works."}, {"pmid": "1", "sentence": "Drug A fails."}, {"pmid": "9", "sentence": "x"}], texts) \
        == ["verified", "not_found", "not_checkable"]
    assert e5.lm.repair_join_period("A title..", "A title. Drug A works.") in ("A title..", "A title.")
    assert e5.lm.verify_sentences([{"pmid": "1", "sentence": e5.lm.repair_join_period("A title.", src)}], texts) == ["verified"]


def test_no_gold_argument_in_feature_side_functions():
    for fn in (e5.record_sentences, e5.thin_records, e5.make_items, e5.shuffle_and_id, e5.build_query, e5.ordered_drug_terms, e5.search_pair):
        params = set(inspect.signature(fn).parameters)
        assert not params & {"gold", "gold_of", "y", "labels", "score"}, fn.__name__


def test_build_pair_index_offline(tmp_path, monkeypatch):
    monkeypatch.setattr(e5, "QUERY_CACHE", tmp_path)
    monkeypatch.setattr(e5, "search_pair", lambda q, g=None, c=tmp_path, k=15: {"hitCount": 2, "records": [
        {"pid": "MED:1", "title": "Aspirin in KD", "abstract": "Aspirin reduced aneurysms. Nothing else.", "year": "2000", "language": "eng"},
        {"pid": "MED:2", "title": "No drug here", "abstract": "Kawasaki disease only.", "year": "2001", "language": "eng"}]})
    doc = {"diseases": {"KD": {"split": "dev", "pairs": {"aspirin": {"query": "q", "pids": ["MED:1", "MED:2"]}, "none": {"query": None, "pids": []}}}}}
    idx = e5.build_pair_index(doc, {"aspirin": ["aspirin", "asa tablet"], "none": []})
    assert idx["KD"]["aspirin"]["records_all"] == ["MED:1", "MED:2"]
    assert idx["KD"]["aspirin"]["records"] == [{"pid": "MED:1", "sents": ["Aspirin in KD.", "Aspirin reduced aneurysms."]}]
    assert idx["KD"]["none"] == {"records_all": [], "records": []}
    json.dumps(idx)
