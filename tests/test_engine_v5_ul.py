"""Offline, synthetic tests of track2.engine_v5_ultralight (no network, no project data)."""
import numpy as np

from track2 import engine_v5 as e5
from track2 import engine_v5_ultralight as ul


def test_surviving_matches_v5_rng():
    pids = [f"MED:{i}" for i in range(15)]
    kept = e5.thin_records(15, 0.03, 4, 7, 3, 2)
    assert ul.surviving_pids(pids, 4, 7, 3, 2, 0.03) == [pids[i] for i in kept]


def test_select_items_by_union_and_ids_unchanged():
    items = [{"id": "ho5-00002", "drug": "a", "disease": "d", "sentence": "s2", "pmid": "2"},
             {"id": "ho5-00001", "drug": "a", "disease": "d", "sentence": "s1", "pmid": "1"},
             {"id": "ho5-00003", "drug": "b", "disease": "d", "sentence": "s3", "pmid": "1"}]
    sel = ul.select_items(items, {("d", "a", "1"), ("d", "a", "2")})
    assert [i["id"] for i in sel] == ["ho5-00001", "ho5-00002"]
    assert all(tuple(ul.strip_item(i)) == ("id", "drug", "disease", "sentence", "pmid") for i in sel)


def test_opus_shards_seeded_30_percent():
    a = ul.opus_shards(13)
    assert a == ul.opus_shards(13) and len(a) == 4 and set(a) <= set(range(1, 14))
    assert len(ul.opus_shards(1)) == 1


def test_kappa_gate_perfect_and_poor():
    s = {str(i): {"relevance": "yes" if i % 2 else "no"} for i in range(10)}
    assert ul.relevance_kappa(s, dict(s))["gate_passed"]
    o = {k: {"relevance": "no" if v["relevance"] == "yes" else "yes"} for k, v in s.items()}
    assert not ul.relevance_kappa(s, o)["gate_passed"]
    assert not ul.relevance_kappa(s, {})["gate_passed"]


def test_item_is_benefit_rules():
    s = {"a": {"relevance": "yes", "direction": "benefit"}, "b": {"relevance": "yes", "direction": "harm"}, "c": {"relevance": "yes", "direction": "benefit"}}
    o = {"c": {"relevance": "no", "direction": "benefit"}}
    assert ul.item_is_benefit(s, o, "a") and not ul.item_is_benefit(s, o, "b") and not ul.item_is_benefit(s, o, "c")


def test_build_B_uses_only_surviving_records():
    b = ul.build_B({"g1": ["p1", "p2"], "g2": ["p3"]}, {"g1": {"p2"}, "g2": set()}, {("g1", "p1"), ("g2", "p3")}, "d")
    assert b == {"g1": 0, "g2": 0}
    assert ul.build_B({"g1": ["p1"]}, {"g1": {"p1"}}, {("g1", "p1")}, "d") == {"g1": 1}


def test_delta_and_combine_default_no_effect():
    po = np.array([0.0, 1.0, 2.0, 3.0])
    d = ul.delta_for(po)
    assert abs(d - 0.5 * np.std(po, ddof=1)) < 1e-12 and ul.delta_for(po[:1]) == 0.0
    out = ul.combine_score(po, [1, 3], [1, 0], d)
    assert out[1] == po[1] + d and out[3] == po[3] and out[0] == 0.0
    assert np.array_equal(ul.combine_score(po, [1], [1], 0.0), po)
    assert ul.DELTA_SD_FRACTION == 0.5 and ul.KAPPA_GATE == 0.6
