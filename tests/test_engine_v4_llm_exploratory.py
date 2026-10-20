"""Offline synthetic tests of engine_v4_llm_exploratory (no network, no cache files, no LLM labels, no model files)."""
import inspect
from pathlib import Path

import numpy as np

from track2 import engine_v4 as v4
from track2 import engine_v4_llm_exploratory as x


def _items():
    base = {"disease": "D", "drug": "aspirin", "pmid": "1"}
    return [dict(base, id="a", sentence="Aspirin improved survival in patients."),
            dict(base, id="b", sentence="Aspirin was toxic and increased mortality."),
            dict(base, id="c", sentence="Aspirin improved outcomes."),
            dict(base, id="d", sentence="Aspirin was studied.")]


def _ann(rel, dr):
    return {i: {"relevance": r, "direction": d} for i, r, d in zip("abcd", rel, dr)}


def test_registry_complete():
    assert v4.check_registry(x.UNVALIDATED) == []
    for k, p in x.UNVALIDATED.items():
        assert p["status"] == "UNVALIDATED" and "conservative_default" in p, k


def test_model_columns_prior_plus_one_feature():
    assert x.MODEL_COLS["PO"] == [4, 6, 7]
    for f in x.FEATURES:
        assert x.MODEL_COLS[x.MODEL_NAME[f]] == [4, 6, 7, x.FEATURE_COLS[f]]
    assert x.REPORT_FRACTIONS == (1.0, 0.1, 0.03, 0.01)


def test_no_gold_argument_and_no_heldout_access():
    for fn in (x.exploratory_label_table, x.exploratory_counts, x.to_model_block):
        bad = [p for p in inspect.signature(fn).parameters if "gold" in p.lower() or p.lower() in ("y", "score")]
        assert not bad, (fn.__name__, bad)
    src = Path(x.__file__).read_text(encoding="utf-8")
    code = "\n".join(l for l in src.splitlines() if not l.lstrip().startswith("#"))
    for forbidden in ("_heldout.json", "run_heldout(", "x.LOCK", "lm.LOCK", "split_of", '"heldout"]'):
        assert forbidden not in code.split('"""', 2)[-1], forbidden
    assert "LOCK" not in code


def test_feature_definitions_on_labels():
    items = _items()
    # a: both yes + benefit ; b: harm ; c: one yes only ; d: unclear
    a = _ann(["yes", "yes", "yes", "no"], ["benefit", "harm", "benefit", "unclear"])
    b = _ann(["yes", "yes", "no", "no"], ["benefit", "harm", "benefit", "unclear"])
    t = x.exploratory_label_table(items, a, b)
    key = lambda s: ("D", "aspirin", "1", s)
    assert t[key(items[0]["sentence"])] == {"F_dir": True, "F_rel_strict": True, "F_rel_any": True}
    assert t[key(items[1]["sentence"])] == {"F_dir": False, "F_rel_strict": False, "F_rel_any": False}      # harm is never counted
    assert t[key(items[2]["sentence"])] == {"F_dir": True, "F_rel_strict": False, "F_rel_any": True}        # either-yes versus both-yes
    assert t[key(items[3]["sentence"])]["F_dir"] is False                                                   # unclear direction


def test_direction_uses_two_of_three_vote():
    items = [dict(_items()[2], id="a")]            # rule extractor says benefit for this sentence
    t = x.exploratory_label_table(items, {"a": {"relevance": "no", "direction": "harm"}}, {"a": {"relevance": "no", "direction": "benefit"}})
    assert t[("D", "aspirin", "1", items[0]["sentence"])]["F_dir"] in (True, False)   # defined (no exception) whichever the rule says
    t2 = x.exploratory_label_table(items, {"a": {"relevance": "no", "direction": "benefit"}}, {"a": {"relevance": "no", "direction": "benefit"}})
    assert t2[("D", "aspirin", "1", items[0]["sentence"])]["F_dir"] is True       # both annotators agree


def test_counts_cap_papers_and_use_kept_only():
    sents = [f"Aspirin improved outcome {i}." for i in range(10)]
    pair_index = {"aspirin": [{"i": i, "pid": f"MED:{i}", "sents": [sents[2 * i], sents[2 * i + 1]]} for i in range(5)]}
    table = {("D", "aspirin", str(i), s): {"F_dir": True, "F_rel_strict": i % 2 == 0, "F_rel_any": True} for i in range(5) for s in pair_index["aspirin"][i]["sents"]}
    C = x.exploratory_counts("D", pair_index, {0, 1, 2, 3, 4}, table, {"aspirin": 0}, 2)
    assert C[0].tolist() == [6, 4, 6] and C[1].tolist() == [0, 0, 0]                   # first 3 papers x 2 sentences
    C2 = x.exploratory_counts("D", pair_index, {3, 4}, table, {"aspirin": 0}, 2)
    assert C2[0].tolist() == [4, 2, 4]                                                  # papers 3 (odd) and 4 (even)
    assert x.exploratory_counts("D", pair_index, set(), table, {"aspirin": 0}, 2).sum() == 0
    assert np.allclose(x.to_model_block(C), np.log1p(C))


def test_coverage_summary():
    r1 = np.array([[1, 0, 2], [0, 0, 0], [3, 1, 1], [0, 0, 0]], dtype=float)
    r2 = np.array([[0, 0, 0], [0, 0, 0], [3, 1, 1], [1, 0, 1]], dtype=float)
    y = np.array([1, 0, 1, 0])
    c = x.coverage_summary([r1, r2], y)
    assert c["F_dir"]["nonzero_pairs_mean"] == 2.0 and c["F_dir"]["nonzero_pairs_min"] == 2 and c["F_dir"]["nonzero_pairs_max"] == 2
    assert c["F_rel_any"]["nonzero_pairs_min"] == 2 and c["F_rel_strict"]["nonzero_pairs_max"] == 1
    assert c["F_dir"]["nonzero_positive_mean"] == 1.5 and c["F_dir"]["nonzero_negative_mean"] == 0.5
    assert abs(c["F_rel_strict"]["nonzero_share"] - 0.25) < 1e-12
    assert c["F_rel_any"]["n_pairs"] == 4


def test_scorers_use_their_columns():
    class H(dict):
        pass
    heads = {m: {"d": v4.fit_head("logit", np.random.default_rng(0).normal(size=(60, len(c))), (np.arange(60) % 2).astype(float), 1.0)}
             for m, c in x.MODEL_COLS.items()}
    sc = x.make_scorers(heads)
    Z = np.random.default_rng(1).normal(size=(5, 11))
    Z2 = Z.copy()
    Z2[:, 8] += 5                                                  # only the F_dir column changes
    assert np.allclose(sc["PO"]("d", Z), sc["PO"]("d", Z2))
    assert not np.allclose(sc["PO_F_dir"]("d", Z), sc["PO_F_dir"]("d", Z2))
    assert np.allclose(sc["PO_F_rel_any"]("d", Z), sc["PO_F_rel_any"]("d", Z2))
