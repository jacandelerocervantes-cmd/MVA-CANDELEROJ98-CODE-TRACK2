"""Offline synthetic tests of engine_v4 (no network, no cache files, no model files)."""
import re
from pathlib import Path

import numpy as np

from track2 import engine_v4 as v4


def _diseases():
    ds = [{"name": f"v2_{i}", "origin": "v2", "group": f"g{i}"} for i in range(10)]
    ds += [{"name": "mate", "origin": "v4_new", "group": "g0"}]             # group-mate of a v2 disease
    ds += [{"name": f"n_{i:02d}", "origin": "v4_new", "group": f"f{i // 2}"} for i in range(20)]   # 10 free groups of 2
    return ds


def test_split_groups_never_straddle_and_are_deterministic():
    s = v4.make_split_v4(_diseases())
    gmap = {d["name"]: d["group"] for d in _diseases()}
    assert not set(s["dev"]) & set(s["heldout"])
    assert not {gmap[d] for d in s["dev"]} & {gmap[d] for d in s["heldout"]}
    assert "mate" in s["dev"]
    assert len(s["heldout"]) == 10 and len(s["dev"]) == 21            # free 20 -> 10 in DEV, 10 held out
    assert s == v4.make_split_v4(list(reversed(_diseases())))
    assert s["sha256_dev_ids"] == v4.hash_ids(s["dev"])


def test_strata_cut_points():
    assert [v4.stratum_of(h) for h in (0, 49, 50, 500, 501)] == ["L", "L", "M", "M", "H"]


def test_group_id_ignores_order():
    assert v4.group_id(["a", "b"]) == v4.group_id(["b", "a"]) != v4.group_id(["a"])


def test_candidate_rules_declared_changes():
    m = {"n_matched": 1, "rare_genetic_identity": True, "orphanet_catalogued": True, "mendelian": True, "genes": ["G1"]}
    u = {"d1": {"molecule_type": "Small molecule"}, "d2": {"molecule_type": "Enzyme"}}
    assert v4.evaluate_candidate_v4("x", {"d1"}, m, u)[0]                     # D1 and D2: one gene, one drug
    assert not v4.evaluate_candidate_v4("x", {"d2"}, m, u)[0]
    assert not v4.evaluate_candidate_v4("x", set(), m, u)[0]
    assert not v4.evaluate_candidate_v4("fanconi anemia", {"d1"}, m, u)[0]
    assert not v4.evaluate_candidate_v4("x", {"d1"}, dict(m, genes=[]), u)[0]


def test_conformal_mondrian_coverage_and_inf_rank():
    rng = np.random.default_rng(1)
    y = (rng.random(20000) < 0.05).astype(float)
    p = np.clip(0.05 + 0.5 * y * rng.random(20000) + 0.1 * rng.random(20000), 0, 1)
    q = v4.conformal_q(p, y, 0.1)
    st = v4.conformal_stats(p, y, q)
    assert st["positive_coverage"] >= 0.89 and st["negative_coverage"] >= 0.89
    # new data from the same law keeps the coverage (exchangeable)
    y2 = (rng.random(20000) < 0.05).astype(float)
    p2 = np.clip(0.05 + 0.5 * y2 * rng.random(20000) + 0.1 * rng.random(20000), 0, 1)
    assert v4.conformal_stats(p2, y2, q)["positive_coverage"] >= 0.87
    # tiny calibration set: rank above n -> threshold 2.0 (always included)
    assert v4.conformal_q(np.array([0.1, 0.2, 0.3, 0.4, 0.5]), np.array([1, 1, 1, 1, 1]), 0.1)["1"] == 2.0


def test_conformal_abstains_on_empty_and_two_label_sets():
    q = {"0": 0.2, "1": 0.7}                     # label 0 kept if p <= 0.2; label 1 kept if p >= 0.3
    p = np.array([0.1, 0.25, 0.9])
    h0, h1 = v4.conformal_sets(p, q)
    assert list(h0) == [True, False, False] and list(h1) == [False, False, True] and True
    q2 = {"0": 0.5, "1": 0.9}
    h0, h1 = v4.conformal_sets(np.array([0.3]), q2)
    assert h0[0] and h1[0]                        # two labels -> abstain
    st = v4.conformal_stats(np.array([0.3, 0.9, 0.0]), np.array([1.0, 1.0, 0.0]), q2)
    assert abs(st["abstention_rate"] - 1 / 3) < 1e-9 and st["two_label_set_rate"] > 0


def test_one_se_inclusion_rule():
    a = np.array([0.7, 0.8, 0.75, 0.72, 0.9])
    b = a - np.array([0.05, 0.06, 0.04, 0.05, 0.05])
    assert v4.one_se_inclusion(a, b)["included"]
    c = a + np.array([0.2, -0.2, 0.1, -0.1, 0.0])
    assert not v4.one_se_inclusion(c, a)["included"]


def test_tier_rule():
    cf = lambda v: [{"type": "Cases/families", "geo": "", "cls": "", "val": v}]  # noqa: E731
    assert v4.tier_from_records(cf(1.0)) == "1 known case"
    assert v4.tier_from_records(cf(4.0)) == "2-5"
    assert v4.tier_from_records(cf(12.0)) == ">5"
    assert v4.tier_from_records([{"type": "Point prevalence", "geo": "", "cls": "1-9 / 100 000", "val": None}]) == ">5"
    assert v4.tier_from_records([{"type": "Point prevalence", "geo": "", "cls": "<1 / 1 000 000", "val": None}]) == "unknown"
    assert v4.tier_from_records([]) == "unknown"
    assert "dose" not in " ".join(v4.DESIGN_BY_TIER.values()).lower().replace("dosing", "")


def test_prevalence_parser():
    xml = (b"<r><Disorder><OrphaCode>5</OrphaCode><PrevalenceList><Prevalence><PrevalenceType><Name>Cases/families</Name></PrevalenceType>"
           b"<ValMoy>3</ValMoy><PrevalenceClass><Name>Not yet documented</Name></PrevalenceClass></Prevalence></PrevalenceList></Disorder></r>")
    p = v4.parse_orphanet_prevalence(xml)
    assert v4.tier_from_records(p["ORPHA:5"]) == "2-5"


def test_drug_sentences_text_only_and_limit():
    papers = [{"pid": "a", "text": "GENE1 is studied. Aspirin helped patients. Aspirin again. Aspirin third. Aspirin fourth."},
              {"pid": "b", "text": "No gene here. Aspirin appears."}]
    from track2 import empirical_literature_control as elc
    dre = elc.drug_regex(["aspirin"])
    out = v4.drug_sentences(papers, ["GENE1"], dre, None, 3)
    assert out["aspirin"].count("Aspirin") == 3 and "No gene" not in out["aspirin"]


def test_head_roundtrip_logit_and_mlp():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(400, 4))
    y = (X[:, 0] + 0.5 * rng.normal(size=400) > 1).astype(float)
    for kind in ("logit", "mlp"):
        h = v4.fit_head(kind, X, y)
        s = v4.head_decision(h, X)
        p = v4.head_params(h, [0, 1, 2, 3])
        assert np.allclose(s, v4.head_decision(p, X))
        from track2 import engine_score as es
        assert es.auc_mw(s, y) > (0.8 if kind == 'logit' else 0.55)


def test_frozen_hash_detects_change():
    f = v4.freeze_models({"heads": {"B0": {"x": [1.0]}}}, {"final_model": "B0"})
    assert v4.thaw_v4(f)["final_model"] == "B0"
    f["model"]["final_model"] = "E"
    try:
        v4.thaw_v4(f)
        raise AssertionError("tampering not detected")
    except RuntimeError:
        pass


def test_group_bootstrap_evaluation_runs_and_groups_share_weight():
    rng = np.random.default_rng(3)
    n = 60
    ds = [f"d{i}" for i in range(8)]
    y, pop, proxy, pred = {}, {}, {}, {}
    for d in ds:
        yy = np.zeros(n)
        yy[rng.choice(n, 3, replace=False)] = 1
        y[d], pop[d], proxy[d] = yy, rng.random(n), rng.random(n)
        s = yy * 1.5 + rng.normal(size=n)
        pred[d] = {"score": s, "platt": v4.sigmoid(s - 3), "prev": np.full(n, 0.05)}
    members = {"g0": ["d0", "d1"], "g1": ["d2"], "g2": ["d3", "d4", "d5"], "g3": ["d6"], "g4": ["d7"]}
    ev, arr = v4.evaluate_groups(members, {"B0": pred, "E": pred}, y, pop, proxy, reference="B0")
    assert ev["n_groups"] == 5 and ev["n_diseases"] == 8 and len(arr["B0"]) == 5
    assert ev["models"]["B0"]["decision"]["verdict"]
    assert ev["models"]["E"]["auc_minus_B0"]["estimate"] == 0.0


def test_unvalidated_registry_complies_and_defaults_apply():
    assert v4.check_registry(v4.UNVALIDATED) == []
    assert v4.check_registry({"x": {"value": 1, "status": "UNVALIDATED"}})
    p = {"value": 0.2, "status": "UNVALIDATED", "conservative_default": 0.1}
    assert v4.resolve(p) == 0.1 and v4.resolve(p, True) == 0.2


def test_feature_function_has_no_gold_argument_and_no_disease_constants():
    import inspect
    assert not any("gold" in p for p in inspect.signature(v4.compute_disease_features).parameters)
    src = Path(v4.__file__).read_text(encoding="utf-8")
    assert not re.search(r"['\"](kidney Wilms|Alzheimer|hemophilia)", src)


def test_evidence_gap_shrinkage_weights():
    bench = {"diseases": [{"name": f"d{i}", "genes": [{"symbol": "G", "mechanism_class_general": "mitochondrial"}],
                           "identifiers": {"orphanet_xref": [], "id_matched": []}} for i in range(5)]}
    aucs = {f"d{i}": 0.6 + 0.05 * i for i in range(5)}
    gap = v4.evidence_gap(bench, aucs, {k: 3 for k in aucs}, 2705, {}, {k: "L" for k in aucs})
    r = gap["rows"]["d0"]
    assert r["tier"] == "unknown" and r["borrowing_pool_size"] == 4 and 0.0 <= r["weight_on_own_auc"] <= 1.0
