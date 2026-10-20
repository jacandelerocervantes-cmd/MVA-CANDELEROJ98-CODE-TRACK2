"""Tests of the multi-evidence integration (synthetic inputs, offline, fast)."""

from __future__ import annotations

import json

import numpy as np

from track2 import multi_evidence as ME


def entry(direction, level, pmid, consequence=None, finding=""):
    return {"direction": direction, "level": level, "pmid": pmid, "consequence": consequence, "finding": finding, "system": "", "sentence": "", "section": "findings", "use": None}


def test_lit_grade_counts_papers_once_and_signs():
    assert ME.lit_grade([])[0] is None
    g, info = ME.lit_grade([entry("may_help", 1, "1"), entry("may_help", 0, "1")])  # same paper: one pair at its highest level
    assert info["n_help_papers"] == 1 and info["best_help_level"] == 1 and g > 0.5
    g2, info2 = ME.lit_grade([entry("may_help", 0, "1"), entry("may_harm", 3, "2")])
    assert g2 < 0.5 and info2["n_harm_papers"] == 1
    g3, _ = ME.lit_grade([entry("unclear", 3, "1"), entry("no_effect", 1, "2")])
    assert g3 == 0.5


def test_mechanism_fit_uses_consequence_confidence_and_cap():
    rules = [{"match_any": ["mtor"], "consequences": ["4"], "context": "x"}, {"match_any": ["tubulin"], "consequences": ["5"], "only_for_list": "y"}]
    base, cons = ME.mech_base([], "mtor inhibitor", rules)
    assert base == ME.CONSEQUENCE_CONFIDENCE["4"] * 0.5 and cons == {"4": 0.5}
    assert ME.mech_base([], "tubulin binder", rules)[0] is None  # rule restricted to a list is not applied
    b2, _ = ME.mech_base([entry("may_help", 1, "1", "tumour")], "", rules)
    assert b2 == 1.0
    b3, _ = ME.mech_base([entry("may_harm", 1, "1", "tumour")], "", rules)
    assert b3 == 0.0
    block = {"primary_gene": "BUB1B", "genes": [
        {"gene": "BUB1B", "tier_weight": 0.0, "bears_on": []},
        {"gene": "IGBP1", "tier_weight": 1.0, "bears_on": [{"mechanism_regex": "mtor|akt", "relevance": 1.0}, {"literature_regex": "ampk", "relevance": 0.5}]},
        {"gene": "CENPE", "tier_weight": 0.5, "bears_on": [{"mechanism_regex": "tubulin", "relevance": 1.0}]}]}
    t, genes = ME.track1_bonus([], "mtor inhibitor", block)
    assert t == 1.0 and genes == ["IGBP1"]
    t2, _ = ME.track1_bonus([entry("may_help", 1, "1", None, "activates AMPK")], "", block)
    assert t2 == 0.5
    assert ME.track1_bonus([], "other", block)[0] == 0.0 and ME.track1_bonus([], "mtor", None)[0] == 0.0
    # the cap raises M by at most cap, never lowers it, and never exceeds 1
    assert abs(ME.mech_grade(0.3, 1.0, 0.1) - 0.4) < 1e-12
    assert ME.mech_grade(0.3, 0.0, 0.1) == 0.3 and ME.mech_grade(0.95, 1.0, 0.1) == 1.0 and ME.mech_grade(None, 0.0, 0.1) is None
    assert ME.mech_grade(0.3, 1.0, 0.0) == 0.3


def test_scores_renormalise_and_penalise_missing():
    G = np.array([[1.0, 1.0], [1.0, np.nan], [np.nan, np.nan]])
    W = np.array([[0.5, 0.5]])
    S = ME.scores(G, W)[0]
    assert abs(S[0] - 1.0) < 1e-12
    assert abs(S[1] - (1.0 * (1 - ME.MISSING_PENALTY * 0.5))) < 1e-12  # mean of the available, times the penalty for half the weight missing
    assert S[2] == 0.0


def test_smaa_acceptability_and_decision():
    rng = np.random.default_rng(1)
    # drug 0 dominates, drug 1 is second, the rest are weak and missing
    G = np.array([[0.9] * 4, [0.6] * 4] + [[0.1, np.nan, np.nan, 0.1] for _ in range(6)])
    W = ME.draw_weights(rng, 2000, 4)
    assert np.allclose(W.sum(1), 1.0)
    res = ME.smaa(G, W)
    assert res["p1"][0] == 1.0 and res["p3"][1] == 1.0
    assert np.all(res["p1"] <= res["p3"]) and np.all(res["p3"] <= res["p5"]) and np.all(res["p5"] <= res["p10"] + 1e-12)
    assert abs(res["p1"].sum() - 1.0) < 1e-9
    assert np.allclose(res["central"][0].sum(), 1.0) and res["confidence"][0] == 1.0 and res["confidence"][2] < 1.0
    names = [f"D{i}" for i in range(8)]
    prop, prim = ME.decide(names, res)
    assert prop[:2] == ["D0", "D1"] and prim[:2] == ["D0", "D1"] and len(prim) <= ME.N_PRIMARY
    assert all(res["p5"][names.index(d)] >= ME.P5_CUT for d in prop)


def test_band_state_and_directions():
    rec = {"L": {"state": "known", "intervals": [[13.0, 18.0]]}, "H_intervals_years": []}
    assert ME.band_state(rec, 0.0, 5.0) == "not_labelled" and ME.band_state(rec, 12.0, 18.0) == "partly"
    rec2 = {"L": {"state": "known", "intervals": [[1.0, 18.0]]}, "H_intervals_years": [[1.0, 2.0]]}
    assert ME.band_state(rec2, 5.0, 12.0) == "labelled" and ME.band_state(rec2, 0.0, 5.0) == "partly"
    assert ME.band_state({"L": {"state": "unknown", "intervals": []}}, 0.0, 5.0) == "verify"
    assert ME.engine_band("partly (from 13 y)") == "partly" and ME.engine_band("age_not_stated") == "labelled" and ME.engine_band(None) is None
    assert ME.direction("LIT", 0.7) == "supportive" and ME.direction("D", 0.0) == "contrary" and ME.direction("A", 0.3) == "neutral" and ME.direction("A", None) == "missing"


def synthetic_case():
    drugs = ["AAA", "BBB", "CCC", "DDD", "EEE", "FFF", "GGG", "HHH"]
    by_name = {d: {"name": d, "mechanisms": [{"mechanism": "mtor inhibitor" if d == "AAA" else "other", "action": ""}], "black_box_warning": 0} for d in drugs}
    entries = {"AAA": {"entries": [entry("may_help", 2, "1", "growth_muscle"), entry("may_help", 1, "2", "growth_muscle")]},
               "BBB": {"entries": [entry("may_help", 1, "3", "tumour")]},
               "CCC": {"entries": [entry("selectively_kills_aneuploid_cells", 1, "4", "tumour")]},
               "DDD": {"entries": [entry("may_harm", 3, "5", "tumour")]}}
    lines = {"A": {d: 0.9 - 0.1 * i for i, d in enumerate(drugs)}, "B1": {"AAA": 0.8, "BBB": 0.6}, "B2": {}, "C": {}, "D": {"AAA": 1.0, "DDD": 0.0},
             "E": {"AAA": 0.9}, "EP": {}, "ENG": {"BBB": 0.95}}
    return {"names": {}, "by_name": by_name, "entries": entries, "rules": [], "gene_block": None, "cap": 0.0, "lines": lines, "tops": {}, "caution_names": {"DDD"},
            "hard_class": {"EEE"}, "hazard": {"CCC"}, "proposal_safety": {}, "pool": drugs, "not_checked": [], "not_eligible": [], "not_small_or_unknown": [], "n_A": 8}


def test_analyse_excludes_hazards_and_reports_stability():
    case = synthetic_case()
    age = {"AAA": {"0-4": "not_labelled", "5-11": "partly", "12-17": "labelled"}, "BBB": {"0-4": "verify", "5-11": "not_labelled", "12-17": "not_labelled"}}
    res = ME.analyse(case, age, n_draws=1500, seed=7, mc_reps=5, mc_draws=300)
    assert set(res["excluded"]) == {"CCC", "EEE"} and "CCC" not in res["drugs"] and "EEE" not in res["drugs"]
    assert res["drugs"]["DDD"]["grades"]["S"] == 0.0 and res["drugs"]["DDD"]["directions"]["LIT"] == "contrary"
    assert "AAA" in res["proposed"] and res["drugs"]["AAA"]["p_top5"] >= 0.5 and res["drugs"]["AAA"]["age_bands"]["12-17"] == "labelled"
    assert res["drugs"]["HHH"]["age_bands"]["0-4"] == "not_assessed"
    assert "level 2 literature" in res["drugs"]["AAA"]["evidence_tier"]
    bands = res["age_bands"]
    assert bands["12-17"]["filter_pass"] == ["AAA"] and bands["0-4"]["filter_pass"] == [] and set(bands["unknown_age_union"]["filter_pass"]) == {"AAA"}
    assert all(d in res["proposed"] for d in bands["5-11"]["proposed"])  # headline band shortlist is a subset of the full-pool shortlist
    s = res["sensitivity"]
    assert set(s["leave_one_criterion_out"]) == set(ME.CRITERIA) and set(s["label_error_monte_carlo"]) == {"10", "20", "30"}
    for v in s["label_error_monte_carlo"].values():
        assert 0 <= v["expected_baseline_retained"] <= v["n_baseline"]
    assert len(s["equal_weights"]["top5"]) == 5 and "track1_gene_block_off" in s and "literature_only_weights" in s
    json.dumps(res, sort_keys=True)  # serialisable


def test_deterministic_and_markdown_wording():
    case = synthetic_case()
    a = ME.analyse(case, {}, n_draws=500, seed=3, mc_reps=2, mc_draws=100)
    b = ME.analyse(case, {}, n_draws=500, seed=3, mc_reps=2, mc_draws=100)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    result = {"label": ME.__doc__.split("\n")[0], "protocol": "p", "seed": 3, "n_draws": 500, "criteria": ME.CRITERIA, "proposal_rule": "r",
              "pool": {"size": 8, "n_after_exclusion": 6, "n_by_source": {"A": 8}, "approval_not_checked": [], "approval_checked_not_eligible": []}, "excluded": a["excluded"],
              "proposed": a["proposed"], "primary": a["primary"], "drugs": a["drugs"], "age_bands": a["age_bands"], "age_label_coverage_n_assessed": {"0-4": 0, "5-11": 0, "12-17": 0},
              "age_reading_disagreements_age_support_vs_engine": [], "sensitivity": a["sensitivity"],
              "comparison": {"old_rule_text": "two lines", "old_rule_candidates": {"x": 0}, "new_proposed": a["proposed"], "nine_report_drugs": {}, "proposed_not_in_report": []},
              "track1_gene_block": {"cap": 0.0, "statement": "hypotheses to test, not findings"}, "data_gaps": ["g"]}
    md = ME.to_markdown(result)
    assert "orders hypotheses" in ME.__doc__ and "Evidence matrix" in md and "Data gaps and limits" in md and "not findings" in md


# ---------------------------------------------------------------- completion (addendum section 11)

def test_completion_line_c_statistic_sign_null_and_missing():
    from track2 import multi_evidence_completion as MC

    rng = np.random.default_rng(0)
    n = 300
    score = rng.normal(size=n)
    Y = {"LETHAL_WHEN_WEAK": score + rng.normal(scale=0.3, size=n),      # low score (weak checkpoint) -> low viability: rho > 0, supportive
         "LETHAL_WHEN_STRONG": -score + rng.normal(scale=0.3, size=n),  # contrary
         "NOISE": rng.normal(size=n)}
    sparse = rng.normal(size=n)
    sparse[: n - 20] = np.nan
    Y["SPARSE"] = sparse                                                  # fewer than the minimum number of lines
    res = MC.line_c_statistics(Y, score, n_perm=200, seed=1)
    assert res["LETHAL_WHEN_WEAK"]["rho"] > 0.8 and res["LETHAL_WHEN_WEAK"]["grade"] > 0.9
    assert res["LETHAL_WHEN_STRONG"]["rho"] < -0.8 and res["LETHAL_WHEN_STRONG"]["grade"] < 0.1
    assert res["NOISE"]["grade"] == 0.5 and res["NOISE"]["q_bh"] >= MC.Q_CUT  # not distinguishable from the null: neutral
    assert res["SPARSE"]["grade"] is None and res["SPARSE"]["n_lines"] == 20
    assert 0 < res["LETHAL_WHEN_WEAK"]["p_perm"] <= 1
    assert res == MC.line_c_statistics(Y, score, n_perm=200, seed=1)  # deterministic
    assert MC.grade_from(0.5, 0.01) == 1.0 and MC.grade_from(-0.5, 0.01) == 0.0 and MC.grade_from(0.5, 0.5) == 0.5


def test_completion_bh_rank_and_approval_reuse():
    from track2 import multi_evidence_completion as MC

    q = MC.bh(np.array([0.01, 0.04, 0.03, 0.5]))
    assert np.all(q >= np.array([0.01, 0.04, 0.03, 0.5])) and np.all(q <= 1) and q[0] < q[3]
    score = np.array([0.1, 0.9, 0.5, 0.7])
    names = ["a", "b", "c", "d"]
    assert MC.rank_desc(score, names, np.array([True, True, False, True])) == {"b": 1, "d": 2, "a": 3}
    calls = []

    def check(n):
        calls.append(n)
        return {"confirmed": n != "NOPE"}

    first = MC.complete_approval(["X", "NOPE"], check)
    assert first["confirmed"] == ["X"] and first["not_confirmed"] == ["NOPE"]
    second = MC.complete_approval(["X", "NOPE", "Y"], check, first)  # answers already stored are not asked again
    assert calls == ["NOPE", "X", "Y"] and second["n_drugs"] == 3


def test_completed_pool_inputs_full_engine_and_cautions():
    rec = {"n_fda_small_molecule": 100, "drugs": {"AAA": {"rank_fda_small_molecule": 1}, "BBB": {"rank_fda_small_molecule": 51}, "CCC": {"rank_fda_small_molecule": None}}}
    grades, top, fda = ME.line_ENG_full(rec)
    assert grades["AAA"] == 1.0 and abs(grades["BBB"] - 0.5) < 1e-12 and "CCC" not in grades and top == ["AAA", "BBB"] and fda == {"AAA": True, "BBB": True}
    case = synthetic_case()
    res = ME.analyse(case, {}, n_draws=500, seed=5, mc_reps=2, mc_draws=100)
    fake = {"label": "x", "proposed": ["HHH"], "drugs": {"HHH": {"grades": {"LIT": None, "M": None}, "n_missing": 7, "confidence_factor": 0.4}}}
    text = " ".join(ME.completed_cautions(fake))
    assert "no verified literature entry" in text and "trivially 1" in text and "not a finding" in text
    fake["proposed"] = ["CALCITRIOL"]
    fake["drugs"] = {"CALCITRIOL": {"grades": {"LIT": None, "M": None}, "n_missing": 4, "confidence_factor": 0.6}}
    assert "nephrocalcinosis" in " ".join(ME.completed_cautions(fake))
    assert res["order"]
