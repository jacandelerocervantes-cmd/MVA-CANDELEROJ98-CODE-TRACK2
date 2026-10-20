"""Tests of the addendum-11 modules (E1 to E15) on small synthetic inputs."""

from __future__ import annotations

import numpy as np

from track2.class_analysis import bh, classes_of, evaluate_classes
from track2.fill_methods_form import answers, plain
from track2.four_lines_review2 import distinct_top, overlap
from track2.phenotype_module import NOT_THE_CHILD, query_terms
from track2.vulnerability_core import hits, scan

PATIENT = {"hpo_terms": [{"id": "HP:0000121"}, {"id": "HP:0004322"}, {"id": "HP:0200067"}]}


def test_parental_term_is_not_a_phenotype_query_term() -> None:
    assert "HP:0200067" in NOT_THE_CHILD
    assert query_terms(PATIENT) == ["HP:0000121", "HP:0004322"]
    assert query_terms(PATIENT, leave_out={"HP:0004322"}) == ["HP:0000121"]
    assert query_terms(PATIENT, only={"HP:0004322"}) == ["HP:0004322"]


def test_distinct_top_merges_nodes_that_map_to_one_drug() -> None:
    ranks = {"A": 1, "B": 2, "C": 3, "D": 4}
    assert distinct_top(ranks, 2) == ["A", "B"]
    # a drug that appears with two ranks keeps the best one (the dict keys are drug ids, so duplicates are merged upstream: best rank wins)
    assert distinct_top({"A": 5, "B": 1}, 1) == ["B"]


def test_overlap_uses_the_common_universe_only() -> None:
    res = overlap({"A": {"x", "y", "z"}, "B": {"y", "w"}}, common={"x", "y", "z"})
    assert res["A_and_B"]["overlap"] == 1 and res["A_and_B"]["size_second"] == 1  # "w" is outside the common universe


def test_bh_is_monotone_and_bounded() -> None:
    q = bh(np.array([0.001, 0.02, 0.5, 0.04]))
    assert (q <= 1).all() and q[0] <= q[1] <= q[3] <= q[2] + 1e-12


def test_class_analysis_finds_an_enriched_class_and_not_a_random_one() -> None:
    rng = np.random.default_rng(0)
    drugs = [f"d{i}" for i in range(200)]
    goodness = {d: float(rng.normal()) for d in drugs}
    degree = {d: float(rng.integers(1, 50)) for d in drugs}
    top = sorted(drugs, key=lambda d: -goodness[d])[:6]
    rand = list(rng.choice(drugs, size=6, replace=False))
    res = evaluate_classes(goodness, degree, {"enriched": set(top), "random": {str(x) for x in rand}}, permutations=500, seed=1)
    by = {r["class"]: r for r in res["classes"]}
    assert by["enriched"]["p"] < 0.01
    assert by["random"]["p"] > by["enriched"]["p"]


def test_class_analysis_skips_tiny_classes_and_unknown_mechanisms() -> None:
    universe = {"drugs": [{"drug_id": "a", "atc": ["L01EX01"], "mechanisms": [{"mechanism": "Unknown"}]},
                          {"drug_id": "b", "atc": ["L01EX02"], "mechanisms": [{"mechanism": "Kinase inhibitor"}]}]}
    cls = classes_of(universe)
    assert "Unknown" not in cls["moa"] and cls["atc4"]["L01EX"] == {"a", "b"}
    res = evaluate_classes({"a": 1.0, "b": 2.0}, {"a": 1.0, "b": 1.0}, cls["atc4"], permutations=50)
    assert res["n_classes_tested"] == 0  # fewer than three members


def test_tp53_covariate_removes_a_confounded_signal() -> None:
    """Sensitivity depends on TP53 and the checkpoint score is correlated with TP53: without the covariate the score looks like a hit, with it the effect disappears."""
    rng = np.random.default_rng(3)
    n = 400
    tp53 = (rng.random(n) < 0.5).astype(float)
    score = tp53 * 1.0 + rng.normal(size=n) * 0.6
    prolif = rng.normal(size=n)
    lineage = np.array(["L"] * n)
    y = 0.8 * (1 - tp53) + rng.normal(size=n) * 0.1  # low where TP53 is mutated; unrelated to the score except through TP53
    Y = np.vstack([y + rng.normal(size=n) * 0.01 for _ in range(3)])
    plain_scan = scan(Y, score, prolif, lineage)
    adjusted = scan(Y, score, prolif, lineage, extra=tp53[:, None])
    assert np.nanmean(np.abs(adjusted["coef"])) < 0.25 * np.nanmean(np.abs(plain_scan["coef"]))
    assert not hits(adjusted).any()


def test_form_answers_are_split_by_question_number_and_stripped_of_markup() -> None:
    md = "### 1. Team name\nIndividual **Submission**\n\n### 2. Approach\nText with `code`.\n\n### 11. Abstract\nBody\n\n(Abstract length: 3 words.)\n"
    a = answers(md)
    assert a[1] == "Individual Submission" and a[2] == "Text with code." and a[11] == "Body"
    assert plain("**b** and *i*") == "b and i"


def test_quote_checker_ignores_whitespace_and_decodes_entities_only() -> None:
    from track2.verify_quotes import squash

    assert squash("mm 3 &lt; 1500") == squash("mm3<1500")
    assert squash("neutrophil  count") != squash("neutrophil counts")


def test_label_levels_cover_the_same_drugs_and_respect_the_config_verdicts() -> None:
    import json
    from pathlib import Path

    from track2.label_levels import ALLOWED, LEVEL

    entries = json.loads(Path("config/pediatric_label_review.json").read_text(encoding="utf-8"))["entries"].values()
    reviewed = {e["name"]: e["verdict"] for e in entries}
    assert len(LEVEL) == 71 and set(LEVEL.values()) <= {1, 2, 3, 4, 5, 6} and set(LEVEL) <= set(reviewed)
    assert LEVEL["DACARBAZINE"] == 5 and LEVEL["PAZOPANIB"] == 6 and LEVEL["TREOSULFAN"] == 2 and LEVEL["CARBOPLATIN"] == 4 and LEVEL["NATALIZUMAB"] == 6
    assert ALLOWED["verify"] == {4} and 3 in ALLOWED["established"] and 3 not in ALLOWED["not_established"]
    for name, level in LEVEL.items():
        if reviewed[name] in ALLOWED:
            assert level in ALLOWED[reviewed[name]], (name, level, reviewed[name])
