"""Offline tests of the empirical age-support function with synthetic labels."""

import pytest

from track2 import empirical_age_support as e


def test_label_support_shapes():
    assert e.label_support("established", 13.0) == {"state": "known", "intervals": [[13.0, 18.0]]}
    assert e.label_support("established", None)["intervals"] == [[0.0, 18.0]]
    assert e.label_support("established", 20.0)["intervals"] == []
    assert e.label_support("not_established", None)["intervals"] == []
    assert e.label_support("verify", 5.0)["state"] == "unknown"


def test_measure_and_bounds():
    L = e.label_support("established", 2.0)
    assert e.support_bounds(L, [])["lower"] == pytest.approx(16 / 18, abs=1e-5)
    assert e.support_bounds(L, [[0.0, 6.0]])["upper"] == pytest.approx(12 / 18, abs=1e-5)
    U = e.label_support("verify", None)
    b = e.support_bounds(U, [[0.0, 6.0]])
    assert b["lower"] == 0.0 and b["upper"] == pytest.approx(12 / 18, abs=1e-5)
    assert e.support_bounds(L, [], 0.0, 3.0)["lower"] == pytest.approx(1 / 3, abs=1e-5)  # 2-3 y inside 0-<3
    assert e.support_bounds(L, [[0, 18]])["upper"] == 0.0


def test_extraction_age_and_weight_verbatim():
    sections = {
        "warnings": "Fatal toxicity occurred in children under 6 years of age. Take with food. Not recommended in patients less than 13 years. "
                    "Safety and effectiveness have been established in patients 2 years and older.",
        "dosage": "Tablets are not recommended for neonates, infants and children up to 20 kg. Adverse reactions in patients aged <2 years were similar.",
    }
    c = e.extract_cautions(sections)
    sents = [x["sentence"] for x in c]
    assert len(c) == 3
    for x in c:
        assert e._norm(x["sentence"]) in e._norm(sections[x["section"]])
    ages = {x["sentence"][:12]: x["age_intervals_years"] for x in c}
    assert [0.0, 6.0] in ages["Fatal toxici"] and [0.0, 13.0] in ages["Not recommen"]
    w = [x for x in c if x["weight_cutoffs"]][0]
    assert w["weight_cutoffs"][0]["cutoff_kg"] == 20.0 and w["weight_flag_not_converted"]
    assert any(0.0 < iv[1] < 0.1 for iv in w["age_intervals_years"])  # neonates -> 28 days
    assert not any("established in patients 2 years" in s for s in sents)
    assert e.h_intervals(c) == [[0.0, 13.0]]


def test_months_and_treatment_duration_ignored():
    c = e.extract_cautions({"warnings": "Not recommended in infants under 6 months of age. Not recommended for use longer than 12 weeks. Avoid in patients less than 4 weeks of therapy."})
    assert len(c) == 1 and c[0]["age_intervals_years"] == [[0.0, 0.5]]


def test_drug_record_and_groups():
    sections = {"pediatric_use": "Use is not recommended in children under 6 years."}
    r = e.drug_record("X", "CHEMBL0", {"verdict": "established", "min_age_years": 2.0, "scope": "s"}, sections, "sid")
    assert r["H_intervals_years"] == [[0.0, 6.0]] and r["n_hazards"] == 1
    assert r["C"]["lower"] == pytest.approx(12 / 18, abs=1e-5)
    assert r["C_by_group"]["0_to_under_3"]["lower"] == 0.0
    assert r["C_by_group"]["3_to_under_12"]["lower"] == pytest.approx(6 / 9, abs=1e-5)
    assert r["C_by_group"]["12_and_older"]["lower"] == 1.0


def test_build_uses_loader_and_sorted_header():
    table = {"farmacos": [{"name": "B", "drug_id": "C2", "etiqueta": {"setid": "s2"}}, {"name": "A", "drug_id": "C1", "etiqueta": {"setid": "s1"}}]}
    config = {"entries": {"C1": {"name": "A", "verdict": "verify"}, "C2": {"name": "B", "verdict": "not_established"}}}
    out = e.build(table, config, loader=lambda s: {"warnings": "Contraindicated in neonates."})
    assert [d["name"] for d in out["drugs"]] == ["A", "B"]
    assert out["drugs"][0]["C"]["lower"] == 0.0 and out["drugs"][0]["C"]["upper"] > 0.99
    assert out["drugs"][1]["C"]["upper"] == 0.0
    assert "NOT a recommendation" in out["header"]["not_a_recommendation"] and "No dose" in out["header"]["not_a_recommendation"]


def test_not_established_does_not_reduce_support_and_weight_flag():
    sections = {"pediatric_use": "Safety and efficacy have not been established in pediatric patients <18 years. The use of tablets for neonates, infants and children to the weight of 20 kg is not recommended.",
                "warnings": "A study in patients 6 to 144 months post-transplantation showed more deaths."}
    r = e.drug_record("X", "C", {"verdict": "established", "min_age_years": 13.0}, sections, "s")
    assert r["C"]["lower"] == pytest.approx(5 / 18, abs=1e-4)  # only the neonate hazard (28 days) is subtracted outside 13-18
    assert r["sensitivity_C_if_not_established_sentences_counted_as_hazards"]["lower"] == 0.0
    assert r["weight_flags"][0]["weight_cutoffs"][0]["cutoff_kg"] == 20.0
    reasons = " ".join(d["reason"] for d in r["dropped"])
    assert "restates" in reasons and "time since" in reasons
