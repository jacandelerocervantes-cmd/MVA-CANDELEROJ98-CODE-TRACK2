from __future__ import annotations

from track2 import age_view
from track2.gates import evaluate_gates
from track2.labels import band_for_age


def test_band_for_age_boundaries() -> None:
    assert [band_for_age(a) for a in (0.0, 0.5, 1.9, 2.0, 11.99, 12.0, 17.9, 18.0, 40.0)] == [
        "neonate", "infant", "infant", "child", "child", "adolescent", "adolescent", "adult", "adult"]


def _age_status(patient_measurements: dict, profile: dict) -> str:
    patient = {"documented_findings": {}, "measurements": patient_measurements}
    gates = evaluate_gates(patient, {"age_profile": profile, "keyword_flags": []}, {})
    return next(g["status"] for g in gates if g["requirement"] == "AGE_LABEL")


def test_age_gate_not_assessable_without_age_and_evaluated_with_it() -> None:
    profile = {"child": "LABELLED_USE", "infant": "NO_LABELLED_USE", "adolescent": "VERIFY_MANUALLY", "adult": "LABELLED"}
    assert _age_status({"age_band": None, "age_years": None}, profile) == "NOT_ASSESSABLE"
    assert _age_status({"age_band": None, "age_years": 7}, profile) == "PASS"
    assert _age_status({"age_band": "infant"}, profile) == "FAIL"
    assert _age_status({"age_band": "adolescent"}, profile) == "NOT_ASSESSABLE"  # needs manual verification, never assumed


def test_view_groups_by_band() -> None:
    cand = [{"drug_id": "D1", "name": "A", "risk_flags": [], "excipients": {"alerts": [{"excipient": "BENZYL ALCOHOL"}]},
             "age_profile": {"neonate": "NO_LABELLED_USE", "infant": "UNKNOWN", "child": "LABELLED_USE", "adolescent": "VERIFY_MANUALLY", "adult": "LABELLED"},
             "layers": {"robustness": {"class": "robust"}}}]
    groups = age_view.by_band(age_view.view(cand, "main"))
    assert groups["child"]["labelled_use"] == ["A"] and groups["neonate"]["no_labelled_use"] == ["A"]
    assert groups["infant"]["unknown_or_verify"] == ["A"] and groups["adult"]["labelled_use"] == ["A"]
