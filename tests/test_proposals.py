from __future__ import annotations

from typing import Any

from track2 import proposals

CONFIG = {"classes": [{"id": "mtor_pathway", "match_any": ["fk506"], "claim_level": "x", "direction_hint": "extrapolation_from_mouse", "note": "n", "pmids": ["1"]},
                      {"id": "microtubule_targeting", "match_any": ["tubulin"], "claim_level": "x", "direction_hint": "unknown", "note": "n", "pmids": ["2"]}],
          "tumour_context_reports": {"applies_to_flag": "NONE", "note": "", "pmids": []}}
ELIGIBLE = {"eligible": True, "matched_name": "A", "marketing_statuses": ["Prescription"], "applications": ["NDA1"], "n_applications": 1}


def _cand(name: str, moa: str) -> dict[str, Any]:
    return {"drug_id": name, "name": name, "mechanism_of_action": [moa], "risk_flags": ["BOXED_WARNING"], "targets": ["P1"],
            "age_profile": {"child": "LABELLED_USE", "adolescent": "LABELLED_USE", "adult": "LABELLED", "neonate": "NO_LABELLED_USE", "infant": "NO_LABELLED_USE"},
            "excipients": {"alerts": [{"excipient": "LACTOSE"}]}}


def test_only_tiers_1_and_2_come_from_the_main_arm() -> None:
    plaus = {"drugs": [{"tier": t, "drug_id": f"D{t}", "name": f"D{t}", "list": "constitutional_defect", "novelty": "x", "approved_indications": [], "target_roles": [],
                        "hypothesis": "h", "L1_robust_proximity": True, "L2_binding_curated_or_measured": True, "L3_mechanistic_rationale": False, "L3_classes": [],
                        "robustness_count": 8, "expanded_best_z": -3, "age_groups": {}, "age_bands": {}, "label_min_age_years": None, "label_caution": None, "label_source": None, "boxed_warning": False,
                        "excipient_alerts": [], "fda": {}} for t in (0, 1, 2, 3)]}
    assert [r["name"] for r in proposals.main_rows(plaus)] == ["D1", "D2"]


def test_extrapolation_keeps_only_rationale_classes_that_are_fda_eligible() -> None:
    arm = [_cand("EVEROLIMUS", "FK506-binding protein 1A inhibitor"), _cand("PACLITAXEL", "Tubulin stabiliser")]
    rows = proposals.extrapolation_rows(arm, CONFIG, {}, lambda n: ELIGIBLE, lambda d: ["Tuberous Sclerosis"], direct_hits={"EVEROLIMUS": True})
    assert [r["name"] for r in rows] == ["EVEROLIMUS"]  # microtubule class is not an extrapolation rationale
    assert rows[0]["arm"] == "extrapolation" and rows[0]["tier"] is None and rows[0]["evidence"]["pmids"] == ["1"]
    none = proposals.extrapolation_rows(arm, CONFIG, {}, lambda n: {**ELIGIBLE, "eligible": False}, lambda d: [])
    assert none == []  # a drug that is not FDA-eligible is never proposed


def test_reviewed_label_sets_the_age_groups() -> None:
    review = {"EVEROLIMUS": {"verdict": "established", "min_age_years": 1.0}}
    rows = proposals.extrapolation_rows([_cand("EVEROLIMUS", "FK506-binding protein 1A inhibitor")], CONFIG, review, lambda n: ELIGIBLE, lambda d: [])
    assert rows[0]["age_groups"]["0_to_under_3"] == "partly (from 1 y)" and rows[0]["age_groups"]["3_to_under_12"] == "labelled"


def test_extrapolation_l2_is_computed_not_assumed() -> None:
    arm = [_cand("EVEROLIMUS", "FK506-binding protein 1A inhibitor")]
    yes = proposals.extrapolation_rows(arm, CONFIG, {}, lambda n: ELIGIBLE, lambda d: [], direct_hits={"EVEROLIMUS": True})
    no = proposals.extrapolation_rows(arm, CONFIG, {}, lambda n: ELIGIBLE, lambda d: [], direct_hits={})
    assert yes[0]["evidence"]["L2_binding"] is True and no[0]["evidence"]["L2_binding"] is False
