from __future__ import annotations

from typing import Any

from track2 import plausibility


def _row(**kw: Any) -> dict[str, Any]:
    base = {"drug_id": "D1", "name": "A", "class": "constitutional_hypothesis", "is_new": True, "robustness_count": 8, "expanded_best_z": -3.0, "curated_best_z": None,
            "n_module_proteins_linked": 2, "promiscuous": False, "oral": True, "parenteral": False, "mechanism_of_action": ["x"],
            "links": {"sub_M": [{"uniprot": "P1", "gene": "MTOR", "source": "measured", "action": None, "median_pchembl": 7.2, "n_measurements": 5, "similarity": None}]}, "literature_context": [{"class": "mtor_pathway"}], "label": {"bands": {"child": "LABELLED_USE"}, "source": "reviewed"}}
    base.update(kw)
    return base


ELIGIBLE = {"eligible": True, "matched_name": "A", "marketing_statuses": ["Prescription"], "applications": ["NDA1"], "n_applications": 1}
NOT_ELIGIBLE = {**ELIGIBLE, "eligible": False}


def test_tier1_needs_l1_l2_l3_and_all_gates() -> None:
    assert plausibility.evaluate(_row(), ELIGIBLE)["tier"] == 1
    assert plausibility.evaluate(_row(literature_context=[]), ELIGIBLE)["tier"] == 2
    assert plausibility.evaluate(_row(robustness_count=3), ELIGIBLE)["tier"] == 3  # L2 only
    assert plausibility.evaluate(_row(), NOT_ELIGIBLE)["tier"] == 0  # a failed gate removes every tier


def test_gates_systemic_route_and_acute_use_class() -> None:
    assert plausibility.evaluate(_row(oral=False, parenteral=False, topical=True), ELIGIBLE)["tier"] == 0  # topical-only is out
    assert plausibility.evaluate(_row(oral=None, parenteral=None, topical=None), ELIGIBLE)["tier"] == 1  # unrecorded route is not read as non-systemic
    assert plausibility.evaluate(_row(**{"class": "excluded_acute_use_anaesthetic_sedative"}), ELIGIBLE)["tier"] == 0


def test_predicted_links_do_not_count_as_l2_and_microtubule_is_not_a_rationale() -> None:
    out = plausibility.evaluate(_row(links={"sub_M": [{"uniprot": "P1", "gene": "MTOR", "source": "predicted", "action": None, "median_pchembl": None, "similarity": 0.8}]}, literature_context=[{"class": "microtubule_targeting"}]), ELIGIBLE)
    assert out["L2_binding_curated_or_measured"] is False and out["L2_predicted_only"] is True and out["L3_mechanistic_rationale"] is False
    assert out["tier"] == 3  # L1 holds, nothing else


def test_antineoplastics_go_to_the_tumour_context_list() -> None:
    assert plausibility.evaluate(_row(**{"class": "antineoplastic_tumour_context"}), ELIGIBLE)["list"] == "tumour_context"


def test_age_band_grouping_and_hypothesis_sentence() -> None:
    out = plausibility.evaluate(_row(), ELIGIBLE, ["Tuberous Sclerosis"])
    assert "binds MTOR" in out["hypothesis"] and "tuberous sclerosis" in out["hypothesis"] and "BUB1B" not in out["hypothesis"]
    groups = plausibility.by_age_band([out])
    assert groups["child"]["labelled_for_this_band"] == ["A"] and groups["infant"]["unknown_or_verify"] == ["A"]


def test_single_measurement_link_does_not_count_as_l2() -> None:
    row = _row(links={"sub_M": [{"uniprot": "P1", "gene": "NEK9", "source": "measured", "action": None, "median_pchembl": 7.4, "n_measurements": 1, "similarity": None}]})
    out = plausibility.evaluate(row, ELIGIBLE)
    assert out["L2_binding_curated_or_measured"] is False and out["L2_single_measurement_only"] is True and out["tier"] == 3
