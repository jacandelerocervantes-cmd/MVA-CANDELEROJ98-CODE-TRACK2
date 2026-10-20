from __future__ import annotations

from track2 import repositioning

CFG = repositioning.load()


def test_proteasome_inhibitor_is_tumour_directed() -> None:
    out = repositioning.proposed_use(["26S proteasome inhibitor"], [], "curated_mechanism_already_near", CFG)
    assert out["rule"] == "proteasome_inhibitor" and out["context"] == "tumour" and out["consequences"] == ["2", "5"]
    assert "does not treat the constitutional defect" in out["new_use"]


def test_adc_names_the_antigen_and_is_not_confused_with_a_plain_microtubule_agent() -> None:
    out = repositioning.proposed_use(["Tubulin inhibitor", "Nectin-4 binding agent"], [], "curated_mechanism_already_near", CFG, "tumour_context")
    assert out["rule"] == "tubulin_toxin_antibody_drug_conjugate" and "Nectin-4" in out["new_use"]
    plain = repositioning.proposed_use(["Tubulin stabiliser"], [], "curated_mechanism_already_near", CFG, "tumour_context")
    assert plain["rule"] == "microtubule_agent_tumour" and "not favoured" in plain["caveat"]


def test_mtor_class_ties_to_growth_and_senescence() -> None:
    out = repositioning.proposed_use(["FK506-binding protein 1A inhibitor"], [], "x", CFG)
    assert out["rule"] == "mtor_pathway" and out["consequences"] == ["3", "4"]


def test_new_drug_with_measured_links_gets_a_hypothesis_with_its_genes_and_unknown_direction() -> None:
    out = repositioning.proposed_use(["Spleen tyrosine kinase inhibitor"], ["NEK9", "PLK4", "AURKA", "X"], "new_via_expanded_links", CFG)
    assert out["rule"] == "measured_links" and "NEK9, PLK4, AURKA" in out["new_use"] and "direction of effect is unknown" in out["new_use"]


def test_no_rule_and_not_new_means_no_use_is_invented() -> None:
    out = repositioning.proposed_use(["Something else"], ["G1"], "curated_mechanism_already_near", CFG)
    assert out["rule"] == "none" and out["consequences"] == []


def test_microtubule_agent_use_depends_on_whether_it_is_an_oncology_drug() -> None:
    onc = repositioning.proposed_use(["Tubulin inhibitor"], [], "x", CFG, "tumour_context")
    non = repositioning.proposed_use(["Tubulin inhibitor"], [], "x", CFG, "constitutional_defect")
    assert onc["rule"] == "microtubule_agent_tumour" and non["rule"] == "microtubule_agent_non_oncology"
    assert non["consequences"] == [] and "No specific repositioned use" in non["new_use"]
