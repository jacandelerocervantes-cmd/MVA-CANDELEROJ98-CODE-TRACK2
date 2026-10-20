"""Offline tests for track2.case_profile (synthetic inputs only)."""

import copy
import json

import pytest

from track2 import case_profile as cp


def minimal():
    return cp.build_profile("synthetic_1", when="2026-01-01")


def test_all_unknown_by_default():
    p = minimal()
    assert p["age_band"] == "unknown" and p["sex"] == "unknown"
    assert p["tumour"] == {"history": "unknown", "current_status": "unknown"}
    assert set(p["findings"].values()) == {"unknown"}
    assert cp.validate_profile(p) == []


def test_full_profile_validates():
    p = cp.build_profile(
        "synthetic_2", age_band="3-<12", sex="female",
        phenotypic_features=[{"hpo_id": "HP:0000252", "status": "present"}, {"hpo_id": "HP:0001250", "status": "absent"}],
        tumour_history="yes", tumour_current_status="remission", renal="present", growth="present", neurological="absent",
        genomic_findings=[{"gene_symbol": "GENE1", "mechanism_class": "spindle-assembly checkpoint",
                           "consequence_classes": ["missense", "loss of function"], "genotype_pattern": "compound heterozygous",
                           "phase_status": "unresolved", "computational_evidence_strength": "moderate"}],
        source_files=["config/patient.json"], when="2026-01-01")
    assert cp.validate_profile(p) == []
    assert p["genomic_findings"][0]["consequence_classes"] == ["loss of function", "missense"]


def test_schema_is_draft_2020_12():
    assert cp.SCHEMA["$schema"].endswith("2020-12/schema")
    cp.Draft202012Validator.check_schema(cp.SCHEMA)


@pytest.mark.parametrize("field,value", [("age_band", "5"), ("sex", "M?"), ("age_band", "adult")])
def test_closed_enums_reject_guesses(field, value):
    p = minimal()
    p[field] = value
    assert cp.validate_profile(p)


def test_unknown_fields_rejected():
    p = minimal()
    p["position"] = 5
    assert cp.validate_profile(p)
    p = minimal()
    p["tumour"]["extra"] = "x"
    assert cp.validate_profile(p)


def test_bad_hpo_status_and_id():
    for bad in ({"hpo_id": "HP:12", "status": "present"}, {"hpo_id": "HP:0000252", "status": "maybe"}):
        p = minimal()
        p["phenotypic_features"] = [bad]
        assert cp.validate_profile(p)


def test_conflicting_status_rejected():
    with pytest.raises(ValueError):
        cp.build_profile("x", phenotypic_features=[{"hpo_id": "HP:0000252", "status": "present"},
                                                   {"hpo_id": "HP:0000252", "status": "absent"}])


GENE_FIELD = ["gene_symbol", "mechanism_class", "consequence_classes", "genotype_pattern", "phase_status",
              "computational_evidence_strength"]
PRIVATE_STRINGS = ["chr7:55181378", "7:55181378", "c.123A>T", "p.Arg893Ser", "A>G", "DP=34", "AF: 0.48", "VAF 0.3",
                   "rs1422532977", "0/1", "read depth 30", "g.12345678"]


@pytest.mark.parametrize("text", PRIVATE_STRINGS)
def test_privacy_no_free_text_field_accepts_coordinates_alleles_or_read_metrics(text):
    # every string-valued slot of the profile must reject a private-looking value
    slots = []
    p = minimal()
    p["genomic_findings"] = [{"gene_symbol": "GENE1", "mechanism_class": "unknown", "consequence_classes": ["unknown"],
                              "genotype_pattern": "unknown", "phase_status": "unknown",
                              "computational_evidence_strength": "unknown"}]
    p["provenance"]["notes"] = ["note"]
    p["phenotypic_features"] = [{"hpo_id": "HP:0000252", "status": "present"}]
    assert cp.validate_profile(p) == []

    def paths(node, path=()):
        if isinstance(node, str):
            slots.append(path)
        elif isinstance(node, dict):
            for k, v in node.items():
                paths(v, path + (k,))
        elif isinstance(node, list):
            for i, v in enumerate(node):
                paths(v, path + (i,))
    paths(p)
    assert slots
    for path in slots:
        q = copy.deepcopy(p)
        ref = q
        for k in path[:-1]:
            ref = ref[k]
        ref[path[-1]] = text
        assert cp.validate_profile(q), f"slot {path} accepted {text!r}"


def test_privacy_scan_catches_key_names_and_nested_values():
    p = minimal()
    p["provenance"]["notes"] = ["depth DP=30"]
    assert any("read metric" in e for e in cp.privacy_violations(p))
    assert cp.privacy_violations({"chr1:100000": 1})


def test_hpo_identifiers_are_not_mistaken_for_positions():
    p = cp.build_profile("x", phenotypic_features=[{"hpo_id": "HP:0000252", "status": "present"}])
    assert cp.privacy_violations(p) == []


def test_child_profile_from_synthetic_files(tmp_path):
    (tmp_path / "config").mkdir()
    (tmp_path / "out").mkdir()
    (tmp_path / "config/patient.json").write_text(json.dumps({
        "hpo_terms": [{"id": "HP:0002859", "label": "x", "group": "cancer_predisposition"},
                      {"id": "HP:0004322", "label": "y", "group": "growth_restriction"}],
        "documented_findings": {"tumour_history": True, "nephrocalcinosis": None},
        "measurements": {"age_band": None, "tumour_status": None, "neurological_findings": None},
        "genome_derived": {"sex": None}}), encoding="utf-8")
    (tmp_path / "out/differential_diagnosis.json").write_text(json.dumps({"genes": {
        "GENE1": {"genotype_type": "compound_het", "c3": "high_impact", "reason": "phase not resolved at gene level"},
        "GENE2": {"genotype_type": None, "c3": "structural_unpaired", "reason": ""}}}), encoding="utf-8")
    (tmp_path / "out/variant_consequence_classes.json").write_text(json.dumps({"genes": {
        "GENE1": {"class_counts": {"missense": 1, "loss of function": 1}}}}), encoding="utf-8")
    (tmp_path / "out/mva_landscape.json").write_text(json.dumps({"genes": [
        {"symbol": "GENE1", "mechanism_class": "spindle-assembly checkpoint"}]}), encoding="utf-8")
    p = cp.build_child_profile(tmp_path, when="2026-01-01")
    assert cp.validate_profile(p) == []
    assert p["age_band"] == "unknown" and p["sex"] == "unknown"
    assert p["tumour"]["history"] == "yes" and p["tumour"]["current_status"] == "unknown"
    assert p["findings"] == {"renal": "unknown", "growth": "present", "neurological": "unknown"}
    g = {x["gene_symbol"]: x for x in p["genomic_findings"]}
    assert g["GENE1"]["genotype_pattern"] == "compound heterozygous" and g["GENE1"]["phase_status"] == "unresolved"
    assert g["GENE2"]["genotype_pattern"] == "structural unpaired" and g["GENE2"]["consequence_classes"] == ["unknown"]
    assert g["GENE1"]["computational_evidence_strength"] == "unknown"
    assert "out/empirical/variant_calibration.json" not in p["provenance"]["source_files"]
    feats = {f["hpo_id"]: f["status"] for f in p["phenotypic_features"]}
    assert feats["HP:0002859"] == "present" and feats["HP:0000252"] == "unknown"


def test_calibration_and_mechanism_summarised_categorically(tmp_path):
    (tmp_path / "config").mkdir()
    (tmp_path / "out/empirical").mkdir(parents=True)
    (tmp_path / "config/patient.json").write_text(json.dumps({"hpo_terms": []}), encoding="utf-8")
    (tmp_path / "out/differential_diagnosis.json").write_text(json.dumps({"genes": {
        "BUB1B": {"genotype_type": "compound_het", "reason": ""}, "GENE2": {"genotype_type": "hemizygous"}}}), encoding="utf-8")
    (tmp_path / "out/variant_consequence_classes.json").write_text(json.dumps({"genes": {}}), encoding="utf-8")
    (tmp_path / "out/mva_landscape.json").write_text(json.dumps({"genes": []}), encoding="utf-8")
    (tmp_path / "out/empirical/variant_calibration.json").write_text(json.dumps({"inputs": {"scores": {"CADD": 24.5}}, "per_predictor": {
        "AlphaMissense": {"strength": "PP3_Moderate", "score": 0.9}, "PolyPhen2": {"strength": "PP3_Supporting"},
        "REVEL": {"strength": "none"}, "CADD": {"strength": "none"}, "SIFT": {"strength": "none"}}}), encoding="utf-8")
    (tmp_path / "out/variant_mechanism.json").write_text(json.dumps({
        "method_constants": {"conserved_fraction_min": 0.75, "buried_sasa_percentile_max": 25.0},
        "variant_a": {"nmd": {"nmd_predicted": True}, "truncated_product_features": {"lost_entirely": [{"description": "Protein kinase"}]}},
        "variant_b": {"conservation": {"all_one2one_orthologs": {"fraction_identical_at_position_of_aligned": 0.9}},
                      "structure": {"sasa_percentile_among_same_aa": 8.0}, "domain_source": {"description": "Protein kinase"}}}), encoding="utf-8")
    p = cp.build_child_profile(tmp_path, when="2026-01-01")
    assert cp.validate_profile(p) == []
    g = {x["gene_symbol"]: x for x in p["genomic_findings"]}
    assert g["BUB1B"]["computational_evidence_strength"] == "discordant across predictors"
    assert {x["predictor"]: x["strength"] for x in g["BUB1B"]["predictor_strengths"]} == {
        "AlphaMissense": "moderate", "PolyPhen-2": "supporting", "REVEL": "none", "CADD": "none", "SIFT": "none"}
    assert g["BUB1B"]["variant_mechanism_statements"] == cp.MECH_STATEMENTS
    assert "predictor_strengths" not in g["GENE2"] and g["GENE2"]["computational_evidence_strength"] == "unknown"
    assert "0.9" not in json.dumps(p) and "24.5" not in json.dumps(p)


def test_mechanism_statements_need_rule_support_and_reject_numbers():
    assert cp.mechanism_statements({"variant_x": {"conservation": {}, "structure": {}}}) == []
    p = cp.build_profile("x", genomic_findings=[{"gene_symbol": "GENE1"}])
    p["genomic_findings"][0]["variant_mechanism_statements"] = ["missense at residue 1002"]
    assert cp.validate_profile(p)
    p["genomic_findings"][0]["variant_mechanism_statements"] = ["chr1:12345"]
    assert cp.validate_profile(p)
