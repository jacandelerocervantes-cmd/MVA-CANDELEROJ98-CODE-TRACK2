"""Tests of labels, patient requirements and the candidate schema (no network)."""

from __future__ import annotations

import json
import re
from pathlib import Path

import jsonschema
import pytest

from track2.genome_covariates import variant_panel
from track2.gates import evaluate_gates
from track2.labels import band_profile, classify_pediatric

ROOT = Path(__file__).resolve().parents[1]


def test_pediatric_text_classification() -> None:
    assert classify_pediatric("")["status"] == "UNKNOWN"
    assert classify_pediatric("The safety and effectiveness of X have not been established in pediatric patients.")["status"] == "NOT_ESTABLISHED"
    established = classify_pediatric("The safety and effectiveness of X have been established in pediatric patients 1 year of age and older.")
    assert established["status"] == "ESTABLISHED" and established["min_age_years"] == 1.0


def test_band_profile_marks_bands_and_never_guesses_when_unknown() -> None:
    profile = band_profile({"status": "ESTABLISHED", "min_age_years": 1.0}, True)
    assert profile["neonate"] == "NO_LABELLED_USE" and profile["infant"] == "LABELLED_USE" and profile["adult"] == "LABELLED"
    assert set(band_profile({"status": "UNKNOWN", "min_age_years": None}, False).values()) == {"UNKNOWN"}
    assert band_profile({"status": "MIXED", "min_age_years": 18.0}, True)["child"] == "VERIFY_MANUALLY"


def test_gates_are_never_pass_or_fail_without_patient_data() -> None:
    patient = json.loads((ROOT / "config" / "patient.json").read_text(encoding="utf-8"))
    label = {"keyword_flags": ["RENAL"], "age_profile": {"adult": "LABELLED"}}
    gates = evaluate_gates(patient, label, {"oral": True, "parenteral": False})
    assert {g["status"] for g in gates} == {"NOT_ASSESSABLE"}  # no label mention of nephrocalcinosis, no patient data
    renal = next(g for g in gates if g["requirement"] == "RENAL_SAFETY")
    assert "label_alert" in renal and renal["needs"] == "renal_function"


def test_patient_config_has_the_eight_terms_grouped_and_no_birth_figures() -> None:
    raw = (ROOT / "config" / "patient.json").read_text(encoding="utf-8")
    patient = json.loads(raw)
    assert len(patient["hpo_terms"]) == 8 and all(t["id"].startswith("HP:") for t in patient["hpo_terms"])
    assert not re.search(r"\d+\s*weeks", raw) and not re.search(r"\d+(\.\d+)?\s*kg", raw)  # no birth figures


def _record() -> dict:
    return {
        "drug_id": "CHEMBL1", "name": "X", "approval": {"source": "ChEMBL", "chembl_version": "ChEMBL_37", "max_phase": 4},
        "targets": ["P1"], "mechanism_scores": [{"mechanism_id": "m", "z": -2.5, "rank": 1, "passes_primary_threshold": True}],
        "combined_score": 2.5, "risk_flags": [], "direction_review": "REQUIRES_REVIEW",
        "age_profile": {"adult": "LABELLED", "neonate": "UNKNOWN", "infant": "UNKNOWN", "child": "UNKNOWN", "adolescent": "UNKNOWN"},
        "gates": [{"requirement": "RENAL_SAFETY", "status": "NOT_ASSESSABLE"}],
        "toxicity": {"available": False, "source": "openFDA FAERS"}, "excipients": {"inactive_ingredients": [], "alerts": []},
        "evidence_for": [{"type": "PubMed", "id": "1", "citation_verified": True}], "evidence_against": [],
    }


def test_schema_accepts_a_valid_record_and_rejects_unapproved_or_invented_status() -> None:
    schema = json.loads((ROOT / "schemas" / "candidate_drug.schema.json").read_text(encoding="utf-8"))
    jsonschema.validate(_record(), schema)
    bad = _record()
    bad["approval"]["max_phase"] = 3                       # not approved
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)
    bad = _record()
    bad["gates"][0]["status"] = "MAYBE"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)
    bad = _record()
    bad["direction_review"] = "OK"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_variant_panel_keeps_rare_coding_variants_only() -> None:
    effects = [
        {"gene": "TPMT", "impact": "MODERATE", "chrom": "chr6", "pos": 1, "ref": "A", "alt": "G", "consequence": "missense_variant", "gt": "0/1"},
        {"gene": "TPMT", "impact": "MODERATE", "chrom": "chr6", "pos": 2, "ref": "A", "alt": "G", "consequence": "missense_variant", "gt": "0/1"},
        {"gene": "OTHER", "impact": "HIGH", "chrom": "chr1", "pos": 3, "ref": "A", "alt": "G", "consequence": "stop_gained", "gt": "0/1"},
    ]
    freqs = {"chr6:1:A:G": {"present": True, "af": 1e-5}, "chr6:2:A:G": {"present": True, "af": 0.2}}
    kept = variant_panel(effects, freqs, ["TPMT"])
    assert [v["pos"] for v in kept] == [1]


def test_nephrocalcinosis_gate_fails_only_with_documented_finding_and_label_evidence() -> None:
    patient = json.loads((ROOT / "config" / "patient.json").read_text(encoding="utf-8"))
    flagged = {"found": True, "keyword_flags": ["NEPHROCALCINOSIS_OR_STONES"], "keyword_excerpts": {"NEPHROCALCINOSIS_OR_STONES": "may cause kidney stones"}}
    gate = next(g for g in evaluate_gates(patient, flagged, {}) if g["requirement"] == "NEPHROCALCINOSIS_CAUTION")
    assert gate["status"] == "FAIL" and "kidney stones" in gate["label_excerpt"]
    silent = next(g for g in evaluate_gates(patient, {"found": True, "keyword_flags": []}, {}) if g["requirement"] == "NEPHROCALCINOSIS_CAUTION")
    assert silent["status"] == "NOT_ASSESSABLE"  # no mention is not proof of safety
    patient["documented_findings"]["nephrocalcinosis"] = False
    assert next(g for g in evaluate_gates(patient, flagged, {}) if g["requirement"] == "NEPHROCALCINOSIS_CAUTION")["status"] == "NOT_ASSESSABLE"


def test_allergy_gate_is_not_assessable_without_history_and_carries_excipient_alerts() -> None:
    patient = json.loads((ROOT / "config" / "patient.json").read_text(encoding="utf-8"))
    label = {"keyword_flags": [], "contraindicated_hypersensitivity": True, "excipient_alerts": [{"excipient": "BENZYL ALCOHOL", "concern": "x"}]}
    gate = next(g for g in evaluate_gates(patient, label, {}) if g["requirement"] == "ALLERGY_HISTORY")
    assert gate["status"] == "NOT_ASSESSABLE" and gate["label_contraindicates_hypersensitivity"] and gate["excipient_alerts"]


def test_excipient_alerts_and_toxicity_groups() -> None:
    from track2.labels import excipient_alerts
    from track2.toxicity import classify_terms
    alerts = excipient_alerts(["BETADEX SULFOBUTYL ETHER SODIUM", "MANNITOL", "BENZYL ALCOHOL"])
    assert {a["excipient"] for a in alerts} == {"BETADEX SULFOBUTYL ETHER SODIUM", "BENZYL ALCOHOL"}
    assert [x["excipient"] for x in excipient_alerts(["ALCOHOL", "POLYVINYL ALCOHOL", "ISOPROPYL ALCOHOL"])] == ["ALCOHOL"]
    shares = classify_terms([{"term": "ACUTE KIDNEY INJURY", "count": 30}, {"term": "NEUTROPENIA", "count": 50}, {"term": "NAUSEA", "count": 20}])
    assert shares["renal"] == 0.3 and shares["haematological"] == 0.5 and shares["hepatic"] == 0.0
