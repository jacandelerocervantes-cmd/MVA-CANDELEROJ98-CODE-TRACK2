"""Per-drug patient requirements with three states: PASS, FAIL and NOT_ASSESSABLE.

A requirement is evaluated only if the datum exists in `config/patient.json` (or was derived from the genome).
If it is missing, the status is NOT_ASSESSABLE and the missing datum is documented, together with the label
alerts that would make it relevant. NOT_ASSESSABLE is never turned into PASS or FAIL by assumption.
"""

from __future__ import annotations

from typing import Any

from track2.labels import band_for_age

GATES = {
    "RENAL_SAFETY": {"needs": "renal_function", "label_flag": "RENAL", "why": "documented nephrocalcinosis"},
    "HEPATIC_SAFETY": {"needs": "hepatic_function", "label_flag": "HEPATIC", "why": "hepatic elimination and hepatotoxicity"},
    "HEMATOLOGIC_SAFETY": {"needs": "blood_counts", "label_flag": "HEMATOLOGIC", "why": "cytopenias and myelosuppression"},
    "AGE_LABEL": {"needs": "age_band", "label_flag": None, "why": "label approved by age"},
    "ADMINISTRATION_ROUTE": {"needs": "can_swallow", "label_flag": None, "why": "oral or parenteral formulation"},
    "TUMOUR_CONTEXT": {"needs": "tumour_status", "label_flag": None, "why": "antitumour indication only if a tumour is present"},
    "ALLERGY_HISTORY": {"needs": "allergies", "label_flag": None, "why": "label hypersensitivity contraindications and excipients"},
}


def _nephrocalcinosis_gate(patient: dict[str, Any], label: dict[str, Any]) -> dict[str, Any]:
    """Nephrocalcinosis is documented, so a label warning about it is real evidence of concern (FAIL).

    The absence of a mention does not prove safety, so without a mention the state stays NOT_ASSESSABLE.
    """
    record = {"requirement": "NEPHROCALCINOSIS_CAUTION", "why_relevant": "nephrocalcinosis is documented in the proband",
              "needs": "documented_findings.nephrocalcinosis"}
    documented = bool(patient.get("documented_findings", {}).get("nephrocalcinosis"))
    flagged = "NEPHROCALCINOSIS_OR_STONES" in label.get("keyword_flags", [])
    if documented and flagged:
        record["status"] = "FAIL"
        record["data_used"] = "documented nephrocalcinosis + label text"
        record["label_excerpt"] = label.get("keyword_excerpts", {}).get("NEPHROCALCINOSIS_OR_STONES", "")
    else:
        record["status"] = "NOT_ASSESSABLE"
        record["data_used"] = None if not label.get("found") else "label has no mention of nephrocalcinosis or kidney stones (not proof of safety)"
    return record


AGE_PASS = {"LABELLED", "LABELLED_USE"}
AGE_FAIL = {"NO_LABELLED_USE"}


def _patient_band(patient: dict[str, Any]) -> str | None:
    """Age band from `age_band` or `age_years` (both are optional; the clinical document gives neither)."""
    m = patient["measurements"]
    if m.get("age_band"):
        return str(m["age_band"])
    return band_for_age(float(m["age_years"])) if m.get("age_years") is not None else None


def _age_gate(record: dict[str, Any], band: str | None, label: dict[str, Any]) -> None:
    profile = label.get("age_profile", {})
    status = profile.get(band) if band else None
    if band is None:
        record["status"], record["data_used"] = "NOT_ASSESSABLE", None
    elif status in AGE_PASS:
        record["status"], record["data_used"] = "PASS", f"age band {band}: label status {status} (heuristic extraction)"
    elif status in AGE_FAIL:
        record["status"], record["data_used"] = "FAIL", f"age band {band}: label status {status} (heuristic extraction)"
    else:
        record["status"], record["data_used"] = "NOT_ASSESSABLE", f"age band {band}: label status {status or 'UNKNOWN'} needs manual verification"


def evaluate_gates(patient: dict[str, Any], label: dict[str, Any], drug: dict[str, Any]) -> list[dict[str, Any]]:
    out = [_nephrocalcinosis_gate(patient, label)]
    for gate_id, rule in GATES.items():
        value = patient["measurements"].get(rule["needs"])
        record = {"requirement": gate_id, "why_relevant": rule["why"], "needs": rule["needs"]}
        flags = label.get("keyword_flags", [])
        if rule["label_flag"] and rule["label_flag"] in flags:
            record["label_alert"] = f"the label mentions {rule['label_flag'].lower()}"
        if gate_id == "ADMINISTRATION_ROUTE":
            record["drug_routes"] = {"oral": drug.get("oral"), "parenteral": drug.get("parenteral")}
        if gate_id == "AGE_LABEL":
            record["label_by_band"] = label.get("age_profile", {})
        if gate_id == "ALLERGY_HISTORY":
            record["label_contraindicates_hypersensitivity"] = bool(label.get("contraindicated_hypersensitivity"))
            record["excipient_alerts"] = label.get("excipient_alerts", [])
        if gate_id == "AGE_LABEL":
            _age_gate(record, _patient_band(patient), label)
        elif value is None:
            record["status"] = "NOT_ASSESSABLE"
            record["data_used"] = None
        else:  # no data today: the real evaluation is implemented once the datum exists
            record["status"] = "NOT_ASSESSABLE"
            record["data_used"] = "datum present but no evaluation rule implemented"
        out.append(record)
    return out
