"""Minimal, privacy-safe case profile (inspired by the GA4GH Phenopacket v2 core).

One small record describes any individual (this child, another person with MVA, a person with another
rare disease). Every field has an explicit "unknown" state; a value that is not documented is never
guessed. Genomic findings are categorical only: gene symbol, mechanism class, consequence class,
genotype pattern, phase status and predicted computational evidence strength. No field can hold a
coordinate, an allele, a genotype call or a read metric: the schema uses closed patterns and the
validator additionally rejects any string that looks like one.
"""

from __future__ import annotations

import json
import re
from datetime import date
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator

SCHEMA_VERSION = "case_profile/1.0"
UNKNOWN = "unknown"

AGE_BANDS = ["0-<3", "3-<12", "12 and older", UNKNOWN]
SEX_VALUES = ["female", "male", "other", UNKNOWN]
STATUS = ["present", "absent", UNKNOWN]
TUMOUR_HISTORY = ["yes", "no", UNKNOWN]
TUMOUR_CURRENT = ["active", "remission", "none", UNKNOWN]
GENOTYPE_PATTERNS = [
    "compound heterozygous", "homozygous", "hemizygous", "structural unpaired",
    "single heterozygous", UNKNOWN,
]
PHASE_STATUS = ["resolved", "unresolved", UNKNOWN]
EVIDENCE_STRENGTH = ["very strong", "strong", "moderate", "supporting", "indeterminate", "discordant across predictors", UNKNOWN]
PREDICTORS = ["AlphaMissense", "PolyPhen-2", "REVEL", "CADD", "SIFT"]
PREDICTOR_STRENGTH = ["moderate", "supporting", "none", UNKNOWN]
_CAL_KEYS = {"AlphaMissense": "AlphaMissense", "PolyPhen2": "PolyPhen-2", "REVEL": "REVEL", "CADD": "CADD", "SIFT": "SIFT"}
_STRENGTH_WORDS = {"PP3_Moderate": "moderate", "PP3_Supporting": "supporting", "none": "none"}
CALIBRATION_FILE = "out/empirical/variant_calibration.json"
MECHANISM_FILE = "out/variant_mechanism.json"
MECH_LOF = "loss of function by predicted nonsense-mediated decay; any truncated product lacks the kinase domain"
MECH_MISSENSE = "missense at a conserved, buried position of the kinase domain"
MECH_STATEMENTS = [MECH_LOF, MECH_MISSENSE]

_CLASS_PATTERN = r"^[a-z][a-z0-9 ,;/()-]{1,79}$"  # lower-case words only: no digits-only tokens, no ':' or '>'
_GENE_PATTERN = r"^[A-Z][A-Z0-9-]{1,14}$"
_HPO_PATTERN = r"^HP:[0-9]{7}$"
_STATEMENT_PATTERN = r"^[a-z][a-z ,;/()-]{1,159}$"  # words only: no digits at all
_FILE_PATTERN = r"^[A-Za-z0-9_./-]{1,120}$"
_DATE_PATTERN = r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$"
_ID_PATTERN = r"^[A-Za-z0-9_-]{1,40}$"

# Strings that look like a coordinate, an allele, a genotype call or a read metric are rejected
# whatever field they are in (defence in depth on top of the closed patterns above).
_FORBIDDEN = [
    (re.compile(r"chr[0-9XYM]+", re.I), "chromosome name"),
    (re.compile(r"\b[0-9XYM]{1,2}\s*[:_-]\s*[0-9]{3,}"), "chromosome:position"),
    (re.compile(r"\b[cgmnpr]\.\s*[0-9*_-]"), "HGVS description"),
    (re.compile(r"\bp\.\(?[A-Z][a-z]{2}[0-9]"), "HGVS protein description"),
    (re.compile(r"[ACGTN]+\s*>\s*[ACGTN]+", re.I), "allele change"),
    (re.compile(r"\b[01][/|][01]\b"), "genotype call"),
    (re.compile(r"\b(DP|AD|AF|GQ|VAF|QUAL|MQ)\s*[=:]", re.I), "read metric"),
    (re.compile(r"\b(read depth|allele fraction|allele frequency|vaf|coverage)\b", re.I), "read metric"),
    (re.compile(r"\brs[0-9]{3,}\b", re.I), "variant identifier"),
    (re.compile(r"(?<![0-9])[0-9]{5,}(?![0-9])"), "long number (possible position)"),
]


def _status_enum() -> dict:
    return {"type": "string", "enum": STATUS}


SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://example.org/track2/case_profile.schema.json",
    "title": "Minimal case profile (Phenopacket v2 core, categorical and privacy-safe)",
    "type": "object",
    "additionalProperties": False,
    "required": ["schema_version", "profile_id", "age_band", "sex", "phenotypic_features", "tumour",
                 "findings", "genomic_findings", "provenance"],
    "properties": {
        "schema_version": {"const": SCHEMA_VERSION},
        "profile_id": {"type": "string", "pattern": _ID_PATTERN},
        "age_band": {"type": "string", "enum": AGE_BANDS},
        "sex": {"type": "string", "enum": SEX_VALUES},
        "phenotypic_features": {
            "type": "array",
            "items": {
                "type": "object", "additionalProperties": False, "required": ["hpo_id", "status"],
                "properties": {"hpo_id": {"type": "string", "pattern": _HPO_PATTERN}, "status": _status_enum()},
            },
        },
        "tumour": {
            "type": "object", "additionalProperties": False, "required": ["history", "current_status"],
            "properties": {
                "history": {"type": "string", "enum": TUMOUR_HISTORY},
                "current_status": {"type": "string", "enum": TUMOUR_CURRENT},
            },
        },
        "findings": {
            "type": "object", "additionalProperties": False, "required": ["renal", "growth", "neurological"],
            "properties": {"renal": _status_enum(), "growth": _status_enum(), "neurological": _status_enum()},
        },
        "genomic_findings": {
            "type": "array",
            "items": {
                "type": "object", "additionalProperties": False,
                "required": ["gene_symbol", "mechanism_class", "consequence_classes", "genotype_pattern",
                             "phase_status", "computational_evidence_strength"],
                "properties": {
                    "gene_symbol": {"type": "string", "pattern": _GENE_PATTERN},
                    "mechanism_class": {"type": "string", "pattern": _CLASS_PATTERN},
                    "consequence_classes": {"type": "array", "items": {"type": "string", "pattern": _CLASS_PATTERN}},
                    "genotype_pattern": {"type": "string", "enum": GENOTYPE_PATTERNS},
                    "phase_status": {"type": "string", "enum": PHASE_STATUS},
                    "computational_evidence_strength": {"type": "string", "enum": EVIDENCE_STRENGTH},
                    "predictor_strengths": {
                        "type": "array",
                        "items": {
                            "type": "object", "additionalProperties": False, "required": ["predictor", "strength"],
                            "properties": {"predictor": {"type": "string", "enum": PREDICTORS},
                                           "strength": {"type": "string", "enum": PREDICTOR_STRENGTH}},
                        },
                    },
                    "variant_mechanism_statements": {"type": "array", "items": {"type": "string", "pattern": _STATEMENT_PATTERN}},
                },
            },
        },
        "provenance": {
            "type": "object", "additionalProperties": False, "required": ["source_files", "date"],
            "properties": {
                "source_files": {"type": "array", "items": {"type": "string", "pattern": _FILE_PATTERN}},
                "date": {"type": "string", "pattern": _DATE_PATTERN},
                "notes": {"type": "array", "items": {"type": "string", "pattern": _CLASS_PATTERN}},
            },
        },
    },
}

_VALIDATOR = Draft202012Validator(SCHEMA)


def _strings(node: Any, path: str = ""):
    if isinstance(node, str):
        yield path, node
    elif isinstance(node, dict):
        for k, v in node.items():
            yield from _strings(k, f"{path}/<key>")
            yield from _strings(v, f"{path}/{k}")
    elif isinstance(node, (list, tuple)):
        for i, v in enumerate(node):
            yield from _strings(v, f"{path}[{i}]")


def privacy_violations(profile: Any) -> list[str]:
    """Strings (and keys) that look like a coordinate, allele, genotype call or read metric."""
    found = []
    for path, text in _strings(profile):
        text = re.sub(r"HP:[0-9]{7}", "HPO", text)  # ontology identifiers are allowed
        for rx, label in _FORBIDDEN:
            if rx.search(text):
                found.append(f"{path}: looks like {label}")
    return found


def validate_profile(profile: Any) -> list[str]:
    """All problems found (schema errors plus privacy violations); an empty list means valid."""
    errors = [f"{'/'.join(map(str, e.absolute_path)) or '<root>'}: {e.message}"
              for e in sorted(_VALIDATOR.iter_errors(profile), key=lambda e: list(map(str, e.absolute_path)))]
    return errors + privacy_violations(profile)


def build_profile(profile_id: str, *, age_band: str = UNKNOWN, sex: str = UNKNOWN,
                  phenotypic_features: list[dict] | None = None, tumour_history: str = UNKNOWN,
                  tumour_current_status: str = UNKNOWN, renal: str = UNKNOWN, growth: str = UNKNOWN,
                  neurological: str = UNKNOWN, genomic_findings: list[dict] | None = None,
                  source_files: list[str] | None = None, when: str | None = None,
                  notes: list[str] | None = None) -> dict:
    """Assemble and validate a profile; every omitted field is the explicit 'unknown'. Raises ValueError."""
    seen: dict[str, str] = {}
    for f in phenotypic_features or []:
        if f["hpo_id"] in seen and seen[f["hpo_id"]] != f["status"]:
            raise ValueError(f"conflicting status for {f['hpo_id']}")
        seen[f["hpo_id"]] = f["status"]
    findings = []
    for g in genomic_findings or []:
        extra = {k: g[k] for k in ("predictor_strengths", "variant_mechanism_statements") if g.get(k)}
        findings.append({
            **extra,
            "gene_symbol": g["gene_symbol"],
            "mechanism_class": g.get("mechanism_class", UNKNOWN),
            "consequence_classes": sorted(set(g.get("consequence_classes") or [UNKNOWN])),
            "genotype_pattern": g.get("genotype_pattern", UNKNOWN),
            "phase_status": g.get("phase_status", UNKNOWN),
            "computational_evidence_strength": g.get("computational_evidence_strength", UNKNOWN),
        })
    prov: dict[str, Any] = {"source_files": list(source_files or []), "date": when or date.today().isoformat()}
    if notes:
        prov["notes"] = list(notes)
    profile = {
        "schema_version": SCHEMA_VERSION,
        "profile_id": profile_id,
        "age_band": age_band,
        "sex": sex,
        "phenotypic_features": [{"hpo_id": k, "status": v} for k, v in sorted(seen.items())],
        "tumour": {"history": tumour_history, "current_status": tumour_current_status},
        "findings": {"renal": renal, "growth": growth, "neurological": neurological},
        "genomic_findings": sorted(findings, key=lambda g: g["gene_symbol"]),
        "provenance": prov,
    }
    problems = validate_profile(profile)
    if problems:
        raise ValueError("invalid profile: " + "; ".join(problems))
    return profile


# ---------------------------------------------------------------------------------------------
# Child profile from existing Track 2 files (categorical content only)
# ---------------------------------------------------------------------------------------------

_GENOTYPE_MAP = {"compound_het": "compound heterozygous", "homozygous": "homozygous",
                 "hemizygous": "hemizygous", "structural_unpaired": "structural unpaired"}
_GROWTH_GROUPS = {"growth_restriction"}

# Panel of HPO terms used for the contrast with published cases (labels from hp.obo); written as
# explicit 'unknown' when the source records do not document them.
CONTRAST_PANEL = [
    "HP:0000252", "HP:0001511", "HP:0001510", "HP:0001263", "HP:0001250", "HP:0001305", "HP:0012443",
    "HP:0002664", "HP:0002859", "HP:0002667", "HP:0001909", "HP:0000121", "HP:0000077",
    "HP:0000518", "HP:0000568", "HP:0000478",
]


def _jload(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def calibration_summary(cal: dict) -> tuple[str, list[dict]]:
    """Categorical per-predictor strengths from the calibration file (scores are never copied)."""
    per = cal.get("per_predictor") or {}
    out = []
    for key, name in _CAL_KEYS.items():
        raw = (per.get(key) or {}).get("strength")
        out.append({"predictor": name, "strength": _STRENGTH_WORDS.get(raw, UNKNOWN) if raw else UNKNOWN})
    levels = {x["strength"] for x in out if x["strength"] != UNKNOWN}
    overall = UNKNOWN if not levels else ("discordant across predictors" if len(levels) > 1 else next(iter(levels)))
    return overall, out


def mechanism_statements(mech: dict) -> list[str]:
    """The two rule-based statements, each emitted only when the file's own rule outputs support it."""
    consts = mech.get("method_constants", {})
    out = []
    for key, v in mech.items():
        if not key.startswith("variant_") or not isinstance(v, dict):
            continue
        nmd = v.get("nmd")
        if nmd and nmd.get("nmd_predicted") is True and any(
                "kinase" in (f.get("description") or "").lower()
                for f in (v.get("truncated_product_features") or {}).get("lost_entirely", [])):
            out.append(MECH_LOF)
        cons, struct = v.get("conservation"), v.get("structure")
        if cons and struct and "kinase" in ((v.get("domain_source") or {}).get("description") or "").lower():
            frac = (cons.get("all_one2one_orthologs") or {}).get("fraction_identical_at_position_of_aligned")
            pct = struct.get("sasa_percentile_among_same_aa")
            if (frac is not None and pct is not None and frac >= consts.get("conserved_fraction_min", 0.75)
                    and pct <= consts.get("buried_sasa_percentile_max", 25.0)):
                out.append(MECH_MISSENSE)
    return sorted(set(out), key=MECH_STATEMENTS.index)


def _tri(value: Any) -> str:
    return "present" if value is True else "absent" if value is False else UNKNOWN


def build_child_profile(root: Path = Path("."), when: str | None = None) -> dict:
    """Profile of the proband from config/patient.json and the Track 2 gene-level outputs."""
    patient = _jload(root / "config/patient.json")
    diff = _jload(root / "out/differential_diagnosis.json")["genes"]
    cons = _jload(root / "out/variant_consequence_classes.json")["genes"]
    land = {g["symbol"]: g for g in _jload(root / "out/mva_landscape.json")["genes"]}
    sources = ["config/patient.json", "out/differential_diagnosis.json", "out/variant_consequence_classes.json",
               "out/mva_landscape.json"]
    calib_path = root / CALIBRATION_FILE
    calib = None
    notes = []
    if calib_path.exists():
        calib = _jload(calib_path)
        sources.append(CALIBRATION_FILE)
    else:
        notes.append("variant calibration file absent so computational evidence strength is unknown")
    mech_path = root / MECHANISM_FILE
    mech_stmts: list[str] = []
    if mech_path.exists():
        mech_stmts = mechanism_statements(_jload(mech_path))
        sources.append(MECHANISM_FILE)

    feats: dict[str, str] = {t["id"]: "present" for t in patient.get("hpo_terms", [])}
    for hpo in CONTRAST_PANEL:
        feats.setdefault(hpo, UNKNOWN)
    groups = {t.get("group") for t in patient.get("hpo_terms", [])}
    docs = patient.get("documented_findings", {})
    meas = patient.get("measurements", {})
    growth = "present" if groups & _GROWTH_GROUPS else UNKNOWN
    renal = _tri(docs.get("nephrocalcinosis")) if docs.get("nephrocalcinosis") is not None else UNKNOWN
    tumour_status = meas.get("tumour_status")

    genes = []
    for sym in sorted(diff):
        d = diff[sym]
        cc = sorted((cons.get(sym) or {}).get("class_counts", {})) or [UNKNOWN]
        pattern = _GENOTYPE_MAP.get(d.get("genotype_type") or "", UNKNOWN)
        if pattern == UNKNOWN and d.get("c3") == "structural_unpaired":
            pattern = "structural unpaired"
        reason = d.get("reason") or ""
        phase = "unresolved" if "phase not resolved" in reason else UNKNOWN
        strength, per_pred, stmts = UNKNOWN, [], []
        if sym == "BUB1B":  # the calibration and mechanism files describe BUB1B only
            if calib:
                strength, per_pred = calibration_summary(calib)
            stmts = mech_stmts
        genes.append({
            "gene_symbol": sym,
            "mechanism_class": (land.get(sym) or {}).get("mechanism_class", UNKNOWN),
            "consequence_classes": cc,
            "genotype_pattern": pattern,
            "phase_status": phase,
            "computational_evidence_strength": strength,
            "predictor_strengths": per_pred,
            "variant_mechanism_statements": stmts,
        })
    return build_profile(
        "child_proband",
        age_band=meas.get("age_band") if meas.get("age_band") in AGE_BANDS else UNKNOWN,
        sex=patient.get("genome_derived", {}).get("sex") if patient.get("genome_derived", {}).get("sex") in SEX_VALUES else UNKNOWN,
        phenotypic_features=[{"hpo_id": k, "status": v} for k, v in feats.items()],
        tumour_history={"present": "yes", "absent": "no"}.get(_tri(docs.get("tumour_history")), UNKNOWN),
        tumour_current_status=tumour_status if tumour_status in TUMOUR_CURRENT else UNKNOWN,
        renal=renal, growth=growth, neurological=UNKNOWN,  # free-text field in the source: not interpreted
        genomic_findings=genes, source_files=sources, when=when, notes=notes,
    )


def main(root: Path = Path(".")) -> Path:
    profile = build_child_profile(root)
    out = root / "out/case_profile_child.json"
    out.write_text(json.dumps(profile, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    return out


if __name__ == "__main__":
    print(main())
