"""Contrast a minimal case profile with the published BUB1B (MVA1) reference set.

Reading by fixed rule only; no inference beyond the rule. Inputs: out/case_profile_child.json (or any
profile of the same schema) and out/case_reference_bub1b.json. Output: out/case_contrast.json.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

NOT_REPORTED = "not reported in the sources"
UNKNOWN = "unknown"


def parse_frequency(token: str | None) -> float | None:
    """'93.5%' -> 0.935; '7/8' -> 0.875; anything else -> None."""
    if not token:
        return None
    m = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)%", token)
    if m:
        return float(m.group(1)) / 100
    m = re.fullmatch(r"([0-9]+)/([0-9]+)", token)
    if m and int(m.group(2)) > 0:
        return int(m.group(1)) / int(m.group(2))
    return None


def child_status(profile: dict, feature: dict) -> str:
    """present if any mapped HPO term is present; absent only if all mapped terms are absent; else unknown."""
    st = {f["hpo_id"]: f["status"] for f in profile.get("phenotypic_features", [])}
    vals = [st.get(h, UNKNOWN) for h in feature["hpo_ids"]]
    if feature["feature"] == "any malignancy" and profile.get("tumour", {}).get("history") == "yes":
        return "present"
    if feature["feature"] == "any malignancy" and profile.get("tumour", {}).get("history") == "no" and "present" not in vals:
        return "absent"
    if "present" in vals:
        return "present"
    if vals and all(v == "absent" for v in vals):
        return "absent"
    return UNKNOWN


def _primary(claims: list[dict]) -> dict | None:
    for c in claims:
        if parse_frequency(c.get("value")) is not None:
            return c
    return claims[0] if claims else None


def _shown(c: dict) -> str:
    return f"{c['value']} (PMID {c['pmid']})" if c.get("value") else f"named without frequency (PMID {c['pmid']})"


def read_row(status: str, claims: list[dict]) -> str:
    """The fixed reading rules."""
    quantified = [c for c in claims if parse_frequency(c.get("value")) is not None]
    p = _primary(claims)
    if p is None:  # nothing in the sources
        return {"present": "present in the child but not reported in published cases",
                "absent": "absent in the child; not reported in the sources",
                UNKNOWN: "not recorded in the child; not reported in the sources"}[status]
    if not quantified:
        named = f"named without a frequency in PMID {p['pmid']}"
        return {"present": f"present in the child; {named}",
                "absent": f"absent in the child; {named}",
                UNKNOWN: f"not recorded in the child; {named}"}[status]
    f = parse_frequency(p["value"])
    unit = "of published cases" if p["value"].endswith("%") else "published cases with BUB1B mutations"
    freq = f"{p['value']} {unit}, PMID {p['pmid']}"
    if status == "present":
        return f"present in the child and reported in {freq}"
    if status == "absent":
        return f"absent in the child; reported in {freq}"
    if f is not None and f > 0.5:
        return f"not recorded in the child although reported in a majority of cases ({freq}): check clinically"
    return f"not recorded in the child; reported in a minority of cases ({freq})"


def genotype_contrast(profile: dict, reference: dict, gene: str = "BUB1B") -> dict:
    """Child's genotype pattern for the disease gene versus the published patterns (class level only)."""
    gf = next((g for g in profile.get("genomic_findings", []) if g["gene_symbol"] == gene), None)
    claims = reference.get("genotype_claims", [])
    pub = [{"pattern": c["pattern"], "pmid": c["pmid"], "sentence": c["sentence"]} for c in claims]
    not_rep = reference.get("genotype_not_reported", {})
    if gf is None:
        return {"gene": gene, "child": None, "published": pub, "not_reported": not_rep,
                "reading": f"no {gene} finding recorded in the child profile"}
    classes = set(gf["consequence_classes"])
    child = {"genotype_pattern": gf["genotype_pattern"], "consequence_classes": sorted(classes),
             "phase_status": gf["phase_status"]}
    if gf["genotype_pattern"] in ("compound heterozygous", "homozygous"):
        biallelic = f"biallelic ({gf['genotype_pattern']})"
    else:
        biallelic = None
    if biallelic and classes == {"loss of function", "missense"}:
        reading = (f"child: {biallelic} with one loss-of-function and one missense class; published: biallelic "
                   "patients carry a missense variant paired with a truncating one (PMID 20516114): pattern "
                   "class matches (loss-of-function class compared with 'truncating' at class level only)")
        if gf["phase_status"] != "resolved":
            reading += "; phase not resolved, so trans configuration is not confirmed"
    elif biallelic and "loss of function" in classes and len(classes) == 1:
        reading = (f"child: {biallelic} with loss-of-function class only; published two-truncating frequency: "
                   f"{NOT_REPORTED}")
    elif biallelic:
        reading = f"child: {biallelic} with classes {sorted(classes)}; pattern not described in the sources beyond biallelic"
    else:
        reading = f"child: {gf['genotype_pattern']}; published patterns are biallelic (PMID 35804254): no match by rule"
    reading += ("; domain of the missense variant is not stored in the profile (categorical level), "
                "so the kinase-domain location reported in PMID 20516114 is not assessed")
    return {"gene": gene, "child": child, "published": pub, "not_reported": not_rep, "reading": reading}


def build_contrast(profile: dict, reference: dict) -> dict:
    rows = []
    for f in reference["features"]:
        status = child_status(profile, f)
        claims = f["claims"]
        rows.append({
            "feature": f["feature"],
            "hpo_ids": f["hpo_ids"],
            "published_frequency": ([{"value": c.get("value"), "pmid": c["pmid"], "population": c["population"],
                                      "sentence": c["sentence"]} for c in claims] if claims else NOT_REPORTED),
            "child": status,
            "reading": read_row(status, claims),
        })
    return {
        "profile_id": profile["profile_id"],
        "disease": reference["disease"],
        "reference_retrieved": reference["retrieved"],
        "rows": rows,
        "genotype": genotype_contrast(profile, reference),
        "rules": [
            "present in the child and reported in X of published cases",
            "not recorded in the child although reported in a majority of cases: check clinically (majority: primary frequency above 50%)",
            "present in the child but not reported in published cases",
            "other combinations are stated literally (absent, minority, named without frequency, not reported)",
        ],
        "limits": "Published frequencies come from reviews of case reports (selection bias); the contrast lists what to check, not a diagnosis.",
    }


def main(root: Path = Path(".")) -> Path:
    profile = json.loads((root / "out/case_profile_child.json").read_text(encoding="utf-8"))
    reference = json.loads((root / "out/case_reference_bub1b.json").read_text(encoding="utf-8"))
    out = root / "out/case_contrast.json"
    out.write_text(json.dumps(build_contrast(profile, reference), indent=1, sort_keys=True, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    return out


if __name__ == "__main__":
    print(main())
