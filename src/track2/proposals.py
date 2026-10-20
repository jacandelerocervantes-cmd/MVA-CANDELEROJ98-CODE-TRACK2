"""The table of proposed repositioning candidates (one row per drug, FDA-approved only), built from the plausibility tiers.

    python -m track2.proposals --variant e0_e1_e2 --out out

Rows come from two places, kept apart by `arm`:
- `main`: Tier 1 and Tier 2 drugs of `out/plausibility_<variant>.json` (FDA-eligible, evidence lines L1 and L2, see `protocol/addendum_plausibility.md`);
- `extrapolation`: drugs of the separate extrapolation arm (`out/arm_extrapolation/candidates.json`) whose mechanism matches a literature class that is a
  mechanistic rationale (mTOR pathway, senescence), also FDA-eligible. They have no proximity tier and are labelled as extrapolations from mouse or
  cancer-cell evidence.
Each row carries what the drug is approved for today, the gene products it acts on and how, the evidence lines, the three age groups
(0-<3, 3-<12, 12 and older) and the cautions shown but not used to exclude. Hypotheses for follow-up; not medical advice and not dosing.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Callable

from track2.age_groups import GROUPS, group_status
from track2.drug_context import approved_indications
from track2.fda_approved import lookup
from track2.literature_context import classify
from track2.pediatric_review import load_review, reviewed_profile

EXTRAPOLATION_CLASSES = {"mtor_pathway", "senescence_axis"}


def core_summary(core: dict[str, Any] | None, drug_id: str, l1: bool) -> dict[str, Any]:
    """L1 without the generic complexes (`protocol/addendum_core.md`): the count of the nine configurations with the core modules and the label."""
    if not core:
        return {"L1_core": None, "label": "not_assessed"}
    row = next((d for d in core["drugs"] if d["drug_id"] == drug_id), None)
    ok = bool(row and row["L1_core"])
    return {"L1_core": ok, "robustness_count_core": row["robustness_count_core"] if row else 0,
            "label": "proximity_by_generic_complex_only" if l1 and not ok else ("proximity_also_without_generic_complexes" if ok else "no_proximity")}


def l4_summary(prism: dict[str, Any] | None, drug_id: str) -> dict[str, Any]:
    """L4 from the PRISM analysis (`out/prism_l4.json`): the call and the strongest row (viability difference, adjusted FDR)."""
    entry = ((prism or {}).get("mean_all") or {}).get("drugs", {}).get(drug_id)
    if not entry:
        return {"call": "not_assessed"}
    rows = entry["rows"]
    best = min(rows, key=lambda t: t["fdr"]) if rows else None
    return {"call": entry["call"], "diff_weak_minus_strong": round(best["diff_weak_minus_strong"], 2) if best else None,
            "fdr_adjusted": round(best["fdr_adjusted"], 4) if best else None, "bub1b_only_call": ((prism or {}).get("BUB1B") or {}).get("drugs", {}).get(drug_id, {}).get("call")}


def main_rows(plaus: dict[str, Any], prism: dict[str, Any] | None = None, core: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    rows = []
    for r in plaus["drugs"]:
        if r["tier"] not in (1, 2):
            continue
        rows.append({
            "arm": "main", "drug_id": r["drug_id"], "name": r["name"], "list": r["list"], "tier": r["tier"], "novelty": r["novelty"],
            "approved_for": r["approved_indications"], "acts_on": [{"gene": x["gene"], "role": x["role"]} for x in r["target_roles"]],
            "hypothesis": r["hypothesis"],
            "evidence": {"L1_robust_proximity": r["L1_robust_proximity"], "L2_binding": r["L2_binding_curated_or_measured"], "L3_rationale": r["L3_mechanistic_rationale"],
                         "L3_classes": r["L3_classes"], "L4": l4_summary(prism, r["drug_id"]), "core": core_summary(core, r["drug_id"], r["L1_robust_proximity"]), "robustness_count_of_9": r["robustness_count"], "expanded_best_z": r["expanded_best_z"]},
            "age_groups": r["age_groups"], "age_bands": r["age_bands"], "label_min_age_years": r["label_min_age_years"], "label_source": r["label_source"],
            "cautions": {"boxed_warning": r["boxed_warning"], "label_caution": r.get("label_caution"), "broad_target_profile": r.get("broad_target_profile", False), "n_targets": r.get("n_targets_expanded"),
                         "excipient_alerts": r["excipient_alerts"]},
            "fda": r["fda"],
        })
    return rows


def extrapolation_rows(arm: list[dict[str, Any]], config: dict[str, Any], review: dict[str, Any], fda_fn: Callable[[str], dict[str, Any]] = lookup,
                       indications_fn: Callable[[str], list[str]] = approved_indications, routes: dict[str, dict[str, Any]] | None = None,
                       prism: dict[str, Any] | None = None, direct_hits: dict[str, bool] | None = None) -> list[dict[str, Any]]:
    rows = []
    for c in arm:
        route = (routes or {}).get(c["drug_id"]) or {}
        if route.get("topical") and not (route.get("oral") or route.get("parenteral")):
            continue  # same gate as the main arm: topical-only drugs are out
        classes = [e for e in classify(c, config) if e["class"] in EXTRAPOLATION_CLASSES]
        if not classes:
            continue
        fda = fda_fn(c["name"] or "")
        if not fda["eligible"]:
            continue
        indications = indications_fn(c["drug_id"])
        entry = review.get(c["drug_id"])
        profile, source = reviewed_profile(c["age_profile"], entry)
        min_age = entry["min_age_years"] if entry and entry["verdict"] == "established" else c["age_profile"].get("min_age_years")
        rows.append({
            "arm": "extrapolation", "drug_id": c["drug_id"], "name": c["name"], "list": "extrapolation", "tier": None, "novelty": "extrapolation_from_mouse_or_cancer_cell_evidence",
            "approved_for": indications, "acts_on": [{"gene": m, "role": f"{m} (curated mechanism)"} for m in c["mechanism_of_action"]],
            "hypothesis": f"{c['name'].title()}: {'; '.join(c['mechanism_of_action'])}. Rationale class: {', '.join(e['class'] for e in classes)} (extrapolation; no drug was tested in the source mouse studies).",
            "evidence": {"L1_robust_proximity": False, "L2_binding": bool((direct_hits or {}).get(c["drug_id"], False)), "L3_rationale": bool(classes), "L3_classes": [e["class"] for e in classes], "L4": l4_summary(prism, c["drug_id"]),
                         "pmids": sorted({p for e in classes for p in e["pmids"]})},
            "age_groups": group_status(profile, min_age, (entry or {}).get("verdict")), "age_bands": {b: profile.get(b, "UNKNOWN") for b in ("neonate", "infant", "child", "adolescent", "adult")},
            "label_min_age_years": min_age, "label_source": source,
            "cautions": {"boxed_warning": any(f == "BOXED_WARNING" for f in c["risk_flags"]), "broad_target_profile": False, "n_targets": len(c.get("targets", [])),
                         "excipient_alerts": [a["excipient"] for a in c["excipients"]["alerts"]]},
            "fda": {k: fda[k] for k in ("matched_name", "marketing_statuses", "applications", "n_applications")},
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=("e0_e1", "e0_e1_e2"), required=True)
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    plaus = json.loads((args.out / f"plausibility_{args.variant}.json").read_text(encoding="utf-8"))
    arm = json.loads((args.out / "arm_extrapolation" / "candidates.json").read_text(encoding="utf-8"))
    config = json.loads(Path("config/literature_context.json").read_text(encoding="utf-8"))
    universe = json.loads(Path("data/universe.json").read_text(encoding="utf-8"))
    routes = {d["drug_id"]: {"oral": d.get("oral"), "parenteral": d.get("parenteral"), "topical": d.get("topical")} for d in universe["drugs"]}
    prism_path = args.out / "prism_l4.json"
    prism = json.loads(prism_path.read_text(encoding="utf-8")) if prism_path.exists() else None
    core_path = args.out / "core_sensitivity.json"
    core = json.loads(core_path.read_text(encoding="utf-8")) if core_path.exists() else None
    direct: dict[str, bool] = {}
    for mid in ("bub1b_down_proteotoxic", "bub1b_down_p53_senescence", "bub1b_down_growth_mtor_ampk"):
        sp = args.out / "comparison_manual" / f"scores_{mid}.json"
        if sp.exists():
            for r in json.loads(sp.read_text(encoding="utf-8"))["rows"]:
                direct[r["drug_id"]] = direct.get(r["drug_id"], False) or bool(r.get("direct_module_targets"))
    rows = main_rows(plaus, prism, core) + extrapolation_rows(arm, config, load_review(), routes=routes, prism=prism, direct_hits=direct)
    result = {"variant": args.variant, "n_main": sum(1 for r in rows if r["arm"] == "main"), "n_extrapolation": sum(1 for r in rows if r["arm"] == "extrapolation"),
              "age_groups": list(GROUPS), "note": "FDA-approved drugs only. Hypotheses for follow-up; not medical advice and not dosing.", "proposals": rows}
    (args.out / "proposals.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(args.out / "proposals.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["arm", "list", "tier", "drug_id", "name", "approved_for", "acts_on", "novelty", *[f"age_{g}" for g in GROUPS], "boxed_warning", "broad_target_profile"])
        for r in rows:
            writer.writerow([r["arm"], r["list"], r["tier"] or "", r["drug_id"], r["name"], "; ".join(r["approved_for"][:4]), "; ".join(x["role"] for x in r["acts_on"][:3]),
                             r["novelty"], *[r["age_groups"][g] for g in GROUPS], r["cautions"]["boxed_warning"], r["cautions"]["broad_target_profile"]])
    print({"n_main": result["n_main"], "n_extrapolation": result["n_extrapolation"]})


if __name__ == "__main__":
    main()
