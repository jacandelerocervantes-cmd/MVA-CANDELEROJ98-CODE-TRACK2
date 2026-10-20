"""Repurposing-plausibility framework: evidence lines, feasibility gates and tiers (rules: `protocol/addendum_plausibility.md`).

    python -m track2.plausibility --variant e0_e1 --out out

Reads `out/novel_hypotheses_<variant>.json` (the pool of drugs in the top 20 of at least 5 configurations plus the new ones). For each
drug: L1 robust proximity, L2 binding to a module protein (curated or measured; predicted is labelled), L3 mechanistic rationale (verified
literature class), L4 not assessed; gates FDA-eligible, systemic route, not an acute-use anaesthetic or sedative; antineoplastics go to a separate
list. Tier 1 = gates + L1 + L2 + L3; Tier 2 = gates + L1 + L2; Tier 3 = gates + L1 or L2. Not a recommendation and not dosing advice.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Callable

from track2.age_groups import GROUPS, group_status, summarise as summarise_groups
from track2.drug_context import approved_indications, hypothesis_sentence, is_tumour_indication, target_roles
from track2.fda_approved import lookup

L1_MIN = 7
MIN_MEASUREMENTS = 2  # amendment after the independent review: single-measurement links are shown but do not count as L2
RATIONALE_EXCLUDED = {"mva_tumour_case_reports", "microtubule_targeting"}  # not a mechanistic rationale / direction is unfavourable-or-unknown
BANDS = ("neonate", "infant", "child", "adolescent", "adult")


def evaluate(row: dict[str, Any], fda: dict[str, Any], indications: list[str] | None = None) -> dict[str, Any]:
    all_links = [x for v in row["links"].values() for x in v]
    sources = {x["source"] for x in all_links}
    l1 = row["robustness_count"] >= L1_MIN
    # A measured link counts only with at least MIN_MEASUREMENTS independent activity records (a single record may come from one screen or one document)
    measured_ok = any(x["source"] == "measured" and (x.get("n_measurements") or 0) >= MIN_MEASUREMENTS for x in all_links)
    l2 = "curated" in sources or measured_ok
    single_only = "measured" in sources and not l2
    l2_pred = "predicted" in sources
    rationale = sorted({e["class"] for e in row.get("literature_context", []) if e["class"] not in RATIONALE_EXCLUDED})
    l3 = bool(rationale)
    topical_only = bool(row.get("topical")) and not (row.get("oral") or row.get("parenteral"))  # ChEMBL route flags are incomplete: only topical-only is excluded
    tumour = row["class"] == "antineoplastic_tumour_context" or is_tumour_indication(indications or [])
    gates = {"fda_eligible": fda["eligible"], "systemic_route": not topical_only,
             "not_acute_use_class": row["class"] != "excluded_acute_use_anaesthetic_sedative"}
    passed = all(gates.values())
    tier = 1 if passed and l1 and l2 and l3 else 2 if passed and l1 and l2 else 3 if passed and (l1 or l2) else 0
    lab = row.get("label") or {}
    roles = target_roles(row["links"])
    return {
        "approved_indications": indications or [], "target_roles": roles, "hypothesis": hypothesis_sentence(row["name"] or "", roles, indications or [], row.get("proximal_targets")),
        "drug_id": row["drug_id"], "name": row["name"], "tier": tier, "list": "tumour_context" if tumour else "constitutional_defect",
        "novelty": "new_via_expanded_links" if row["is_new"] else "curated_mechanism_already_near",
        "L1_robust_proximity": l1, "L2_binding_curated_or_measured": l2, "L2_predicted_only": l2_pred and not l2, "L2_single_measurement_only": single_only, "L3_mechanistic_rationale": l3,
        "L3_classes": rationale, "L4_independent_functional_data": "not_assessed", "gates": gates,
        "fda": {k: fda[k] for k in ("matched_name", "marketing_statuses", "applications", "n_applications")},
        "robustness_count": row["robustness_count"], "expanded_best_z": row["expanded_best_z"], "curated_best_z": row["curated_best_z"],
        "n_module_proteins_linked": row["n_module_proteins_linked"], "promiscuous": row["promiscuous"], "n_targets_expanded": row.get("n_targets_expanded"), "broad_target_profile": row.get("broad_target_profile", False), "proximal_targets": row.get("proximal_targets", []), "boxed_warning": row.get("black_box_warning", False),
        "age_bands": {b: (lab.get("bands") or {}).get(b, "UNKNOWN") for b in BANDS}, "label_source": lab.get("source"),
        "age_groups": group_status(lab.get("bands") or {}, lab.get("min_age_years"), lab.get("verdict")) if lab else {g: "verify" for g in GROUPS},
        "label_min_age_years": lab.get("min_age_years"), "label_caution": lab.get("caution"),
        "excipient_alerts": lab.get("excipient_alerts", []), "mechanism_of_action": row["mechanism_of_action"], "links": row["links"],
    }


def by_age_band(rows: list[dict[str, Any]]) -> dict[str, dict[str, list[str]]]:
    """For each band, the Tier 1 and Tier 2 drugs grouped by what the US label says for that band."""
    out: dict[str, dict[str, list[str]]] = {}
    for band in BANDS:
        groups: dict[str, list[str]] = {"labelled_for_this_band": [], "no_labelled_use": [], "unknown_or_verify": []}
        for r in rows:
            if r["tier"] not in (1, 2):
                continue
            status = r["age_bands"][band]
            key = "labelled_for_this_band" if status in ("LABELLED", "LABELLED_USE") else "no_labelled_use" if status == "NO_LABELLED_USE" else "unknown_or_verify"
            groups[key].append(r["name"])
        out[band] = {k: sorted(v) for k, v in groups.items()}
    return out


def run(pool: dict[str, Any], lookup_fn: Callable[[str], dict[str, Any]] = lookup,
        indications_fn: Callable[[str], list[str]] = approved_indications) -> dict[str, Any]:
    rows = [evaluate(r, lookup_fn(r["name"] or ""), indications_fn(r["drug_id"])) for r in pool["drugs"]]
    rows.sort(key=lambda r: (r["tier"] == 0, r["tier"], -r["robustness_count"], r["expanded_best_z"], r["drug_id"]))
    tiers = {f"tier_{t}": sum(1 for r in rows if r["tier"] == t) for t in (1, 2, 3)}
    by_list = {lst: {f"tier_{t}": sum(1 for r in rows if r["tier"] == t and r["list"] == lst) for t in (1, 2, 3)} for lst in ("constitutional_defect", "tumour_context")}
    return {"summary": {**pool["summary"], **tiers, "by_list": by_list, "n_not_fda_eligible": sum(1 for r in rows if not r["gates"]["fda_eligible"]),
                        "L4": "not assessed in this version"}, "by_age_band": by_age_band(rows), "by_age_group": summarise_groups([r for r in rows if r["tier"] in (1, 2, 3)]), "drugs": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=("e0_e1", "e0_e1_e2"), required=True)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--tag", default="", help="suffix of the pool and output files, for example _wide")
    args = parser.parse_args()
    pool = json.loads((args.out / f"novel_hypotheses_{args.variant}{args.tag}.json").read_text(encoding="utf-8"))
    result = run(pool)
    (args.out / f"plausibility_{args.variant}{args.tag}.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(args.out / f"plausibility_{args.variant}{args.tag}.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["tier", "list", "drug_id", "name", "approved_indications", "genes_acted_on", "novelty", "robustness_count", "L1", "L2", "L2_predicted_only", "L3", "fda_eligible", "boxed_warning", "promiscuous", *[f"age_{b}" for b in BANDS], *[f"group_{g}" for g in GROUPS]])
        for r in result["drugs"]:
            writer.writerow([r["tier"], r["list"], r["drug_id"], r["name"], "; ".join(r["approved_indications"]), "; ".join(x["role"] for x in r["target_roles"][:3]), r["novelty"], r["robustness_count"], r["L1_robust_proximity"], r["L2_binding_curated_or_measured"],
                             r["L2_predicted_only"], r["L3_mechanistic_rationale"], r["gates"]["fda_eligible"], r["boxed_warning"], r["promiscuous"], *[r["age_bands"][b] for b in BANDS], *[r["age_groups"][g] for g in GROUPS]])
    print(result["summary"])


if __name__ == "__main__":
    main()
