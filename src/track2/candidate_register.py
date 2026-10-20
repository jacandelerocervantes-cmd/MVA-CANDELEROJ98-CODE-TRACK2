"""The complete, justified register of candidates (rule: `protocol/addendum_wide.md`).

    python -m track2.candidate_register --out out

One row per FDA-approved drug that was evaluated and reached at least one evidence line in any arm of the analysis: the strict plausibility tiers, the wider ring, the
extrapolation arm and the paediatric search. Each row says which lines it satisfies and which it fails, its PRISM call, what it is approved for today, which gene products
it acts on, the three age groups and the cautions, and a one-sentence justification generated from those fields. The report's tables are views of this register.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from track2.age_groups import GROUPS
from track2 import repositioning
from track2.proposals import core_summary, l4_summary

L4_TEXT = {"supports_resistance_when_SAC_weak": "PRISM: cells with a weak checkpoint are more resistant", "supports_sensitivity_when_SAC_weak": "PRISM: cells with a weak checkpoint are more sensitive",
           "no_association": "PRISM: no association with checkpoint strength", "lineage_confounded": "PRISM: association explained by cell lineage", "mixed": "PRISM: mixed",
           "not_in_PRISM": "not in PRISM", "not_assessed": "PRISM not assessed", "not_testable": "PRISM not testable"}


def justification(r: dict[str, Any]) -> str:
    ok = [n for n, v in (("L1 robust proximity", r["L1"]), ("L2 binding to a module protein", r["L2"]), ("L3 literature rationale", r["L3"])) if v]
    no = [n for n, v in (("L1 robust proximity", r["L1"]), ("L2 binding to a module protein", r["L2"]), ("L3 literature rationale", r["L3"])) if not v]
    text = f"Satisfies {', '.join(ok) or 'no evidence line'}; does not satisfy {', '.join(no) or 'none'}; {L4_TEXT.get(r['L4_call'], r['L4_call'])}."
    if r.get("core", {}).get("label") == "proximity_by_generic_complex_only":
        text += " Its proximity comes only from membership of generic complexes (proteasome or tubulin): it disappears when they are removed from the modules."
    if r["novelty"] == "new_via_expanded_links":
        text += " New: close to a module only thanks to measured or predicted links."
    if r["broad_target_profile"]:
        text += " Broad target profile: closeness may reflect breadth."
    return text


def build(out: Path) -> dict[str, Any]:
    prism = json.loads((out / "prism_l4.json").read_text(encoding="utf-8")) if (out / "prism_l4.json").exists() else None
    core = json.loads((out / "core_sensitivity.json").read_text(encoding="utf-8")) if (out / "core_sensitivity.json").exists() else None
    rows: dict[str, dict[str, Any]] = {}

    def add(source: str, r: dict[str, Any], tier_key: str | None) -> None:
        row = rows.setdefault(r["drug_id"], {"drug_id": r["drug_id"], "name": r["name"], "sources": [], "tier_strict": None, "tier_wide": None})
        row["sources"].append(source)
        if tier_key:
            row[tier_key] = r["tier"] or None
        row.setdefault("list", r["list"])
        for k in ("L1", "L2", "L3"):
            row.setdefault(k, False)
        row["L1"] = row["L1"] or r.get("L1_robust_proximity", False)
        row["L2"] = row["L2"] or r.get("L2_binding_curated_or_measured", False)
        row["L3"] = row["L3"] or r.get("L3_mechanistic_rationale", False)
        row.setdefault("novelty", r.get("novelty", ""))
        for k in ("approved_indications", "target_roles", "age_groups", "age_bands", "label_caution", "boxed_warning", "broad_target_profile", "gates", "excipient_alerts", "mechanism_of_action", "proximal_targets"):
            if k in r and k not in row:
                row[k] = r[k]
        row.setdefault("fda_eligible", (r.get("gates") or {}).get("fda_eligible"))

    for name, tag, key in (("plausibility_e0_e1_e2.json", "strict", "tier_strict"), ("plausibility_e0_e1_e2_wide.json", "wide_ring", "tier_wide")):
        path = out / name
        if path.exists():
            for r in json.loads(path.read_text(encoding="utf-8"))["drugs"]:
                if all(r["gates"].values()) and (r["L1_robust_proximity"] or r["L2_binding_curated_or_measured"] or r["L3_mechanistic_rationale"]):
                    add(tag, r, key)
    props = out / "proposals.json"
    if props.exists():
        for r in json.loads(props.read_text(encoding="utf-8"))["proposals"]:
            if r["arm"] == "extrapolation":
                ev = r["evidence"]
                add("extrapolation_arm", {**r, "L1_robust_proximity": False, "L2_binding_curated_or_measured": ev["L2_binding"], "L3_mechanistic_rationale": ev["L3_rationale"],
                                          "gates": {"fda_eligible": True}, "approved_indications": r["approved_for"], "target_roles": r["acts_on"], "boxed_warning": r["cautions"]["boxed_warning"],
                                          "broad_target_profile": False, "excipient_alerts": r["cautions"]["excipient_alerts"], "label_caution": r["cautions"].get("label_caution"),
                                          "mechanism_of_action": [x["gene"] for x in r["acts_on"]]}, None)
    ped = out / "pediatric_arm.json"
    if ped.exists():
        for r in json.loads(ped.read_text(encoding="utf-8"))["drugs"]:
            if r["headline"] and r["drug_id"] in rows:
                rows[r["drug_id"]]["sources"].append("paediatric_search")
    out_rows = []
    config = repositioning.load()
    for d, row in rows.items():
        row["repositioning"] = repositioning.proposed_use(row.get("mechanism_of_action") or [], [p["gene"] for p in row.get("proximal_targets") or []], row.get("novelty", ""), config, row.get("list"))
        row["L4_call"] = l4_summary(prism, d)["call"]
        row["core"] = core_summary(core, d, row["L1"])
        row["sources"] = sorted(set(row["sources"]))
        row["age_groups"] = row.get("age_groups") or {g: "verify" for g in GROUPS}
        row["age_bands"] = row.get("age_bands") or {b: "UNKNOWN" for b in ("neonate", "infant", "child", "adolescent", "adult")}
        row["justification"] = justification(row)
        out_rows.append(row)
    out_rows.sort(key=lambda r: (r["tier_strict"] or 9, r["tier_wide"] or 9, r["name"]))
    return {"n_drugs": len(out_rows), "by_source": {s: sum(1 for r in out_rows if s in r["sources"]) for s in ("strict", "wide_ring", "extrapolation_arm", "paediatric_search")},
            "note": "FDA-approved drugs only. Hypotheses for follow-up; not medical advice and not dosing.", "drugs": out_rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    reg = build(args.out)
    (args.out / "candidate_register.json").write_text(json.dumps(reg, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(args.out / "candidate_register.csv", "w", encoding="utf-8", newline="") as handle:
        w = csv.writer(handle, lineterminator="\n")
        w.writerow(["drug_id", "name", "list", "tier_strict", "tier_wide", "sources", "L1", "L1_core", "L2", "L3", "L4_PRISM", "approved_for", "proposed_repositioned_use", "caveat", "how_to_validate_or_refute", "acts_on", *[f"age_{g}" for g in GROUPS], "boxed_warning", "justification"])
        for r in reg["drugs"]:
            w.writerow([r["drug_id"], r["name"], r["list"], r["tier_strict"] or "", r["tier_wide"] or "", ";".join(r["sources"]), r["L1"], (r.get("core") or {}).get("L1_core", ""), r["L2"], r["L3"], r["L4_call"],
                        "; ".join((r.get("approved_indications") or [])[:4]), r["repositioning"]["new_use"], r["repositioning"]["caveat"], r["repositioning"].get("validation", ""), "; ".join((x.get("role") or x.get("gene", "")) for x in (r.get("target_roles") or [])[:3]),
                        *[r["age_groups"][g] for g in GROUPS], r.get("boxed_warning", ""), r["justification"]])
    print(reg["n_drugs"], "drugs;", reg["by_source"])


if __name__ == "__main__":
    main()
