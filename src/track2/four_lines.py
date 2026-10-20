"""Combine the four lines of evidence into candidates (rule: `protocol/addendum_four_lines.md`, section "Combining the lines").

    python -m track2.four_lines --out out

Inputs (all in `out/`): kg_validation.json and kg_query.json (A), phenotype_validation.json, phenotype_proximity.json and core_sensitivity.json (B), vulnerability.json (C),
signed_effect.json (D). A line counts only if its own validation passed. A drug is a *candidate* if it is FDA-eligible (openFDA), not topical-only, not an acute-use anaesthetic
or sedative (ATC N01, N05C), and supported by at least two of the usable lines A, B, C, or (amendment 3) by one usable line plus a literature rationale accepted in `config/literature_rationale.json`. Tier A: candidate and D usable and favourable; tier B: candidate and D unclear or not usable
(exploratory); *avoid*: candidate and D usable and unfavourable. Tier C: supported by one usable line only (listed, not proposed). Writes `out/four_lines.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from scipy.stats import hypergeom

from track2.fda_approved import lookup

TOP_A = 50


def load(out: Path, name: str) -> dict[str, Any]:
    return json.loads((out / name).read_text(encoding="utf-8"))


def line_status(out: Path) -> dict[str, dict[str, Any]]:
    kg, ph, vu, sg = load(out, "kg_validation.json"), load(out, "phenotype_validation.json"), load(out, "vulnerability.json"), load(out, "signed_effect.json")
    return {
        "A": {"usable": bool(kg["validated"]), "evidence": f"mean AUC {kg['summary']['mean_auc']:.3f}, 95% interval {kg['summary']['ci95_low']:.3f}-{kg['summary']['ci95_high']:.3f}, {kg['n_evaluated']} diseases"},
        "B": {"usable": bool(ph["primary"]["validated"]), "evidence": f"mean AUC {ph['primary']['summary']['mean_auc']:.3f}, 95% interval {ph['primary']['summary']['ci95_low']:.3f}-{ph['primary']['summary']['ci95_high']:.3f}, {ph['primary']['n_evaluable']} diseases"},
        "C": {"usable": bool(vu["line_c_used"]), "evidence": f"observed {vu['observed_total']} hits, null 95th percentile {vu['null']['total_p95']:.1f}"},
        "D": {"usable": bool(sg["validation"]["used"]), "evidence": f"{sg['validation']['n_recovered']} of 6 textbook perturbations recovered"},
    }


def support_sets(out: Path) -> dict[str, set[str]]:
    a = {d for d, r in load(out, "kg_query.json")["rank_combined"].items() if r <= TOP_A}
    pheno = {d["drug_id"] for d in load(out, "phenotype_proximity.json")["drugs"] if d["robust"]}
    core = {d["drug_id"] for d in load(out, "core_sensitivity.json")["drugs"] if d["L1_core"]}
    vu = load(out, "vulnerability.json")
    return {"A": a, "B": pheno | core, "B_phenotype_module": pheno, "B_mitosis_module_without_generic_complexes": core, "C": set(vu["supported_drugs"])}


def load_rationale(path: Path) -> set[str]:
    """Names of the drugs whose literature rationale was accepted by the author (`config/literature_rationale.json`)."""
    return {x["drug"].upper() for x in json.loads(path.read_text(encoding="utf-8"))["accepted"]}


def gates(drug: dict[str, Any], fda: dict[str, Any]) -> dict[str, bool]:
    atc = drug.get("atc") or []
    topical_only = bool(drug.get("topical")) and not (drug.get("oral") or drug.get("parenteral"))
    return {"fda_eligible": bool(fda["eligible"]), "systemic_route": not topical_only,
            "not_acute_use_class": not any(a.startswith(("N01", "N05C")) for a in atc)}


def combine(universe: dict[str, Any], status: dict[str, dict[str, Any]], sets: dict[str, set[str]], d_calls: dict[str, dict[str, Any]],
            fda_fn=lookup, literature: set[str] | None = None) -> list[dict[str, Any]]:
    by_id = {d["drug_id"]: d for d in universe["drugs"]}
    usable = [k for k in ("A", "B", "C") if status[k]["usable"]]
    rows = []
    for drug_id in sorted(set().union(*(sets[k] for k in usable)) if usable else set()):
        drug = by_id[drug_id]
        lines = [k for k in usable if drug_id in sets[k]]
        g = gates(drug, fda_fn(drug["name"]))
        d_call = (d_calls.get(drug_id) or {}).get("call", "unclear")
        passed = all(g.values())
        with_literature = len(lines) == 1 and drug["name"].upper() in (literature or set())  # amendment 3: one usable line plus a literature rationale
        evidence = "two_or_more_lines" if len(lines) >= 2 else "one_line_plus_literature" if with_literature else "one_line"
        candidate = passed and (len(lines) >= 2 or with_literature)
        if candidate:
            tier = ("A" if d_call == "strengthens_checkpoint" else "avoid" if d_call == "weakens_checkpoint" else "B") if status["D"]["usable"] else "B"
        else:
            tier = "C" if passed else "not_eligible"
        rows.append({"drug_id": drug_id, "name": drug["name"], "supported_by": lines, "n_lines": len(lines), "evidence": evidence, "gates": g, "candidate": candidate, "tier": tier,
                     "list": "tumour_context" if any(a.startswith("L01") for a in drug.get("atc") or []) else "constitutional_defect",
                     "D_direction_exploratory" if not status["D"]["usable"] else "D_direction": d_call,
                     "B_modules": [m for m in ("B_phenotype_module", "B_mitosis_module_without_generic_complexes") if drug_id in sets[m]] if "B" in lines else []})
    return sorted(rows, key=lambda r: (-r["n_lines"], r["tier"], r["name"]))


def overlap_statistics(sets: dict[str, set[str]], ranked: set[str]) -> dict[str, Any]:
    """Hypergeometric test of the overlap between lines among the drugs that line A ranked (the common universe)."""
    out = {}
    a, b, c = sets["A"] & ranked, sets["B"] & ranked, sets["C"] & ranked
    for name, (x, y) in {"A_and_B": (a, b), "A_and_C": (a, c), "B_and_C": (b, c)}.items():
        k = len(x & y)
        out[name] = {"overlap": k, "size_first": len(x), "size_second": len(y), "universe": len(ranked), "p_hypergeometric": float(hypergeom.sf(k - 1, len(ranked), len(x), len(y))) if k else 1.0}
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    args = parser.parse_args()
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    status, sets = line_status(args.out), support_sets(args.out)
    d_calls = load(args.out, "signed_effect.json")["drugs"]
    rationale = load_rationale(Path("config/literature_rationale.json"))
    rows = combine(universe, status, sets, d_calls, literature=rationale)
    ranked = set(load(args.out, "kg_query.json")["rank_combined"])
    result = {"lines": status, "n_usable_lines": sum(1 for k in ("A", "B", "C") if status[k]["usable"]), "support_set_sizes": {k: len(v) for k, v in sets.items()},
              "overlap_vs_chance": overlap_statistics(sets, ranked), "candidates": [r for r in rows if r["candidate"]],
              "tier_c_one_line": [r for r in rows if r["tier"] == "C"], "not_eligible": [r["name"] for r in rows if r["tier"] == "not_eligible"]}
    result["negative_result"] = result["n_usable_lines"] < 2
    (args.out / "four_lines.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print({k: v["usable"] for k, v in status.items()}, "| candidates:", [(r["name"], r["supported_by"], r["tier"]) for r in result["candidates"]])
    print("overlap vs chance:", {k: (v["overlap"], round(v["p_hypergeometric"], 3)) for k, v in result["overlap_vs_chance"].items()})


if __name__ == "__main__":
    main()
