"""Combine the lines of evidence under addendum 11 (`protocol/addendum_review2.md`: E1, E3, E4, E7, E9, E14).

    python -m track2.four_lines_review2 --out out

Rules applied (all fixed before the runs):
- A is usable if kg_validation passed and the bootstrap lower bound of its excess over the degree baseline is > 0 (E3). A drug is supported by A if it is among the 50 best *distinct* drugs of
  the patient query (PrimeKG nodes that map to one drug are merged, best rank kept) and it appears in the top 50 of fewer than 5 % of the 500 degree-matched random queries (E4); a drug whose
  specificity was not measured does not count and is listed.
- B is one line. B1 (phenotype module, validated without the causal gene and with excess over degree > 0, E2/E3) and B2 (mitosis module of BUB1B alone, W3, links with >= 2 records, E7/E8; usable only
  if `module_rule_validation` passed, E1). B supports a drug through whichever of B1 or B2 is usable.
- C is usable if the adjusted total exceeds the 95th percentile of its 1,000-set null (E6).
- D is exploratory: it gives a direction note and no tier (E7).
- A candidate is an approved (E14: FDA-eligible and ChEMBL max_phase 4, the universe), systemic, non-N01/N05C drug supported by >= 2 usable lines among A, B, C (E9). One line plus literature is listed apart.
- The overlap between lines is tested on the drugs that both lines can score (E9).
B2 passed its pre-specified validation but not the check that leaves out the drugs acting on the causal genes (not testable); the result is therefore reported twice: with B2 counted as
pre-specified and with B2 held as exploratory. Writes `out/four_lines_review2.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from scipy.stats import hypergeom

from track2.fda_approved import lookup
from track2.four_lines import gates, load_rationale

TOP_A, GENERIC_FRACTION = 50, 0.05


def load(out: Path, name: str) -> Any:
    return json.loads((out / name).read_text(encoding="utf-8"))


def distinct_top(rank_combined: dict[str, int], n: int) -> list[str]:
    """Drug ids by best rank: several graph nodes map to one drug, so the top n has n distinct drugs."""
    best: dict[str, int] = {}
    for drug, r in rank_combined.items():
        best[drug] = min(r, best.get(drug, r))
    return sorted(best, key=lambda d: (best[d], d))[:n]


def line_status(out: Path) -> dict[str, dict[str, Any]]:
    kg, extra = load(out, "kg_validation.json"), load(out, "kg_extra.json")["results"]["disease_only"]
    e2, mrv = load(out, "phenotype_validation_e2.json"), load(out, "module_rule_validation.json")["evaluation"]
    mrv_x = load(out, "module_rule_validation.json").get("evaluation_without_drugs_acting_on_causal_genes")
    vu = load(out, "vulnerability_e6.json") if (out / "vulnerability_e6.json").exists() else load(out, "c1000_e6/vulnerability_e6.json")
    sg = load(out, "signed_effect.json")
    a_ok = bool(kg["validated"]) and extra["excess_ci95"][0] > 0
    b1 = e2["primary"]
    b1_ok = bool(b1["validated"] and b1["excess_condition_met"])
    b2_ok = bool(mrv["validated"])
    return {
        "A": {"usable": a_ok, "evidence": f"mean AUC {kg['summary']['mean_auc']:.3f} ({kg['summary']['ci95_low']:.3f}-{kg['summary']['ci95_high']:.3f}), {kg['n_evaluated']} diseases; "
                                          f"excess over degree {extra['mean_excess_over_degree']:.3f} (lower bound {extra['excess_ci95'][0]:.3f})"},
        "B1_phenotype_module": {"usable": b1_ok, "evidence": f"without causal gene: mean AUC {b1['summary']['mean_auc']:.3f} ({b1['summary']['ci95_low']:.3f}-{b1['summary']['ci95_high']:.3f}), "
                                                           f"{b1['n_evaluable']} diseases; excess over degree {b1['excess_over_degree']['mean_excess_over_degree']:.3f} "
                                                           f"(lower bound {b1['excess_over_degree']['ci95_low']:.3f}); excess condition met: {b1['excess_condition_met']}"},
        "B2_mitosis_module": {"usable": b2_ok, "evidence": f"module rule: mean AUC {mrv['summary']['mean_auc']:.3f} ({mrv['summary']['ci95_low']:.3f}-{mrv['summary']['ci95_high']:.3f}), {mrv['n_evaluable']} diseases",
                              "check_without_drugs_acting_on_causal_genes": None if not mrv_x else
                              {"n_evaluable": mrv_x["n_evaluable"], "test_is_valid": mrv_x["test_is_valid"], "mean_auc": mrv_x["summary"]["mean_auc"],
                               "excess_lower_bound": mrv_x["excess_over_degree"]["ci95_low"]}},
        "C": {"usable": bool(vu["line_c_used"]), "evidence": f"adjusted total {vu['observed_total']} against null 95th percentile {vu['null']['total_p95']:.1f} ({vu['null']['n_sets']} sets)"},
        "D": {"usable": False, "evidence": f"exploratory: {sg['validation']['n_recovered']} of 6 controls recovered, two by construction; no tier assigned"},
    }


def support_sets(out: Path, universe: dict[str, Any]) -> tuple[dict[str, set[str]], dict[str, Any]]:
    by_name = {d["name"].upper(): d["drug_id"] for d in universe["drugs"]}
    names = {d["drug_id"]: d["name"] for d in universe["drugs"]}
    q, nq = load(out, "kg_query_e5.json"), load(out, "kg_nullquery_e5.json")
    top = distinct_top(q["rank_combined"], TOP_A)
    frac = {by_name[n.upper()]: f for n, f in nq["fraction_of_random_queries_with_drug_in_top50"].items() if n.upper() in by_name}
    a_specific = {d for d in top if d in frac and frac[d] < GENERIC_FRACTION}
    a_generic = {d for d in top if d in frac and frac[d] >= GENERIC_FRACTION}
    a_unmeasured = {d for d in top if d not in frac}
    b1 = {d["drug_id"] for d in load(out, "phenotype_proximity_e5.json")["drugs"] if d["robust"]}
    w3 = load(out, "track1_scenarios.json")["variants"]["two_record_links_no_predicted"]["scheme_robust"]["W3_pathogenic_gene_only"]
    b2 = {by_name[n.upper()] for n in w3 if n.upper() in by_name}
    vu = load(out, "vulnerability_e6.json") if (out / "vulnerability_e6.json").exists() else load(out, "c1000_e6/vulnerability_e6.json")
    sets = {"A": a_specific, "B1": b1, "B2": b2, "C": set(vu["supported_drugs"])}
    detail = {"A_top50_distinct": [names[d] for d in top], "A_generic_not_counted": sorted(names[d] for d in a_generic),
              "A_specificity_not_measured_not_counted": sorted(names[d] for d in a_unmeasured), "B1_robust": sorted(names[d] for d in b1), "B2_W3_two_record": sorted(names[d] for d in b2),
              "random_query_fraction_of_A_specific": {names[d]: round(frac[d], 3) for d in sorted(a_specific, key=lambda x: frac[x])}}
    return sets, detail


def combine(universe: dict[str, Any], status: dict[str, dict[str, Any]], sets: dict[str, set[str]], d_calls: dict[str, Any], b2_counts: bool,
            fda_fn=lookup, literature: set[str] | None = None) -> list[dict[str, Any]]:
    by_id = {d["drug_id"]: d for d in universe["drugs"]}
    b_parts = [k for k in ("B1", "B2") if (status[f"{k}_{'phenotype_module' if k == 'B1' else 'mitosis_module'}"]["usable"] and (k == "B1" or b2_counts))]
    b_set = set().union(*(sets[k] for k in b_parts)) if b_parts else set()
    lines = {"A": (status["A"]["usable"], sets["A"]), "B": (bool(b_parts), b_set), "C": (status["C"]["usable"], sets["C"])}
    usable = {k: s for k, (u, s) in lines.items() if u}
    rows = []
    for drug_id in sorted(set().union(*usable.values()) if usable else set()):
        drug = by_id[drug_id]
        supported = [k for k, s in usable.items() if drug_id in s]
        g = gates(drug, fda_fn(drug["name"]))
        passed = all(g.values())
        with_lit = len(supported) == 1 and drug["name"].upper() in (literature or set())
        candidate = passed and len(supported) >= 2
        rows.append({"drug_id": drug_id, "name": drug["name"], "supported_by": supported, "b_through": [p for p in b_parts if drug_id in sets[p]] if "B" in supported else [],
                     "gates": g, "candidate": candidate, "one_line_plus_literature": bool(passed and with_lit),
                     "direction_note_exploratory": (d_calls.get(drug_id) or {}).get("call", "unclear"), "atc": drug.get("atc") or []})
    return sorted(rows, key=lambda r: (-len(r["supported_by"]), r["name"]))


def overlap(sets_by_line: dict[str, set[str]], common: set[str]) -> dict[str, Any]:
    out = {}
    names = sorted(sets_by_line)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            x, y = sets_by_line[a] & common, sets_by_line[b] & common
            k = len(x & y)
            out[f"{a}_and_{b}"] = {"overlap": k, "size_first": len(x), "size_second": len(y), "universe": len(common),
                                   "p_hypergeometric": float(hypergeom.sf(k - 1, len(common), len(x), len(y))) if k else 1.0}
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    args = parser.parse_args()
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    status = line_status(args.out)
    sets, detail = support_sets(args.out, universe)
    d_calls = load(args.out, "signed_effect.json")["drugs"]
    rationale = load_rationale(Path("config/literature_rationale.json"))
    names = {d["drug_id"]: d["name"] for d in universe["drugs"]}
    readings = {}
    for label, counts in (("b2_counted_as_prespecified", True), ("b2_held_exploratory", False)):
        rows = combine(universe, status, sets, d_calls, counts, literature=rationale)
        readings[label] = {"candidates": [r for r in rows if r["candidate"]], "one_line_only": [r["name"] for r in rows if not r["candidate"] and len(r["supported_by"]) == 1 and all(r["gates"].values())],
                           "one_line_plus_literature": [r["name"] for r in rows if r["one_line_plus_literature"]], "n_candidates": sum(1 for r in rows if r["candidate"]),
                           "not_eligible": sorted(r["name"] for r in rows if not all(r["gates"].values()))}
    scored_a = set(load(args.out, "kg_query_e5.json")["rank_combined"])
    scored_b = set(load(args.out, "phenotype_threshold_scores_e5.json")["min2_terms"]["700"])
    common = scored_a & scored_b
    result = {"addendum": "protocol/addendum_review2.md", "lines": status, "support_set_sizes": {k: len(v) for k, v in sets.items()}, "support_detail": detail,
              "n_usable_lines_with_b2_as_prespecified": sum(1 for k in ("A", "C") if status[k]["usable"]) + int(status["B1_phenotype_module"]["usable"] or status["B2_mitosis_module"]["usable"]),
              "n_usable_lines_with_b2_exploratory": sum(1 for k in ("A", "C") if status[k]["usable"]) + int(status["B1_phenotype_module"]["usable"]),
              "readings": readings,
              "overlap_vs_chance_common_universe": overlap({"A": sets["A"], "B1": sets["B1"], "B2": sets["B2"], "C": sets["C"]}, common),
              "scorable_in_both_A_and_B": len(common),
              "sensitivity_single_record_links_and_original_query": "see out/four_lines.json (addendum 10), kept as sensitivity"}
    result["negative_result"] = all(v["n_candidates"] == 0 for v in readings.values())
    (args.out / "four_lines_review2.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    for label, v in readings.items():
        print(label, "| candidates:", [(r["name"], r["supported_by"]) for r in v["candidates"]], "| one line only:", len(v["one_line_only"]))
    print({k: v["usable"] for k, v in status.items()}, "| negative result:", result["negative_result"])


if __name__ == "__main__":
    main()
