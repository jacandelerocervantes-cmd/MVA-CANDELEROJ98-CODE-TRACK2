"""Forking-path ledger (addendum 11, E12): every variant of the analysis that was run, with the one declared as primary for each line.

    python -m track2.forking_ledger --out out

The entries are written by hand from `protocol/` (addenda 10 and 11, `deviations.md`); the figures are read from the result files so that they cannot drift. Everything not marked
`primary` is a sensitivity analysis and is never presented as a finding. Writes `out/forking_ledger.json`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _load(out: Path, name: str) -> Any:
    return json.loads((out / name).read_text(encoding="utf-8"))


def ledger(out: Path) -> list[dict[str, Any]]:
    kg, mrv, e2, orig = _load(out, "kg_validation.json"), _load(out, "module_rule_validation.json"), _load(out, "phenotype_validation_e2.json"), _load(out, "phenotype_validation.json")
    c200, c1000 = _load(out, "vulnerability.json"), _load(out, "vulnerability_e6.json")
    fl_old, fl_new = _load(out, "four_lines.json"), _load(out, "four_lines_review2.json")
    return [
        {"line": "A", "variant": "PrimeKG personalised PageRank, graph with drug_effect edges (first run)", "status": "superseded", "when": "addendum 10",
         "result": "mean AUC 0.821 (0.801-0.839); ranking dominated by drugs that cause the phenotypes", "why_changed": "post hoc: side-effect edges removed after the ranking was read (amendment 2)"},
        {"line": "A", "variant": "graph without drug-disease and drug_effect edges; patient query with 8 HPO terms (7 with a node)", "status": "sensitivity", "when": "addendum 10",
         "result": f"validation mean AUC {kg['summary']['mean_auc']:.3f}; patient ranking in out/kg_query.json"},
        {"line": "A", "variant": "same graph; patient query without the parental term HP:0200067 (6 phenotype nodes), specificity against 500 degree-matched random queries", "status": "primary",
         "when": "addendum 11, E4/E5", "result": "43 of 46 drugs of the real top 50 are in the top 50 of >= 5 % of random queries; fostamatinib is in the top 5 of 98.8 % of them"},
        {"line": "A", "variant": "term ablation (each of 6 terms left out), core-only and non-specific-only queries", "status": "sensitivity", "when": "addendum 11, E5",
         "result": "Jaccard of the top 50 against the primary query between 0.89 and 1.00"},
        {"line": "B1", "variant": "phenotype module, curated + measured links with single-record links, validation on 19 diseases", "status": "superseded", "when": "addendum 10, amendment 1",
         "result": f"mean AUC {orig['primary']['summary']['mean_auc']:.3f}", "why_changed": "validation replaced before it ran (the 30-disease benchmark could not be run); thresholds unchanged"},
        {"line": "B1", "variant": "validation without the causal genes of each disease and without the drugs that act on them; excess over the degree baseline required", "status": "primary",
         "when": "addendum 11, E2/E3", "result": f"mean AUC {e2['primary']['summary']['mean_auc']:.3f}, excess lower bound {e2['primary']['excess_over_degree']['ci95_low']:.3f}: not usable"},
        {"line": "B1", "variant": "module defined by >= 3 terms, and by the 100 best terms by information content (sensitivity of the nine configurations)", "status": "sensitivity", "when": "addendum 10",
         "result": f"top-100 variant mean AUC {e2['secondary_top100']['summary']['mean_auc']:.3f} in the primary validation, excess condition met: {e2['secondary_top100']['excess_condition_met']}"},
        {"line": "B2", "variant": "mitosis modules W1/W2/W3, nine configurations, single-record links", "status": "superseded", "when": "addendum 10",
         "result": "fostamatinib robust 9/9 (single-record links)", "why_changed": "post hoc: two-record rule applied after fostamatinib's support was read (E8)"},
        {"line": "B2", "variant": "BUB1B-only module (W3), links with >= 2 records, no predicted links; module rule validated on diseases with a known causal gene", "status": "primary", "when": "addendum 11, E1/E7/E8",
         "result": f"module rule validated: mean AUC {mrv['evaluation']['summary']['mean_auc']:.3f}, {mrv['evaluation']['n_evaluable']} diseases; without the drugs acting on the causal genes "
                   f"{mrv['evaluation_without_drugs_acting_on_causal_genes']['n_evaluable']} evaluable diseases (not testable)"},
        {"line": "B2", "variant": "W1 (Track 1 weights) and W2 (equal weights)", "status": "sensitivity", "when": "addendum 11, E7", "result": "pathway-context scenarios; see out/track1_scenarios.json"},
        {"line": "C", "variant": "lineage + proliferation covariates, 200 null sets", "status": "superseded", "when": "addendum 10",
         "result": f"{c200['observed_total']} hits against p95 {c200['null']['total_p95']:.2f}: not usable", "why_changed": "TP53 confounding (6 of 19 PRISM hits were MDM2 inhibitors) and a null that was too small"},
        {"line": "C", "variant": "lineage + proliferation + TP53 + MSI + CIN, 1,000 null sets", "status": "primary", "when": "addendum 11, E6",
         "result": f"{c1000['observed_total']} hits against p95 {c1000['null']['total_p95']:.2f}: not usable; no drug supported"},
        {"line": "C", "variant": "SAC score orthogonalised to proliferation; PRISM and CRISPR totals separate; collapsed by MOA class", "status": "sensitivity", "when": "addendum 11, E6",
         "result": f"orthogonalised hits {c1000['covariates']['orthogonalised_score_hits']}; VIF {c1000['covariates']['vif_of_score']:.2f}"},
        {"line": "D", "variant": "OmniPath signed paths to checkpoint-strength nodes", "status": "exploratory", "when": "addendum 10", "result": "4 of 6 controls recovered, two by construction; no tier assigned"},
        {"line": "combination", "variant": "at least two usable lines among A, B, C; B1 and B2 count as one line", "status": "primary", "when": "addendum 11, E1/E9",
         "result": f"candidates: {fl_new['readings']['b2_counted_as_prespecified']['n_candidates']} (B2 as pre-specified), {fl_new['readings']['b2_held_exploratory']['n_candidates']} (B2 exploratory)"},
        {"line": "combination", "variant": "addendum 10 rule with single-record links and the 8-term query", "status": "superseded", "when": "addendum 10",
         "result": f"candidates: {[c['name'] for c in fl_old['candidates']]}", "why_changed": "its B support disappears with two-record links and its A support is generic"},
        {"line": "combination", "variant": "amendment 3: one usable line plus a literature rationale", "status": "sensitivity", "when": "addendum 10, amendment 3 (post hoc)", "result": "nobody accepted; kept as a labelled list only"},
        {"line": "class", "variant": "ChEMBL mechanism-of-action text and ATC level 4; 10,000 degree-matched random classes; FDR < 0.10", "status": "primary", "when": "addendum 11, E10",
         "result": "see out/class_analysis.json"},
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    entries = ledger(args.out)
    result = {"note": "Every variant of the computational lines that was run. The literature searches (two rounds, 419 and 256 queries) and the selection and editorial rules of E16 are not entries of this ledger; they are described in protocol/addendum_e16.md. Only entries marked primary are findings; the rest are sensitivity analyses or superseded.", "entries": entries,
              "n_entries": len(entries), "n_primary": sum(1 for e in entries if e["status"] == "primary")}
    (args.out / "forking_ledger.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(len(entries), "entries;", result["n_primary"], "primary")


if __name__ == "__main__":
    main()
