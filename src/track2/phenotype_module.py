"""Phenotype module of line B (`protocol/addendum_four_lines.md`): genes annotated to the patient's HPO terms.

Rule, fixed before any data were read:
1. A gene is annotated to a patient term if it carries the term or any descendant of it (annotations are propagated to ancestors).
2. The module is the set of genes annotated to at least `MIN_TERMS` (2) of the patient terms.
3. If that gives fewer than `MIN_GENES` (10) genes, the module is the `TOP_N` (100) genes with the highest score, where a gene's
   score is the sum of the information content of the patient terms it is annotated to (genes with a zero score are never added).
Terms with no annotation in the release are reported, never silently dropped.
"""

from __future__ import annotations

from typing import Any

from track2.hpo import information_content, propagate

MIN_TERMS = 2
MIN_GENES = 10
TOP_N = 100
# Addendum 11, E5: terms that describe the parents and not the child are not part of any phenotype query of lines A and B.
NOT_THE_CHILD = {"HP:0200067"}  # recurrent spontaneous abortion (parental history)
CORE_TERMS = {"HP:0002859", "HP:0000121", "HP:0003202"}  # rhabdomyosarcoma, nephrocalcinosis, skeletal muscle atrophy
NONSPECIFIC_TERMS = {"HP:0004322", "HP:0001508", "HP:0001518", "HP:0001622"}  # short stature, failure to thrive, small for gestational age, premature birth


def query_terms(patient: dict[str, Any], leave_out: set[str] | None = None, only: set[str] | None = None) -> list[str]:
    """HPO terms of the child used in the phenotype queries: all of the patient's terms except those of `NOT_THE_CHILD` and `leave_out`; restricted to `only` if given."""
    terms = [t["id"] for t in patient["hpo_terms"] if t["id"] not in NOT_THE_CHILD and t["id"] not in (leave_out or set())]
    return [t for t in terms if only is None or t in only]


def build_phenotype_module(terms: list[str], annotations: dict[str, set[str]], parents: dict[str, set[str]],
                           alt: dict[str, str] | None = None, weights: dict[str, float] | None = None) -> dict[str, Any]:
    """`weights` (term -> weight, default 1) are used by the age scenarios; the count rule uses the terms only."""
    alt = alt or {}
    terms = sorted({alt.get(t, t) for t in terms})
    propagated = propagate(annotations, parents)
    ic = information_content(propagated)
    w = {t: (weights or {}).get(t, 1.0) for t in terms}
    unannotated = [t for t in terms if t not in ic]
    per_gene: dict[str, dict[str, Any]] = {}
    for gene, have in propagated.items():
        hit = [t for t in terms if t in have]
        if hit:
            per_gene[gene] = {"terms": hit, "n_terms": len(hit), "score": sum(ic[t] * w[t] for t in hit)}
    chosen = sorted(g for g, v in per_gene.items() if v["n_terms"] >= MIN_TERMS)
    ranked = sorted(per_gene, key=lambda g: (-per_gene[g]["score"], g))
    top_n = sorted(ranked[:TOP_N])
    rule = f"annotated_to_at_least_{MIN_TERMS}_terms"
    if len(chosen) < MIN_GENES:
        chosen, rule = top_n, f"top_{TOP_N}_by_information_content"
    # `variants` are the alternative definitions of the module used for the nine configurations of line B (phenotype_proximity.py)
    variants = {"min2_terms": sorted(g for g, v in per_gene.items() if v["n_terms"] >= 2),
                "min3_terms": sorted(g for g, v in per_gene.items() if v["n_terms"] >= 3),
                f"top{TOP_N}_information_content": top_n}
    # scores are rounded only in the written detail (ranking above uses the unrounded values): the last digits of math.log differ between operating systems (about 1e-15), which
    # would otherwise make the file differ byte for byte between Windows and Linux
    detail = {g: {**per_gene[g], "score": round(per_gene[g]["score"], 10)} for g in sorted(set(chosen) | set(top_n))}
    return {"terms": terms, "terms_without_annotation": unannotated, "rule": rule, "genes": chosen, "variants": variants,
            "n_genes_annotated_to_any_term": len(per_gene), "gene_detail": detail}


def main() -> None:
    import argparse
    import json
    from pathlib import Path

    from track2.hpo import parse_gene_annotations, parse_obo

    parser = argparse.ArgumentParser(description="Build the phenotype module of line B from the patient's HPO terms (addendum 11, E5: without the parental term).")
    parser.add_argument("--patient", type=Path, default=Path("config/patient.json"))
    parser.add_argument("--sources", type=Path, default=Path("data/sources/hpo_ontology"))
    parser.add_argument("--out", type=Path, default=Path("out/phenotype_module_e5.json"))
    parser.add_argument("--leave-out", nargs="*", default=[], help="HPO ids left out (ablation)")
    parser.add_argument("--only", choices=["core", "nonspecific"], default=None, help="restrict to the core MVA terms or to the non-specific constitutional terms")
    args = parser.parse_args()
    patient = json.loads(args.patient.read_text(encoding="utf-8"))
    only = {"core": CORE_TERMS, "nonspecific": NONSPECIFIC_TERMS}.get(args.only)
    terms = query_terms(patient, set(args.leave_out), only)
    parents, alt = parse_obo(args.sources / "hp.obo")
    module = build_phenotype_module(terms, parse_gene_annotations(args.sources / "genes_to_phenotype.txt", alt), parents, alt)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(module, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("terms:", module["terms"], "| rule:", module["rule"], "| genes:", len(module["genes"]), "| variants:", {k: len(v) for k, v in module["variants"].items()})


if __name__ == "__main__":
    main()
