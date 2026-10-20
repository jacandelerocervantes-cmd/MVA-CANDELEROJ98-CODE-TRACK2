"""Validation of the module rule used for the Track 1 genes (protocol/addendum_four_lines.md, "Three further analyses", item 3).

    python -m track2.module_rule_validation --out out [--workers 4]

For the diseases of the validation of line B, the causal genes are the genes annotated to the disease's OMIM entries in `genes_to_phenotype.txt`. Each gene gets a module with the same
rule used for the Track 1 genes (`auto_modules.build_mechanism`), the proximity of all scorable drugs is computed at STRING 700 with the two-record measured links and no predicted links,
the per-gene scores are combined with equal weights (`consensus.combined`), and the ROC AUC of the disease's known approved drugs is computed as in the other validations.
Valid with >= 15 evaluable diseases; validated if the mean AUC >= 0.60 with a lower bound > 0.50. Writes `out/module_rule_validation.json`.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from track2.augmented_proximity import augmented_universe
from track2.consensus import combined
from track2.core_n2 import restrict_links

THRESHOLD, ITERATIONS = 700, 1000
_STATE: dict[str, Any] = {}


def causal_genes(g2p: Path, omim_ids: set[str]) -> dict[str, set[str]]:
    """OMIM id -> gene symbols annotated to it in genes_to_phenotype.txt (columns gene_symbol and disease_id)."""
    lines = g2p.read_text(encoding="utf-8").splitlines()
    header = lines[0].lstrip("#").split("\t")
    ig, idis = header.index("gene_symbol"), header.index("disease_id")
    out: dict[str, set[str]] = {o: set() for o in omim_ids}
    for line in lines[1:]:
        c = line.split("\t")
        if len(c) > max(ig, idis) and c[idis] in out and c[ig] not in ("", "-"):
            out[c[idis]].add(c[ig])
    return out


def combine_equal(z_by_gene: dict[str, dict[str, float]]) -> dict[str, float]:
    """Per-gene z maps -> one pseudo-z per drug (negative of the equal-weight sum of max(0, -z)), so that `-z` is the score used by the AUC code."""
    scores = combined(z_by_gene, {g: 1.0 for g in z_by_gene})
    return {d: -v for d, v in scores.items()}


def _init(universe: dict[str, Any]) -> None:
    from track2.string_net import load_interactome

    _STATE["universe"], _STATE["net"] = universe, load_interactome(THRESHOLD)


def _gene_z(gene: str) -> tuple[str, dict[str, float] | None]:
    from track2.auto_modules import build_mechanism, modellable
    from track2.run_pipeline import run_mechanism

    mech = build_mechanism(gene, 1.0, _STATE["net"], {})
    if not modellable(mech, _STATE["net"]):
        return gene, None
    rows = run_mechanism(mech, _STATE["universe"], THRESHOLD, ITERATIONS)["rows"]
    return gene, {r["drug_id"]: r["proximity"]["z"] for r in rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--sources", type=Path, default=Path("data/sources"))
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    from track2 import ot_client
    from track2.benchmark import drug_mean_degree
    from track2.phenotype_validation import evaluate
    from track2.run_pipeline import load_universe
    from track2.string_net import load_interactome

    val = json.loads((args.out / "phenotype_validation.json").read_text(encoding="utf-8"))
    diseases = {d: m["omim"] for d, m in val["modules"].items() if d in val["primary"]["per_disease"]}
    genes_of = causal_genes(args.sources / "hpo_ontology/genes_to_phenotype.txt", {o for v in diseases.values() for o in v})
    disease_genes = {d: sorted(set().union(*(genes_of[o] for o in omims))) for d, omims in diseases.items()}
    universe = load_universe()
    aug = augmented_universe(universe, restrict_links(json.loads(Path("data/dti_e1_all.json").read_text(encoding="utf-8"))["edges"]), None)
    degrees = drug_mean_degree(universe, load_interactome(THRESHOLD))
    known, _ = ot_client.approved_indications(sorted(degrees))
    all_genes = sorted({g for v in disease_genes.values() for g in v})
    print(len(diseases), "diseases;", len(all_genes), "causal genes to model", flush=True)
    z_gene: dict[str, dict[str, float]] = {}
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init, initargs=(aug,)) as ex:
        for gene, z in ex.map(_gene_z, all_genes):
            if z is not None:
                z_gene[gene] = z
            print(gene, "modelled" if z else "no module of >= 10 proteins", flush=True)
    z_by, usable = {}, {}
    for d, genes in disease_genes.items():
        have = [g for g in genes if g in z_gene]
        if have:
            z_by[(d, "rule")] = combine_equal({g: z_gene[g] for g in have})
            usable[d] = have
    # Circularity check (addendum 11, E1/E2): the known drugs of a disease include those that act on its causal gene, which a module built around that gene recovers almost by construction.
    # The decision stays on the pre-specified evaluation; this one leaves those drugs out of the positives and the negatives.
    net = load_interactome(THRESHOLD)
    by_drug_nodes = {d["drug_id"]: set(net.nodes_for_uniprot(d["human_targets"])) for d in universe["drugs"] if d["human_targets"]}
    exclude = {}
    for d in usable:
        c_nodes = set(net.nodes_for_symbols(disease_genes[d]))
        exclude[d] = {drug for drug, nodes in by_drug_nodes.items() if nodes & c_nodes}
    result = {"n_diseases": len(diseases), "n_with_modellable_causal_gene": len(usable), "causal_genes_modelled": usable,
              "diseases_without_causal_gene_module": sorted(set(diseases) - set(usable)), "evaluation": evaluate(z_by, known, degrees, "rule"),
              "evaluation_without_drugs_acting_on_causal_genes": evaluate(z_by, known, degrees, "rule", exclude),
              "n_drugs_excluded_per_disease": {d: len(v) for d, v in exclude.items()}}
    (args.out / "module_rule_validation.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    e = result["evaluation"]
    s = e["summary"]
    print("evaluable", e["n_evaluable"], "| valid test:", e["test_is_valid"], "| mean AUC", s and round(s["mean_auc"], 3), "| CI low", s and round(s["ci95_low"], 3), "| validated:", e["validated"])
    x = result["evaluation_without_drugs_acting_on_causal_genes"]
    xs = x["summary"]
    print("without drugs acting on the causal genes: evaluable", x["n_evaluable"], "| mean AUC", xs and round(xs["mean_auc"], 3), "| CI low", xs and round(xs["ci95_low"], 3), "| validated:", x["validated"],
          "| excess over degree", x["excess_over_degree"])


if __name__ == "__main__":
    main()
