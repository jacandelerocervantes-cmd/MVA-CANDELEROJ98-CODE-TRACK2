"""Validation of line B (rule: `protocol/addendum_four_lines.md`, amendment 1).

    python -m track2.phenotype_validation --out out [--workers 6]

For rare diseases that have OMIM cross-references in MONDO, HPO annotations, at least 3 known approved scorable drugs and a module of at least 10
nodes, builds a phenotype module from the disease's own HPO terms with the rule used for the patient, computes the proximity of all scorable drugs at STRING 700
and the ROC AUC of the disease's known approved drugs. Valid only with >= 15 evaluable diseases; line B is validated if mean AUC >= 0.60 and the
lower bound of the bootstrap interval > 0.50 for the patient's rule. Writes `out/phenotype_validation.json`.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

from track2.benchmark_stats import DiseaseData, auc_from_scores, drug_bins, summarize
from track2.hpo import parse_gene_annotations, parse_obo
from track2.phenotype_module import build_phenotype_module

MIN_KNOWN = 3
MIN_HPO_TERMS = 5
MIN_NODES = 10
MAX_DISEASES = 60
MIN_VALID = 15
AUC_MIN, CI_LOWER_MIN = 0.60, 0.50
THRESHOLD = 700
PRIMARY, SECONDARY = "min2_terms", "top100_information_content"
_STATE: dict[str, Any] = {}


def parse_mondo_omim(path: Path) -> dict[str, set[str]]:
    """MONDO id (as Open Targets writes it, MONDO_0000001) -> OMIM ids from the term's own `xref:` lines."""
    out: dict[str, set[str]] = {}
    current, obsolete, omim = None, False, set()

    def close() -> None:
        if current and not obsolete and omim:  # `is_obsolete` can come after the xrefs, so a term is stored when it ends
            out[current] = omim

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("["):
            close()
            current, obsolete, omim = None, False, set()
        elif line.startswith("id: MONDO:"):
            current = line[4:].replace(":", "_")
        elif line == "is_obsolete: true":
            obsolete = True
        elif line.startswith("xref: OMIM:"):
            omim.add(line[6:].split()[0].split("{")[0])
    close()
    return out


def parse_hpoa(path: Path, alt: dict[str, str]) -> dict[str, set[str]]:
    """OMIM id -> HPO terms (phenotypic abnormality aspect, not negated)."""
    out: dict[str, set[str]] = defaultdict(set)
    header: list[str] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if header is None:
            header = cols
            continue
        row = dict(zip(header, cols))
        if row["database_id"].startswith("OMIM:") and row.get("qualifier") != "NOT" and row.get("aspect") == "P":
            out[row["database_id"]].add(alt.get(row["hpo_id"], row["hpo_id"]))
    return dict(out)


def select_diseases(known: dict[str, set[str]], scorable: set[str], mondo: dict[str, set[str]], hpoa: dict[str, set[str]]) -> list[dict[str, Any]]:
    """Ascending MONDO id; terms = union of the HPO terms of the disease's OMIM entries (module size is checked later)."""
    out = []
    for disease in sorted(known):
        omim = sorted(mondo.get(disease, ()))
        terms = sorted(set().union(*(hpoa.get(o, set()) for o in omim))) if omim else []
        if len(terms) >= MIN_HPO_TERMS and len(known[disease] & scorable) >= MIN_KNOWN:
            out.append({"disease": disease, "omim": omim, "terms": terms})
    return out


def _init(universe: dict[str, Any]) -> None:
    from track2.string_net import load_interactome

    _STATE["universe"], _STATE["net"] = universe, load_interactome(THRESHOLD)


def _score(job: tuple[str, str, list[str]]) -> tuple[str, str, dict[str, float] | None]:
    from track2.benchmark import ITERATIONS
    from track2.proximity import build_null
    from track2.scoring import score_drugs

    disease, variant, genes = job
    net = _STATE["net"]
    nodes = np.array(net.nodes_for_symbols(genes), dtype=np.int64)
    if len(nodes) < MIN_NODES:
        return disease, variant, None
    null = build_null(net.adjacency, net.degree, nodes, repeats=ITERATIONS)
    rows = score_drugs(_STATE["universe"]["drugs"], net, nodes, set(), null)
    return disease, variant, {r["drug_id"]: r["proximity"]["z"] for r in rows}


def evaluate(z_by: dict[tuple[str, str], dict[str, float] | None], known: dict[str, set[str]], degrees: dict[str, float], variant: str,
             exclude: dict[str, set[str]] | None = None) -> dict[str, Any]:
    """AUC per disease and summary. `exclude[disease]` are drugs left out of both the positives and the negatives (addendum 11, E2: drugs that act on a causal gene).
    Also reports the excess of each AUC over the degree-only baseline (score = mean target degree of the drug) with a bootstrap interval over diseases (E3)."""
    from track2.benchmark_stats import bootstrap_mean_ci

    aucs, items, per, excess = [], [], {}, []
    for (disease, v), z in sorted(z_by.items()):
        if v != variant or z is None:
            continue
        ids = sorted((set(z) & set(degrees)) - (exclude or {}).get(disease, set()))
        positive = np.array([i in known[disease] for i in ids])
        if positive.sum() < MIN_KNOWN:
            continue
        score = -np.array([z[i] for i in ids])
        deg = np.array([degrees[i] for i in ids])
        auc = auc_from_scores(score, positive)
        aucs.append(auc)
        excess.append(auc - auc_from_scores(deg, positive))
        items.append(DiseaseData(score, positive, drug_bins(deg)))
        per[disease] = {"auc": round(auc, 4), "auc_degree_only": round(auc - excess[-1], 4), "n_known": int(positive.sum()), "n_scorable": len(ids)}
    valid = len(aucs) >= MIN_VALID
    summary = summarize(aucs, items) if len(aucs) >= 2 else None
    ex = None
    if len(excess) >= 2:
        low, high = bootstrap_mean_ci(np.array(excess))
        ex = {"mean_excess_over_degree": float(np.mean(excess)), "ci95_low": low, "ci95_high": high}
    return {"n_evaluable": len(aucs), "test_is_valid": valid, "summary": summary, "per_disease": per, "excess_over_degree": ex,
            "validated": bool(valid and summary and summary["mean_auc"] >= AUC_MIN and summary["ci95_low"] > CI_LOWER_MIN),
            "excess_condition_met": bool(ex and ex["ci95_low"] > 0)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--sources", type=Path, default=Path("data/sources"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--leave-causal", action="store_true",
                        help="addendum 11, E2: remove the causal genes of each disease from its module and the drugs that act on them from the ground truth; writes phenotype_validation_e2.json")
    args = parser.parse_args()
    from track2 import ot_client
    from track2.benchmark import drug_mean_degree
    from track2.run_pipeline import load_universe
    from track2.string_net import load_interactome

    parents, alt = parse_obo(args.sources / "hpo_ontology/hp.obo")
    annotations = parse_gene_annotations(args.sources / "hpo_ontology/genes_to_phenotype.txt", alt)
    mondo = parse_mondo_omim(args.sources / "mondo/mondo.obo")
    hpoa = parse_hpoa(args.sources / "hpo_annotations/phenotype.hpoa", alt)
    universe = load_universe()
    net = load_interactome(THRESHOLD)
    degrees = drug_mean_degree(universe, net)
    scorable = set(degrees)
    known, _ = ot_client.approved_indications(sorted(scorable))
    candidates = select_diseases(known, scorable, mondo, hpoa)
    print("candidate diseases (MONDO+OMIM, >=5 HPO terms, >=3 known drugs):", len(candidates), flush=True)
    causal: dict[str, set[str]] = {}
    exclude: dict[str, set[str]] = {}
    if args.leave_causal:
        from track2.module_rule_validation import causal_genes

        genes_of = causal_genes(args.sources / "hpo_ontology/genes_to_phenotype.txt", {o for d in candidates for o in d["omim"]})
        by_drug_nodes = {d["drug_id"]: set(net.nodes_for_uniprot(d["human_targets"])) for d in universe["drugs"] if d["human_targets"]}
        for d in candidates:
            causal[d["disease"]] = set().union(*(genes_of.get(o, set()) for o in d["omim"]))
            c_nodes = set(net.nodes_for_symbols(sorted(causal[d["disease"]])))
            exclude[d["disease"]] = {drug for drug, nodes in by_drug_nodes.items() if nodes & c_nodes}
    jobs, modules = [], {}
    for d in candidates:
        m = build_phenotype_module(d["terms"], annotations, parents, alt)
        if args.leave_causal:  # the causal genes leave every definition of the module
            m["variants"] = {k: sorted(set(v) - causal[d["disease"]]) for k, v in m["variants"].items()}
        modules[d["disease"]] = {"omim": d["omim"], "n_terms": len(d["terms"]), "rule": m["rule"], "sizes": {k: len(v) for k, v in m["variants"].items()},
                                 **({"causal_genes": sorted(causal[d["disease"]]), "n_drugs_excluded": len(exclude[d["disease"]])} if args.leave_causal else {})}
        for variant in (PRIMARY, SECONDARY):
            if len(net.nodes_for_symbols(m["variants"][variant])) >= MIN_NODES:
                jobs.append((d["disease"], variant, m["variants"][variant]))
    eligible = sorted({j[0] for j in jobs if j[1] == PRIMARY})[:MAX_DISEASES]
    jobs = [j for j in jobs if j[0] in eligible]
    print("diseases with a module of >= 10 nodes (primary rule):", len(eligible), "| jobs:", len(jobs), flush=True)
    z_by: dict[tuple[str, str], dict[str, float] | None] = {}
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init, initargs=(universe,)) as ex:
        for disease, variant, z in ex.map(_score, jobs):
            z_by[(disease, variant)] = z
            print(disease, variant, "ok" if z else "module too small", flush=True)
    result = {"rule": {"min_known": MIN_KNOWN, "min_hpo_terms": MIN_HPO_TERMS, "min_valid_diseases": MIN_VALID, "auc_min": AUC_MIN, "ci_lower_min": CI_LOWER_MIN,
                       "string_score": THRESHOLD, "primary": PRIMARY},
              "n_candidates": len(candidates), "modules": {d: modules[d] for d in eligible},
              "primary": evaluate(z_by, known, degrees, PRIMARY, exclude), "secondary_top100": evaluate(z_by, known, degrees, SECONDARY, exclude)}
    p = result["primary"]
    result["e2_leave_causal_gene"] = bool(args.leave_causal)
    usable = p["validated"] and (p["excess_condition_met"] if args.leave_causal else True)
    result["verdict"] = ("line B validated" if usable else "line B not validated" if p["test_is_valid"]
                         else f"not testable: {p['n_evaluable']} evaluable diseases (< {MIN_VALID}); line B stays exploratory")
    (args.out / ("phenotype_validation_e2.json" if args.leave_causal else "phenotype_validation.json")).write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(result["verdict"], "| mean AUC:", p["summary"] and round(p["summary"]["mean_auc"], 3), "| CI low:", p["summary"] and round(p["summary"]["ci95_low"], 3))


if __name__ == "__main__":
    main()
