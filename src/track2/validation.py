"""Validation according to `protocol/frozen_protocol.md` (P1 to P4). The parameters and rules are not touched here.

    python -m track2.validation --out out
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr

from track2.modules import build_module, load_mechanism
from track2.proximity import build_null, degree_bins, min_distance_to_set, proximity
from track2.run_pipeline import ITERATIONS, load_universe, run_mechanism
from track2.scoring import PROXIMITY_Z_PRIMARY
from track2.string_net import load_interactome

AUC_MIN = 0.70
N_NEGATIVES = 1000
N_BOOTSTRAP = 1000
MECHANISM = Path("config/mechanisms/bub1b_sac_lof.json")
FANCONI = Path("config/mechanisms/fanconi_anemia_demo.json")


def auc(positive: np.ndarray, negative: np.ndarray) -> float:
    """Mann-Whitney AUC: probability that a positive has a higher score than a negative."""
    combined = np.concatenate([positive, negative])
    ranks = combined.argsort().argsort() + 1.0
    # average of ranks in ties
    for value in np.unique(combined):
        mask = combined == value
        if mask.sum() > 1:
            ranks[mask] = ranks[mask].mean()
    u = ranks[: len(positive)].sum() - len(positive) * (len(positive) + 1) / 2
    return float(u / (len(positive) * len(negative)))


def bootstrap_auc(positive: np.ndarray, negative: np.ndarray, repeats: int = N_BOOTSTRAP, seed: int = 11) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = [auc(rng.choice(positive, len(positive)), rng.choice(negative, len(negative))) for _ in range(repeats)]
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def p1_leave_one_out(interactome, module_nodes: np.ndarray, max_seeds: int | None = None) -> dict[str, Any]:
    """z of each seed against the module without it, versus random genes outside the module (same degree)."""
    rng = np.random.default_rng(101)
    bin_of, bins = degree_bins(interactome.degree)
    z_pos = []
    positives = module_nodes
    if max_seeds and len(module_nodes) > max_seeds:  # addendum_robustness.md: fixed subsample for large modules
        positives = np.sort(np.random.default_rng(101).choice(module_nodes, max_seeds, replace=False))
    for node in positives:
        rest = module_nodes[module_nodes != node]
        null = build_null(interactome.adjacency, interactome.degree, rest, ITERATIONS)
        z_pos.append(proximity(min_distance_to_set(interactome.adjacency, rest), np.array([node]), null)["z"])
    # negatives: nodes outside the module with the same degree distribution as the seeds
    outside = np.setdiff1d(np.arange(len(interactome.proteins)), module_nodes)
    negatives = []
    for node in rng.choice(module_nodes, N_NEGATIVES):
        pool = np.intersect1d(bins[bin_of[node]], outside)
        negatives.append(int(rng.choice(pool)) if pool.size else int(rng.choice(outside)))
    full_null = build_null(interactome.adjacency, interactome.degree, module_nodes, ITERATIONS)
    full_dist = min_distance_to_set(interactome.adjacency, module_nodes)
    z_neg = [proximity(full_dist, np.array([n]), full_null, seed_offset=i)["z"] for i, n in enumerate(negatives)]
    pos, neg = -np.array(z_pos), -np.array(z_neg)
    low, high = bootstrap_auc(pos, neg)
    value = auc(pos, neg)
    return {
        "n_seeds": int(len(module_nodes)), "n_seeds_tested": int(len(positives)), "fraction_seeds_z_le_-2": float((np.array(z_pos) <= PROXIMITY_Z_PRIMARY).mean()),
        "auc": value, "auc_ci95": [low, high], "decision": "SUPPORTED_AT_PATHWAY_LEVEL" if value >= AUC_MIN and low > 0.5 else "NOT_SUPPORTED",
        "median_z_seeds": float(np.median(z_pos)), "median_z_negatives": float(np.median(z_neg)),
    }


def _jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if sa | sb else 0.0


def p3_robustness(mechanism: dict[str, Any], universe: dict[str, Any]) -> dict[str, Any]:
    ref = run_mechanism(mechanism, universe, 700)
    top_ref = [r["name"] for r in ref["rows"][:20]]
    z_ref = {r["drug_id"]: r["proximity"]["z"] for r in ref["rows"]}
    out: dict[str, Any] = {"string_thresholds": {}}
    for score in (400, 900):
        other = run_mechanism(mechanism, universe, score)
        common = [d for d in z_ref if d in {r["drug_id"] for r in other["rows"]}]
        zo = {r["drug_id"]: r["proximity"]["z"] for r in other["rows"]}
        out["string_thresholds"][str(score)] = {
            "jaccard_top20": _jaccard(top_ref, [r["name"] for r in other["rows"][:20]]),
            "spearman_z": float(spearmanr([z_ref[d] for d in common], [zo[d] for d in common]).statistic), "n_common": len(common),
        }
    rng = np.random.default_rng(303)
    jaccards = []
    interactome = load_interactome(700)
    module = build_module(mechanism, interactome)
    from track2.scoring import score_drugs
    from track2.modules import reactome_uniprot
    uni = set()
    for p in mechanism["seeds"].get("reactome_pathways", []):
        uni.update(reactome_uniprot(p))
    for rep in range(20):
        keep = rng.choice(module["nodes"], int(round(0.8 * len(module["nodes"]))), replace=False)
        null = build_null(interactome.adjacency, interactome.degree, keep, ITERATIONS)
        rows = score_drugs(universe["drugs"], interactome, keep, uni, null)
        jaccards.append(_jaccard(top_ref, [r["name"] for r in rows[:20]]))
    out["drop_20pct_seeds"] = {"mean_jaccard_top20": float(np.mean(jaccards)), "min": float(np.min(jaccards)), "max": float(np.max(jaccards))}
    out["verdict"] = "STABLE" if out["drop_20pct_seeds"]["mean_jaccard_top20"] >= 0.5 else "UNSTABLE"
    # P4: degree bias
    degree_by_drug = {}
    drugs = {d["drug_id"]: d for d in universe["drugs"]}
    for r in ref["rows"]:
        nodes = interactome.nodes_for_uniprot(drugs[r["drug_id"]]["human_targets"])
        degree_by_drug[r["drug_id"]] = float(np.mean(interactome.degree[nodes]))
    ids = [r["drug_id"] for r in ref["rows"]]
    out["p4_degree_bias_spearman"] = float(spearmanr([degree_by_drug[i] for i in ids], [z_ref[i] for i in ids]).statistic)
    return out


def p2_fanconi(universe: dict[str, Any]) -> dict[str, Any]:
    result = run_mechanism(load_mechanism(FANCONI), universe, 700)
    rows = {r["name"]: r for r in result["rows"]}
    n = len(result["rows"])
    names = {d["name"]: d for d in universe["drugs"] if d["name"]}
    out = {}
    for control in ("METFORMIN", "OXYMETHOLONE", "CISPLATIN", "MITOMYCIN"):
        if control in rows:
            out[control] = {"scorable": True, "z": rows[control]["proximity"]["z"], "rank": rows[control]["rank"], "of": n,
                            "percentile": 100.0 * (1 - (rows[control]["rank"] - 1) / n)}
        else:
            known = control in names
            out[control] = {"scorable": False, "reason": "no protein target with curated mechanism in ChEMBL" if known else "not in the approved universe"}
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    mechanism = load_mechanism(MECHANISM)
    universe = load_universe()
    interactome = load_interactome(700)
    module = build_module(mechanism, interactome)
    results = {
        "protocol": "protocol/frozen_protocol.md",
        "P1_leave_one_out": p1_leave_one_out(interactome, module["nodes"]),
        "P2_fanconi_controls": p2_fanconi(universe),
        "P3_P4_robustness": p3_robustness(mechanism, universe),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "validation_results.json").write_text(json.dumps(results, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(results, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
