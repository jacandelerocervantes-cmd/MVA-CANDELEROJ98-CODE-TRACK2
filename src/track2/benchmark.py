"""Multi-disease validation benchmark (pre-registered in `protocol/addendum_benchmark.md`).

    python -m track2.benchmark --out out

Disease selection, module, ground truth, metrics and decision rule are fixed by the addendum; nothing is tuned here.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import rankdata

from track2 import http_cache, ot_client
from track2.benchmark_report import write_reports
from track2.benchmark_stats import DiseaseData, auc_from_scores, drug_bins, percentile_ranks, summarize
from track2.proximity import build_null
from track2.run_pipeline import ITERATIONS, load_universe
from track2.scoring import score_drugs
from track2.string_net import Interactome, load_interactome

THRESHOLDS = (700, 400, 900)
MAIN = 700
MIN_SCORE = 0.5
MIN_MODULE_NODES = 10
MIN_TARGETS = 20
MIN_KNOWN_SELECTION = 5
MIN_KNOWN_SENSITIVITY = 3
N_DISEASES = 30
AUC_MIN = 0.60
CI_LOWER_MIN = 0.50
EXCLUDED_PREFIXES = ("HP_",)


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def scorable_ids(universe: dict[str, Any], interactome: Interactome) -> set[str]:
    return {d["drug_id"] for d in universe["drugs"] if interactome.nodes_for_uniprot(d["human_targets"])}


def select_diseases(universe: dict[str, Any], interactome: Interactome, limit: int = N_DISEASES) -> list[dict[str, Any]]:
    """Deterministic rule of section (a): first `limit` qualifying diseases by ascending id at threshold 700."""
    scorable = scorable_ids(universe, interactome)
    by_disease, names = ot_client.approved_indications(sorted(scorable))
    selected: list[dict[str, Any]] = []
    for disease_id in sorted(by_disease):
        if disease_id.startswith(EXCLUDED_PREFIXES) or len(by_disease[disease_id]) < MIN_KNOWN_SELECTION:
            continue
        symbols = ot_client.associated_symbols(disease_id, MIN_SCORE)
        n_mapped = len(interactome.nodes_for_symbols(symbols))
        if len(symbols) >= MIN_TARGETS and n_mapped >= MIN_MODULE_NODES:
            selected.append({"efo_id": disease_id, "name": names[disease_id], "symbols": symbols})
            log(f"selected {len(selected)}: {disease_id} {names[disease_id]} ({len(symbols)} targets, {n_mapped} mapped)")
            if len(selected) == limit:
                break
    return selected


def known_drugs(universe: dict[str, Any]) -> dict[str, set[str]]:
    by_disease, _ = ot_client.approved_indications([d["drug_id"] for d in universe["drugs"] if d["human_targets"]])
    return by_disease


def disease_z(disease: dict[str, Any], interactome: Interactome, universe: dict[str, Any], threshold: int) -> dict[str, float] | None:
    """z of every scorable drug against the disease module; the null is built once per (disease, threshold)."""
    nodes = np.array(interactome.nodes_for_symbols(disease["symbols"]), dtype=np.int64)
    if len(nodes) < MIN_MODULE_NODES:
        return None
    path = http_cache.CACHE_DIR / f"bench_{disease['efo_id']}_{threshold}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))["z"]
    null = build_null(interactome.adjacency, interactome.degree, nodes, repeats=ITERATIONS)
    rows = score_drugs(universe["drugs"], interactome, nodes, set(), null)
    z = {r["drug_id"]: r["proximity"]["z"] for r in rows}
    path.write_text(json.dumps({"z": z, "module_nodes": len(nodes)}, sort_keys=True), encoding="utf-8", newline="\n")
    return z


def drug_mean_degree(universe: dict[str, Any], interactome: Interactome) -> dict[str, float]:
    out: dict[str, float] = {}
    for drug in universe["drugs"]:
        nodes = interactome.nodes_for_uniprot(drug["human_targets"])
        if nodes:
            out[drug["drug_id"]] = float(interactome.degree[nodes].mean())
    return out


def _metrics(score: np.ndarray, positive: np.ndarray) -> dict[str, Any]:
    pct = percentile_ranks(score, positive)
    return {
        "n_scorable": int(len(score)), "n_known": int(positive.sum()), "auc": auc_from_scores(score, positive),
        "mean_percentile": float(pct.mean()), "median_percentile": float(np.median(pct)),
    }


def evaluate_threshold(
    diseases: list[dict[str, Any]], interactome: Interactome, universe: dict[str, Any], known: dict[str, set[str]], threshold: int,
) -> dict[str, dict[str, Any]]:
    """Per disease: z map plus metrics and permutation inputs; diseases that are not evaluable carry a reason."""
    degrees = drug_mean_degree(universe, interactome)
    out: dict[str, dict[str, Any]] = {}
    for disease in diseases:
        started = time.time()
        z = disease_z(disease, interactome, universe, threshold)
        if z is None:
            out[disease["efo_id"]] = {"evaluable": False, "reason": f"module < {MIN_MODULE_NODES} mapped nodes"}
            continue
        ids = sorted(z)
        positive = np.array([i in known[disease["efo_id"]] for i in ids])
        n_known = int(positive.sum())
        min_known = MIN_KNOWN_SELECTION if threshold == MAIN else MIN_KNOWN_SENSITIVITY
        if n_known < min_known:
            out[disease["efo_id"]] = {"evaluable": False, "reason": f"{n_known} known scorable drugs < {min_known}"}
            continue
        score = -np.array([z[i] for i in ids])
        data = DiseaseData(score, positive, drug_bins(np.array([degrees[i] for i in ids])))
        out[disease["efo_id"]] = {
            "evaluable": True, "module_nodes": len(interactome.nodes_for_symbols(disease["symbols"])), "z": z, "data": data,
            "metrics": _metrics(score, positive),
        }
        log(f"threshold {threshold} {disease['efo_id']} {disease['name']}: AUC {out[disease['efo_id']]['metrics']['auc']:.3f} ({time.time() - started:.0f}s)")
    return out


def consensus(per_threshold: dict[int, dict[str, dict[str, Any]]], diseases: list[dict[str, Any]], known: dict[str, set[str]],
              universe: dict[str, Any], interactome_main: Interactome) -> dict[str, dict[str, Any]]:
    """Mean z-rank across the three thresholds, over drugs scorable at all three and diseases evaluable at all three."""
    degrees = drug_mean_degree(universe, interactome_main)
    out: dict[str, dict[str, Any]] = {}
    for disease in diseases:
        eid = disease["efo_id"]
        parts = [per_threshold[t][eid] for t in THRESHOLDS]
        if not all(p["evaluable"] for p in parts):
            out[eid] = {"evaluable": False, "reason": "not evaluable at all three thresholds"}
            continue
        ids = sorted(set.intersection(*(set(p["z"]) for p in parts)) & set(degrees))
        ranks = np.mean([rankdata([p["z"][i] for i in ids]) for p in parts], axis=0)
        positive = np.array([i in known[eid] for i in ids])
        if positive.sum() < MIN_KNOWN_SENSITIVITY:
            out[eid] = {"evaluable": False, "reason": f"{int(positive.sum())} known drugs scorable at all thresholds"}
            continue
        data = DiseaseData(-ranks, positive, drug_bins(np.array([degrees[i] for i in ids])))
        out[eid] = {"evaluable": True, "data": data, "metrics": _metrics(-ranks, positive)}
    return out


def summarize_threshold(results: dict[str, dict[str, Any]]) -> dict[str, Any] | None:
    ok = [r for r in results.values() if r["evaluable"]]
    if len(ok) < 2:
        return None
    return summarize([r["metrics"]["auc"] for r in ok], [r["data"] for r in ok])


def decide(summary_main: dict[str, Any] | None) -> dict[str, Any]:
    supported = bool(summary_main and summary_main["mean_auc"] >= AUC_MIN and summary_main["ci95_low"] > CI_LOWER_MIN)
    return {
        "rule": f"mean AUC >= {AUC_MIN} and 95% bootstrap CI lower bound > {CI_LOWER_MIN} at STRING threshold {MAIN}",
        "supported_across_diseases": supported,
        "verdict": "supported across diseases" if supported else "not supported across diseases (ranking stays exploratory)",
    }


def run(out_dir: Path) -> dict[str, Any]:
    universe = load_universe()
    interactomes: dict[int, Interactome] = {MAIN: load_interactome(MAIN)}
    log("selecting diseases")
    diseases = select_diseases(universe, interactomes[MAIN])
    known = known_drugs(universe)
    per_threshold: dict[int, dict[str, dict[str, Any]]] = {}
    for threshold in THRESHOLDS:
        interactomes.setdefault(threshold, load_interactome(threshold))
        per_threshold[threshold] = evaluate_threshold(diseases, interactomes[threshold], universe, known, threshold)
        if threshold != MAIN:
            del interactomes[threshold]
    per_threshold_all: dict[Any, dict[str, dict[str, Any]]] = dict(per_threshold)
    per_threshold_all["consensus"] = consensus(per_threshold, diseases, known, universe, interactomes[MAIN])
    summaries = {str(k): summarize_threshold(v) for k, v in per_threshold_all.items()}
    decision = decide(summaries[str(MAIN)])
    log(f"decision: {decision['verdict']}")
    write_reports(out_dir, diseases, per_threshold_all, summaries, decision, universe, ot_client.EXCLUDED_DATASOURCES)
    return decision


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    run(args.out)


if __name__ == "__main__":
    sys.exit(main())
