"""Statistics of the multi-disease benchmark: AUC, percentile ranks, bootstrap over diseases, degree-bin permutation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import rankdata

N_BOOTSTRAP = 10000
BOOTSTRAP_SEED = 7
N_PERMUTATIONS = 1000
PERMUTATION_SEED = 13
N_DRUG_BINS = 10


@dataclass
class DiseaseData:
    """Scores for one disease: higher score = predicted more likely indicated (score = -z or -mean rank)."""

    score: np.ndarray        # (n_drugs,)
    positive: np.ndarray     # (n_drugs,) bool
    drug_bin: np.ndarray     # (n_drugs,) int, degree-matched bin of each drug


def auc_from_scores(score: np.ndarray, positive: np.ndarray) -> float:
    """Mann-Whitney AUC, ties count 1/2."""
    n_pos = int(positive.sum())
    n_neg = len(score) - n_pos
    ranks = rankdata(score)
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def percentile_ranks(score: np.ndarray, positive: np.ndarray) -> np.ndarray:
    """Rank of each positive from the top (best = 1) divided by the number of drugs; lower is better."""
    top_rank = rankdata(-score)
    return top_rank[positive] / len(score)


def drug_bins(mean_degree: np.ndarray, n_bins: int = N_DRUG_BINS) -> np.ndarray:
    """Equal-count bins of the drugs' mean target degree."""
    order = np.argsort(mean_degree, kind="stable")
    bins = np.empty(len(mean_degree), dtype=np.int64)
    for b, members in enumerate(np.array_split(order, min(n_bins, len(order)))):
        bins[members] = b
    return bins


def bootstrap_mean_ci(values: np.ndarray, repeats: int = N_BOOTSTRAP, seed: int = BOOTSTRAP_SEED) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(values), size=(repeats, len(values)))
    means = values[draws].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def permutation_test(items: list[DiseaseData], repeats: int = N_PERMUTATIONS, seed: int = PERMUTATION_SEED) -> dict[str, float]:
    """Shuffle positive labels inside degree-matched drug bins; statistic = mean AUC across diseases."""
    rng = np.random.default_rng(seed)
    observed = float(np.mean([auc_from_scores(d.score, d.positive) for d in items]))
    ranks = [rankdata(d.score) for d in items]
    members = [[np.flatnonzero(d.drug_bin == b) for b in np.unique(d.drug_bin)] for d in items]
    n_pos_in_bin = [[int(d.positive[m].sum()) for m in ms] for d, ms in zip(items, members)]
    null = np.empty(repeats)
    for r in range(repeats):
        aucs = []
        for d, rk, ms, ks in zip(items, ranks, members, n_pos_in_bin):
            rank_sum = 0.0
            for m, k in zip(ms, ks):
                if k:
                    rank_sum += rk[rng.choice(m, size=k, replace=False)].sum()
            n_pos = int(d.positive.sum())
            aucs.append((rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * (len(d.score) - n_pos)))
        null[r] = np.mean(aucs)
    return {
        "observed_mean_auc": observed, "null_mean": float(null.mean()), "null_sd": float(null.std()),
        "p_one_sided": float((1 + np.sum(null >= observed)) / (1 + repeats)), "permutations": repeats,
    }


def summarize(aucs: list[float], items: list[DiseaseData]) -> dict[str, object]:
    """Cross-disease summary: mean AUC, bootstrap CI over diseases, permutation test."""
    values = np.array(aucs)
    low, high = bootstrap_mean_ci(values)
    return {
        "n_diseases": len(aucs), "mean_auc": float(values.mean()), "median_auc": float(np.median(values)),
        "ci95_low": low, "ci95_high": high, "n_auc_above_0_5": int((values > 0.5).sum()),
        "permutation": permutation_test(items),
    }
