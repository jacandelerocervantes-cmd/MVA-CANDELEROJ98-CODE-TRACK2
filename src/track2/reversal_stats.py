"""Aggregation of SigCom LINCS two-sided results per drug, empirical null and FDR (pure functions, no network)."""

from __future__ import annotations

from typing import Any

import numpy as np


def per_drug(results: list[dict[str, Any]], drug_of: dict[str, list[str]]) -> dict[str, dict[str, Any]]:
    """Best reversal / mimic per drug. z-sum < 0 means the signature moves against the query (reversal).

    `drug_of` maps a signature UUID to the approved drug ids it belongs to. The API FDR of a signature is
    conservatively taken as the larger of its up-gene and down-gene values.
    """
    out: dict[str, dict[str, Any]] = {}
    for row in results:
        for drug_id in drug_of.get(row["uuid"], []):
            rec = out.setdefault(drug_id, {"rev": 0.0, "mim": 0.0, "api_fdr_rev": 1.0})
            z = float(row["z-sum"])
            if z < 0 and -z > rec["rev"]:
                rec["rev"] = -z
                rec["api_fdr_rev"] = max(float(row["fdr-up"]), float(row["fdr-down"]))
            elif z > rec["mim"]:
                rec["mim"] = z
    return out


def empirical_p(observed: np.ndarray, null: np.ndarray) -> np.ndarray:
    """P(null >= observed) with the +1 correction; `null` is the pooled drug-level null."""
    ordered = np.sort(null)
    ge = len(ordered) - np.searchsorted(ordered, observed, side="left")
    return (1.0 + ge) / (1.0 + len(ordered))


def bh_fdr(p: np.ndarray) -> np.ndarray:
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def direction(fdr_rev: float, fdr_mim: float, alpha: float = 0.05) -> str:
    if fdr_rev <= alpha and fdr_rev <= fdr_mim:
        return "reverses"
    if fdr_mim <= alpha:
        return "mimics"
    return "none"


def null_summary(null_rev: np.ndarray, per_perm_false_hits: list[int], n_drugs: int) -> dict[str, Any]:
    q = np.quantile(null_rev, [0.5, 0.95, 0.99, 0.999])
    return {
        "n_null_values": int(len(null_rev)), "n_permutations": len(per_perm_false_hits), "n_drugs_scored": n_drugs,
        "null_reversal_score_quantiles": {"p50": round(float(q[0]), 4), "p95": round(float(q[1]), 4),
                                          "p99": round(float(q[2]), 4), "p99.9": round(float(q[3]), 4)},
        "null_max_reversal_score": round(float(null_rev.max()), 4),
        "false_reverses_calls_per_random_set_leave_one_out": per_perm_false_hits,
        "mean_false_calls_per_random_set": round(float(np.mean(per_perm_false_hits)), 3) if per_perm_false_hits else None,
    }
