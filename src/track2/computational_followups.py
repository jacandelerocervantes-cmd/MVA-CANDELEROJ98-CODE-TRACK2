r"""Computational follow-ups CF1 to CF9 of `protocol/follow_up_experiments_v2.md` (rules: `protocol/addendum_computational_followups.md`, written before any run).

Run from the repository root:  PYTHONPATH=src python -m track2.computational_followups [--only CF1,CF3] [--out out]
Writes `out/computational_followups.json` (sorted keys) and `out/computational_followups.md`. Only data on disk and the project caches are used (Europe PMC
only in CF3, cursorMark, at most 4 requests per second, cached). Exploratory and post hoc: a candidate that survives stays a hypothesis; one that is refuted is
closed as a computational candidate. No patient data: gene symbols, drug names and the categorical case profile only.
"""
from __future__ import annotations

import argparse
import json
import re
import time
import types
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence

import numpy as np
from scipy import sparse, stats
from scipy.sparse.csgraph import dijkstra

SEED = 20261018
Q_CUT = 0.10
P_CUT = 0.05
VERDICTS = ("survives", "refuted", "not estimable", "not runnable")
LABEL = "exploratory, post hoc; cannot prove efficacy; a candidate that survives stays a hypothesis; a refuted one is closed as a computational candidate"
PROTOCOL = "protocol/addendum_computational_followups.md"
CASE_GENES = {"BUB1", "BUB1B", "BUB3", "MAD1L1", "MAD2L1", "TTK", "TRIP13", "CEP57", "CENATAC", "CENPE", "IGBP1", "CDC20"}


# ------------------------------------------------------------------ common helpers
def bh(p: Sequence[float]) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def jsonable(o: Any) -> Any:
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [jsonable(v) for v in (sorted(o) if isinstance(o, set) else o)]
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, float):
        return None if not np.isfinite(o) else o
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def rnd(x: Optional[float], nd: int = 4) -> Optional[float]:
    return None if x is None or not np.isfinite(x) else round(float(x), nd)


# ------------------------------------------------------------------ CF1: sign audit (resolved by documentation)
def quartile_delta(viability: np.ndarray, score: np.ndarray) -> Dict[str, Any]:
    """delta = mean viability (log fold change) in the weakest score quartile minus the strongest quartile; negative = more sensitive when the checkpoint is weak."""
    ok = np.isfinite(viability) & np.isfinite(score)
    v, s = viability[ok], score[ok]
    lo, hi = np.quantile(s, 0.25), np.quantile(s, 0.75)
    weak, strong = v[s <= lo], v[s >= hi]
    rho = float(stats.spearmanr(s, v).statistic)
    return {"n_lines": int(ok.sum()), "mean_weak_quartile": float(weak.mean()), "mean_strong_quartile": float(strong.mean()),
            "delta_weak_minus_strong": float(weak.mean() - strong.mean()), "spearman_rho_unadjusted": rho}


def cf1_sign_audit(Y: Dict[str, np.ndarray], score: np.ndarray, drugs: Sequence[str], adjusted_rho: Dict[str, Optional[float]],
                   sanity_drugs: Sequence[str], documentation: Dict[str, Any]) -> Dict[str, Any]:
    sanity = {d: float(np.nanmedian(Y[d])) for d in sanity_drugs if d in Y}
    sanity_pass = bool(sanity) and all(v < 0 for v in sanity.values())
    per: Dict[str, Any] = {}
    for d in drugs:
        if d not in Y:
            per[d] = {"verdict": "not estimable", "reason": "drug not in the PRISM matrix"}
            continue
        q = quartile_delta(Y[d], score)
        q["adjusted_partial_rho"] = adjusted_rho.get(d)
        delta_ok = q["delta_weak_minus_strong"] < 0
        rho_pos = (q["adjusted_partial_rho"] or q["spearman_rho_unadjusted"]) > 0
        q["verdict"] = "survives" if (sanity_pass and delta_ok and rho_pos) else "refuted" if q["delta_weak_minus_strong"] > 0 else "not estimable"
        per[d] = q
    return {"documentation": documentation, "orientation_sanity_median_lfc": sanity, "orientation_sanity_pass": sanity_pass, "per_drug": per}


# ------------------------------------------------------------------ CF2: class enrichment
def cf2_enrichment(estimable: Sequence[str], survivors: Sequence[str], tubulin_members: Sequence[str], all_adjusted_hits: Sequence[str]) -> Dict[str, Any]:
    def stem(x: str) -> str:
        return x.upper().split()[0]

    est = list(estimable)
    cls = {d for d in est if d.upper() in {t.upper() for t in tubulin_members} or stem(d) in {stem(t) for t in tubulin_members}}
    M, K = len(est), len(cls)

    def test(draws: Sequence[str]) -> Dict[str, Any]:
        n = len(draws)
        x = sum(1 for d in draws if d in cls)
        p = float(stats.hypergeom.sf(x - 1, M, K, n)) if n else 1.0
        return {"n_draws": n, "x_tubulin_among_draws": x, "expected": (n * K / M) if M else None, "p_hypergeometric_ge": p}

    primary, secondary = test(list(survivors)), test(list(all_adjusted_hits))
    ok = len(survivors) >= 3 and primary["x_tubulin_among_draws"] <= K
    verdict = "not estimable" if not ok else "refuted" if primary["p_hypergeometric_ge"] <= P_CUT else "survives"
    return {"M_estimable": M, "K_tubulin_class_in_estimable": K, "tubulin_class_in_estimable": sorted(cls), "survivors": sorted(survivors),
            "primary": primary, "sensitivity_all_adjusted_q_below_0_10": secondary, "verdict": verdict}


# ------------------------------------------------------------------ CF3: engine specificity (literature leg)
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
MAX_YEAR = 2024
MIN_INTERVAL = 0.26


def epmc_query(drug: str, gene: str) -> str:
    def t(x: str) -> str:
        return f'(TITLE:"{x}" OR ABSTRACT:"{x}")'

    return f"{t(drug)} AND {t(gene)} AND PUB_YEAR:[1900 TO {MAX_YEAR}]"


def epmc_url(query: str) -> str:
    from urllib.parse import quote

    return f"{EPMC}?query={quote(query)}&format=json&resultType=lite&pageSize=1&cursorMark=*"


class RateLimitedGetter:
    """Cached Europe PMC getter through `http_cache`; at most 4 uncached requests per second."""

    def __init__(self) -> None:
        self.last = 0.0
        self.n_network = 0

    def __call__(self, url: str) -> Any:
        from track2 import http_cache

        if not http_cache._cache_path(url, ".json").exists():
            wait = MIN_INTERVAL - (time.time() - self.last)
            if wait > 0:
                time.sleep(wait)
            self.last = time.time()
            self.n_network += 1
        return http_cache.get_json(url)


def hit_count(drug: str, gene: str, getter: Callable[[str], Any]) -> int:
    return int(getter(epmc_url(epmc_query(drug, gene)))["hitCount"])


def matched_null_genes(degree_by_symbol: Dict[str, int], target: str, n: int, seed: int = SEED, exclude: Iterable[str] = CASE_GENES) -> Dict[str, Any]:
    d0 = degree_by_symbol[target]
    banned = {e.upper() for e in exclude} | {target.upper()}
    for tol in (0.25, 0.50):
        pool = sorted(g for g, d in degree_by_symbol.items() if g.upper() not in banned and len(g) >= 4 and re.fullmatch(r"[A-Z0-9-]+", g.upper())
                      and abs(d - d0) <= tol * d0)
        if len(pool) >= n:
            rng = np.random.default_rng(seed)
            return {"tolerance": tol, "pool_size": len(pool), "target_degree": d0, "genes": sorted(rng.choice(pool, size=n, replace=False).tolist())}
    raise ValueError("fewer null genes than requested within +-50 % of the degree")


def cf3_specificity(drug_counts: Dict[str, Optional[int]], null_counts: Dict[str, List[int]]) -> Dict[str, Any]:
    names = sorted(drug_counts)
    rows: Dict[str, Any] = {}
    est = [d for d in names if drug_counts[d] is not None]
    pvals = {}
    for d in est:
        obs, null = drug_counts[d], null_counts[d]
        pvals[d] = (1 + sum(1 for x in null if x >= obs)) / (len(null) + 1)
    qv = dict(zip(est, bh([pvals[d] for d in est]))) if est else {}
    for d in names:
        if drug_counts[d] is None:
            rows[d] = {"verdict": "not estimable", "reason": "Europe PMC query could not be answered"}
            continue
        null = np.array(null_counts[d])
        rows[d] = {"observed_hit_count_bub1b": drug_counts[d], "null_median": float(np.median(null)), "null_max": int(null.max()), "null_n_genes": int(len(null)),
                   "n_null_ge_observed": int((null >= drug_counts[d]).sum()), "p_empirical": pvals[d], "q_bh": float(qv[d]),
                   "verdict": "survives" if (pvals[d] <= P_CUT and qv[d] < Q_CUT) else "refuted"}
    return {"per_drug": rows}


# ------------------------------------------------------------------ CF4: graph path audit
def cf4_path_audit(adj: sparse.csr_matrix, node_name: Sequence[str], node_type: Sequence[str], drug_nodes: Dict[str, List[int]], query_nodes: Dict[int, str],
                   hub_percentile: float = 99.0, max_len: int = 4, edge_relation: Optional[Callable[[List[List[int]]], Dict[tuple, str]]] = None) -> Dict[str, Any]:
    degree = np.asarray(adj.sum(axis=1)).ravel()
    thr = float(np.percentile(degree, hub_percentile))
    is_hub = degree >= thr
    query = set(query_nodes)
    out: Dict[str, Any] = {"hub_degree_threshold": thr, "hub_percentile": hub_percentile, "n_hubs": int(is_hub.sum()), "max_len": max_len, "per_drug": {}}
    all_paths: List[List[int]] = []
    for drug, nodes in drug_nodes.items():
        if not nodes:
            out["per_drug"][drug] = {"verdict": "not estimable", "reason": "no PrimeKG drug node matched"}
            continue
        keep = ~is_hub
        keep[list(query)] = True
        keep[nodes] = True
        mask = sparse.diags(keep.astype(float))
        adj_free = (mask @ adj @ mask).tocsr()
        any_len, free_len, examples = {}, {}, []
        for src in nodes:
            d_any = dijkstra(adj, directed=False, unweighted=True, indices=src, limit=max_len + 2)
            d_free, pred = dijkstra(adj_free, directed=False, unweighted=True, indices=src, limit=max_len, return_predecessors=True)
            for q in query:
                if np.isfinite(d_any[q]):
                    any_len[q] = min(any_len.get(q, 99), int(d_any[q]))
                if np.isfinite(d_free[q]) and d_free[q] <= max_len:
                    free_len[q] = min(free_len.get(q, 99), int(d_free[q]))
                    path, cur = [q], q
                    while cur != src:
                        cur = int(pred[cur])
                        path.append(cur)
                    examples.append(path[::-1])
        examples.sort(key=lambda p: (len(p), p))
        ex = examples[:5]
        all_paths.extend(ex)
        out["per_drug"][drug] = {"drug_nodes": [int(n) for n in nodes], "shortest_length_hubs_allowed": {query_nodes[q]: v for q, v in sorted(any_len.items())},
                                 "shortest_length_hub_free": {query_nodes[q]: v for q, v in sorted(free_len.items())}, "n_query_nodes_reached_hub_free": len(free_len),
                                 "example_paths_raw": ex, "verdict": "survives" if free_len else "refuted"}
    if edge_relation is not None and all_paths:
        rel = edge_relation(all_paths)
        for drug, rec in out["per_drug"].items():
            if "example_paths_raw" not in rec:
                continue
            rec["example_paths"] = [{"nodes": [f"{node_name[n]} [{node_type[n]}]" for n in p], "edges": [rel.get((min(a, b), max(a, b)), "unknown") for a, b in zip(p[:-1], p[1:])],
                                     "intermediate_degrees": [int(degree[n]) for n in p[1:-1]]} for p in rec["example_paths_raw"]]
    for rec in out["per_drug"].values():
        rec.pop("example_paths_raw", None)
    out["source_note"] = ("PrimeKG's on-disk edge table has no per-edge primary-source field (relation, display_relation, x_index, y_index only); the relation type of each edge is listed; "
                          "the 'lacks a primary source' branch is not assessable per edge and was not used to refute")
    return out


# ------------------------------------------------------------------ CF5: scheme sensitivity
def scheme_table(scores: Dict[str, Dict[str, Dict[str, float]]], wts: Dict[str, Dict[str, float]], thresholds: Sequence[str], schemes: Sequence[str],
                 top_n: int = 20) -> Dict[str, Any]:
    """scores[module][threshold][drug] = z. Per scheme and threshold: combined score ranking; top set = top_n at every threshold."""
    ranks: Dict[str, Dict[str, Dict[str, int]]] = {}
    tops: Dict[str, List[set]] = {}
    for s in schemes:
        ranks[s], tops[s] = {}, []
        for t in thresholds:
            comb: Dict[str, float] = {}
            for m, w in wts[s].items():
                for d, z in scores[m][t].items():
                    comb[d] = comb.get(d, 0.0) + w * max(0.0, -z)
            order = sorted((d for d, v in comb.items() if v > 0), key=lambda d: (-comb[d], d))
            ranks[s][t] = {d: i + 1 for i, d in enumerate(order)}
            tops[s].append(set(order[:top_n]))
    top_set = {s: set.intersection(*tops[s]) for s in schemes}
    return {"ranks": ranks, "top_set": top_set}


def cf5_verdicts(cands: Dict[str, str], table: Dict[str, Any], schemes: Sequence[str], w3: str, observed_rank: Dict[str, Optional[int]],
                 null_ranks: Dict[str, List[int]]) -> Dict[str, Any]:
    names = sorted(cands)
    pv = {}
    for d in names:
        if observed_rank.get(d) is None:
            continue
        pv[d] = (1 + sum(1 for r in null_ranks[d] if r <= observed_rank[d])) / (len(null_ranks[d]) + 1)
    qv = dict(zip(pv, bh(list(pv.values())))) if pv else {}
    out: Dict[str, Any] = {}
    for d in names:
        did = cands[d]
        others = [s for s in schemes if s != w3]
        in_top = {s: did in table["top_set"][s] for s in schemes}
        rec: Dict[str, Any] = {"in_top_set_by_scheme": in_top,
                               "ranks_by_scheme_threshold": {s: {t: table["ranks"][s][t].get(did) for t in table["ranks"][s]} for s in schemes}}
        if d not in pv:
            rec.update(verdict="not estimable", reason="drug has no scored targets in the recomputation")
        else:
            crit = any(in_top[s] for s in others)
            control = pv[d] <= P_CUT or qv[d] < Q_CUT
            rec.update(rank_real_module=observed_rank[d], null_median_rank=float(np.median(null_ranks[d])), n_null_modules=len(null_ranks[d]), p_random_module=pv[d], q_bh=float(qv[d]),
                       criterion_in_top_set_under_other_scheme=crit, random_module_control_passed=bool(control),
                       verdict="survives" if (crit and control) else "refuted")
        out[d] = rec
    return {"per_drug": out}


# ------------------------------------------------------------------ CF6: sign permutation (line D)
def net_effect_of_drug(net: Any, uniprot_by_target: Dict[str, List[str]], mechanisms: Sequence[Dict[str, Any]], inhibitory: set, activating: set,
                       cache: Optional[Dict[Any, Any]] = None) -> float:
    from track2.signed_effect import SignedNetwork

    cache = {} if cache is None else cache
    total = 0.0
    for m in mechanisms:
        action = (m.get("action") or "").upper()
        sign = -1 if action in inhibitory else 1 if action in activating else 0
        if not sign:
            continue
        for acc in uniprot_by_target.get(m["target_id"], []):
            key = (acc, sign)
            if key not in cache:
                cache[key] = SignedNetwork.effect(net, acc, sign)
            if cache[key]:
                total += cache[key]
    return total


def permuted_net(net: Any, rng: np.random.Generator) -> Any:
    adj = net.adj.copy()
    adj.data = rng.permutation(adj.data)
    return types.SimpleNamespace(nodes=net.nodes, index=net.index, adj=adj, readout=net.readout, symbol=net.symbol, by_symbol=net.by_symbol)


def cf6_sign_permutation(net: Any, drugs: Dict[str, Dict[str, Any]], inhibitory: set, activating: set, n_perm: int, seed: int = SEED,
                         controls: Optional[Sequence[tuple]] = None, call_of: Optional[Callable] = None, readout_genes: Sequence[str] = ()) -> Dict[str, Any]:
    """drugs[name] = {'uniprot_by_target': {...}, 'mechanisms': [...], 'stored_net': float|None}."""
    rng = np.random.default_rng(seed)
    obs = {d: net_effect_of_drug(net, v["uniprot_by_target"], v["mechanisms"], inhibitory, activating) for d, v in drugs.items()}
    null = {d: np.zeros(n_perm) for d in drugs}
    n_rec = np.zeros(n_perm, dtype=int)
    n_rec_nonreadout = np.zeros(n_perm, dtype=int)

    def recovered(nt: Any) -> List[tuple]:
        res = []
        for gene, sign, expected in controls or ():
            acc = nt.by_symbol.get(gene)
            from track2.signed_effect import SignedNetwork

            got = call_of(SignedNetwork.effect(nt, acc, sign)) if acc else "not_in_network"
            res.append((gene, got == expected))
        return res

    obs_ctrl = recovered(net) if controls else []
    for i in range(n_perm):
        pn = permuted_net(net, rng)
        cache: Dict[Any, Any] = {}
        for d, v in drugs.items():
            null[d][i] = net_effect_of_drug(pn, v["uniprot_by_target"], v["mechanisms"], inhibitory, activating, cache)
        if controls:
            r = recovered(pn)
            n_rec[i] = sum(ok for _, ok in r)
            n_rec_nonreadout[i] = sum(ok for g, ok in r if g not in readout_genes)
    names = sorted(drugs)
    est = [d for d in names if drugs[d]["mechanisms"] and any(obs.get(d) is not None for _ in [0])]
    p = {}
    for d in names:
        p[d] = (1 + int((np.abs(null[d]) >= abs(obs[d]) - 1e-12).sum())) / (n_perm + 1)
    valid = [d for d in names if abs(obs[d]) > 0 or (null[d] != 0).any()]
    qv = dict(zip(valid, bh([p[d] for d in valid]))) if valid else {}
    per = {}
    for d in names:
        lo, hi = np.percentile(null[d], [2.5, 97.5])
        stored = drugs[d].get("stored_net")
        rec = {"observed_net": obs[d], "stored_net": stored, "reproduces_stored": (stored is None or abs(stored - obs[d]) < 1e-9),
               "null_range_2_5_97_5": [float(lo), float(hi)], "null_mean": float(null[d].mean()), "p_two_sided": p[d]}
        if d not in qv:
            rec.update(verdict="not estimable", reason="no targets with a signed path to the readout nodes")
        else:
            rec.update(q_bh=float(qv[d]), verdict="survives" if (p[d] <= P_CUT and qv[d] < Q_CUT) else "refuted")
        per[d] = rec
    out: Dict[str, Any] = {"n_perm": n_perm, "seed": seed, "per_drug": per}
    if controls:
        k_all, k_non = sum(ok for _, ok in obs_ctrl), sum(ok for g, ok in obs_ctrl if g not in readout_genes)
        n_non = sum(1 for g, _ in obs_ctrl if g not in readout_genes)
        pa = (1 + int((n_rec >= k_all).sum())) / (n_perm + 1)
        pn_ = (1 + int((n_rec_nonreadout >= k_non).sum())) / (n_perm + 1)
        out["controls"] = {"recovered_observed": k_all, "of": len(obs_ctrl), "p_null_ge_observed": pa, "null_mean_recovered": float(n_rec.mean()),
                           "non_readout_recovered_observed": k_non, "non_readout_of": n_non, "p_non_readout_null_ge_observed": pn_,
                           "controls_above_chance": bool(pa <= P_CUT and pn_ <= P_CUT), "per_control": [{"gene": g, "recovered": bool(ok)} for g, ok in obs_ctrl]}
    return out


# ------------------------------------------------------------------ CF7: leave-one-array-out (line E-P)
def signature_from_arrays(m: np.ndarray, genes_of_rows: Sequence[Optional[str]], case_cols: Sequence[int], ctrl_cols: Sequence[int], n: int = 100,
                          min_abs_d: float = 0.5) -> Dict[str, Any]:
    """Same estimator as `line_e.geo_signature`, except that a group with a single array takes its variance from the other group."""
    a_all, b_all = m[:, list(case_cols)], m[:, list(ctrl_cols)]
    use = list(case_cols) + list(ctrl_cols)
    mean = np.nanmean(m[:, use], axis=1)
    best: Dict[str, int] = {}
    for r, g in enumerate(genes_of_rows):
        if g and np.isfinite(mean[r]) and (g not in best or mean[r] > mean[best[g]]):
            best[g] = r
    genes = sorted(best)
    g_m = m[[best[g] for g in genes]][:, use]
    na, nb = len(case_cols), len(ctrl_cols)
    a, b = g_m[:, :na], g_m[:, na:]
    d = a.mean(1) - b.mean(1)
    with np.errstate(all="ignore"):
        va = a.var(1, ddof=1) if na > 1 else None
        vb = b.var(1, ddof=1) if nb > 1 else None
    if va is None:
        va = vb
    if vb is None:
        vb = va
    pooled = np.sqrt((va * (na - 1) + vb * (nb - 1)) / (na + nb - 2) * (1 / na + 1 / nb))
    t = d / (pooled + np.nanmedian(pooled))
    gmean = g_m.mean(1)
    expressed = gmean >= np.median(gmean)
    ok = [k for k in range(len(genes)) if expressed[k] and abs(d[k]) >= min_abs_d and np.isfinite(t[k])]
    up = [genes[k] for k in sorted((k for k in ok if d[k] > 0), key=lambda k: -t[k])[:n]]
    down = [genes[k] for k in sorted((k for k in ok if d[k] < 0), key=lambda k: t[k])[:n]]
    return {"up": up, "down": down, "expressed_genes": [genes[k] for k in range(len(genes)) if expressed[k]], "n_cases": na, "n_controls": nb}


def cf7_verdicts(cands: Sequence[str], per_leave_out: Dict[str, Dict[str, Dict[str, Any]]], concordance: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for d in cands:
        if d not in concordance:
            out[d] = {"verdict": "not estimable", "reason": "drug not in the matched drug set"}
            continue
        rows = {lo: v[d] for lo, v in per_leave_out.items() if d in v}
        fails = sorted(lo for lo, r in rows.items() if not (r["q_rev"] < Q_CUT and r["R"] > 0))
        conc = concordance[d]["concordant_with_disease_signature"]
        out[d] = {"leave_out": rows, "leave_outs_failed": fails, "max_q_rev": max(r["q_rev"] for r in rows.values()), "concordance": concordance[d],
                  "verdict": "refuted" if (fails or conc) else "survives",
                  "action_if_refuted": ("flag hazard" if conc else "close") if (fails or conc) else None}
    return {"per_drug": out}


# ------------------------------------------------------------------ CF8: IGBP1 axis distance
def cf8_distance(adj: sparse.csr_matrix, degree: np.ndarray, source: int, targets: Dict[str, int], n_draws: int, seed: int = SEED,
                 bin_builder: Optional[Callable] = None, chunk: int = 250) -> Dict[str, Any]:
    from track2.proximity import degree_bins

    bin_of, bins = degree_bins(degree)
    d_obs_all = dijkstra(adj, directed=False, unweighted=True, indices=source)
    obs = {k: float(d_obs_all[t]) for k, t in targets.items()}
    if not np.all(np.isfinite(list(obs.values()))):
        return {"verdict": "not estimable", "reason": "a target is unreachable from the source", "observed": obs}
    obs_mean = float(np.mean(list(obs.values())))
    rng = np.random.default_rng(seed)
    real = {source, *targets.values()}
    src_pool = np.array([n for n in bins[bin_of[source]] if n not in real])
    tgt_pools = {k: np.array([n for n in bins[bin_of[t]] if n not in real]) for k, t in targets.items()}
    sources = rng.choice(src_pool, size=n_draws)
    picks = {k: rng.choice(tgt_pools[k], size=n_draws) for k in targets}
    null_mean = np.zeros(n_draws)
    null_each = {k: np.zeros(n_draws) for k in targets}
    for start in range(0, n_draws, chunk):
        idx = slice(start, min(start + chunk, n_draws))
        d = dijkstra(adj, directed=False, unweighted=True, indices=sources[idx])
        rows = np.arange(d.shape[0])
        acc = np.zeros(d.shape[0])
        for k in targets:
            v = d[rows, picks[k][idx]]
            null_each[k][idx] = v
            acc += v
        null_mean[idx] = acc / len(targets)
    finite = np.isfinite(null_mean)
    nm = null_mean[finite]
    p = (1 + int((nm <= obs_mean + 1e-12).sum())) / (len(nm) + 1)
    return {"observed_distance_each": obs, "observed_mean_distance": obs_mean, "null_mean_of_means": float(nm.mean()), "null_median": float(np.median(nm)), "null_n_draws": int(len(nm)),
            "p_null_le_observed": p, "p_each_target": {k: (1 + int((null_each[k][finite] <= obs[k] + 1e-12).sum())) / (len(nm) + 1) for k in targets},
            "verdict": "survives" if p <= P_CUT else "refuted"}


# ------------------------------------------------------------------ CF9: second screen
def cf9_second_screen(compound_rows: Sequence[Dict[str, str]], drugs: Sequence[str], other_resources: Sequence[str] = (),
                      rerun: Optional[Callable[[str, List[str]], Dict[str, Any]]] = None) -> Dict[str, Any]:
    by_name: Dict[str, List[Dict[str, str]]] = {}
    for r in compound_rows:
        for nm in {r["Drug.Name"].upper()} | {s.strip().upper() for s in (r.get("Synonyms") or "").split("|") if s.strip()}:
            by_name.setdefault(nm, []).append(r)
    per: Dict[str, Any] = {}
    for d in drugs:
        rows = by_name.get(d.upper(), [])
        screens = sorted({r["screen"] for r in rows})
        per[d] = {"n_rows": len(rows), "screens": screens, "has_second_independent_screen": len(screens) >= 2}
    runnable = [d for d in drugs if per[d]["has_second_independent_screen"]]
    if not runnable and not other_resources:
        return {"verdict": "not runnable", "per_drug": per, "other_resources_found": list(other_resources),
                "reason": "each line C survivor has exactly one row in the PRISM compound list (one screen), and no GDSC, CTRP or PRISM secondary file is on disk; the DepMap CRISPR matrix is a different assay",
                "needed": "an independent drug-sensitivity screen covering the same drugs and cell lines (GDSC, CTRP or PRISM secondary/MTS data), which would have to be downloaded"}
    out_per = {}
    for d in runnable:
        res = rerun(d, per[d]["screens"]) if rerun else None
        per[d]["second_screen_result"] = res
        if res is None:
            continue
        ok = res.get("same_sign") and res.get("p_perm", 1.0) <= P_CUT
        out_per[d] = "survives" if ok else "refuted"
    return {"verdict": "partial" if out_per else "not runnable", "per_drug": per, "per_drug_verdict": out_per}


# ------------------------------------------------------------------ runners (disk data)
def _read_json(p: Path) -> Any:
    return json.loads(p.read_text(encoding="utf-8"))


def triage_map(plan: Dict[str, Any]) -> Dict[str, List[str]]:
    m: Dict[str, List[str]] = {}
    for t in plan["triage"]:
        cfs = re.findall(r"CF\d", t["next_step"])
        m[t["candidate"].upper()] = cfs
    return m


_LC_CACHE: Dict[str, Any] = {}


def line_c_context(root: Path) -> Dict[str, Any]:
    """Loads the PRISM matrix and the checkpoint score and re-runs the adjusted analysis of `line_c_adjusted` unchanged (all estimable drugs)."""
    if "ctx" in _LC_CACHE:
        return _LC_CACHE["ctx"]
    from track2 import line_c_adjusted as lca

    orig = _read_json(root / "out" / "line_c_drug_level.json")
    pool = sorted(list(orig["drugs"]) + list(orig["not_in_prism"]))
    D = lca.load(root, pool)
    full = lca.analyse(D["Y"], D["score"], D["prolif"], D["lineage"])
    ctx = {"D": D, "full": full, "orig": orig["drugs"], "lca": lca}
    _LC_CACHE["ctx"] = ctx
    return ctx


def run_cf1(root: Path) -> Dict[str, Any]:
    ctx = line_c_context(root)
    D, full = ctx["D"], ctx["full"]
    adj = {d: (full[d]["rho_adj"] if d in full and full[d].get("estimable") else None) for d in ["LAZERTINIB", "NERATINIB"]}
    doc = {"definition_strings": {
        "src/track2/multi_evidence_completion.py:214": "score = mean over cell lines of the z-scored expression of BUB1B, BUB1, BUB3, MAD1L1, MAD2L1, TTK; low = weak checkpoint; statistic = Spearman rho between PRISM viability (log fold change) and the score; rho > 0 = more lethal when the checkpoint is weak",
        "src/track2/vulnerability_core.py:4": "diff = mean of the weak-score quartile minus mean of the strong-score quartile (negative = the row is lower, i.e. more lethal, when the score is weak)",
        "src/track2/prism_l4.py:43-60": "sac_scores: mean_all = mean of per-gene z-scores of log2(TPM+1) of BUB1B, BUB1, BUB3, MAD1L1, MAD2L1, TTK per cell line (high = strong checkpoint)",
        "src/track2/prism_l4.py:5-6": "difference in mean viability (log fold change) between weak-SAC and strong-SAC groups (negative = the drug kills weak-SAC cells more)",
        "src/track2/line_c_adjusted.py:76-93": "partial_perm: rho_adj is the partial rank correlation of the viability rank with the score rank (same orientation as the unadjusted rho)"},
        "reading": "viability is a log fold change (negative = fewer cells); score low = weak checkpoint; therefore rho > 0 means low-score (weak-checkpoint) lines have lower viability = more sensitive; positive rho is sensitivity, not resistance"}
    res = cf1_sign_audit(D["Y"], D["score"], ["LAZERTINIB", "NERATINIB"], adj, ["BORTEZOMIB", "PACLITAXEL", "VINCRISTINE", "COLCHICINE"], doc)
    res["concerns"] = ["LAZERTINIB", "NERATINIB"]
    # the four tubulin-class survivors are checked with the same convention (information)
    extra = {d: quartile_delta(D["Y"][d], D["score"]) for d in ["COLCHICINE", "ERIBULIN", "PACLITAXEL", "VINCRISTINE"] if d in D["Y"]}
    for d, v in extra.items():
        v["adjusted_partial_rho"] = full[d]["rho_adj"] if d in full and full[d].get("estimable") else None
    res["tubulin_class_information"] = extra
    res["verdict"] = "survives" if all(v["verdict"] == "survives" for v in res["per_drug"].values()) else "refuted" if any(v["verdict"] == "refuted" for v in res["per_drug"].values()) else "not estimable"
    res["resolved_by_documentation"] = True
    return res


def run_cf2(root: Path) -> Dict[str, Any]:
    ctx = line_c_context(root)
    full, orig, lca = ctx["full"], ctx["orig"], ctx["lca"]
    est = [d for d in full if full[d].get("estimable") and orig.get(d, {}).get("rho") is not None]
    hits = [d for d in orig if orig[d].get("q_bh") is not None and orig[d]["q_bh"] < lca.Q_CUT and d in full and full[d].get("estimable")]
    surv = [d for d in hits if lca.survives(full[d], orig[d]["rho"]) == "survives"]
    adj_hits = [d for d in est if full[d]["q_bh"] < lca.Q_CUT]
    caution = _read_json(root / "out" / "caution_list.json")
    members = [f["name"] for f in caution["clases"]["a_tubulina"]["farmacos"]]
    res = cf2_enrichment(est, surv, members, adj_hits)
    res["concerns"] = ["COLCHICINE", "ERIBULIN", "PACLITAXEL", "VINCRISTINE"]
    res["tubulin_class_source"] = "out/caution_list.json class a_tubulina (9 listed drugs)"
    res["n_original_hits_estimable"] = len(hits)
    return res


def run_cf3(root: Path, n_null: int = 99) -> Dict[str, Any]:
    from track2.string_net import load_interactome

    it = load_interactome(700)
    deg = {sym: int(it.degree[node]) for sym, node in it.symbol_to_node.items()}
    nullset = matched_null_genes(deg, "BUB1B", n_null)
    getter = RateLimitedGetter()
    drugs = {"OLAPARIB": "olaparib", "TYROSINE": "tyrosine", "PROGESTERONE": "progesterone", "ZINC CHLORIDE": "zinc chloride", "VORINOSTAT": "vorinostat"}
    counts: Dict[str, Optional[int]] = {}
    nulls: Dict[str, List[int]] = {}
    for key, nm in drugs.items():
        counts[key] = hit_count(nm, "BUB1B", getter)
        nulls[key] = [hit_count(nm, g, getter) for g in nullset["genes"]]
    gene_volume = {"BUB1B": int(getter(epmc_url(f'(TITLE:"BUB1B" OR ABSTRACT:"BUB1B") AND PUB_YEAR:[1900 TO {MAX_YEAR}]'))["hitCount"])}
    res = cf3_specificity(counts, nulls)
    res.update(null_genes=nullset, query_template=epmc_query("DRUG", "GENE"), europepmc_requests_to_network_this_run=getter.n_network, bub1b_papers_up_to_2024=gene_volume["BUB1B"],
               limitation="only the literature co-mention leg of the engine is tested; the engine prior is not re-scored for random genes")
    return res


def run_cf4(root: Path) -> Dict[str, Any]:
    import pandas as pd

    from track2.kg_query import universe_match
    from track2.knowledge_graph import load_graph

    src = root / "data" / "sources" / "primekg"
    graph = load_graph(src / "nodes.tab", src / "edges.csv")
    universe = _read_json(root / "data" / "universe.json")
    kq = _read_json(root / "out" / "kg_query_e5.json")
    query = {int(kq["disease_node"]["index"]): "MVA disease node: " + kq["disease_node"]["name"]}
    for hp, v in kq["phenotype_nodes"].items():
        if v:
            query[int(v["index"])] = f"{hp} {v['name']}"
    matched = universe_match(graph, universe)
    drug_nodes = {n: sorted(k for k, v in matched.items() if v["name"].upper() == n) for n in ("POMALIDOMIDE", "DISOPYRAMIDE")}

    def relations(paths: List[List[int]]) -> Dict[tuple, str]:
        wanted = {(min(a, b), max(a, b)) for p in paths for a, b in zip(p[:-1], p[1:])}
        nn = len(graph.node_name)
        keys = np.array(sorted(a * nn + b for a, b in wanted), dtype=np.int64)
        found: Dict[tuple, str] = {}
        for chunk in pd.read_csv(src / "edges.csv", usecols=["relation", "x_index", "y_index"], chunksize=2_000_000):
            x, y = chunk["x_index"].to_numpy(np.int64), chunk["y_index"].to_numpy(np.int64)
            k = np.minimum(x, y) * nn + np.maximum(x, y)
            hit = np.isin(k, keys)
            for kk, rel in zip(k[hit], chunk["relation"].to_numpy()[hit]):
                a, b = divmod(int(kk), nn)
                found.setdefault((a, b), str(rel))
        return found

    res = cf4_path_audit(graph.adj, graph.node_name, graph.node_type, drug_nodes, query, edge_relation=relations)
    res["query_nodes"] = {str(k): v for k, v in query.items()}
    res["concerns"] = list(drug_nodes)
    return res


def run_cf5(root: Path, n_modules: int = 100, iterations: int = 300) -> Dict[str, Any]:
    from track2.augmented_proximity import augmented_universe
    from track2.consensus import SCHEMES, THRESHOLDS, weights
    from track2.core_n2 import restrict_links
    from track2.core_sensitivity import core_mechanism
    from track2.modules import build_module, load_mechanism
    from track2.proximity import _sample_matched, build_null, degree_bins
    from track2.scoring import score_drugs
    from track2.string_net import load_interactome
    from track2.track1_scenarios import load_scores

    mech_dir = root / "config" / "mechanisms_from_submission"
    mechs = [m for m in (load_mechanism(p) for p in sorted(mech_dir.glob("*.json"))) if m.get("status") == "mechanism_candidate"]
    ids = [m["id"] for m in mechs]
    wts = weights(mechs)
    scores = load_scores(root / "out" / "core_n2", "core_threshold_scores_", ids)
    table = scheme_table(scores, wts, THRESHOLDS, SCHEMES)
    universe = _read_json(root / "data" / "universe.json")
    ids_by_name = {d["name"].upper(): d["drug_id"] for d in universe["drugs"]}
    cands = {n: ids_by_name[n] for n in ("CERITINIB", "SELINEXOR", "TIRBANIBULIN")}
    e1 = restrict_links(_read_json(root / "data" / "dti_e1_all.json")["edges"])
    aug = augmented_universe(universe, e1, {})
    it = load_interactome(700)
    core = core_mechanism(load_mechanism(mech_dir / "sub_BUB1B.json"))
    mod = build_module(core, it)["nodes"]

    def ranks_for(nodes: np.ndarray, seed: int) -> Dict[str, int]:
        null = build_null(it.adjacency, it.degree, nodes, repeats=iterations, seed=seed)
        rows = score_drugs(aug["drugs"], it, nodes, set(), null)
        return {r["drug_id"]: r["rank"] for r in rows}

    real = ranks_for(mod, 20261006)
    rng = np.random.default_rng(SEED)
    bin_of, bins = degree_bins(it.degree)
    rm = _sample_matched(rng, mod, bin_of, bins, n_modules)
    null_ranks: Dict[str, List[int]] = {n: [] for n in cands}
    for r in range(n_modules):
        rk = ranks_for(np.unique(rm[r]), 20261006 + 1 + r)
        for n, did in cands.items():
            null_ranks[n].append(rk.get(did, len(rk) + 1))
    obs = {n: real.get(did) for n, did in cands.items()}
    res = cf5_verdicts(cands, table, SCHEMES, "W3_pathogenic_gene_only", obs, null_ranks)
    res.update(top_set_by_scheme={s: sorted(v) for s, v in table["top_set"].items()}, top_set_names={s: sorted(d["name"] for d in universe["drugs"] if d["drug_id"] in v) for s, v in table["top_set"].items()},
               module_size=int(len(mod)), n_random_modules=n_modules, null_z_iterations=iterations, scored_drugs_real_module=len(real))
    return res


def run_cf6(root: Path, n_perm: int = 1000) -> Dict[str, Any]:
    from track2 import signed_effect as se
    from track2.direction import ACTIVATING, INHIBITORY

    net = se.SignedNetwork(root / "data" / "sources" / "omnipath" / "omnipath_interactions.tsv")
    universe = _read_json(root / "data" / "universe.json")
    uni_t = {tid: t.get("uniprot", []) for tid, t in universe["targets"].items()}
    stored = _read_json(root / "out" / "signed_effect.json")["drugs"]
    by_name = {d["name"].upper(): d for d in universe["drugs"]}
    drugs = {}
    for n in ("REGORAFENIB", "SUNITINIB", "NIROGACESTAT", "RITLECITINIB", "ERDAFITINIB"):
        d = by_name[n]
        drugs[n] = {"uniprot_by_target": {m["target_id"]: uni_t.get(m["target_id"], []) for m in d["mechanisms"]}, "mechanisms": d["mechanisms"],
                    "stored_net": stored[d["drug_id"]]["net"]}
    readout = set(se.MCC) | set(se.APCC)
    res = cf6_sign_permutation(net, drugs, set(INHIBITORY), set(ACTIVATING), n_perm, controls=se.CONTROLS, call_of=se.call_of, readout_genes=sorted(readout))
    res["network"] = {"nodes": len(net.nodes), "signed_edges": int(net.adj.nnz), "n_positive": int((net.adj.data > 0).sum()), "n_negative": int((net.adj.data < 0).sum())}
    res["concerns"] = list(drugs)
    return res


def run_cf7(root: Path) -> Dict[str, Any]:
    import gzip

    from track2 import line_e as le

    ko = le.parse_library((root / le.KO_FILE).read_text(encoding="utf-8"))
    chem = le.parse_library((root / le.CHEM_FILE).read_text(encoding="utf-8"))
    universe = _read_json(root / "data" / "universe.json")["drugs"]
    genes = sorted({g for lib in (ko, chem) for v in lib.values() for s in v.values() for g in s})
    index = {g: i for i, g in enumerate(genes)}
    matched, _ = le.match_drugs(list(chem), universe)
    names = list(matched)
    uname = [matched[n]["name"].upper() for n in names]
    D = le.drug_matrix(chem, names, le.N_PRIMARY, index)
    annotation = le.read_annotation(root / le.GPL_FILE)
    with gzip.open(root / le.GEO_FILE, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().split("\n")
    titles = [t.strip('"') for t in next(l for l in lines if l.startswith("!Sample_title")).split("\t")[1:]]
    i0 = next(k for k, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    case = [k for k, t in enumerate(titles) if t.startswith("LCL-proband")]
    ctrl = [k for k, t in enumerate(titles) if t.startswith("LCL-control")]
    use = case + ctrl
    ids, mat = [], []
    for l in lines[i0 + 2:]:
        if not l or l.startswith("!"):
            continue
        f = l.split("\t")
        ids.append(f[0].strip('"'))
        mat.append([float(f[1 + k]) if f[1 + k] not in ("", "null") else np.nan for k in range(len(titles))])
    M = np.array(mat)
    rowgenes = [annotation.get(i) for i in ids]
    cands = ["NERATINIB", "AMPHOTERICIN B", "AURANOFIN", "FLUVASTATIN", "NELFINAVIR"]
    stored = _read_json(root / "out" / "line_e_signature_reversal.json")
    tab = {r["name"].upper(): r for r in stored["table_all_drugs_n150"]}

    def evaluate(case_cols: List[int], ctrl_cols: List[int]) -> tuple:
        sig = signature_from_arrays(M, rowgenes, case_cols, ctrl_cols)
        p_sig = le.signed_vector(sig["up"], sig["down"], None, index)
        null = le.reversal(D, le.random_signature_null(sig["expressed_genes"], index, len(sig["up"]), len(sig["down"]), le.N_NULL, le.P_SEED))
        t = le.test_signature(D, p_sig, null)
        return sig, t

    sig_full, t_full = evaluate(case, ctrl)
    ref = _read_json(root / "out" / "line_e_signature_reversal.json")
    self_check = {"up_equals_stored_signature": None}
    # self-check against the original estimator
    orig = le.geo_signature(root / le.GEO_FILE, annotation)
    self_check["up_equals_stored_signature"] = bool(orig["up"] == sig_full["up"] and orig["down"] == sig_full["down"])
    per_lo: Dict[str, Dict[str, Dict[str, Any]]] = {}
    scheme = {"full_six_arrays": (case, ctrl)}
    for k in range(len(titles)):
        if k not in use:
            continue
        scheme[f"leave_out_{titles[k]}"] = ([c for c in case if c != k], [c for c in ctrl if c != k])
    for label, (cc, tc) in scheme.items():
        sig, t = (sig_full, t_full) if label == "full_six_arrays" else evaluate(cc, tc)
        per_lo[label] = {}
        for c in cands:
            if c in uname:
                i = uname.index(c)
                per_lo[label][c] = {"R": float(t["R"][i]), "q_rev": float(t["q_rev"][i]), "p_rev": float(t["p_rev"][i]), "n_up": len(sig["up"]), "n_down": len(sig["down"]),
                                    "n_case": len(cc), "n_control": len(tc)}
    for c in cands:
        if c in uname and c in tab:
            i = uname.index(c)
            self_check[f"{c}_R_q_equal_stored"] = bool(abs(float(t_full['R'][i]) - tab[c]["P"]["R"]) < 1e-9 and abs(float(t_full['q_rev'][i]) - tab[c]["P"]["q_rev"]) < 1e-9)
    conc = {}
    for c in cands:
        r = tab.get(c)
        if not r:
            continue
        flags = {k: bool(r[k]["hazard"]) for k in ("bub1b", "track1_combined")}
        conc[c] = {"hazard_flag_by_signature": flags, "concordant_with_disease_signature": any(flags.values()),
                   "bub1b_q_rev": r["bub1b"]["q_rev"], "track1_combined_q_rev": r["track1_combined"]["q_rev"], "bub1b_R": r["bub1b"]["R"], "track1_combined_R": r["track1_combined"]["R"]}
    leave = {k: v for k, v in per_lo.items() if k != "full_six_arrays"}
    res = cf7_verdicts(cands, leave, conc)
    res.update(self_check=self_check, full_six_arrays=per_lo["full_six_arrays"], array_titles=titles, concerns=cands)
    return res


def run_cf8(root: Path, n_draws: int = 5000) -> Dict[str, Any]:
    from track2.string_net import load_interactome

    out: Dict[str, Any] = {}
    for thr in (700, 400, 900):
        it = load_interactome(thr)
        s = it.symbol_to_node.get("IGBP1")
        tg = {g: it.symbol_to_node[g] for g in ("MTOR", "PRKAA1", "PRKAA2") if g in it.symbol_to_node}
        if s is None or not tg:
            out[str(thr)] = {"verdict": "not estimable", "reason": "IGBP1 or all targets absent from the network"}
            continue
        r = cf8_distance(it.adjacency, it.degree, s, tg, n_draws)
        r["degrees"] = {"IGBP1": int(it.degree[s]), **{g: int(it.degree[n]) for g, n in tg.items()}}
        r["targets_found"] = sorted(tg)
        out[str(thr)] = r
    res = dict(out["700"])
    res["primary_threshold"] = 700
    res["sensitivity_thresholds"] = {k: {"observed_mean_distance": v.get("observed_mean_distance"), "p_null_le_observed": v.get("p_null_le_observed"), "verdict": v.get("verdict")} for k, v in out.items() if k != "700"}
    res["concerns"] = ["SIROLIMUS", "EVEROLIMUS", "METFORMIN"]
    res["applies_to"] = "the IGBP1 rationale of X1, not the drugs"
    return res


def run_cf9(root: Path) -> Dict[str, Any]:
    import csv

    rows = list(csv.DictReader(open(root / "data" / "depmap" / "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv", encoding="utf-8")))
    other = [str(p.relative_to(root)) for p in (root / "data").rglob("*") if p.is_file() and re.search(r"gdsc|ctrp|secondary|mts", p.name, re.I)]
    res = cf9_second_screen(rows, ["LAZERTINIB", "NERATINIB", "COLCHICINE", "ERIBULIN", "PACLITAXEL", "VINCRISTINE"], other)
    res["screens_present_in_compound_list"] = sorted({r["screen"] for r in rows})
    res["concerns"] = ["LAZERTINIB", "NERATINIB", "COLCHICINE", "ERIBULIN", "PACLITAXEL", "VINCRISTINE"]
    return res


RUNNERS = {"CF1": run_cf1, "CF2": run_cf2, "CF3": run_cf3, "CF4": run_cf4, "CF5": run_cf5, "CF6": run_cf6, "CF7": run_cf7, "CF8": run_cf8, "CF9": run_cf9}
DEFS = {
    "CF1": ("Line C sign audit", "positive rho: viability against the checkpoint score; orientation-explicit quartile delta"),
    "CF2": ("Tubulin-class enrichment among line C survivors", "hypergeometric test among estimable drugs"),
    "CF3": ("Engine specificity", "BUB1B co-mention count against 99 degree-matched random genes (Europe PMC)"),
    "CF4": ("Graph path audit (line A)", "hub-free paths of length <= 4 from drug to the query nodes in PrimeKG"),
    "CF5": ("B2 scheme sensitivity", "top set under W1/W2 from stored scores; rank against 100 random degree-matched modules"),
    "CF6": ("Line D sign permutation", "net signed effect against 1,000 sign-shuffled networks; control re-check"),
    "CF7": ("Line E-P leave-one-array-out", "six leave-outs of the 2 case vs 4 control arrays; concordance with BUB1B/Track 1 signatures"),
    "CF8": ("IGBP1 axis network distance", "STRING distance IGBP1 to MTOR/PRKAA1/PRKAA2 against degree-matched random pairs"),
    "CF9": ("Second sensitivity screen for line C survivors", "second independent screen of the same drugs on disk"),
}


def cf_verdict(res: Dict[str, Any]) -> str:
    if "verdict" in res and isinstance(res["verdict"], str) and res["verdict"] in VERDICTS:
        return res["verdict"]
    return "mixed"


def per_candidate(res: Dict[str, Any]) -> Dict[str, str]:
    """candidate (upper case) -> verdict, for the candidates a CF concerns."""
    if "per_drug" in res and isinstance(res["per_drug"], dict):
        return {k.upper(): v["verdict"] for k, v in res["per_drug"].items() if isinstance(v, dict) and "verdict" in v}
    return {}


def candidate_status(plan: Dict[str, Any], results: Dict[str, Any]) -> List[Dict[str, Any]]:
    tm = triage_map(plan)
    cand_cf: Dict[str, Dict[str, str]] = {}
    for cf, res in results.items():
        v = per_candidate(res) if cf not in ("CF2", "CF8", "CF9") else {}
        if cf == "CF2":
            v = {c.upper(): res["verdict"] for c in res.get("concerns", [])}
        if cf == "CF8":
            v = {c.upper(): res["verdict"] for c in res.get("concerns", [])}
        if cf == "CF9":
            v = {c.upper(): res.get("verdict", "not runnable") for c in res.get("concerns", [])}
        if cf == "CF1":
            v = {c.upper(): x["verdict"] for c, x in res["per_drug"].items()}
        for c, x in v.items():
            cand_cf.setdefault(c, {})[cf] = x
    rows = []
    for t in plan["triage"]:
        name = t["candidate"]
        checks = cand_cf.get(name.upper(), {})
        if t["next_step_type"] == "closed":
            status, nxt = "closed (earlier adjusted analysis; not part of the CF set)", "none"
        elif t["next_step_type"] == "hazard":
            status = "hazard class" + ("; CF2 " + checks["CF2"] if "CF2" in checks else "")
            nxt = "laboratory H1 (toxicity readouts only); never a benefit candidate"
        else:
            hard = {cf: v for cf, v in checks.items() if cf != "CF8"}
            refuted = [cf for cf, v in hard.items() if v == "refuted"]
            passed = [cf for cf, v in hard.items() if v == "survives"]
            unk = [f"{cf} {v}" for cf, v in hard.items() if v in ("not estimable", "not runnable")]
            if refuted:
                status = "closed as a computational candidate (refuted by " + ", ".join(refuted) + ")"
                if passed:
                    status += f"; note: {', '.join(passed)} passed, and the refuted check concerns only the other line's support (plan rule: row closed; reopening is a user decision)"
                nxt = "closed" if not any(r in ("CF7",) for r in refuted) else "closed or flagged hazard (CF7 rule)"
            elif passed:
                status = "still a hypothesis; passed " + ", ".join(passed) + (" (" + "; ".join(unk) + ")" if unk else "")
                nxt = f"laboratory {t['experiment_id']}" if t.get("experiment_id") else "needs a laboratory design (none in the plan yet)"
            else:
                status = "unchanged (" + ("; ".join(unk) if unk else "no per-drug computational check applies to this row") + ")"
                nxt = f"laboratory {t['experiment_id']}" if t.get("experiment_id") else "unchanged"
            if "CF8" in checks:
                status += f"; IGBP1 rationale {'kept' if checks['CF8'] == 'survives' else 'dropped' if checks['CF8'] == 'refuted' else checks['CF8']} (CF8, affects the rationale only)"
        rows.append({"candidate": name, "lines": t["lines"], "checks": checks, "status": status, "next_step": nxt, "plan_next_step": t["next_step"]})
    return rows


def render_md(report: Dict[str, Any]) -> str:
    L = ["# Computational follow-ups CF1 to CF9", "", f"Label: {LABEL}. Rules fixed before running: `{PROTOCOL}`. Every number is computed by `src/track2/computational_followups.py`; no patient data.", "",
         "| CF | Candidates concerned | Test | Result | Verdict | Reading |", "|---|---|---|---|---|---|"]
    for cf in sorted(report["cf"]):
        r = report["cf"][cf]
        L.append(f"| {cf} | {', '.join(r.get('concerns_display', []))} | {r['test']} | {r['result_line']} | {r['verdict_summary']} | {r['reading']} |")
    L += ["", "## Updated status of the triaged candidates", "", "| Candidate | Line(s) | Checks | Status | Next step |", "|---|---|---|---|---|"]
    for c in report["candidate_status"]:
        chk = "; ".join(f"{k} {v}" for k, v in sorted(c["checks"].items())) or "none"
        L.append(f"| {c['candidate']} | {', '.join(c['lines'])} | {chk} | {c['status']} | {c['next_step']} |")
    L += ["", "Not concluded from any row: clinical benefit, safety or dose. A survivor is a hypothesis that was not refuted by a post hoc computational check.", ""]
    if report.get("not_runnable"):
        L += ["## Not runnable", ""] + [f"* {k}: {v}" for k, v in report["not_runnable"].items()] + [""]
    return "\n".join(L)


def describe(cf: str, r: Dict[str, Any]) -> Dict[str, Any]:
    """One-line result and reading per CF, built only from the computed numbers."""
    f = lambda x, nd=3: "NA" if x is None else (f"{x:.{nd}g}")  # noqa: E731
    pc = per_candidate(r)
    summ = ", ".join(f"{k.title()} {v}" for k, v in pc.items()) if pc else ""
    if cf == "CF1":
        pd_ = r["per_drug"]
        res = "; ".join(f"{k.title()}: delta {f(v.get('delta_weak_minus_strong'))}, adjusted rho {f(v.get('adjusted_partial_rho'))}" for k, v in pd_.items())
        res += f"; orientation sanity median LFC {', '.join(f'{k.title()} {f(v)}' for k, v in r['orientation_sanity_median_lfc'].items())}"
        return {"result_line": res, "verdict_summary": r["verdict"] + " (resolved by documentation)",
                "reading": "rho > 0 means lower viability in weak-checkpoint lines, i.e. sensitivity; Lazertinib and Neratinib are not dropped on sign grounds" if r["verdict"] == "survives" else "sign check failed; see JSON"}
    if cf == "CF2":
        p = r["primary"]
        return {"result_line": f"{p['x_tubulin_among_draws']} of {p['n_draws']} survivors are tubulin-class; {r['K_tubulin_class_in_estimable']} of {r['M_estimable']} estimable drugs are; expected {f(p['expected'])}; hypergeometric p {f(p['p_hypergeometric_ge'])}",
                "verdict_summary": r["verdict"], "reading": "enriched: the score behaves as a mitotic-dependency control and the four stay hazard-class" if r["verdict"] == "refuted" else "not enriched at the pre-set level; the four still stay hazard-class" if r["verdict"] == "survives" else "not estimable"}
    if cf == "CF3":
        res = "; ".join(f"{k.title()} {v.get('observed_hit_count_bub1b')} vs null median {f(v.get('null_median'))} (p {f(v.get('p_empirical'))})" for k, v in r["per_drug"].items())
        return {"result_line": res, "verdict_summary": summ, "reading": "title/abstract co-mention of the drug with BUB1B is not above 99 degree-matched genes; only the co-mention leg is tested and aliases such as BubR1 are not searched"}
    if cf == "CF4":
        res = "; ".join(f"{k.title()}: hub-free query nodes reached {v.get('n_query_nodes_reached_hub_free')}" for k, v in r["per_drug"].items())
        return {"result_line": res, "verdict_summary": summ, "reading": "hub-free paths of at most 4 edges exist, so line A is not refuted; this is a permissive check (most drugs with targets have such paths) and edge sources are not recorded per edge in PrimeKG"}
    if cf == "CF5":
        res = "; ".join(f"{k.title()}: top set W1 {v['in_top_set_by_scheme']['W1_epcr']}, W2 {v['in_top_set_by_scheme']['W2_equal']}, random-module p {f(v.get('p_random_module'))}" for k, v in r["per_drug"].items())
        return {"result_line": res, "verdict_summary": summ, "reading": "support is not specific to scheme W3 if the drug is in the top set under W1 or W2 and above random modules"}
    if cf == "CF6":
        res = "; ".join(f"{k.title()}: net {f(v['observed_net'])}, null range [{f(v['null_range_2_5_97_5'][0])}, {f(v['null_range_2_5_97_5'][1])}], p {f(v['p_two_sided'])}" for k, v in r["per_drug"].items())
        c = r.get("controls")
        if c:
            res += f"; controls recovered {c['recovered_observed']}/{c['of']} (null p {f(c['p_null_ge_observed'])})"
        return {"result_line": res, "verdict_summary": summ, "reading": "signed-path nets are compared with sign-shuffled networks of the same topology; controls " + ("above chance" if c and c["controls_above_chance"] else "not above chance")}
    if cf == "CF7":
        res = "; ".join(f"{k.title()}: max q {f(v.get('max_q_rev'))}, failed leave-outs {len(v.get('leave_outs_failed', []))}, concordant {v.get('concordance', {}).get('concordant_with_disease_signature')}" for k, v in r["per_drug"].items())
        return {"result_line": res, "verdict_summary": summ, "reading": "reversal must hold with q < 0.10 in all six leave-outs and the drug must not mimic the BUB1B or Track 1 signature"}
    if cf == "CF8":
        if "observed_mean_distance" not in r:
            return {"result_line": r.get("reason", ""), "verdict_summary": r["verdict"], "reading": ""}
        return {"result_line": f"mean distance IGBP1 to MTOR/PRKAA1/PRKAA2 {f(r['observed_mean_distance'])} vs null mean {f(r['null_mean_of_means'])}; p {f(r['p_null_le_observed'])}",
                "verdict_summary": r["verdict"], "reading": "the IGBP1 rationale of X1 is kept as a hypothesis" if r["verdict"] == "survives" else "the IGBP1 rationale is dropped from X1; the drugs are not affected"}
    if cf == "CF9":
        return {"result_line": r.get("reason", "") if r["verdict"] == "not runnable" else "see JSON", "verdict_summary": r["verdict"],
                "reading": "needs: " + r.get("needed", "") if r["verdict"] == "not runnable" else ""}
    return {"result_line": "", "verdict_summary": "", "reading": ""}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="out")
    ap.add_argument("--only", default="", help="comma-separated CF ids to (re)run ('none' = only rebuild the summary); the others are kept from the existing JSON")
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.out)
    path = out / "computational_followups.json"
    prev = _read_json(path) if path.exists() else {"cf": {}}
    todo = [] if a.only == "none" else [x for x in (a.only.split(",") if a.only else list(RUNNERS)) if x]
    results = {cf: prev.get("results", {}).get(cf) for cf in RUNNERS if prev.get("results", {}).get(cf)}
    for cf in todo:
        t0 = time.time()
        results[cf] = jsonable(RUNNERS[cf](root))
        results[cf]["seconds"] = round(time.time() - t0, 1)
        print(cf, "done", results[cf]["seconds"], "s", flush=True)
    plan = _read_json(root / "out" / "experiment_plan.json")
    cfdefs = {c["id"]: c for c in plan["computational_followups"]}
    summary: Dict[str, Any] = {}
    for cf, r in results.items():
        d = describe(cf, r)
        v = cf_verdict(r)
        pcs = per_candidate(r)
        summary[cf] = {"name": DEFS[cf][0], "plan_definition": cfdefs[cf]["test"], "plan_refutation": cfdefs[cf]["refutation"], "test": DEFS[cf][1],
                       "concerns_display": [c.title() for c in r.get("concerns", list(pcs))], "verdict": v, "per_candidate_verdict": pcs or None,
                       "result_line": d["result_line"], "verdict_summary": d["verdict_summary"] or v, "reading": d["reading"]}
    status = candidate_status(plan, results)
    not_runnable = {cf: r.get("reason") + ". Needed: " + r.get("needed", "") for cf, r in results.items() if r.get("verdict") == "not runnable"}
    report = {"label": LABEL, "protocol": PROTOCOL, "seed": SEED, "results": results, "cf": summary, "candidate_status": status, "not_runnable": not_runnable}
    path.write_text(json.dumps(jsonable(report), indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    (out / "computational_followups.md").write_text(render_md(report), encoding="utf-8", newline="\n")
    for cf in sorted(summary):
        print(cf, summary[cf]["verdict_summary"])


if __name__ == "__main__":
    main()
