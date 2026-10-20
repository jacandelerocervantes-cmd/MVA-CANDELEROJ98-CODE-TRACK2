"""Completion of three data gaps of the multi-evidence integration (`protocol/addendum_multi_evidence.md`, section 11).

    PYTHONPATH=src python -m track2.multi_evidence_completion approval|engine|linec [--out out]

approval: the project's FDA approval check for candidates without an approval record (openFDA only; stored).
engine:   score and rank of every drug of the universe under the frozen B0 engine (no refit; Europe PMC caches only; stored).
linec:    drug-level association between PRISM sensitivity and the checkpoint-weakness score of the cell lines, with a permutation null and BH (stored).

No existing module is modified. Gene symbols, drug names and categorical outputs only. Hypotheses, not evidence of efficacy.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy import stats

SEED = 20261018
N_PERM = 1000
MIN_LINES = 100
Q_CUT = 0.10
RHO_SCALE = 0.3


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


# ---------------------------------------------------------------- approval

def complete_approval(names: list[str], check: Callable[[str], dict[str, Any]], prior: dict[str, Any] | None = None) -> dict[str, Any]:
    """Run the approval check for each name (answers already stored in `prior` are reused, so a rerun asks nothing twice)."""
    rows = dict((prior or {}).get("drugs", {}))
    for n in sorted(set(names)):
        if n not in rows:
            rows[n] = check(n)
    confirmed = sorted(n for n, r in rows.items() if r["confirmed"])
    return {"source": "Drugs@FDA through openFDA (https://api.fda.gov/drug/drugsfda.json)", "rule": "an NDA, BLA or ANDA record with a product marketed as Prescription or Over-the-counter",
            "code": "track2.fda_approval_check.check (unmodified)", "n_drugs": len(rows), "n_confirmed": len(confirmed), "confirmed": confirmed,
            "not_confirmed": sorted(set(rows) - set(confirmed)), "drugs": rows}


# ---------------------------------------------------------------- engine

def rank_desc(score: np.ndarray, names: list[str], mask: np.ndarray) -> dict[str, int]:
    """Rank (1 = best) of the masked drugs by decreasing score, ties by name."""
    idx = [i for i in range(len(names)) if mask[i]]
    idx.sort(key=lambda i: (-score[i], names[i]))
    return {names[i]: k + 1 for k, i in enumerate(idx)}


def engine_full_ranking(root: Path) -> dict[str, Any]:
    """Frozen B0 score of every drug of the universe, offline. Mirrors the scoring steps of `engine_apply.run` (which is not modified and not called,
    because it writes its own outputs); a Europe PMC cache miss raises."""
    from track2 import engine_apply as ea
    from track2 import empirical_literature_control as elc
    from track2 import engine_benchmark as eb
    from track2 import engine_case_study as cs
    from track2 import engine_score as es
    from track2 import engine_score_v2 as v2
    from track2 import engine_v4 as v4

    def no_net(*a: Any, **k: Any) -> Any:
        raise RuntimeError("Europe PMC cache miss (this step uses caches only)")

    landscape = json.loads((root / ea.DISEASES["mva"]["landscape"]).read_text(encoding="utf-8"))
    cur, _ = ea.curated_genes(landscape)
    symbols = [g["symbol"] for g in cur]
    model = v4.thaw_v4(json.loads((root / "out" / "engine" / "engine_v4_frozen.json").read_text(encoding="utf-8")))  # raises if the model hash does not match
    drugs = eb.load_universe(root / "data" / "universe.json")
    universe = eb.universe_index(drugs)
    names = sorted(universe)
    dre = elc.drug_regex(names)
    gate = es.name_gate(names)
    guard = v2.RequestGuard(no_net, None)
    cdir = root / "data" / "cache_engine"
    targets, ppi = es.load_prior_graph(root, cdir / "primekg_prior_edges.json")
    pop, tot = es.fetch_popularity(names, None, cdir / "epmc_popularity_all.json", guard)
    papers = es.fetch_papers(symbols, ea.MAX_YEAR, cdir / "epmc", guard)
    lit = es.literature_counts(papers, symbols, dre, gate)
    prior = es.prior_features(symbols, names, targets, ppi)
    X = es.disease_matrix(lit, names, pop, tot, prior)
    sc = np.asarray(v4.predict_frozen_v4(model, "B0", X)["score"], dtype=float)
    # filter chain by label cache only (no service): FDA label small molecules, as in the stored chain
    from track2 import engine_benchmark_v2 as b2mod
    label_only = {n: cs.has_fda_label(n, b2mod.fetch_labels(n, getter=ea._no_network)) for n in names}
    return {"names": names, "score": sc, "universe": universe, "label_only": label_only, "cs": cs}


def full_ranking_record(root: Path, stored_apply: dict[str, Any], lookup: Callable[[str], dict[str, Any]] | None = None) -> dict[str, Any]:
    """Ranks for every drug: among all drugs, and among FDA small molecules (label match, plus Drugs@FDA for the top 100 without a label match, as in the
    engine application). The stored top 20 and the nine report positions of `apply_mva.json` are used as a consistency check."""
    from track2 import fda_approved

    e = engine_full_ranking(root)
    names, sc, universe, label_only, cs = e["names"], e["score"], e["universe"], e["label_only"], e["cs"]
    lookup_fn = lookup or fda_approved.lookup
    order_all = sorted(range(len(names)), key=lambda i: (-sc[i], names[i]))
    check = sorted({names[i] for i in order_all[:cs.RESCUE_TOP] if not label_only[names[i]]})
    fda: dict[str, dict[str, Any]] = {}
    for n in check:
        try:
            fda[n] = lookup_fn(universe[n]["name"])
        except Exception as err:  # noqa: BLE001 - a service error is never read as approved
            fda[n] = {"eligible": None, "error": repr(err)[:120]}
    has_label = {n: label_only[n] or bool((fda.get(n) or {}).get("eligible")) for n in names}
    chain = cs.filter_chain(names, universe, has_label)
    sm = set(chain["small_molecule_names"])
    r_all = rank_desc(sc, names, np.ones(len(names), dtype=bool))
    r_sm = rank_desc(sc, names, np.array([n in sm for n in names]))
    stored_top = [(r["name"].upper(), r["position"]) for r in stored_apply.get("ranking_top20_fda_small_molecule", [])]
    mine_top = sorted(((n.upper(), r) for n, r in r_sm.items() if r <= 20), key=lambda x: x[1])
    consistent = [t for t in stored_top] == [(n, r) for n, r in mine_top][:len(stored_top)] if stored_top else None
    return {"label": "frozen B0 engine, exploratory; a prior from literature co-mention and graph targets; not a candidate list and not evidence of efficacy",
            "model": "B0 (frozen; no refit)", "n_universe": len(names), "n_fda_small_molecule": len(r_sm),
            "consistent_with_stored_top20": consistent, "n_label_lookups_online_or_cached": len(check),
            "drugs": {n.upper(): {"score": round(float(sc[i]), 6), "rank_all": r_all[n], "rank_fda_small_molecule": r_sm.get(n), "fda_small_molecule": n in sm}
                      for i, n in enumerate(names)}}


# ---------------------------------------------------------------- line C at drug level

def spearman_perm(y: np.ndarray, score: np.ndarray, perms: np.ndarray) -> tuple[float, float, int]:
    """Spearman rho of y with score over the lines where both exist, and the two-sided permutation p from `perms` (rows = permutations of the line index)."""
    ok = ~np.isnan(y) & ~np.isnan(score)
    n = int(ok.sum())
    if n < 3:
        return float("nan"), float("nan"), n
    yr = stats.rankdata(y[ok])
    yr = (yr - yr.mean()) / yr.std()
    s_perm = score[perms][:, ok]  # permuted score on the lines where y exists
    valid = ~np.isnan(s_perm).any(axis=1)
    sr = stats.rankdata(s_perm[valid], axis=1)
    sr = (sr - sr.mean(axis=1, keepdims=True)) / sr.std(axis=1, keepdims=True)
    rho_perm = (sr * yr).mean(axis=1)
    s0 = stats.rankdata(score[ok])
    s0 = (s0 - s0.mean()) / s0.std()
    rho = float((s0 * yr).mean())
    p = (1 + int((np.abs(rho_perm) >= abs(rho) - 1e-12).sum())) / (1 + len(rho_perm))
    return rho, p, n


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    adj = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


def grade_from(rho: float, q: float) -> float:
    if q >= Q_CUT:
        return 0.5
    return 0.5 + 0.5 * max(-1.0, min(1.0, rho / RHO_SCALE))


def line_c_statistics(Y: dict[str, np.ndarray], score: np.ndarray, n_perm: int = N_PERM, seed: int = SEED, min_lines: int = MIN_LINES) -> dict[str, Any]:
    """Y: drug -> viability vector over the same lines as `score` (NaN allowed). Missing when fewer than `min_lines` lines have both values."""
    rng = np.random.default_rng(seed)
    perms = np.array([rng.permutation(len(score)) for _ in range(n_perm)])
    rows: dict[str, dict[str, Any]] = {}
    for d in sorted(Y):
        rho, p, n = spearman_perm(Y[d], score, perms)
        rows[d] = {"n_lines": n, "rho": None if np.isnan(rho) else round(rho, 4), "p_perm": None if np.isnan(p) else p}
    tested = [d for d in rows if rows[d]["n_lines"] >= min_lines and rows[d]["rho"] is not None]
    if tested:
        q = bh(np.array([rows[d]["p_perm"] for d in tested]))
        for d, qv in zip(tested, q):
            rows[d]["q_bh"] = float(qv)
            rows[d]["grade"] = round(grade_from(rows[d]["rho"], float(qv)), 4)
    for d in rows:
        rows[d].setdefault("grade", None)
        if d not in tested:
            rows[d]["note"] = "fewer than the minimum number of cell lines"
    return rows


def line_c_drug_level(root: Path, pool_names: list[str]) -> dict[str, Any]:
    from track2 import prism_l4

    data = root / "data" / "depmap"
    for f in ("OmicsExpressionProteinCodingGenesTPMLogp1.csv", "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv", "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"):
        if not (data / f).exists():
            raise FileNotFoundError(f"{f} is not on disk; nothing is downloaded")
    scores = prism_l4.sac_scores(data / "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    keys = prism_l4.compound_keys(data / "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")
    matches = prism_l4.match_rows({n: n for n in pool_names}, keys)
    with open(data / "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        col = {c: i + 1 for i, c in enumerate(header[1:]) if c in scores}
        lines = list(col)
        wanted = {r for rows in matches.values() for r in rows}
        mat = {row[0]: np.array([float(row[col[c]]) if row[col[c]] else np.nan for c in lines]) for row in reader if row[0] in wanted}
    score = np.array([scores[c]["mean_all"] for c in lines])
    prolif = np.array([scores[c]["prolif"] for c in lines])
    Y = {}
    for d, rows in matches.items():
        with np.errstate(all="ignore"):
            Y[d] = np.nanmean(np.vstack([mat[r] for r in rows]), axis=0)
    res = line_c_statistics(Y, score)
    for d in res:
        res[d]["n_prism_rows"] = len(matches[d])
    not_in = sorted(set(pool_names) - set(matches))
    n_sig = sum(1 for d in res if res[d].get("q_bh") is not None and res[d]["q_bh"] < Q_CUT)
    return {"label": "exploratory, correlational, from cancer cell lines; not validated; a hypothesis signal only",
            "definition": {"score": "mean over cell lines of the z-scored expression of BUB1B, BUB1, BUB3, MAD1L1, MAD2L1, TTK; low = weak checkpoint", "statistic": "Spearman rho between PRISM viability (log fold change) and the score; rho > 0 = more lethal when the checkpoint is weak",
                           "null": f"{N_PERM} permutations of the score across cell lines, seed {SEED}; BH across tested drugs", "grade": f"q < {Q_CUT}: 0.5 + 0.5 clip(rho / {RHO_SCALE}, -1, 1); else 0.5; missing if not in PRISM or fewer than {MIN_LINES} lines",
                           "not_adjusted_for": "proliferation and lineage"},
            "n_cell_lines": len(lines), "spearman_score_vs_proliferation": round(float(stats.spearmanr(score, prolif, nan_policy="omit").statistic), 3),
            "n_pool_drugs": len(pool_names), "n_in_prism": len(matches), "n_tested": sum(1 for d in res if res[d].get("q_bh") is not None), "n_significant_q_below_cut": n_sig,
            "not_in_prism": not_in, "drugs": res}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=("approval", "engine", "linec"))
    ap.add_argument("--out", type=Path, default=Path("out"))
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--pool-file", type=Path, default=None, help="linec: JSON list of pool drug names (default: pool of the completed inputs)")
    args = ap.parse_args()
    out = args.out
    if args.stage == "approval":
        from track2.fda_approval_check import check
        first = json.loads((out / "multi_evidence_mva.json").read_text(encoding="utf-8"))
        names = first["pool"]["approval_not_checked"]
        target = out / "fda_approval_check_pool_completion.json"
        prior = json.loads(target.read_text(encoding="utf-8")) if target.exists() else None
        res = complete_approval(names, check, prior)
        write_json(target, res)
        print({"asked": len(names), "confirmed": res["n_confirmed"], "not_confirmed": len(res["not_confirmed"])})
    elif args.stage == "engine":
        stored = json.loads((out / "engine" / "apply_mva.json").read_text(encoding="utf-8"))
        res = full_ranking_record(args.root, stored)
        write_json(out / "engine" / "apply_mva_full_ranking.json", res)
        print({"n_universe": res["n_universe"], "n_fda_small_molecule": res["n_fda_small_molecule"], "consistent_with_stored_top20": res["consistent_with_stored_top20"]})
    else:
        pool = json.loads(args.pool_file.read_text(encoding="utf-8"))
        res = line_c_drug_level(args.root, pool)
        write_json(out / "line_c_drug_level.json", res)
        print({k: res[k] for k in ("n_pool_drugs", "n_in_prism", "n_tested", "n_significant_q_below_cut", "spearman_score_vs_proliferation")})


if __name__ == "__main__":
    main()
