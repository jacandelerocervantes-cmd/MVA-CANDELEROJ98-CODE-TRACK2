r"""Engine v4 literature-dilution experiment (protocol/addendum_engine_v4_dilution.md, written before this module ran).

Run (PYTHONPATH=src), in this order:
    python -m track2.engine_v4_dilution features   thinned-paper features + drug texts of every benchmark disease (cache only, no labels)
    python -m track2.engine_v4_dilution embed      frozen BiomedBERT / SapBERT vectors of the NEW texts only (isolated `python -I` worker)
    python -m track2.engine_v4_dilution dev        frozen B0 / E on diluted DEV diseases + Europe PMC mapping of the project diseases
    python -m track2.engine_v4_dilution heldout    frozen B0 / E on diluted HELD-OUT diseases, once (lock file); completes the JSON

Reuses engine_v4 / engine_score functions unchanged; edits nothing. Label: post hoc development, not a pre-registered result.
"""
from __future__ import annotations

import json
import math
import pickle
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from track2 import engine_benchmark as v1
from track2 import engine_score as es
from track2 import engine_score_v2 as v2s
from track2 import engine_v3 as v3
from track2 import engine_v4 as v4

ROOT = v4.ROOT
SEED = 20261015
LABEL = "post hoc development (engine v4 literature dilution), not a pre-registered primary result"
OUT = v4.OUT
DIL_DIR = v4.CACHE / "v4_dilution"
FEATS_FILE = DIL_DIR / "dilution_features.pkl"
NEW_TEXTS = DIL_DIR / "new_texts.json"
PROTOCOL = ROOT / "protocol" / "addendum_engine_v4_dilution.md"
OUT_FILE = OUT / "validation_v4_dilution.json"
LOCK = OUT / "heldout_v4_dilution_scored.lock"
FRACTIONS = [1.0, 0.3, 0.1, 0.03, 0.01]
N_REPS = 5
FLOOR = 5
CAP = es.MAX_PAPERS
N_GRID = 200
FROZEN_FILE_SHA256 = "09e3a8e8f5a32ea79d455d69439a2edd86a96058fc85cd5a62be1c13b398e673"
MODELS = ("B0", "E")
PROJECT_SETS = {"Fanconi anaemia (project case genes ATR, FANCA)": ["ATR", "FANCA"],
                "MVA gene set (BUB1B, CEP57, TRIP13)": ["BUB1B", "CEP57", "TRIP13"]}

UNVALIDATED: Dict[str, Dict[str, Any]] = {
    "fractions": {"value": FRACTIONS, "status": "UNVALIDATED", "conservative_default": FRACTIONS, "note": "fixed by the task"},
    "floor_papers": {"value": FLOOR, "status": "UNVALIDATED", "conservative_default": FLOOR, "note": "fixed by the task"},
    "replicates": {"value": N_REPS, "status": "UNVALIDATED", "conservative_default": N_REPS, "note": "fixed by the task"},
    "gap_curve_form": {"value": "OLS gap ~ log10(n_papers)", "status": "UNVALIDATED", "conservative_default": "OLS gap ~ log10(n_papers)",
                       "note": "pooled summary, not a per-disease law"},
    "break_even_definition": {"value": "smallest grid N with bootstrap lower bound > 0 for all grid values >= N (grid 200 points, 5..300)",
                              "status": "UNVALIDATED", "conservative_default": "same", "note": "stated in protocol section 3"},
    "dev_weights_in_sample": {"value": "frozen weights fitted on DEV; DEV scores optimistic", "status": "UNVALIDATED",
                              "conservative_default": "HELD-OUT is the primary reading", "note": "protocol section 1"},
}


# ============================================================================= design
def n_sub(n: int, f: float) -> int:
    """Subsample size: min(n, max(floor, ceil(f * n))); f = 1 keeps everything."""
    if f >= 1.0:
        return int(n)
    return int(min(n, max(FLOOR, math.ceil(f * n - 1e-9))))


def subsample_indices(n: int, f: float, rank: int, frac_index: int, rep: int, seed: int = SEED) -> np.ndarray:
    """Sorted indices (original retrieval order kept) of a without-replacement subsample; deterministic in (seed, rank, fraction, rep)."""
    k = n_sub(n, f)
    if k >= n:
        return np.arange(n)
    rng = np.random.default_rng([seed, rank, frac_index, rep])
    return np.sort(rng.choice(n, size=k, replace=False))


def sample_keys() -> List[Tuple[int, int]]:
    """(frac_index, rep) of every sample: f = 1.0 has one deterministic replicate, the others N_REPS."""
    return [(0, 0)] + [(fi, r) for fi in range(1, len(FRACTIONS)) for r in range(N_REPS)]


def key_of(fi: int, rep: int) -> Tuple[int, int]:
    return (0, 0) if fi == 0 else (fi, rep)


def cache_only_getter(url: str, params: dict) -> dict:
    raise RuntimeError("cache miss: this stage never calls Europe PMC")


# ============================================================================= features
def thinned_features(symbols: Sequence[str], papers: Sequence[dict], rank: int, W: Dict[str, Any]) -> Dict[Tuple[int, int], Dict[str, Any]]:
    """Features of ONE disease for every (fraction, replicate): (n_drugs, 8) matrix + drug texts. No gold argument exists."""
    names, dre, gate, pop = W["names"], W["dre"], W["gate"], W["pop"]
    prior = es.prior_features(symbols, names, W["targets"], W["ppi"])
    out: Dict[Tuple[int, int], Dict[str, Any]] = {}
    for fi, rep in sample_keys():
        idx = subsample_indices(len(papers), FRACTIONS[fi], rank, fi, rep)
        sub = [papers[i] for i in idx]
        lit = es.literature_counts(sub, symbols, dre, gate)
        X = es.disease_matrix(lit, names, pop, W["total"], prior)
        out[(fi, rep)] = {"X": X, "texts": v4.drug_sentences(sub, symbols, dre, gate), "n_sub": len(sub),
                          "n_papers_with_gene": lit["n_papers_with_gene"]}
    return out


def _job(args: Tuple[str, Sequence[str], int]) -> Tuple[str, Dict[str, Any]]:
    name, symbols, rank = args
    papers = es.fetch_papers(symbols, es.YEAR_PRIMARY, v4.CACHE / "epmc", cache_only_getter)
    return name, {"samples": thinned_features(symbols, papers, rank, v4._W), "n_retrieved": len(papers)}


def run_features(root: Path = ROOT, workers: int = 6) -> dict:
    bench = v4.load_benchmark(root)
    drugs = v1.load_universe(root / "data" / "universe.json")
    names = sorted(v1.universe_index(drugs))
    ranked = sorted(d["name"] for d in bench["diseases"])
    rank = {n: i for i, n in enumerate(ranked)}
    DIL_DIR.mkdir(parents=True, exist_ok=True)
    done: Dict[str, Any] = {}
    if FEATS_FILE.exists():
        with FEATS_FILE.open("rb") as fh:  # project-own cache written by this module
            done = pickle.load(fh)
    jobs = []
    for dis in bench["diseases"]:
        if dis["name"] not in done:
            jobs.append((dis["name"], sorted({g["symbol"] for g in dis["genes"]}), rank[dis["name"]]))
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers, initializer=v4._init_worker, initargs=(names,)) as ex:
        for i, (name, res) in enumerate(ex.map(_job, jobs, chunksize=1)):
            done[name] = res
            if (i + 1) % 10 == 0:
                print(json.dumps({"features": i + 1, "of": len(jobs), "seconds": round(time.time() - t0)}), flush=True)
    with FEATS_FILE.open("wb") as fh:
        pickle.dump(done, fh)
    print(json.dumps({"features_done": len(done), "new": len(jobs)}), flush=True)
    return done


# ============================================================================= embeddings of the new texts
def run_embed(root: Path = ROOT) -> dict:
    with FEATS_FILE.open("rb") as fh:
        feats = pickle.load(fh)
    old = set(json.loads((v4.EMB_DIR / "texts.json").read_text(encoding="utf-8")))
    new: Dict[str, int] = {}
    for d in sorted(feats):
        for s in feats[d]["samples"].values():
            for t in s["texts"].values():
                if t not in old:
                    new.setdefault(t, len(new))
    lst = [t for t, _ in sorted(new.items(), key=lambda kv: kv[1])]
    NEW_TEXTS.write_text(json.dumps(lst), encoding="utf-8")
    stats: Dict[str, Any] = {"n_new_texts": len(lst), "n_old_texts": len(old)}
    worker = ROOT / "src" / "track2" / "embed_worker.py"
    for key, spec in v4.MODEL_SPECS.items():
        dst = DIL_DIR / f"{key}_new.npy"
        if dst.exists():
            stats[key] = {"cached": True}
            continue
        cmd = [sys.executable, "-I", str(worker), str(v4.MODELS_DIR / key), spec["pool"], str(NEW_TEXTS), str(dst), str(v4.MAX_TOK), "32"]
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(DIL_DIR))
        if r.returncode != 0:
            raise RuntimeError(f"embedding worker failed for {key}: {r.stderr[-800:]}")
        stats[key] = json.loads(r.stdout.strip().splitlines()[-1])
        print(json.dumps({key: stats[key]}), flush=True)
    v4.write_json(DIL_DIR / "embed_stats_new.json", stats)
    return stats


def load_all_embeddings() -> Tuple[Dict[str, int], Dict[str, np.ndarray]]:
    """v4 embeddings plus the embeddings of the new diluted texts, as one index."""
    tindex, emb = v4.load_embeddings()
    new = json.loads(NEW_TEXTS.read_text(encoding="utf-8"))
    merged = dict(tindex)
    for j, t in enumerate(new):
        merged[t] = len(tindex) + j
    out = {k: np.vstack([emb[k], np.load(DIL_DIR / f"{k}_new.npy")]) for k in emb}
    return merged, out


# ============================================================================= scoring (frozen heads)
def sample_Z(name: str, symbols: Sequence[str], sample: Dict[str, Any], names: Sequence[str], tindex, emb) -> np.ndarray:
    f = {"symbols": symbols, "drug_texts": sample["texts"]}
    return np.column_stack([sample["X"], v4.emb_block(name, f, names, tindex, emb)])


def frozen_heads(frozen_path: Path) -> Dict[str, Any]:
    fr = json.loads(frozen_path.read_text(encoding="utf-8"))
    m = v4.thaw_v4(fr)
    return {"model": m, "file_sha256": v4.sha256_file(frozen_path), "sha256_of_model": fr["sha256_of_model"]}


def metrics_one(score: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    return es.auc_mw(score, y), es.precision_at_k(score, y, 10), es.recall_at_k(score, y, 20)


def score_diseases(diseases: Sequence[str], feats: Dict[str, Any], heads: Dict[str, Any], y: Dict[str, np.ndarray],
                   pop: Dict[str, np.ndarray], symbols: Dict[str, Sequence[str]], names: Sequence[str], tindex, emb) -> Dict[str, Any]:
    """Per disease, per (fraction, replicate): AUC / P@10 / R@20 of frozen B0 and E, plus popularity AUC/P@10/R@20 (replicate-free)
    and n_sub. Returns arrays indexed [disease][model][fraction] -> list over replicates."""
    res: Dict[str, Any] = {}
    for d in diseases:
        r: Dict[str, Any] = {"n_sub": {}, "metrics": {m: {} for m in MODELS}}
        for fi, rep in sample_keys():
            s = feats[d]["samples"][(fi, rep)]
            Z = sample_Z(d, symbols[d], s, names, tindex, emb)
            r["n_sub"].setdefault(fi, []).append(s["n_sub"])
            for m in MODELS:
                h = heads["heads"][m]
                sc = v4.head_decision(h, Z[:, h["cols"]])
                r["metrics"][m].setdefault(fi, []).append(metrics_one(sc, y[d]))
        r["popularity"] = metrics_one(pop[d], y[d])
        res[d] = r
    return res


# ============================================================================= aggregation
def per_fraction_arrays(res: Dict[str, Any], dn: Sequence[str], model: str, fi: int) -> Dict[str, np.ndarray]:
    """Per-disease mean over replicates for fraction fi: auc, p10, r20 (NaN-safe)."""
    a = np.array([[np.nanmean([x[k] for x in res[d]["metrics"][model][fi]]) for k in range(3)] for d in dn])
    return {"auc": a[:, 0], "p10": a[:, 1], "r20": a[:, 2]}


def boot_idx(n_units: int, seed: int = SEED, n_boot: int = es.N_BOOT) -> np.ndarray:
    return np.random.default_rng(seed).integers(0, n_units, size=(n_boot, n_units))


def table_for_set(res: Dict[str, Any], members: Dict[str, List[str]], y: Dict[str, np.ndarray], disease_level: bool = False) -> Dict[str, Any]:
    """Per-fraction table: macro AUC / gap vs popularity (CI), P@10, R@20, mean n_sub, replicate SD of the macro gap."""
    if disease_level:
        members = {d: [d] for g in members for d in members[g]}
    gs = sorted(members)
    dn = sorted(d for g in gs for d in members[g])
    idx = boot_idx(len(gs))
    pop = {k: np.array([res[d]["popularity"][i] for d in dn]) for i, k in enumerate(("auc", "p10", "r20"))}
    pop_g = {k: v4.agg_groups(v, dn, members) for k, v in pop.items()}
    rnd = es.random_per_disease(y, dn)
    rnd_g = {k: v4.agg_groups(v, dn, members) for k, v in rnd.items()}
    out: Dict[str, Any] = {"n_diseases": len(dn), "n_groups": len(gs), "random": {k: es.macro_with_ci(v, idx) for k, v in rnd_g.items()},
                           "popularity": {k: es.macro_with_ci(v, idx) for k, v in pop_g.items()}, "fractions": {}}
    for fi, f in enumerate(FRACTIONS):
        row: Dict[str, Any] = {"fraction": f}
        ns = [x for d in dn for x in res[d]["n_sub"][fi]]
        row["n_papers_mean"], row["n_papers_median"] = float(np.mean(ns)), float(np.median(ns))
        row["n_papers_min"], row["n_papers_max"] = int(np.min(ns)), int(np.max(ns))
        for m in MODELS:
            arr = per_fraction_arrays(res, dn, m, fi)
            g = {k: v4.agg_groups(v, dn, members) for k, v in arr.items()}
            row[m] = {"auc": es.macro_with_ci(g["auc"], idx), "auc_gap_vs_popularity": es.diff_with_ci(g["auc"], pop_g["auc"], idx),
                      "precision_at_10": es.macro_with_ci(g["p10"], idx), "recall_at_20": es.macro_with_ci(g["r20"], idx),
                      "precision_at_10_minus_popularity": es.diff_with_ci(g["p10"], pop_g["p10"], idx)}
            nrep = len(res[dn[0]]["metrics"][m][fi])
            reps = []
            for r in range(nrep):
                a = np.array([res[d]["metrics"][m][fi][r][0] for d in dn])
                reps.append(float(np.nanmean(v4.agg_groups(a, dn, members) - pop_g["auc"])))
            row[m]["replicate_gap_sd"] = float(np.std(reps, ddof=1)) if len(reps) > 1 else 0.0
            row[m]["n_replicates"] = nrep
        out["fractions"][str(f)] = row
    return out


# ============================================================================= gap versus log10(papers) and break-even
def fit_gap_curve(res: Dict[str, Any], members: Dict[str, List[str]], model: str, seed: int = SEED, n_boot: int = es.N_BOOT) -> Dict[str, Any]:
    """OLS gap = a + b * log10(n_sub) over all (disease, fraction, replicate) observations with a cluster bootstrap over disease
    groups; the pointwise lower bound on a log grid gives the break-even literature volume."""
    gs = sorted(members)
    xs, ys, cl = [], [], []
    for gi, g in enumerate(gs):
        for d in members[g]:
            pa = res[d]["popularity"][0]
            for fi in range(len(FRACTIONS)):
                for k, n in enumerate(res[d]["n_sub"][fi]):
                    a = res[d]["metrics"][model][fi][k][0]
                    if not (math.isnan(a) or math.isnan(pa)):
                        xs.append(math.log10(n))
                        ys.append(a - pa)
                        cl.append(gi)
    x, yv, cl = np.array(xs), np.array(ys), np.array(cl)
    G = len(gs)
    stats = np.zeros((G, 5))
    for gi in range(G):
        m = cl == gi
        stats[gi] = [m.sum(), x[m].sum(), yv[m].sum(), (x[m] ** 2).sum(), (x[m] * yv[m]).sum()]

    def ols(t: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        n, sx, sy, sxx, sxy = t.T if t.ndim == 2 else t
        den = n * sxx - sx ** 2
        with np.errstate(all="ignore"):
            b = np.where(np.abs(den) > 1e-12, (n * sxy - sx * sy) / den, np.nan)
            a = (sy - b * sx) / n
        return a, b

    a0, b0 = ols(stats.sum(0))
    a0, b0 = float(a0), float(b0)
    idx = boot_idx(G, seed, n_boot)
    counts = np.stack([np.bincount(r, minlength=G) for r in idx]).astype(float)
    ab, bb = ols(counts @ stats)
    grid_n = np.logspace(math.log10(FLOOR), math.log10(CAP), N_GRID)
    gx = np.log10(grid_n)
    fit = a0 + b0 * gx
    boots = ab[:, None] + bb[:, None] * gx[None, :]
    ok = ~np.isnan(boots).any(1)
    lo = np.percentile(boots[ok], 2.5, axis=0)
    hi = np.percentile(boots[ok], 97.5, axis=0)
    return {"model": model, "n_observations": int(len(x)), "intercept": a0, "slope_per_log10_papers": b0,
            "slope_ci95": [float(np.percentile(bb[ok], 2.5)), float(np.percentile(bb[ok], 97.5))],
            "grid_n_papers": grid_n.tolist(), "fit": fit.tolist(), "lower95": lo.tolist(), "upper95": hi.tolist(),
            "break_even": break_even(grid_n, lo, fit)}


def break_even(grid_n: np.ndarray, lower: np.ndarray, fit: np.ndarray) -> Dict[str, Any]:
    """Break-even literature volume: smallest grid N whose lower bound is > 0 for every grid value >= N (protocol section 3)."""
    pos = np.asarray(lower) > 0
    point = np.asarray(fit) > 0
    suffix = np.flip(np.cumprod(np.flip(pos).astype(int))).astype(bool)      # True where lower > 0 from here to the end
    res: Dict[str, Any] = {}
    if suffix.all():
        res.update(value=None, statement="lower bound above 0 over the whole grid: break-even at or below the smallest simulated volume "
                                          f"({FLOOR} papers)")
    elif not suffix.any():
        res.update(value=None, statement=f"not reached: the lower bound is not above 0 at {CAP} papers")
    else:
        i = int(np.argmax(suffix))
        res.update(value=float(grid_n[i]), statement="lower bound above 0 from this volume upward")
    sp = np.flip(np.cumprod(np.flip(point).astype(int))).astype(bool)
    res["point_estimate_crossing"] = (float(grid_n[int(np.argmax(sp))]) if sp.any() and not sp.all() else
                                      ("fitted gap above 0 over the whole grid" if sp.all() else "fitted gap not above 0 at the cap"))
    return res


def placement(n_papers: int, be: Dict[str, Any]) -> str:
    v = be.get("value")
    if v is None:
        return "break-even below the simulated range" if be["statement"].startswith("lower bound above 0 over") else "break-even not reached"
    return ("above" if n_papers >= v else "below") + f" the break-even volume ({v:.1f} papers)" + \
           (" (at the retrieval cap of 300)" if n_papers >= CAP else "")


# ============================================================================= mapping of the project diseases
def map_projects(feats_v4: Dict[str, Any], getter: Optional[Callable[[str, dict], dict]] = None) -> Dict[str, Any]:
    """Papers retrieved (cap 300) and hit count for the project gene sets (Europe PMC only if not cached) and for the 6 benchmark
    diseases with the fewest retrieved papers (from the v4 feature cache)."""
    getter = getter or v2s.epmc_getter
    hits = json.loads(v4.HITS_FILE.read_text(encoding="utf-8")) if v4.HITS_FILE.exists() else {}
    out: Dict[str, Any] = {"project": {}, "fewest_benchmark_diseases": []}
    cnt = v4.ReqCounter(getter)
    for label, syms in PROJECT_SETS.items():
        papers = es.fetch_papers(sorted(syms), es.YEAR_PRIMARY, v4.CACHE / "epmc", cnt)
        h = v4.hit_count(sorted(syms), cnt, hits)
        out["project"][label] = {"genes": sorted(syms), "n_retrieved": len(papers), "hit_count": h, "at_retrieval_cap": len(papers) >= CAP}
    v4.HITS_FILE.write_text(json.dumps(hits), encoding="utf-8")
    out["europe_pmc_requests_made"] = cnt.n
    rows = sorted((v["n_retrieved"], k) for k, v in feats_v4.items())[:6]
    for n, k in rows:
        out["fewest_benchmark_diseases"].append({"disease": k, "n_retrieved": int(n), "hit_count": int(feats_v4[k]["hit_count"]),
                                                 "genes": feats_v4[k]["symbols"]})
    return out


def place_mapping(mapping: Dict[str, Any], be_by_set: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"project": {}, "fewest_benchmark_diseases": []}
    for label, v in mapping["project"].items():
        out["project"][label] = dict(v, relative_to_break_even={s: placement(v["n_retrieved"], be) for s, be in be_by_set.items()})
    for v in mapping["fewest_benchmark_diseases"]:
        out["fewest_benchmark_diseases"].append(dict(v, relative_to_break_even={s: placement(v["n_retrieved"], be) for s, be in be_by_set.items()}))
    out["europe_pmc_requests_made"] = mapping.get("europe_pmc_requests_made")
    return out


# ============================================================================= commands
def _context(root: Path, subset: Sequence[str], read_labels: bool) -> Dict[str, Any]:
    bench = v4.load_benchmark(root)
    drugs = v1.load_universe(root / "data" / "universe.json")
    universe = v1.universe_index(drugs)
    names = sorted(universe)
    index = {n: i for i, n in enumerate(names)}
    with v4.FEATURES_FILE.open("rb") as fh:
        f4 = pickle.load(fh)
    with FEATS_FILE.open("rb") as fh:
        fd = pickle.load(fh)
    tindex, emb = load_all_embeddings()
    by_name = {d["name"]: d for d in bench["diseases"]}
    y: Dict[str, np.ndarray] = {}
    if read_labels:
        for d in subset:
            yy = np.zeros(len(names))
            for g in by_name[d]["gold1"]:
                if g in index:
                    yy[index[g]] = 1
            y[d] = yy
    return {"bench": bench, "names": names, "f4": f4, "fd": fd, "tindex": tindex, "emb": emb, "y": y,
            "pop": {d: f4[d]["popularity"] for d in subset}, "symbols": {d: f4[d]["symbols"] for d in subset}}


def _verify(root: Path, bench: dict) -> Dict[str, Any]:
    fz = frozen_heads(OUT / "engine_v4_frozen.json")
    m = fz["model"]
    split = bench["split"]
    chk = {"frozen_file_sha256": fz["file_sha256"], "frozen_file_sha256_expected": FROZEN_FILE_SHA256,
           "frozen_model_hash_ok": True, "dev_ids_hash_ok": v4.hash_ids(split["dev"]) == m["sha256_dev_ids"],
           "heldout_ids_hash_ok": v4.hash_ids(split["heldout"]) == m["sha256_heldout_ids"],
           "frozen_file_sha256_ok": fz["file_sha256"] == FROZEN_FILE_SHA256}
    if not (chk["dev_ids_hash_ok"] and chk["heldout_ids_hash_ok"] and chk["frozen_file_sha256_ok"]):
        raise RuntimeError("frozen file or split hashes differ from the frozen record: " + json.dumps(chk))
    fz["check"] = chk
    return fz


def _f1_check(ctx: Dict[str, Any], diseases: Sequence[str]) -> float:
    """f = 1.0 must reproduce the v4 feature matrix and the v4 drug texts exactly (max abs difference)."""
    worst = 0.0
    for d in diseases:
        s = ctx["fd"][d]["samples"][(0, 0)]
        worst = max(worst, float(np.abs(s["X"] - ctx["f4"][d]["X"]).max()))
        if s["texts"] != ctx["f4"][d]["drug_texts"]:
            worst = max(worst, 1.0)
    return worst


def _block(root: Path, which: str) -> Dict[str, Any]:
    bench = v4.load_benchmark(root)
    subset = bench["split"][which]
    ctx = _context(root, subset, read_labels=True)
    fz = _verify(root, bench)
    members = v4.members_of(bench, subset)
    t0 = time.time()
    res = score_diseases(subset, ctx["fd"], fz["model"], ctx["y"], ctx["pop"], ctx["symbols"], ctx["names"], ctx["tindex"], ctx["emb"])
    block: Dict[str, Any] = {"n_diseases": len(subset), "n_groups": len(members), "diseases": sorted(subset),
                             "f1_reproduces_v4_features_max_abs_diff": _f1_check(ctx, subset),
                             "table_group_bootstrap": table_for_set(res, members, ctx["y"]),
                             "table_disease_level_bootstrap_sensitivity": table_for_set(res, members, ctx["y"], disease_level=True),
                             "gap_curve": {m: fit_gap_curve(res, members, m) for m in MODELS},
                             "n_retrieved_by_disease": {d: ctx["fd"][d]["n_retrieved"] for d in subset},
                             "scoring_seconds": round(time.time() - t0, 1), "frozen_check": fz["check"]}
    return block


def run_dev(root: Path = ROOT) -> dict:
    out = {"label": LABEL, "protocol": str(PROTOCOL.relative_to(ROOT)), "protocol_sha256_at_dev_run": v4.sha256_file(PROTOCOL),
           "settings": {"fractions": FRACTIONS, "replicates": N_REPS, "floor_papers": FLOOR, "seed": SEED, "n_boot": es.N_BOOT,
                        "grid_points": N_GRID, "retrieval_cap": CAP, "models": list(MODELS)},
           "unvalidated_registry": UNVALIDATED,
           "status": "DEV part written; HELD-OUT not scored yet",
           "dev_note": "frozen weights were fitted on DEV: DEV scores are in-sample for the weights (optimistic); HELD-OUT is the primary reading"}
    out["dev"] = _block(root, "dev")
    with v4.FEATURES_FILE.open("rb") as fh:
        f4 = pickle.load(fh)
    try:
        out["mapping_raw"] = map_projects(f4)
    except Exception as e:  # noqa: BLE001 - reported, never hidden
        out["mapping_raw"] = {"error": repr(e)[:300], "project": {}, "fewest_benchmark_diseases": []}
    v4.write_json(OUT_FILE, out)
    be = {m: out["dev"]["gap_curve"][m]["break_even"] for m in MODELS}
    print(json.dumps({"dev_done": True, "break_even": be}, sort_keys=True), flush=True)
    return out


def run_heldout(root: Path = ROOT) -> dict:
    if LOCK.exists():
        raise RuntimeError("HELD-OUT dilution already scored once (lock file exists); refusing a second scoring")
    out = json.loads(OUT_FILE.read_text(encoding="utf-8"))
    if "dev" not in out:
        raise RuntimeError("run the dev command first")
    bench = v4.load_benchmark(root)
    fz = _verify(root, bench)
    LOCK.write_text(json.dumps({"scored_at": "first and only held-out dilution scoring", "frozen_file_sha256": fz["file_sha256"],
                                "protocol_sha256": v4.sha256_file(PROTOCOL)}) + "\n", encoding="utf-8")
    out["heldout"] = _block(root, "heldout")
    be_ho = {m: out["heldout"]["gap_curve"][m]["break_even"] for m in MODELS}
    be_dev = {m: out["dev"]["gap_curve"][m]["break_even"] for m in MODELS}
    out["break_even_literature_volume"] = {"primary_B0_heldout": be_ho["B0"], "secondary_E_heldout": be_ho["E"],
                                           "sensitivity_B0_dev_in_sample_weights": be_dev["B0"], "sensitivity_E_dev_in_sample_weights": be_dev["E"]}
    be_sets = {"B0_heldout": be_ho["B0"], "B0_dev": be_dev["B0"], "E_heldout": be_ho["E"]}
    out["mapping"] = place_mapping(out.pop("mapping_raw"), be_sets)
    out["protocol_sha256_at_heldout_run"] = v4.sha256_file(PROTOCOL)
    out["status"] = "complete: DEV and HELD-OUT; HELD-OUT scored once for frozen B0 and E"
    v4.write_json(OUT_FILE, out)
    print(json.dumps({"heldout_done": True, "break_even": be_ho}, sort_keys=True, default=v3._default), flush=True)
    return out


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = {"features": run_features, "embed": run_embed, "dev": run_dev, "heldout": run_heldout}.get(cmd)
    if fn is None:
        raise SystemExit("usage: python -m track2.engine_v4_dilution features|embed|dev|heldout")
    fn()


if __name__ == "__main__":
    main()
