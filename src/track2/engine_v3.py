r"""Evidence engine v3 - post hoc development with a sealed held-out set (protocol/addendum_engine_v3.md).

Run (PYTHONPATH=src):
    python -m track2.engine_v3 dev        development on DEV only, freeze, writes out/engine/engine_v3_frozen.json,
                                          out/engine/engine_v1_weights_frozen.json and out/engine/validation_v3.json (heldout = null)
    python -m track2.engine_v3 heldout    scores the 14 HELD-OUT diseases exactly once with the frozen files (lock file)

Reuses engine_score (v1) functions and engine_score_v2.build_features unchanged. Label: post hoc development, not a
pre-registered primary result. No network is used: every retrieval comes from data/cache_engine.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
import re
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression

from track2 import engine_benchmark as eb
from track2 import engine_score as es
from track2 import engine_score_v2 as v2

ROOT = Path(__file__).resolve().parents[2]
SEED = es.SEED
LABEL = "post hoc development (engine v3), not a pre-registered primary result"
C_GRID = (0.01, 0.03, 0.1, 0.3, 1.0)
ALPHA_F3 = 20.0
SEED_INNER, SEED_FINAL = SEED + 1, SEED + 2
WILMS = "kidney Wilms tumor"
WILMS_DRUGS = ("doxorubicin", "dactinomycin", "vincristine", "tamoxifen", "progesterone", "danazol")

F1_NAMES = ["f1_tfidf_papers", "f1_tfidf_benefit", "f1_tfidf_human", "f1_pmi_benefit"]
F2_NAMES = ["f2_class_gene_papers", "f2_class_gene_pmi", "f2_class_keyword_papers"]
F3_NAME = "f3_atc_class_rate"
ALL_NAMES = list(es.FEATURE_NAMES) + F1_NAMES + F2_NAMES          # 8 + 4 + 3 = 15 stored columns; F3 is added in the fit
FAMILY_COLS = {"F1": list(range(8, 12)), "F2": list(range(12, 15))}
CONFIGS: List[Tuple[str, ...]] = [(), ("F1",), ("F2",), ("F3",), ("F1", "F2"), ("F1", "F3"), ("F2", "F3"), ("F1", "F2", "F3")]

CLASS_KEYWORDS = {
    "lysosomal_vesicular_trafficking": ["lysosom", "vesicl", "endosom", "autophag", "traffick", "storage", "exocyt", "endocyt"],
    "ion_channel_transporter": ["channel", "transporter", "sodium", "potassium", "calcium", "chloride", "electrolyte"],
    "metabolic_enzyme": ["enzyme", "metabol", "substrate", "synthase", "dehydrogenase", "hydroxylase", "cofactor"],
    "receptor_signalling": ["receptor", "signalling", "signaling", "kinase", "agonist", "antagonist", "inhibitor", "pathway"],
    "transcription_chromatin_regulation": ["transcription", "chromatin", "histone", "methylat", "epigenet", "gene expression"],
    "structural_cytoskeletal": ["collagen", "cytoskelet", "actin", "myosin", "muscle", "fibrosis", "cardiomyopath"],
    "mitochondrial": ["mitochond", "oxidative phosphorylation", "respiratory chain", "oxidative stress"],
    "dna_repair_replication": ["DNA repair", "DNA damage", "replication", "genotoxic", "crosslink", "PARP"],
    "cell_cycle_mitotic_checkpoint": ["cell cycle", "mitosis", "mitotic", "checkpoint", "proliferat", "apoptosis"],
    "rna_processing_translation": ["RNA", "splic", "translation", "ribosom"],
    "extracellular_secreted": ["secreted", "extracellular", "matrix", "cytokine", "complement", "growth factor", "hormone"],
}


def _no_network(url: str, params: dict) -> dict:
    raise RuntimeError("engine_v3 uses caches only; a network request was attempted: " + url)


# ============================================================================= split
def make_split(diseases: Sequence[dict], seed: int = SEED) -> Dict[str, List[str]]:
    """Protocol section 1: DEV = all v1 + Wilms + floor(n/2) of the other v2_new (sorted, seeded shuffle); HELD-OUT = the rest."""
    v1 = sorted(d["name"] for d in diseases if d["origin"] == "v1")
    new = sorted(d["name"] for d in diseases if d["origin"] == "v2_new" and d["name"] != WILMS)
    order = list(new)
    random.Random(seed).shuffle(order)
    k = len(order) // 2
    dev = sorted(v1 + [WILMS] + order[:k])
    held = sorted(order[k:])
    return {"dev": dev, "heldout": held}


# ============================================================================= features
def cfg_name(cfg: Sequence[str]) -> str:
    return "base" if not cfg else "+".join(cfg)


def cfg_cols(cfg: Sequence[str]) -> List[int]:
    cols = list(range(8))
    for f in cfg:
        cols += FAMILY_COLS.get(f, [])
    return cols


def cfg_features(cfg: Sequence[str]) -> List[str]:
    return [ALL_NAMES[c] for c in cfg_cols(cfg)] + ([F3_NAME] if "F3" in cfg else [])


def complexity(cfg: Sequence[str]) -> Tuple[int, int]:
    return len(cfg), len(cfg_features(cfg))


def keyword_regex(cls: Optional[str]) -> Optional["re.Pattern[str]"]:
    kws = CLASS_KEYWORDS.get(cls or "")
    if not kws:
        return None
    return re.compile(r"\b(?:" + "|".join(re.escape(k) for k in kws) + ")", re.I)


def disease_class_genes(genes: Sequence[dict]) -> Tuple[Optional[str], List[str]]:
    """Plurality mechanism_class_general ignoring 'other' (ties alphabetical) and the symbols of that class."""
    c = Counter(g.get("mechanism_class_general") or "other" for g in genes if (g.get("mechanism_class_general") or "other") != "other")
    if not c:
        return None, []
    top = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    return top, sorted({g["symbol"] for g in genes if g.get("mechanism_class_general") == top})


def class_counts(papers: Sequence[dict], all_symbols: Sequence[str], class_symbols: Sequence[str],
                 kw_rx: Optional["re.Pattern[str]"], dre: "re.Pattern[str]", gate: Optional[set]) -> Dict[str, Any]:
    """Per drug: papers with the drug and a class gene (cg) and papers with a sentence holding the drug and a class keyword (kw)."""
    from track2 import empirical_literature_control as elc
    g_all = es.gene_regex(all_symbols)
    g_cls = es.gene_regex(class_symbols) if class_symbols else None
    cg: Dict[str, set] = {}
    kw: Dict[str, set] = {}
    n_cg = 0
    for p in papers:
        text = p["text"]
        if not g_all.search(text):
            continue
        has_cg = bool(g_cls and g_cls.search(text))
        n_cg += has_cg
        if not has_cg and kw_rx is None:
            continue
        for sent in elc.split_sentences(text):
            if gate is not None and not (set(re.findall(r"\w+", sent.lower())) & gate):
                continue
            has_kw = bool(kw_rx and kw_rx.search(sent))
            if not has_cg and not has_kw:
                continue
            for name in {m.group(0).lower() for m in dre.finditer(sent)}:
                if has_cg:
                    cg.setdefault(name, set()).add(p["pid"])
                if has_kw:
                    kw.setdefault(name, set()).add(p["pid"])
    return {"cg": {k: len(v) for k, v in cg.items()}, "kw": {k: len(v) for k, v in kw.items()}, "n_cg": n_cg}


def f1_block(c: np.ndarray, nb: np.ndarray, nh: np.ndarray, pop: np.ndarray, n_gene: float, total: float) -> np.ndarray:
    idf = np.log(total / (pop + 1.0))
    expected = n_gene * pop / max(total, 1.0)
    rho = float(nb.sum() / c.sum()) if c.sum() > 0 else 0.0
    return np.column_stack([np.log1p(c) * idf, np.log1p(nb) * idf, np.log1p(nh) * idf,
                            np.log((nb + 0.5) / (expected * rho + 0.5))])


def f2_block(names: Sequence[str], cc: Dict[str, Any], pop: np.ndarray, total: float, has_class: bool) -> np.ndarray:
    n = len(names)
    if not has_class:
        return np.zeros((n, 3))
    cg = np.array([cc["cg"].get(x, 0) for x in names], dtype=float)
    kw = np.array([cc["kw"].get(x, 0) for x in names], dtype=float)
    e_cg = cc["n_cg"] * pop / max(total, 1.0)
    return np.column_stack([np.log1p(cg), np.log((cg + 0.5) / (e_cg + 0.5)), np.log1p(kw)])


def build_Z(root: Path, bench2: dict, v1_bench: dict, subset: Sequence[str]) -> Dict[str, Any]:
    """Features of the given diseases from caches only: base (v1 pipeline of engine_score_v2), F1, F2."""
    from track2 import empirical_literature_control as elc
    sub = {"diseases": [d for d in bench2["diseases"] if d["name"] in set(subset)]}
    feats = v2.build_features(root, sub, v1_bench, _no_network)
    names, classes, per = feats["names"], feats["classes"], feats["per"]
    drugs = eb.load_universe(root / "data" / "universe.json")
    universe = eb.universe_index(drugs)
    dre = elc.drug_regex(names)
    gate = es.name_gate(names)
    _, total = es.fetch_popularity(names, None, root / "data" / "cache_engine" / "epmc_popularity_all.json", _no_network)
    total = float(total)
    cache = root / "data" / "cache_engine" / "epmc"
    Z: Dict[str, np.ndarray] = {}
    y: Dict[str, np.ndarray] = {}
    pop: Dict[str, np.ndarray] = {}
    proxy: Dict[str, np.ndarray] = {}
    info: Dict[str, Any] = {}
    for dis in sub["diseases"]:
        d = dis["name"]
        P = per[d]["primary"]
        X = P["X"]
        c = np.rint(np.expm1(X[:, 0]))
        nb = np.rint(np.expm1(X[:, 1]))
        nh = np.rint(np.expm1(X[:, 2]))
        p = P["popularity"]
        f1 = f1_block(c, nb, nh, p, P["n_papers_with_gene"], total)
        cls, csyms = disease_class_genes(dis["genes"])
        symbols = per[d]["symbols"]
        papers = es.fetch_papers(symbols, es.YEAR_PRIMARY, cache, _no_network)
        cc = class_counts(papers, symbols, csyms, keyword_regex(cls), dre, gate)
        f2 = f2_block(names, cc, p, total, cls is not None)
        Z[d] = np.column_stack([X[:, :8], f1, f2])
        y[d], pop[d], proxy[d] = per[d]["gold1"], p, P["proxy"]
        info[d] = {"mechanism_class": cls, "n_class_genes": len(csyms), "n_class_gene_papers": cc["n_cg"],
                   "n_gold1": int(per[d]["gold1"].sum())}
    cls_names, cls_id = np.unique(classes, return_inverse=True)
    return {"Z": Z, "y": y, "pop": pop, "proxy": proxy, "names": names, "classes": classes, "cls_id": cls_id,
            "n_classes": len(cls_names), "info": info, "total": total}


# ============================================================================= model
class F3Ref:
    """ATC class rate from reference diseases; leave-one-disease-out tables for training rows."""

    def __init__(self, ds: Sequence[str], y: Dict[str, np.ndarray], cls_id: np.ndarray, k: int):
        self.g = {d: np.bincount(cls_id, weights=y[d], minlength=k) for d in ds}
        self.m = np.bincount(cls_id, minlength=k).astype(float)
        self.G = sum(self.g.values()) if self.g else np.zeros(k)
        self.n = len(self.g)

    def table(self, exclude: Optional[str] = None) -> np.ndarray:
        G = self.G - (self.g[exclude] if exclude is not None else 0.0)
        n = self.n - (1 if exclude is not None else 0)
        M = self.m * n
        if n <= 0 or M.sum() <= 0 or G.sum() <= 0:
            return np.zeros(len(self.m))
        p0 = G.sum() / M.sum()
        return np.log(((G + ALPHA_F3 * p0) / (M + ALPHA_F3)) / p0)


def sigmoid(z: np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(-np.clip(z, -35, 35)))


def fit_model(cfg: Sequence[str], C: float, train: Sequence[str], Z: Dict[str, np.ndarray], y: Dict[str, np.ndarray],
              cls_id: np.ndarray, k: int, keep: Optional[np.ndarray] = None) -> Dict[str, Any]:
    cols = cfg_cols(cfg)
    use3 = "F3" in cfg
    ref = F3Ref(train, y, cls_id, k) if use3 else None
    keep = np.ones(len(cls_id), dtype=bool) if keep is None else keep
    Xs, ys = [], []
    for d in train:
        X = Z[d][:, cols]
        if use3:
            X = np.column_stack([X, ref.table(d)[cls_id]])
        Xs.append(X[keep])
        ys.append(y[d][keep])
    X, yy = np.vstack(Xs), np.concatenate(ys)
    mu, sd = X.mean(0), X.std(0)
    sd[sd < 1e-12] = 1.0
    clf = LogisticRegression(C=C, max_iter=3000)
    clf.fit((X - mu) / sd, yy)
    return {"config": list(cfg), "C": C, "cols": cols, "features": cfg_features(cfg), "mu": mu, "sd": sd,
            "coef": clf.coef_[0].copy(), "intercept": float(clf.intercept_[0]),
            "f3_table": ref.table(None) if use3 else None, "prev": float(yy.mean())}


def score_model(m: Dict[str, Any], Zd: np.ndarray, cls_id: np.ndarray) -> np.ndarray:
    X = Zd[:, m["cols"]]
    if m["f3_table"] is not None:
        X = np.column_stack([X, np.asarray(m["f3_table"])[cls_id]])
    return ((X - np.asarray(m["mu"])) / np.asarray(m["sd"])) @ np.asarray(m["coef"]) + m["intercept"]


def pick_C(grid: Sequence[float], aucs: Dict[float, Dict[str, float]]) -> Tuple[float, Dict[str, Any]]:
    """Highest mean AUC, then the smallest C whose mean >= best - SE (SE = sd/sqrt(n) of the best)."""
    means = {C: float(np.mean(list(aucs[C].values()))) for C in grid}
    best = max(sorted(grid), key=lambda C: means[C])
    v = np.array(list(aucs[best].values()))
    se = float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 1 else 0.0
    chosen = min(C for C in grid if means[C] >= means[best] - se)
    return chosen, {"mean_auc_by_C": {str(C): means[C] for C in grid}, "best_C": best, "se_of_best": se}


def select_C(cfg: Sequence[str], train: Sequence[str], Z, y, classes, cls_id, k, grid: Sequence[float], seed: int):
    inner = es.make_folds(list(train), seed)
    aucs: Dict[float, Dict[str, float]] = {C: {} for C in grid}
    oof: Dict[float, Tuple[List[np.ndarray], List[np.ndarray]]] = {C: ([], []) for C in grid}
    for ifold in inner:
        itrain = [d for d in train if d not in ifold]
        if not itrain:
            continue
        ikeep = es._class_mask(classes, es._gold_classes(ifold, y, classes))
        for C in grid:
            m = fit_model(cfg, C, itrain, Z, y, cls_id, k, ikeep)
            for d in ifold:
                s = score_model(m, Z[d], cls_id)
                aucs[C][d] = es.auc_mw(s, y[d])
                oof[C][0].append(s)
                oof[C][1].append(y[d])
    C, diag = pick_C(grid, aucs)
    return C, oof, diag


def outer_cv(cfg: Sequence[str], dev: Sequence[str], Z, y, classes, cls_id, k, grid: Sequence[float] = C_GRID):
    """DEV grouped CV with nested C selection. Returns out-of-fold predictions per disease and the C chosen per fold."""
    dev = sorted(dev)
    folds = es.make_folds(dev, SEED)
    res: Dict[str, Dict[str, np.ndarray]] = {}
    chosen = []
    for fold in folds:
        train = [d for d in dev if d not in fold]
        C, oof, diag = select_C(cfg, train, Z, y, classes, cls_id, k, grid, SEED_INNER)
        s_in, y_in = np.concatenate(oof[C][0]), np.concatenate(oof[C][1])
        a, b = es.logit_fit(s_in, y_in)
        keep = es._class_mask(classes, es._gold_classes(fold, y, classes))
        m = fit_model(cfg, C, train, Z, y, cls_id, k, keep)
        for d in fold:
            s = score_model(m, Z[d], cls_id)
            res[d] = {"score": s, "platt": sigmoid(a + b * s), "prev": np.full(len(s), m["prev"])}
        chosen.append(C)
    return res, chosen, folds


def _outer_job(args):
    cfg, tag, grid, dev, Z, y, classes, cls_id, k = args
    res, chosen, _ = outer_cv(cfg, dev, Z, y, classes, cls_id, k, grid)
    return tag, res, chosen


def select_config(stats: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Protocol section 4: simplest configuration with mean >= best - SE of the best."""
    best = max(stats, key=lambda n: (stats[n]["mean_auc"], -stats[n]["complexity"][0], -stats[n]["complexity"][1]))
    thr = stats[best]["mean_auc"] - stats[best]["se"]
    elig = [n for n in stats if stats[n]["mean_auc"] >= thr]
    chosen = sorted(elig, key=lambda n: (stats[n]["complexity"][0], stats[n]["complexity"][1], -stats[n]["mean_auc"], n))[0]
    return {"best_by_mean": best, "threshold": thr, "eligible": sorted(elig), "selected": chosen}


# ============================================================================= evaluation
def evaluate(diseases: Sequence[str], preds: Dict[str, Dict[str, Dict[str, np.ndarray]]], y, pop, proxy,
             reference: Optional[str] = None) -> Dict[str, Any]:
    """Baselines + models with disease-level bootstrap and the three v1 criteria (es.decide) for each model."""
    diseases = sorted(diseases)
    idx = np.random.default_rng(SEED).integers(0, len(diseases), size=(es.N_BOOT, len(diseases)))
    pop_m, prx_m = es.per_disease_metrics(pop, y, diseases), es.per_disease_metrics(proxy, y, diseases)
    rnd_m = es.random_per_disease(y, diseases)
    base = {n: {"auc": es.macro_with_ci(m["auc"], idx), "precision_at_10": es.macro_with_ci(m["p10"], idx),
                "recall_at_20": es.macro_with_ci(m["r20"], idx)}
            for n, m in (("popularity", pop_m), ("proxy", prx_m), ("random", rnd_m))}
    best_name = max(("popularity", "proxy", "random"), key=lambda n: base[n]["precision_at_10"]["estimate"] or -1)
    best_m = {"popularity": pop_m, "proxy": prx_m, "random": rnd_m}[best_name]
    out: Dict[str, Any] = {"diseases": diseases, "n_diseases": len(diseases), "baselines": base,
                           "best_baseline_by_precision_at_10": best_name, "models": {}}
    ms: Dict[str, Dict[str, np.ndarray]] = {}
    for name, pred in preds.items():
        m = es.per_disease_metrics({d: pred[d]["score"] for d in diseases}, y, diseases)
        ms[name] = m
        e = {"auc": es.macro_with_ci(m["auc"], idx), "precision_at_10": es.macro_with_ci(m["p10"], idx),
             "recall_at_20": es.macro_with_ci(m["r20"], idx),
             "auc_minus_popularity": es.diff_with_ci(m["auc"], pop_m["auc"], idx),
             "auc_minus_proxy": es.diff_with_ci(m["auc"], prx_m["auc"], idx),
             "auc_minus_random": es.diff_with_ci(m["auc"], rnd_m["auc"], idx),
             "precision_at_10_minus_best_baseline": es.diff_with_ci(m["p10"], best_m["p10"], idx),
             "calibration_platt": es.calibration_summary(pred, y, diseases, idx, "platt"),
             "per_disease_auc": {d: es._r(v) for d, v in zip(diseases, m["auc"])}}
        e["decision"] = es.decide(e)
        out["models"][name] = e
    if reference is not None:
        for name in preds:
            if name != reference:
                out["models"][name]["auc_minus_" + reference] = es.diff_with_ci(ms[name]["auc"], ms[reference]["auc"], idx)
                out["models"][name]["precision_at_10_minus_" + reference] = es.diff_with_ci(ms[name]["p10"], ms[reference]["p10"], idx)
    return out


def paired_vs(a_auc: np.ndarray, b_auc: np.ndarray) -> Dict[str, Any]:
    idx = np.random.default_rng(SEED).integers(0, len(a_auc), size=(es.N_BOOT, len(a_auc)))
    return es.diff_with_ci(a_auc, b_auc, idx)


# ============================================================================= freezing
def canonical_hash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), default=_default).encode("utf-8")).hexdigest()


def _default(o: Any) -> Any:
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, tuple)):
        return sorted(o)
    raise TypeError(type(o))


def freeze(model: Dict[str, Any], platt: Tuple[float, float], classes_unique: Sequence[str], extra: Dict[str, Any]) -> Dict[str, Any]:
    body = {"config": model["config"], "features": model["features"], "C": model["C"], "cols": model["cols"],
            "mu": [float(v) for v in model["mu"]], "sd": [float(v) for v in model["sd"]], "coef": [float(v) for v in model["coef"]],
            "intercept": model["intercept"], "prev_dev": model["prev"], "platt_a": float(platt[0]), "platt_b": float(platt[1]),
            "f3_table": None if model["f3_table"] is None else [float(v) for v in model["f3_table"]],
            "classes_unique": list(classes_unique)}
    body.update(extra)
    return {"model": body, "sha256_of_model": canonical_hash(body)}


def thaw(frozen: Dict[str, Any]) -> Dict[str, Any]:
    if canonical_hash(frozen["model"]) != frozen["sha256_of_model"]:
        raise RuntimeError("frozen model hash mismatch")
    m = dict(frozen["model"])
    m["f3_table"] = None if m["f3_table"] is None else np.array(m["f3_table"])
    m["prev"] = m["prev_dev"]
    return m


def predict_frozen(m: Dict[str, Any], Zd: np.ndarray, cls_id: np.ndarray) -> Dict[str, np.ndarray]:
    s = score_model(m, Zd, cls_id)
    return {"score": s, "platt": sigmoid(m["platt_a"] + m["platt_b"] * s), "prev": np.full(len(s), m["prev_dev"])}


def rank_table(score: np.ndarray, names: Sequence[str], drugs: Sequence[str]) -> Dict[str, Any]:
    out = {}
    for dname in drugs:
        i = names.index(dname) if dname in names else None
        out[dname] = None if i is None else {"rank": int(1 + (score > score[i]).sum()), "score": es._r(score[i], 6)}
    return out


# ============================================================================= drivers
def _load(root: Path):
    b2 = json.loads((root / "out" / "engine" / "benchmark_v2.json").read_text(encoding="utf-8"))
    b1 = json.loads((root / "out" / "engine" / "benchmark.json").read_text(encoding="utf-8"))
    return b2, b1


def run_dev(root: Path = ROOT, workers: int = 6) -> Dict[str, Any]:
    out_dir = root / "out" / "engine"
    proto = root / "protocol" / "addendum_engine_v3.md"
    b2, b1 = _load(root)
    split = make_split(b2["diseases"])
    dev = split["dev"]
    print("DEV", len(dev), "HELD-OUT", len(split["heldout"]), file=sys.stderr, flush=True)
    D = build_Z(root, b2, b1, dev)
    Z, y, classes, cls_id, k, names = D["Z"], D["y"], D["classes"], D["cls_id"], D["n_classes"], D["names"]
    cls_unique = [str(c) for c in np.unique(classes)]

    # sanity: v1 base features equal the case-study cache (read-only)
    check = {"case_study_cache_available": False}
    pkl = root / "data" / "cache_engine" / "case_study_features.pkl"
    if pkl.exists():
        import pickle
        with pkl.open("rb") as fh:  # project-own cache written by engine_case_study
            blob = pickle.load(fh)
        perc = blob["per"]
        same = [bool(np.array_equal(perc[d]["primary"]["X"], Z[d][:, :8])) for d in dev if d in perc]
        check = {"case_study_cache_available": True, "n_compared": len(same), "all_identical": bool(all(same))}

    # ---------- candidate grid on DEV (nested C selection inside each outer training set)
    jobs = [(cfg, cfg_name(cfg), C_GRID, dev, Z, y, classes, cls_id, k) for cfg in CONFIGS]
    grid_res: Dict[str, Any] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for tag, res, chosen in ex.map(_outer_job, jobs):
            grid_res[tag] = (res, chosen)
            print("grid done", tag, file=sys.stderr, flush=True)
    stats: Dict[str, Dict[str, Any]] = {}
    for cfg in CONFIGS:
        tag = cfg_name(cfg)
        res, chosen = grid_res[tag]
        auc = np.array([es.auc_mw(res[d]["score"], y[d]) for d in dev])
        stats[tag] = {"config": list(cfg), "features": cfg_features(cfg), "complexity": list(complexity(cfg)),
                      "mean_auc": float(auc.mean()), "se": float(auc.std(ddof=1) / math.sqrt(len(auc))),
                      "C_chosen_per_outer_fold": chosen, "per_disease_auc": {d: es._r(v) for d, v in zip(dev, auc)},
                      "auc_ci95": es.macro_with_ci(auc, np.random.default_rng(SEED).integers(0, len(dev), size=(es.N_BOOT, len(dev))))["ci95"]}
    sel = select_config(stats)
    sel_cfg = tuple(stats[sel["selected"]]["config"])
    sel_tag = sel["selected"]
    print("selected", sel, file=sys.stderr, flush=True)

    # ---------- ablations on DEV: drop-one family, fixed C = 1
    abl_jobs = []
    for f in sel_cfg:
        sub = tuple(x for x in sel_cfg if x != f)
        abl_jobs.append((sub, "drop_" + f, C_GRID, dev, Z, y, classes, cls_id, k))
    abl_jobs.append((sel_cfg, "fixed_C_1", (1.0,), dev, Z, y, classes, cls_id, k))
    abl_res: Dict[str, Any] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for tag, res, chosen in ex.map(_outer_job, abl_jobs):
            abl_res[tag] = (res, chosen)
    sel_auc = np.array([es.auc_mw(grid_res[sel_tag][0][d]["score"], y[d]) for d in dev])
    ablation: Dict[str, Any] = {"selected": sel_tag, "selected_mean_auc": float(sel_auc.mean()), "drop_one_and_fixed_C": {}}
    for tag, (res, chosen) in sorted(abl_res.items()):
        a = np.array([es.auc_mw(res[d]["score"], y[d]) for d in dev])
        ablation["drop_one_and_fixed_C"][tag] = {"mean_auc": float(a.mean()), "C_chosen_per_outer_fold": chosen,
                                                  "selected_minus_this": paired_vs(sel_auc, a),
                                                  "features": cfg_features(tuple(x for x in sel_cfg if "drop_" + x != tag)) if tag.startswith("drop_") else cfg_features(sel_cfg)}
    ablation["add_one_from_grid"] = {t: {"mean_auc": stats[t]["mean_auc"], "minus_base": paired_vs(
        np.array([stats[t]["per_disease_auc"][d] for d in dev]), np.array([stats["base"]["per_disease_auc"][d] for d in dev]))}
        for t in ("base", "F1", "F2", "F3")}

    # ---------- v1 reference on DEV (fixed C = 1, arm A) and full evaluation of the selected configuration
    Xa = {d: Z[d][:, :8] for d in dev}
    v1_folds = es.make_folds(dev)
    v1_oof = es.oof_predictions(Xa, y, classes, v1_folds)
    sel_oof = grid_res[sel_tag][0]
    dev_eval = evaluate(dev, {"v1_weights_C1": v1_oof, "v3_selected": sel_oof}, y, D["pop"], D["proxy"], reference="v1_weights_C1")
    all_oof_eval = {t: grid_res[t][0] for t in grid_res}

    # ---------- final C, refit on all DEV, Platt from pooled outer OOF, freeze
    C_fin, _, diag = select_C(sel_cfg, dev, Z, y, classes, cls_id, k, C_GRID, SEED_FINAL)
    model = fit_model(sel_cfg, C_fin, dev, Z, y, cls_id, k, None)
    s_all = np.concatenate([sel_oof[d]["score"] for d in dev])
    y_all = np.concatenate([y[d] for d in dev])
    platt = es.logit_fit(s_all, y_all)
    extra = {"label": LABEL, "protocol": "protocol/addendum_engine_v3.md", "dev": dev, "heldout": split["heldout"],
             "selection": sel, "final_C_selection": diag, "n_universe_drugs": len(names)}
    frozen_v3 = freeze(model, platt, cls_unique, extra)
    m1 = fit_model((), 1.0, dev, Z, y, cls_id, k, None)
    s1 = np.concatenate([v1_oof[d]["score"] for d in dev])
    platt1 = es.logit_fit(s1, y_all)
    frozen_v1 = freeze(m1, platt1, cls_unique, {"label": "engine v1 weights: 8 base features, C = 1, fitted on DEV", "dev": dev,
                                                "heldout": split["heldout"]})
    (out_dir / "engine_v3_frozen.json").write_text(json.dumps(frozen_v3, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out_dir / "engine_v1_weights_frozen.json").write_text(json.dumps(frozen_v1, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    # ---------- dry run of the frozen scoring path on DEV diseases (smoke test of the held-out code path; no held-out data)
    t3, t1 = thaw(frozen_v3), thaw(frozen_v1)
    dry = {d: {"v3": predict_frozen(t3, Z[d], cls_id), "v1": predict_frozen(t1, Z[d], cls_id)} for d in dev[:3]}
    dry_ok = all(np.isfinite(v["v3"]["score"]).all() and np.isfinite(v["v1"]["score"]).all() for v in dry.values())

    # ---------- Wilms diagnostic: strict refit without Wilms
    wtrain = [d for d in dev if d != WILMS]
    wm3 = fit_model(sel_cfg, C_fin, wtrain, Z, y, cls_id, k, None)
    wm1 = fit_model((), 1.0, wtrain, Z, y, cls_id, k, None)
    s3, s1w = score_model(wm3, Z[WILMS], cls_id), score_model(wm1, Z[WILMS], cls_id)
    order3 = [names[i] for i in np.argsort(-s3, kind="stable")[:20]]
    wilms = {"label": "diagnostic, not evidence (Wilms is a DEV disease inspected before the protocol)",
             "training_diseases": len(wtrain), "gold1": sorted(names[i] for i in np.flatnonzero(y[WILMS])),
             "ranks_v3": rank_table(s3, names, WILMS_DRUGS), "ranks_v1_features_C1": rank_table(s1w, names, WILMS_DRUGS),
             "auc_v3": es._r(es.auc_mw(s3, y[WILMS])), "auc_v1_features_C1": es._r(es.auc_mw(s1w, y[WILMS])),
             "top20_v3": order3, "top20_v1": [names[i] for i in np.argsort(-s1w, kind="stable")[:20]],
             "n_universe_drugs": len(names)}
    cs = root / "out" / "engine" / "case_wilms.json"
    if cs.exists():
        wilms["case_wilms_json_sha256"] = es.sha256(cs)

    out = {
        "label": LABEL, "protocol": "protocol/addendum_engine_v3.md", "protocol_sha256_at_dev_run": es.sha256(proto),
        "status": "v1 (supported, 17 diseases) and v2 results are unchanged; this is post hoc development. HELD-OUT scored once, see 'heldout'.",
        "settings": {"seed": SEED, "n_bootstrap": es.N_BOOT, "C_grid": list(C_GRID), "alpha_f3": ALPHA_F3, "inner_seed": SEED_INNER,
                     "final_seed": SEED_FINAL, "thresholds": es.THRESHOLDS, "feature_names": ALL_NAMES + [F3_NAME],
                     "class_keywords": CLASS_KEYWORDS, "n_dev": len(dev), "n_heldout": len(split["heldout"]),
                     "new_network_requests": 0},
        "split": split, "feature_cache_check": check,
        "disease_info_dev": D["info"],
        "dev_cv_candidates": stats, "dev_selection": sel,
        "selected_configuration": {"name": sel_tag, "features": cfg_features(sel_cfg), "C": C_fin, "final_C_selection": diag,
                                   "coefficients": dict(zip(model["features"], [float(v) for v in model["coef"]])),
                                   "intercept": model["intercept"], "platt": {"a": float(platt[0]), "b": float(platt[1])}},
        "frozen": {"engine_v3_sha256": frozen_v3["sha256_of_model"], "engine_v1_weights_sha256": frozen_v1["sha256_of_model"],
                   "files": ["out/engine/engine_v3_frozen.json", "out/engine/engine_v1_weights_frozen.json"],
                   "dry_run_on_dev_finite": bool(dry_ok)},
        "dev_evaluation_selected_vs_v1": dev_eval,
        "dev_ablation": ablation,
        "wilms_diagnostic": wilms,
        "heldout": None,
    }
    (out_dir / "validation_v3.json").write_text(json.dumps(out, indent=2, sort_keys=True, ensure_ascii=False, default=_default) + "\n", encoding="utf-8")
    print(json.dumps({"selected": sel, "C": C_fin, "v3_sha": frozen_v3["sha256_of_model"], "v1_sha": frozen_v1["sha256_of_model"]}), flush=True)
    return out


def run_heldout(root: Path = ROOT) -> Dict[str, Any]:
    out_dir = root / "out" / "engine"
    lock = out_dir / "heldout_v3_scored.lock"
    if lock.exists():
        raise RuntimeError("HELD-OUT was already scored once (lock file exists); refusing a second scoring")
    f3 = json.loads((out_dir / "engine_v3_frozen.json").read_text(encoding="utf-8"))
    f1 = json.loads((out_dir / "engine_v1_weights_frozen.json").read_text(encoding="utf-8"))
    m3, m1 = thaw(f3), thaw(f1)                     # verifies both hashes
    b2, b1 = _load(root)
    split = make_split(b2["diseases"])
    if split["heldout"] != m3["heldout"] or split["dev"] != m3["dev"]:
        raise RuntimeError("split differs from the frozen file")
    held = split["heldout"]
    D = build_Z(root, b2, b1, held)                 # features only; no score yet
    lock.write_text(json.dumps({"scored_at": "first and only held-out scoring", "v3_sha256": f3["sha256_of_model"],
                                "v1_weights_sha256": f1["sha256_of_model"]}) + "\n", encoding="utf-8")
    Z, y, cls_id = D["Z"], D["y"], D["cls_id"]
    if [str(c) for c in np.unique(D["classes"])] != m3["classes_unique"]:
        raise RuntimeError("drug class vocabulary differs from the frozen file")
    preds = {"v1_weights": {d: predict_frozen(m1, Z[d], cls_id) for d in held},
             "v3": {d: predict_frozen(m3, Z[d], cls_id) for d in held}}
    ev = evaluate(held, preds, y, D["pop"], D["proxy"], reference="v1_weights")
    v3 = ev["models"]["v3"]
    diff = v3["auc_minus_v1_weights"]
    improved = bool(diff["ci95"][0] is not None and diff["ci95"][0] > 0)
    heldout = {"label": LABEL + "; HELD-OUT scored once", "diseases": held, "n_diseases": len(held),
               "frozen_v3_sha256": f3["sha256_of_model"], "frozen_v1_weights_sha256": f1["sha256_of_model"],
               "evaluation": ev, "disease_info": D["info"],
               "improvement_rule": "v3 improved over v1-weights only if the 95 % interval of the paired AUC difference (v3 minus v1-weights) lies above 0",
               "improvement_demonstrated": improved, "paired_auc_v3_minus_v1_weights": diff,
               "criteria_v3": v3["decision"], "criteria_v1_weights": ev["models"]["v1_weights"]["decision"]}
    path = out_dir / "validation_v3.json"
    out = json.loads(path.read_text(encoding="utf-8"))
    out["heldout"] = heldout
    out["status"] = out["status"].replace("HELD-OUT scored once, see 'heldout'.", "HELD-OUT scored once (see 'heldout').")
    path.write_text(json.dumps(out, indent=2, sort_keys=True, ensure_ascii=False, default=_default) + "\n", encoding="utf-8")
    print(json.dumps({"improvement_demonstrated": improved, "paired": diff, "criteria_v3": v3["decision"],
                      "criteria_v1": ev["models"]["v1_weights"]["decision"]}, default=_default, sort_keys=True), flush=True)
    return out


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in ("dev", "heldout"):
        raise SystemExit("usage: python -m track2.engine_v3 dev|heldout")
    (run_dev if sys.argv[1] == "dev" else run_heldout)()


if __name__ == "__main__":
    main()
