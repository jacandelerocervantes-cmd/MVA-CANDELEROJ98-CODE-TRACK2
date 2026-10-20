r"""Literature-grounded evidence engine: scorer and validation exactly as pre-registered in protocol/addendum_engine.md.

Run:  PYTHONPATH=src python -m track2.engine_score      (writes out/engine/validation.json, sorted keys)

What is implemented (nothing was tuned after seeing a result; all constants are fixed below):
  * Query = disease gene symbols only (never drug names): TITLE_ABS:"SYMBOL" OR ..., up to 300 abstracts, restricted by
    PUB_YEAR (<= 2024 primary, <= 2018 temporal).
  * Drug mentions: whole-word, case-insensitive match of the names of the approved drugs of data/universe.json
    (names of >= 4 characters, the drug_regex of the earlier proxy).
  * Features per (disease, drug): independent papers that mention the drug AND (case-sensitively) a disease gene symbol;
    papers with a sentence classed `benefit` by engine_extractors.rule_direction; papers with a human-patient / human-cell
    sentence (rule_system); number of `harm` sentences; Europe PMC record count of the drug (popularity); excess of the
    observed co-mention over the popularity expectation, log((c + 0.5) / (E + 0.5)), E = papers_with_gene * pop / N_epmc;
    PrimeKG prior (drug_protein targets among the disease genes, and disease genes with a PPI neighbour targeted by the
    drug; `indication` edges are never used, they define gold 1) - NOT used in the temporal analysis.
  * Model: L2 logistic regression (C = 1, standardised features) with out-of-fold predictions; folds by disease (seed
    20261015, 5 folds, fewer if < 10 diseases); drugs of one ATC level-3 class are never split between train and test of a
    fold; the Platt calibration map (protocol section 4.2) is fitted on inner out-of-fold scores of the training diseases
    only. An isotonic map fitted the same way is reported as a supplementary, non-decisive variant.
  * Metrics with disease-level bootstrap (2000 resamples, seed 20261015, percentile 95 % intervals, paired differences).

WHAT THE FEATURES CAN AND CANNOT CAPTURE: they are co-mention counts and rule-extracted cues (lexicon direction, system and
design cues); they are NOT the AI-assisted direction labels of the report, do not read evidence level, dose, species
fidelity or rules R1/R3/R4, and cannot separate a drug studied as a comparator from one that helps. The 'excess over
popularity' removes only the part of co-mention predictable from how often the drug is cited at all.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
from scipy.stats import rankdata
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression

from track2 import empirical_literature_control as elc
from track2 import engine_benchmark as eb
from track2.engine_extractors import rule_direction, rule_system

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261015
N_BOOT = 2000
MAX_PAPERS = 300
YEAR_PRIMARY = 2024
YEAR_TEMPORAL = 2018
MIN_PAPERS_SENSITIVITY = 10
SMOOTH = 0.5
C_REG = 1.0
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CACHE = ROOT / "data" / "cache_engine"
ARMS = {"A": (True, True), "B": (False, True), "C": (True, False), "D": (False, False)}  # (prior, popularity correction)
HUMAN_SYSTEMS = ("human patients", "human cells")
THRESHOLDS = {"auc_diff_min": 0.05, "calibration_slope": [0.8, 1.2]}


# ============================================================================= Europe PMC access
class _Limiter:
    def __init__(self, per_sec: float = 4.0):
        self.gap, self.next, self.lock = 1.0 / per_sec, 0.0, threading.Lock()

    def wait(self) -> None:
        with self.lock:
            now = time.monotonic()
            t = max(now, self.next)
            self.next = t + self.gap
        if t > now:
            time.sleep(t - now)


_LIMITER = _Limiter()


def default_getter(url: str, params: dict) -> dict:
    import requests
    for attempt in range(7):
        try:
            _LIMITER.wait()
            r = requests.get(url, params=params, timeout=60, headers={"User-Agent": "hselith-track2/0.3 (research)"})
            r.raise_for_status()
            return r.json()
        except Exception:  # noqa: BLE001 - SSL / transient errors: exponential backoff
            if attempt == 6:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def build_gene_query(symbols: Sequence[str], max_year: int) -> str:
    terms = " OR ".join(f'TITLE_ABS:"{s}"' for s in sorted(set(symbols)))
    return f"({terms}) AND PUB_YEAR:[1800 TO {max_year}] AND HAS_ABSTRACT:y"


def fetch_papers(symbols: Sequence[str], max_year: int, cache_dir: Path = CACHE / "epmc",
                 getter: Optional[Callable[[str, dict], dict]] = None, max_n: int = MAX_PAPERS) -> List[dict]:
    """Up to max_n papers {pid, year, text(title + abstract)} for the gene-symbol query; cached on disk."""
    query = build_gene_query(symbols, max_year)
    cache_dir.mkdir(parents=True, exist_ok=True)
    f = cache_dir / (hashlib.sha256(f"{query}|{max_n}".encode()).hexdigest()[:24] + ".json")
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    getter = getter or default_getter
    papers: List[dict] = []
    seen: Set[str] = set()
    cursor = "*"
    while len(papers) < max_n:
        data = getter(EPMC, {"query": query, "format": "json", "resultType": "core", "pageSize": 100, "cursorMark": cursor})
        items = (data.get("resultList") or {}).get("result") or []
        nxt = data.get("nextCursorMark")
        n_before = len(papers)
        for it in items:
            pid = f"{it.get('source', '')}:{it.get('id') or it.get('pmid') or ''}"
            if not it.get("abstractText") or pid in seen:
                continue
            seen.add(pid)
            try:
                year = int(it.get("pubYear"))
            except (TypeError, ValueError):
                year = None
            papers.append({"pid": pid, "year": year, "text": f"{it.get('title') or ''}. {it['abstractText']}"})
        if len(items) < 100 or not nxt or nxt == cursor or (len(papers) == n_before and not items):
            break
        cursor = nxt
    papers = papers[:max_n]
    f.write_text(json.dumps(papers), encoding="utf-8")
    return papers


def _phrase(name: str) -> str:
    return '"' + name.replace('"', " ").strip() + '"'


def fetch_popularity(names: Sequence[str], max_year: Optional[int], cache_file: Path,
                     getter: Optional[Callable[[str, dict], dict]] = None, workers: int = 4) -> Tuple[Dict[str, int], int]:
    """Europe PMC hit count of the exact phrase of each drug name (all years, or PUB_YEAR <= max_year) and the total
    number of records for the same period. Cached incrementally in cache_file."""
    getter = getter or default_getter
    cache: Dict[str, int] = json.loads(cache_file.read_text(encoding="utf-8")) if cache_file.exists() else {}
    suffix = f" AND PUB_YEAR:[1800 TO {max_year}]" if max_year else ""
    lock = threading.Lock()

    def hits(q: str) -> int:
        d = getter(EPMC, {"query": q, "format": "json", "resultType": "lite", "pageSize": 1})
        return int(d.get("hitCount", 0))

    def work(name: str) -> None:
        if name in cache:
            return
        n = hits(_phrase(name) + suffix)
        with lock:
            cache[name] = n
            if len(cache) % 200 == 0:
                cache_file.write_text(json.dumps(cache), encoding="utf-8")

    todo = [n for n in names if n not in cache]
    if "__total__" not in cache:
        cache["__total__"] = hits(f"PUB_YEAR:[1800 TO {max_year or 2100}]")
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(work, todo))
    cache_file.write_text(json.dumps(cache), encoding="utf-8")
    total = cache["__total__"]
    return {n: cache[n] for n in names}, total


# ============================================================================= PrimeKG prior
def load_prior_graph(root: Path = ROOT, cache_file: Path = CACHE / "primekg_prior_edges.json") -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """drug (lower-case name) -> target gene symbols (drug_protein) and gene -> PPI neighbours. No indication edges."""
    if cache_file.exists():
        raw = json.loads(cache_file.read_text(encoding="utf-8"))
        return {k: set(v) for k, v in raw["targets"].items()}, {k: set(v) for k, v in raw["ppi"].items()}
    import pandas as pd
    d = root / "data" / "sources" / "primekg"
    nodes = pd.read_csv(d / "nodes.tab")
    ntype = dict(zip(nodes["node_index"], nodes["node_type"]))
    nname = dict(zip(nodes["node_index"], nodes["node_name"]))
    targets: Dict[str, Set[str]] = {}
    ppi: Dict[str, Set[str]] = {}
    for ch in pd.read_csv(d / "edges.csv", chunksize=2_000_000, usecols=["relation", "x_index", "y_index"]):
        ch = ch[ch["relation"].isin(["drug_protein", "protein_protein"])]
        for rel, x, y in zip(ch["relation"], ch["x_index"], ch["y_index"]):
            if rel == "drug_protein":
                for a, b in ((x, y), (y, x)):
                    if ntype[a] == "drug" and ntype[b] == "gene/protein":
                        targets.setdefault(str(nname[a]).lower(), set()).add(str(nname[b]))
            else:
                ppi.setdefault(str(nname[x]), set()).add(str(nname[y]))
                ppi.setdefault(str(nname[y]), set()).add(str(nname[x]))
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps({"targets": {k: sorted(v) for k, v in targets.items()},
                                      "ppi": {k: sorted(v) for k, v in ppi.items()}}), encoding="utf-8")
    return targets, ppi


def prior_features(genes: Sequence[str], drug_names: Sequence[str], targets: Dict[str, Set[str]],
                   ppi: Dict[str, Set[str]]) -> np.ndarray:
    """(n_drugs, 2): log1p(#disease genes directly targeted), log1p(#disease genes with a PPI neighbour targeted)."""
    gs = set(genes)
    out = np.zeros((len(drug_names), 2))
    for i, n in enumerate(drug_names):
        t = targets.get(n, set())
        if not t:
            continue
        out[i, 0] = math.log1p(len(t & gs))
        out[i, 1] = math.log1p(sum(1 for g in gs if ppi.get(g, set()) & t))
    return out


# ============================================================================= literature features
def gene_regex(symbols: Sequence[str]) -> "re.Pattern[str]":
    return re.compile("|".join(r"(?<![\w])" + re.escape(s) + r"(?![\w])" for s in sorted(set(symbols), key=len, reverse=True)))


def name_gate(names: Sequence[str]) -> Set[str]:
    """First word-run of every matchable drug name; a sentence without any of them cannot contain a drug (speed only)."""
    out: Set[str] = set()
    for n in names:
        if len(n) >= 4:
            m = re.search(r"\w+", n.lower())
            if m:
                out.add(m.group(0))
    return out


def literature_counts(papers: Sequence[dict], symbols: Sequence[str], dre: "re.Pattern[str]",
                      gate: Optional[Set[str]] = None) -> Dict[str, Any]:
    """Per drug (lower-case name): n_papers, n_benefit, n_human, n_harm; plus the number of papers with a gene symbol."""
    grx = gene_regex(symbols)
    agg: Dict[str, Dict[str, Set[str]]] = {}
    harm: Dict[str, int] = {}
    n_gene = 0
    for p in papers:
        text = p["text"]
        if not grx.search(text):
            continue
        n_gene += 1
        for sent in elc.split_sentences(text):
            if gate is not None and not (set(re.findall(r"\w+", sent.lower())) & gate):
                continue
            for name in {m.group(0).lower() for m in dre.finditer(sent)}:
                a = agg.setdefault(name, {"papers": set(), "benefit": set(), "human": set()})
                a["papers"].add(p["pid"])
                direction = rule_direction(sent, name)
                if direction == "benefit":
                    a["benefit"].add(p["pid"])
                elif direction == "harm":
                    harm[name] = harm.get(name, 0) + 1
                if rule_system(sent) in HUMAN_SYSTEMS:
                    a["human"].add(p["pid"])
    counts = {n: {"n_papers": len(a["papers"]), "n_benefit": len(a["benefit"]), "n_human": len(a["human"]),
                  "n_harm": harm.get(n, 0)} for n, a in agg.items()}
    return {"counts": counts, "n_papers_with_gene": n_gene}


FEATURE_NAMES = ["log1p_n_papers", "log1p_n_benefit", "log1p_n_human", "log1p_n_harm", "log_popularity", "excess_over_popularity",
                 "prior_direct_targets", "prior_ppi_targets"]


def disease_matrix(lit: Dict[str, Any], names: Sequence[str], pop: Dict[str, int], total: int,
                   prior: Optional[np.ndarray]) -> np.ndarray:
    """(n_drugs, 8) matrix; prior columns are zero when prior is None (they are dropped by the arm anyway)."""
    n_g = lit["n_papers_with_gene"]
    X = np.zeros((len(names), len(FEATURE_NAMES)))
    for i, n in enumerate(names):
        c = lit["counts"].get(n, {})
        npap = c.get("n_papers", 0)
        X[i, 0] = math.log1p(npap)
        X[i, 1] = math.log1p(c.get("n_benefit", 0))
        X[i, 2] = math.log1p(c.get("n_human", 0))
        X[i, 3] = math.log1p(c.get("n_harm", 0))
        X[i, 4] = math.log1p(pop.get(n, 0))
        expected = n_g * pop.get(n, 0) / max(total, 1)
        X[i, 5] = math.log((npap + SMOOTH) / (expected + SMOOTH))
    if prior is not None:
        X[:, 6:8] = prior
    return X


def arm_columns(use_prior: bool, use_pop: bool) -> List[int]:
    cols = [0, 1, 2, 3]
    if use_pop:
        cols += [4, 5]
    if use_prior:
        cols += [6, 7]
    return cols


# ============================================================================= metrics
def auc_mw(score: np.ndarray, y: np.ndarray) -> float:
    pos = int(y.sum())
    neg = len(y) - pos
    if pos == 0 or neg == 0:
        return float("nan")
    r = rankdata(score)
    return float((r[y == 1].sum() - pos * (pos + 1) / 2) / (pos * neg))


def expected_hits_at_k(score: np.ndarray, y: np.ndarray, k: int) -> float:
    """Expected number of positives in the top k when ties are broken uniformly at random (exact, group-wise)."""
    uniq, inv = np.unique(score, return_inverse=True)
    size = np.bincount(inv)
    npos = np.bincount(inv, weights=y)
    left, hits = k, 0.0
    for g in range(len(uniq) - 1, -1, -1):
        take = min(left, int(size[g]))
        hits += take * npos[g] / size[g]
        left -= take
        if left <= 0:
            break
    return hits


def precision_at_k(score: np.ndarray, y: np.ndarray, k: int = 10) -> float:
    return expected_hits_at_k(score, y, k) / k if y.sum() > 0 else float("nan")


def recall_at_k(score: np.ndarray, y: np.ndarray, k: int = 20) -> float:
    return expected_hits_at_k(score, y, k) / y.sum() if y.sum() > 0 else float("nan")


def random_metrics(y: np.ndarray) -> Tuple[float, float, float]:
    """Exact expectations for uniformly random scores."""
    n, p = len(y), int(y.sum())
    if p == 0:
        return float("nan"), float("nan"), float("nan")
    auc = 0.5 if p < n else float("nan")
    return auc, p / n, min(20, n) / n


def logit_fit(x: np.ndarray, y: np.ndarray, w: Optional[np.ndarray] = None, iters: int = 50, ridge: float = 1e-9) -> Tuple[float, float]:
    """Weighted unpenalised logistic regression y ~ a + b x (Newton). Returns (a, b)."""
    w = np.ones_like(x, dtype=float) if w is None else w.astype(float)
    a, b = 0.0, 0.0
    for _ in range(iters):
        z = a + b * x
        p = 1 / (1 + np.exp(-np.clip(z, -35, 35)))
        g0, g1 = (w * (y - p)).sum(), (w * (y - p) * x).sum()
        s = w * p * (1 - p)
        h00, h01, h11 = s.sum() + ridge, (s * x).sum(), (s * x * x).sum() + ridge
        det = h00 * h11 - h01 * h01
        if det <= 0:
            break
        da, db = (h11 * g0 - h01 * g1) / det, (h00 * g1 - h01 * g0) / det
        a, b = a + da, b + db
        if abs(da) + abs(db) < 1e-9:
            break
    return float(a), float(b)


def _logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


# ============================================================================= grouped cross-validation
def make_folds(diseases: Sequence[str], seed: int = SEED) -> List[List[str]]:
    n = len(diseases)
    k = 5 if n >= 10 else max(2, min(5, n // 2))
    if n < 2:
        return [list(diseases)]
    order = sorted(diseases)
    random.Random(seed).shuffle(order)
    return [sorted(order[i::k]) for i in range(k)]


def drug_classes(universe: Dict[str, dict], names: Sequence[str]) -> List[str]:
    out = []
    for n in names:
        c = eb.atc_class(universe[n].get("atc") or [])
        out.append(c if c else f"nc:{n}")
    return out


def _fit_scorer(Xs: List[np.ndarray], ys: List[np.ndarray], keep: List[np.ndarray]):
    X = np.vstack([x[m] for x, m in zip(Xs, keep)])
    y = np.concatenate([yy[m] for yy, m in zip(ys, keep)])
    mu, sd = X.mean(0), X.std(0)
    sd[sd < 1e-12] = 1.0
    clf = LogisticRegression(C=C_REG, max_iter=2000)
    clf.fit((X - mu) / sd, y)
    return lambda Z: clf.decision_function((Z - mu) / sd), float(y.mean())


def _class_mask(classes: np.ndarray, bad: Set[str]) -> np.ndarray:
    return ~np.isin(classes, list(bad)) if bad else np.ones(len(classes), dtype=bool)


def _gold_classes(diseases: Sequence[str], ys: Dict[str, np.ndarray], classes: np.ndarray) -> Set[str]:
    return {c for d in diseases for c in classes[ys[d] == 1]}


def oof_predictions(X: Dict[str, np.ndarray], y: Dict[str, np.ndarray], classes: np.ndarray,
                    folds: List[List[str]]) -> Dict[str, Dict[str, np.ndarray]]:
    """Out-of-fold raw score, Platt probability, isotonic probability and training prevalence, per held-out disease."""
    res: Dict[str, Dict[str, np.ndarray]] = {}
    all_d = sorted(d for f in folds for d in f)
    for fold in folds:
        train = [d for d in all_d if d not in fold]
        bad = _gold_classes(fold, y, classes)
        keep = _class_mask(classes, bad)
        # inner out-of-fold raw scores of the training diseases (calibration is learned only from them)
        inner_scores, inner_y = [], []
        inner = make_folds(train, SEED + 1) if len(train) >= 2 else [train]
        for ifold in inner:
            itrain = [d for d in train if d not in ifold]
            if not itrain:
                continue
            ikeep = _class_mask(classes, _gold_classes(ifold, y, classes))
            score, _ = _fit_scorer([X[d] for d in itrain], [y[d] for d in itrain], [ikeep] * len(itrain))
            for d in ifold:
                inner_scores.append(score(X[d]))
                inner_y.append(y[d])
        s_in, y_in = np.concatenate(inner_scores), np.concatenate(inner_y)
        a, b = logit_fit(s_in, y_in)
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0).fit(s_in, y_in)
        score, prev = _fit_scorer([X[d] for d in train], [y[d] for d in train], [keep] * len(train))
        for d in fold:
            s = score(X[d])
            res[d] = {"score": s, "platt": 1 / (1 + np.exp(-np.clip(a + b * s, -35, 35))),
                      "isotonic": iso.predict(s), "prev": np.full(len(s), prev)}
    return res


# ============================================================================= bootstrap summaries
def _ci(samples: np.ndarray) -> Tuple[Optional[float], Optional[float]]:
    s = samples[~np.isnan(samples)]
    if len(s) == 0:
        return None, None
    return float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5))


def _r(x: Optional[float], nd: int = 6) -> Optional[float]:
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else round(float(x), nd)


def macro_with_ci(vals: np.ndarray, idx: np.ndarray) -> Dict[str, Any]:
    ok = ~np.isnan(vals)
    est = float(np.nanmean(vals)) if ok.any() else None
    with np.errstate(all="ignore"):
        boots = np.array([np.nanmean(v) if (~np.isnan(v)).any() else np.nan for v in vals[idx]])
    lo, hi = _ci(boots)
    return {"estimate": _r(est), "ci95": [_r(lo), _r(hi)], "n_diseases": int(ok.sum())}


def diff_with_ci(a: np.ndarray, b: np.ndarray, idx: np.ndarray) -> Dict[str, Any]:
    d = a - b
    out = macro_with_ci(d, idx)
    out["p_boot_le_0"] = None
    ok = ~np.isnan(d)
    if ok.any():
        with np.errstate(all="ignore"):
            boots = np.array([np.nanmean(v) if (~np.isnan(v)).any() else np.nan for v in d[idx]])
        boots = boots[~np.isnan(boots)]
        out["p_boot_le_0"] = _r(float((boots <= 0).mean())) if len(boots) else None
    return out


def calibration_summary(pred: Dict[str, Dict[str, np.ndarray]], y: Dict[str, np.ndarray], diseases: Sequence[str],
                        idx: np.ndarray, key: str = "platt") -> Dict[str, Any]:
    """Pooled slope (logit of the calibrated probability) and Brier score with disease-level bootstrap."""
    D = len(diseases)
    p = np.concatenate([pred[d][key] for d in diseases])
    prev = np.concatenate([pred[d]["prev"] for d in diseases])
    yy = np.concatenate([y[d] for d in diseases]).astype(float)
    di = np.concatenate([np.full(len(y[d]), i) for i, d in enumerate(diseases)])
    x = _logit(p)
    slope = logit_fit(x, yy)[1]
    brier, brier_c = float(((p - yy) ** 2).mean()), float(((prev - yy) ** 2).mean())
    cnt = np.stack([np.bincount(r, minlength=D) for r in idx]).astype(float)
    slopes, bs, bcs = [], [], []
    for c in cnt:
        w = c[di]
        if w.sum() == 0 or (w * yy).sum() == 0:
            continue
        slopes.append(logit_fit(x, yy, w)[1])
        bs.append(float((w * (p - yy) ** 2).sum() / w.sum()))
        bcs.append(float((w * (prev - yy) ** 2).sum() / w.sum()))
    ci = lambda v: [_r(c) for c in _ci(np.array(v))]  # noqa: E731
    return {"calibration_slope": {"estimate": _r(slope), "ci95": ci(slopes)},
            "brier": {"estimate": _r(brier), "ci95": ci(bs)},
            "brier_constant_prevalence": {"estimate": _r(brier_c), "ci95": ci(bcs)},
            "n_pairs": int(len(yy)), "n_positive_pairs": int(yy.sum())}


# ============================================================================= analysis of one setting
def per_disease_metrics(score_by_d: Dict[str, np.ndarray], y: Dict[str, np.ndarray], diseases: Sequence[str],
                        valid_mask: Optional[Dict[str, np.ndarray]] = None) -> Dict[str, np.ndarray]:
    auc, p10, r20 = [], [], []
    for d in diseases:
        s, yy = score_by_d[d], y[d]
        if valid_mask is not None:
            m = valid_mask[d]
            s, yy = s[m], yy[m]
        auc.append(auc_mw(s, yy))
        p10.append(precision_at_k(s, yy, 10))
        r20.append(recall_at_k(s, yy, 20))
    return {"auc": np.array(auc), "p10": np.array(p10), "r20": np.array(r20)}


def random_per_disease(y: Dict[str, np.ndarray], diseases: Sequence[str], valid_mask: Optional[Dict[str, np.ndarray]] = None) -> Dict[str, np.ndarray]:
    rows = [random_metrics(y[d][valid_mask[d]] if valid_mask is not None else y[d]) for d in diseases]
    return {k: np.array([r[i] for r in rows]) for i, k in enumerate(("auc", "p10", "r20"))}


def analyse(diseases: List[str], X: Dict[str, np.ndarray], y: Dict[str, np.ndarray], base_pop: Dict[str, np.ndarray],
            base_proxy: Dict[str, np.ndarray], classes: np.ndarray, arms: Dict[str, Tuple[bool, bool]],
            drop_mask: Optional[Dict[str, np.ndarray]] = None, with_calibration: bool = True) -> Dict[str, Any]:
    """Cross-validate every arm and compute all metrics for one (benchmark subset, gold, literature period).

    drop_mask: per disease boolean (True = keep) applied at evaluation only (sensitivity: removing gold-1 drugs)."""
    diseases = sorted(diseases)
    folds = make_folds(diseases)
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, len(diseases), size=(N_BOOT, len(diseases)))
    pop_m = per_disease_metrics(base_pop, y, diseases, drop_mask)
    prx_m = per_disease_metrics(base_proxy, y, diseases, drop_mask)
    rnd_m = random_per_disease(y, diseases, drop_mask)
    baselines = {name: {"auc": macro_with_ci(m["auc"], idx), "precision_at_10": macro_with_ci(m["p10"], idx),
                        "recall_at_20": macro_with_ci(m["r20"], idx)}
                 for name, m in (("popularity", pop_m), ("proxy", prx_m), ("random", rnd_m))}
    best_name = max(("popularity", "proxy", "random"), key=lambda n: baselines[n]["precision_at_10"]["estimate"] or -1)
    best_m = {"popularity": pop_m, "proxy": prx_m, "random": rnd_m}[best_name]
    result: Dict[str, Any] = {"diseases": diseases, "n_folds": len(folds), "folds": folds, "baselines": baselines,
                              "best_baseline_by_precision_at_10": best_name, "arms": {}}
    for arm, (use_prior, use_pop) in arms.items():
        cols = arm_columns(use_prior, use_pop)
        Xa = {d: X[d][:, cols] for d in diseases}
        pred = oof_predictions(Xa, y, classes, folds)
        m = per_disease_metrics({d: pred[d]["score"] for d in diseases}, y, diseases, drop_mask)
        out: Dict[str, Any] = {
            "features": [FEATURE_NAMES[c] for c in cols],
            "auc": macro_with_ci(m["auc"], idx),
            "precision_at_10": macro_with_ci(m["p10"], idx),
            "recall_at_20": macro_with_ci(m["r20"], idx),
            "auc_minus_popularity": diff_with_ci(m["auc"], pop_m["auc"], idx),
            "auc_minus_proxy": diff_with_ci(m["auc"], prx_m["auc"], idx),
            "auc_minus_random": diff_with_ci(m["auc"], rnd_m["auc"], idx),
            "precision_at_10_minus_best_baseline": diff_with_ci(m["p10"], best_m["p10"], idx),
            "per_disease_auc": {d: _r(v) for d, v in zip(diseases, m["auc"])},
        }
        if with_calibration:
            out["calibration_platt"] = calibration_summary(pred, y, diseases, idx, "platt")
            out["calibration_isotonic_supplementary"] = calibration_summary(pred, y, diseases, idx, "isotonic")
        result["arms"][arm] = out
    return result


def decide(arm_a: Dict[str, Any]) -> Dict[str, Any]:
    """The three pre-registered success criteria (protocol section 4), evaluated on arm A, whatever the result."""
    d = arm_a["auc_minus_popularity"]
    lo = d["ci95"][0]
    c1 = bool(d["estimate"] is not None and lo is not None and d["estimate"] >= THRESHOLDS["auc_diff_min"] and lo > 0)
    sl = arm_a["calibration_platt"]["calibration_slope"]["estimate"]
    c2 = bool(sl is not None and THRESHOLDS["calibration_slope"][0] <= sl <= THRESHOLDS["calibration_slope"][1])
    p = arm_a["precision_at_10_minus_best_baseline"]
    c3 = bool(p["estimate"] is not None and p["estimate"] > 0)
    return {
        "criterion_1_auc_minus_popularity_ge_0.05_and_ci_lower_gt_0": {"met": c1, "estimate": d["estimate"], "ci95": d["ci95"]},
        "criterion_2_calibration_slope_in_0.8_1.2": {"met": c2, "estimate": sl},
        "criterion_3_precision_at_10_above_best_baseline": {"met": c3, "difference_estimate": p["estimate"]},
        "verdict": "supported" if (c1 and c2 and c3) else "not supported (ranking stays exploratory)",
    }


# ============================================================================= driver
def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(root: Path = ROOT, getter: Optional[Callable[[str, dict], dict]] = None, out_path: Optional[Path] = None) -> Dict[str, Any]:
    bench_path = root / "out" / "engine" / "benchmark.json"
    bench = json.loads(bench_path.read_text(encoding="utf-8"))
    drugs = eb.load_universe(root / "data" / "universe.json")
    universe = eb.universe_index(drugs)
    names = sorted(universe)
    dre = elc.drug_regex(names)
    gate = name_gate(names)
    index = {n: i for i, n in enumerate(names)}
    classes = np.array(drug_classes(universe, names))
    targets, ppi = load_prior_graph(root, root / "data" / "cache_engine" / "primekg_prior_edges.json")
    pop_all, tot_all = fetch_popularity(names, None, root / "data" / "cache_engine" / "epmc_popularity_all.json", getter)
    pop_18, tot_18 = fetch_popularity(names, YEAR_TEMPORAL, root / "data" / "cache_engine" / "epmc_popularity_2018.json", getter)
    cache = root / "data" / "cache_engine" / "epmc"

    per: Dict[str, Dict[str, Any]] = {}
    for dis in bench["diseases"]:
        name = dis["name"]
        symbols = sorted({g["symbol"] for g in dis["genes"]})
        mech = sorted(set(symbols) | {t for g in dis["genes"] if g["mechanism_class"] != "other" for t in g["matched_terms"]})
        rec: Dict[str, Any] = {"symbols": symbols, "mech_terms": mech}
        gold1 = np.zeros(len(names)); gold2 = np.zeros(len(names))
        for g in dis["gold1"]:
            if g["name"].lower() in index:
                gold1[index[g["name"].lower()]] = 1
        for g in dis["gold2"]:
            if g["name"].lower() in index:
                gold2[index[g["name"].lower()]] = 1
        rec["gold1"], rec["gold2"] = gold1, gold2
        for period, year, pop, tot in (("primary", YEAR_PRIMARY, pop_all, tot_all), ("temporal", YEAR_TEMPORAL, pop_18, tot_18)):
            papers = fetch_papers(symbols, year, cache, getter)
            lit = literature_counts(papers, symbols, dre, gate)
            prior = prior_features(symbols, names, targets, ppi) if period == "primary" else None
            proxy_counts = elc.score_abstracts([(p["pid"], p["text"]) for p in papers], mech, dre)
            rec[period] = {
                "n_papers_retrieved": len(papers), "n_papers_with_gene": lit["n_papers_with_gene"],
                "X": disease_matrix(lit, names, pop, tot, prior),
                "popularity": np.array([pop.get(n, 0) for n in names], dtype=float),
                "proxy": np.array([proxy_counts.get(n, 0) for n in names], dtype=float),
            }
        per[name] = rec

    def run_set(sel: List[str]) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        y1 = {d: per[d]["gold1"] for d in sel}
        out["primary_gold1"] = analyse(sel, {d: per[d]["primary"]["X"] for d in sel}, y1,
                                       {d: per[d]["primary"]["popularity"] for d in sel},
                                       {d: per[d]["primary"]["proxy"] for d in sel}, classes, ARMS)
        out["primary_gold1"]["decision_arm_A"] = decide(out["primary_gold1"]["arms"]["A"])
        sel2 = [d for d in sel if per[d]["gold2"].sum() >= 1]
        temporal_arms = {k: v for k, v in ARMS.items() if not v[0]}
        y2 = {d: per[d]["gold2"] for d in sel2}
        tX = {d: per[d]["temporal"]["X"] for d in sel2}
        tp = {d: per[d]["temporal"]["popularity"] for d in sel2}
        tq = {d: per[d]["temporal"]["proxy"] for d in sel2}
        out["temporal_gold2"] = analyse(sel2, tX, y2, tp, tq, classes, temporal_arms)
        out["temporal_gold2"]["note"] = ("Literature up to 2018 only; knowledge-graph prior excluded, so only arms B and D are valid; "
                                         "arms A and C are not computed. Diseases with no gold 2 drug are excluded.")
        out["temporal_gold2"]["diseases_excluded_empty_gold2"] = sorted(set(sel) - set(sel2))
        keep = {d: per[d]["gold1"] == 0 for d in sel2}
        out["temporal_gold2_excluding_gold1_drugs"] = analyse(
            [d for d in sel2 if (per[d]["gold2"][keep[d]].sum() >= 1)], tX, y2, tp, tq, classes, temporal_arms, drop_mask=keep) \
            if sel2 else None
        return out

    all_d = sorted(per)
    main_res = run_set(all_d)
    big = [d for d in all_d if per[d]["primary"]["n_papers_retrieved"] >= MIN_PAPERS_SENSITIVITY]
    sens = run_set(big) if len(big) >= 4 else None

    bench_sha = sha256(bench_path)
    out = {
        "benchmark_sha256": bench_sha,
        "benchmark_sha256_matches_addendum": bench_sha == "48c53a059a3f35218a3fa57efaf26d8d28d4002d86dde536ab4432211dde9a74",
        "addendum_sha256": sha256(root / "protocol" / "addendum_engine.md"),
        "protocol": "protocol/addendum_engine.md (pre-registered 2026-10-15); thresholds and endpoint unchanged",
        "settings": {"seed": SEED, "n_bootstrap": N_BOOT, "max_papers": MAX_PAPERS, "year_primary_max": YEAR_PRIMARY,
                     "year_temporal_max": YEAR_TEMPORAL, "logistic_C": C_REG, "smoothing": SMOOTH,
                     "n_universe_drugs": len(names), "n_benchmark_diseases": len(all_d),
                     "arms": {k: {"knowledge_graph_prior": v[0], "popularity_correction": v[1]} for k, v in ARMS.items()},
                     "thresholds": THRESHOLDS, "feature_names": FEATURE_NAMES},
        "per_disease_retrieval": {d: {"symbols": per[d]["symbols"], "n_gold1": int(per[d]["gold1"].sum()),
                                      "n_gold2": int(per[d]["gold2"].sum()),
                                      "n_papers_retrieved_primary": per[d]["primary"]["n_papers_retrieved"],
                                      "n_papers_with_gene_primary": per[d]["primary"]["n_papers_with_gene"],
                                      "n_papers_retrieved_temporal": per[d]["temporal"]["n_papers_retrieved"],
                                      "n_papers_with_gene_temporal": per[d]["temporal"]["n_papers_with_gene"]} for d in all_d},
        "main": main_res,
        "sensitivity_min_10_retrieved_papers": {
            "label": "SENSITIVITY ANALYSIS (not the pre-registered primary analysis): diseases with >= 10 retrieved primary papers; "
                     "the whole cross-validation is re-run on this subset", "diseases": big, "results": sens},
        "deviations_and_clarifications": [
            "Calibration map: the protocol (section 4.2) specifies a logistic (Platt) map; this is the decisive one. An isotonic map "
            "(requested in the task brief) is reported as supplementary only, never used for the decision.",
            "Platt/isotonic maps are fitted on inner out-of-fold scores of the training diseases of each outer fold (nested grouped CV).",
            "Random baseline: exact expectations (AUC 0.5, precision@k = positives/N, recall@20 = 20/N), no sampling.",
            "Popularity counts all records up to the retrieval date (primary) or PUB_YEAR <= 2018 (temporal); the primary literature "
            "is restricted to <= 2024, a negligible mismatch.",
            "Popularity expectation uses all Europe PMC records as the denominator while co-mention is measured in abstracts; "
            "the learned weight absorbs the constant.",
            "Gene presence in a paper is a case-sensitive whole-word match of the HGNC symbol in title or abstract.",
        ],
        "what_features_can_and_cannot_capture": {
            "can": ["how many independent abstracts mention the drug together with a disease gene",
                    "whether such a paper has a sentence the lexicon rules read as benefit or harm",
                    "whether a drug sentence has a human patient or human cell cue",
                    "co-mention above what the drug's overall popularity predicts", "network proximity of drug targets to disease genes (non-temporal only)"],
            "cannot": ["the AI-assisted direction labels of the report (these are not used)",
                       "negation or comparator roles beyond the rule lexicon; a sentence 'benefit' may concern another entity",
                       "evidence level, dose, species fidelity, rules R1/R3/R4", "brand names, salts or synonyms outside the universe name",
                       "drugs with names shorter than 4 characters (not matched in text)", "full texts (abstracts only, max 300 by relevance)"]},
        "known_limits": ["Benchmark of 17 diseases: bootstrap intervals are wide.", "Gold 1 is incomplete and partly circular.",
                         "Gold 2 is a name match in intervention fields and may include comparators.",
                         "PrimeKG-based gold may reflect post-2018 knowledge; the prior is excluded from the temporal analysis."],
    }
    out_path = out_path or root / "out" / "engine" / "validation.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    return out


def main() -> None:
    res = run()
    d = res["main"]["primary_gold1"]
    print(json.dumps({"decision_arm_A": d["decision_arm_A"], "baselines": d["baselines"]}, indent=1, sort_keys=True)[:4000])


if __name__ == "__main__":
    main()
