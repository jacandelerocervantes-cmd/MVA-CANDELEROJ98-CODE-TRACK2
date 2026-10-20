r"""Engine v4, LLM-read evidence at low literature volume (protocol/addendum_engine_v4_llm.md, written before this module ran).

Run (PYTHONPATH=src), in this order:
    python -m track2.engine_v4_llm pairs     candidate pairs and annotated sentences (label-free, offline, counts printed)
    python -m track2.engine_v4_llm build     (reuses the pairs)  Europe PMC verbatim verification, annotation shards + manifest (no LLM is called)
    -- the author dispatches the blind annotators and places out/engine/llm_labels/<annotator>_<split>.json --
    python -m track2.engine_v4_llm dev       agreement + kappa gate, nested grouped CV on DEV, freeze of the refit heads
    python -m track2.engine_v4_llm heldout   scores the sealed HELD-OUT exactly once (lock file) and evaluates the decision rule

Reuses engine_v4 / engine_v4_dilution / engine_v4_ablation / engine_score / verify_literature / engine_extractors functions unchanged; edits nothing.
Label: post hoc development, not a pre-registered primary result. The feature-construction functions have no gold argument.
"""
from __future__ import annotations

import hashlib
import json
import math
import pickle
import random
import re
import sys
import time
import urllib.parse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from track2 import empirical_literature_control as elc
from track2 import engine_benchmark as v1
from track2 import engine_extractors as ex
from track2 import engine_llm_agreement as ag
from track2 import engine_score as es
from track2 import engine_v3 as v3
from track2 import engine_v4 as v4
from track2 import engine_v4_ablation as ab
from track2 import engine_v4_dilution as dl
from track2 import verify_literature as vl

ROOT = v4.ROOT
SEED = 20261015
LABEL = "post hoc development (engine v4, LLM-read evidence), not a pre-registered primary result"
OUT = v4.OUT
PROTOCOL = ROOT / "protocol" / "addendum_engine_v4_llm.md"
INPUT_DIR = OUT / "llm_input"            # shards + manifest: the only files given to the annotators
BUILD_DIR = OUT / "llm_build"            # private: candidate pairs with provenance (contains gold flags), never shipped to annotators
LABELS_DIR = OUT / "llm_labels"          # annotator outputs, placed by the author
CACHE_DIR = v4.CACHE / "v4_llm"
PAIR_INDEX_FILE = CACHE_DIR / "pair_index.json"
SOURCE_TEXTS_FILE = CACHE_DIR / "source_texts.json"
CAND_FILE = BUILD_DIR / "candidate_pairs.json"
OUT_FILE = OUT / "validation_v4_llm.json"
FROZEN_FILE = OUT / "engine_v4_llm_frozen.json"
LOCK = OUT / "heldout_v4_llm_scored.lock"

TOP_K = 30
MAX_PAPERS_PER_PAIR = 3
MAX_SENT_PER_PAPER = 2
SHARD_SIZE = 100
KAPPA_GATE = 0.6
C_GRID = (0.1, 1.0, 10.0)
INNER_FOLDS = 3
BATCH_IDS = 50
ANNOTATORS = ("opus", "sonnet")
ITEM_FIELDS = ("id", "drug", "disease", "sentence", "pmid")      # the only fields an item may carry
SPLIT_PREFIX = {"dev": "dev", "heldout": "ho"}
LLM_NAMES = ["llm_log1p_rel_benefit", "llm_log1p_rel_harm", "llm_log1p_rel_robust"]
PRIOR_COLS = ab.PRIOR_COLS                                         # [4, 6, 7]
LLM_COLS = [8, 9, 10]                                              # appended after the 8 B0 columns
MODEL_COLS = {"PL": PRIOR_COLS + LLM_COLS, "PO": PRIOR_COLS, "B0L": list(range(11))}
MODELS = ("PL", "PO", "B0L", "B0", "PO_frozen")
PAIRS = [("PL", "PO"), ("B0L", "B0"), ("B0L", "PO_frozen"), ("PL", "PO_frozen"), ("B0", "PO_frozen")]
PRIMARY_FRACTIONS = (0.03, 0.01)
VOCAB_DIR = ex.LABELS
VOCAB_REL = ["yes", "no", "unclear"]

UNVALIDATED: Dict[str, Dict[str, Any]] = {
    "top_k_candidates": {"value": TOP_K, "status": "UNVALIDATED", "conservative_default": TOP_K, "note": "budget control, fixed by the task"},
    "papers_per_pair": {"value": MAX_PAPERS_PER_PAIR, "status": "UNVALIDATED", "conservative_default": MAX_PAPERS_PER_PAIR, "note": "budget control, fixed by the task"},
    "sentences_per_paper": {"value": MAX_SENT_PER_PAPER, "status": "UNVALIDATED", "conservative_default": MAX_SENT_PER_PAPER, "note": "budget control, fixed by the task"},
    "relevance_rule": {"value": "both annotators say yes", "status": "UNVALIDATED", "conservative_default": "both annotators say yes", "note": "strict conjunction"},
    "direction_rule": {"value": "2-of-3 vote of opus, sonnet, rule; unresolved excluded", "status": "UNVALIDATED",
                       "conservative_default": "2-of-3 vote of opus, sonnet, rule; unresolved excluded", "note": "reliability flag of engine_llm_agreement"},
    "feature_definitions": {"value": LLM_NAMES, "status": "UNVALIDATED", "conservative_default": LLM_NAMES, "note": "protocol section 3"},
    "c_grid": {"value": list(C_GRID), "status": "UNVALIDATED", "conservative_default": [1.0], "note": "nested grouped CV, inner 3-fold on DEV"},
    "kappa_gate": {"value": KAPPA_GATE, "status": "UNVALIDATED", "conservative_default": KAPPA_GATE, "note": "fixed by the task"},
    "evaluation_universe": {"value": "candidate set P_d", "status": "UNVALIDATED", "conservative_default": "candidate set P_d",
                            "note": "avoids the selection asymmetry between gold and non-gold annotation"},
    "training_features": {"value": "f = 1.0 capped counts", "status": "UNVALIDATED", "conservative_default": "f = 1.0 capped counts", "note": "as the ablation refit"},
    "verification_rate_limit": {"value": "4 requests/s", "status": "UNVALIDATED", "conservative_default": "4 requests/s", "note": "Europe PMC limit given by the task"},
}


# ============================================================================= candidate pairs
def candidate_indices(score: np.ndarray, gold_idx: Iterable[int], top_k: int = TOP_K) -> Tuple[List[int], Dict[int, Dict[str, bool]]]:
    """P_d = the top_k drugs by score (stable order) union the gold indices; returns sorted indices and top-flag provenance."""
    top = [int(i) for i in np.argsort(-np.asarray(score), kind="stable")[:top_k]]
    gold = {int(i) for i in gold_idx}
    cand = sorted(set(top) | gold)
    return cand, {i: {"top": i in top} for i in cand}


# ============================================================================= sentences of a paper and thinning reuse
def paper_drug_sentences(papers: Sequence[dict], symbols: Sequence[str], dre: "re.Pattern[str]", gate: Optional[Set[str]],
                         cand: Set[str], k_sent: int = MAX_SENT_PER_PAPER) -> List[Dict[str, List[str]]]:
    """Per paper (retrieval order): {drug: first k_sent sentences that name the drug}, only for candidate drugs and only for papers that
    contain a disease gene symbol (the v1 paper gate). Text only; no label argument."""
    grx = es.gene_regex(symbols)
    out: List[Dict[str, List[str]]] = []
    for p in papers:
        m: Dict[str, List[str]] = {}
        if grx.search(p["text"]):
            for s in elc.split_sentences(p["text"]):
                if gate is not None and not (set(re.findall(r"\w+", s.lower())) & gate):
                    continue
                for nm in sorted({x.group(0).lower() for x in dre.finditer(s)}):
                    if nm in cand:
                        lst = m.setdefault(nm, [])
                        if len(lst) < k_sent:
                            lst.append(s)
        out.append(m)
    return out


def pair_papers(per_paper: Sequence[Dict[str, List[str]]], kept: Iterable[int], k_papers: int = MAX_PAPERS_PER_PAIR) -> Dict[str, List[int]]:
    """Per drug: the first k_papers kept papers (kept order = retrieval order) that mention it."""
    out: Dict[str, List[int]] = {}
    for i in sorted(kept):
        for nm in per_paper[i]:
            lst = out.setdefault(nm, [])
            if len(lst) < k_papers:
                lst.append(i)
    return out


def union_pair_papers(per_paper: Sequence[Dict[str, List[str]]], rank: int) -> Dict[str, List[int]]:
    """Union over the 21 samples (f = 1.0 once, the other fractions x 5 replicates; the dilution/ablation subsamples) of the papers
    that each sample annotates for a drug. Sentences depend only on (paper, drug), so annotating this union once serves every sample."""
    n = len(per_paper)
    acc: Dict[str, Set[int]] = {}
    for fi, rep in dl.sample_keys():
        kept = dl.subsample_indices(n, dl.FRACTIONS[fi], rank, fi, rep)
        for nm, lst in pair_papers(per_paper, kept).items():
            acc.setdefault(nm, set()).update(lst)
    return {nm: sorted(s) for nm, s in acc.items()}


def pmid_of(pid: str) -> str:
    return pid.split(":", 1)[1] if pid.startswith("MED:") else pid


def pid_of(pmid: str) -> str:
    return pmid if ":" in pmid else f"MED:{pmid}"


def build_pair_index(per_paper: Sequence[Dict[str, List[str]]], papers: Sequence[dict], rank: int) -> Dict[str, List[dict]]:
    """drug -> [{i, pid, sents}] for the union papers; this is everything the feature code needs later (offline)."""
    union = union_pair_papers(per_paper, rank)
    return {nm: [{"i": i, "pid": papers[i]["pid"], "sents": per_paper[i][nm]} for i in lst] for nm, lst in sorted(union.items())}


# ============================================================================= verbatim verification (use 3)
def _source_query(source: str, ids: Sequence[str]) -> str:
    return "(" + " OR ".join(f"EXT_ID:{i}" for i in ids) + f") AND SRC:{source}"


def fetch_source_texts(pids: Sequence[str], getter=None, cache: Optional[Dict[str, Optional[str]]] = None) -> Dict[str, Optional[str]]:
    """Title + abstract exactly as `verify_literature.abstract_text` composes them, fetched from Europe PMC by external id in batches
    (rate-limited by the engine limiter, 4 req/s). A pid not returned is None (= not checkable). Results are cached by the caller."""
    getter = getter or es.default_getter
    cache = {} if cache is None else cache
    todo = sorted({p for p in pids if p not in cache})
    by_src: Dict[str, List[str]] = {}
    for p in todo:
        src, _, ident = p.partition(":")
        by_src.setdefault(src, []).append(ident)
    for src, ids in sorted(by_src.items()):
        for k in range(0, len(ids), BATCH_IDS):
            batch = ids[k:k + BATCH_IDS]
            data = getter(es.EPMC, {"query": _source_query(src, batch), "format": "json", "resultType": "core", "pageSize": 100})
            got = {}
            for it in (data.get("resultList") or {}).get("result") or []:
                got[f"{it.get('source', src)}:{it.get('id') or it.get('pmid') or ''}"] = " ".join([it.get("title") or "", it.get("abstractText") or ""])
            for ident in batch:
                cache[f"{src}:{ident}"] = got.get(f"{src}:{ident}")
    return {p: cache.get(p) for p in pids}


def verify_sentences(items: Sequence[dict], texts: Dict[str, Optional[str]]) -> List[str]:
    """Status per item: verified (the squashed sentence is a substring of the squashed source text), not_found, or not_checkable
    (no source text). Same rule as `verify_literature.check` (squash imported unchanged)."""
    sq: Dict[str, str] = {}
    out = []
    for it in items:
        pid = pid_of(it["pmid"])
        src = texts.get(pid)
        if not src:
            out.append("not_checkable")
            continue
        if pid not in sq:
            sq[pid] = vl.squash(src)
        out.append("verified" if vl.squash(it["sentence"]) in sq[pid] else "not_found")
    return out


def repair_join_period(sentence: str, source_text: Optional[str]) -> str:
    """The engine's retrieval stores f"{title}. {abstract}", so a title sentence carries one period that Europe PMC's text does not have
    (a title ending in "." becomes ".." and a title without a final period gains one). If the sentence is not a substring of the source
    but is once its single last period is removed, return the sentence without it; otherwise return it unchanged. The result is then
    checked by the strict rule of `verify_sentences`; nothing else is normalised."""
    if not source_text:
        return sentence
    sq = vl.squash(source_text)
    if vl.squash(sentence) in sq or not sentence.endswith("."):
        return sentence
    cut = sentence[:-1]
    return cut if vl.squash(cut) in sq else sentence


def repair_pair_index(index_by_d: Dict[str, Any], texts: Dict[str, Optional[str]]) -> int:
    """Apply repair_join_period to every sentence of the pair index in place (so that feature keys and item text stay identical); returns
    the number of sentences changed."""
    n = 0
    for v in index_by_d.values():
        for entries in v["pairs"].values():
            for e in entries:
                new = [repair_join_period(x, texts.get(e["pid"])) for x in e["sents"]]
                n += sum(a != b for a, b in zip(new, e["sents"]))
                e["sents"] = new
    return n


# ============================================================================= items and shards
def make_items(disease: str, pair_index: Dict[str, List[dict]]) -> List[dict]:
    """Annotation items of one disease (internal records: disease, drug, pmid, sentence). Deterministic order."""
    out = []
    for drug in sorted(pair_index):
        for e in pair_index[drug]:
            for s in e["sents"]:
                out.append({"disease": disease, "drug": drug, "pmid": pmid_of(e["pid"]), "sentence": s})
    return out


def dedupe(items: Sequence[dict]) -> List[dict]:
    seen, out = set(), []
    for it in items:
        k = (it["disease"], it["drug"], it["pmid"], it["sentence"])
        if k not in seen:
            seen.add(k)
            out.append(it)
    return out


def shuffle_and_id(items: Sequence[dict], split: str, seed: int = SEED) -> List[dict]:
    """Deterministic shuffle, then opaque sequential ids; items carry exactly ITEM_FIELDS (no gold flag, no score)."""
    lst = sorted(items, key=lambda it: (it["disease"], it["drug"], it["pmid"], it["sentence"]))
    random.Random(seed + (0 if split == "dev" else 1)).shuffle(lst)
    return [{"id": f"{SPLIT_PREFIX[split]}-{n:05d}", "drug": it["drug"], "disease": it["disease"], "sentence": it["sentence"], "pmid": it["pmid"]}
            for n, it in enumerate(lst)]


def shard(items: Sequence[dict], size: int = SHARD_SIZE) -> List[List[dict]]:
    return [list(items[k:k + size]) for k in range(0, len(items), size)]


def estimate_tokens(items: Sequence[dict]) -> Dict[str, int]:
    """Crude estimate: 1 token per 4 characters of the JSON text of the items (input) and ~24 tokens of output per item per task
    (direction + relevance, two tasks per item)."""
    text = json.dumps(list(items), ensure_ascii=False)
    return {"input_tokens_items": math.ceil(len(text) / 4), "output_tokens_two_tasks": 24 * len(items) * 2}


def assert_item_fields(items: Sequence[dict]) -> None:
    for it in items:
        if tuple(it) != ITEM_FIELDS:
            raise AssertionError(f"item fields {tuple(it)} differ from {ITEM_FIELDS}")


# ============================================================================= labels, consensus, features
def load_annotations(path: Path, items: Sequence[dict]) -> Dict[str, Dict[str, Any]]:
    """{id: {direction, relevance, ...}} from an annotator file; ids must be exactly the shard ids and labels in the vocabularies."""
    doc = json.loads(path.read_text(encoding="utf-8"))
    rows = {r["id"]: r for r in doc["labels"]}
    if sorted(rows) != sorted(it["id"] for it in items):
        raise ValueError(f"{path.name}: ids differ from the shard ids")
    for r in rows.values():
        if r["direction"] not in VOCAB_DIR or r["relevance"] not in VOCAB_REL:
            raise ValueError(f"{path.name}: unknown label in {r['id']}")
    return rows


def consensus(rel_a: str, rel_b: str, dir_a: str, dir_b: str, sentence: str, drug: str) -> Dict[str, Any]:
    """relevant = both annotators say yes; direction = label shared by at least 2 of (a, b, rule); reliability = robust/majority/unresolved."""
    rule = ex.rule_direction(sentence, drug)
    votes = [dir_a, dir_b, rule]
    flag = ag.reliability(votes)
    top = max(set(votes), key=lambda v: (votes.count(v), -VOCAB_DIR.index(v)))
    return {"relevant": rel_a == "yes" and rel_b == "yes", "direction": None if flag == "unresolved" else top, "reliability": flag}


def label_table(items: Sequence[dict], ann_a: Dict[str, Dict[str, Any]], ann_b: Dict[str, Dict[str, Any]]) -> Dict[Tuple[str, str, str, str], Dict[str, Any]]:
    """(disease, drug, pmid, sentence) -> consensus labels. Only the sentences that were annotated (= verified) are present."""
    out = {}
    for it in items:
        a, b = ann_a[it["id"]], ann_b[it["id"]]
        out[(it["disease"], it["drug"], it["pmid"], it["sentence"])] = consensus(a["relevance"], b["relevance"], a["direction"], b["direction"], it["sentence"], it["drug"])
    return out


def llm_features(disease: str, pair_index: Dict[str, List[dict]], kept: Set[int], labels: Dict[Tuple[str, str, str, str], Dict[str, Any]],
                 drug_row: Dict[str, int], n_rows: int) -> np.ndarray:
    """(n_rows, 3) LLM block of one disease for one sample: log1p of the counts of relevant-benefit, relevant-harm and robust counted sentences
    among the first MAX_PAPERS_PER_PAIR kept papers (MAX_SENT_PER_PAPER sentences each) of every pair. Rows of unannotated pairs are 0.
    There is no gold argument: the function cannot see gold membership."""
    L = np.zeros((n_rows, 3))
    for drug, entries in pair_index.items():
        r = drug_row.get(drug)
        if r is None:
            continue
        nb = nh = nr = 0
        for e in [e for e in entries if e["i"] in kept][:MAX_PAPERS_PER_PAIR]:
            for s in e["sents"]:
                lab = labels.get((disease, drug, pmid_of(e["pid"]), s))
                if lab is None or not lab["relevant"] or lab["direction"] not in ("benefit", "harm"):
                    continue
                nb += lab["direction"] == "benefit"
                nh += lab["direction"] == "harm"
                nr += lab["reliability"] == "robust"
        L[r] = [math.log1p(nb), math.log1p(nh), math.log1p(nr)]
    return L


# ============================================================================= agreement (quality gate)
def agreement_generic(ref: Sequence[str], hyp: Sequence[str], groups: Sequence[str], n_boot: int = 2000, seed: int = SEED) -> Dict[str, Any]:
    """Accuracy and Cohen's kappa for any label set, with percentile 95 % bootstrap intervals resampling papers (groups)."""
    r, h, g = np.array(ref), np.array(hyp), np.array(groups)
    papers = sorted(set(groups))
    members = [np.where(g == p)[0] for p in papers]
    rng = np.random.default_rng(seed)
    ks, accs = [], []
    for _ in range(n_boot):
        idx = np.concatenate([members[i] for i in rng.integers(0, len(papers), len(papers))])
        ks.append(ex.kappa(r[idx], h[idx]))
        accs.append(float(np.mean(r[idx] == h[idx])))

    def ci(v):
        a = np.array(v, dtype=float)
        a = a[~np.isnan(a)]
        return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))] if len(a) else [None, None]
    return {"n": int(len(r)), "n_papers": len(papers), "accuracy": {"estimate": float(np.mean(r == h)), "ci95": ci(accs)},
            "kappa": {"estimate": float(ex.kappa(r, h)), "ci95": ci(ks)}}


def relevance_gate(kappa_estimate: Optional[float], gate: float = KAPPA_GATE) -> bool:
    """True if the relevance labels are reliable enough to be used (point estimate >= gate); NaN or None fails."""
    return bool(kappa_estimate is not None and kappa_estimate == kappa_estimate and kappa_estimate >= gate)


# ============================================================================= nested grouped CV and heads (candidate-set rows)
def k_group_folds(members: Dict[str, List[str]], k: int, seed: int) -> List[List[str]]:
    keys = sorted(members)
    random.Random(seed).shuffle(keys)
    k = max(2, min(k, len(keys)))
    return [sorted(d for g in keys[i::k] for d in members[g]) for i in range(k)]


def _rows(diseases: Sequence[str], Z: Dict[str, np.ndarray], y: Dict[str, np.ndarray], cls: Dict[str, np.ndarray], bad: Set[str],
          cols: Sequence[int]) -> Tuple[np.ndarray, np.ndarray]:
    """Stack candidate-set rows of the diseases; rows whose ATC class is in `bad` (the held-out fold's gold classes) are masked (v1 rule)."""
    X, yy = [], []
    for d in diseases:
        keep = ~np.isin(cls[d], list(bad)) if bad else np.ones(len(cls[d]), dtype=bool)
        X.append(Z[d][keep][:, cols])
        yy.append(y[d][keep])
    return np.vstack(X), np.concatenate(yy)


def _macro_auc(head: Dict[str, Any], diseases: Sequence[str], Z, y, cols) -> float:
    vals = [es.auc_mw(v4.head_decision(head, Z[d][:, cols]), y[d]) for d in diseases]
    vals = [v for v in vals if v == v]
    return float(np.mean(vals)) if vals else float("nan")


def gold_classes(diseases: Sequence[str], y: Dict[str, np.ndarray], cls: Dict[str, np.ndarray]) -> Set[str]:
    return {c for d in diseases for c in cls[d][y[d] == 1]}


def select_C(train: Sequence[str], members: Dict[str, List[str]], Z, y, cls, cols: Sequence[int], seed: int = SEED + 2,
             grid: Sequence[float] = C_GRID, k: int = INNER_FOLDS) -> Tuple[float, Dict[str, float]]:
    """Inner grouped k-fold on the training diseases: the C with the highest mean macro AUC (ties: the smaller C, the more regularised)."""
    tm = {g: [d for d in ds if d in set(train)] for g, ds in members.items()}
    tm = {g: ds for g, ds in tm.items() if ds}
    folds = k_group_folds(tm, k, seed)
    score = {}
    for C in grid:
        aucs = []
        for fold in folds:
            tr = [d for d in train if d not in fold]
            X, yy = _rows(tr, Z, y, cls, gold_classes(fold, y, cls), cols)
            h = v4.fit_head("logit", X, yy, C)
            aucs.append(_macro_auc(h, fold, Z, y, cols))
        score[str(C)] = float(np.nanmean(aucs))
    best = max(grid, key=lambda C: (round(score[str(C)], 12), -C))
    return float(best), score


def nested_oof_heads(dev: Sequence[str], members: Dict[str, List[str]], Z, y, cls, cols: Sequence[int]) -> Tuple[Dict[str, Dict[str, Any]], List[float]]:
    """Outer grouped 5-fold (v4 folds): per disease, the head of its outer fold, fitted with the C chosen by the inner CV of that fold's training set."""
    heads: Dict[str, Dict[str, Any]] = {}
    chosen = []
    for fold in v4.group_folds(members):
        train = [d for d in dev if d not in fold]
        C, _ = select_C(train, members, Z, y, cls, cols)
        X, yy = _rows(train, Z, y, cls, gold_classes(fold, y, cls), cols)
        h = v4.fit_head("logit", X, yy, C)
        h["C"] = C
        chosen.append(C)
        for d in fold:
            heads[d] = h
    return heads, chosen


def fit_final(dev: Sequence[str], members: Dict[str, List[str]], Z, y, cls, cols: Sequence[int]) -> Dict[str, Any]:
    """Head on all DEV diseases (no masking), C chosen by the inner grouped CV over all DEV."""
    C, scores = select_C(dev, members, Z, y, cls, cols)
    X, yy = _rows(dev, Z, y, cls, set(), cols)
    h = v4.fit_head("logit", X, yy, C)
    h["C"], h["inner_cv_macro_auc"] = C, scores
    return v4.head_params(h, cols)


# ============================================================================= scoring on candidate sets
def score_models_pd(diseases: Sequence[str], samples: Dict[str, Dict[Tuple[int, int], np.ndarray]], n_sub: Dict[str, Dict[Tuple[int, int], int]],
                    scorers: Dict[str, Any], y: Dict[str, np.ndarray], pop: Dict[str, np.ndarray]) -> Dict[str, Any]:
    """Layout of engine_v4_ablation.score_models, restricted to candidate rows. samples[d][(fi, rep)] is the (n_cand, 11) matrix;
    scorers[model](d, Zc) returns scores for those rows."""
    res: Dict[str, Any] = {}
    for d in diseases:
        r: Dict[str, Any] = {"n_sub": {}, "metrics": {m: {} for m in scorers}}
        for fi, rep in dl.sample_keys():
            Zc = samples[d][(fi, rep)]
            r["n_sub"].setdefault(fi, []).append(n_sub[d][(fi, rep)])
            for m, fn in scorers.items():
                r["metrics"][m].setdefault(fi, []).append(dl.metrics_one(fn(d, Zc), y[d]))
        r["popularity"] = dl.metrics_one(pop[d], y[d])
        res[d] = r
    return res


def decide(tab: Dict[str, Any], gate_ok: bool) -> Dict[str, Any]:
    """Protocol section 7: supported only if the lower bound of the PL minus PO paired AUC difference is above 0 at both fractions."""
    if not gate_ok:
        return {"hypothesis": "not testable with these labels (relevance kappa below the gate)", "supported": False}
    per = {}
    for f in PRIMARY_FRACTIONS:
        d = tab["fractions"][str(f)]["paired"]["PL_minus_PO"]["auc"]
        per[str(f)] = {"difference": d["estimate"], "ci95": d["ci95"], "lower_above_0": bool(d["ci95"][0] is not None and d["ci95"][0] > 0)}
    n_ok = sum(v["lower_above_0"] for v in per.values())
    return {"per_fraction": per, "supported": n_ok == len(PRIMARY_FRACTIONS),
            "hypothesis": "supported" if n_ok == len(PRIMARY_FRACTIONS) else ("partial, not supported" if n_ok else "not supported")}


def coverage(L: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    """Share of candidate pairs with at least one counted relevant sentence (all, positives, negatives)."""
    has = L[:, :2].sum(1) > 0
    return {"all": float(has.mean()), "positive": float(has[y == 1].mean()) if (y == 1).any() else float("nan"),
            "negative": float(has[y == 0].mean()) if (y == 0).any() else float("nan")}


# ============================================================================= build (label-free)
_P: Dict[str, Any] = {}


def _init_pair_worker(names: List[str]) -> None:
    _P.update({"dre": elc.drug_regex(names), "gate": es.name_gate(names)})


def _pair_job(args: Tuple[str, Sequence[str], int, List[str]]) -> Tuple[str, Dict[str, Any]]:
    name, symbols, rank, cand = args
    papers = es.fetch_papers(symbols, es.YEAR_PRIMARY, v4.CACHE / "epmc", dl.cache_only_getter)
    per_paper = paper_drug_sentences(papers, symbols, _P["dre"], _P["gate"], set(cand))
    return name, {"n_retrieved": len(papers), "pairs": build_pair_index(per_paper, papers, rank)}


def run_pairs(root: Path = ROOT, workers: int = 6) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Label-free stage: candidate sets P_d (top-30 by frozen B0 union gold) and the annotated sentences of every pair (union over the 21
    samples). Writes the private candidate file and the pair index; prints counts. No network."""
    import pickle as pk
    bench = v4.load_benchmark(root)
    fz = dl._verify(root, bench)
    head = fz["model"]["heads"]["B0"]
    drugs = v1.load_universe(root / "data" / "universe.json")
    names = sorted(v1.universe_index(drugs))
    index = {n: i for i, n in enumerate(names)}
    ranked = sorted(d["name"] for d in bench["diseases"])
    rank = {n: i for i, n in enumerate(ranked)}
    split_of = {d: s for s in ("dev", "heldout") for d in bench["split"][s]}
    b2 = {d["name"]: d for d in json.loads((root / "out" / "engine" / "benchmark_v2.json").read_text(encoding="utf-8"))["diseases"]}
    with v4.FEATURES_FILE.open("rb") as fh:
        f4 = pk.load(fh)
    cands: Dict[str, Any] = {}
    jobs = []
    for dis in bench["diseases"]:
        d = dis["name"]
        sc = v4.head_decision(head, f4[d]["X"][:, head["cols"]])
        g1 = {index[g] for g in dis["gold1"] if g in index}
        g3 = {index[g] for g in dis["gold3"] if g in index}
        g2 = {index[x["name"]] for x in b2[d]["gold2"] if x["name"] in index} if d in b2 else set()
        cand, prov = candidate_indices(sc, g1 | g2 | g3)
        cands[d] = {"split": split_of[d], "gold2_available": d in b2,
                    "candidates": {names[i]: {"top": prov[i]["top"], "g1": i in g1, "g2": i in g2, "g3": i in g3} for i in cand}}
        jobs.append((d, sorted({g["symbol"] for g in dis["genes"]}), rank[d], [names[i] for i in cand]))
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    index_by_d: Dict[str, Any] = {}
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_pair_worker, initargs=(names,)) as pool:
        for i, (name, res) in enumerate(pool.map(_pair_job, jobs, chunksize=1)):
            index_by_d[name] = res
            if (i + 1) % 10 == 0:
                print(json.dumps({"pair_index": i + 1, "of": len(jobs), "seconds": round(time.time() - t0)}), flush=True)
    v4.write_json(PAIR_INDEX_FILE, index_by_d)
    v4.write_json(CAND_FILE, cands)
    n_items = {s: sum(len(make_items(d, index_by_d[d]["pairs"])) for d in index_by_d if split_of[d] == s) for s in ("dev", "heldout")}
    print(json.dumps({"pairs_done": True, "candidate_pairs": {s: sum(len(cands[d]["candidates"]) for d in cands if cands[d]["split"] == s) for s in ("dev", "heldout")},
                      "items_before_dedupe_and_verification": n_items, "seconds": round(time.time() - t0)}), flush=True)
    return bench, cands, index_by_d


def run_build(root: Path = ROOT, workers: int = 6, getter=None) -> dict:
    """Pairs (reused from disk if present), verbatim verification against Europe PMC, shards and manifest. Calls no LLM."""
    if PAIR_INDEX_FILE.exists() and CAND_FILE.exists():
        bench = v4.load_benchmark(root)
        fz = dl._verify(root, bench)
        cands = json.loads(CAND_FILE.read_text(encoding="utf-8"))
        index_by_d = json.loads(PAIR_INDEX_FILE.read_text(encoding="utf-8"))
        split_of = {d: s for s in ("dev", "heldout") for d in bench["split"][s]}
    else:
        bench, cands, index_by_d = run_pairs(root, workers)
        fz = dl._verify(root, bench)
        split_of = {d: s for s in ("dev", "heldout") for d in bench["split"][s]}
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    # source texts from Europe PMC, join-period repair of the title sentences, then items and strict verbatim verification
    pids = sorted({e["pid"] for v in index_by_d.values() for entries in v["pairs"].values() for e in entries})
    cache = json.loads(SOURCE_TEXTS_FILE.read_text(encoding="utf-8")) if SOURCE_TEXTS_FILE.exists() else {}
    t1 = time.time()
    texts = fetch_source_texts(pids, getter, cache)
    v4.write_json(SOURCE_TEXTS_FILE, cache)
    verify_seconds = time.time() - t1
    strict = {s: 0 for s in ("dev", "heldout")}
    for d in index_by_d:
        strict[split_of[d]] += sum(st == "verified" for st in verify_sentences(dedupe(make_items(d, index_by_d[d]["pairs"])), texts))
    n_repaired = repair_pair_index(index_by_d, texts)
    v4.write_json(PAIR_INDEX_FILE, index_by_d)
    internal: Dict[str, List[dict]] = {"dev": [], "heldout": []}
    for d in sorted(index_by_d):
        internal[split_of[d]].extend(make_items(d, index_by_d[d]["pairs"]))
    internal = {s: dedupe(v) for s, v in internal.items()}
    man: Dict[str, Any] = {"label": LABEL, "protocol": str(PROTOCOL.relative_to(root)), "protocol_sha256": v4.sha256_file(PROTOCOL),
                           "guidelines": {"direction": "protocol/annotation_guideline.md", "relevance": "protocol/relevance_guideline.md"},
                           "guideline_sha256": {"direction": v4.sha256_file(root / "protocol" / "annotation_guideline.md"),
                                                "relevance": v4.sha256_file(root / "protocol" / "relevance_guideline.md")},
                           "guideline_tokens": math.ceil((len((root / "protocol" / "annotation_guideline.md").read_text(encoding="utf-8"))
                                                          + len((root / "protocol" / "relevance_guideline.md").read_text(encoding="utf-8"))) / 4),
                           "token_estimate_method": "1 token per 4 characters of item JSON and guidelines; 24 output tokens per item per task; one shard, one annotator, one pass over both tasks",
                           "item_fields": list(ITEM_FIELDS), "shard_size": SHARD_SIZE, "seed": SEED, "top_k": TOP_K,
                           "max_papers_per_pair": MAX_PAPERS_PER_PAIR, "max_sentences_per_paper": MAX_SENT_PER_PAPER,
                           "candidate_file_sha256": v4.sha256_file(CAND_FILE), "frozen_check": fz["check"],
                           "verification": {"source": "Europe PMC (EXT_ID batches), rule of track2.verify_literature.squash", "seconds": round(verify_seconds, 1),
                                            "n_papers_requested": len(pids), "n_papers_with_source_text": sum(1 for p in pids if texts.get(p)),
                                            "join_period_repairs": n_repaired, "verified_strict_before_repair": strict},
                           "splits": {}}
    for old in INPUT_DIR.glob("*_shard_*.json"):
        old.unlink()
    pid_sets: Dict[str, Set[str]] = {}
    for s in ("dev", "heldout"):
        its = internal[s]
        status = verify_sentences(its, texts)
        kept = [it for it, st in zip(its, status) if st == "verified"]
        items = shuffle_and_id(kept, s)
        assert_item_fields(items)
        shards = shard(items)
        files = []
        for n, sh in enumerate(shards, 1):
            fn = INPUT_DIR / f"{s}_shard_{n:03d}.json"
            v4.write_json(fn, {"split": s, "shard": n, "n_items": len(sh), "items": sh})
            tk = estimate_tokens(sh)
            files.append({"file": fn.name, "n_items": len(sh), "sha256": v4.sha256_file(fn), **tk,
                          "total_tokens_one_pass_both_tasks": tk["input_tokens_items"] + man["guideline_tokens"] + tk["output_tokens_two_tasks"]})
        ds = [d for d in cands if cands[d]["split"] == s]
        pid_sets[s] = {pid_of(it["pmid"]) for it in items}
        man["splits"][s] = {
            "n_diseases": len(ds), "n_candidate_pairs": sum(len(cands[d]["candidates"]) for d in ds),
            "n_pairs_with_sentences": sum(len(index_by_d[d]["pairs"]) for d in ds),
            "n_items_before_verification": len(its), "verification_status_counts": {k: status.count(k) for k in ("verified", "not_found", "not_checkable")},
            "n_items": len(items), "n_papers": len(pid_sets[s]), "n_shards": len(shards), "shards": files,
            "tokens_per_shard_mean": float(np.mean([f["total_tokens_one_pass_both_tasks"] for f in files])) if files else 0.0,
            "tokens_per_shard_max": max((f["total_tokens_one_pass_both_tasks"] for f in files), default=0)}
    man["papers_in_both_splits"] = len(pid_sets["dev"] & pid_sets["heldout"])
    man["n_items_total"] = sum(man["splits"][s]["n_items"] for s in man["splits"])
    man["n_shards_total"] = sum(man["splits"][s]["n_shards"] for s in man["splits"])
    v4.write_json(INPUT_DIR / "manifest.json", man)
    print(json.dumps({"items": man["n_items_total"], "shards": man["n_shards_total"],
                      "per_split": {s: [man["splits"][s]["n_items"], man["splits"][s]["n_shards"]] for s in man["splits"]}}), flush=True)
    return man


# ============================================================================= DEV and HELD-OUT stages (run by the author after annotation)
def _load_split_items(split: str) -> List[dict]:
    items: List[dict] = []
    for fn in sorted(INPUT_DIR.glob(f"{split}_shard_*.json")):
        items.extend(json.loads(fn.read_text(encoding="utf-8"))["items"])
    return items


def _annotations(split: str, items: Sequence[dict]) -> Dict[str, Dict[str, Dict[str, Any]]]:
    return {a: load_annotations(LABELS_DIR / f"{a}_{split}.json", items) for a in ANNOTATORS}


def agreement_block(items: Sequence[dict], ann: Dict[str, Dict[str, Dict[str, Any]]]) -> Dict[str, Any]:
    groups = [it["pmid"] for it in items]
    a, b = ann[ANNOTATORS[0]], ann[ANNOTATORS[1]]
    out = {"relevance": agreement_generic([a[i["id"]]["relevance"] for i in items], [b[i["id"]]["relevance"] for i in items], groups),
           "direction": agreement_generic([a[i["id"]]["direction"] for i in items], [b[i["id"]]["direction"] for i in items], groups),
           "label_counts": {n: {"relevance": {c: sum(1 for i in items if ann[n][i["id"]]["relevance"] == c) for c in VOCAB_REL},
                                "direction": {c: sum(1 for i in items if ann[n][i["id"]]["direction"] == c) for c in VOCAB_DIR}} for n in ANNOTATORS}}
    lab = label_table(items, a, b)
    out["reliability_counts"] = {f: sum(1 for v in lab.values() if v["reliability"] == f) for f in ("robust", "majority", "unresolved")}
    return out


def _context(split_diseases: Sequence[str], with_labels: bool) -> Dict[str, Any]:
    """Universe, candidate sets, features at all 21 samples (label-free); gold labels only if with_labels."""
    bench = v4.load_benchmark(ROOT)
    drugs = v1.load_universe(ROOT / "data" / "universe.json")
    universe = v1.universe_index(drugs)
    names = sorted(universe)
    index = {n: i for i, n in enumerate(names)}
    classes = np.array(es.drug_classes(universe, names))
    cands = json.loads(CAND_FILE.read_text(encoding="utf-8"))
    pair_idx = json.loads(PAIR_INDEX_FILE.read_text(encoding="utf-8"))
    with v4.FEATURES_FILE.open("rb") as fh:
        f4 = pickle.load(fh)
    with dl.FEATS_FILE.open("rb") as fh:
        fd = pickle.load(fh)
    by_name = {d["name"]: d for d in bench["diseases"]}
    ctx: Dict[str, Any] = {"bench": bench, "names": names, "index": index, "f4": f4, "fd": fd, "pair_idx": pair_idx, "cands": {}, "cls": {}, "pop": {}, "y": {}, "classes": classes}
    for d in split_diseases:
        c = [index[n] for n in cands[d]["candidates"]]
        ctx["cands"][d] = np.array(sorted(c), dtype=int)
        ctx["cls"][d] = classes[ctx["cands"][d]]
        ctx["pop"][d] = f4[d]["popularity"][ctx["cands"][d]]
        if with_labels:
            yy = np.zeros(len(names))
            for g in by_name[d]["gold1"]:
                if g in index:
                    yy[index[g]] = 1
            ctx["y"][d] = yy
    return ctx


def sample_matrices(ctx: Dict[str, Any], diseases: Sequence[str], labels) -> Tuple[Dict[str, Dict[Tuple[int, int], np.ndarray]], Dict[str, Dict[Tuple[int, int], int]], Dict[str, Dict[Tuple[int, int], np.ndarray]]]:
    """Per disease and sample: the (n_cand, 11) matrix [8 thinned B0 columns | 3 LLM columns], n_sub, and the LLM block (for coverage)."""
    Zs, ns, Ls = {}, {}, {}
    for d in diseases:
        cand = ctx["cands"][d]
        row_of = {ctx["names"][i]: r for r, i in enumerate(cand)}
        n_ret = ctx["pair_idx"][d]["n_retrieved"]
        Zs[d], ns[d], Ls[d] = {}, {}, {}
        for fi, rep in dl.sample_keys():
            kept = set(dl.subsample_indices(n_ret, dl.FRACTIONS[fi], ctx["rank"][d], fi, rep).tolist())
            L = llm_features(d, ctx["pair_idx"][d]["pairs"], kept, labels, row_of, len(cand))
            X = ctx["fd"][d]["samples"][(fi, rep)]["X"][cand]
            Zs[d][(fi, rep)], Ls[d][(fi, rep)] = np.hstack([X, L]), L
            ns[d][(fi, rep)] = ctx["fd"][d]["samples"][(fi, rep)]["n_sub"]
    return Zs, ns, Ls


def _prepare(split: str, diseases: Sequence[str], with_labels: bool) -> Dict[str, Any]:
    items = _load_split_items(split)
    ann = _annotations(split, items)
    ctx = _context(diseases, with_labels)
    ctx["rank"] = {n: i for i, n in enumerate(sorted(d["name"] for d in ctx["bench"]["diseases"]))}
    labels = label_table(items, ann[ANNOTATORS[0]], ann[ANNOTATORS[1]])
    Zs, ns, Ls = sample_matrices(ctx, diseases, labels)
    return {"items": items, "ann": ann, "ctx": ctx, "labels": labels, "Zs": Zs, "ns": ns, "Ls": Ls}


def _scorers(heads: Dict[str, Dict[str, Dict[str, Any]]], b0_head: Dict[str, Any], mean: np.ndarray) -> Dict[str, Any]:
    """heads[model][disease] = fitted head. B0 and PO_frozen use the frozen B0 weights (ablation scorers) on the first 8 columns."""
    fb, fp = ab.frozen_scorer("B0", b0_head, mean), ab.frozen_scorer("PRIOR_ONLY", b0_head, mean)
    sc = {m: (lambda d, Zc, m=m: v4.head_decision(heads[m][d], Zc[:, MODEL_COLS[m]])) for m in ("PL", "PO", "B0L")}
    sc["B0"] = lambda d, Zc: fb(d, Zc[:, :8])
    sc["PO_frozen"] = lambda d, Zc: fp(d, Zc[:, :8])
    return sc


def run_dev(root: Path = ROOT) -> dict:
    t0 = time.time()
    bench = v4.load_benchmark(root)
    dev = list(bench["split"]["dev"])
    P = _prepare("dev", dev, True)
    ctx = P["ctx"]
    agr = agreement_block(P["items"], P["ann"])
    gate = relevance_gate(agr["relevance"]["kappa"]["estimate"])
    out: Dict[str, Any] = {"label": LABEL, "protocol": str(PROTOCOL.relative_to(root)), "protocol_sha256_at_dev_run": v4.sha256_file(PROTOCOL),
                           "unvalidated_registry": UNVALIDATED, "dev": {"agreement": agr, "relevance_gate": {"kappa_point": agr["relevance"]["kappa"]["estimate"], "threshold": KAPPA_GATE, "passed": gate}}}
    if not gate:
        out["status"] = "relevance kappa below the gate: labels unreliable, LLM features not used, no HELD-OUT scoring"
        v4.write_json(OUT_FILE, out)
        print(json.dumps(out["dev"]["relevance_gate"]), flush=True)
        return out
    y = {d: ctx["y"][d][ctx["cands"][d]] for d in dev}
    members = v4.members_of(bench, dev)
    Z1 = {d: P["Zs"][d][(0, 0)] for d in dev}
    heads = {m: nested_oof_heads(dev, members, Z1, y, ctx["cls"], MODEL_COLS[m]) for m in ("PL", "PO", "B0L")}
    fz = dl._verify(root, bench)
    b0 = fz["model"]["heads"]["B0"]
    mean = ab.group_means(ctx["f4"], dev)
    res = score_models_pd(dev, P["Zs"], P["ns"], _scorers({m: h[0] for m, h in heads.items()}, b0, mean), y, ctx["pop"])
    tab = ab.table(res, members, y, MODELS, PAIRS)
    cov = {str(dl.FRACTIONS[fi]): {k: float(np.mean([coverage(P["Ls"][d][(fi, r)], y[d])[k] for d in dev for r in ([0] if fi == 0 else range(dl.N_REPS))]))
                                   for k in ("all", "positive", "negative")} for fi in range(len(dl.FRACTIONS))}
    final = {m: fit_final(dev, members, Z1, y, ctx["cls"], MODEL_COLS[m]) for m in ("PL", "PO", "B0L")}
    frozen = v4.freeze_models({"heads": final, "llm_feature_names": LLM_NAMES, "model_cols": MODEL_COLS},
                              {"sha256_dev_ids": v4.hash_ids(bench["split"]["dev"]), "sha256_heldout_ids": v4.hash_ids(bench["split"]["heldout"]),
                               "candidate_file_sha256": v4.sha256_file(CAND_FILE), "protocol_sha256": v4.sha256_file(PROTOCOL)})
    v4.write_json(FROZEN_FILE, frozen)
    out["dev"].update({"n_diseases": len(dev), "n_groups": len(members), "table_group_bootstrap": tab,
                       "table_disease_level_bootstrap_sensitivity": ab.table(res, members, y, MODELS, PAIRS, disease_level=True),
                       "chosen_C_per_outer_fold": {m: h[1] for m, h in heads.items()}, "coverage_of_llm_features": cov,
                       "frozen_final_heads": {m: {"C": h["C"], "coef": h["coef"], "intercept": h["intercept"]} for m, h in final.items()},
                       "reading": "DEV out-of-fold refit heads (nested grouped CV); evaluation on candidate sets P_d; frozen B0 / PO_frozen are in-sample"})
    out["frozen"] = {"file": str(FROZEN_FILE.relative_to(root)), "sha256_of_model": frozen["sha256_of_model"], "file_sha256": v4.sha256_file(FROZEN_FILE)}
    out["status"] = "DEV part written; HELD-OUT not scored yet"
    out["dev"]["seconds"] = round(time.time() - t0, 1)
    v4.write_json(OUT_FILE, out)
    print(json.dumps({"dev_done": True}), flush=True)
    return out


def run_heldout(root: Path = ROOT) -> dict:
    if LOCK.exists():
        raise RuntimeError("HELD-OUT LLM experiment already scored once (lock file exists); refusing a second scoring")
    out = json.loads(OUT_FILE.read_text(encoding="utf-8"))
    if "table_group_bootstrap" not in out.get("dev", {}):
        raise RuntimeError("run the dev command first (and it must pass the relevance gate)")
    bench = v4.load_benchmark(root)
    ho = list(bench["split"]["heldout"])
    fr = json.loads(FROZEN_FILE.read_text(encoding="utf-8"))
    m = v4.thaw_v4(fr)
    if m["sha256_dev_ids"] != v4.hash_ids(bench["split"]["dev"]) or m["sha256_heldout_ids"] != v4.hash_ids(ho) or fr["model"]["candidate_file_sha256"] != v4.sha256_file(CAND_FILE):
        raise RuntimeError("frozen LLM model or split / candidate hashes differ from the freeze")
    LOCK.write_text(json.dumps({"scored_at": "first and only held-out LLM scoring", "frozen_file_sha256": v4.sha256_file(FROZEN_FILE),
                                "protocol_sha256": v4.sha256_file(PROTOCOL)}) + "\n", encoding="utf-8")
    P = _prepare("heldout", ho, True)                 # HELD-OUT gold is read only here, after the lock exists
    ctx = P["ctx"]
    agr = agreement_block(P["items"], P["ann"])
    y = {d: ctx["y"][d][ctx["cands"][d]] for d in ho}
    members = v4.members_of(bench, ho)
    heads = {mm: {d: m["heads"][mm] for d in ho} for mm in ("PL", "PO", "B0L")}
    fz = dl._verify(root, bench)
    dev = list(bench["split"]["dev"])
    mean = ab.group_means(ctx["f4"], dev)
    res = score_models_pd(ho, P["Zs"], P["ns"], _scorers(heads, fz["model"]["heads"]["B0"], mean), y, ctx["pop"])
    tab = ab.table(res, members, y, MODELS, PAIRS)
    out["heldout"] = {"n_diseases": len(ho), "n_groups": len(members), "agreement_report_only": agr, "table_group_bootstrap": tab,
                      "table_disease_level_bootstrap_sensitivity": ab.table(res, members, y, MODELS, PAIRS, disease_level=True)}
    out["decision_rule_heldout"] = decide(tab, True)
    out["decision_rule_heldout_disease_level_sensitivity"] = decide(out["heldout"]["table_disease_level_bootstrap_sensitivity"], True)
    out["protocol_sha256_at_heldout_run"] = v4.sha256_file(PROTOCOL)
    out["status"] = "complete: DEV and HELD-OUT; HELD-OUT scored once"
    v4.write_json(OUT_FILE, out)
    print(json.dumps({"heldout_done": True, "decision": out["decision_rule_heldout"]["hypothesis"]}), flush=True)
    return out


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = {"pairs": run_pairs, "build": run_build, "dev": run_dev, "heldout": run_heldout}.get(cmd)
    if fn is None:
        raise SystemExit("usage: python -m track2.engine_v4_llm pairs|build|dev|heldout")
    fn()


if __name__ == "__main__":
    main()
