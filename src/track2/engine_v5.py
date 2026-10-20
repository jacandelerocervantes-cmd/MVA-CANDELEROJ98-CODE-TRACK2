r"""Engine v5, per-pair retrieval of drug-and-disease evidence (protocol/addendum_engine_v5.md, written before this module ran).

Run (PYTHONPATH=src, python without -I):
    python -m track2.engine_v5 retrieve    one cached Europe PMC query per candidate (disease, drug) pair; records the queries, hit counts, record ids
    python -m track2.engine_v5 build       sentences, verbatim verification, shards of <= 100 items and the manifest (no LLM is called)

Label-free: no annotation, no score, no model. P_d is the v4 candidate file (gold decides which pairs enter, as in v4); gold is never used in a
query, an item or a feature. Reuses engine_v4 / engine_v4_llm / engine_v4_llm_audit functions unchanged; edits nothing.
"""
from __future__ import annotations

import datetime
import json
import math
import random
import re
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from track2 import empirical_literature_control as elc
from track2 import engine_benchmark as v1
from track2 import engine_score as es
from track2 import engine_v4 as v4
from track2 import engine_v4_dilution as dl
from track2 import engine_v4_llm as lm
from track2 import engine_v4_llm_audit as au

ROOT = v4.ROOT
SEED = 20261015
LABEL = "post hoc development (engine v5, per-pair retrieval), extends the pre-registered v4 programme; not a pre-registered primary result"
OUT = v4.OUT
PROTOCOL = ROOT / "protocol" / "addendum_engine_v5.md"
DIR_GUIDELINE = ROOT / "protocol" / "annotation_guideline.md"
REL_GUIDELINE = ROOT / "protocol" / "relevance_guideline_v2.md"
GUIDELINE_SHA256 = {"direction": "f56abc1ff93ec8e2a06d8eabb062c318e51bba53cd1c0e738a6a6092fea2e45d",
                    "relevance_v2": "5f55313ea0046b795f32e5fd28baa7330b249c8c9b16c8ce3e51d452808be5d9"}   # frozen in the protocol
CACHE_DIR = v4.CACHE / "v5"
QUERY_CACHE = CACHE_DIR / "query"
SOURCE_TEXTS_FILE = CACHE_DIR / "source_texts.json"
DRUG_TERMS_FILE = CACHE_DIR / "drug_terms.json"
PAIR_INDEX_FILE = CACHE_DIR / "pair_index.json"
BUILD_DIR = OUT / "llm_build_v5"                  # private (retrieval index, provenance); never given to annotators
RETRIEVAL_FILE = BUILD_DIR / "retrieval.json"
INPUT_DIR = OUT / "llm_input_v5"
CAND_FILE = lm.CAND_FILE                          # the v4 candidate sets P_d (reused, hashed in the manifest)

K_RECORDS = 15
MAX_DRUG_TERMS_QUERY = 8
MAX_DIS_TERMS_QUERY = au.MAX_DIS_TERMS            # 6
MAX_SENT_PER_RECORD = 3
MIN_TERM_LEN = 4
SCARCITY_MAX_RECORDS = 2
SHARD_SIZE = 100
FRACTIONS = tuple(dl.FRACTIONS)
N_REPS = dl.N_REPS
ITEM_FIELDS = lm.ITEM_FIELDS
SPLIT_PREFIX = {"dev": "dev5", "heldout": "ho5"}
SRC_PRIORITY = ["salt_stripped", "hyphen_space", "chembl_synonyms", "ctgov_other_names"]


def _uv(value: Any, default: Any, note: str) -> Dict[str, Any]:
    return {"value": value, "status": "UNVALIDATED", "conservative_default": default, "note": note}


UNVALIDATED: Dict[str, Dict[str, Any]] = {
    "k_records_per_pair": _uv(K_RECORDS, 10, "records retrieved per pair; 10 = the audit's direct-query probe; K is never changed after a result"),
    "drug_terms_in_query": _uv(MAX_DRUG_TERMS_QUERY, 1, "name first, then salt-stripped, hyphen/space, ChEMBL, ClinicalTrials.gov terms; matching uses all terms"),
    "disease_terms_in_query": _uv(MAX_DIS_TERMS_QUERY, 1, "name plus the shortest name_set synonyms"),
    "sentences_per_record": _uv(MAX_SENT_PER_RECORD, MAX_SENT_PER_RECORD, "first sentences that name the drug"),
    "min_term_length": _uv(MIN_TERM_LEN, MIN_TERM_LEN, "as the v1 drug regex"),
    "bernoulli_thinning": _uv("each record kept with probability f", "each record kept with probability f", "keeps f = 0.03 and 0.01 distinct"),
    "natural_scarcity_max_records": _uv(SCARCITY_MAX_RECORDS, SCARCITY_MAX_RECORDS, "stratum of pairs with 0-2 direct records at f = 1.0"),
    "top_k_candidates": _uv(lm.TOP_K, lm.TOP_K, "inherited from v4"),
    "relevance_rule": _uv("both annotators say yes", "both annotators say yes", "inherited from v4"),
    "direction_rule": _uv("2-of-3 vote of opus, sonnet, rule; unresolved excluded", "2-of-3 vote of opus, sonnet, rule; unresolved excluded", "inherited from v4"),
    "c_grid": _uv(list(lm.C_GRID), [1.0], "nested grouped CV, inner 3-fold on DEV"),
    "kappa_gate": _uv(lm.KAPPA_GATE, lm.KAPPA_GATE, "DEV relevance kappa point estimate"),
    "evaluation_universe": _uv("candidate set P_d", "candidate set P_d", "as v4 LLM"),
    "requests_per_second": _uv(4, 4, "Europe PMC limit given by the task"),
}


# ============================================================================= drug and disease terms
def ordered_drug_terms(name: str, syn: Dict[str, Set[str]]) -> List[str]:
    """[name] then the synonym terms (>= MIN_TERM_LEN characters, no quotes) by source priority; inside a source: no digit first, shorter, alphabetical.
    Each term appears once (lower case). `syn` is {term: {sources}} for this drug."""
    low = name.lower().strip()
    out: List[str] = [low] if len(low) >= MIN_TERM_LEN else []
    seen = set(out)
    for src in SRC_PRIORITY:
        cand = sorted({t for t, ss in syn.items() if src in ss and len(t) >= MIN_TERM_LEN and '"' not in t},
                      key=lambda t: (any(c.isdigit() for c in t), len(t), t))
        # a term already taken from a higher-priority source is skipped
        for t in cand:
            if t not in seen:
                seen.add(t)
                out.append(t)
    return out


def build_query(drug_terms: Sequence[str], dis_terms: Sequence[str]) -> str:
    """(drug clause) AND (disease clause) AND HAS_ABSTRACT:y; no year limit (protocol section 1)."""
    return f"{au.clause(drug_terms)} AND {au.clause(dis_terms)} AND HAS_ABSTRACT:y"


def search_pair(query: str, getter: Optional[Callable[[str, dict], dict]] = None, cache_dir: Path = QUERY_CACHE, k: int = K_RECORDS) -> Dict[str, Any]:
    """One cached Europe PMC call (relevance order, cursorMark '*', pageSize k), via the audit's cached search."""
    return au.epmc_search(query, k, "*", getter, cache_dir)


# ============================================================================= sentences, thinning
def record_sentences(text: str, rx: "re.Pattern[str]", k: int = MAX_SENT_PER_RECORD) -> List[str]:
    """The first k sentences of the record text that name the drug by any of its terms (no label or gold argument)."""
    out = []
    for s in elc.split_sentences(text):
        if rx.search(s):
            out.append(s)
            if len(out) >= k:
                break
    return out


def thin_records(n_records: int, f: float, disease_rank: int, drug_index: int, frac_index: int, rep: int, seed: int = SEED) -> List[int]:
    """Indices (retrieval order) of the records kept by one sample: each record independently with probability f (f = 1.0 keeps all)."""
    if f >= 1.0:
        return list(range(n_records))
    rng = np.random.default_rng([seed, disease_rank, drug_index, frac_index, rep])
    return [i for i, u in enumerate(rng.random(n_records)) if u < f]


# ============================================================================= items
def make_items(disease: str, pairs: Dict[str, Dict[str, Any]]) -> List[dict]:
    out = []
    for drug in sorted(pairs):
        for rec in pairs[drug]["records"]:
            for s in rec["sents"]:
                out.append({"disease": disease, "drug": drug, "pmid": lm.pmid_of(rec["pid"]), "sentence": s})
    return out


def shuffle_and_id(items: Sequence[dict], split: str, seed: int = SEED + 5) -> List[dict]:
    """Deterministic shuffle, opaque sequential ids (dev5- / ho5-); items carry exactly ITEM_FIELDS."""
    lst = sorted(items, key=lambda it: (it["disease"], it["drug"], it["pmid"], it["sentence"]))
    random.Random(seed + (0 if split == "dev" else 1)).shuffle(lst)
    return [{"id": f"{SPLIT_PREFIX[split]}-{n:05d}", "drug": it["drug"], "disease": it["disease"], "sentence": it["sentence"], "pmid": it["pmid"]}
            for n, it in enumerate(lst)]


def pair_stats(pairs_by_disease: Dict[str, Dict[str, Dict[str, Any]]], diseases: Sequence[str], gold_of: Optional[Dict[Tuple[str, str], Dict[str, bool]]] = None) -> Dict[str, Any]:
    """Counts over the pairs of the given diseases: pairs, pairs with >= 1 direct record, records, verified-sentence pairs; by gold membership if gold_of is given."""
    def agg(keys: List[Tuple[str, str]]) -> Dict[str, Any]:
        recs = [len(pairs_by_disease[d][g]["records_all"]) for d, g in keys]
        sent = [sum(len(r["sents"]) for r in pairs_by_disease[d][g]["records"]) for d, g in keys]
        n = len(keys)
        return {"n_pairs": n, "n_pairs_with_ge1_direct_record": sum(r >= 1 for r in recs), "n_direct_records": sum(recs),
                "n_pairs_with_ge1_verified_sentence": sum(s >= 1 for s in sent), "n_verified_sentences": sum(sent),
                "n_pairs_scarce_0_to_2_records": sum(r <= SCARCITY_MAX_RECORDS for r in recs),
                "share_with_ge1_direct_record": (sum(r >= 1 for r in recs) / n) if n else None}
    keys = [(d, g) for d in diseases for g in pairs_by_disease[d]]
    out: Dict[str, Any] = {"all": agg(keys)}
    if gold_of is not None:
        out["gold1"] = agg([k for k in keys if gold_of[k]["g1"]])
        out["not_gold1"] = agg([k for k in keys if not gold_of[k]["g1"]])
        out["any_gold"] = agg([k for k in keys if gold_of[k]["g1"] or gold_of[k]["g2"] or gold_of[k]["g3"]])
        out["no_gold"] = agg([k for k in keys if not (gold_of[k]["g1"] or gold_of[k]["g2"] or gold_of[k]["g3"])])
    return out


def assert_item_fields(items: Sequence[dict]) -> None:
    lm.assert_item_fields(items)


# ============================================================================= retrieve
def _load_synonyms(names: Sequence[str]) -> Dict[str, Dict[str, Set[str]]]:
    return au.build_synonyms(names, sorted((v4.CACHE / "ctgov").glob("*.json")) + sorted((v4.CACHE / "ctgov_v2").glob("*.json")),
                             sorted((ROOT / "data" / "cache").glob("*.json")))


def _context(root: Path = ROOT) -> Dict[str, Any]:
    bench = v4.load_benchmark(root)
    drugs = v1.load_universe(root / "data" / "universe.json")
    names = sorted(v1.universe_index(drugs))
    cands = json.loads(CAND_FILE.read_text(encoding="utf-8"))
    rank = {n: i for i, n in enumerate(sorted(d["name"] for d in bench["diseases"]))}
    dis = {d["name"]: list(d["identifiers"].get("name_set", [])) for d in bench["diseases"]}
    return {"bench": bench, "names": names, "index": {n: i for i, n in enumerate(names)}, "cands": cands, "rank": rank, "dis": dis}


def run_retrieve(root: Path = ROOT, workers: int = 4, getter=None) -> Dict[str, Any]:
    t0 = time.time()
    ctx = _context(root)
    names, cands, dis = ctx["names"], ctx["cands"], ctx["dis"]
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    syn = _load_synonyms(names)
    print(json.dumps({"synonyms_built": len(syn), "seconds": round(time.time() - t0)}), flush=True)
    all_terms: Dict[str, List[str]] = {n: ordered_drug_terms(n, syn.get(n.lower(), {})) for n in names}
    v4.write_json(DRUG_TERMS_FILE, {"terms": all_terms, "sources": {n: {t: sorted(s) for t, s in syn.get(n.lower(), {}).items()} for n in names if n.lower() in syn}})
    jobs = []
    for d in sorted(cands):
        dt = au.query_disease_terms(d, dis[d], MAX_DIS_TERMS_QUERY)
        for g in sorted(cands[d]["candidates"]):
            jobs.append((d, g, dt))
    todo = [j for j in jobs if all_terms.get(j[1])]
    print(json.dumps({"pairs": len(jobs), "pairs_with_drug_terms": len(todo)}), flush=True)

    def one(job: Tuple[str, str, List[str]]) -> Tuple[str, str, Dict[str, Any]]:
        d, g, dt = job
        q = build_query(all_terms[g][:MAX_DRUG_TERMS_QUERY], dt)
        page = search_pair(q, getter)
        recs = [r for r in page["records"] if r["abstract"]]
        seen, pids = set(), []
        for r in recs:
            if r["pid"] not in seen:
                seen.add(r["pid"])
                pids.append(r["pid"])
        return d, g, {"query": q, "drug_terms_query": all_terms[g][:MAX_DRUG_TERMS_QUERY], "n_drug_terms_total": len(all_terms[g]), "disease_terms": dt,
                      "hit_count": page["hitCount"], "n_records": len(pids), "pids": pids}

    result: Dict[str, Any] = {d: {"split": cands[d]["split"], "pairs": {}} for d in cands}
    for d, g, dt in jobs:
        if not all_terms.get(g):
            result[d]["pairs"][g] = {"query": None, "drug_terms_query": [], "n_drug_terms_total": 0, "disease_terms": dt, "hit_count": 0, "n_records": 0, "pids": []}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, (d, g, r) in enumerate(pool.map(one, todo)):
            result[d]["pairs"][g] = r
            if (i + 1) % 250 == 0:
                print(json.dumps({"queries_done": i + 1, "of": len(todo), "seconds": round(time.time() - t0)}), flush=True)
    doc = {"label": LABEL, "protocol": str(PROTOCOL.relative_to(root)), "retrieved_on": datetime.date.today().isoformat(), "k_records": K_RECORDS,
           "candidate_file_sha256": v4.sha256_file(CAND_FILE), "diseases": result}
    v4.write_json(RETRIEVAL_FILE, doc)
    n_by = {s: {"pairs": 0, "pairs_ge1_record": 0, "records": 0} for s in ("dev", "heldout")}
    for d, v in result.items():
        for r in v["pairs"].values():
            n = n_by[v["split"]]
            n["pairs"] += 1
            n["pairs_ge1_record"] += r["n_records"] >= 1
            n["records"] += r["n_records"]
    print(json.dumps({"retrieve_done": True, "per_split": n_by, "seconds": round(time.time() - t0)}), flush=True)
    return doc


# ============================================================================= build
def build_pair_index(doc: Dict[str, Any], all_terms: Dict[str, List[str]], getter=None) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """disease -> drug -> {records_all: [pid in order], records: [{pid, sents}]}; sentences from the cached query responses (a cache miss is an error)."""
    miss = getter or dl.cache_only_getter
    out: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for d in sorted(doc["diseases"]):
        out[d] = {}
        for g, r in sorted(doc["diseases"][d]["pairs"].items()):
            if not r["query"]:
                out[d][g] = {"records_all": [], "records": []}
                continue
            page = search_pair(r["query"], miss)
            rx = elc.term_regex(all_terms[g])
            by_pid = {x["pid"]: x for x in page["records"] if x["abstract"]}
            recs = []
            for pid in r["pids"]:
                x = by_pid[pid]
                sents = record_sentences(f"{x['title']}. {x['abstract']}", rx)
                if sents:
                    recs.append({"pid": pid, "sents": sents})
            out[d][g] = {"records_all": list(r["pids"]), "records": recs}
    return out


def run_build(root: Path = ROOT, getter=None) -> Dict[str, Any]:
    t0 = time.time()
    doc = json.loads(RETRIEVAL_FILE.read_text(encoding="utf-8"))
    if doc["candidate_file_sha256"] != v4.sha256_file(CAND_FILE):
        raise RuntimeError("candidate file differs from the one used by the retrieval")
    all_terms = json.loads(DRUG_TERMS_FILE.read_text(encoding="utf-8"))["terms"]
    split_of = {d: v["split"] for d, v in doc["diseases"].items()}
    cands = json.loads(CAND_FILE.read_text(encoding="utf-8"))
    idx = build_pair_index(doc, all_terms)
    n_before = {s: sum(len(make_items(d, idx[d])) for d in idx if split_of[d] == s) for s in ("dev", "heldout")}
    pids = sorted({r["pid"] for d in idx for p in idx[d].values() for r in p["records"]})
    cache = json.loads(SOURCE_TEXTS_FILE.read_text(encoding="utf-8")) if SOURCE_TEXTS_FILE.exists() else {}
    t1 = time.time()
    texts = lm.fetch_source_texts(pids, getter, cache)
    v4.write_json(SOURCE_TEXTS_FILE, cache)
    verify_seconds = time.time() - t1
    # join-period repair, then strict verification; unverified sentences are removed from the index
    n_repaired = 0
    for d in idx:
        for p in idx[d].values():
            for r in p["records"]:
                new = [lm.repair_join_period(s, texts.get(r["pid"])) for s in r["sents"]]
                n_repaired += sum(a != b for a, b in zip(new, r["sents"]))
                r["sents"] = new
    status_counts = {s: Counter() for s in ("dev", "heldout")}
    for d in idx:
        for g, p in idx[d].items():
            kept = []
            for r in p["records"]:
                st = lm.verify_sentences([{"pmid": lm.pmid_of(r["pid"]), "sentence": s} for s in r["sents"]], texts)
                status_counts[split_of[d]].update(st)
                ok = [s for s, x in zip(r["sents"], st) if x == "verified"]
                if ok:
                    kept.append({"pid": r["pid"], "sents": ok})
            p["records"] = kept
    v4.write_json(PAIR_INDEX_FILE, idx)
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    for old in INPUT_DIR.glob("*_shard_*.json"):
        old.unlink()
    sha = {"direction": v4.sha256_file(DIR_GUIDELINE), "relevance_v2": v4.sha256_file(REL_GUIDELINE)}
    if sha != GUIDELINE_SHA256:
        raise RuntimeError(f"guideline hashes differ from the protocol: {sha}")
    guideline_tokens = math.ceil((len(DIR_GUIDELINE.read_text(encoding="utf-8")) + len(REL_GUIDELINE.read_text(encoding="utf-8"))) / 4)
    gold_of = {(d, g): {k: bool(v[k]) for k in ("g1", "g2", "g3")} for d in cands if cands[d]["split"] == "dev" for g, v in cands[d]["candidates"].items()}
    man: Dict[str, Any] = {
        "label": LABEL, "protocol": str(PROTOCOL.relative_to(root)), "protocol_sha256": v4.sha256_file(PROTOCOL),
        "guidelines": {"direction": "protocol/annotation_guideline.md", "relevance": "protocol/relevance_guideline_v2.md"}, "guideline_sha256": sha,
        "guideline_tokens": guideline_tokens,
        "token_estimate_method": "1 token per 4 characters of item JSON and guidelines; 24 output tokens per item per task; one shard, one annotator, one pass over both tasks",
        "item_fields": list(ITEM_FIELDS), "shard_size": SHARD_SIZE, "seed": SEED, "k_records": K_RECORDS, "max_sentences_per_record": MAX_SENT_PER_RECORD,
        "candidate_file_sha256": v4.sha256_file(CAND_FILE), "retrieval_file_sha256": v4.sha256_file(RETRIEVAL_FILE), "retrieved_on": doc["retrieved_on"],
        "unvalidated": UNVALIDATED,
        "verification": {"source": "Europe PMC (EXT_ID batches), rule of track2.verify_literature.squash via engine_v4_llm", "seconds": round(verify_seconds, 1),
                         "n_papers_requested": len(pids), "n_papers_with_source_text": sum(1 for p in pids if texts.get(p)), "join_period_repairs": n_repaired},
        "splits": {}}
    pid_sets: Dict[str, Set[str]] = {}
    for s in ("dev", "heldout"):
        ds = [d for d in sorted(idx) if split_of[d] == s]
        its = lm.dedupe([it for d in ds for it in make_items(d, idx[d])])
        items = shuffle_and_id(its, s)
        assert_item_fields(items)
        shards = lm.shard(items, SHARD_SIZE)
        files = []
        for n, sh in enumerate(shards, 1):
            fn = INPUT_DIR / f"{s}_shard_{n:03d}.json"
            v4.write_json(fn, {"split": s, "shard": n, "n_items": len(sh), "items": sh})
            tk = lm.estimate_tokens(sh)
            files.append({"file": fn.name, "n_items": len(sh), "sha256": v4.sha256_file(fn), **tk,
                          "total_tokens_one_pass_both_tasks": tk["input_tokens_items"] + guideline_tokens + tk["output_tokens_two_tasks"]})
        pid_sets[s] = {lm.pid_of(it["pmid"]) for it in items}
        st = pair_stats(idx, ds, gold_of if s == "dev" else None)
        tot = sum(f["total_tokens_one_pass_both_tasks"] for f in files)
        man["splits"][s] = {"n_diseases": len(ds), "pair_stats": st, "n_items_before_verification": n_before[s], "n_items_after_dedupe_and_verification": len(items),
                            "verification_status_counts": {k: status_counts[s][k] for k in ("verified", "not_found", "not_checkable")},
                            "n_items": len(items), "n_papers": len(pid_sets[s]), "n_shards": len(shards), "shards": files,
                            "tokens_one_annotator": tot, "tokens_two_annotators": 2 * tot,
                            "tokens_per_shard_mean": float(np.mean([f["total_tokens_one_pass_both_tasks"] for f in files])) if files else 0.0,
                            "tokens_per_shard_max": max((f["total_tokens_one_pass_both_tasks"] for f in files), default=0)}
    man["papers_in_both_splits"] = len(pid_sets["dev"] & pid_sets["heldout"])
    man["n_items_total"] = sum(man["splits"][s]["n_items"] for s in man["splits"])
    man["n_shards_total"] = sum(man["splits"][s]["n_shards"] for s in man["splits"])
    man["tokens_two_annotators_total"] = sum(man["splits"][s]["tokens_two_annotators"] for s in man["splits"])
    man["seconds"] = round(time.time() - t0, 1)
    v4.write_json(INPUT_DIR / "manifest.json", man)
    print(json.dumps({"build_done": True, "items": man["n_items_total"], "shards": man["n_shards_total"],
                      "per_split": {s: [man["splits"][s]["n_items"], man["splits"][s]["n_shards"]] for s in man["splits"]}}), flush=True)
    return man


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = {"retrieve": run_retrieve, "build": run_build}.get(cmd)
    if fn is None:
        raise SystemExit("usage: python -m track2.engine_v5 retrieve|build")
    fn()


if __name__ == "__main__":
    main()
