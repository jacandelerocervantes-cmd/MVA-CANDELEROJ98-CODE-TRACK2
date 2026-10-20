r"""Evidence engine v4, phase 1 - scale-up, literature-volume strata, embeddings feature, abstention
(protocol/addendum_engine_v4.md, written before this module ran).

Run (PYTHONPATH=src), in this order:
    python -m track2.engine_v4 benchmark   eligibility + groups + split -> out/engine/benchmark_v4.json (no literature, no score)
    python -m track2.engine_v4 download    the two authorised Hugging Face models -> data/models/ (git-ignored), revisions recorded
    python -m track2.engine_v4 features    Europe PMC retrieval + v1 features + embedding texts of every disease -> data/cache_engine/
    python -m track2.engine_v4 embed       frozen BiomedBERT / SapBERT vectors of the texts (isolated `python -I` subprocess)
    python -m track2.engine_v4 dev         DEV only: out-of-fold B0 / E / M, 1-SE inclusion, conformal, freeze, validation_v4.json
    python -m track2.engine_v4 heldout     scores the sealed HELD-OUT exactly once (lock file) and completes validation_v4.json

Reuses engine_score (v1), engine_benchmark_v2 and engine_v3 functions unchanged. Label: post hoc development, not a
pre-registered primary result. engine_score / engine_v3 / engine_benchmark_v2 are imported, never edited.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier

from track2 import empirical_literature_control as elc
from track2 import engine_benchmark as v1
from track2 import engine_benchmark_v2 as b2
from track2 import engine_score as es
from track2 import engine_score_v2 as v2s
from track2 import engine_v3 as v3

ROOT = Path(__file__).resolve().parents[2]
SEED = es.SEED
LABEL = "post hoc development (engine v4, phase 1), not a pre-registered primary result"
OUT = ROOT / "out" / "engine"
CACHE = es.CACHE
FEATURES_FILE = CACHE / "v4_features.pkl"
HITS_FILE = CACHE / "v4_hitcount.json"
EMB_DIR = CACHE / "v4_embeddings"
MODELS_DIR = ROOT / "data" / "models"
PROTOCOL = ROOT / "protocol" / "addendum_engine_v4.md"
ORPHA_PREV_URL = "https://www.orphadata.com/data/xml/en_product9_prev.xml"
ALPHA = 0.1
MIN_GENES_V4, MIN_GOLD_V4 = 1, 1
N_SENT, MAX_TOK = 3, 128
STRATA = ("L", "M", "H")
MIN_GROUPS_STRATUM = 5
MODEL_SPECS = {
    "biomedbert": {"repo": "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext", "pool": "mean",
                   "files": ["config.json", "tokenizer_config.json", "vocab.txt", "pytorch_model.bin"]},
    "sapbert": {"repo": "cambridgeltl/SapBERT-from-PubMedBERT-fulltext", "pool": "cls",
                "files": ["config.json", "tokenizer_config.json", "special_tokens_map.json", "vocab.txt", "model.safetensors"]},
}
B0_COLS = list(range(8))
E_COLS = list(range(11))                       # 8 v1 features + 3 embedding features
EMB_NAMES = ["emb_cos_biomedbert", "emb_cos_sapbert", "emb_has_text"]
MLP_PARAMS = {"hidden": 16, "alpha": 0.01, "max_iter": 100}

UNVALIDATED: Dict[str, Dict[str, Any]] = {
    "conformal_alpha": {"value": ALPHA, "status": "UNVALIDATED", "conservative_default": ALPHA, "note": "fixed by the task"},
    "stratum_cut_points": {"value": [50, 500], "status": "UNVALIDATED", "conservative_default": [50, 500], "note": "fixed by the task"},
    "embedding_sentences": {"value": N_SENT, "status": "UNVALIDATED", "conservative_default": N_SENT, "note": "not tuned"},
    "embedding_max_tokens": {"value": MAX_TOK, "status": "UNVALIDATED", "conservative_default": MAX_TOK, "note": "not tuned"},
    "mlp_hyperparameters": {"value": MLP_PARAMS, "status": "UNVALIDATED", "conservative_default": MLP_PARAMS,
                            "note": "ablation only, never tuned, never selectable"},
    "prevalence_tier_rule": {"value": "Cases/families value, else class >= 1-9 / 1 000 000 -> >5, else unknown",
                             "status": "UNVALIDATED", "conservative_default": "unknown", "note": "rule fixed in protocol section 6"},
    "inclusion_margin_1se": {"value": 1.0, "status": "UNVALIDATED", "conservative_default": 1.0, "note": "fixed by the task"},
    "coverage_validity_margin": {"value": "0.90 minus bootstrap SE", "status": "UNVALIDATED", "conservative_default": "0.90",
                                 "note": "validity margin for conformal inclusion"},
}


def resolve(p: Dict[str, Any], experimental: bool = False) -> Any:
    """unvalidated-gate: an UNVALIDATED value resolves to its conservative default unless experimental mode is on."""
    return p["value"] if p["status"] == "VALIDATED" or experimental else p["conservative_default"]


def check_registry(reg: Dict[str, Dict[str, Any]]) -> List[str]:
    bad = []
    for k, p in reg.items():
        if p.get("status") not in ("VALIDATED", "UNVALIDATED"):
            bad.append(k + ": no status")
        elif p["status"] == "UNVALIDATED" and "conservative_default" not in p:
            bad.append(k + ": no conservative_default")
        elif p["status"] == "VALIDATED" and not p.get("evidence"):
            bad.append(k + ": VALIDATED without evidence")
    return bad


def hash_ids(ids: Sequence[str]) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, ids))).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, sort_keys=True, ensure_ascii=False, default=v3._default) + "\n", encoding="utf-8", newline="\n")


def append_protocol(text: str) -> None:
    with PROTOCOL.open("a", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


# ============================================================================= benchmark v4
def evaluate_candidate_v4(name: str, gold1_names: Set[str], match: Optional[Dict[str, Any]], uidx: Dict[str, dict],
                          curated: bool = False) -> Tuple[bool, str]:
    """The v2 checks in the v2 order with the declared changes D1 (>= 1 gene) and D2 (>= 1 gold 1 drug)."""
    if v1.APPLICATION_CASES.search(name):
        return False, "application case (MVA / Fanconi anaemia) excluded by design"
    if match is None or match["n_matched"] == 0:
        return False, "no HPO or Orphanet disorder matched by identifier, name or synonym"
    if match["n_matched"] > v1.MAX_MATCHED:
        return False, f"umbrella term: {match['n_matched']} matched disorders (> {v1.MAX_MATCHED})"
    if not match["rare_genetic_identity"] and not (curated and match["orphanet_catalogued"] and match["mendelian"]):
        return False, "not a rare genetic disease (no Orphanet-catalogued identity disorder with Mendelian/mitochondrial/sex-linked inheritance)"
    if len(match["genes"]) < MIN_GENES_V4:
        return False, f"fewer than {MIN_GENES_V4} genes ({len(match['genes'])})"
    gold = {g for g in gold1_names if g in uidx}
    if len(gold) < MIN_GOLD_V4:
        return False, f"fewer than {MIN_GOLD_V4} approved universe drugs indicated ({len(gold)})"
    bad = sorted(g for g in gold if uidx[g].get("molecule_type") in v1.BIOLOGIC_THERAPY_TYPES)
    if bad:
        return False, "gold drug is an enzyme replacement / gene / cell therapy"
    return True, "eligible"


def group_key(gold_names: Sequence[str]) -> str:
    return "|".join(sorted(gold_names))


def group_id(gold_names: Sequence[str]) -> str:
    return "G" + hashlib.sha1(group_key(gold_names).encode("utf-8")).hexdigest()[:10]


def make_split_v4(diseases: Sequence[dict], seed: int = SEED) -> Dict[str, Any]:
    """Protocol 1.3. diseases: dicts with name, origin ('v2' | 'v4_new'), group (gold 1 set key).

    DEV = v2 diseases + new diseases whose group contains a v2 disease + groups of the free pool taken in seeded order while
    the number of free diseases in DEV is below floor(n_free / 2); HELD-OUT = the other free groups."""
    by_group: Dict[str, List[str]] = defaultdict(list)
    for d in diseases:
        by_group[d["group"]].append(d["name"])
    v2_groups = {d["group"] for d in diseases if d["origin"] == "v2"}
    dev = sorted(n for g in v2_groups for n in by_group[g])
    free = sorted(g for g in by_group if g not in v2_groups)
    n_free = sum(len(by_group[g]) for g in free)
    order = list(free)
    random.Random(seed).shuffle(order)
    target = n_free // 2
    taken = 0
    held: List[str] = []
    for g in order:
        if taken < target:
            dev += by_group[g]
            taken += len(by_group[g])
        else:
            held += by_group[g]
    dev, held = sorted(dev), sorted(held)
    return {"dev": dev, "heldout": held, "n_free_diseases": n_free, "n_free_groups": len(free),
            "sha256_dev_ids": hash_ids(dev), "sha256_heldout_ids": hash_ids(held)}


def cache_only_annotate(symbol: str) -> dict:
    """UniProt annotation from the existing cache only (v4 sends no UniProt request); missing -> exception -> class 'other'."""
    from track2 import disease_landscape as dl
    old_cache, old_fetch = dl.CACHE, dl._fetch

    def only(url: str, retries: int = 6, timeout: int = 180) -> bytes:
        path = dl.CACHE / (hashlib.sha256(url.encode()).hexdigest()[:24] + ".bin")
        if path.exists():
            return path.read_bytes()
        raise RuntimeError("annotation not cached (v4 makes no UniProt request)")

    dl.CACHE, dl._fetch = CACHE / "landscape", only
    try:
        return dl.uniprot_annotation(symbol)
    finally:
        dl.CACHE, dl._fetch = old_cache, old_fetch


def gene_record_v4(symbol: str, annotate: Optional[Callable[[str], dict]] = None) -> dict:
    from track2 import disease_landscape as dl
    from track2.mechanism_general import classify_general
    try:
        ann = (annotate or cache_only_annotate)(symbol)
    except Exception as e:  # noqa: BLE001 - keep the gene, class 'other', failure recorded
        ann = {"terms": [], "accession": None, "error": repr(e)[:120]}
    narrow, matched = dl.classify_terms(ann["terms"])
    general, gmatched = classify_general(ann["terms"])
    rec = {"symbol": symbol, "uniprot_accession": ann.get("accession"), "mechanism_class": narrow, "matched_terms": list(matched),
           "mechanism_class_general": general, "general_matched_terms": gmatched[:12], "n_terms": len(ann["terms"])}
    if ann.get("error"):
        rec["annotation_error"] = ann["error"]
    return rec


def build_benchmark(root: Path = ROOT, hpo=None, orph=None, orph_names=None, mondo=None, mondo_version=None, primekg_ids=None,
                    universe=None, gold1=None, v2_bench=None, curated=None, labels_by_drug=None,
                    annotate: Optional[Callable[[str], dict]] = None) -> Dict[str, Any]:
    drugs = universe if universe is not None else v1.load_universe(root / "data" / "universe.json")
    uidx = v1.universe_index(drugs)
    gold1 = gold1 if gold1 is not None else v1.primekg_gold1(root, set(uidx))
    v2b = v2_bench if v2_bench is not None else json.loads((root / "out" / "engine" / "benchmark_v2.json").read_text(encoding="utf-8"))
    if hpo is None:
        src = root / "data" / "sources"
        hpo = v1.load_hpo(next(src.rglob("phenotype.hpoa")), next(src.rglob("genes_to_phenotype.txt")), next(src.rglob("hp.obo")))
    if orph is None:
        orph = v1.fetch_orphanet()
    if orph_names is None:
        orph_names = b2.parse_orphanet_names(v1._cached_bytes(b2.ORPHADATA_NAMES))
    if mondo is None:
        mondo, mondo_version = b2.load_mondo(root / "data" / "sources" / "mondo" / "mondo.obo")
    if primekg_ids is None:
        primekg_ids = b2.load_primekg_disease_ids(root / "data" / "sources" / "primekg" / "nodes.tab")
    if curated is None:
        curated = {c[0] for c in elc.CANDIDATES}
    index = b2.DisorderIndex(hpo, orph)
    v2_by_name = {d["name"]: d for d in v2b["diseases"]}
    log: List[dict] = []
    eligible: List[str] = []
    matches: Dict[str, dict] = {}
    for name in sorted(gold1):
        if name in v2_by_name:
            log.append({"disease": name, "eligible": False, "retained_from_v2": True, "reason": "retained from benchmark_v2"})
            continue
        m = b2.match_disease_v2(name, primekg_ids.get(name, []), mondo, index, orph_names)
        ok, why = evaluate_candidate_v4(name, gold1[name], m, uidx, name in curated)
        matches[name] = m
        if ok:
            eligible.append(name)
        log.append({"disease": name, "eligible": ok, "retained_from_v2": False, "reason": why})
    if labels_by_drug is None:
        labels_by_drug = {d: b2.fetch_labels(d) for d in sorted(uidx)}

    def drug_name_list(items: Sequence[dict]) -> List[str]:
        return sorted(g["name"].lower() for g in items)

    diseases: List[dict] = []
    for n, d in sorted(v2_by_name.items()):
        g1 = drug_name_list(d["gold1"])
        diseases.append({"name": n, "origin": "v2", "origin_v2": d["origin"], "group": group_id(g1), "gold1": g1,
                         "gold3": drug_name_list(d["gold3"]), "genes": d["genes"], "identifiers": d["identifiers"],
                         "inheritance": d.get("inheritance", [])})
    for n in sorted(eligible):
        m = matches[n]
        g1 = sorted(g for g in gold1[n] if g in uidx)
        g3 = sorted(d for d in b2.gold3_for_disease(sorted(set(m["name_set"]) | {b2.clean_name(n)}), labels_by_drug) if d in uidx)
        diseases.append({"name": n, "origin": "v4_new", "origin_v2": None, "group": group_id(g1), "gold1": g1, "gold3": g3,
                         "genes": [gene_record_v4(s, annotate) for s in m["genes"]],
                         "identifiers": {"primekg_name": n, "mondo": m["mondo"], "omim_xref": m["omim"], "orphanet_xref": m["orphanet_xref"],
                                         "hpo": m["hpo"], "orphanet": m["orphanet"], "id_matched": m["id_matched"], "name_set": m["name_set"]},
                         "inheritance": m["inheritance"]})
    diseases.sort(key=lambda d: (d["origin"] != "v2", d["name"]))
    split = make_split_v4(diseases)
    gsize = Counter(d["group"] for d in diseases)
    reasons = Counter(re.sub(r"\s*[:(].*$", "", e["reason"]).strip() for e in log)
    n_new = len(eligible)
    dropped = {k: v for k, v in sorted(reasons.items()) if k not in ("retained from benchmark_v2",)}
    n_dropped = sum(v for k, v in dropped.items() if k != "eligible")
    return {
        "version": "v4", "seed": SEED, "rules": {"min_genes": MIN_GENES_V4, "min_gold1": MIN_GOLD_V4, "max_matched": v1.MAX_MATCHED,
                                                 "duplicate_gold1": "grouped (D3), not dropped", "declared_changes": ["D1", "D2", "D3"]},
        "n_benchmark": len(diseases), "n_v2_retained": len(v2_by_name), "n_v4_new": n_new, "n_groups": len(gsize),
        "n_pool_evaluated": len(log), "n_dropped_total": n_dropped, "target_150_met": len(diseases) >= 150,
        "retained_and_dropped_by_reason": {"retained_from_benchmark_v2": len(v2_by_name), "eligible_new": n_new,
                                           "dropped": {k: v for k, v in dropped.items() if k != "eligible"}},
        "groups_with_more_than_one_disease": {g: n for g, n in sorted(gsize.items()) if n > 1},
        "split": split, "n_dev": len(split["dev"]), "n_heldout": len(split["heldout"]),
        "heldout_at_least_40": len(split["heldout"]) >= 40,
        "sources": {"mondo_data_version": mondo_version, "n_drugs_with_labels": sum(1 for v in labels_by_drug.values() if v)},
        "diseases": diseases, "selection_log": sorted(log, key=lambda e: e["disease"]),
    }


def load_benchmark(root: Path = ROOT) -> dict:
    return json.loads((root / "out" / "engine" / "benchmark_v4.json").read_text(encoding="utf-8"))


def run_benchmark(root: Path = ROOT) -> dict:
    t = time.time()
    b = build_benchmark(root)
    b["build_seconds"] = round(time.time() - t, 1)
    write_json(root / "out" / "engine" / "benchmark_v4.json", b)
    s = b["split"]
    print(json.dumps({k: b[k] for k in ("n_benchmark", "n_v2_retained", "n_v4_new", "n_groups", "n_dev", "n_heldout",
                                        "target_150_met", "heldout_at_least_40", "retained_and_dropped_by_reason")}, sort_keys=True), flush=True)
    append_protocol(
        "\n### Freeze record 1 - benchmark and split (written by `engine_v4 benchmark`, before any feature or score exists)\n\n"
        f"- `out/engine/benchmark_v4.json` SHA-256 `{sha256_file(root / 'out' / 'engine' / 'benchmark_v4.json')}`; "
        f"{b['n_benchmark']} diseases ({b['n_v2_retained']} from benchmark_v2 + {b['n_v4_new']} new) in {b['n_groups']} groups; "
        f"target of 150 met: {b['target_150_met']}.\n"
        f"- DEV {b['n_dev']} diseases, SHA-256 of the sorted id list `{s['sha256_dev_ids']}`.\n"
        f"- HELD-OUT {b['n_heldout']} diseases, SHA-256 of the sorted id list `{s['sha256_heldout_ids']}` "
        f"(at least 40: {b['heldout_at_least_40']}).\n"
        f"- Retained and dropped by first failing rule: `{json.dumps(b['retained_and_dropped_by_reason'], sort_keys=True)}`.\n")
    return b


# ============================================================================= models download
def download_models() -> dict:
    """The one authorised download: two Hugging Face models, weights and tokenizer files only, pinned to the current commit."""
    from huggingface_hub import HfApi, snapshot_download
    api = HfApi()
    rec = {}
    for key, spec in MODEL_SPECS.items():
        info = api.model_info(spec["repo"])
        target = MODELS_DIR / key
        target.mkdir(parents=True, exist_ok=True)
        snapshot_download(repo_id=spec["repo"], revision=info.sha, local_dir=str(target), allow_patterns=spec["files"])
        files = {f: {"sha256": sha256_file(target / f), "bytes": (target / f).stat().st_size} for f in spec["files"]}
        rec[key] = {"repo": spec["repo"], "revision_sha": info.sha, "pool": spec["pool"], "files": files}
        print(json.dumps({key: info.sha}), flush=True)
    write_json(MODELS_DIR / "revisions.json", rec)
    return rec


# ============================================================================= features
class ReqCounter:
    def __init__(self, getter: Callable[[str, dict], dict]):
        self.getter, self.n = getter, 0

    def __call__(self, url: str, params: dict) -> dict:
        self.n += 1
        return self.getter(url, params)


def hit_count(symbols: Sequence[str], getter: Callable[[str, dict], dict], cache: Dict[str, int]) -> int:
    q = es.build_gene_query(symbols, es.YEAR_PRIMARY)
    key = hashlib.sha256(q.encode()).hexdigest()[:24]
    if key not in cache:
        d = getter(es.EPMC, {"query": q, "format": "json", "resultType": "lite", "pageSize": 1})
        cache[key] = int(d.get("hitCount", 0))
    return cache[key]


def stratum_of(hits: int) -> str:
    """L < 50, M 50..500, H > 500 (Europe PMC hit count of the primary retrieval query)."""
    return "L" if hits < 50 else ("M" if hits <= 500 else "H")


def disease_text(name: str, symbols: Sequence[str]) -> str:
    return f"{b2.clean_name(name)}. Genes: {', '.join(sorted(set(symbols)))}"


def drug_sentences(papers: Sequence[dict], symbols: Sequence[str], dre: "re.Pattern[str]", gate: Optional[set], k: int = N_SENT) -> Dict[str, str]:
    """Per drug: the first k sentences (paper order, then sentence order) of the retrieved papers that contain a disease gene
    symbol (paper level, as the v1 paper count) and mention the drug (sentence level, v1 drug match). Text only; no label."""
    grx = es.gene_regex(symbols)
    acc: Dict[str, List[str]] = {}
    for p in papers:
        if not grx.search(p["text"]):
            continue
        for sent in elc.split_sentences(p["text"]):
            if gate is not None and not (set(re.findall(r"\w+", sent.lower())) & gate):
                continue
            for name in sorted({m.group(0).lower() for m in dre.finditer(sent)}):
                lst = acc.setdefault(name, [])
                if len(lst) < k:
                    lst.append(sent)
    return {n: " ".join(s) for n, s in acc.items()}


_W: Dict[str, Any] = {}


def _init_worker(names: List[str]) -> None:
    drugs = v1.load_universe(ROOT / "data" / "universe.json")
    universe = v1.universe_index(drugs)
    pop, total = es.fetch_popularity(names, None, CACHE / "epmc_popularity_all.json", None)  # fully cached
    targets, ppi = es.load_prior_graph(ROOT, CACHE / "primekg_prior_edges.json")
    _W.update({"names": names, "dre": elc.drug_regex(names), "gate": es.name_gate(names), "pop": pop, "total": total,
               "targets": targets, "ppi": ppi, "universe": universe})


def compute_disease_features(symbols: Sequence[str], mech_terms: Sequence[str], papers: Sequence[dict]) -> Dict[str, Any]:
    """Features of ONE disease from its retrieved papers: no gold argument exists (leakage audit L2)."""
    names, dre, gate, pop = _W["names"], _W["dre"], _W["gate"], _W["pop"]
    t0 = time.time()
    lit = es.literature_counts(papers, symbols, dre, gate)
    prior = es.prior_features(symbols, names, _W["targets"], _W["ppi"])
    X = es.disease_matrix(lit, names, pop, _W["total"], prior)
    proxy_counts = elc.score_abstracts([(p["pid"], p["text"]) for p in papers], mech_terms, dre)
    texts = drug_sentences(papers, symbols, dre, gate)
    return {"X": X, "popularity": np.array([pop.get(n, 0) for n in names], dtype=float),
            "proxy": np.array([proxy_counts.get(n, 0) for n in names], dtype=float), "n_papers_with_gene": lit["n_papers_with_gene"],
            "drug_texts": texts, "feature_seconds": time.time() - t0}


def _feat_job(args: Tuple[str, Sequence[str], Sequence[str], Sequence[dict]]) -> Tuple[str, Dict[str, Any]]:
    name, symbols, mech, papers = args
    return name, compute_disease_features(symbols, mech, papers)


def proxy_terms(dis: dict, v1_terms: Dict[str, dict]) -> List[str]:
    """Proxy-baseline mechanism terms, as engine_score_v2.build_features: v1 diseases from the v1 benchmark, others by the narrow
    rule on the cached annotation (cache only)."""
    symbols = sorted({g["symbol"] for g in dis["genes"]})
    if dis["name"] in v1_terms and dis.get("origin_v2") == "v1":
        return sorted(set(symbols) | {t for g in v1_terms[dis["name"]]["genes"] if g["mechanism_class"] != "other" for t in g["matched_terms"]})
    extra: Set[str] = set()
    for g in dis["genes"]:
        if g["mechanism_class"] != "other":
            try:
                from track2 import disease_landscape as dl
                extra |= set(dl.classify_terms(cache_only_annotate(g["symbol"])["terms"])[1])
            except Exception:  # noqa: BLE001 - symbols only
                pass
    return sorted(set(symbols) | extra)


def run_features(root: Path = ROOT, workers: int = 6) -> dict:
    bench = load_benchmark(root)
    v1b = json.loads((root / "out" / "engine" / "benchmark.json").read_text(encoding="utf-8"))
    v1_terms = {d["name"]: d for d in v1b["diseases"]}
    drugs = v1.load_universe(root / "data" / "universe.json")
    names = sorted(v1.universe_index(drugs))
    done: Dict[str, Any] = {}
    if FEATURES_FILE.exists():
        import pickle
        with FEATURES_FILE.open("rb") as fh:  # project-own cache written by this module
            done = pickle.load(fh)
    hits: Dict[str, int] = json.loads(HITS_FILE.read_text(encoding="utf-8")) if HITS_FILE.exists() else {}
    jobs, meta = [], {}
    n_total = len(bench["diseases"])
    for i, dis in enumerate(bench["diseases"]):
        d = dis["name"]
        if d in done:
            continue
        symbols = sorted({g["symbol"] for g in dis["genes"]})
        counter = ReqCounter(v2s.epmc_getter)
        t0 = time.time()
        papers = es.fetch_papers(symbols, es.YEAR_PRIMARY, CACHE / "epmc", counter)
        h = hit_count(symbols, counter, hits)
        retr_s = time.time() - t0
        meta[d] = {"symbols": symbols, "n_retrieved": len(papers), "hit_count": h, "retrieval_requests": counter.n,
                   "retrieval_seconds": round(retr_s, 3), "pids": [p["pid"] for p in papers]}
        jobs.append((d, symbols, proxy_terms(dis, v1_terms), papers))
        if (i + 1) % 10 == 0:
            HITS_FILE.write_text(json.dumps(hits), encoding="utf-8")
            print(json.dumps({"retrieved": i + 1, "of": n_total}), flush=True)
    HITS_FILE.write_text(json.dumps(hits), encoding="utf-8")
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(names,)) as ex:
        for name, res in ex.map(_feat_job, jobs, chunksize=2):
            res.update(meta[name])
            done[name] = res
    import pickle
    with FEATURES_FILE.open("wb") as fh:
        pickle.dump(done, fh)
    # Orphanet prevalence file (evidence-gap block), cached under the landscape cache
    try:
        v1._cached_bytes(ORPHA_PREV_URL)
    except Exception as e:  # noqa: BLE001
        print("orphanet prevalence not available:", repr(e)[:100], flush=True)
    print(json.dumps({"features": len(done), "new": len(jobs)}), flush=True)
    return done


# ============================================================================= embeddings
def run_embed(root: Path = ROOT) -> dict:
    import pickle
    bench = load_benchmark(root)
    with FEATURES_FILE.open("rb") as fh:
        feats = pickle.load(fh)
    rev = json.loads((MODELS_DIR / "revisions.json").read_text(encoding="utf-8"))
    texts: Dict[str, int] = {}
    for dis in bench["diseases"]:
        f = feats[dis["name"]]
        for t in [disease_text(dis["name"], f["symbols"])] + [f["drug_texts"][k] for k in sorted(f["drug_texts"])]:
            texts.setdefault(t, len(texts))
    lst = [t for t, _ in sorted(texts.items(), key=lambda kv: kv[1])]
    EMB_DIR.mkdir(parents=True, exist_ok=True)
    src = EMB_DIR / "texts.json"
    src.write_text(json.dumps(lst), encoding="utf-8")
    stats = {"n_unique_texts": len(lst)}
    worker = ROOT / "src" / "track2" / "embed_worker.py"
    for key, spec in MODEL_SPECS.items():
        dst = EMB_DIR / f"{key}.npy"
        if dst.exists():
            stats[key] = {"cached": True}
            continue
        cmd = [sys.executable, "-I", str(worker), str(MODELS_DIR / key), spec["pool"], str(src), str(dst), str(MAX_TOK), "32"]
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(EMB_DIR))
        if r.returncode != 0:
            raise RuntimeError(f"embedding worker failed for {key}: {r.stderr[-800:]}")
        stats[key] = json.loads(r.stdout.strip().splitlines()[-1])
        print(json.dumps({key: stats[key]}), flush=True)
    stats["revisions"] = {k: v["revision_sha"] for k, v in rev.items()}
    write_json(EMB_DIR / "embed_stats.json", stats)
    return stats


def load_embeddings() -> Tuple[Dict[str, int], Dict[str, np.ndarray]]:
    lst = json.loads((EMB_DIR / "texts.json").read_text(encoding="utf-8"))
    return {t: i for i, t in enumerate(lst)}, {k: np.load(EMB_DIR / f"{k}.npy") for k in MODEL_SPECS}


def emb_block(name: str, f: Dict[str, Any], names: Sequence[str], tindex: Dict[str, int], emb: Dict[str, np.ndarray]) -> np.ndarray:
    """(n_drugs, 3): cos(BiomedBERT), cos(SapBERT), has_text. Disease text vs drug text; pairs without sentences are 0 / 0 / 0."""
    out = np.zeros((len(names), 3))
    dvec = {k: emb[k][tindex[disease_text(name, f["symbols"])]] for k in emb}
    for i, n in enumerate(names):
        t = f["drug_texts"].get(n)
        if t is None:
            continue
        j = tindex[t]
        out[i, 0] = float(emb["biomedbert"][j] @ dvec["biomedbert"])
        out[i, 1] = float(emb["sapbert"][j] @ dvec["sapbert"])
        out[i, 2] = 1.0
    return out


# ============================================================================= data assembly
def assemble(root: Path, bench: dict, subset: Sequence[str]) -> Dict[str, Any]:
    """Z (8 v1 + 3 embedding columns), labels, baselines, strata for the given diseases. Labels of a disease are read only
    for the diseases asked for."""
    import pickle
    with FEATURES_FILE.open("rb") as fh:
        feats = pickle.load(fh)
    drugs = v1.load_universe(root / "data" / "universe.json")
    universe = v1.universe_index(drugs)
    names = sorted(universe)
    index = {n: i for i, n in enumerate(names)}
    classes = np.array(es.drug_classes(universe, names))
    tindex, emb = load_embeddings()
    by_name = {d["name"]: d for d in bench["diseases"]}
    Z, y, pop, proxy, info = {}, {}, {}, {}, {}
    for d in sorted(subset):
        f = feats[d]
        Z[d] = np.column_stack([f["X"], emb_block(d, f, names, tindex, emb)])
        yy = np.zeros(len(names))
        for g in by_name[d]["gold1"]:
            if g in index:
                yy[index[g]] = 1
        y[d] = yy
        pop[d], proxy[d] = f["popularity"], f["proxy"]
        info[d] = {"hit_count": f["hit_count"], "stratum": stratum_of(f["hit_count"]), "n_retrieved": f["n_retrieved"],
                   "n_papers_with_gene": f["n_papers_with_gene"], "n_gold1": int(yy.sum()), "group": by_name[d]["group"],
                   "n_pairs_with_text": int(Z[d][:, 10].sum())}
    return {"Z": Z, "y": y, "pop": pop, "proxy": proxy, "names": names, "classes": classes, "info": info, "feats": feats}


def members_of(bench: dict, diseases: Sequence[str]) -> Dict[str, List[str]]:
    gmap = {d["name"]: d["group"] for d in bench["diseases"]}
    m: Dict[str, List[str]] = defaultdict(list)
    for d in sorted(diseases):
        m[gmap[d]].append(d)
    return dict(m)


def group_folds(members: Dict[str, List[str]], seed: int = SEED) -> List[List[str]]:
    """5 folds over groups (v1 fold maker on the group keys); a group is never split across folds."""
    keys = sorted(members)
    return [sorted(d for g in f for d in members[g]) for f in es.make_folds(keys, seed)]


# ============================================================================= heads
def sigmoid(z: np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(-np.clip(z, -35, 35)))


def fit_head(kind: str, X: np.ndarray, y: np.ndarray, C: float = 1.0) -> Dict[str, Any]:
    mu, sd = X.mean(0), X.std(0)
    sd[sd < 1e-12] = 1.0
    Xs = (X - mu) / sd
    if kind == "logit":
        clf = LogisticRegression(C=C, max_iter=3000).fit(Xs, y)
        return {"kind": "logit", "mu": mu, "sd": sd, "coef": clf.coef_[0].copy(), "intercept": float(clf.intercept_[0]), "prev": float(y.mean())}
    clf = MLPClassifier(hidden_layer_sizes=(MLP_PARAMS["hidden"],), alpha=MLP_PARAMS["alpha"], max_iter=MLP_PARAMS["max_iter"],
                        random_state=SEED).fit(Xs, y)
    return {"kind": "mlp", "mu": mu, "sd": sd, "W0": clf.coefs_[0].copy(), "b0": clf.intercepts_[0].copy(),
            "W1": clf.coefs_[1][:, 0].copy(), "b1": float(clf.intercepts_[1][0]), "prev": float(y.mean())}


def head_decision(h: Dict[str, Any], X: np.ndarray) -> np.ndarray:
    Xs = (X - np.asarray(h["mu"])) / np.asarray(h["sd"])
    if h["kind"] == "logit":
        return Xs @ np.asarray(h["coef"]) + h["intercept"]
    hid = np.maximum(Xs @ np.asarray(h["W0"]) + np.asarray(h["b0"]), 0.0)
    return hid @ np.asarray(h["W1"]) + h["b1"]


def head_params(h: Dict[str, Any], cols: Sequence[int]) -> Dict[str, Any]:
    out = {k: (np.asarray(v).tolist() if isinstance(v, np.ndarray) else v) for k, v in h.items()}
    out["cols"] = list(cols)
    return out


def fit_rows(kind: str, cols: Sequence[int], train: Sequence[str], Z, y, keep: Optional[np.ndarray]) -> Dict[str, Any]:
    keep = np.ones(len(next(iter(y.values()))), dtype=bool) if keep is None else keep
    X = np.vstack([Z[d][keep][:, cols] for d in train])
    yy = np.concatenate([y[d][keep] for d in train])
    return fit_head(kind, X, yy)


def oof_head(kind: str, cols: Sequence[int], dev: Sequence[str], members: Dict[str, List[str]], Z, y, classes, nested_platt: bool,
             seed: int = SEED) -> Dict[str, Dict[str, np.ndarray]]:
    """Out-of-fold scores per disease: group folds, ATC-class masking of training rows (v1), optional nested Platt map fitted on
    inner out-of-fold scores of the training groups only."""
    folds = group_folds(members, seed)
    res: Dict[str, Dict[str, np.ndarray]] = {}
    for fold in folds:
        train = [d for d in dev if d not in fold]
        keep = es._class_mask(classes, es._gold_classes(fold, y, classes))
        a = b = None
        if nested_platt:
            tm = {g: [d for d in ds if d in set(train)] for g, ds in members.items()}
            tm = {g: ds for g, ds in tm.items() if ds}
            s_in, y_in = [], []
            for ifold in group_folds(tm, SEED + 1):
                itrain = [d for d in train if d not in ifold]
                ikeep = es._class_mask(classes, es._gold_classes(ifold, y, classes))
                h = fit_rows(kind, cols, itrain, Z, y, ikeep)
                for d in ifold:
                    s_in.append(head_decision(h, Z[d][:, cols]))
                    y_in.append(y[d])
            a, b = es.logit_fit(np.concatenate(s_in), np.concatenate(y_in))
        h = fit_rows(kind, cols, train, Z, y, keep)
        for d in fold:
            s = head_decision(h, Z[d][:, cols])
            res[d] = {"score": s, "prev": np.full(len(s), h["prev"])}
            if nested_platt:
                res[d]["platt"] = sigmoid(a + b * s)
    return res


def add_platt(pred: Dict[str, Dict[str, np.ndarray]], dev: Sequence[str], y) -> Tuple[float, float]:
    """Pooled Platt map on the out-of-fold scores of the given diseases; adds 'platt_pooled' to every entry."""
    a, b = es.logit_fit(np.concatenate([pred[d]["score"] for d in dev]), np.concatenate([y[d] for d in dev]))
    for d in dev:
        pred[d]["platt_pooled"] = sigmoid(a + b * pred[d]["score"])
    return float(a), float(b)


# ============================================================================= evaluation with groups
def agg_groups(vals: np.ndarray, dnames: Sequence[str], members: Dict[str, List[str]]) -> np.ndarray:
    pos = {d: i for i, d in enumerate(dnames)}
    out = []
    for g in sorted(members):
        v = vals[[pos[d] for d in members[g]]]
        out.append(float(np.nanmean(v)) if (~np.isnan(v)).any() else float("nan"))
    return np.array(out)


def evaluate_groups(members: Dict[str, List[str]], preds: Dict[str, Dict[str, Dict[str, np.ndarray]]], y, pop, proxy,
                    reference: Optional[str] = None, platt_key: str = "platt") -> Tuple[Dict[str, Any], Dict[str, np.ndarray]]:
    """Baselines + models, three v1 criteria (es.decide) with a bootstrap that resamples GROUPS (a group = one observation, the
    mean over its member diseases). members = {group: [diseases]} restricted to the set being evaluated."""
    gs = sorted(members)
    dn = sorted(d for g in gs for d in members[g])
    out: Dict[str, Any] = {"n_diseases": len(dn), "n_groups": len(gs), "models": {}}
    if not gs:
        return out, {}
    idx = np.random.default_rng(SEED).integers(0, len(gs), size=(es.N_BOOT, len(gs)))

    def grp(m: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        return {k: agg_groups(v, dn, members) for k, v in m.items()}

    pop_m = grp(es.per_disease_metrics(pop, y, dn))
    prx_m = grp(es.per_disease_metrics(proxy, y, dn))
    rnd_m = grp(es.random_per_disease(y, dn))
    base = {n: {"auc": es.macro_with_ci(m["auc"], idx), "precision_at_10": es.macro_with_ci(m["p10"], idx),
                "recall_at_20": es.macro_with_ci(m["r20"], idx)} for n, m in (("popularity", pop_m), ("proxy", prx_m), ("random", rnd_m))}
    best_name = max(("popularity", "proxy", "random"), key=lambda n: base[n]["precision_at_10"]["estimate"] or -1)
    best_m = {"popularity": pop_m, "proxy": prx_m, "random": rnd_m}[best_name]
    out["baselines"], out["best_baseline_by_precision_at_10"] = base, best_name
    arrays: Dict[str, np.ndarray] = {"popularity": pop_m["auc"]}
    ms: Dict[str, Dict[str, np.ndarray]] = {}
    for name, pred in preds.items():
        raw = es.per_disease_metrics({d: pred[d]["score"] for d in dn}, y, dn)
        m = grp(raw)
        ms[name] = m
        arrays[name] = m["auc"]
        e: Dict[str, Any] = {"auc": es.macro_with_ci(m["auc"], idx), "precision_at_10": es.macro_with_ci(m["p10"], idx),
                             "recall_at_20": es.macro_with_ci(m["r20"], idx),
                             "per_disease_auc": {d: es._r(v) for d, v in zip(dn, raw["auc"])},
                             "auc_minus_popularity": es.diff_with_ci(m["auc"], pop_m["auc"], idx),
                             "auc_minus_proxy": es.diff_with_ci(m["auc"], prx_m["auc"], idx),
                             "auc_minus_random": es.diff_with_ci(m["auc"], rnd_m["auc"], idx),
                             "precision_at_10_minus_best_baseline": es.diff_with_ci(m["p10"], best_m["p10"], idx)}
        try:
            pg = {g: {"platt": np.concatenate([pred[d][platt_key] for d in members[g]]),
                      "prev": np.concatenate([pred[d]["prev"] for d in members[g]])} for g in gs}
            yg = {g: np.concatenate([y[d] for d in members[g]]) for g in gs}
            e["calibration_platt"] = es.calibration_summary(pg, yg, gs, idx, "platt")
            e["decision"] = es.decide(e)
        except Exception as ex:  # noqa: BLE001 - tiny strata
            e["calibration_platt"], e["decision"] = None, {"error": repr(ex)[:120]}
        out["models"][name] = e
    if reference is not None and reference in ms:
        for name in preds:
            if name != reference:
                out["models"][name]["auc_minus_" + reference] = es.diff_with_ci(ms[name]["auc"], ms[reference]["auc"], idx)
                out["models"][name]["precision_at_10_minus_" + reference] = es.diff_with_ci(ms[name]["p10"], ms[reference]["p10"], idx)
    return out, arrays


def stratified(members_all: Dict[str, List[str]], info: Dict[str, dict], preds, y, pop, proxy, reference: Optional[str],
               platt_key: str = "platt") -> Dict[str, Any]:
    res: Dict[str, Any] = {}
    for s in STRATA:
        mem = {g: [d for d in ds if info[d]["stratum"] == s] for g, ds in members_all.items()}
        mem = {g: ds for g, ds in mem.items() if ds}
        ev, _ = evaluate_groups(mem, preds, y, pop, proxy, reference, platt_key)
        ev["descriptive_only_n_too_small"] = len(mem) < MIN_GROUPS_STRATUM
        res[s] = ev
    return res


def paired_group_diff(a: np.ndarray, b: np.ndarray) -> Dict[str, Any]:
    d = a - b
    d = d[~np.isnan(d)]
    n = len(d)
    se = float(d.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
    return {"mean": float(d.mean()) if n else None, "se": se, "n_groups": n, "noise_threshold_2se": 2 * se, "mde_2p8se": 2.8 * se}


def one_se_inclusion(a: np.ndarray, b: np.ndarray) -> Dict[str, Any]:
    """Protocol 4.2: candidate a enters only if the mean paired difference exceeds 1 SE of the paired differences and 0."""
    p = paired_group_diff(a, b)
    p["included"] = bool(p["mean"] is not None and p["mean"] > 0 and p["mean"] > p["se"])
    return p


# ============================================================================= conformal abstention
def conformal_q(p: np.ndarray, y: np.ndarray, alpha: float = ALPHA) -> Dict[str, float]:
    """Mondrian split conformal thresholds: label 1 nonconformity 1 - p, label 0 nonconformity p; rank ceil((n+1)(1-alpha));
    2.0 (always included) when the rank exceeds n."""
    q = {}
    for lab in (0, 1):
        a = np.sort((1 - p[y == 1]) if lab == 1 else p[y == 0])
        k = math.ceil((len(a) + 1) * (1 - alpha))
        q[str(lab)] = float(a[k - 1]) if 1 <= k <= len(a) else 2.0
    return q


def conformal_sets(p: np.ndarray, q: Dict[str, float]) -> Tuple[np.ndarray, np.ndarray]:
    return p <= q["0"], (1 - p) <= q["1"]


def conformal_stats(p: np.ndarray, y: np.ndarray, q: Dict[str, float], top_mask: Optional[np.ndarray] = None) -> Dict[str, Any]:
    h0, h1 = conformal_sets(p, q)
    size = h0.astype(int) + h1.astype(int)
    cov = np.where(y == 1, h1, h0)
    out = {"n_pairs": int(len(p)), "abstention_rate": float((size != 1).mean()), "empty_set_rate": float((size == 0).mean()),
           "two_label_set_rate": float((size == 2).mean()), "marginal_coverage": float(cov.mean()),
           "positive_coverage": float(cov[y == 1].mean()) if (y == 1).any() else None,
           "negative_coverage": float(cov[y == 0].mean()) if (y == 0).any() else None, "n_positive_pairs": int(y.sum())}
    if top_mask is not None and top_mask.any():
        out["abstention_rate_top10"] = float((size[top_mask] != 1).mean())
        out["confident_positive_in_top10"] = float(((size == 1) & h1)[top_mask].mean())
    return out


def top10_mask(score: np.ndarray) -> np.ndarray:
    m = np.zeros(len(score), dtype=bool)
    m[np.argsort(-score, kind="stable")[:10]] = True
    return m


def conformal_report(members: Dict[str, List[str]], info, pred, y, key: str, q: Dict[str, float]) -> Dict[str, Any]:
    ds = sorted(d for g in members for d in members[g])

    def block(sel: Sequence[str]) -> Optional[Dict[str, Any]]:
        if not sel:
            return None
        p = np.concatenate([pred[d][key] for d in sel])
        yy = np.concatenate([y[d] for d in sel])
        top = np.concatenate([top10_mask(pred[d]["score"]) for d in sel])
        return conformal_stats(p, yy, q, top)

    res = {"overall": block(ds)}
    for s in STRATA:
        res[s] = block([d for d in ds if info[d]["stratum"] == s])
    return res


def dev_conformal_validity(members: Dict[str, List[str]], pred, y, key: str, seed: int = SEED) -> Dict[str, Any]:
    """Cross-conformal check on DEV: thresholds from the out-of-fold pairs of the other folds, positive coverage on the fold."""
    folds = group_folds(members, seed)
    gmap = {d: g for g, ds in members.items() for d in ds}
    cov, npos, dis, grp = [], [], [], []
    for fold in folds:
        other = [d for d in gmap if d not in fold]
        q = conformal_q(np.concatenate([pred[d][key] for d in other]), np.concatenate([y[d] for d in other]))
        for d in fold:
            p, yy = pred[d][key], y[d]
            _, h1 = conformal_sets(p, q)
            cov.append(int(h1[yy == 1].sum()))
            npos.append(int((yy == 1).sum()))
            grp.append(gmap[d])
    gs = sorted(set(grp))
    c = np.array([sum(c_ for c_, g in zip(cov, grp) if g == k) for k in gs], dtype=float)
    n = np.array([sum(n_ for n_, g in zip(npos, grp) if g == k) for k in gs], dtype=float)
    est = float(c.sum() / n.sum())
    idx = np.random.default_rng(SEED).integers(0, len(gs), size=(es.N_BOOT, len(gs)))
    boots = np.array([c[r].sum() / max(n[r].sum(), 1) for r in idx])
    se = float(boots.std(ddof=1))
    return {"positive_coverage": est, "bootstrap_se": se, "threshold_0.9_minus_se": 0.9 - se, "valid": bool(est >= 0.9 - se)}


# ============================================================================= freezing
def freeze_models(models: Dict[str, Any], extra: Dict[str, Any]) -> Dict[str, Any]:
    body = dict(models)
    body.update(extra)
    return {"model": body, "sha256_of_model": v3.canonical_hash(body)}


def thaw_v4(frozen: Dict[str, Any]) -> Dict[str, Any]:
    if v3.canonical_hash(frozen["model"]) != frozen["sha256_of_model"]:
        raise RuntimeError("frozen v4 model hash mismatch")
    return frozen["model"]


def predict_frozen_v4(m: Dict[str, Any], key: str, Zd: np.ndarray) -> Dict[str, np.ndarray]:
    h = m["heads"][key]
    s = head_decision(h, Zd[:, h["cols"]])
    pl = m["platt"][key]
    return {"score": s, "platt": sigmoid(pl["a"] + pl["b"] * s), "prev": np.full(len(s), h["prev"])}


# ============================================================================= evidence gap (hierarchical borrowing)
def parse_orphanet_prevalence(xml: bytes) -> Dict[str, List[dict]]:
    """Orphadata en_product9_prev: 'ORPHA:code' -> [{type, geo, cls, val}]."""
    out: Dict[str, List[dict]] = {}
    for d in ET.fromstring(xml).iter("Disorder"):
        code = (d.findtext("OrphaCode") or "").strip()
        recs = []
        for p in d.iter("Prevalence"):
            val = (p.findtext("ValMoy") or "").strip()
            try:
                v = float(val) if val else None
            except ValueError:
                v = None
            recs.append({"type": (p.findtext("./PrevalenceType/Name") or "").strip(),
                         "geo": (p.findtext("./PrevalenceGeographic/Name") or "").strip(),
                         "cls": (p.findtext("./PrevalenceClass/Name") or "").strip(), "val": v})
        if code:
            out[f"ORPHA:{code}"] = recs
    return out


_CLASS_AT_LEAST_1_PER_MILLION = re.compile(r"^(1-9 / 1 000 000|1-9 / 100 000|1-5 / 10 000|6-9 / 10 000|>1 / 1000)$")


def tier_from_records(recs: Sequence[dict]) -> str:
    """Case-count tier (protocol 6): Cases/families value if present; else a class >= 1-9 / 1 000 000 implies > 5; else unknown."""
    vals = [r["val"] for r in recs if r["type"].lower().startswith("cases/families") and r["val"] is not None and r["val"] > 0]
    if vals:
        v = max(vals)
        return "1 known case" if round(v) <= 1 else ("2-5" if round(v) <= 5 else ">5")
    if any(_CLASS_AT_LEAST_1_PER_MILLION.match(r["cls"]) for r in recs):
        return ">5"
    return "unknown"


TIER_ORDER = {"unknown": 0, "1 known case": 1, "2-5": 2, ">5": 3}
DESIGN_BY_TIER = {
    "1 known case": "Single-patient (N-of-1) crossover design with a pre-specified outcome, analysed with a hierarchical model that borrows from the mechanism class.",
    "2-5": "Series of N-of-1 trials analysed jointly with a Bayesian hierarchical model that borrows across the mechanism class.",
    ">5": "Multi-centre single-arm study against a natural-history control, with class-level borrowing as an informative prior.",
    "unknown": "Ascertain the number of patients through registries before choosing a design.",
}


def disease_tier(dis: dict, prev: Dict[str, List[dict]]) -> Tuple[str, List[str]]:
    ident = dis["identifiers"]
    codes = [f"ORPHA:{c}" for c in ident.get("orphanet_xref", [])] or [c for c in ident.get("id_matched", []) if c.startswith("ORPHA:")]
    tiers = [tier_from_records(prev[c]) for c in codes if c in prev]
    return (max(tiers, key=lambda t: TIER_ORDER[t]) if tiers else "unknown"), codes


def hanley_mcneil_var(auc: float, n_pos: int, n_neg: int) -> float:
    q1, q2 = auc / (2 - auc), 2 * auc * auc / (1 + auc)
    return (auc * (1 - auc) + (n_pos - 1) * (q1 - auc * auc) + (n_neg - 1) * (q2 - auc * auc)) / (n_pos * n_neg)


def dersimonian_laird(y: np.ndarray, v: np.ndarray) -> Dict[str, float]:
    w = 1 / v
    mu_fixed = float((w * y).sum() / w.sum())
    q = float((w * (y - mu_fixed) ** 2).sum())
    k = len(y)
    denom = float(w.sum() - (w ** 2).sum() / w.sum())
    tau2 = max(0.0, (q - (k - 1)) / denom) if denom > 0 and k > 1 else 0.0
    wr = 1 / (v + tau2)
    return {"tau2": tau2, "mu": float((wr * y).sum() / wr.sum()), "k": k}


def evidence_gap(bench: dict, aucs: Dict[str, float], n_pos: Dict[str, int], n_drugs: int, prev: Dict[str, List[dict]],
                 strata: Dict[str, str]) -> Dict[str, Any]:
    """Per disease: tier, ORPHA codes, borrowing pool (same plurality general mechanism class), empirical-Bayes weight of the
    disease's own AUC, and the fixed design text. Computed, never fed back into a model; no dose or drug advice."""
    cls_of = {d["name"]: v3.disease_class_genes(d["genes"])[0] for d in bench["diseases"]}
    pools: Dict[str, List[str]] = defaultdict(list)
    for d, c in cls_of.items():
        if c is not None and d in aucs and not math.isnan(aucs[d]):
            pools[c].append(d)
    dl: Dict[str, Dict[str, float]] = {}
    for c, ds in pools.items():
        if len(ds) >= 3:
            y = np.array([aucs[d] for d in ds])
            v = np.array([max(hanley_mcneil_var(min(max(aucs[d], 1e-3), 1 - 1e-3), max(n_pos[d], 1), n_drugs - n_pos[d]), 1e-8) for d in ds])
            dl[c] = dersimonian_laird(y, v)
    rows = {}
    for dis in bench["diseases"]:
        n = dis["name"]
        tier, codes = disease_tier(dis, prev)
        c = cls_of[n]
        row = {"tier": tier, "orpha_codes": codes, "stratum": strata.get(n), "mechanism_class": c,
               "borrowing_pool_size": max(len(pools.get(c, [])) - (1 if n in pools.get(c, []) else 0), 0) if c else 0,
               "design": DESIGN_BY_TIER[tier]}
        if c in dl and n in aucs and not math.isnan(aucs[n]):
            v = max(hanley_mcneil_var(min(max(aucs[n], 1e-3), 1 - 1e-3), max(n_pos[n], 1), n_drugs - n_pos[n]), 1e-8)
            w = dl[c]["tau2"] / (dl[c]["tau2"] + v)
            row.update({"own_auc": aucs[n], "class_mean_auc": dl[c]["mu"], "weight_on_own_auc": w,
                        "borrowed_auc": dl[c]["mu"] + w * (aucs[n] - dl[c]["mu"]), "class_tau2": dl[c]["tau2"]})
        rows[n] = row
    return {"rows": rows, "tier_counts": dict(Counter(r["tier"] for r in rows.values())),
            "classes": {c: dl[c] for c in sorted(dl)}, "note": "Hypotheses for study design only: no dose, no drug recommendation."}


# ============================================================================= leakage audit
def leakage_audit(root: Path, bench: dict, feats: Dict[str, Any]) -> Dict[str, Any]:
    import inspect
    out: Dict[str, Any] = {}
    prior = json.loads((CACHE / "primekg_prior_edges.json").read_text(encoding="utf-8"))
    pk_names = set(json.loads((CACHE / "primekg_indications.json").read_text(encoding="utf-8")))
    pk_names |= {d["name"] for d in bench["diseases"]}
    keys = set(prior["targets"]) | set(prior["ppi"])
    vals = {x for v in prior["targets"].values() for x in v} | {x for v in prior["ppi"].values() for x in v}
    src = inspect.getsource(es.load_prior_graph)
    out["L1_prior_graph"] = {"relations_filter_in_source": '["drug_protein", "protein_protein"]' in src,
                             "disease_names_checked": len(pk_names), "disease_names_found_as_key_or_value": len(pk_names & (keys | vals)),
                             "relations_in_filter": re.findall(r"isin\(\[([^\]]*)\]", src),
                             "indication_relations_read": any("indication" in r for r in re.findall(r"isin\(\[([^\]]*)\]", src)),
                             "passed": bool(not (pk_names & (keys | vals)) and '["drug_protein", "protein_protein"]' in src and re.findall(r"isin\(\[([^\]]*)\]", src) == ['"drug_protein", "protein_protein"'])}
    g3 = {d["name"]: d["gold3"] for d in bench["diseases"]}
    g1 = {d["name"]: d["gold1"] for d in bench["diseases"]}
    params = list(inspect.signature(compute_disease_features).parameters)
    out["L2_gold3_and_texts"] = {"compute_disease_features_parameters": params, "gold_argument_present": any("gold" in p for p in params),
                                 "n_diseases_with_gold3": sum(1 for v in g3.values() if v),
                                 "n_gold3_pairs": sum(len(v) for v in g3.values()),
                                 "n_gold3_pairs_also_gold1": sum(len(set(g3[d]) & set(g1[d])) for d in g3),
                                 "gold3_used_as_feature_or_training_label": False,
                                 "embedding_texts_source": "retrieved abstract sentences only", "passed": not any("gold" in p for p in params)}
    split = bench["split"]
    gmap = {d["name"]: d["group"] for d in bench["diseases"]}
    shared_groups = {gmap[d] for d in split["dev"]} & {gmap[d] for d in split["heldout"]}
    out["L3_group_split"] = {"groups_in_both_dev_and_heldout": len(shared_groups), "passed": not shared_groups}
    dev_genes = {g["symbol"] for d in bench["diseases"] if d["name"] in set(split["dev"]) for g in d["genes"]}
    dev_pids: Set[str] = set()
    for d in split["dev"]:
        dev_pids |= set(feats[d]["pids"])
    ho = [d for d in bench["diseases"] if d["name"] in set(split["heldout"])]
    shared_gene = [d["name"] for d in ho if {g["symbol"] for g in d["genes"]} & dev_genes]
    fr = [len(set(feats[d["name"]]["pids"]) & dev_pids) / max(len(feats[d["name"]]["pids"]), 1) for d in ho]
    out["L5_overlap_heldout_vs_dev"] = {"heldout_diseases_sharing_a_gene_with_dev": len(shared_gene), "n_heldout": len(ho),
                                        "fraction_diseases_sharing_gene": len(shared_gene) / max(len(ho), 1),
                                        "mean_fraction_of_heldout_papers_also_retrieved_for_a_dev_disease": float(np.mean(fr)) if fr else None,
                                        "note": "reported, not removed: features are literature counts; labels of DEV are not used by held-out features"}
    out["L6_atc_class_masking"] = "training rows exclude drugs whose ATC level-3 class holds a gold 1 drug of any test-fold disease (es._class_mask)"
    return out


# ============================================================================= drivers
def _fold_check(members: Dict[str, List[str]], folds: List[List[str]]) -> bool:
    where = {}
    for i, f in enumerate(folds):
        for d in f:
            where[d] = i
    return all(len({where[d] for d in ds}) == 1 for ds in members.values())


def run_dev(root: Path = ROOT) -> dict:
    t_start = time.time()
    bench = load_benchmark(root)
    split = bench["split"]
    if hash_ids(split["dev"]) != split["sha256_dev_ids"] or hash_ids(split["heldout"]) != split["sha256_heldout_ids"]:
        raise RuntimeError("split hashes do not match benchmark_v4.json")
    dev = split["dev"]
    members = members_of(bench, dev)
    D = assemble(root, bench, dev)           # labels and features of DEV diseases only
    Z, y, pop, proxy, classes, names, info = D["Z"], D["y"], D["pop"], D["proxy"], D["classes"], D["names"], D["info"]
    folds = group_folds(members)
    assert _fold_check(members, folds), "a group straddles folds"
    # sanity: v1 base features equal the earlier case-study cache for the diseases it holds (read-only, project-own file)
    check = {"case_study_cache_available": False}
    pkl = CACHE / "case_study_features.pkl"
    if pkl.exists():
        import pickle
        with pkl.open("rb") as fh:
            blob = pickle.load(fh)["per"]
        same = [bool(np.array_equal(blob[d]["primary"]["X"], Z[d][:, :8])) for d in dev if d in blob]
        check = {"case_study_cache_available": True, "n_compared": len(same), "all_identical": bool(all(same))}
    t_fit = time.time()
    pred = {"B0": oof_head("logit", B0_COLS, dev, members, Z, y, classes, True),
            "E": oof_head("logit", E_COLS, dev, members, Z, y, classes, True),
            "M": oof_head("mlp", E_COLS, dev, members, Z, y, classes, False)}
    platt = {k: add_platt(pred[k], dev, y) for k in pred}
    ev_all, arr = evaluate_groups(members, {"B0": pred["B0"], "E": pred["E"], "M": pred["M"]}, y, pop, proxy, reference="B0", platt_key="platt_pooled")
    ev_nested, _ = evaluate_groups(members, {"B0": pred["B0"], "E": pred["E"]}, y, pop, proxy, reference="B0", platt_key="platt")
    incl_E = one_se_inclusion(arr["E"], arr["B0"])
    incl_M = {"eligible": False, "reason": "ablation only (protocol section 3)", "paired_vs_B0": paired_group_diff(arr["M"], arr["B0"])}
    final = "E" if incl_E["included"] else "B0"
    strat = stratified(members, info, {"B0": pred["B0"], "E": pred["E"], "M": pred["M"]}, y, pop, proxy, "B0", "platt_pooled")
    # ---- conformal on the final model (frozen Platt applied to pooled OOF scores)
    p_dev = np.concatenate([pred[final][d]["platt_pooled"] for d in dev])
    y_dev = np.concatenate([y[d] for d in dev])
    q = conformal_q(p_dev, y_dev)
    conf_dev = conformal_report(members, info, pred[final], y, "platt_pooled", q)
    conf_valid = dev_conformal_validity(members, pred[final], y, "platt_pooled")
    # ---- refit on all DEV (no masking), freeze
    heads = {}
    for key, kind, cols in (("B0", "logit", B0_COLS), ("E", "logit", E_COLS), ("M", "mlp", E_COLS)):
        h = head_params(fit_rows(kind, cols, dev, Z, y, None), cols)
        heads[key] = h
    rev = json.loads((MODELS_DIR / "revisions.json").read_text(encoding="utf-8"))
    extra = {"label": LABEL, "protocol": "protocol/addendum_engine_v4.md", "final_model": final, "dev": dev, "heldout": split["heldout"],
             "sha256_dev_ids": split["sha256_dev_ids"], "sha256_heldout_ids": split["sha256_heldout_ids"],
             "platt": {k: {"a": platt[k][0], "b": platt[k][1]} for k in platt},
             "conformal": {"alpha": ALPHA, "q": q, "calibration": "pooled DEV out-of-fold probabilities of the final model (frozen Platt)",
                           "model": final},
             "embedding_models": {k: v["revision_sha"] for k, v in rev.items()}, "feature_names": es.FEATURE_NAMES + EMB_NAMES,
             "n_universe_drugs": len(names)}
    frozen = freeze_models({"heads": heads}, extra)
    write_json(OUT / "engine_v4_frozen.json", frozen)
    fh = sha256_file(OUT / "engine_v4_frozen.json")
    # ---- dry run of the frozen path on three DEV diseases (smoke test; no held-out data)
    mdl = thaw_v4(json.loads((OUT / "engine_v4_frozen.json").read_text(encoding="utf-8")))
    dry_ok = all(np.isfinite(predict_frozen_v4(mdl, final, Z[d])["score"]).all() for d in dev[:3])
    out = {
        "label": LABEL, "protocol": "protocol/addendum_engine_v4.md", "protocol_sha256_at_dev_run": es.sha256(PROTOCOL),
        "benchmark_sha256": es.sha256(root / "out" / "engine" / "benchmark_v4.json"),
        "status": "v1 (supported, 17 diseases), v2 and v3 results are unchanged; this is post hoc development. HELD-OUT not yet scored.",
        "settings": {"seed": SEED, "n_bootstrap": es.N_BOOT, "alpha": ALPHA, "strata": {"L": "<50", "M": "50-500", "H": ">500"},
                     "n_dev_diseases": len(dev), "n_dev_groups": len(members), "n_heldout": len(split["heldout"]),
                     "folds": folds, "feature_names": es.FEATURE_NAMES + EMB_NAMES, "thresholds": es.THRESHOLDS,
                     "mlp": MLP_PARAMS, "n_sentences": N_SENT, "max_tokens": MAX_TOK},
        "unvalidated_registry": UNVALIDATED, "feature_cache_check": check,
        "benchmark_counts": {k: bench[k] for k in ("n_benchmark", "n_v2_retained", "n_v4_new", "n_groups", "n_dev", "n_heldout",
                                                   "target_150_met", "heldout_at_least_40", "retained_and_dropped_by_reason")},
        "dev": {"per_disease_auc_final_model": {d: es._r(es.auc_mw(pred[final][d]["score"], y[d])) for d in dev},
                "disease_info": info, "evaluation_overall_pooled_platt": ev_all, "evaluation_overall_nested_platt": ev_nested,
                "by_stratum": strat, "stratum_sizes": {s: {"diseases": sum(1 for d in dev if info[d]["stratum"] == s),
                                                           "groups": len({info[d]["group"] for d in dev if info[d]["stratum"] == s})} for s in STRATA}},
        "inclusion": {"E_vs_B0_1se": incl_E, "M": incl_M, "final_model": final,
                      "K_conformal": {"dev_cross_conformal_validity": conf_valid, "included": bool(conf_valid["valid"]),
                                      "dev_in_sample_report": conf_dev, "thresholds": q}},
        "frozen": {"file": "out/engine/engine_v4_frozen.json", "model_body_sha256": frozen["sha256_of_model"], "file_sha256": fh,
                   "final_model": final, "dry_run_finite": bool(dry_ok)},
        "timings": {"dev_total_seconds": round(time.time() - t_start, 1), "dev_fit_cv_seconds": round(time.time() - t_fit, 1)},
        "heldout": None,
    }
    write_json(OUT / "validation_v4.json", out)
    append_protocol(
        "\n### Freeze record 2 - model (written by `engine_v4 dev`, before any HELD-OUT score exists)\n\n"
        f"- Final model by the 1-SE rule on DEV: `{final}` (E minus B0 mean paired group AUC {incl_E['mean']:.4f}, SE {incl_E['se']:.4f}, "
        f"included: {incl_E['included']}); MLP is ablation only.\n"
        f"- Conformal K (alpha {ALPHA}): DEV cross-conformal positive coverage {conf_valid['positive_coverage']:.3f} "
        f"(0.9 minus SE = {conf_valid['threshold_0.9_minus_se']:.3f}), valid: {conf_valid['valid']}.\n"
        f"- `out/engine/engine_v4_frozen.json`: SHA-256 of the model body `{frozen['sha256_of_model']}`; file SHA-256 `{fh}`.\n"
        f"- Embedding model revisions: `{json.dumps({k: v['revision_sha'] for k, v in rev.items()})}`.\n"
        "- Clarification applied in code: the drug text takes sentences that mention the drug from retrieved papers that contain a "
        "disease gene symbol (paper level, as the v1 paper count).\n"
        "- The held-out command verifies this hash and the split hashes, writes `out/engine/heldout_v4_scored.lock` and refuses a second run.\n")
    print(json.dumps({"final": final, "E_vs_B0": incl_E, "conformal_valid": conf_valid, "frozen_sha": frozen["sha256_of_model"]},
                     default=v3._default, sort_keys=True), flush=True)
    return out


def run_heldout(root: Path = ROOT) -> dict:
    lock = OUT / "heldout_v4_scored.lock"
    if lock.exists():
        raise RuntimeError("HELD-OUT v4 was already scored once (lock file exists); refusing a second scoring")
    bench = load_benchmark(root)
    split = bench["split"]
    frozen = json.loads((OUT / "engine_v4_frozen.json").read_text(encoding="utf-8"))
    m = thaw_v4(frozen)
    if hash_ids(split["heldout"]) != m["sha256_heldout_ids"] or hash_ids(split["dev"]) != m["sha256_dev_ids"]:
        raise RuntimeError("split differs from the frozen file")
    v1f = json.loads((OUT / "engine_v1_weights_frozen.json").read_text(encoding="utf-8"))
    m1 = v3.thaw(v1f)
    held, dev = split["heldout"], split["dev"]
    D = assemble(root, bench, held)          # features only so far; labels are read here, after the freeze
    lock.write_text(json.dumps({"scored_at": "first and only held-out scoring", "frozen_model_sha256": frozen["sha256_of_model"],
                                "v1_weights_sha256": v1f["sha256_of_model"]}) + "\n", encoding="utf-8")
    Z, y, pop, proxy, info = D["Z"], D["y"], D["pop"], D["proxy"], D["info"]
    members = members_of(bench, held)
    final = m["final_model"]
    t0 = time.time()
    preds = {k: {d: predict_frozen_v4(m, k, Z[d]) for d in held} for k in ("B0", "E", "M")}
    cls_id = np.zeros(len(D["names"]), dtype=int)          # unused by the v1-weights base model (no F3 table)
    preds["v1_weights_v3freeze"] = {d: v3.predict_frozen(m1, Z[d][:, :8], cls_id) for d in held}
    score_seconds = time.time() - t0
    ev, arr = evaluate_groups(members, preds, y, pop, proxy, reference="B0")
    strat = stratified(members, info, preds, y, pop, proxy, "B0")
    ev_dis, _ = evaluate_groups({d: [d] for d in held}, preds, y, pop, proxy, reference="B0")
    q = m["conformal"]["q"]
    conf = conformal_report(members, info, preds[final], y, "platt", q)
    conf_sets = {}
    for s in STRATA:
        mem = {g: [d for d in ds if info[d]["stratum"] == s] for g, ds in members.items()}
        conf_sets[s] = len(sum(mem.values(), []))
    # ---- evidence-gap block (DEV out-of-fold AUCs stored by the DEV run + frozen HELD-OUT AUCs)
    out = json.loads((OUT / "validation_v4.json").read_text(encoding="utf-8"))
    n_drugs = len(D["names"])
    gold_n = {d["name"]: len(d["gold1"]) for d in bench["diseases"]}
    auc_all = {d: (v if v is not None else float("nan")) for d, v in out["dev"]["per_disease_auc_final_model"].items()}
    npos_all = {d: gold_n[d] for d in dev}
    strata_all = {d: out["dev"]["disease_info"][d]["stratum"] for d in dev}
    for d in held:
        auc_all[d] = es.auc_mw(preds[final][d]["score"], y[d])
        npos_all[d] = int(y[d].sum())
        strata_all[d] = info[d]["stratum"]
    try:
        prev = parse_orphanet_prevalence(v1._cached_bytes(ORPHA_PREV_URL))
    except Exception as e:  # noqa: BLE001
        prev = {}
        print("orphanet prevalence unavailable:", repr(e)[:100], flush=True)
    gap = evidence_gap(bench, auc_all, npos_all, n_drugs, prev, strata_all)
    gap["orphanet_prevalence_records_available"] = len(prev)
    audit = leakage_audit(root, bench, D["feats"])
    # ---- cost per disease
    cost = cost_report(D["feats"], root, bench, score_seconds, len(held))
    improvement = {k: ev["models"][k].get("auc_minus_B0") for k in ("E", "M", "v1_weights_v3freeze")}
    heldout = {"label": LABEL + "; HELD-OUT scored once", "n_diseases": len(held), "n_groups": len(members), "diseases": held,
               "frozen_model_sha256": frozen["sha256_of_model"], "final_model": final, "evaluation_group_bootstrap": ev,
               "evaluation_disease_level_bootstrap_sensitivity": ev_dis, "by_stratum": strat,
               "stratum_sizes_diseases": conf_sets, "disease_info": info,
               "criteria_final_model": ev["models"][final]["decision"], "criteria_all_models": {k: v.get("decision") for k, v in ev["models"].items()},
               "paired_vs_B0": improvement,
               "conformal": {"alpha": ALPHA, "thresholds": q, "report": conf, "applied_to_model": final}}
    out["heldout"] = heldout
    out["status"] = "v1, v2 and v3 results are unchanged; post hoc development. HELD-OUT scored once (see 'heldout')."
    out["evidence_gap_hierarchical_borrowing"] = gap
    out["leakage_audit"] = audit
    out["cost_per_disease"] = cost
    out["predicted_direction_check"] = direction_check(strat, out["dev"]["by_stratum"], conf)
    write_json(OUT / "validation_v4.json", out)
    print(json.dumps({"final": final, "criteria": heldout["criteria_final_model"], "n": len(held)}, default=v3._default, sort_keys=True), flush=True)
    return out


def direction_check(ho_strat: Dict[str, Any], dev_strat: Dict[str, Any], conf: Dict[str, Any]) -> Dict[str, Any]:
    """Computed comparison with the directions stated before scoring (P1, P4); the model key is read from the first model."""
    def gap(strat: Dict[str, Any], s: str, model: str) -> Optional[float]:
        m = strat.get(s, {}).get("models", {}).get(model)
        return None if not m else m["auc_minus_popularity"]["estimate"]
    res = {}
    for name, strat in (("dev", dev_strat), ("heldout", ho_strat)):
        g = {s: gap(strat, s, "B0") for s in STRATA}
        res[name] = {"auc_minus_popularity_by_stratum": g,
                     "P1_L_smaller_than_H": None if g["L"] is None or g["H"] is None else bool(g["L"] < g["H"]),
                     "P2_ordered_L_lt_M_le_H": None if None in g.values() else bool(g["L"] < g["M"] <= g["H"])}
    ab = {s: (conf.get(s) or {}).get("abstention_rate") for s in STRATA}
    res["P4_abstention_by_stratum_heldout"] = ab
    res["P4_L_higher_than_H"] = None if ab["L"] is None or ab["H"] is None else bool(ab["L"] > ab["H"])
    return res


def cost_report(feats: Dict[str, Any], root: Path, bench: dict, score_seconds: float, n_scored: int) -> Dict[str, Any]:
    emb = json.loads((EMB_DIR / "embed_stats.json").read_text(encoding="utf-8")) if (EMB_DIR / "embed_stats.json").exists() else {}
    emb_sec = sum(v.get("seconds", 0.0) for k, v in emb.items() if isinstance(v, dict))
    n_tx = max(emb.get("n_unique_texts", 1), 1)
    rows = {}
    for dis in bench["diseases"]:
        f = feats[dis["name"]]
        cold = f["retrieval_requests"] > 0
        est_req = 1 + math.ceil(f["n_retrieved"] / 100)
        n_txt = 1 + len(f["drug_texts"])
        rows[dis["name"]] = {"retrieval_requests": f["retrieval_requests"] if cold else est_req,
                             "requests_kind": "measured (cold)" if cold else "estimated (cached earlier)",
                             "retrieval_seconds": f["retrieval_seconds"] if cold else None, "feature_seconds": round(f["feature_seconds"], 3),
                             "embedding_seconds_share": round(emb_sec * n_txt / n_tx, 3), "n_texts": n_txt}
    cold_rows = [r for r in rows.values() if r["requests_kind"].startswith("measured")]

    def stat(v: Sequence[float]) -> Dict[str, Optional[float]]:
        v = [x for x in v if x is not None]
        return {"mean": float(np.mean(v)), "median": float(np.median(v)), "max": float(np.max(v))} if v else {"mean": None, "median": None, "max": None}
    total_sec = [(r["retrieval_seconds"] or 0) + r["feature_seconds"] + r["embedding_seconds_share"] for r in rows.values()]
    return {"n_diseases": len(rows), "n_measured_cold": len(cold_rows), "requests_per_disease": stat([r["retrieval_requests"] for r in rows.values()]),
            "requests_per_disease_measured_cold": stat([r["retrieval_requests"] for r in cold_rows]),
            "retrieval_seconds_cold": stat([r["retrieval_seconds"] for r in cold_rows]), "feature_seconds": stat([r["feature_seconds"] for r in rows.values()]),
            "embedding_seconds_share": stat([r["embedding_seconds_share"] for r in rows.values()]),
            "scoring_seconds_per_disease_frozen_models": score_seconds / max(n_scored, 1),
            "total_seconds_per_disease_cold_estimate": stat(total_sec),
            "extrapolation_1000_diseases_hours": float(np.mean(total_sec) * 1000 / 3600), "embedding_total_seconds": emb_sec,
            "note": "Retrieval is rate-limited to 4 requests/s; cached diseases' retrieval seconds were not measured.",
            "per_disease": rows}


def main() -> None:
    cmds = {"benchmark": run_benchmark, "download": download_models, "features": run_features, "embed": run_embed,
            "dev": run_dev, "heldout": run_heldout}
    if len(sys.argv) != 2 or sys.argv[1] not in cmds:
        raise SystemExit("usage: python -m track2.engine_v4 " + "|".join(cmds))
    cmds[sys.argv[1]]()


if __name__ == "__main__":
    main()
