r"""Benchmark v2 analysis of the evidence engine - 'v2 extension, not pre-registered primary result'.

Run:  PYTHONPATH=src python -m track2.engine_score_v2     (writes out/engine/validation_v2.json, sorted keys)

Protocol: protocol/addendum_engine_v2.md. All v1 thresholds, endpoints, baselines, bootstrap, grouped CV, ablations and
reporting rules (protocol/addendum_engine.md) are unchanged and are reused by importing track2.engine_score (not edited).
The v1-only result in out/engine/validation.json stays the primary pre-registered result; nothing here overturns it.

What is computed
  * Golds on benchmark v2 (45 diseases = 17 v1 + 28 new): gold 1 (endpoint gold), gold 3 (openFDA labels) and the union
    gold 1 or 3, each for (i) all 45 diseases and (ii) the 28 new diseases alone (held-out). For a gold other than gold 1 a
    disease enters only if its gold is non-empty (an empty gold gives no AUC; the cross-validation uses the same diseases).
  * Arm A (full), ablations B, C, D, baselines popularity / proxy / random, 95 % disease-level bootstrap, the three
    pre-registered criteria evaluated verbatim on arm A (engine_score.decide) for every set.
  * Temporal analysis with the corrected gold 2 (EXPERIMENTAL arms), literature <= 2018, no prior (arms B and D).
  * AUC by general mechanism class (rule in MECH_RULE below), and a post hoc sensitivity analysis without 'common' diseases.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import threading
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from track2 import empirical_literature_control as elc
from track2 import engine_benchmark as eb
from track2 import engine_score as es

ROOT = Path(__file__).resolve().parents[2]
LABEL = "v2 extension, not pre-registered primary result"
MAX_REQUESTS = 4000          # progress guard: hard cap on new network requests of one run
GOLDS = ("gold1", "gold3", "union")

# ----- post hoc sensitivity: 'common' (non-rare) diseases. LABELLED POST HOC, defined from the v2 hand-back + a documented rule.
COMMON_EXPLICIT = ("small cell lung carcinoma", "ovarian cancer", "systemic lupus erythematosus (disease)",
                   "Parkinson disease", "Alzheimer disease")
COMMON_INHERITANCE_RE = re.compile(r"non-mendelian|polygenic|multigenic|multifactorial", re.I)
COMMON_RULE_TEXT = ("Post hoc sensitivity (not pre-registered). Rule S1 = the five diseases named in the v2 hand-back "
                    "(small cell lung carcinoma, ovarian cancer, systemic lupus erythematosus, Parkinson disease, Alzheimer "
                    "disease). Rule S2 = S1 plus every disease of the v2 file whose Orphanet/HPO inheritance annotation contains "
                    "a complex-disease mode (Non-Mendelian, Polygenic, Multigenic/multifactorial). 'Sporadic' alone does not "
                    "flag a disease. The v2 file itself carries no 'common' flag.")
MECH_RULE = ("A disease belongs to a general mechanism class when strictly more than 50 % of its genes carry that "
             "mechanism_class_general (genes of class 'other' count in the denominator); a disease with no such majority is "
             "'no_majority_class' and is not assigned. Plurality assignment (largest class, ties by alphabetical order) is "
             "reported as a supplementary variant. Class AUC = macro mean over assigned diseases of the disease AUC of arm A.")


# ============================================================================= request guard
class RequestGuard:
    """Wraps a getter: counts new network requests, writes a progress file, aborts after MAX_REQUESTS."""

    def __init__(self, getter: Callable[[str, dict], dict], progress_file: Optional[Path] = None, limit: int = MAX_REQUESTS,
                 every: int = 25):
        self.getter, self.file, self.limit, self.every = getter, progress_file, limit, every
        self.n, self.lock = 0, threading.Lock()

    def __call__(self, url: str, params: dict) -> dict:
        with self.lock:
            self.n += 1
            n = self.n
        if n > self.limit:
            raise RuntimeError(f"request guard: more than {self.limit} new requests, aborting")
        if self.file is not None and (n == 1 or n % self.every == 0):
            self.file.parent.mkdir(parents=True, exist_ok=True)
            self.file.write_text(json.dumps({"requests": n, "limit": self.limit}), encoding="utf-8")
        return self.getter(url, params)


def epmc_getter(url: str, params: dict) -> dict:
    """Europe PMC GET (POST for very long gene queries), <= 4 req/s through the engine_score limiter, retries."""
    import requests
    import time
    for attempt in range(7):
        try:
            es._LIMITER.wait()
            hdr = {"User-Agent": "hselith-track2/0.3 (research)"}
            if len(params.get("query", "")) > 1500:
                r = requests.post(url, data=params, timeout=90, headers=hdr)
            else:
                r = requests.get(url, params=params, timeout=90, headers=hdr)
            r.raise_for_status()
            return r.json()
        except Exception:  # noqa: BLE001
            if attempt == 6:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


# ============================================================================= gold and class helpers
def gold_vectors(dis: dict, index: Dict[str, int], n: int) -> Dict[str, np.ndarray]:
    def vec(items: Sequence[dict]) -> np.ndarray:
        v = np.zeros(n)
        for g in items:
            i = index.get(g["name"].lower())
            if i is not None:
                v[i] = 1
        return v
    g1, g2, g3 = vec(dis["gold1"]), vec(dis["gold2"]), vec(dis["gold3"])
    return {"gold1": g1, "gold2": g2, "gold3": g3, "union": np.maximum(g1, g3)}


def disease_mech_class(genes: Sequence[dict], key: str = "mechanism_class_general") -> Dict[str, Optional[str]]:
    """Majority (> 50 %) class and plurality class of a disease."""
    from collections import Counter
    c = Counter(g.get(key) or "other" for g in genes)
    if not c:
        return {"majority": None, "plurality": None}
    top, k = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0]
    return {"majority": top if k * 2 > sum(c.values()) else None, "plurality": top}


def common_flags(diseases: Sequence[dict]) -> Dict[str, Dict[str, bool]]:
    out = {}
    for d in diseases:
        explicit = d["name"] in COMMON_EXPLICIT
        rule = any(COMMON_INHERITANCE_RE.search(t) for t in d.get("inheritance", []))
        out[d["name"]] = {"S1": explicit, "S2": explicit or rule}
    return out


def class_auc_table(auc_by_d: Dict[str, float], base_by_d: Dict[str, float], mech: Dict[str, Dict[str, Optional[str]]],
                    variant: str = "majority", n_boot: int = es.N_BOOT) -> Dict[str, Any]:
    """Macro AUC of the engine (and of popularity) per class, disease-level bootstrap within class (percentile 95 %)."""
    by: Dict[str, List[str]] = {}
    for d, v in auc_by_d.items():
        c = mech.get(d, {}).get(variant)
        if c is not None and v is not None and not np.isnan(v):
            by.setdefault(c, []).append(d)
    out = {}
    for c in sorted(by):
        ds = sorted(by[c])
        a = np.array([auc_by_d[d] for d in ds])
        b = np.array([base_by_d[d] for d in ds])
        rng = np.random.default_rng(es.SEED)
        idx = rng.integers(0, len(ds), size=(n_boot, len(ds)))
        out[c] = {"n_diseases": len(ds), "diseases": ds,
                  "engine_auc": es.macro_with_ci(a, idx), "popularity_auc": es.macro_with_ci(b, idx),
                  "engine_minus_popularity": es.diff_with_ci(a, b, idx)}
    return out


# ============================================================================= data assembly
def build_features(root: Path, bench2: dict, v1_bench: dict, getter: Optional[Callable[[str, dict], dict]] = None) -> Dict[str, Any]:
    drugs = eb.load_universe(root / "data" / "universe.json")
    universe = eb.universe_index(drugs)
    names = sorted(universe)
    index = {n: i for i, n in enumerate(names)}
    dre = elc.drug_regex(names)
    gate = es.name_gate(names)
    classes = np.array(es.drug_classes(universe, names))
    targets, ppi = es.load_prior_graph(root, root / "data" / "cache_engine" / "primekg_prior_edges.json")
    pop_all, tot_all = es.fetch_popularity(names, None, root / "data" / "cache_engine" / "epmc_popularity_all.json", getter)
    pop_18, tot_18 = es.fetch_popularity(names, es.YEAR_TEMPORAL, root / "data" / "cache_engine" / "epmc_popularity_2018.json", getter)
    cache = root / "data" / "cache_engine" / "epmc"
    v1_terms = {d["name"]: d for d in v1_bench["diseases"]}

    per: Dict[str, Dict[str, Any]] = {}
    for dis in bench2["diseases"]:
        name = dis["name"]
        symbols = sorted({g["symbol"] for g in dis["genes"]})
        if name in v1_terms:   # identical proxy mechanism terms as v1
            mech = sorted(set(symbols) | {t for g in v1_terms[name]["genes"] if g["mechanism_class"] != "other" for t in g["matched_terms"]})
            mech_src = "v1 benchmark matched_terms"
        else:                  # narrow v1 rule re-applied to the cached UniProt annotation
            extra: set = set()
            for g in dis["genes"]:
                if g["mechanism_class"] != "other":
                    try:
                        from track2 import disease_landscape as dl
                        from track2 import engine_benchmark_v2 as ebv2
                        extra |= set(dl.classify_terms(ebv2._default_annotate(g["symbol"])["terms"])[1])
                    except Exception:  # noqa: BLE001 - symbols only
                        pass
            mech = sorted(set(symbols) | extra)
            mech_src = "narrow rule on cached annotation"
        rec: Dict[str, Any] = {"symbols": symbols, "mech_terms": mech, "mech_terms_source": mech_src}
        rec.update(gold_vectors(dis, index, len(names)))
        for period, year, pop, tot in (("primary", es.YEAR_PRIMARY, pop_all, tot_all), ("temporal", es.YEAR_TEMPORAL, pop_18, tot_18)):
            papers = es.fetch_papers(symbols, year, cache, getter)
            lit = es.literature_counts(papers, symbols, dre, gate)
            prior = es.prior_features(symbols, names, targets, ppi) if period == "primary" else None
            proxy_counts = elc.score_abstracts([(p["pid"], p["text"]) for p in papers], mech, dre)
            rec[period] = {"n_papers_retrieved": len(papers), "n_papers_with_gene": lit["n_papers_with_gene"],
                           "X": es.disease_matrix(lit, names, pop, tot, prior),
                           "popularity": np.array([pop.get(n, 0) for n in names], dtype=float),
                           "proxy": np.array([proxy_counts.get(n, 0) for n in names], dtype=float)}
        per[name] = rec
    return {"per": per, "classes": classes, "names": names}


# ============================================================================= analysis of one set
def analyse_set(per: Dict[str, Dict[str, Any]], classes: np.ndarray, diseases: Sequence[str], gold: str) -> Dict[str, Any]:
    """Primary-period analysis (all four arms, three criteria on arm A) for one gold and one disease subset."""
    sel = sorted(d for d in diseases if per[d][gold].sum() >= 1)
    y = {d: per[d][gold] for d in sel}
    res = es.analyse(sel, {d: per[d]["primary"]["X"] for d in sel}, y, {d: per[d]["primary"]["popularity"] for d in sel},
                     {d: per[d]["primary"]["proxy"] for d in sel}, classes, es.ARMS)
    res["decision_arm_A"] = es.decide(res["arms"]["A"])
    res["diseases_excluded_empty_gold"] = sorted(set(diseases) - set(sel))
    res["label"] = LABEL
    # per-disease baseline AUCs (popularity) for the class table
    res["per_disease_popularity_auc"] = {d: es._r(es.auc_mw(per[d]["primary"]["popularity"], y[d])) for d in sel}
    return res


def analyse_temporal(per: Dict[str, Dict[str, Any]], classes: np.ndarray, diseases: Sequence[str]) -> Dict[str, Any]:
    sel2 = sorted(d for d in diseases if per[d]["gold2"].sum() >= 1)
    arms = {k: v for k, v in es.ARMS.items() if not v[0]}
    y2 = {d: per[d]["gold2"] for d in sel2}
    tX = {d: per[d]["temporal"]["X"] for d in sel2}
    tp = {d: per[d]["temporal"]["popularity"] for d in sel2}
    tq = {d: per[d]["temporal"]["proxy"] for d in sel2}
    out: Dict[str, Any] = {"label": LABEL, "gold": "corrected gold 2 (EXPERIMENTAL arms, start > 2018-12-31)",
                           "note": "Literature <= 2018, knowledge-graph prior excluded: only arms B and D are valid.",
                           "diseases_excluded_empty_gold2": sorted(set(diseases) - set(sel2))}
    out["all_gold2"] = es.analyse(sel2, tX, y2, tp, tq, classes, arms) if len(sel2) >= 2 else None
    keep = {d: per[d]["gold1"] == 0 for d in sel2}
    sel3 = [d for d in sel2 if per[d]["gold2"][keep[d]].sum() >= 1]
    out["excluding_gold1_drugs"] = es.analyse(sel3, tX, y2, tp, tq, classes, arms, drop_mask=keep) if len(sel3) >= 2 else None
    return out


def _job(args: Tuple[str, str, Any, Any, Any, Any]) -> Tuple[str, Dict[str, Any]]:
    key, kind, per, classes, diseases, gold = args
    if kind == "temporal":
        return key, analyse_temporal(per, classes, diseases)
    return key, analyse_set(per, classes, diseases, gold)


# ============================================================================= driver
def run(root: Path = ROOT, getter: Optional[Callable[[str, dict], dict]] = None, out_path: Optional[Path] = None,
        workers: int = 6) -> Dict[str, Any]:
    b2_path = root / "out" / "engine" / "benchmark_v2.json"
    b1_path = root / "out" / "engine" / "benchmark.json"
    bench2 = json.loads(b2_path.read_text(encoding="utf-8"))
    v1_bench = json.loads(b1_path.read_text(encoding="utf-8"))
    guard = RequestGuard(getter or epmc_getter, root / "data" / "cache_engine" / "v2_progress.json")
    feats = build_features(root, bench2, v1_bench, guard)
    per, classes = feats["per"], feats["classes"]
    all_d = sorted(per)
    origin = {d["name"]: d["origin"] for d in bench2["diseases"]}
    new_d = sorted(d for d in all_d if origin[d] == "v2_new")
    v1_d = sorted(d for d in all_d if origin[d] == "v1")
    flags = common_flags(bench2["diseases"])
    mech = {d["name"]: disease_mech_class(d["genes"]) for d in bench2["diseases"]}
    subsets = {"all_45": all_d, "new_28_heldout": new_d}

    jobs: List[Tuple[str, str, Any, Any, Any, Any]] = []
    for g in GOLDS:
        for sk, sv in subsets.items():
            jobs.append((f"{g}|{sk}", "set", per, classes, sv, g))
    jobs.append(("temporal|all_45", "temporal", per, classes, all_d, None))
    jobs.append(("temporal|new_28_heldout", "temporal", per, classes, new_d, None))
    jobs.append(("replication_v1_gold1|v1_17", "set", per, classes, v1_d, "gold1"))
    for rule in ("S1", "S2"):
        for g in ("gold1", "union"):
            for sk, sv in subsets.items():
                kept = [d for d in sv if not flags[d][rule]]
                jobs.append((f"sens_{rule}_{g}|{sk}", "set", per, classes, kept, g))

    results: Dict[str, Any] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for key, res in ex.map(_job, jobs):
            results[key] = res
            print("done", key, file=sys.stderr, flush=True)

    def slim(r: Dict[str, Any]) -> Dict[str, Any]:
        return r

    sets: Dict[str, Any] = {}
    for g in GOLDS:
        sets[g] = {sk: results[f"{g}|{sk}"] for sk in subsets}
    # per-mechanism-class AUC (arm A, all 45 and, separately, the 28 new) for each gold
    class_tables: Dict[str, Any] = {}
    for g in GOLDS:
        for sk in subsets:
            r = results[f"{g}|{sk}"]
            aucs = {d: v for d, v in r["arms"]["A"]["per_disease_auc"].items()}
            pops = {d: v for d, v in r["per_disease_popularity_auc"].items()}
            class_tables[f"{g}|{sk}"] = {
                "rule": MECH_RULE,
                "majority": class_auc_table(aucs, pops, mech, "majority"),
                "plurality_supplementary": class_auc_table(aucs, pops, mech, "plurality"),
                "diseases_without_majority_class": sorted(d for d in aucs if mech[d]["majority"] is None)}
    sens = {"label": "POST HOC sensitivity analysis, " + LABEL, "rule": COMMON_RULE_TEXT,
            "excluded_S1": sorted(d for d in all_d if flags[d]["S1"]), "excluded_S2": sorted(d for d in all_d if flags[d]["S2"]),
            "results": {k[len("sens_"):]: v for k, v in results.items() if k.startswith("sens_")}}

    # reproduction check against the v1 result (read-only)
    repl = results["replication_v1_gold1|v1_17"]
    v1val_path = root / "out" / "engine" / "validation.json"
    check = {"v2_pipeline_v1_17_gold1_armA_auc": repl["arms"]["A"]["auc"]["estimate"]}
    if v1val_path.exists():
        v1val = json.loads(v1val_path.read_text(encoding="utf-8"))
        check["validation_json_v1_armA_auc"] = v1val["main"]["primary_gold1"]["arms"]["A"]["auc"]["estimate"]
        check["validation_json_v1_armA_decision"] = v1val["main"]["primary_gold1"]["decision_arm_A"]["verdict"]
        check["validation_json_sha256"] = es.sha256(v1val_path)
        check["identical_auc"] = check["v2_pipeline_v1_17_gold1_armA_auc"] == check["validation_json_v1_armA_auc"]

    out = {
        "label": LABEL,
        "status": "v1-only result (out/engine/validation.json, 17 diseases, gold 1) remains the primary pre-registered result; "
                  "this file is a documented extension and does not change that decision",
        "protocol": "protocol/addendum_engine_v2.md; thresholds and endpoints of protocol/addendum_engine.md unchanged",
        "benchmark_v2_sha256": es.sha256(b2_path), "addendum_v2_sha256": es.sha256(root / "protocol" / "addendum_engine_v2.md"),
        "settings": {"seed": es.SEED, "n_bootstrap": es.N_BOOT, "thresholds": es.THRESHOLDS, "n_diseases_all": len(all_d),
                     "n_diseases_v1": len(v1_d), "n_diseases_new_heldout": len(new_d), "new_diseases": new_d,
                     "arms": {k: {"knowledge_graph_prior": v[0], "popularity_correction": v[1]} for k, v in es.ARMS.items()},
                     "feature_names": es.FEATURE_NAMES, "request_guard_limit": MAX_REQUESTS,
                     "new_network_requests_this_run": guard.n},
        "gold_definitions": bench2["gold_rules"],
        "gold_sizes_summary": bench2["summary"],
        "gold_handling": "Gold 1 is the endpoint gold. For gold 3 and the union a disease enters an analysis only if that gold "
                         "is non-empty (diseases_excluded_empty_gold lists the others); the cross-validation uses the same diseases.",
        "per_disease_retrieval": {d: {"origin": origin[d], "n_gold1": int(per[d]["gold1"].sum()), "n_gold2": int(per[d]["gold2"].sum()),
                                      "n_gold3": int(per[d]["gold3"].sum()), "n_union": int(per[d]["union"].sum()),
                                      "n_symbols": len(per[d]["symbols"]),
                                      "mech_terms_source": per[d]["mech_terms_source"],
                                      "n_papers_retrieved_primary": per[d]["primary"]["n_papers_retrieved"],
                                      "n_papers_retrieved_temporal": per[d]["temporal"]["n_papers_retrieved"],
                                      "mechanism_class_majority": mech[d]["majority"], "mechanism_class_plurality": mech[d]["plurality"],
                                      "flag_common_S1": flags[d]["S1"], "flag_common_S2": flags[d]["S2"]} for d in all_d},
        "primary_sets": sets,
        "temporal_gold2_corrected": {sk: results[f"temporal|{sk}"] for sk in subsets},
        "auc_by_general_mechanism_class": class_tables,
        "sensitivity_excluding_common_diseases_post_hoc": sens,
        "reproduction_check_v1": check,
        "deviations_and_clarifications": [
            "The v1 pipeline (engine_score) is imported unchanged; features of the 17 v1 diseases are recomputed from the same caches.",
            "Gold 3 and union analyses keep only diseases with a non-empty gold; this is a clarification, v1 had no empty gold 1.",
            "Proxy baseline mechanism terms: v1 diseases use the v1 benchmark matched_terms; new diseases use the narrow v1 rule on the cached "
            "UniProt annotation plus gene symbols.",
            "Europe PMC POST is used for gene queries longer than 1500 characters (same endpoint, same parameters).",
            "The 'common disease' exclusion and the per-class AUC table are post hoc/descriptive and change no decision.",
            "The temporal analysis for the v2 corrected gold 2 has no pre-registered threshold; it is descriptive.",
        ],
        "known_limits": [
            "Gold 1 is incomplete and partly circular; gold 3 is a text match without negation handling.",
            "Several new diseases are common or non-rare (multi-subtype umbrella terms); see the post hoc sensitivity analysis.",
            "Bootstrap intervals are disease-level; small subsets (held-out, non-empty gold 3) give wide intervals.",
        ],
    }
    out_path = out_path or root / "out" / "engine" / "validation_v2.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True, ensure_ascii=False, default=_default) + "\n", encoding="utf-8")
    return out


def _default(o: Any) -> Any:
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, tuple)):
        return sorted(o)
    raise TypeError(type(o))


def main() -> None:
    res = run()
    for g in GOLDS:
        for sk, r in res["primary_sets"][g].items():
            print(g, sk, json.dumps(r["decision_arm_A"], sort_keys=True))


if __name__ == "__main__":
    main()
