r"""Worked case of the evidence engine on a v2 disease: kidney Wilms tumor.

LABEL: worked case on a v2 disease, not part of the pre-registered primary validation.

Run:  PYTHONPATH=src python -m track2.engine_case_study     (from the repository root; writes out/engine/case_wilms.json
      with sorted keys and out/engine/case_wilms.md)

Purpose: show the full output format of the engine on a paediatric rare disease whose answer is known (approved treatment),
so that the same output can later be produced for a disease without one (MVA). Nothing here modifies the engine, the
benchmarks, the protocols or any result file; engine_score / engine_score_v2 are imported and reused unchanged.

Ranking WITHOUT the disease's own gold (headline)
  The engine (arm A, protocol addendum_engine.md section 7) is refitted on the 44 OTHER benchmark-v2 diseases and applied to
  Wilms tumor. The labels of Wilms tumor (gold 1, 2, 3) are not read at all by this refit: no training row, no class
  exclusion, no calibration. A second run reuses the grouped cross-validation fold in which the disease is held out (same
  folds as out/engine/validation_v2.json, seed 20261015); that fold excludes from training the drugs of the pharmacological
  classes of the held-out diseases' gold 1, which is a protocol rule that reads gold 1 only to REMOVE training rows. It is
  reported as a check and must reproduce the per-disease AUC of validation_v2.json.
Temporal check: literature <= 2018 only, no knowledge-graph prior (arm B), labels = corrected gold 2 of the other diseases
  (as in the v2 temporal analysis; that protocol trains on post-2018 trial labels of OTHER diseases, a documented limit).

Nothing here is dosing or a treatment recommendation. The drugs are oncology drugs with serious toxicity.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from track2 import empirical_literature_control as elc
from track2 import engine_benchmark as eb
from track2 import engine_score as es
from track2 import engine_score_v2 as v2

ROOT = Path(__file__).resolve().parents[2]
TARGET = "kidney Wilms tumor"
SEED = es.SEED
LABEL = "worked case on a v2 disease, not part of the pre-registered primary validation"
TOP_N = 20
RESCUE_TOP = 100
CLASS_TOP = 5
KS = (10, 20)
BANDS = {"0_4": (0.0, 5.0), "5_11": (5.0, 12.0), "12_17": (12.0, 18.0)}  # [lo, hi) in years; unknown age is never assumed
BAND_TEXT = {"0_4": "0-4 y", "5_11": "5-11 y", "12_17": "12-17 y"}
MITOTIC = "cell_cycle_mitotic_checkpoint"
MOL_SMALL = "Small molecule"
BIOLOGIC_TYPES = {"Protein", "Antibody", "Antibody drug conjugate", "Vaccine component", "Enzyme", "Gene", "Cell", "Oligosaccharide"}
NOT_READ = "read time not recorded (deterministic output)"
HAZARD_CATEGORIES = ["boxed_warning", "toxicidad_medular", "segundos_tumores", "embrio_fetal"]
GLOBAL_CAUTION = ("Hypotheses, not advice. No dosing and no treatment recommendation. Every drug below that is an oncology drug "
                  "(ATC L01) carries serious toxicity (chemotherapy caution): bone-marrow suppression, organ toxicity, "
                  "second malignancy and embryo-fetal harm are label-documented hazards, not features.")
RUNGS = [("L1", "drug co-mentioned with a disease gene in retrieved abstracts"),
         ("L2", "at least one such paper with a sentence the rule lexicon reads as benefit"),
         ("L3", "at least one such paper with a human-patient or human-cell cue"),
         ("L4", "drug is an intervention of an EXPERIMENTAL arm of a registered phase 2+ study of the disease started after 2018"),
         ("L5", "current FDA label names the disease (openFDA indications text)")]


# ============================================================================= small helpers
def _r(x: Any, nd: int = 6) -> Optional[float]:
    return es._r(x, nd) if x is not None else None


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _name_rx(name: str) -> "re.Pattern[str]":
    return re.compile(r"(?<![\w])" + re.escape(name.lower()) + r"(?![\w])")


# ============================================================================= filter chain
def has_fda_label(name: str, labels: Sequence[dict]) -> bool:
    """openFDA label exists: at least one retrieved current label whose generic_name field names the drug as whole words."""
    rx = _name_rx(name)
    return any(rx.search(" ".join(l.get("generic_name") or []).lower()) for l in labels)


def molecule_group(mtype: Optional[str]) -> str:
    if mtype == MOL_SMALL:
        return "small_molecule"
    if mtype == "Oligonucleotide":
        return "oligonucleotide"
    if mtype in BIOLOGIC_TYPES:
        return "biologic"
    return "unknown_type"


def filter_chain(names: Sequence[str], universe: Dict[str, dict], has_label: Dict[str, bool]) -> Dict[str, Any]:
    """Universe -> FDA label -> small molecule. Counts and the dropped-by-type breakdown, all computed."""
    s1 = [n for n in names if has_label.get(n)]
    by_type: Dict[str, int] = {}
    groups = {"biologic": 0, "oligonucleotide": 0, "unknown_type": 0, "small_molecule": 0}
    for n in s1:
        t = universe[n].get("molecule_type") or "None"
        by_type[t] = by_type.get(t, 0) + 1
        groups[molecule_group(universe[n].get("molecule_type"))] += 1
    s2 = [n for n in s1 if molecule_group(universe[n].get("molecule_type")) == "small_molecule"]
    return {"n_universe": len(names), "n_fda_label": len(s1), "n_without_fda_label": len(names) - len(s1),
            "n_small_molecule": len(s2), "dropped_by_type_after_fda": {k: groups[k] for k in ("biologic", "oligonucleotide", "unknown_type")},
            "fda_label_by_molecule_type": dict(sorted(by_type.items())),
            "fda_label_names": s1, "small_molecule_names": s2}


# ============================================================================= ranking and concordance
def rank_arrays(score: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Competition ranks of descending scores: rank_min = 1 + #strictly higher, rank_max = #higher or equal."""
    s = np.sort(score)
    n = len(score)
    rmin = n - np.searchsorted(s, score, side="right") + 1
    rmax = n - np.searchsorted(s, score, side="left")
    return rmin.astype(int), rmax.astype(int)


def rank_in_mask(score: np.ndarray, mask: np.ndarray, names: Sequence[str]) -> Dict[str, int]:
    """name -> rank_min among the candidates of `mask`."""
    idx = [i for i in range(len(names)) if mask[i]]
    rmin, _ = rank_arrays(score[idx])
    return {names[i]: int(rmin[j]) for j, i in enumerate(idx)}


def top_table(score: np.ndarray, names: Sequence[str], mask: np.ndarray, n: int = TOP_N) -> List[dict]:
    """Top n of the candidates in mask by (-score, name): position, rank_min/max within the mask, score."""
    idx = [i for i in range(len(names)) if mask[i]]
    sub = score[idx]
    rmin, rmax = rank_arrays(sub)
    order = sorted(range(len(idx)), key=lambda j: (-sub[j], names[idx[j]]))
    return [{"position": k + 1, "name": names[idx[j]], "score": _r(sub[j]), "rank_min": int(rmin[j]), "rank_max": int(rmax[j]),
             "tie_group_size": int(rmax[j] - rmin[j] + 1), "_i": int(idx[j])} for k, j in enumerate(order[:n])]


def concordance(score: np.ndarray, names: Sequence[str], mask: np.ndarray, gold: np.ndarray, ks: Sequence[int] = KS) -> Dict[str, Any]:
    """Where the gold drugs rank among the candidates of `mask`. Hit at k is certain only if rank_max <= k (ties count against)."""
    idx = [i for i in range(len(names)) if mask[i]]
    sub, yy = score[idx], gold[idx]
    rmin, rmax = rank_arrays(sub)
    n = len(idx)
    drugs = []
    for j, i in enumerate(idx):
        if gold[i] == 1:
            drugs.append({"name": names[i], "rank_min": int(rmin[j]), "rank_max": int(rmax[j]),
                          "percentile_top": _r(100.0 * rmin[j] / n, 3), "tie_group_size": int(rmax[j] - rmin[j] + 1)})
    drugs.sort(key=lambda d: (d["rank_min"], d["name"]))
    absent = sorted(names[i] for i in range(len(names)) if gold[i] == 1 and not mask[i])
    out: Dict[str, Any] = {"n_candidates": n, "n_gold_total": int(gold.sum()), "n_gold_in_candidates": len(drugs), "gold_not_in_candidates": absent,
                           "drugs": drugs, "auc": _r(es.auc_mw(sub, yy)) if 0 < yy.sum() < n else None}
    for k in ks:
        out[f"hits_at_{k}_certain"] = sum(1 for d in drugs if d["rank_max"] <= k)
        out[f"expected_hits_at_{k}_tie_averaged"] = _r(es.expected_hits_at_k(sub, yy, k), 4) if yy.sum() > 0 else None
    return out


# ============================================================================= held-out refit
def heldout_scores(X: Dict[str, np.ndarray], y: Dict[str, np.ndarray], classes: np.ndarray, test_fold: Sequence[str],
                   train: Sequence[str], exclude_gold_classes: bool) -> Dict[str, Dict[str, np.ndarray]]:
    """Same nested procedure as engine_score.oof_predictions for ONE outer fold, with an explicit train set.

    exclude_gold_classes=True reproduces the protocol fold (training rows of drugs in the ATC classes of the test fold's gold 1 are
    removed; y of the test fold is read). False is the strict refit: y of the test diseases is never read."""
    train = list(train)
    bad = es._gold_classes(test_fold, y, classes) if exclude_gold_classes else set()
    keep = es._class_mask(classes, bad)
    inner_scores, inner_y = [], []
    inner = es.make_folds(train, es.SEED + 1) if len(train) >= 2 else [train]
    for ifold in inner:
        itrain = [d for d in train if d not in ifold]
        if not itrain:
            continue
        ikeep = es._class_mask(classes, es._gold_classes(ifold, y, classes))
        score, _ = es._fit_scorer([X[d] for d in itrain], [y[d] for d in itrain], [ikeep] * len(itrain))
        for d in ifold:
            inner_scores.append(score(X[d]))
            inner_y.append(y[d])
    s_in, y_in = np.concatenate(inner_scores), np.concatenate(inner_y)
    a, b = es.logit_fit(s_in, y_in)
    score, prev = es._fit_scorer([X[d] for d in train], [y[d] for d in train], [keep] * len(train))
    out = {}
    for d in test_fold:
        s = score(X[d])
        out[d] = {"score": s, "platt": 1 / (1 + np.exp(-np.clip(a + b * s, -35, 35))), "prev": np.full(len(s), prev)}
    return out


def compute_rankings(per: Dict[str, Dict[str, Any]], classes: np.ndarray, target: str) -> Dict[str, Any]:
    """Headline (strict refit) and check (protocol fold) rankings, primary arm A and temporal arm B."""
    res: Dict[str, Any] = {}
    # ---------------- primary: arm A, gold 1 endpoint set (every v2 disease has gold 1)
    sel = sorted(d for d in per if per[d]["gold1"].sum() >= 1)
    cols = es.arm_columns(True, True)
    X = {d: per[d]["primary"]["X"][:, cols] for d in sel}
    y = {d: per[d]["gold1"] for d in sel}
    others = [d for d in sel if d != target]
    y_blank = {d: y[d] for d in others}                      # the target's labels are not even passed
    strict = heldout_scores(X, y_blank, classes, [target], others, exclude_gold_classes=False)[target]
    folds = es.make_folds(sel)
    fold = next(f for f in folds if target in f)
    chk = heldout_scores(X, y, classes, fold, [d for d in sel if d not in fold], exclude_gold_classes=True)[target]
    res["primary"] = {"strict": strict, "fold": chk, "folds": folds, "fold_of_target": fold, "n_train_strict": len(others),
                      "auc_fold_vs_gold1": es.auc_mw(chk["score"], per[target]["gold1"]),
                      "auc_strict_vs_gold1": es.auc_mw(strict["score"], per[target]["gold1"])}
    # ---------------- temporal: arm B (no prior), literature <= 2018, labels = corrected gold 2
    sel2 = sorted(d for d in per if per[d]["gold2"].sum() >= 1)
    if target in sel2 and len(sel2) >= 3:
        cols_t = es.arm_columns(False, True)
        Xt = {d: per[d]["temporal"]["X"][:, cols_t] for d in sel2}
        yt = {d: per[d]["gold2"] for d in sel2}
        others2 = [d for d in sel2 if d != target]
        strict_t = heldout_scores(Xt, {d: yt[d] for d in others2}, classes, [target], others2, exclude_gold_classes=False)[target]
        folds_t = es.make_folds(sel2)
        fold_t = next(f for f in folds_t if target in f)
        chk_t = heldout_scores(Xt, yt, classes, fold_t, [d for d in sel2 if d not in fold_t], exclude_gold_classes=True)[target]
        res["temporal"] = {"strict": strict_t, "fold": chk_t, "folds": folds_t, "fold_of_target": fold_t, "n_train_strict": len(others2),
                           "auc_fold_vs_gold2": es.auc_mw(chk_t["score"], per[target]["gold2"]),
                           "auc_strict_vs_gold2": es.auc_mw(strict_t["score"], per[target]["gold2"])}
    else:
        res["temporal"] = None
    return res


# ============================================================================= age bands from label status
def band_status(status: str, min_age: Optional[float]) -> Dict[str, str]:
    """Label support per band 0-4, 5-11, 12-17. status: ESTABLISHED / NOT_ESTABLISHED / anything else (-> verify).

    Established with a stated minimum age m: labelled if m <= band start, `partly (from m y)` if the band start < m < band end,
    not_labelled if m >= band end. Established without a stated age: age_not_stated. The age of a patient is never assumed."""
    out = {}
    for b, (lo, hi) in BANDS.items():
        if status == "NOT_ESTABLISHED":
            out[b] = "not_labelled"
        elif status == "ESTABLISHED":
            if min_age is None:
                out[b] = "age_not_stated"
            elif min_age <= lo:
                out[b] = "labelled"
            elif min_age < hi:
                out[b] = f"partly (from {min_age:g} y)"
            else:
                out[b] = "not_labelled"
        else:
            out[b] = "verify"
    return out


def status_class(s: str) -> str:
    return "partly" if s.startswith("partly") else s


def label_support(drug: dict, fetch: Callable[[str], Tuple[Optional[dict], List[str]]], review: dict, row_for: Callable) -> Dict[str, Any]:
    """Label-based age support of one drug, reusing age_label_table.row_for (heuristic classifier + manual review in config)."""
    lab, failures = fetch(drug["name"])
    row = row_for(drug, ["engine_case_study"], lab, failures, review, NOT_READ)
    manual = row.get("revision_manual")
    if row.get("etiqueta") is None:
        status, min_age, basis = "UNREADABLE", None, "no label readable; nothing is inferred"
    elif manual and manual.get("verdict") == "established":
        status, min_age, basis = "ESTABLISHED", manual.get("min_age_years"), "manual review (config/pediatric_label_review.json)"
    elif manual and manual.get("verdict") == "not_established":
        status, min_age, basis = "NOT_ESTABLISHED", None, "manual review (config/pediatric_label_review.json)"
    elif manual:
        status, min_age, basis = "VERIFY", None, "manual review verdict: " + str(manual.get("verdict"))
    else:
        auto = (row.get("clasificador_automatico") or {})
        status, min_age, basis = auto.get("estado") or "UNKNOWN", auto.get("edad_minima_anios"), "heuristic classifier of labels.classify_pediatric (not hand-verified)"
        if not (lab or {}).get("sections", {}).get("pediatric_use"):
            status = "UNKNOWN"
            basis += "; no readable paediatric-use section"
    bands = band_status(status, min_age)
    quotes = []
    if lab:
        from track2 import caution_list as cl
        ref = {"source": lab["source"], "setid": lab["setid"], "title": lab["title"]}
        quotes = cl.extract_quotes(lab["sections"], [c for c in HAZARD_CATEGORIES if c in cl.CATEGORIES], NOT_READ, ref)
    return {"name": drug["name"], "drug_id": drug["drug_id"], "label_status": status, "minimum_age_years": min_age, "basis": basis,
            "bands": bands, "manual_review_required": bool(row.get("manual_review_required", True)),
            "scope_of_manual_review": (manual or {}).get("scope"),
            "label": row.get("etiqueta"), "reading_failures": failures or None,
            "verbatim_pediatric_use": row.get("frase_literal_sustento"),
            "verbatim_pediatric_warnings": row.get("advertencias_pediatricas") or [],
            "verbatim_hazard_quotes": [{"category": q["categoria"], "section": q["seccion"], "quote": q["texto_literal"], "truncated": q["truncada"],
                                         "set_id": q["setid"]} for q in quotes],
            "three_group_status_existing_pipeline": row.get("estado_por_grupo")}


# ============================================================================= evidence rungs and gaps
def evidence_rungs(lit: Dict[str, int], in_gold2: bool, in_gold3: bool) -> Dict[str, bool]:
    return {"L1": lit.get("n_papers", 0) >= 1, "L2": lit.get("n_benefit", 0) >= 1, "L3": lit.get("n_human", 0) >= 1,
            "L4": bool(in_gold2), "L5": bool(in_gold3)}


def evidence_gap(name: str, atc: Sequence[str], rungs: Dict[str, bool], bands: Dict[str, str], trials: Sequence[str]) -> Dict[str, Any]:
    reached = [r for r, _ in RUNGS if rungs.get(r)]
    top = reached[-1] if reached else "L0"
    missing = [{"rung": r, "meaning": t} for r, t in RUNGS if not rungs.get(r)]
    oncology = any(a.startswith("L01") for a in atc)
    gaps = [b for b, s in bands.items() if status_class(s) in ("verify", "not_labelled", "age_not_stated")]
    closers: List[str] = []
    if not rungs.get("L4") and not rungs.get("L5"):
        closers.append("a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, "
                       "with pre-specified response and toxicity endpoints and an independent safety committee")
    if rungs.get("L4") and not rungs.get("L5"):
        closers.append("completion and reporting of the registered trial(s) " + ", ".join(trials) + "; a label supplement for the disease needs the sponsor's "
                       "regulatory submission, which is not something this analysis can supply")
    if rungs.get("L5") and gaps:
        closers.append("an age-stratified label supplement (paediatric study data for the bands without labelled support)")
    if not rungs.get("L5") and not oncology and not rungs.get("L4"):
        closers.append("a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review")
    if rungs.get("L3") and not rungs.get("L4"):
        closers.append("a registry of exposure and outcome inside an existing paediatric oncology cooperative-group protocol, to collect the "
                       "human data that abstracts cannot provide")
    if oncology and not closers:
        closers.append("a registry of exposure and outcome inside an existing paediatric oncology protocol")
    return {"name": name, "evidence_level_reached": top, "rungs": rungs, "missing": missing, "age_band_gaps": gaps,
            "single_study_or_data_that_would_close_the_gap": closers[0] if closers else "none identified beyond the above",
            "other_options": closers[1:], "oncology_drug_chemotherapy_caution": oncology}


# ============================================================================= data assembly
def build_disease_features(dis: dict, names: Sequence[str], index: Dict[str, int], dre: "re.Pattern[str]", gate: Set[str],
                           targets: dict, ppi: dict, pops: Dict[str, Tuple[dict, int]], cache: Path, getter: Optional[Callable]) -> Dict[str, Any]:
    symbols = sorted({g["symbol"] for g in dis["genes"]})
    rec: Dict[str, Any] = {"symbols": symbols}
    rec.update(v2.gold_vectors(dis, index, len(names)))
    for period, year in (("primary", es.YEAR_PRIMARY), ("temporal", es.YEAR_TEMPORAL)):
        pop, tot = pops[period]
        papers = es.fetch_papers(symbols, year, cache, getter)
        lit = es.literature_counts(papers, symbols, dre, gate)
        prior = es.prior_features(symbols, names, targets, ppi) if period == "primary" else None
        rec[period] = {"X": es.disease_matrix(lit, names, pop, tot, prior), "lit": lit, "n_papers_retrieved": len(papers),
                       "papers": papers if dis.get("name") == TARGET else None}
    return rec


def class_comention(papers: Sequence[dict], genes: Sequence[str], dre: "re.Pattern[str]", gate: Set[str]) -> Dict[str, Dict[str, int]]:
    """Per drug: independent papers that mention the drug and (case-sensitively) a gene of the class, with benefit / human papers."""
    return es.literature_counts(papers, genes, dre, gate)["counts"]


def class_top_drugs(counts: Dict[str, Dict[str, int]], allowed: Set[str], n: int = CLASS_TOP) -> List[dict]:
    rows = [{"name": k, "n_papers": v["n_papers"], "n_benefit": v["n_benefit"], "n_human": v["n_human"]} for k, v in counts.items()
            if k in allowed and v["n_papers"] >= 1]
    rows.sort(key=lambda r: (-r["n_papers"], -r["n_benefit"], r["name"]))
    return rows[:n]


def age_table(class_rows: Dict[str, List[dict]], supports: Dict[str, dict]) -> Dict[str, Dict[str, Dict[str, List[str]]]]:
    """class -> band -> status class -> drug names (only drugs whose label support was read)."""
    out: Dict[str, Dict[str, Dict[str, List[str]]]] = {}
    for c, rows in class_rows.items():
        out[c] = {}
        for b in BANDS:
            cell: Dict[str, List[str]] = {}
            for r in rows:
                s = supports.get(r["name"])
                st = status_class(s["bands"][b]) if s else "no_label_read"
                cell.setdefault(st, []).append(r["name"])
            out[c][b] = {k: sorted(v) for k, v in sorted(cell.items())}
    return out


# ============================================================================= driver
def run(root: Path = ROOT, getter: Optional[Callable] = None, fetch: Optional[Callable] = None, lookup: Optional[Callable] = None,
        out_json: Optional[Path] = None, out_md: Optional[Path] = None, target: str = TARGET) -> Dict[str, Any]:
    from track2 import age_label_table as alt
    from track2 import caution_list as cl
    from track2 import engine_benchmark_v2 as b2mod
    from track2 import fda_approved
    from track2.pediatric_review import load_review

    bench2_path = root / "out" / "engine" / "benchmark_v2.json"
    val_path = root / "out" / "engine" / "validation_v2.json"
    bench2 = json.loads(bench2_path.read_text(encoding="utf-8"))
    val = json.loads(val_path.read_text(encoding="utf-8"))
    dis = next(d for d in bench2["diseases"] if d["name"] == target)

    drugs = eb.load_universe(root / "data" / "universe.json")
    universe = eb.universe_index(drugs)
    names = sorted(universe)
    index = {n: i for i, n in enumerate(names)}
    dre = elc.drug_regex(names)
    gate = es.name_gate(names)
    classes = np.array(es.drug_classes(universe, names))
    guard = v2.RequestGuard(getter or v2.epmc_getter, None)
    cdir = root / "data" / "cache_engine"
    targets, ppi = es.load_prior_graph(root, cdir / "primekg_prior_edges.json")
    pops = {"primary": es.fetch_popularity(names, None, cdir / "epmc_popularity_all.json", guard),
            "temporal": es.fetch_popularity(names, es.YEAR_TEMPORAL, cdir / "epmc_popularity_2018.json", guard)}
    key = sha256_file(bench2_path)[:16] + sha256_file(root / "data" / "universe.json")[:16] + "|" + target
    fcache = cdir / "case_study_features.pkl"
    per = None
    if fcache.exists():
        import pickle
        with fcache.open("rb") as fh:   # own cache file written by this module
            blob = pickle.load(fh)
        per = blob["per"] if blob.get("key") == key else None
    if per is None:
        per = {d["name"]: build_disease_features(d, names, index, dre, gate, targets, ppi, pops, cdir / "epmc", guard) for d in bench2["diseases"]}
        import pickle
        with fcache.open("wb") as fh:
            pickle.dump({"key": key, "per": per}, fh)
    W = per[target]

    # ------------------------------------------------------------------ 1. rankings
    rk = compute_rankings(per, classes, target)
    chk_auc = val["primary_sets"]["gold1"]["all_45"]["arms"]["A"]["per_disease_auc"].get(target)
    chk_auc_t = (((val["temporal_gold2_corrected"]["all_45"].get("all_gold2") or {}).get("arms") or {}).get("B") or {}).get("per_disease_auc", {}).get(target)
    stored_folds = val["primary_sets"]["gold1"]["all_45"]["folds"]
    stored_folds_t = ((val["temporal_gold2_corrected"]["all_45"].get("all_gold2") or {}).get("folds"))
    repro = {"primary_fold_auc_recomputed": _r(rk["primary"]["auc_fold_vs_gold1"]), "primary_fold_auc_in_validation_v2": chk_auc,
             "primary_fold_auc_identical": chk_auc is not None and _r(rk["primary"]["auc_fold_vs_gold1"]) == chk_auc,
             "primary_folds_identical_to_validation_v2": rk["primary"]["folds"] == stored_folds,
             "primary_strict_refit_auc_vs_gold1": _r(rk["primary"]["auc_strict_vs_gold1"])}
    if rk["temporal"]:
        repro.update({"temporal_fold_auc_recomputed": _r(rk["temporal"]["auc_fold_vs_gold2"]), "temporal_fold_auc_in_validation_v2_armB": chk_auc_t,
                      "temporal_fold_auc_identical": chk_auc_t is not None and _r(rk["temporal"]["auc_fold_vs_gold2"]) == chk_auc_t,
                      "temporal_folds_identical_to_validation_v2": rk["temporal"]["folds"] == stored_folds_t,
                      "temporal_strict_refit_auc_vs_gold2": _r(rk["temporal"]["auc_strict_vs_gold2"])})

    # ------------------------------------------------------------------ 2. filter chain (FDA label or Drugs@FDA, small molecule)
    sc = rk["primary"]["strict"]["score"]
    pl = rk["primary"]["strict"]["platt"]
    lc = W["primary"]["lit"]["counts"]
    gold = {"gold1": W["gold1"], "gold3": W["gold3"], "gold2": W["gold2"]}
    gold_names = {k: sorted(names[i] for i in np.flatnonzero(v)) for k, v in gold.items()}
    label_only = {n: has_fda_label(n, b2mod.fetch_labels(n, getter=_no_network)) for n in names}
    lookup_fn = lookup or fda_approved.lookup
    # The cached openFDA label query misses some approved drugs (for example brand-only substance names). For the drugs where this matters
    # (no label match AND (top RESCUE_TOP of the universe ranking OR a known drug of any gold)) Drugs@FDA is asked; an error is never read as eligible.
    order_all = sorted(range(len(names)), key=lambda i: (-sc[i], names[i]))
    check_names = sorted({names[i] for i in order_all[:RESCUE_TOP] if not label_only[names[i]]} |
                         {n for v in gold_names.values() for n in v if not label_only[n]})
    rescue: Dict[str, dict] = {}
    for n in check_names:
        try:
            rescue[n] = lookup_fn(universe[n]["name"])
        except Exception as e:  # noqa: BLE001
            rescue[n] = {"eligible": None, "error": repr(e)[:200]}
    has_label = {n: label_only[n] or bool((rescue.get(n) or {}).get("eligible")) for n in names}
    chain = filter_chain(names, universe, has_label)
    chain["n_fda_label_only_rule"] = sum(label_only.values())
    chain["n_added_by_drugs_at_fda_check"] = sum(1 for n in names if has_label[n] and not label_only[n])
    chain["n_checked_at_drugs_at_fda"] = len(check_names)
    chain["added_by_drugs_at_fda_check"] = sorted(n for n in names if has_label[n] and not label_only[n])
    m_all = np.array([has_label[n] for n in names])
    sm_names_set = set(chain["small_molecule_names"])
    m_sm = np.array([n in sm_names_set for n in names])
    dropped_gold = {k: [{"name": n, "has_fda_label_or_drugs_at_fda": has_label[n], "molecule_type": universe[n].get("molecule_type"),
                         "group": molecule_group(universe[n].get("molecule_type"))} for n in v if not (has_label[n] and molecule_group(universe[n].get("molecule_type")) == "small_molecule")]
                    for k, v in gold_names.items()}

    def decorate(rows: List[dict], pl_arr: np.ndarray) -> List[dict]:
        for r in rows:
            i = r.pop("_i")
            n = r["name"]
            r.update({"calibrated_probability_platt": _r(pl_arr[i], 6), "molecule_type": universe[n].get("molecule_type"), "fda_label": has_label[n],
                      "atc": universe[n].get("atc") or [], "first_approval": universe[n].get("first_approval"),
                      "in_gold1": bool(gold["gold1"][i]), "in_gold2": bool(gold["gold2"][i]), "in_gold3": bool(gold["gold3"][i]),
                      "n_papers_with_gene": lc.get(n, {}).get("n_papers", 0), "n_benefit": lc.get(n, {}).get("n_benefit", 0),
                      "n_human": lc.get(n, {}).get("n_human", 0), "n_harm_sentences": lc.get(n, {}).get("n_harm", 0)})
        return rows

    all_mask = np.ones(len(names), dtype=bool)
    top_universe = decorate(top_table(sc, names, all_mask), pl)
    top_filtered = decorate(top_table(sc, names, m_sm), pl)
    for r in top_universe:
        r["passes_fda_label"] = r["fda_label"]
        r["passes_small_molecule"] = bool(m_sm[index[r["name"]]])

    concord = {k: {"primary_literature_to_2024": {"approved_fda_label": concordance(sc, names, m_all, v),
                                                  "approved_small_molecule": concordance(sc, names, m_sm, v)}}
               for k, v in gold.items()}
    concord["strict_refit_vs_protocol_fold_gold1"] = {
        "fold_approved_fda_label": concordance(rk["primary"]["fold"]["score"], names, m_all, gold["gold1"]),
        "fold_approved_small_molecule": concordance(rk["primary"]["fold"]["score"], names, m_sm, gold["gold1"])}

    # ------------------------------------------------------------------ 3. temporal check
    temporal: Optional[Dict[str, Any]] = None
    if rk["temporal"]:
        st = rk["temporal"]["strict"]["score"]
        lct = W["temporal"]["lit"]["counts"]
        later = np.array([(universe[n].get("first_approval") or 0) > es.YEAR_TEMPORAL for n in names])
        tt = []
        for tm, mname in ((m_all, "approved_fda_label"), (m_sm, "approved_small_molecule")):
            rows = top_table(st, names, tm)
            for r in rows:
                i = r.pop("_i")
                n = r["name"]
                r.update({"in_gold1": bool(gold["gold1"][i]), "in_gold2_later_trialled": bool(gold["gold2"][i]),
                          "first_approval": universe[n].get("first_approval"), "first_approved_after_2018": bool(later[i]),
                          "n_papers_with_gene_to_2018": lct.get(n, {}).get("n_papers", 0)})
            tt.append((mname, rows))
        temporal = {"literature": "Europe PMC records with PUB_YEAR <= 2018 only; knowledge-graph prior excluded (arm B)",
                    "labels_for_training": "corrected gold 2 of the other diseases (v2 temporal protocol; post-2018 trials of OTHER diseases)",
                    "top20": {m: rows for m, rows in tt},
                    "gold2_later_trialled": {mname: concordance(st, names, tm, gold["gold2"]) for tm, mname in ((m_all, "approved_fda_label"), (m_sm, "approved_small_molecule"))},
                    "gold1_known_drugs": {mname: concordance(st, names, tm, gold["gold1"]) for tm, mname in ((m_all, "approved_fda_label"), (m_sm, "approved_small_molecule"))},
                    "later_approved_first_approval_after_2018": {
                        mname: [{"name": names[i], "first_approval": universe[names[i]].get("first_approval"), "rank_min": rank_in_mask(st, tm, names)[names[i]],
                                 "in_gold2": bool(gold["gold2"][i])} for i in np.flatnonzero(later & tm)
                                if lct.get(names[i], {}).get("n_papers", 0) > 0 or gold["gold2"][i]]
                        for tm, mname in ((m_all, "approved_fda_label"), (m_sm, "approved_small_molecule"))},
                    "n_universe_drugs_first_approved_after_2018": int(later.sum())}

    # ------------------------------------------------------------------ 4. label support (top filtered + class tops + gold + mitotic mechanism)
    genes_by_class: Dict[str, List[str]] = {}
    for g in dis["genes"]:
        genes_by_class.setdefault(g["mechanism_class_general"], []).append(g["symbol"])
    sm_set = set(chain["small_molecule_names"])
    papers = W["primary"]["papers"] or []
    class_rows: Dict[str, List[dict]] = {}
    class_blocks: Dict[str, Any] = {}
    for c in sorted(genes_by_class):
        cnt = class_comention(papers, genes_by_class[c], dre, gate)
        n_gene_papers = es.literature_counts(papers, genes_by_class[c], dre, gate)["n_papers_with_gene"]
        rows = class_top_drugs(cnt, sm_set)
        rmin, _ = rank_arrays(sc[[i for i in range(len(names)) if m_sm[i]]])
        pos = {names[i]: k for k, i in enumerate([i for i in range(len(names)) if m_sm[i]])}
        for r in rows:
            r["engine_rank_small_molecule_min"] = int(rmin[pos[r["name"]]])
            r["in_gold1"] = bool(gold["gold1"][index[r["name"]]])
            r["in_gold2"] = bool(gold["gold2"][index[r["name"]]])
        class_rows[c] = rows
        class_blocks[c] = {"genes": sorted(genes_by_class[c]), "n_retrieved_papers_mentioning_class_gene": n_gene_papers,
                           "top_drugs_by_comention_small_molecule_fda": rows,
                           "scope_note": "counts are independent papers among the retrieved abstracts (<= 300, query = all 14 disease genes) that name the drug and a gene of this class"}
    mito_mech = [m for m in cl.classify_universe(drugs)["b_checkpoint_mitotico"]]
    mito_rows = []
    for m in mito_mech:
        n = m["name"].lower()
        i = index.get(n)
        if i is None:
            continue
        mito_rows.append({"name": n, "mechanism_note": m["motivo"], "has_fda_label": has_label[n], "molecule_type": universe[n].get("molecule_type"),
                          "small_molecule_fda": bool(m_sm[i]), "engine_rank_all_min": int(rank_arrays(sc)[0][i]),
                          "n_papers_with_wilms_gene": lc.get(n, {}).get("n_papers", 0)})
    mito_rows.sort(key=lambda r: (r["engine_rank_all_min"], r["name"]))
    rmin_all = rank_arrays(sc)[0]
    tub_rows = []
    for m in cl.classify_universe(drugs)["a_tubulina"]:
        n = m["name"].lower()
        i = index.get(n)
        if i is None:
            continue
        tub_rows.append({"name": n, "criterion": m["motivo"], "has_fda_label": has_label[n], "molecule_type": universe[n].get("molecule_type"),
                         "small_molecule_fda": bool(m_sm[i]), "engine_rank_all_min": int(rmin_all[i]),
                         "n_papers_with_wilms_gene": lc.get(n, {}).get("n_papers", 0), "in_gold1": bool(gold["gold1"][i]),
                         "in_gold2": bool(gold["gold2"][i])})
    tub_rows.sort(key=lambda r: (r["engine_rank_all_min"], r["name"]))
    known = sorted({n for k in ("gold1", "gold3") for n in gold_names[k] if has_label[n]})

    need: List[str] = [r["name"] for r in top_filtered]
    for c in class_rows:
        need += [r["name"] for r in class_rows[c]]
    need += known
    need += [r["name"] for r in mito_rows if r["small_molecule_fda"]]
    need += [r["name"] for r in tub_rows if r["small_molecule_fda"]]
    need = sorted(set(need))
    review = load_review(root / "config" / "pediatric_label_review.json")
    setids: Dict[str, str] = {}
    for e in review.values():
        m = re.search(r"setid=([0-9a-f-]{36})", e.get("source") or "")
        if m:
            setids[e["name"].upper()] = m.group(1)

    def default_fetch(name: str) -> Tuple[Optional[dict], List[str]]:
        if name.upper() in setids:
            lab, failures = cl.fetch_by_setid(setids[name.upper()])
            if lab:
                return lab, failures
        return cl.fetch_label(name)

    fetch_fn = fetch or default_fetch
    supports = {n: label_support(universe[n], fetch_fn, review, alt.row_for) for n in need}
    approvals = {}
    for n in sorted({r["name"] for r in top_filtered} | {n for k in ("gold1", "gold3") for n in gold_names[k]}):
        if n in rescue:
            approvals[n] = rescue[n]
            continue
        try:
            approvals[n] = lookup_fn(universe[n]["name"])
        except Exception as e:  # noqa: BLE001 - a service error is never read as eligible / not eligible
            approvals[n] = {"eligible": None, "error": repr(e)[:200]}

    # ------------------------------------------------------------------ 5. evidence gaps
    trials_by_drug = {g["name"].lower(): sorted({t["nct"] for t in g.get("trials", [])}) for g in dis["gold2"]}
    gaps = []
    for r in top_filtered:
        n = r["name"]
        rungs = evidence_rungs(lc.get(n, {}), r["in_gold2"], r["in_gold3"])
        g = evidence_gap(n, universe[n].get("atc") or [], rungs, supports[n]["bands"], trials_by_drug.get(n, []))
        g["position_in_filtered_ranking"] = r["position"]
        g["hazards_label_documented"] = supports[n]["verbatim_hazard_quotes"] or "no hazard sentence retrieved from the label sections read: not read as safe"
        g["trials_gold2"] = trials_by_drug.get(n, [])
        gaps.append(g)

    result = {
        "label": LABEL, "disease": target, "seed": SEED, "caution": GLOBAL_CAUTION,
        "inputs": {"benchmark_v2_sha256": sha256_file(bench2_path), "validation_v2_sha256": sha256_file(val_path),
                   "universe_sha256": sha256_file(root / "data" / "universe.json"), "chembl_version": json.loads((root / "data" / "universe.json").read_text(encoding="utf-8"))["chembl_version"],
                   "n_diseases_in_benchmark_v2": len(per), "gold_sizes": dis["gold_sizes"], "genes": [g["symbol"] for g in dis["genes"]],
                   "n_papers_retrieved_primary": W["primary"]["n_papers_retrieved"], "n_papers_with_gene_primary": W["primary"]["lit"]["n_papers_with_gene"],
                   "gold_names": gold_names},
        "method": {"headline_ranking": "arm A (prior + popularity correction), refit on the 44 other benchmark-v2 diseases (strict: no label of the target is read), "
                                       "literature <= 2024, applied to the target for all %d universe drugs" % len(names),
                   "check_ranking": "protocol grouped-CV fold in which the disease is held out (same folds as validation_v2.json)",
                   "fold_of_target": rk["primary"]["fold_of_target"],
                   "folds_primary": rk["primary"]["folds"],
                   "temporal_fold_of_target": rk["temporal"]["fold_of_target"] if rk["temporal"] else None,
                   "approval_rule": "FDA-approved = a current openFDA drug label exists whose generic_name field names the drug as whole words "
                                    "(cached labels of the gold 3 retrieval, up to 5 per drug); Drugs@FDA confirmation (marketing status Prescription) "
                                    "by track2.fda_approved.lookup is added for the displayed candidates",
                   "small_molecule_rule": "ChEMBL molecule_type == 'Small molecule' in data/universe.json (as recorded; no override)",
                   "age_bands": {BAND_TEXT[b]: list(v) for b, v in BANDS.items()},
                   "age_rule": "label text only; established use with a stated minimum age m is labelled for a band if m <= band start, partly if inside the band, "
                               "not labelled if m >= band end; established without an age = age_not_stated; unknown/mixed = verify; the patient's age is never assumed",
                   "evidence_rungs": [{"rung": r, "meaning": t} for r, t in RUNGS]},
        "reproduction_checks": repro,
        "ranking_top20_universe": top_universe,
        "ranking_top20_fda_small_molecule": top_filtered,
        "filter_chain": {k: v for k, v in chain.items() if k not in ("fda_label_names", "small_molecule_names")},
        "gold_drugs_dropped_by_filters": dropped_gold,
        "drugs_at_fda_confirmation": approvals,
        "concordance": concord,
        "temporal_check": temporal,
        "mechanism_classes": class_blocks,
        "mitotic_checkpoint_class": {"class": MITOTIC, "genes": sorted(genes_by_class.get(MITOTIC, [])),
                                     "literature_block": class_blocks.get(MITOTIC),
                                     "universe_drugs_with_mitotic_checkpoint_mechanism": mito_rows,
                                     "n_universe_drugs_with_mitotic_checkpoint_mechanism": len(mito_rows),
                                     "tubulin_binders_acting_on_the_mitotic_spindle": tub_rows,
                                     "note": "mechanism list = caution_list class b (ChEMBL mechanism: Aurora, PLK, TTK/MPS1, BUB1, MAD2, mitotic checkpoint, spindle assembly)"},
        "mechanism_x_age_band": age_table(class_rows, supports),
        "known_drugs_label_support": {n: supports[n] for n in known},
        "label_support": supports,
        "evidence_gaps": gaps,
        "limits": ["Label support is label text read by a heuristic classifier unless a manual review exists in config; manual_review_required marks the former.",
                   "Co-mention counts are among at most 300 retrieved abstracts for 14 genes together; they are not disease-specific effect evidence.",
                   "ChEMBL molecule type is used as recorded (for example dactinomycin is typed Protein and is therefore dropped by the small-molecule filter).",
                   "openFDA label existence uses current labels by generic-name field and may include combination or homeopathic products that name the substance.",
                   "Wilms tumour is a v2 disease; this is a worked case, not part of the pre-registered primary validation, and no threshold is evaluated here."],
    }
    out_json = out_json or root / "out" / "engine" / "case_wilms.json"
    out_md = out_md or root / "out" / "engine" / "case_wilms.md"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False, default=v2._default) + "\n", encoding="utf-8")
    out_md.write_text(render_md(result), encoding="utf-8")
    return result


def _no_network(url: str, params: dict) -> dict:
    raise RuntimeError("openFDA label cache miss (the case study reads the cached gold 3 retrieval only)")


# ============================================================================= markdown
def _q(t: Optional[str], n: int = 260) -> str:
    if not t:
        return ""
    t = t.replace("|", "/").replace("\n", " ")
    return t if len(t) <= n else t[:n] + " [...]"


def render_md(r: Dict[str, Any]) -> str:
    L: List[str] = []
    a = L.append
    a(f"# Worked case: {r['disease']}")
    a("")
    a(f"**{r['label']}.**")
    a("")
    a(f"> {r['caution']}")
    a("")
    a(f"Seed {r['seed']}. Inputs: benchmark_v2 sha256 `{r['inputs']['benchmark_v2_sha256'][:16]}`, universe sha256 `{r['inputs']['universe_sha256'][:16]}` ({r['inputs']['chembl_version']}). "
      f"Genes ({len(r['inputs']['genes'])}): {', '.join(r['inputs']['genes'])}. Known drugs: gold 1 = {', '.join(r['inputs']['gold_names']['gold1'])}; "
      f"gold 3 = {', '.join(r['inputs']['gold_names']['gold3'])}; corrected gold 2 has {len(r['inputs']['gold_names']['gold2'])} drugs.")
    a("")
    a("## 1. Engine ranking without the disease's own gold")
    a("")
    a(f"Headline: {r['method']['headline_ranking']}. Check: {r['method']['check_ranking']} (the fold of {len(r['method']['fold_of_target'])} diseases that contains the target).")
    rc = r["reproduction_checks"]
    a(f"Reproduction of validation_v2.json: fold AUC vs gold 1 = {rc['primary_fold_auc_recomputed']} (stored {rc['primary_fold_auc_in_validation_v2']}; identical: {rc['primary_fold_auc_identical']}; "
      f"folds identical: {rc['primary_folds_identical_to_validation_v2']}). Strict refit AUC vs gold 1 = {rc['primary_strict_refit_auc_vs_gold1']}.")
    a("")
    a("Top 20 of all universe drugs (score = logistic decision value; probability = Platt-calibrated):")
    a("")
    a("| # | drug | score | probability | FDA label | type | gold 1 | gold 2 | gold 3 | papers with gene | benefit | human |")
    a("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for x in r["ranking_top20_universe"]:
        a(f"| {x['position']} | {x['name']} | {x['score']} | {x['calibrated_probability_platt']} | {'yes' if x['fda_label'] else 'no'} | {x['molecule_type']} | "
          f"{'x' if x['in_gold1'] else ''} | {'x' if x['in_gold2'] else ''} | {'x' if x['in_gold3'] else ''} | {x['n_papers_with_gene']} | {x['n_benefit']} | {x['n_human']} |")
    a("")
    a("## 2. Filter chain and label-based age support")
    a("")
    fc = r["filter_chain"]
    a(f"Universe {fc['n_universe']} -> FDA approved {fc['n_fda_label']} (label rule alone: {fc['n_fda_label_only_rule']}; added by the Drugs@FDA check of {fc['n_checked_at_drugs_at_fda']} drugs: "
      f"{', '.join(fc['added_by_drugs_at_fda_check']) or 'none'}; not approved by either rule: {fc['n_without_fda_label']}) -> small molecule {fc['n_small_molecule']}. "
      f"Dropped by type after the FDA step: {fc['dropped_by_type_after_fda']['biologic']} biologics, {fc['dropped_by_type_after_fda']['oligonucleotide']} oligonucleotides, "
      f"{fc['dropped_by_type_after_fda']['unknown_type']} of unknown type.")
    a("")
    for k, v in r["gold_drugs_dropped_by_filters"].items():
        if v:
            a(f"- Known drugs of {k} dropped by the filters: " + "; ".join(f"{d['name']} (label: {d['has_fda_label_or_drugs_at_fda']}, ChEMBL type {d['molecule_type']})" for d in v) + ".")
    a("")
    a("Top 20 after the filters (FDA label, small molecule):")
    a("")
    a("| # | drug | score | Drugs@FDA | 0-4 y | 5-11 y | 12-17 y | basis |")
    a("|---|---|---|---|---|---|---|---|")
    for x in r["ranking_top20_fda_small_molecule"]:
        s = r["label_support"][x["name"]]
        ap = r["drugs_at_fda_confirmation"].get(x["name"], {})
        a(f"| {x['position']} | {x['name']} | {x['score']} | {ap.get('eligible')} | {s['bands']['0_4']} | {s['bands']['5_11']} | {s['bands']['12_17']} | {_q(s['basis'], 60)} |")
    a("")
    a("Verbatim label sentences (section pediatric_use unless stated; not dosing information):")
    a("")
    for x in r["ranking_top20_fda_small_molecule"]:
        s = r["label_support"][x["name"]]
        v = s["verbatim_pediatric_use"]
        if v:
            a(f"- **{x['name']}** ({v['fuente']}, set id {v['setid']}, section {v['seccion']}): \"{_q(v['texto_literal'], 420)}\"")
        else:
            a(f"- **{x['name']}**: no readable paediatric-use sentence ({s['basis']}); status {s['label_status']}")
    a("")
    a("Known drugs (gold 1 / gold 3) with a current FDA label, label-based age support (any molecule type; not filtered):")
    a("")
    a("| drug | 0-4 y | 5-11 y | 12-17 y | basis | verbatim paediatric-use sentence |")
    a("|---|---|---|---|---|---|")
    for n, s_ in r["known_drugs_label_support"].items():
        v = s_["verbatim_pediatric_use"]
        a(f"| {n} | {s_['bands']['0_4']} | {s_['bands']['5_11']} | {s_['bands']['12_17']} | {_q(s_['basis'], 70)} | "
          + (f"\"{_q(v['texto_literal'], 300)}\" ({v['fuente']}, set id {v['setid']}, section {v['seccion']})" if v else "none readable") + " |")
    if any(s_["scope_of_manual_review"] for s_ in r["known_drugs_label_support"].values()):
        a("")
        for n, s_ in r["known_drugs_label_support"].items():
            if s_["scope_of_manual_review"]:
                a(f"- Manual-review scope for {n}: \"{_q(s_['scope_of_manual_review'], 300)}\"")
    a("")
    a("## 3. Concordance with the known drugs")
    a("")
    for g in ("gold1", "gold3", "gold2"):
        a(f"### {g}")
        for sub, label in (("approved_fda_label", "among approved (FDA label) candidates"), ("approved_small_molecule", "among approved small molecules")):
            c = r["concordance"][g]["primary_literature_to_2024"][sub]
            a("")
            a(f"{label}: {c['n_candidates']} candidates, {c['n_gold_in_candidates']} of {c['n_gold_total']} known drugs present, AUC {c['auc']}; certain hits at 10: {c['hits_at_10_certain']}, at 20: {c['hits_at_20_certain']} "
              f"(tie-averaged expected {c['expected_hits_at_10_tie_averaged']} / {c['expected_hits_at_20_tie_averaged']}).")
            if c["drugs"]:
                a("")
                a("| drug | rank (min-max) | top % |")
                a("|---|---|---|")
                for d in c["drugs"][:25]:
                    a(f"| {d['name']} | {d['rank_min']}-{d['rank_max']} | {d['percentile_top']} |")
            if c["gold_not_in_candidates"]:
                a("")
                a("Not in the candidate set: " + ", ".join(c["gold_not_in_candidates"]))
        a("")
    t = r["temporal_check"]
    a("### Temporal check (literature <= 2018, no prior)")
    a("")
    if t:
        a(f"{t['literature']}. Labels: {t['labels_for_training']}. Drugs first approved after 2018 in the universe: {t['n_universe_drugs_first_approved_after_2018']}.")
        a("")
        for m in ("approved_fda_label", "approved_small_molecule"):
            g2 = t["gold2_later_trialled"][m]
            g1 = t["gold1_known_drugs"][m]
            a(f"- {m}: later-trialled (gold 2) drugs present {g2['n_gold_in_candidates']}, AUC {g2['auc']}, certain hits at 10/20: {g2['hits_at_10_certain']}/{g2['hits_at_20_certain']}; "
              f"gold 1 drugs: " + ", ".join(f"{d['name']} #{d['rank_min']}" for d in g1["drugs"]) + ".")
        a("")
        a("Top 20 with literature up to 2018 (approved small molecules):")
        a("")
        a("| # | drug | score | gold 1 | later trialled (gold 2) | first approval | approved after 2018 |")
        a("|---|---|---|---|---|---|---|")
        for x in t["top20"]["approved_small_molecule"]:
            a(f"| {x['position']} | {x['name']} | {x['score']} | {'x' if x['in_gold1'] else ''} | {'x' if x['in_gold2_later_trialled'] else ''} | {x['first_approval'] if x['first_approval'] else 'n/a'} | {'x' if x['first_approved_after_2018'] else ''} |")
        later = t["later_approved_first_approval_after_2018"]["approved_small_molecule"]
        if later:
            a("")
            a("Drugs first approved after 2018 with a literature signal or a later trial (approved small molecules): " + ", ".join(
                f"{d['name']} (#{d['rank_min']}{', trialled' if d['in_gold2'] else ''})" for d in sorted(later, key=lambda d: (d['rank_min'] or 10**9, d['name']))[:20]) + ".")
    else:
        a("Not computed (disease has no corrected gold 2).")
    a("")
    a("## 4. Mechanism classes of the disease genes")
    a("")
    for c, b in r["mechanism_classes"].items():
        a(f"### {c}: {', '.join(b['genes'])}")
        a("")
        a(f"Retrieved papers naming a gene of this class: {b['n_retrieved_papers_mentioning_class_gene']}.")
        if b["top_drugs_by_comention_small_molecule_fda"]:
            a("")
            a("| drug | papers (class genes) | benefit | human | engine rank (small molecules) | gold 1 | gold 2 |")
            a("|---|---|---|---|---|---|---|")
            for d in b["top_drugs_by_comention_small_molecule_fda"]:
                a(f"| {d['name']} | {d['n_papers']} | {d['n_benefit']} | {d['n_human']} | {d['engine_rank_small_molecule_min']} | {'x' if d['in_gold1'] else ''} | {'x' if d['in_gold2'] else ''} |")
        else:
            a("")
            a("No approved small molecule is co-mentioned with a gene of this class in the retrieved abstracts.")
        a("")
    mc = r["mitotic_checkpoint_class"]
    a(f"### Mitotic-checkpoint class (explicit; the class relevant to BUB1B in the MVA case)")
    a("")
    a(f"Genes of {r['disease']} in `{MITOTIC}`: {', '.join(mc['genes']) or 'none'}. Universe drugs whose ChEMBL mechanism is a mitotic-checkpoint, Aurora, PLK or TTK mechanism ({mc['note']}):")
    a("")
    if mc["universe_drugs_with_mitotic_checkpoint_mechanism"]:
        a("| drug | type | small molecule + FDA | Wilms engine rank (all) | papers with a Wilms gene |")
        a("|---|---|---|---|---|")
        for d in mc["universe_drugs_with_mitotic_checkpoint_mechanism"]:
            a(f"| {d['name']} | {d['molecule_type']} | {d['small_molecule_fda']} | {d['engine_rank_all_min']} | {d['n_papers_with_wilms_gene']} |")
    else:
        a(f"None: the approved-drug universe contains {mc['n_universe_drugs_with_mitotic_checkpoint_mechanism']} drugs with such a mechanism (no approved mitotic-checkpoint, Aurora, PLK or TTK inhibitor). "
          "The approved drugs that act on the mitotic spindle are the tubulin binders:")
    a("")
    a("Tubulin binders (caution_list class a; they act on the mitotic spindle):")
    a("")
    a("| drug | type | small molecule + FDA | Wilms engine rank (all) | papers with a Wilms gene | gold 1 | gold 2 |")
    a("|---|---|---|---|---|---|---|")
    for d in mc["tubulin_binders_acting_on_the_mitotic_spindle"]:
        sup = r["label_support"].get(d["name"])
        a(f"| {d['name']} | {d['molecule_type']} | {d['small_molecule_fda']} | {d['engine_rank_all_min']} | {d['n_papers_with_wilms_gene']} | {'x' if d['in_gold1'] else ''} | {'x' if d['in_gold2'] else ''} |"
          + (f" (label bands 0-4/5-11/12-17: {sup['bands']['0_4']} / {sup['bands']['5_11']} / {sup['bands']['12_17']})" if sup else ""))
    a("")
    a("### Mechanism x age-band table (top co-mentioned approved small molecules per class; label support)")
    a("")
    a("| class | 0-4 y | 5-11 y | 12-17 y |")
    a("|---|---|---|---|")
    for c, bands in r["mechanism_x_age_band"].items():
        cells = []
        for b in BANDS:
            cells.append("; ".join(f"{st}: {', '.join(ns)}" for st, ns in bands[b].items()) or "-")
        a(f"| {c} | " + " | ".join(cells) + " |")
    a("")
    a("## 5. Evidence gaps (no dosing, no treatment recommendation)")
    a("")
    a("Evidence rungs: " + "; ".join(f"{x['rung']} {x['meaning']}" for x in r["method"]["evidence_rungs"]) + ".")
    a("")
    for g in r["evidence_gaps"]:
        a(f"### #{g['position_in_filtered_ranking']} {g['name']} - level reached {g['evidence_level_reached']}")
        a("")
        a("- Missing: " + ("; ".join(f"{m['rung']}" for m in g["missing"]) or "none") + (f". Age bands without labelled support: {', '.join(BAND_TEXT[b] for b in g['age_band_gaps'])}" if g["age_band_gaps"] else "") + ".")
        a(f"- Single study or data that would close the gap: {g['single_study_or_data_that_would_close_the_gap']}.")
        for o in g["other_options"]:
            a(f"- Other option: {o}.")
        hz = g["hazards_label_documented"]
        if isinstance(hz, str):
            a(f"- Hazards: {hz}.")
        else:
            for h in hz[:3]:
                a(f"- Label sentence ({h['section']}{'' if h['section'] == h['category'] else ', ' + h['category']}): \"{_q(h['quote'], 300)}\"")
        if g["oncology_drug_chemotherapy_caution"]:
            a("- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.")
    a("")
    a("## Limits")
    a("")
    for x in r["limits"]:
        a(f"- {x}")
    a("")
    return "\n".join(L)


def main() -> None:
    res = run()
    print(json.dumps({"top5": [x["name"] for x in res["ranking_top20_universe"][:5]], "repro": res["reproduction_checks"]}, indent=1, sort_keys=True))


if __name__ == "__main__":
    sys.exit(main())
