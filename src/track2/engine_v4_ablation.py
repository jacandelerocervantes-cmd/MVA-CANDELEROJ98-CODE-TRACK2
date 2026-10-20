r"""Engine v4 prior / literature ablation of the dilution result (protocol/addendum_engine_v4_ablation.md, written before this module ran).

Run (PYTHONPATH=src), in this order (cache only, no network, no embeddings):
    python -m track2.engine_v4_ablation dev       frozen ablations + refit models (grouped CV) on DEV
    python -m track2.engine_v4_ablation heldout   the same on HELD-OUT, once (lock file); evaluates the decision rules

Reuses engine_v4 / engine_v4_dilution / engine_score functions unchanged; edits nothing. Label: post hoc development.
"""
from __future__ import annotations

import json
import pickle
import sys
import time
from typing import Any, Callable, Dict, List, Sequence, Tuple

import numpy as np

from track2 import engine_benchmark as v1
from track2 import engine_score as es
from track2 import engine_v4 as v4
from track2 import engine_v4_dilution as dl

ROOT = v4.ROOT
SEED = dl.SEED
LABEL = "post hoc development (engine v4 prior/literature ablation), not a pre-registered primary result"
OUT = v4.OUT
PROTOCOL = ROOT / "protocol" / "addendum_engine_v4_ablation.md"
OUT_FILE = OUT / "validation_v4_ablation.json"
LOCK = OUT / "heldout_v4_ablation_scored.lock"
DIL_JSON = OUT / "validation_v4_dilution.json"
FRACTIONS = dl.FRACTIONS
LIT_COLS = [0, 1, 2, 3, 5]          # paper-derived (thinned in the dilution experiment)
PRIOR_COLS = [4, 6, 7]              # log_popularity + PrimeKG prior (fixed in the dilution experiment)
MODELS = ("B0", "PRIOR_ONLY", "LIT_ONLY", "PRIOR_ONLY_refit", "LIT_ONLY_refit")
REFIT_COLS = {"PRIOR_ONLY_refit": PRIOR_COLS, "LIT_ONLY_refit": LIT_COLS}
PAIRS = [("B0", "PRIOR_ONLY"), ("B0", "LIT_ONLY"), ("LIT_ONLY", "PRIOR_ONLY"),
         ("B0", "PRIOR_ONLY_refit"), ("B0", "LIT_ONLY_refit"), ("LIT_ONLY_refit", "PRIOR_ONLY_refit")]

UNVALIDATED: Dict[str, Dict[str, Any]] = {
    "literature_group_columns": {"value": LIT_COLS, "status": "UNVALIDATED", "conservative_default": LIT_COLS, "note": "protocol section 1"},
    "prior_group_columns": {"value": PRIOR_COLS, "status": "UNVALIDATED", "conservative_default": PRIOR_COLS, "note": "protocol section 1"},
    "removal_operation": {"value": "replace by the DEV training mean", "status": "UNVALIDATED",
                          "conservative_default": "replace by the DEV training mean", "note": "equals zero after standardisation"},
    "h_lit_rule": {"value": "LIT_ONLY frozen, HELD-OUT AUC gap CI lower bound > 0 at both f=0.03 and f=0.01", "status": "UNVALIDATED",
                   "conservative_default": "same", "note": "protocol section 4"},
    "h_prior_rule": {"value": "PRIOR_ONLY frozen gap at f=0.01 inside the 95% CI of the B0 gap, HELD-OUT", "status": "UNVALIDATED",
                     "conservative_default": "same", "note": "protocol section 4"},
    "refit_models_sensitivity_only": {"value": "refit models do not enter the decision rules", "status": "UNVALIDATED",
                                      "conservative_default": "same", "note": "protocol section 4"},
}


# ============================================================================= scorers
def group_means(f4: Dict[str, Any], dev: Sequence[str]) -> np.ndarray:
    """DEV training mean of every B0 column over all (DEV disease, drug) rows at f = 1.0 (no class masking)."""
    return np.vstack([f4[d]["X"] for d in dev]).mean(0)


def mean_replaced(X: np.ndarray, cols: Sequence[int], mean: np.ndarray) -> np.ndarray:
    out = X.copy()
    out[:, list(cols)] = mean[list(cols)]
    return out


def frozen_scorer(model: str, head: Dict[str, Any], mean: np.ndarray) -> Callable[[str, np.ndarray], np.ndarray]:
    """Frozen B0 weights; PRIOR_ONLY sets the literature columns, LIT_ONLY the prior columns, to the DEV mean."""
    cols = {"B0": [], "PRIOR_ONLY": LIT_COLS, "LIT_ONLY": PRIOR_COLS}[model]

    def f(d: str, X: np.ndarray) -> np.ndarray:
        return v4.head_decision(head, mean_replaced(X, cols, mean)[:, head["cols"]])
    return f


def refit_scorer(cols: Sequence[int], heads: Dict[str, Dict[str, Any]]) -> Callable[[str, np.ndarray], np.ndarray]:
    """heads: disease -> fitted head (out-of-fold head of its fold on DEV; the all-DEV head for every HELD-OUT disease)."""
    def f(d: str, X: np.ndarray) -> np.ndarray:
        return v4.head_decision(heads[d], X[:, list(cols)])
    return f


def fit_refit_heads(dev: Sequence[str], members: Dict[str, List[str]], Z: Dict[str, np.ndarray], y: Dict[str, np.ndarray],
                    classes: np.ndarray, cols: Sequence[int]) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
    """Grouped 5-fold: per-disease out-of-fold head (training rows ATC-masked as in v1) and the all-DEV head (no masking)."""
    folds = v4.group_folds(members)
    oof: Dict[str, Dict[str, Any]] = {}
    for fold in folds:
        train = [d for d in dev if d not in fold]
        keep = es._class_mask(classes, es._gold_classes(fold, y, classes))
        h = v4.fit_rows("logit", cols, train, Z, y, keep)
        for d in fold:
            oof[d] = h
    return oof, v4.fit_rows("logit", cols, dev, Z, y, None)


def oof_f1_check(dev, members, Z, y, classes, cols, oof_heads) -> float:
    """Maximum |difference| between our f=1.0 out-of-fold scores and engine_v4.oof_head(nested_platt=True)."""
    ref = v4.oof_head("logit", list(cols), dev, members, Z, y, classes, True)
    return float(max(np.abs(ref[d]["score"] - v4.head_decision(oof_heads[d], Z[d][:, list(cols)])).max() for d in dev))


# ============================================================================= scoring
def score_models(diseases: Sequence[str], fd: Dict[str, Any], scorers: Dict[str, Callable[[str, np.ndarray], np.ndarray]],
                 y: Dict[str, np.ndarray], pop: Dict[str, np.ndarray]) -> Dict[str, Any]:
    """Same result layout as engine_v4_dilution.score_diseases: res[d] = {n_sub, metrics[model][fi] -> [(auc,p10,r20)], popularity}."""
    res: Dict[str, Any] = {}
    for d in diseases:
        r: Dict[str, Any] = {"n_sub": {}, "metrics": {m: {} for m in scorers}}
        for fi, rep in dl.sample_keys():
            s = fd[d]["samples"][(fi, rep)]
            r["n_sub"].setdefault(fi, []).append(s["n_sub"])
            for m, fn in scorers.items():
                r["metrics"][m].setdefault(fi, []).append(dl.metrics_one(fn(d, s["X"]), y[d]))
        r["popularity"] = dl.metrics_one(pop[d], y[d])
        res[d] = r
    return res


def table(res: Dict[str, Any], members: Dict[str, List[str]], y: Dict[str, np.ndarray], models: Sequence[str],
          pairs: Sequence[Tuple[str, str]], disease_level: bool = False) -> Dict[str, Any]:
    """Per fraction and model: macro AUC, gap vs popularity (CI), P@10, R@20; paired differences with CI (same bootstrap indices)."""
    if disease_level:
        members = {d: [d] for g in members for d in members[g]}
    gs = sorted(members)
    dn = sorted(d for g in gs for d in members[g])
    idx = dl.boot_idx(len(gs))
    pop = {k: np.array([res[d]["popularity"][i] for d in dn]) for i, k in enumerate(("auc", "p10", "r20"))}
    pop_g = {k: v4.agg_groups(v, dn, members) for k, v in pop.items()}
    out: Dict[str, Any] = {"n_diseases": len(dn), "n_groups": len(gs), "popularity": {k: es.macro_with_ci(v, idx) for k, v in pop_g.items()},
                           "fractions": {}}
    for fi, f in enumerate(FRACTIONS):
        row: Dict[str, Any] = {"fraction": f}
        ns = [x for d in dn for x in res[d]["n_sub"][fi]]
        row["n_papers_mean"] = float(np.mean(ns))
        garr: Dict[str, Dict[str, np.ndarray]] = {}
        for m in models:
            arr = dl.per_fraction_arrays(res, dn, m, fi)
            g = {k: v4.agg_groups(v, dn, members) for k, v in arr.items()}
            garr[m] = g
            row[m] = {"auc": es.macro_with_ci(g["auc"], idx), "auc_gap_vs_popularity": es.diff_with_ci(g["auc"], pop_g["auc"], idx),
                      "precision_at_10": es.macro_with_ci(g["p10"], idx), "recall_at_20": es.macro_with_ci(g["r20"], idx),
                      "precision_at_10_minus_popularity": es.diff_with_ci(g["p10"], pop_g["p10"], idx),
                      "n_replicates": len(res[dn[0]]["metrics"][m][fi])}
        row["paired"] = {f"{a}_minus_{b}": {"auc": es.diff_with_ci(garr[a]["auc"], garr[b]["auc"], idx),
                                            "precision_at_10": es.diff_with_ci(garr[a]["p10"], garr[b]["p10"], idx)}
                         for a, b in pairs if a in garr and b in garr}
        out["fractions"][str(f)] = row
    return out


# ============================================================================= decision rules (protocol section 4)
def decide(tab: Dict[str, Any]) -> Dict[str, Any]:
    """H-LIT / H-PRIOR on a HELD-OUT table (group bootstrap). Pure function of the table."""
    fr = tab["fractions"]

    def gap(m: str, f: float) -> Dict[str, Any]:
        return fr[str(f)][m]["auc_gap_vs_popularity"]
    lit = {str(f): {"gap": gap("LIT_ONLY", f)["estimate"], "ci95": gap("LIT_ONLY", f)["ci95"],
                    "lower_above_0": bool(gap("LIT_ONLY", f)["ci95"][0] is not None and gap("LIT_ONLY", f)["ci95"][0] > 0)} for f in (0.03, 0.01)}
    h_lit = all(v["lower_above_0"] for v in lit.values())
    status_lit = "holds" if h_lit else ("partial" if any(v["lower_above_0"] for v in lit.values()) else "does not hold")
    p, b = gap("PRIOR_ONLY", 0.01), gap("B0", 0.01)
    inside = bool(b["ci95"][0] <= p["estimate"] <= b["ci95"][1])
    pd = fr["0.01"]["paired"]["B0_minus_PRIOR_ONLY"]["auc"]
    h_prior = inside
    combos = {(True, True): "both: the prior reproduces B0's low-volume gap, yet the literature features alone also keep a gap above 0",
              (True, False): "H-PRIOR only: the low-volume gap is explained by the prior; no evidence that the literature contributes",
              (False, True): "H-LIT only: the literature features contribute; the prior alone does not reproduce B0's gap",
              (False, False): "neither: inconclusive"}
    return {"H_LIT_literature_contributes_at_low_volume": {"holds": h_lit, "status": status_lit, "literature_only_frozen": lit},
            "H_PRIOR_gap_explained_by_prior": {"holds": h_prior, "prior_only_gap_f0.01": p["estimate"], "prior_only_ci95": p["ci95"],
                                               "b0_gap_f0.01": b["estimate"], "b0_ci95": b["ci95"], "prior_only_inside_b0_ci": inside,
                                               "corroboration_paired_B0_minus_PRIOR_ONLY_f0.01": pd,
                                               "paired_ci_contains_0": bool(pd["ci95"][0] <= 0 <= pd["ci95"][1])},
            "outcome": combos[(h_prior, h_lit)]}


# ============================================================================= context and blocks
def _context(root, which_labels: Sequence[str]) -> Dict[str, Any]:
    bench = v4.load_benchmark(root)
    drugs = v1.load_universe(root / "data" / "universe.json")
    universe = v1.universe_index(drugs)
    names = sorted(universe)
    index = {n: i for i, n in enumerate(names)}
    with v4.FEATURES_FILE.open("rb") as fh:
        f4 = pickle.load(fh)
    with dl.FEATS_FILE.open("rb") as fh:
        fd = pickle.load(fh)
    by_name = {d["name"]: d for d in bench["diseases"]}
    y: Dict[str, np.ndarray] = {}
    for d in which_labels:
        yy = np.zeros(len(names))
        for g in by_name[d]["gold1"]:
            if g in index:
                yy[index[g]] = 1
        y[d] = yy
    return {"bench": bench, "names": names, "classes": np.array(es.drug_classes(universe, names)), "f4": f4, "fd": fd, "y": y}


def _reproduction_check(tab: Dict[str, Any], which: str) -> Dict[str, Any]:
    """B0 macro AUC and gap per fraction against validation_v4_dilution.json (same frozen model, same features)."""
    if not DIL_JSON.exists():
        return {"available": False}
    ref = json.loads(DIL_JSON.read_text(encoding="utf-8"))[which]["table_group_bootstrap"]["fractions"]
    worst = 0.0
    for k, row in tab["fractions"].items():
        for key in ("auc", "auc_gap_vs_popularity"):
            worst = max(worst, abs(row["B0"][key]["estimate"] - ref[k]["B0"][key]["estimate"]))
    return {"available": True, "max_abs_diff_B0_auc_and_gap_vs_dilution_json": worst}


def _heads_and_scorers(ctx, dev, which: str):
    fz = dl._verify(ROOT, ctx["bench"])
    m = fz["model"]
    head = m["heads"]["B0"]
    mean = group_means(ctx["f4"], dev)
    Zdev = {d: ctx["f4"][d]["X"] for d in dev}
    members = v4.members_of(ctx["bench"], dev)
    scorers: Dict[str, Any] = {"B0": frozen_scorer("B0", head, mean), "PRIOR_ONLY": frozen_scorer("PRIOR_ONLY", head, mean),
                               "LIT_ONLY": frozen_scorer("LIT_ONLY", head, mean)}
    checks: Dict[str, Any] = {"dev_mean_vs_frozen_mu_max_abs_diff": float(np.abs(mean - np.asarray(head["mu"])).max()),
                              "frozen_check": fz["check"]}
    return fz, head, mean, Zdev, members, scorers, checks


def run_dev(root=ROOT) -> dict:
    t0 = time.time()
    bench = v4.load_benchmark(root)
    dev = bench["split"]["dev"]
    ctx = _context(root, dev)
    fz, head, mean, Zdev, members, scorers, checks = _heads_and_scorers(ctx, dev, "dev")
    refit_info: Dict[str, Any] = {}
    for name, cols in REFIT_COLS.items():
        oof, full = fit_refit_heads(dev, members, Zdev, ctx["y"], ctx["classes"], cols)
        scorers[name] = refit_scorer(cols, oof)
        refit_info[name] = {"columns": cols, "oof_f1_scores_vs_engine_v4_oof_head_max_abs_diff": oof_f1_check(dev, members, Zdev, ctx["y"], ctx["classes"], cols, oof),
                            "all_dev_coef": [float(c) for c in full["coef"]], "all_dev_intercept": full["intercept"]}
    pop = {d: ctx["f4"][d]["popularity"] for d in dev}
    res = score_models(dev, ctx["fd"], scorers, ctx["y"], pop)
    tab = table(res, members, ctx["y"], MODELS, PAIRS)
    prior_aucs = [tab["fractions"][str(f)]["PRIOR_ONLY"]["auc"]["estimate"] for f in FRACTIONS]
    block = {"n_diseases": len(dev), "n_groups": len(members), "table_group_bootstrap": tab,
             "table_disease_level_bootstrap_sensitivity": table(res, members, ctx["y"], MODELS, PAIRS, disease_level=True),
             "gap_curve": {m: dl.fit_gap_curve(res, members, m) for m in ("B0", "PRIOR_ONLY", "LIT_ONLY", "PRIOR_ONLY_refit", "LIT_ONLY_refit")},
             "checks": dict(checks, prior_only_auc_range_over_fractions=float(max(prior_aucs) - min(prior_aucs)),
                            b0_reproduction=_reproduction_check(tab, "dev"), refit=refit_info),
             "scoring_seconds": round(time.time() - t0, 1),
             "reading": "DEV: frozen B0 and frozen ablations are in-sample for the weights; refit models are out-of-fold"}
    out = {"label": LABEL, "protocol": str(PROTOCOL.relative_to(ROOT)), "protocol_sha256_at_dev_run": v4.sha256_file(PROTOCOL),
           "settings": {"fractions": FRACTIONS, "replicates": dl.N_REPS, "floor_papers": dl.FLOOR, "seed": SEED, "n_boot": es.N_BOOT,
                        "models": list(MODELS), "literature_columns": LIT_COLS, "prior_columns": PRIOR_COLS,
                        "feature_names": es.FEATURE_NAMES},
           "unvalidated_registry": UNVALIDATED, "status": "DEV part written; HELD-OUT not scored yet", "dev": block}
    v4.write_json(OUT_FILE, out)
    print(json.dumps({"dev_done": True, "gap_f0.01": {m: tab["fractions"]["0.01"][m]["auc_gap_vs_popularity"]["estimate"] for m in MODELS}},
                     sort_keys=True), flush=True)
    return out


def run_heldout(root=ROOT) -> dict:
    if LOCK.exists():
        raise RuntimeError("HELD-OUT ablation already scored once (lock file exists); refusing a second scoring")
    out = json.loads(OUT_FILE.read_text(encoding="utf-8"))
    if "dev" not in out:
        raise RuntimeError("run the dev command first")
    bench = v4.load_benchmark(root)
    dev, ho = bench["split"]["dev"], bench["split"]["heldout"]
    fz = dl._verify(root, bench)
    LOCK.write_text(json.dumps({"scored_at": "first and only held-out ablation scoring", "frozen_file_sha256": fz["file_sha256"],
                                "protocol_sha256": v4.sha256_file(PROTOCOL)}) + "\n", encoding="utf-8")
    ctx = _context(root, list(dev) + list(ho))            # HELD-OUT labels are read only here, after the lock exists
    fz, head, mean, Zdev, members_dev, scorers, checks = _heads_and_scorers(ctx, dev, "heldout")
    refit_info: Dict[str, Any] = {}
    for name, cols in REFIT_COLS.items():
        _, full = fit_refit_heads(dev, members_dev, Zdev, ctx["y"], ctx["classes"], cols)
        scorers[name] = refit_scorer(cols, {d: full for d in ho})
        refit_info[name] = {"columns": cols, "all_dev_coef": [float(c) for c in full["coef"]], "all_dev_intercept": full["intercept"]}
    members = v4.members_of(bench, ho)
    pop = {d: ctx["f4"][d]["popularity"] for d in ho}
    res = score_models(ho, ctx["fd"], scorers, ctx["y"], pop)
    tab = table(res, members, ctx["y"], MODELS, PAIRS)
    prior_aucs = [tab["fractions"][str(f)]["PRIOR_ONLY"]["auc"]["estimate"] for f in FRACTIONS]
    out["heldout"] = {"n_diseases": len(ho), "n_groups": len(members), "table_group_bootstrap": tab,
                      "table_disease_level_bootstrap_sensitivity": table(res, members, ctx["y"], MODELS, PAIRS, disease_level=True),
                      "gap_curve": {m: dl.fit_gap_curve(res, members, m) for m in ("B0", "PRIOR_ONLY", "LIT_ONLY", "PRIOR_ONLY_refit", "LIT_ONLY_refit")},
                      "checks": dict(checks, prior_only_auc_range_over_fractions=float(max(prior_aucs) - min(prior_aucs)),
                                     b0_reproduction=_reproduction_check(tab, "heldout"), refit=refit_info)}
    out["decision_rules_heldout"] = decide(tab)
    sens = {}
    for m in ("PRIOR_ONLY_refit", "LIT_ONLY_refit"):
        sens[m] = {str(f): tab["fractions"][str(f)][m]["auc_gap_vs_popularity"] for f in FRACTIONS}
    out["refit_sensitivity_heldout_gap"] = sens
    out["protocol_sha256_at_heldout_run"] = v4.sha256_file(PROTOCOL)
    out["status"] = "complete: DEV and HELD-OUT; HELD-OUT scored once"
    v4.write_json(OUT_FILE, out)
    print(json.dumps({"heldout_done": True, "decision": out["decision_rules_heldout"]["outcome"]}, sort_keys=True), flush=True)
    return out


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = {"dev": run_dev, "heldout": run_heldout}.get(cmd)
    if fn is None:
        raise SystemExit("usage: python -m track2.engine_v4_ablation dev|heldout")
    fn()


if __name__ == "__main__":
    main()
