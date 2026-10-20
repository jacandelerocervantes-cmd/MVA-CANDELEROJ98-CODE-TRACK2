r"""EXPLORATORY, DEV-ONLY, POST HOC analysis of LLM-read evidence beyond the PrimeKG prior
(protocol/addendum_engine_v4_llm_exploratory.md, written before this module ran).

Run (PYTHONPATH=src):   python -m track2.engine_v4_llm_exploratory

It cannot confirm anything and does not change the pre-registered outcome (relevance kappa 0.565 < 0.6: that experiment is closed).
It never reads HELD-OUT labels, HELD-OUT gold or HELD-OUT annotation files, never runs engine_v4_llm heldout and creates no lock file.
Reuses engine_v4_llm / engine_v4 / engine_v4_dilution / engine_v4_ablation functions unchanged; edits nothing.
"""
from __future__ import annotations

import json
import math
import pickle
import time
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from track2 import engine_benchmark as v1
from track2 import engine_score as es
from track2 import engine_v4 as v4
from track2 import engine_v4_ablation as ab
from track2 import engine_v4_dilution as dl
from track2 import engine_v4_llm as lm

ROOT = v4.ROOT
SEED = 20261015
LABEL = "EXPLORATORY, DEV-ONLY, POST HOC (engine v4, LLM-read evidence beyond the prior); cannot confirm anything; does not change the pre-registered outcome"
PROTOCOL = ROOT / "protocol" / "addendum_engine_v4_llm_exploratory.md"
OUT_FILE = v4.OUT / "exploratory_v4_llm_dev.json"
FEATURES = ("F_dir", "F_rel_strict", "F_rel_any")
FEATURE_COLS = {"F_dir": 8, "F_rel_strict": 9, "F_rel_any": 10}       # appended after the 8 B0 columns
MODEL_NAME = {f: "PO_" + f for f in FEATURES}
MODEL_COLS = {"PO": list(lm.PRIOR_COLS), **{MODEL_NAME[f]: list(lm.PRIOR_COLS) + [FEATURE_COLS[f]] for f in FEATURES}}
PAIRS = [(MODEL_NAME[f], "PO") for f in FEATURES]
REPORT_FRACTIONS = (1.0, 0.1, 0.03, 0.01)
MAX_PAPERS = lm.MAX_PAPERS_PER_PAIR
MAX_REPS_NOTE = "f = 1.0 one deterministic sample, other fractions 5 replicates"

UNVALIDATED: Dict[str, Dict[str, Any]] = {
    "feature_definitions": {"value": list(FEATURES), "status": "UNVALIDATED", "conservative_default": list(FEATURES), "note": "exploratory protocol section 1"},
    "either_annotator_rule": {"value": "relevance yes by at least one annotator (F_rel_any)", "status": "UNVALIDATED",
                              "conservative_default": "relevance yes by both annotators (F_rel_strict)", "note": "relevance labels failed the 0.6 kappa gate"},
    "direction_rule": {"value": "2-of-3 vote of opus, sonnet, rule; unresolved excluded", "status": "UNVALIDATED",
                       "conservative_default": "2-of-3 vote of opus, sonnet, rule; unresolved excluded", "note": "as the parent"},
    "one_feature_per_model": {"value": True, "status": "UNVALIDATED", "conservative_default": True, "note": "avoids collinearity between the three counts"},
    "report_fractions": {"value": list(REPORT_FRACTIONS), "status": "UNVALIDATED", "conservative_default": list(REPORT_FRACTIONS), "note": "exploratory protocol section 3"},
    "inherited_constants": {"value": "top 30, 3 papers, 2 sentences, C grid, inner 3-fold, evaluation on P_d", "status": "UNVALIDATED",
                            "conservative_default": "as engine_v4_llm.UNVALIDATED", "note": "parent registry"},
}


# ============================================================================= labels and features (no gold argument)
def exploratory_label_table(items: Sequence[dict], ann_a: Dict[str, Dict[str, Any]], ann_b: Dict[str, Dict[str, Any]]) -> Dict[Tuple[str, str, str, str], Dict[str, bool]]:
    """(disease, drug, pmid, sentence) -> which of the three features count the sentence. Direction = parent consensus
    (2 of 3 of opus, sonnet, rule; unresolved never counted)."""
    out = {}
    for it in items:
        a, b = ann_a[it["id"]], ann_b[it["id"]]
        c = lm.consensus(a["relevance"], b["relevance"], a["direction"], b["direction"], it["sentence"], it["drug"])
        ben = c["direction"] == "benefit"
        out[(it["disease"], it["drug"], it["pmid"], it["sentence"])] = {
            "F_dir": ben,
            "F_rel_strict": ben and a["relevance"] == "yes" and b["relevance"] == "yes",
            "F_rel_any": ben and (a["relevance"] == "yes" or b["relevance"] == "yes")}
    return out


def exploratory_counts(disease: str, pair_index: Dict[str, List[dict]], kept: Set[int], table: Dict[Tuple[str, str, str, str], Dict[str, bool]],
                       drug_row: Dict[str, int], n_rows: int) -> np.ndarray:
    """(n_rows, 3) RAW counts [F_dir, F_rel_strict, F_rel_any] over the first MAX_PAPERS kept papers of every pair (2 sentences per paper).
    Unannotated pairs are 0. No gold argument."""
    C = np.zeros((n_rows, len(FEATURES)))
    for drug, entries in pair_index.items():
        r = drug_row.get(drug)
        if r is None:
            continue
        for e in [e for e in entries if e["i"] in kept][:MAX_PAPERS]:
            for s in e["sents"]:
                lab = table.get((disease, drug, lm.pmid_of(e["pid"]), s))
                if lab is None:
                    continue
                for j, f in enumerate(FEATURES):
                    C[r, j] += bool(lab[f])
    return C


def to_model_block(counts: np.ndarray) -> np.ndarray:
    return np.log1p(counts)


def coverage_summary(reps: Sequence[np.ndarray], y: np.ndarray) -> Dict[str, Any]:
    """reps: list over replicates of (n_pairs, 3) raw counts for all DEV candidate pairs stacked; y: (n_pairs,) gold-1 flag.
    Per feature: number and share of pairs with non-zero value, mean over replicates (min, max), overall / positives / negatives."""
    out: Dict[str, Any] = {}
    y = np.asarray(y)
    for j, f in enumerate(FEATURES):
        nz = np.array([(r[:, j] > 0) for r in reps])                      # (n_reps, n_pairs)
        n_all, n_pos, n_neg = nz.sum(1), nz[:, y == 1].sum(1), nz[:, y == 0].sum(1)
        out[f] = {"n_pairs": int(len(y)), "n_pairs_positive": int((y == 1).sum()), "n_pairs_negative": int((y == 0).sum()),
                  "nonzero_pairs_mean": float(n_all.mean()), "nonzero_pairs_min": int(n_all.min()), "nonzero_pairs_max": int(n_all.max()),
                  "nonzero_share": float(n_all.mean() / len(y)),
                  "nonzero_positive_mean": float(n_pos.mean()), "nonzero_positive_share": float(n_pos.mean() / max(1, (y == 1).sum())),
                  "nonzero_negative_mean": float(n_neg.mean()), "nonzero_negative_share": float(n_neg.mean() / max(1, (y == 0).sum())),
                  "mean_count_given_nonzero": float(np.mean([r[r[:, j] > 0, j].mean() if (r[:, j] > 0).any() else 0.0 for r in reps]))}
    return out


# ============================================================================= DEV-only context
def dev_context(bench: dict, dev: Sequence[str]) -> Dict[str, Any]:
    """Like engine_v4_llm._context but keeps only DEV diseases from every multi-split file at once, and reads gold only for DEV diseases."""
    dev_set = set(dev)
    drugs = v1.load_universe(ROOT / "data" / "universe.json")
    universe = v1.universe_index(drugs)
    names = sorted(universe)
    index = {n: i for i, n in enumerate(names)}
    classes = np.array(es.drug_classes(universe, names))
    cands_all = json.loads(lm.CAND_FILE.read_text(encoding="utf-8"))
    cands = {d: v for d, v in cands_all.items() if d in dev_set and v["split"] == "dev"}
    del cands_all
    pair_all = json.loads(lm.PAIR_INDEX_FILE.read_text(encoding="utf-8"))
    pair_idx = {d: v for d, v in pair_all.items() if d in dev_set}
    del pair_all
    with v4.FEATURES_FILE.open("rb") as fh:
        f4 = pickle.load(fh)
    with dl.FEATS_FILE.open("rb") as fh:
        fd = pickle.load(fh)
    by_name = {d["name"]: d for d in bench["diseases"] if d["name"] in dev_set}
    rank = {n: i for i, n in enumerate(sorted(d["name"] for d in bench["diseases"]))}
    ctx: Dict[str, Any] = {"names": names, "f4": {d: f4[d] for d in dev}, "fd": {d: fd[d] for d in dev}, "pair_idx": pair_idx, "rank": rank,
                           "cands": {}, "cls": {}, "pop": {}, "y": {}}
    for d in dev:
        cand = np.array(sorted(index[n] for n in cands[d]["candidates"]), dtype=int)
        ctx["cands"][d] = cand
        ctx["cls"][d] = classes[cand]
        ctx["pop"][d] = ctx["f4"][d]["popularity"][cand]
        yy = np.zeros(len(names))
        for g in by_name[d]["gold1"]:
            if g in index:
                yy[index[g]] = 1
        ctx["y"][d] = yy[cand]
    return ctx


def sample_matrices(ctx: Dict[str, Any], dev: Sequence[str], table) -> Tuple[Dict, Dict, Dict]:
    """Per disease and sample: (n_cand, 11) matrix [8 thinned B0 columns | log1p of the 3 features], n_sub, and the raw counts."""
    Zs, ns, Cs = {}, {}, {}
    for d in dev:
        cand = ctx["cands"][d]
        row_of = {ctx["names"][i]: r for r, i in enumerate(cand)}
        n_ret = ctx["pair_idx"][d]["n_retrieved"]
        Zs[d], ns[d], Cs[d] = {}, {}, {}
        for fi, rep in dl.sample_keys():
            kept = set(dl.subsample_indices(n_ret, dl.FRACTIONS[fi], ctx["rank"][d], fi, rep).tolist())
            C = exploratory_counts(d, ctx["pair_idx"][d]["pairs"], kept, table, row_of, len(cand))
            X = ctx["fd"][d]["samples"][(fi, rep)]["X"][cand]
            Zs[d][(fi, rep)], Cs[d][(fi, rep)] = np.hstack([X, to_model_block(C)]), C
            ns[d][(fi, rep)] = ctx["fd"][d]["samples"][(fi, rep)]["n_sub"]
    return Zs, ns, Cs


def make_scorers(heads: Dict[str, Dict[str, Dict[str, Any]]]) -> Dict[str, Any]:
    return {m: (lambda d, Zc, m=m: v4.head_decision(heads[m][d], Zc[:, MODEL_COLS[m]])) for m in MODEL_COLS}


def reduce_table(tab: Dict[str, Any]) -> Dict[str, Any]:
    """Keep only the reported fractions and the paired differences of each feature model over PO."""
    out = {"n_diseases": tab["n_diseases"], "n_groups": tab["n_groups"], "fractions": {}}
    for f in REPORT_FRACTIONS:
        row = tab["fractions"][str(f)]
        out["fractions"][str(f)] = {"n_papers_mean": row["n_papers_mean"], "PO_auc": row["PO"]["auc"],
                                    "paired_minus_PO": {m: row["paired"][m + "_minus_PO"] for m in (MODEL_NAME[x] for x in FEATURES)},
                                    "model_auc": {m: row[m]["auc"] for m in (MODEL_NAME[x] for x in FEATURES)}}
    return out


def run(root=ROOT) -> dict:
    t0 = time.time()
    bench = v4.load_benchmark(root)
    dev = list(bench["split"]["dev"])
    items = lm._load_split_items("dev")
    ann = lm._annotations("dev", items)
    agr = json.loads((lm.LABELS_DIR / "agreement_dev.json").read_text(encoding="utf-8"))
    table = exploratory_label_table(items, ann[lm.ANNOTATORS[0]], ann[lm.ANNOTATORS[1]])
    ctx = dev_context(bench, dev)
    Zs, ns, Cs = sample_matrices(ctx, dev, table)
    y = ctx["y"]
    members = v4.members_of(bench, dev)
    Z1 = {d: Zs[d][(0, 0)] for d in dev}
    heads, chosen = {}, {}
    for m, cols in MODEL_COLS.items():
        heads[m], chosen[m] = lm.nested_oof_heads(dev, members, Z1, y, ctx["cls"], cols)
    res = lm.score_models_pd(dev, Zs, ns, make_scorers(heads), y, ctx["pop"])
    tab = ab.table(res, members, y, tuple(MODEL_COLS), PAIRS)
    tab_dl = ab.table(res, members, y, tuple(MODEL_COLS), PAIRS, disease_level=True)
    y_all = np.concatenate([y[d] for d in dev])
    cov = {}
    for fi, f in enumerate(dl.FRACTIONS):
        if f not in REPORT_FRACTIONS:
            continue
        reps = range(1) if fi == 0 else range(dl.N_REPS)
        cov[str(f)] = coverage_summary([np.vstack([Cs[d][(fi, r)] for d in dev]) for r in reps], y_all)
    coef = {}
    for m, cols in MODEL_COLS.items():
        coef[m] = {"columns": cols, "chosen_C_per_outer_fold": chosen[m]}
    out = {"label": LABEL, "status": "exploratory, DEV only; cannot confirm anything; the pre-registered experiment stays closed (relevance kappa below 0.6)",
           "protocol": str(PROTOCOL.relative_to(root)), "protocol_sha256": v4.sha256_file(PROTOCOL),
           "heldout_touched": False, "lock_file_created": False,
           "unvalidated_registry": UNVALIDATED,
           "agreement_dev_from_file": {"direction_kappa": agr["direction"]["kappa"], "relevance_kappa": agr["relevance"]["kappa"], "gate": lm.KAPPA_GATE},
           "n_items_dev": len(items), "n_dev_diseases": len(dev), "n_dev_candidate_pairs": int(len(y_all)), "n_dev_positive_pairs": int(y_all.sum()),
           "replicates": MAX_REPS_NOTE, "seed": SEED, "bootstrap": {"n": es.N_BOOT, "unit": "engine_v4 group", "seed": dl.SEED},
           "models": coef,
           "feature_totals_f1": {f: int(sum(Cs[d][(0, 0)][:, j].sum() for d in dev)) for j, f in enumerate(FEATURES)},
           "paired_auc_group_bootstrap": reduce_table(tab), "paired_auc_disease_bootstrap_sensitivity": reduce_table(tab_dl),
           "coverage_dev_candidate_pairs": cov,
           "reading": "paired macro AUC of (PO + one feature) minus PO on DEV out-of-fold heads, candidate sets P_d; descriptive intervals, no multiplicity correction",
           "seconds": round(time.time() - t0, 1)}
    v4.write_json(OUT_FILE, out)
    print(json.dumps({"done": True, "out": str(OUT_FILE)}), flush=True)
    return out


def main() -> None:
    run()


if __name__ == "__main__":
    main()
