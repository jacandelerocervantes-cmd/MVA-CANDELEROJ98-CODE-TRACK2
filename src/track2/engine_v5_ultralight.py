r"""Engine v5 ULTRALIGHT (protocol/addendum_engine_v5_ultralight.md, written before this module ran).

Run (PYTHONPATH=src, python without -I):
    python -m track2.engine_v5_ultralight select   items of HELD-OUT records that survive the thinning at f = 0.03 or 0.01 in >= 1 of 5 replicates;
                                                   Sonnet shards, seeded 30 % Opus shards, manifest (no LLM is called, no score, no gold)
    python -m track2.engine_v5_ultralight score    needs the labels; kappa gate, fixed-form feature B, single scoring with a lock file

Fixed a-priori feature, no training: score = PO_score + delta * B. Reuses engine_v5 / engine_v4* functions unchanged; edits nothing.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence, Set, Tuple

import numpy as np

from track2 import engine_extractors as ex
from track2 import engine_score as es
from track2 import engine_v4 as v4
from track2 import engine_v4_ablation as ab
from track2 import engine_v4_dilution as dl
from track2 import engine_v4_llm as lm
from track2 import engine_v5 as e5

ROOT = v4.ROOT
OUT = v4.OUT
SEED = e5.SEED
LABEL = "post hoc development (engine v5 ULTRALIGHT), extends the pre-registered v4 programme; not a pre-registered primary result"
PROTOCOL = ROOT / "protocol" / "addendum_engine_v5_ultralight.md"
FULL_PROTOCOL = ROOT / "protocol" / "addendum_engine_v5.md"
INPUT_DIR = OUT / "llm_input_v5_ul"
LABELS_DIR = OUT / "llm_labels_v5_ul"
OUT_FILE = OUT / "validation_v5_ul.json"
LOCK = OUT / "heldout_v5_ul_scored.lock"
FRACTION_INDEX = {0.03: 3, 0.01: 4}
FRACTIONS_UL = (0.03, 0.01)
N_REPS = e5.N_REPS
SHARD_SIZE = e5.SHARD_SIZE
OPUS_SHARE = 0.30
OPUS_SEED = SEED + 7
DELTA_SD_FRACTION = 0.5            # UNVALIDATED; conservative default 0 (no effect)
KAPPA_GATE = lm.KAPPA_GATE         # 0.6
COVERAGE_MIN = 0.01
TOKENS_PER_ITEM = 155


# ============================================================================= selection
def surviving_pids(pids: Sequence[str], rank: int, drug_index: int, frac_index: int, rep: int, f: float) -> List[str]:
    """Retrieved records (in retrieval order) kept by the v5 Bernoulli thinning of this sample."""
    return [pids[i] for i in e5.thin_records(len(pids), f, rank, drug_index, frac_index, rep)]


def survivors_by_sample(retrieval_pairs: Dict[str, Dict[str, Any]], rank: int, drug_index: Dict[str, int]) -> Dict[Tuple[float, int], Dict[str, Set[str]]]:
    """(f, rep) -> {drug: set of surviving pids} for one disease."""
    out: Dict[Tuple[float, int], Dict[str, Set[str]]] = {}
    for f in FRACTIONS_UL:
        for rep in range(N_REPS):
            out[(f, rep)] = {g: set(surviving_pids(r["pids"], rank, drug_index[g], FRACTION_INDEX[f], rep, f)) for g, r in retrieval_pairs.items()}
    return out


def select_items(heldout_items: Sequence[dict], union_keys: Set[Tuple[str, str, str]]) -> List[dict]:
    """HELD-OUT v5 items whose (disease, drug, pmid) is in the union of surviving records; ordered by id (ids unchanged)."""
    return sorted((it for it in heldout_items if (it["disease"], it["drug"], it["pmid"]) in union_keys), key=lambda it: it["id"])


def opus_shards(n_shards: int, share: float = OPUS_SHARE, seed: int = OPUS_SEED) -> List[int]:
    """Seeded random subset of shard numbers (1-based), round(share * n) of them (at least 1)."""
    k = max(1, int(round(share * n_shards)))
    return sorted(random.Random(seed).sample(range(1, n_shards + 1), k))


def strip_item(it: dict) -> dict:
    return {k: it[k] for k in lm.ITEM_FIELDS}


def compute_selection(root: Path = ROOT) -> Dict[str, Any]:
    ctx = e5._context(root)
    doc = json.loads(e5.RETRIEVAL_FILE.read_text(encoding="utf-8"))
    ho_dis = [d for d in sorted(doc["diseases"]) if doc["diseases"][d]["split"] == "heldout"]
    union: Set[Tuple[str, str, str]] = set()
    per_pair_union: Dict[Tuple[str, str], Set[str]] = {}
    for d in ho_dis:
        sv = survivors_by_sample(doc["diseases"][d]["pairs"], ctx["rank"][d], ctx["index"])
        for g in doc["diseases"][d]["pairs"]:
            pids = set().union(*(sv[k][g] for k in sv))
            per_pair_union[(d, g)] = pids
            union.update((d, g, lm.pmid_of(p)) for p in pids)
    items: List[dict] = []
    for p in sorted(e5.INPUT_DIR.glob("heldout_shard_*.json")):
        items.extend(json.loads(p.read_text(encoding="utf-8"))["items"])
    sel = select_items(items, union)
    return {"ctx": ctx, "doc": doc, "ho_diseases": ho_dis, "n_ho_items": len(items), "n_ho_pairs": len(per_pair_union), "union": union,
            "per_pair_union": per_pair_union, "selected": sel}


def run_select(root: Path = ROOT) -> Dict[str, Any]:
    s = compute_selection(root)
    sel = s["selected"]
    e5.assert_item_fields([strip_item(it) for it in sel])
    shards = lm.shard(sel, SHARD_SIZE)
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    for old in list(INPUT_DIR.glob("ul_*.json")) + list(INPUT_DIR.glob("manifest.json")):
        old.unlink()
    files = []
    for n, sh in enumerate(shards, 1):
        fn = INPUT_DIR / f"ul_sonnet_{n:03d}.json"
        v4.write_json(fn, {"split": "heldout", "shard": n, "n_items": len(sh), "items": [strip_item(it) for it in sh]})
        files.append({"file": fn.name, "n_items": len(sh), "sha256": v4.sha256_file(fn)})
    op = opus_shards(len(shards))
    opus_files = []
    for n in op:
        src = INPUT_DIR / f"ul_sonnet_{n:03d}.json"
        fn = INPUT_DIR / f"ul_opus_{n:03d}.json"
        fn.write_bytes(src.read_bytes())
        opus_files.append({"file": fn.name, "n_items": files[n - 1]["n_items"], "sha256": v4.sha256_file(fn)})
    n_opus_items = sum(f["n_items"] for f in opus_files)
    n_pairs = s["n_ho_pairs"]
    pairs_with_items = len({(it["disease"], it["drug"]) for it in sel})
    guideline_tokens = lm_guideline_tokens()
    man = {"label": LABEL, "protocol": str(PROTOCOL.relative_to(root)), "protocol_sha256": v4.sha256_file(PROTOCOL),
           "full_protocol_sha256": v4.sha256_file(FULL_PROTOCOL), "seed": SEED, "opus_seed": OPUS_SEED, "opus_share": OPUS_SHARE,
           "fractions": list(FRACTIONS_UL), "n_replicates": N_REPS, "shard_size": SHARD_SIZE,
           "guidelines": {"direction": "protocol/annotation_guideline.md", "relevance": "protocol/relevance_guideline_v2.md"}, "guideline_sha256": dict(e5.GUIDELINE_SHA256),
           "candidate_file_sha256": v4.sha256_file(e5.CAND_FILE), "retrieval_file_sha256": v4.sha256_file(e5.RETRIEVAL_FILE),
           "item_fields": list(lm.ITEM_FIELDS), "n_heldout_items_v5": s["n_ho_items"], "n_heldout_pairs": n_pairs, "n_heldout_diseases": len(s["ho_diseases"]),
           "n_items_selected": len(sel), "share_of_heldout_items": len(sel) / s["n_ho_items"], "n_pairs_with_selected_items": pairs_with_items,
           "n_sonnet_shards": len(shards), "n_opus_shards": len(op), "n_opus_items": n_opus_items, "opus_shard_numbers": op,
           "items_per_heldout_pair": len(sel) / n_pairs, "items_per_pair_with_items": len(sel) / pairs_with_items,
           "tokens_per_item": TOKENS_PER_ITEM, "token_estimate_sonnet": TOKENS_PER_ITEM * len(sel), "token_estimate_opus": TOKENS_PER_ITEM * n_opus_items,
           "token_estimate_both": TOKENS_PER_ITEM * (len(sel) + n_opus_items), "token_estimate_both_if_opus_on_all": TOKENS_PER_ITEM * 2 * len(sel),
           "tokens_per_heldout_pair_both": TOKENS_PER_ITEM * (len(sel) + n_opus_items) / n_pairs,
           "guideline_tokens_per_shard_pass": guideline_tokens,
           "sonnet_shards": files, "opus_shards": opus_files}
    v4.write_json(INPUT_DIR / "manifest.json", man)
    print(json.dumps({k: man[k] for k in ("n_items_selected", "n_sonnet_shards", "n_opus_shards", "n_opus_items", "items_per_heldout_pair",
                                          "token_estimate_both", "protocol_sha256")}), flush=True)
    return man


def lm_guideline_tokens() -> int:
    import math
    return math.ceil((len(e5.DIR_GUIDELINE.read_text(encoding="utf-8")) + len(e5.REL_GUIDELINE.read_text(encoding="utf-8"))) / 4)


# ============================================================================= labels, gate, feature
def load_label_set(directory: Path, prefix: str, input_dir: Path, shard_numbers: Sequence[int]) -> Tuple[Dict[str, dict], Dict[str, dict]]:
    """(items by id, labels by id) for the given shard numbers; the ids must be exactly the shard ids."""
    items: Dict[str, dict] = {}
    labs: Dict[str, dict] = {}
    for n in shard_numbers:
        sh = json.loads((input_dir / f"ul_sonnet_{n:03d}.json").read_text(encoding="utf-8"))["items"]
        rows = lm.load_annotations(directory / f"{prefix}_{n:03d}.json", sh)
        items.update({it["id"]: it for it in sh})
        labs.update(rows)
    return items, labs


def relevance_kappa(sonnet: Dict[str, dict], opus: Dict[str, dict]) -> Dict[str, Any]:
    ids = sorted(set(sonnet) & set(opus))
    a = np.array([sonnet[i]["relevance"] for i in ids])
    b = np.array([opus[i]["relevance"] for i in ids])
    k = float(ex.kappa(a, b)) if len(ids) else float("nan")
    return {"n_overlap": len(ids), "kappa": k, "accuracy": float(np.mean(a == b)) if len(ids) else None, "gate": KAPPA_GATE, "gate_passed": lm.relevance_gate(k)}


def item_is_benefit(sonnet: Dict[str, dict], opus: Dict[str, dict], iid: str) -> bool:
    """Relevance yes and direction benefit by Sonnet, and Opus (if it labelled the item) does not say relevance != yes."""
    s = sonnet[iid]
    if s["relevance"] != "yes" or s["direction"] != "benefit":
        return False
    return not (iid in opus and opus[iid]["relevance"] != "yes")


def build_B(pair_pids: Dict[str, List[str]], surviving: Dict[str, Set[str]], benefit_keys: Set[Tuple[str, str]], disease: str) -> Dict[str, int]:
    """drug -> B in {0, 1}; benefit_keys = {(drug, pid)} of benefit items. B = 1 if a surviving record carries a benefit item."""
    return {g: int(any((g, p) in benefit_keys for p in surviving.get(g, ()))) for g in pair_pids}


def delta_for(po_candidates: np.ndarray, fraction: float = DELTA_SD_FRACTION) -> float:
    """delta = fraction * SD of the PO scores among the candidate pairs of the disease (0 if fewer than 2 candidates)."""
    return float(fraction * np.std(po_candidates, ddof=1)) if len(po_candidates) > 1 else 0.0


def combine_score(po: np.ndarray, rows: Sequence[int], b: Sequence[int], delta: float) -> np.ndarray:
    out = np.array(po, dtype=float).copy()
    for r, x in zip(rows, b):
        out[r] += delta * x
    return out


# ============================================================================= score
def run_score(root: Path = ROOT) -> Dict[str, Any]:
    if LOCK.exists():
        raise RuntimeError("ultralight HELD-OUT scoring already done once (lock file exists); refusing a second scoring")
    man = json.loads((INPUT_DIR / "manifest.json").read_text(encoding="utf-8"))
    if man["protocol_sha256"] != v4.sha256_file(PROTOCOL):
        raise RuntimeError("protocol changed after the select stage")
    n_s, op = man["n_sonnet_shards"], man["opus_shard_numbers"]
    s_items, s_lab = load_label_set(LABELS_DIR, "sonnet", INPUT_DIR, range(1, n_s + 1))
    o_items, o_lab = load_label_set(LABELS_DIR, "opus", INPUT_DIR, op)
    gate = relevance_kappa(s_lab, o_lab)
    out: Dict[str, Any] = {"label": LABEL, "protocol": str(PROTOCOL.relative_to(root)), "protocol_sha256": v4.sha256_file(PROTOCOL),
                           "delta_sd_fraction": {"value": DELTA_SD_FRACTION, "status": "UNVALIDATED", "conservative_default": 0.0},
                           "quality_gate": gate, "manifest_items": man["n_items_selected"]}
    if not gate["gate_passed"]:
        out["status"] = "not testable: Opus-vs-Sonnet relevance kappa below the gate; B not used; no scoring"
        v4.write_json(OUT_FILE, out)
        print(json.dumps({"gate_passed": False, "kappa": gate["kappa"]}), flush=True)
        return out
    bench = v4.load_benchmark(root)
    dev, ho = bench["split"]["dev"], bench["split"]["heldout"]
    fz = dl._verify(root, bench)
    LOCK.write_text(json.dumps({"scored_at": "first and only ultralight HELD-OUT scoring", "frozen_file_sha256": fz["file_sha256"],
                                "protocol_sha256": v4.sha256_file(PROTOCOL)}) + "\n", encoding="utf-8")
    ctx = ab._context(root, list(ho))                      # HELD-OUT gold is read only here, after the lock exists
    fz, head, mean, Zdev, members_dev, scorers, _ = ab._heads_and_scorers(ctx, dev, "heldout")
    po_fn = scorers["PRIOR_ONLY"]
    e5ctx = e5._context(root)
    doc = json.loads(e5.RETRIEVAL_FILE.read_text(encoding="utf-8"))
    benefit = {(it["disease"], it["drug"], lm.pid_of(it["pmid"])) for i, it in s_items.items() if item_is_benefit(s_lab, o_lab, i)}
    members = v4.members_of(bench, ho)
    dn = sorted(ho)
    res: Dict[str, Dict[Tuple[float, int], Dict[str, float]]] = {d: {} for d in dn}
    cov: Dict[Tuple[float, int], List[int]] = {(f, r): [] for f in FRACTIONS_UL for r in range(N_REPS)}
    for d in dn:
        pairs = doc["diseases"][d]["pairs"]
        rows = {g: e5ctx["index"][g] for g in pairs}
        sv = survivors_by_sample(pairs, e5ctx["rank"][d], e5ctx["index"])
        bk = {(g, p) for (dd, g, p) in benefit if dd == d}
        for f in FRACTIONS_UL:
            fi = FRACTION_INDEX[f]
            for rep in range(N_REPS):
                po = po_fn(d, ctx["fd"][d]["samples"][(fi, rep)]["X"])
                gs = sorted(pairs)
                b = build_B({g: pairs[g]["pids"] for g in gs}, sv[(f, rep)], bk, d)
                delta = delta_for(po[[rows[g] for g in gs]])
                sc = combine_score(po, [rows[g] for g in gs], [b[g] for g in gs], delta)
                res[d][(f, rep)] = {"po": es.auc_mw(po, ctx["y"][d]), "ul": es.auc_mw(sc, ctx["y"][d]), "n_pairs": len(gs), "n_b1": int(sum(b.values())), "delta": delta}
                cov[(f, rep)].append(int(sum(b.values())))
    idx = dl.boot_idx(len(members))
    table: Dict[str, Any] = {}
    for f in FRACTIONS_UL:
        po_m = np.array([[np.mean([res[d][(f, r)]["po"] for r in range(N_REPS)])] for d in dn])[:, 0]
        ul_m = np.array([np.mean([res[d][(f, r)]["ul"] for r in range(N_REPS)]) for d in dn])
        g_po, g_ul = v4.agg_groups(po_m, dn, members), v4.agg_groups(ul_m, dn, members)
        tot_pairs = sum(res[d][(f, 0)]["n_pairs"] for d in dn)
        coverage = float(np.mean([sum(cov[(f, r)]) / tot_pairs for r in range(N_REPS)]))
        diff = es.diff_with_ci(g_ul, g_po, idx)
        lo = diff["ci95"][0]
        informative = coverage >= COVERAGE_MIN
        table[str(f)] = {"auc_prior_only": es.macro_with_ci(g_po, idx), "auc_score": es.macro_with_ci(g_ul, idx), "auc_difference": diff,
                         "coverage_share_pairs_B1": coverage, "informative": informative,
                         "supported": bool(informative and lo is not None and lo > 0)}
    supported = all(t["supported"] for t in table.values())
    out.update({"n_diseases": len(dn), "n_groups": len(members), "table": table,
                "outcome": "supported" if supported else ("not informative" if not all(t["informative"] for t in table.values()) else "not supported"),
                "status": "complete: HELD-OUT scored once; cannot change the pre-registered v4 outcome"})
    v4.write_json(OUT_FILE, out)
    print(json.dumps({"score_done": True, "outcome": out["outcome"]}), flush=True)
    return out


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = {"select": run_select, "score": run_score}.get(cmd)
    if fn is None:
        raise SystemExit("usage: python -m track2.engine_v5_ultralight select|score")
    fn()


if __name__ == "__main__":
    main()
