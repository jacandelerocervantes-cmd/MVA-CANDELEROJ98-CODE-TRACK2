"""Multi-evidence integration of all lines into a ranked list of proposed candidates (hypotheses, low evidence).

    PYTHONPATH=src python -m track2.multi_evidence --out out [--universe data/universe.json] [--gene-block config/multi_evidence_track1_gene_block.json]

Rules are fixed in `protocol/addendum_multi_evidence.md` (post hoc). The integration orders hypotheses; it does not create evidence.
No clinical recommendation, no dose, no claim of efficacy. Only gene symbols, drug names and categorical outputs are read or written.

The case is defined by its output directory (`--out`): every input is read from there under the file names of the MVA case; a missing
input makes its criterion missing for all drugs. Another case (e.g. `out_case_fanconi_anemia_run1`) can be run by changing `--out`.

Criteria (grade in [0, 1] or missing): M mechanism fit, LIT literature, A, B1, B2, C, D, E (BUB1B signature), EP (patient signature),
ENG (engine rank), S (label safety). Age-band label support is a filter and a column, never a criterion.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import re
from pathlib import Path
from typing import Any, Callable

import numpy as np

CRITERIA = ["M", "LIT", "A", "B1", "B2", "C", "D", "E", "EP", "ENG", "S"]
LINES = ["A", "B1", "B2", "C", "D", "E", "EP", "ENG"]
SIGNED = {"LIT", "D", "E", "EP", "S"}
SUPPORT_CUT = 0.9  # a line "supports" a drug when its grade is at least this (rank lines: top 10 %)
MISSING_PENALTY = 0.5
N_DRAWS = 10000
SEED = 20261018
TOP_N = 50
P5_CUT = 0.5
N_PRIMARY = 3
BANDS = {"0-4": (0.0, 5.0), "5-11": (5.0, 12.0), "12-17": (12.0, 18.0)}
BAND_PASS = {"labelled", "partly"}
# Confidence of the five consequences of reduced BubR1 function = (highest verified evidence level of the consequence itself + 1) / 4,
# from the table of report section 3.1 (levels 2, 0, 0, 0, 3). `out/differential_diagnosis.json` holds the gene-level differential, not these values.
CONSEQUENCE_CONFIDENCE = {"1": 0.75, "2": 0.25, "3": 0.25, "4": 0.25, "5": 1.0}
TAG_TO_CONSEQUENCE = {"growth_muscle": "4", "proteotoxic_energy": "2", "senescence_p53": "3", "aneuploidy_tumour": "5", "tumour": "5"}
HARD_EXCLUDE_CLASSES = ("e_hormona_crecimiento",)  # label contraindicates active malignancy; the record documents a rhabdomyosarcoma
MC_RATES = (0.10, 0.20, 0.30)
MC_REPS = 200
MC_DRAWS = 2000


# ---------------------------------------------------------------- small helpers

def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def up(name: str) -> str:
    return re.sub(r"\s+", " ", (name or "").upper().strip())


def clip01(x: float) -> float:
    return max(0.0, min(1.0, x))


# ---------------------------------------------------------------- literature layer: LIT and M

def lit_grade(entries: list[dict[str, Any]]) -> tuple[float | None, dict[str, Any]]:
    """LIT = 0.5 + 0.5 tanh(net / 2); net = sum over (paper, direction) pairs of (level + 1) / 4, help minus harm. None when no entry."""
    if not entries:
        return None, {"best_help_level": None, "n_help_papers": 0, "n_harm_papers": 0, "net": 0.0}
    best: dict[tuple[Any, str], int] = {}
    for i, r in enumerate(entries):
        if r.get("direction") in ("may_help", "may_harm"):
            key = (r.get("pmid") or f"row{i}", r["direction"])
            best[key] = max(best.get(key, -1), int(r.get("level") or 0))
    helps = [lv for (_, d), lv in best.items() if d == "may_help"]
    harms = [lv for (_, d), lv in best.items() if d == "may_harm"]
    net = sum((lv + 1) / 4 for lv in helps) - sum((lv + 1) / 4 for lv in harms)
    return 0.5 + 0.5 * math.tanh(net / 2), {"best_help_level": max(helps) if helps else None, "n_help_papers": len(helps), "n_harm_papers": len(harms), "net": round(net, 3)}


def rule_consequences(mech_text: str, rules: list[dict[str, Any]]) -> set[str]:
    out: set[str] = set()
    for r in rules:
        if r.get("match_any") and not r.get("only_for_list") and not r.get("match_all_groups") and r.get("consequences"):
            if any(k in mech_text for k in r["match_any"]):
                out.update(str(c) for c in r["consequences"])
    return out


def mech_base(entries: list[dict[str, Any]], mech_text: str, rules: list[dict[str, Any]], conf: dict[str, float] = CONSEQUENCE_CONFIDENCE) -> tuple[float | None, dict[str, float]]:
    """M(BUB1B) = max over addressed consequences of confidence x support (1 'may help' entry, 0.5 rule-only or unclear, 0 harm only)."""
    support: dict[str, float] = {}
    for r in entries:
        c = TAG_TO_CONSEQUENCE.get(r.get("consequence") or "")
        if not c:
            continue
        s = 1.0 if r.get("direction") in ("may_help", "selectively_kills_aneuploid_cells") else (0.0 if r.get("direction") == "may_harm" else 0.5)
        support[c] = max(support.get(c, s), s)
    for c in rule_consequences(mech_text, rules):
        support[c] = max(support.get(c, 0.0), 0.5)
    if not support:
        return None, {}
    return max(conf[c] * s for c, s in support.items() if c in conf), {c: support[c] for c in sorted(support)}


def track1_bonus(entries: list[dict[str, Any]], mech_text: str, block: dict[str, Any] | None) -> tuple[float, list[str]]:
    """T = max over non-primary genes of tier weight x relevance, in [0, 1]; also the genes that bear on the drug."""
    if not block:
        return 0.0, []
    text = " ".join(f"{r.get('finding') or ''} {r.get('system') or ''} {r.get('sentence') or ''}" for r in entries).lower()
    tags = {r.get("consequence") for r in entries}
    best, genes = 0.0, []
    for g in block["genes"]:
        if g["gene"] == block.get("primary_gene") or g["tier_weight"] <= 0:
            continue
        rel = 0.0
        for b in g.get("bears_on", []):
            hit = (b.get("mechanism_regex") and re.search(b["mechanism_regex"], mech_text)) or (b.get("literature_regex") and re.search(b["literature_regex"], text)) \
                or (b.get("literature_tags") and tags & set(b["literature_tags"]))
            if hit:
                rel = max(rel, float(b["relevance"]))
        if rel > 0:
            genes.append(g["gene"])
            best = max(best, g["tier_weight"] * rel)
    return best, genes


def mech_grade(base: float | None, bonus: float, cap: float) -> float | None:
    if base is None and bonus <= 0:
        return None
    return min(1.0, (base or 0.0) + cap * bonus)


# ---------------------------------------------------------------- line grades (each returns {NAME: grade or None}, extra info)

def rank_grade(rank: float, n: int) -> float:
    return clip01(1 - (rank - 1) / n)


def line_A(out: Path, names: dict[str, str]) -> tuple[dict[str, float], list[str], int]:
    q, nq = load_json(out / "kg_query_e5.json"), load_json(out / "kg_nullquery_e5.json")
    if not q:
        return {}, [], 0
    best: dict[str, int] = {}
    for did, r in q["rank_combined"].items():
        if did in names:
            best[names[did]] = min(r, best.get(names[did], r))
    n = len(best)
    frac = {up(k): v for k, v in (nq or {}).get("fraction_of_random_queries_with_drug_in_top50", {}).items()}
    order = sorted(best, key=lambda d: (best[d], d))
    top = order[:TOP_N]
    grades = {}
    for d in best:
        g = rank_grade(best[d], n)
        if d in frac:
            g *= 1 - frac[d]
        elif d in top:
            g *= 0.5
        grades[d] = g
    return grades, top, n


def line_B1(out: Path, names: dict[str, str]) -> tuple[dict[str, float], list[str]]:
    p, t = load_json(out / "phenotype_proximity_e5.json"), load_json(out / "phenotype_threshold_scores_e5.json")
    if not p:
        return {}, []
    scored = {names[d] for d in ((t or {}).get("min2_terms", {}).get("700", {})) if d in names}
    conf = int((p.get("rule") or {}).get("configurations", 9))
    listed = {up(d["name"]): d["robustness_count"] / conf for d in p["drugs"]}
    grades = {d: 0.0 for d in scored}
    grades.update(listed)
    top = sorted(listed, key=lambda d: (-listed[d], d))[:TOP_N]
    return grades, top


def line_B2(out: Path) -> dict[str, float]:
    s = load_json(out / "track1_scenarios.json")
    if not s:
        return {}
    schemes = s["variants"]["two_record_links_no_predicted"]["scheme_robust"]
    counts: dict[str, int] = {}
    for lst in schemes.values():
        for n in lst:
            counts[up(n)] = counts.get(up(n), 0) + 1
    return {d: c / len(schemes) for d, c in counts.items()}


def line_C(out: Path, names: dict[str, str]) -> dict[str, float]:
    v = load_json(out / "vulnerability_e6.json")
    sup = (v or {}).get("supported_drugs") or {}
    if not v or not sup:
        return {}  # no drug-level result: missing for all
    keys = sup if isinstance(sup, list) else list(sup)
    return {up(names.get(k, k)): 1.0 for k in keys}


def line_D(out: Path) -> tuple[dict[str, float], list[str]]:
    s = load_json(out / "signed_effect.json")
    if not s:
        return {}, []
    grades, strong = {}, []
    for rec in s["drugs"].values():
        n = up(rec["name"])
        if rec["call"] == "strengthens_checkpoint":
            grades[n] = 1.0
            strong.append((rec.get("net", 0), n))
        elif rec["call"] == "weakens_checkpoint":
            grades[n] = 0.0
    return grades, [n for _, n in sorted(strong, key=lambda x: (-x[0], x[1]))[:TOP_N]]


def line_E(out: Path) -> tuple[dict[str, float], dict[str, float], list[str], list[str]]:
    e = load_json(out / "line_e_signature_reversal.json")
    if not e:
        return {}, {}, [], []
    res: dict[str, dict[str, float]] = {"bub1b": {}, "P": {}}
    tops: dict[str, list[tuple[float, str]]] = {"bub1b": [], "P": []}
    for r in e["table_all_drugs_n150"]:
        for key in ("bub1b", "P"):
            s = r.get(key)
            if not s or s.get("q_rev") is None:
                continue
            n = up(r["name"])
            if s.get("hazard"):
                res[key][n] = 0.0
            else:
                res[key][n] = clip01(1 - s["q_rev"])
                tops[key].append((s["p_rev"], n))
    return res["bub1b"], res["P"], [n for _, n in sorted(tops["bub1b"])[:TOP_N]], [n for _, n in sorted(tops["P"])[:TOP_N]]


def line_ENG(out: Path) -> tuple[dict[str, float], list[str], dict[str, bool]]:
    g = load_json(out / "engine" / "apply_mva.json")
    if not g:
        return {}, [], {}
    ranks: dict[str, int] = {}
    n = 1293
    fda: dict[str, bool] = {}
    for r in g.get("ranking_top20_fda_small_molecule", []):
        ranks[up(r["name"])] = int(r["position"])
        fda[up(r["name"])] = True
    for k, v in g.get("nine_report_drugs_positions", {}).items():
        n = int(v.get("of_fda_small_molecule") or n)
        if v.get("rank_fda_small_molecule"):
            ranks[up(k)] = int(v["rank_fda_small_molecule"])
        fda[up(k)] = bool(v.get("fda_label"))
    return {d: rank_grade(r, n) for d, r in ranks.items()}, sorted(ranks, key=lambda d: (ranks[d], d))[:TOP_N], fda


def line_ENG_full(rec: dict[str, Any]) -> tuple[dict[str, float], list[str], dict[str, bool]]:
    """ENG from the stored full ranking (`out/engine/apply_mva_full_ranking.json`): every FDA small molecule has a rank."""
    drugs = {d: v for d, v in rec["drugs"].items() if v.get("rank_fda_small_molecule")}
    n = int(rec.get("n_fda_small_molecule") or len(drugs))
    ranks = {d: int(v["rank_fda_small_molecule"]) for d, v in drugs.items()}
    return {d: rank_grade(r, n) for d, r in ranks.items()}, sorted(ranks, key=lambda d: (ranks[d], d))[:TOP_N], {d: True for d in ranks}


# ---------------------------------------------------------------- approval, safety, age

class CacheMiss(Exception):
    """Raised on a cache miss; not a LookupError, because `fda_approved.lookup` reads LookupError as "no record"."""


def cache_only_getter(cache_dir: Path) -> Callable[[str], Any]:
    import hashlib

    def get(url: str) -> Any:
        for d in (cache_dir, cache_dir.parent / "cache_engine", cache_dir.parent / "cache_empirical"):
            path = d / f"{hashlib.sha256(url.encode()).hexdigest()[:24]}.json"
            if path.exists():
                return json.loads(path.read_text(encoding="utf-8"))
        raise CacheMiss(url)
    return get


def approval_by_cache(name: str, cache_dir: Path) -> bool | None:
    """True/False from the local openFDA cache through `fda_approved.lookup`; None when the cache has no answer (offline)."""
    from track2.fda_approved import lookup
    try:
        return bool(lookup(name, get=cache_only_getter(cache_dir))["eligible"])
    except CacheMiss:
        return None


def safety_grade(name: str, universe_rec: dict[str, Any] | None, proposal_safety: dict[str, Any], caution_names: set[str]) -> float:
    if name in caution_names:
        return 0.0
    if (universe_rec or {}).get("black_box_warning") or name in {up(k) for k in (proposal_safety or {}).get("entries", {})}:
        return 0.5
    return 1.0


def _measure(intervals: list[list[float]], lo: float, hi: float) -> list[tuple[float, float]]:
    out = []
    for a, b in intervals:
        a2, b2 = max(a, lo), min(b, hi)
        if b2 > a2:
            out.append((a2, b2))
    return out


def band_state(rec: dict[str, Any], lo: float, hi: float) -> str:
    """Label state in the band [lo, hi) from the reviewed label intervals L and the age-specific cautions H (empirical age support)."""
    if rec["L"]["state"] != "known":
        return "verify"
    grid = np.arange(lo, hi, 0.05) + 0.025
    ok = np.zeros(len(grid), dtype=bool)
    for a, b in rec["L"]["intervals"]:
        ok |= (grid >= a) & (grid < b)
    for a, b in rec.get("H_intervals_years", []):
        ok &= ~((grid >= a) & (grid < b))
    frac = float(ok.mean())
    return "labelled" if frac >= 0.999 else ("partly" if frac > 0 else "not_labelled")


def engine_band(state: str | None) -> str | None:
    if state is None:
        return None
    return "partly" if state.startswith("partly") else ("labelled" if state == "age_not_stated" else state)


def age_table(out: Path) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    """Band states per drug from the age-support file, falling back to the engine reading; and the disagreements between the two."""
    a = load_json(out / "empirical" / "age_support.json")
    g = load_json(out / "engine" / "apply_mva.json") or {}
    table: dict[str, dict[str, str]] = {}
    for rec in (a or {}).get("drugs", []):
        table[up(rec["name"])] = {b: band_state(rec, lo, hi) for b, (lo, hi) in BANDS.items()}
    eng: dict[str, dict[str, str]] = {}
    for k, v in g.get("label_support", {}).items():
        eng[up(k)] = {"0-4": engine_band(v["bands"].get("0_4")), "5-11": engine_band(v["bands"].get("5_11")), "12-17": engine_band(v["bands"].get("12_17"))}
    diffs = []
    for d, bands in eng.items():
        if d in table:
            for b in BANDS:
                if bands[b] != table[d][b]:
                    diffs.append({"drug": d, "band": b, "age_support": table[d][b], "engine": bands[b] or "none"})
        else:
            table[d] = {b: bands[b] or "not_assessed" for b in BANDS}
    return table, sorted(diffs, key=lambda x: (x["drug"], x["band"]))


# ---------------------------------------------------------------- SMAA

def draw_weights(rng: np.random.Generator, n: int, k: int) -> np.ndarray:
    return rng.dirichlet(np.ones(k), size=n)


def scores(G: np.ndarray, W: np.ndarray) -> np.ndarray:
    """score = weighted mean over available criteria x (1 - penalty x missing weight share); G has NaN for missing; W rows sum to 1."""
    avail = ~np.isnan(G)
    num = W @ np.nan_to_num(G).T
    den = W @ avail.T.astype(float)
    mean = np.divide(num, den, out=np.zeros_like(num), where=den > 1e-12)
    return mean * (1 - MISSING_PENALTY * np.clip(1 - den, 0, 1))


def smaa(G: np.ndarray, W: np.ndarray) -> dict[str, np.ndarray]:
    """Rank acceptability and central weights; ordinal ranks, ties broken by row order (rows are sorted by drug name)."""
    S = scores(G, W)
    order = np.argsort(-S, axis=1, kind="stable")
    ranks = np.empty_like(order)
    rows = np.arange(S.shape[0])[:, None]
    ranks[rows, order] = np.arange(1, S.shape[1] + 1)[None, :]
    n = G.shape[0]
    avail = ~np.isnan(G)
    den = W @ avail.T.astype(float)
    central = np.zeros((n, W.shape[1]))
    best_rank = ranks.min(axis=0)
    for i in range(n):
        sel = ranks[:, i] == 1 if (ranks[:, i] == 1).any() else ranks[:, i] == best_rank[i]
        central[i] = W[sel].mean(axis=0)
    return {"p1": (ranks == 1).mean(0), "p3": (ranks <= 3).mean(0), "p5": (ranks <= 5).mean(0), "p10": (ranks <= 10).mean(0), "central": central,
            "confidence": (1 - np.clip(1 - den, 0, 1)).mean(0), "mean_rank": ranks.mean(0), "best_rank": best_rank}


def decide(names: list[str], res: dict[str, np.ndarray]) -> tuple[list[str], list[str]]:
    """Proposed candidates: P(top 5) >= 0.5; primary: the top 3 of them by P(top 5), ties by P(top 3) then name."""
    prop = [i for i in range(len(names)) if res["p5"][i] >= P5_CUT]
    prop.sort(key=lambda i: (-res["p5"][i], -res["p3"][i], names[i]))
    return [names[i] for i in prop], [names[i] for i in prop[:N_PRIMARY]]


# ---------------------------------------------------------------- assembling the case

def build_case(out: Path, universe: dict[str, Any], rules_cfg: dict[str, Any] | None, gene_block: dict[str, Any] | None, cache_dir: Path | None = None,
               approval_fn: Callable[[str], bool | None] | None = None, cap: float | None = None, extra_approved: set[str] | None = None,
               engine_full: dict[str, Any] | None = None) -> dict[str, Any]:
    from track2 import literature_summary as LS

    by_name = {up(d["name"]): d for d in universe["drugs"]}
    names = {d["drug_id"]: up(d["name"]) for d in universe["drugs"]}
    cap = float((gene_block or {}).get("cap", 0.0)) if cap is None else cap
    lit_dir = out / "literature_hypotheses"
    entries = LS.collect(lit_dir) if (lit_dir / "verification.json").exists() else {}
    entries = {up(k): v for k, v in entries.items()}
    rules = (rules_cfg or {}).get("rules", [])
    A, topA, nA = line_A(out, names)
    B1, topB1 = line_B1(out, names)
    B2 = line_B2(out)
    C = line_C(out, names)
    D, topD = line_D(out)
    E, EP, topE, topEP = line_E(out)
    ENG, topENG, eng_fda = line_ENG(out)
    if engine_full:
        ENG, topENG, eng_fda = line_ENG_full(engine_full)
    tops = {"A": topA, "B1": topB1, "B2": sorted(B2), "D": topD, "E": topE, "EP": topEP, "ENG": topENG}
    caution = load_json(out / "caution_list.json") or {}
    caution_names = {up(f["name"]) for c in caution.get("clases", {}).values() for f in c.get("farmacos", [])}
    hard_class = {up(f["name"]) for k in HARD_EXCLUDE_CLASSES for f in caution.get("clases", {}).get(k, {}).get("farmacos", [])}
    hazard = {n for n, v in entries.items() if any(r["direction"] == "selectively_kills_aneuploid_cells" and r["level"] >= 1 for r in v["entries"])}
    proposal_safety = load_json(Path(out).parent / "config" / "proposal_safety.json") or load_json(Path("config/proposal_safety.json")) or {}
    approved_files = {d for d, v in (load_json(out / "fda_approval_check.json") or {}).get("drugs", {}).items() if v.get("confirmed")} | set(extra_approved or ())
    e_json = load_json(out / "line_e_signature_reversal.json") or {}
    approved_e = {up(r["name"]) for r in e_json.get("table_all_drugs_n150", []) if r.get("fda_label")}
    approved_eng = {d for d, ok in eng_fda.items() if ok}

    candidates = set(entries) | {d for lst in tops.values() for d in lst}
    pool, not_checked, not_eligible, not_small = [], [], [], []
    for d in sorted(candidates):
        rec = by_name.get(d)
        if not rec or rec.get("molecule_type") != "Small molecule":
            not_small.append(d)
            continue
        ok = d in approved_files or d in approved_e or d in approved_eng
        if not ok and approval_fn is not None:
            r = approval_fn(d)
            ok = bool(r)
            if r is False:
                not_eligible.append(d)
                continue
        if ok:
            pool.append(d)
        else:
            not_checked.append(d)
    return {"names": names, "by_name": by_name, "entries": entries, "rules": rules, "gene_block": gene_block, "cap": cap, "lines": {"A": A, "B1": B1, "B2": B2, "C": C, "D": D, "E": E, "EP": EP, "ENG": ENG},
            "tops": tops, "caution_names": caution_names, "hard_class": hard_class, "hazard": hazard, "proposal_safety": proposal_safety, "pool": pool, "not_checked": not_checked, "not_eligible": not_eligible,
            "not_small_or_unknown": not_small, "n_A": nA}


def mech_text(rec: dict[str, Any] | None) -> str:
    return " ".join(f"{m.get('mechanism', '')} {m.get('action', '')}" for m in (rec or {}).get("mechanisms", [])).lower()


def drug_rows(case: dict[str, Any], entries: dict[str, dict[str, list[dict[str, Any]]]] | None = None, cap: float | None = None) -> dict[str, dict[str, Any]]:
    """Criteria grades per pool drug. `entries` can be replaced (label-error Monte Carlo) so that LIT and M are recomputed."""
    ent = case["entries"] if entries is None else entries
    cap = case["cap"] if cap is None else cap
    rows = {}
    for d in case["pool"]:
        rec = case["by_name"].get(d)
        e = ent.get(d, {"entries": []})["entries"]
        txt = mech_text(rec)
        lit, litinfo = lit_grade(e)
        base, cons = mech_base(e, txt, case["rules"])
        bonus, genes = track1_bonus(e, txt, case["gene_block"])
        g = {"M": mech_grade(base, bonus, cap), "LIT": lit}
        for k in ("A", "B1", "B2", "C", "D", "E", "EP", "ENG"):
            g[k] = case["lines"][k].get(d)
        g["S"] = safety_grade(d, rec, case["proposal_safety"], case["caution_names"])
        rows[d] = {"grades": g, "lit": litinfo, "consequences": cons, "track1_genes": genes, "track1_bonus": round(bonus, 3), "mech_base": base}
    return rows


def matrix(rows: dict[str, dict[str, Any]], order: list[str], criteria: list[str] = CRITERIA) -> np.ndarray:
    return np.array([[np.nan if rows[d]["grades"][c] is None else rows[d]["grades"][c] for c in criteria] for d in order], dtype=float)


def direction(c: str, g: float | None) -> str:
    if g is None:
        return "missing"
    if c in SIGNED:
        return "supportive" if g > 0.5 else ("contrary" if g < 0.5 else "neutral")
    return "supportive" if g >= SUPPORT_CUT else "neutral"


def tier_text(row: dict[str, Any]) -> str:
    g, lit = row["grades"], row["lit"]
    sup = sum(1 for c in LINES if g[c] is not None and g[c] >= SUPPORT_CUT)
    con = sum(1 for c in ("D", "E", "EP") if g[c] is not None and g[c] == 0.0)
    if g["LIT"] is None:
        head = "no verified literature entry"
    elif lit["best_help_level"] is None:
        head = f"literature without a supportive finding ({lit['n_harm_papers']} harm papers)"
    else:
        head = f"level {lit['best_help_level']} literature ({lit['n_help_papers']} help, {lit['n_harm_papers']} harm papers)"
    return f"{head}; {sup} lines supporting, {con} contrary"


def bands_for(d: str, age: dict[str, dict[str, str]]) -> dict[str, str]:
    return age.get(d, {b: "not_assessed" for b in BANDS})


# ---------------------------------------------------------------- the analysis

def analyse(case: dict[str, Any], age: dict[str, dict[str, str]], n_draws: int = N_DRAWS, seed: int = SEED, mc_reps: int = MC_REPS, mc_draws: int = MC_DRAWS) -> dict[str, Any]:
    excluded: dict[str, str] = {}
    for d in case["pool"]:
        if d in case["hazard"]:
            excluded[d] = "hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)"
        elif d in case["hard_class"]:
            excluded[d] = "label contraindicates the hypothesised use (growth-hormone/IGF-1 class: active malignancy; the record documents a rhabdomyosarcoma)"
    rows = drug_rows(case)
    order = sorted(d for d in rows if d not in excluded)
    G = matrix(rows, order)
    rng = np.random.default_rng(seed)
    W = draw_weights(rng, n_draws, len(CRITERIA))
    base = smaa(G, W)
    proposed, primary = decide(order, base)
    drugs = {}
    for i, d in enumerate(order):
        r = rows[d]
        drugs[d] = {"grades": {c: (None if r["grades"][c] is None else round(r["grades"][c], 4)) for c in CRITERIA},
                    "directions": {c: direction(c, r["grades"][c]) for c in CRITERIA},
                    "p_rank1": round(float(base["p1"][i]), 4), "p_top3": round(float(base["p3"][i]), 4), "p_top5": round(float(base["p5"][i]), 4), "p_top10": round(float(base["p10"][i]), 4),
                    "mean_rank": round(float(base["mean_rank"][i]), 2), "best_rank": int(base["best_rank"][i]),
                    "central_weight": {c: round(float(base["central"][i][j]), 3) for j, c in enumerate(CRITERIA)}, "confidence_factor": round(float(base["confidence"][i]), 3),
                    "evidence_tier": tier_text(r), "literature": r["lit"], "consequences_addressed": r["consequences"], "track1_genes_bearing": r["track1_genes"],
                    "track1_bonus": r["track1_bonus"], "age_bands": bands_for(d, age), "proposed": d in proposed, "primary": d in primary,
                    "n_missing": int(np.isnan(G[i]).sum())}
    # age bands: same weights, ranking restricted to the drugs that pass the band filter
    age_out = {}
    for b in list(BANDS) + ["unknown_age_union"]:
        keep = [d for d in order if (any(bands_for(d, age)[x] in BAND_PASS for x in BANDS) if b == "unknown_age_union" else bands_for(d, age)[b] in BAND_PASS)]
        res = smaa(matrix(rows, keep), W) if keep else None
        prop, prim = decide(keep, res) if res else ([], [])
        glob = [d for d in proposed if d in keep]
        age_out[b] = {"filter_pass": keep, "proposed": glob, "primary": [d for d in primary if d in keep], "restricted_pool_proposed": prop, "restricted_pool_primary": prim,
                      "p_top5": {d: round(float(res["p5"][i]), 4) for i, d in enumerate(keep)} if res else {},
                      "outside_filter": {d: (bands_for(d, age)[b] if b != "unknown_age_union" else "/".join(bands_for(d, age)[x] for x in BANDS)) for d in order if d not in keep and d in proposed},
                      "verify_or_unassessed_in_global_top10": [d for d in sorted(order, key=lambda x: (-drugs[x]["p_top5"], x))[:10] if d not in keep]}
    sens = sensitivity(case, rows, order, excluded, proposed, base, age, n_draws, seed, mc_reps, mc_draws)
    return {"order": order, "drugs": drugs, "proposed": proposed, "primary": primary, "excluded": excluded, "age_bands": age_out, "sensitivity": sens}


def jaccard(a: set[str], b: set[str]) -> float:
    return 1.0 if not a and not b else len(a & b) / len(a | b)


def sensitivity(case: dict[str, Any], rows: dict[str, dict[str, Any]], order: list[str], excluded: dict[str, str], proposed: list[str], base: dict[str, np.ndarray],
                age: dict[str, dict[str, str]], n_draws: int, seed: int, mc_reps: int, mc_draws: int) -> dict[str, Any]:
    G = matrix(rows, order)
    base_set = set(proposed)
    out: dict[str, Any] = {}

    def summary(prop: list[str], res: dict[str, np.ndarray] | None = None) -> dict[str, Any]:
        s = {"proposed": prop, "jaccard_with_baseline": round(jaccard(set(prop), base_set), 3), "baseline_retained": sorted(set(prop) & base_set)}
        if res is not None:
            top = sorted(range(len(order)), key=lambda i: (-res["p5"][i], order[i]))[:5]
            s["top5_by_p_top5"] = [{"drug": order[i], "p_top5": round(float(res["p5"][i]), 3)} for i in top]
        return s

    loo = {}
    for k, c in enumerate(CRITERIA):
        cols = [j for j in range(len(CRITERIA)) if j != k]
        res = smaa(G[:, cols], draw_weights(np.random.default_rng(seed + 1 + k), n_draws, len(cols)))
        loo[c] = summary(decide(order, res)[0], res)
    out["leave_one_criterion_out"] = loo
    # equal weights: one vector
    eq = scores(G, np.full((1, len(CRITERIA)), 1 / len(CRITERIA)))[0]
    eq_order = sorted(range(len(order)), key=lambda i: (-eq[i], order[i]))
    top5 = eq_order[:5]
    ties = sorted({order[i] for i in range(len(order)) if i not in top5 and abs(eq[i] - eq[top5[-1]]) < 1e-12})
    out["equal_weights"] = {"top5": [{"drug": order[i], "score": round(float(eq[i]), 4)} for i in top5], "ties_with_fifth": ties,
                            "jaccard_top5_with_baseline_proposed": round(jaccard({order[i] for i in top5}, base_set), 3), "top5_in_baseline_proposed": sorted({order[i] for i in top5} & base_set)}
    # literature-only weights: SMAA over LIT and M only
    cols = [CRITERIA.index("LIT"), CRITERIA.index("M")]
    res = smaa(G[:, cols], draw_weights(np.random.default_rng(seed + 50), n_draws, 2))
    out["literature_only_weights"] = summary(decide(order, res)[0], res)
    # Track 1 non-BUB1B genes off (cap 0)
    rows0 = drug_rows(case, cap=0.0)
    res = smaa(matrix(rows0, order), draw_weights(np.random.default_rng(seed), n_draws, len(CRITERIA)))
    prop0 = decide(order, res)[0]
    out["track1_gene_block_off"] = {**summary(prop0, res), "p_top5_change": {d: round(float(base["p5"][i] - res["p5"][i]), 4) for i, d in enumerate(order) if abs(base["p5"][i] - res["p5"][i]) > 1e-9},
                                    "drugs_with_gene_bonus": {d: {"genes": rows[d]["track1_genes"], "bonus": rows[d]["track1_bonus"]} for d in order if rows[d]["track1_bonus"] > 0}}
    # Monte Carlo on label errors
    from track2.empirical_label_robustness import corrupt
    mc = {}
    for rate in MC_RATES:
        rng = random.Random(f"{seed}-{rate}")
        stay = {d: 0 for d in proposed}
        sizes, jac = [], []
        for rep in range(mc_reps):
            ent = corrupt({d: v for d, v in case["entries"].items()}, rate, rng)
            r2 = drug_rows(case, entries=ent)
            res = smaa(matrix(r2, order), draw_weights(np.random.default_rng(seed * 7 + rep + int(rate * 100)), mc_draws, len(CRITERIA)))
            prop = set(decide(order, res)[0])
            for d in prop & base_set:
                stay[d] += 1
            sizes.append(len(prop))
            jac.append(jaccard(prop, base_set))
        mc[f"{round(rate * 100)}"] = {"replicates": mc_reps, "draws_per_replicate": mc_draws, "retained_fraction_by_baseline_candidate": {d: round(stay[d] / mc_reps, 3) for d in sorted(stay)},
                                      "expected_baseline_retained": round(sum(stay.values()) / mc_reps, 3), "n_baseline": len(base_set), "mean_set_size": round(float(np.mean(sizes)), 2),
                                      "mean_jaccard": round(float(np.mean(jac)), 3)}
    out["label_error_monte_carlo"] = mc
    return out


NINE = {"SIROLIMUS": "A constitutional", "EVEROLIMUS": "A constitutional (sibling)", "METFORMIN": "A constitutional", "PHENYLBUTANOIC ACID": "A constitutional (sodium phenylbutyrate)",
        "TRAMETINIB": "B tumour context", "SELUMETINIB": "B tumour context (sibling)", "BORTEZOMIB": "C hazard", "CHLOROQUINE": "C hazard", "BOSUTINIB": "C hazard"}


def comparison(case: dict[str, Any], result: dict[str, Any], out: Path) -> dict[str, Any]:
    fl = load_json(out / "four_lines_review2.json") or {}
    old = {k: v.get("n_candidates") for k, v in (fl.get("readings") or {}).items()}
    order = result["order"]
    ranked = sorted(order, key=lambda d: (-result["drugs"][d]["p_top5"], -result["drugs"][d]["p_top3"], d))
    nine = {}
    for d, grp in NINE.items():
        if d in result["excluded"]:
            nine[d] = {"group": grp, "status": "hard-excluded", "reason": result["excluded"][d]}
        elif d in result["drugs"]:
            r = result["drugs"][d]
            nine[d] = {"group": grp, "status": "proposed" if r["proposed"] else "in pool, not proposed", "rank_by_p_top5": ranked.index(d) + 1, "of": len(ranked), "p_top5": r["p_top5"], "p_top10": r["p_top10"]}
        else:
            nine[d] = {"group": grp, "status": "not in the pool"}
    return {"old_rule_candidates": old, "old_rule_text": "at least two usable computational lines", "new_proposed": result["proposed"], "nine_report_drugs": nine,
            "proposed_not_in_report": [d for d in result["proposed"] if d not in NINE]}


def run(out: Path, universe_path: Path, rationale_path: Path, gene_block_path: Path | None, cache_dir: Path | None, n_draws: int = N_DRAWS, mc_reps: int = MC_REPS, mc_draws: int = MC_DRAWS) -> dict[str, Any]:
    universe = load_json(universe_path)
    case = build_case(out, universe, load_json(rationale_path), load_json(gene_block_path) if gene_block_path else None,
                      approval_fn=(lambda n: approval_by_cache(n, cache_dir)) if cache_dir else None)
    age, diffs = age_table(out)
    result = analyse(case, age, n_draws=n_draws, mc_reps=mc_reps, mc_draws=mc_draws)
    comp = comparison(case, result, out)
    coverage = {b: sum(1 for d in result["order"] if bands_for(d, age)[b] != "not_assessed") for b in BANDS}
    gb = case["gene_block"] or {}
    return {
        "label": "hypotheses for follow-up; low evidence; the integration orders hypotheses and creates no evidence; no clinical recommendation, dose or claim of efficacy",
        "protocol": "protocol/addendum_multi_evidence.md (post hoc)", "seed": SEED, "n_draws": n_draws, "criteria": CRITERIA, "missing_penalty": MISSING_PENALTY, "proposal_rule": f"P(top 5) >= {P5_CUT}, not hard-excluded; primary = top {N_PRIMARY}",
        "pool": {"size": len(case["pool"]), "n_after_exclusion": len(result["order"]), "drugs": sorted(case["pool"]), "approval_not_checked": sorted(case["not_checked"]), "approval_checked_not_eligible": sorted(case["not_eligible"]), "n_not_small_molecule_or_not_in_universe": len(case["not_small_or_unknown"]),
                 "n_by_source": {k: len(set(v) & set(case["pool"])) for k, v in case["tops"].items()} | {"literature": len(set(case["entries"]) & set(case["pool"]))}},
        "excluded": result["excluded"], "proposed": result["proposed"], "primary": result["primary"], "drugs": result["drugs"],
        "age_bands": result["age_bands"], "age_label_coverage_n_assessed": coverage, "age_reading_disagreements_age_support_vs_engine": diffs,
        "sensitivity": result["sensitivity"], "comparison": comp,
        "track1_gene_block": {"cap": case["cap"], "genes": [{k: g[k] for k in ("gene", "role", "tier", "tier_weight")} for g in gb.get("genes", [])], "tests": gb.get("tests", {}),
                              "statement": "These are hypotheses to test (karyotype or single-cell sequencing for the megabase events; function tests for IGBP1), not findings."},
        "consequence_confidence": CONSEQUENCE_CONFIDENCE,
        "data_gaps": ["engine ranking stored only for the top 20 and the nine report drugs: ENG is missing for the other pool drugs and 'top 50 of the engine' is the stored top 20",
                      "line C gives no drug-level result: missing for every drug",
                      "confidences of the five consequences are taken from report section 3.1 (out/differential_diagnosis.json holds the gene-level differential only)",
                      "age-band label support is assessed for the drugs of the label table only; the rest are 'not_assessed'"],
    }


# ---------------------------------------------------------------- markdown

def fmt(x: float | None) -> str:
    return "-" if x is None else f"{x:.2f}"


def to_markdown(r: dict[str, Any]) -> str:
    L = [f"# Multi-evidence integration of proposed candidates (case output: `out`)", "", f"*{r['label']}.*", "", f"Protocol: `{r['protocol']}`. Seed {r['seed']}, {r['n_draws']} weight vectors over {len(r['criteria'])} criteria; rule: {r['proposal_rule']}.", ""]
    p = r["pool"]
    L += ["## 1. Pool", "", f"{p['size']} approved small molecules entered the pool ({p['n_after_exclusion']} after hard exclusions). Members by source (a drug can be in several): "
          + ", ".join(f"{k} {v}" for k, v in sorted(p["n_by_source"].items())) + f". {len(p['approval_not_checked'])} candidates had no approval record in the project outputs or the local openFDA cache (offline: not in the pool, so the pool is incomplete for them) and {len(p['approval_checked_not_eligible'])} were checked and are not eligible.", ""]
    if r["excluded"]:
        L += ["Hard exclusions: " + "; ".join(f"{d} ({why})" for d, why in sorted(r["excluded"].items())) + ".", ""]
    ds = r["drugs"]
    L += ["## 2. Proposed candidates (hypotheses, low evidence)", ""]
    if not r["proposed"]:
        L += ["No drug reached P(top 5) >= 0.5: zero proposed candidates.", ""]
    else:
        L += ["| Drug | Status | P(top 5) | P(top 3) | P(rank 1) | P(top 10) | Confidence | Evidence tier |", "|---|---|--:|--:|--:|--:|--:|---|"]
        for d in r["proposed"]:
            x = ds[d]
            L.append(f"| {d.title()} | {'PRIMARY' if x['primary'] else 'proposed'} | {x['p_top5']:.3f} | {x['p_top3']:.3f} | {x['p_rank1']:.3f} | {x['p_top10']:.3f} | {x['confidence_factor']:.2f} | {x['evidence_tier']} |")
        L += ["", "### Evidence matrix of the shortlist (grade in [0, 1]; `-` missing; arrows: + supportive, = neutral, x contrary)", "",
              "| Drug | " + " | ".join(r["criteria"]) + " | Age 0-4 | 5-11 | 12-17 |", "|---|" + "--:|" * len(r["criteria"]) + "---|---|---|"]
        mark = {"supportive": "+", "neutral": "=", "contrary": "x", "missing": ""}
        for d in r["proposed"]:
            x = ds[d]
            cells = [f"{fmt(x['grades'][c])}{mark[x['directions'][c]]}" if x["grades"][c] is not None else "-" for c in r["criteria"]]
            L.append(f"| {d.title()} | " + " | ".join(cells) + " | " + " | ".join(x["age_bands"][b] for b in BANDS) + " |")
        L += ["", "Central weight vectors (weights that put the drug first, or at its best rank):", ""]
        for d in r["proposed"]:
            cw = ds[d]["central_weight"]
            L.append(f"- {d.title()}: " + ", ".join(f"{c} {cw[c]:.2f}" for c in r["criteria"]))
        L += [""]
    top = sorted(ds, key=lambda d: (-ds[d]["p_top5"], -ds[d]["p_top3"], d))[:15]
    L += ["## 3. Rank acceptability, top 15 by P(top 5)", "", "| Drug | P(top 5) | P(top 3) | P(rank 1) | P(top 10) | Mean rank | Missing criteria | Proposed | Evidence tier | Age 0-4 / 5-11 / 12-17 |", "|---|--:|--:|--:|--:|--:|--:|---|---|---|"]
    for d in top:
        x = ds[d]
        L.append(f"| {d.title()} | {x['p_top5']:.3f} | {x['p_top3']:.3f} | {x['p_rank1']:.3f} | {x['p_top10']:.3f} | {x['mean_rank']:.1f} | {x['n_missing']} | {'yes' if x['proposed'] else 'no'} | {x['evidence_tier']} | {' / '.join(x['age_bands'][b] for b in BANDS)} |")
    L += ["", "## 4. Age bands (label support is a filter and a column, not evidence of efficacy)", "",
          "Number of pool drugs with a label reading per band: " + ", ".join(f"{b} {n}" for b, n in r["age_label_coverage_n_assessed"].items()) + ". Unknown age is shown as the union; the age is never assumed.", ""]
    L += ["Headline reading (conservative): a drug is shown for a band only if it is proposed in the full pool and passes the band filter. The restricted-pool reading re-ranks only the drugs that pass the filter; with a small pool the top-5 test is permissive (declared in the addendum, section 10) and it is shown for information only.", ""]
    for b, a in r["age_bands"].items():
        L.append(f"- **{b}**: {len(a['filter_pass'])} pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: " + (", ".join(f"{d.title()}" for d in a["proposed"]) or "none")
                 + ". Primary passing: " + (", ".join(d.title() for d in a["primary"]) or "none") + ". "
                 + ("Proposed but outside the filter: " + ", ".join(f"{d.title()} [{s}]" for d, s in sorted(a["outside_filter"].items())) + ". " if a["outside_filter"] else "")
                 + "Restricted-pool reading (information only): " + (", ".join(f"{d.title()} ({a['p_top5'][d]:.2f})" for d in a["restricted_pool_proposed"]) or "none") + ".")
    if r["age_reading_disagreements_age_support_vs_engine"]:
        L += ["", "Disagreements between the age-support reading and the engine's own band reading: " + "; ".join(f"{x['drug'].title()} {x['band']} ({x['age_support']} vs {x['engine']})" for x in r["age_reading_disagreements_age_support_vs_engine"]) + "."]
    s = r["sensitivity"]
    L += ["", "## 5. Sensitivity", "", f"Baseline proposed set: {', '.join(d.title() for d in r['proposed']) or 'none'}.", "", "| Variant | Proposed set | Jaccard | Baseline retained |", "|---|---|--:|--:|"]
    for c, v in s["leave_one_criterion_out"].items():
        L.append(f"| leave out {c} | {', '.join(d.title() for d in v['proposed']) or 'none'} | {v['jaccard_with_baseline']:.2f} | {len(v['baseline_retained'])} |")
    v = s["literature_only_weights"]
    L.append(f"| literature-only weights (LIT, M) | {', '.join(d.title() for d in v['proposed']) or 'none'} | {v['jaccard_with_baseline']:.2f} | {len(v['baseline_retained'])} |")
    v = s["track1_gene_block_off"]
    L.append(f"| Track 1 non-BUB1B genes off | {', '.join(d.title() for d in v['proposed']) or 'none'} | {v['jaccard_with_baseline']:.2f} | {len(v['baseline_retained'])} |")
    e = s["equal_weights"]
    L += ["", "Equal weights (one vector), top 5: " + ", ".join(f"{x['drug'].title()} ({x['score']:.3f})" for x in e["top5"]) + (f"; ties with the fifth: {', '.join(t.title() for t in e['ties_with_fifth'])}" if e["ties_with_fifth"] else "")
          + f". Overlap with the baseline proposed set (Jaccard) {e['jaccard_top5_with_baseline_proposed']:.2f}.", ""]
    g = s["track1_gene_block_off"]
    L += ["Track 1 gene block (cap " + f"{r['track1_gene_block']['cap']:.2f}" + "): drugs with a bonus: " + (", ".join(f"{d.title()} ({'/'.join(x['genes'])}, T={x['bonus']})" for d, x in g["drugs_with_gene_bonus"].items() if d in ds) or "none")
          + ". Change in P(top 5) when the genes are off: " + (", ".join(f"{d.title()} {dv:+.3f}" for d, dv in sorted(g["p_top5_change"].items(), key=lambda x: -abs(x[1]))[:8]) or "none") + ". "
          + r["track1_gene_block"]["statement"], "",
          "Monte Carlo on label errors (levels or directions of verified literature entries corrupted; pool and exclusions fixed):", "", "| Error rate | Expected baseline candidates retained | Mean set size | Mean Jaccard | Retained fraction by candidate |", "|---|--:|--:|--:|---|"]
    for k, v in s["label_error_monte_carlo"].items():
        L.append(f"| {k} % | {v['expected_baseline_retained']:.2f} of {v['n_baseline']} | {v['mean_set_size']:.2f} | {v['mean_jaccard']:.2f} | " + ", ".join(f"{d.title()} {f:.2f}" for d, f in v["retained_fraction_by_baseline_candidate"].items()) + " |")
    c = r["comparison"]
    L += ["", "## 6. Comparison with the binary rule and with the nine drugs of the report", "", f"The earlier rule ({c['old_rule_text']}) gave candidates: {json.dumps(c['old_rule_candidates'], sort_keys=True)}. The new rule proposes {len(c['new_proposed'])}.", "",
          "| Report drug | Group | Status | Rank by P(top 5) | P(top 5) |", "|---|---|---|--:|--:|"]
    for d, v in c["nine_report_drugs"].items():
        L.append(f"| {d.title()} | {v['group']} | {v['status']} | {v.get('rank_by_p_top5', '-')} | {v.get('p_top5', '-')} |")
    L += ["", "Proposed here but not among the nine: " + (", ".join(d.title() for d in c["proposed_not_in_report"]) or "none") + ".", "",
          "## 7. Data gaps and limits", ""] + [f"- {x}" for x in r["data_gaps"]]
    L += ["- The integration orders hypotheses; it does not create evidence. Pool membership comes from the same lines that are integrated, literature labels were assigned by AI search tasks and not reviewed by a specialist, and weights are uniform over a simplex, not elicited.",
          "- Line grades are positions in lists (A, B1, B2, ENG) or exploratory signals (D, E, EP); the patient signature comes from one adult cell line. A rank-based grade is not an effect size.",
          "- Label support by age describes the US label, not benefit. No clinical recommendation, no dosing, no claim of efficacy.", ""]
    return "\n".join(L)


def pools_summary(r: dict[str, Any]) -> dict[str, Any]:
    ds = r["drugs"]
    ranked = sorted(ds, key=lambda d: (-ds[d]["p_top5"], -ds[d]["p_top3"], d))
    s = r["sensitivity"]
    return {"pool_size": r["pool"]["size"], "proposed": r["proposed"], "primary": r["primary"], "top_by_p_top5": [(d, ds[d]["p_top5"]) for d in ranked[:8]],
            "leave_one_out_proposed": {c: v["proposed"] for c, v in s["leave_one_criterion_out"].items()},
            "mc_retained": {k: (v["expected_baseline_retained"], v["n_baseline"], v["retained_fraction_by_baseline_candidate"]) for k, v in s["label_error_monte_carlo"].items()},
            "equal_weights_top5": [x["drug"] for x in s["equal_weights"]["top5"]], "literature_only_proposed": s["literature_only_weights"]["proposed"],
            "track1_off_proposed": s["track1_gene_block_off"]["proposed"], "age_bands": {b: a["proposed"] for b, a in r["age_bands"].items()},
            "nine": {d: (v["status"], v.get("rank_by_p_top5"), v.get("p_top5")) for d, v in r["comparison"]["nine_report_drugs"].items()}}


def run_completed(out: Path, universe_path: Path, rationale_path: Path, gene_block_path: Path | None, root: Path, n_draws: int = N_DRAWS, mc_reps: int = MC_REPS, mc_draws: int = MC_DRAWS) -> dict[str, Any]:
    """Second run with the three completed gaps: approval check of the previously unverified candidates, the full engine ranking, and drug-level line C."""
    from track2 import multi_evidence_completion as MC

    original = load_json(out / "multi_evidence_mva.json")
    appr = load_json(out / "fda_approval_check_pool_completion.json") or {}
    eng = load_json(out / "engine" / "apply_mva_full_ranking.json")
    universe = load_json(universe_path)
    case = build_case(out, universe, load_json(rationale_path), load_json(gene_block_path) if gene_block_path else None, extra_approved=set(appr.get("confirmed", [])), engine_full=eng)
    linec = MC.line_c_drug_level(root, sorted(case["pool"]))  # data on disk only
    MC.write_json(out / "line_c_drug_level.json", linec)
    case["lines"]["C"] = {d: v["grade"] for d, v in linec["drugs"].items() if v.get("grade") is not None}
    age, diffs = age_table(out)
    result = analyse(case, age, n_draws=n_draws, mc_reps=mc_reps, mc_draws=mc_draws)
    comp = comparison(case, result, out)
    gb = case["gene_block"] or {}
    coverage = {b: sum(1 for d in result["order"] if bands_for(d, age)[b] != "not_assessed") for b in BANDS}
    completed = {
        "label": "hypotheses for follow-up; low evidence; the integration orders hypotheses and creates no evidence; no clinical recommendation, dose or claim of efficacy",
        "protocol": "protocol/addendum_multi_evidence.md (post hoc; section 11 completion)", "seed": SEED, "n_draws": n_draws, "criteria": CRITERIA, "missing_penalty": MISSING_PENALTY,
        "proposal_rule": f"P(top 5) >= {P5_CUT}, not hard-excluded; primary = top {N_PRIMARY}",
        "pool": {"size": len(case["pool"]), "n_after_exclusion": len(result["order"]), "drugs": sorted(case["pool"]), "approval_not_checked": sorted(case["not_checked"]), "approval_checked_not_eligible": sorted(case["not_eligible"]),
                 "n_not_small_molecule_or_not_in_universe": len(case["not_small_or_unknown"]),
                 "n_by_source": {k: len(set(v) & set(case["pool"])) for k, v in case["tops"].items()} | {"literature": len(set(case["entries"]) & set(case["pool"]))}},
        "excluded": result["excluded"], "proposed": result["proposed"], "primary": result["primary"], "drugs": result["drugs"], "age_bands": result["age_bands"],
        "age_label_coverage_n_assessed": coverage, "age_reading_disagreements_age_support_vs_engine": diffs, "sensitivity": result["sensitivity"], "comparison": comp,
        "track1_gene_block": {"cap": case["cap"], "genes": [{k: g[k] for k in ("gene", "role", "tier", "tier_weight")} for g in gb.get("genes", [])], "tests": gb.get("tests", {}),
                              "statement": "These are hypotheses to test (karyotype or single-cell sequencing for the megabase events; function tests for IGBP1), not findings."},
        "consequence_confidence": CONSEQUENCE_CONFIDENCE,
        "completion": {"approval_file": "out/fda_approval_check_pool_completion.json", "n_asked": len(appr.get("drugs", {})), "n_confirmed": appr.get("n_confirmed"), "not_confirmed": appr.get("not_confirmed", []),
                       "engine_file": "out/engine/apply_mva_full_ranking.json", "engine_consistent_with_stored_top20": (eng or {}).get("consistent_with_stored_top20"),
                       "line_c_file": "out/line_c_drug_level.json", "line_c_n_in_prism": linec["n_in_prism"], "line_c_n_tested": linec["n_tested"], "line_c_n_significant": linec["n_significant_q_below_cut"],
                       "line_c_spearman_score_vs_proliferation": linec["spearman_score_vs_proliferation"]},
        "data_gaps": ["28 of the 152 candidates were not confirmed by the approval check (listed in the completion file); the check uses name matching and some may be approved under another name",
                      "line C is correlational, from cancer cell lines, not adjusted for proliferation or lineage, and restricted to pool drugs",
                      "confidences of the five consequences are taken from report section 3.1 (out/differential_diagnosis.json holds the gene-level differential only)",
                      "age-band label support is assessed for the drugs of the label table only; the rest are 'not_assessed'"],
    }
    completed["side_by_side"] = {"original": pools_summary(original), "completed": pools_summary(completed)}
    completed["trametinib_survives"] = "TRAMETINIB" in completed["proposed"]
    completed["new_proposed_vs_original"] = sorted(set(completed["proposed"]) - set(original["proposed"]))
    return completed


def completed_cautions(c: dict[str, Any]) -> list[str]:
    """Plain-language cautions about the proposed drugs of the completed run (generated from the data)."""
    out = []
    for d in c["proposed"]:
        x = c["drugs"][d]
        if x["grades"]["LIT"] is None and x["grades"]["M"] is None:
            out.append(f"- {d.title()} has no verified literature entry and no mechanism fit (LIT and M missing, {x['n_missing']} of {len(CRITERIA)} criteria missing, confidence factor {x['confidence_factor']:.2f}): it is proposed only because of list-position and exploratory lines, "
                       "so it is the weakest kind of proposal, and its stability under label errors (Monte Carlo) is trivially 1 because it has no literature labels to corrupt. It is a hypothesis to check against the label and the literature, not a finding.")
    if "CALCITRIOL" in c["proposed"]:
        out.append("- Calcitriol, flag to verify (not part of the pre-specified safety criterion, which found no hazard flag in the project files): the record of the child documents nephrocalcinosis (report section 3.1) and vitamin D analogues act on calcium handling, so the label must be read for hypercalcaemia and kidney warnings before any follow-up.")
    return out


def to_markdown_completed(c: dict[str, Any], original: dict[str, Any]) -> str:
    o, n = c["side_by_side"]["original"], c["side_by_side"]["completed"]

    def t(xs: list[str]) -> str:
        return ", ".join(x.title() for x in xs) or "none"

    tram = c["drugs"].get("TRAMETINIB", {}).get("p_top5")
    L = ["# Multi-evidence integration: original and completed runs side by side", "", f"*{c['label']}.*", "",
         "Protocol: `protocol/addendum_multi_evidence.md`, section 11 (completion fixed before the rerun). The original outputs are unchanged; the proposal rule, weights, SMAA settings and sensitivity analyses are the same in both runs.", "",
         "## Answer", "",
         f"- Trametinib {'survives: it is still a proposed candidate' if c['trametinib_survives'] else 'does not survive: it is not a proposed candidate in the completed run'} "
         f"(P(top 5) {tram:.3f} in the completed run; {original['drugs']['TRAMETINIB']['p_top5']:.3f} in the original).",
         f"- New proposed drugs in the completed run: {t(c['new_proposed_vs_original'])}.",
         *completed_cautions(c), ""]
    L += [
         "## What was completed", "",
         f"- Approval: the project's FDA check was run for {c['completion']['n_asked']} candidates without a record; {c['completion']['n_confirmed']} were confirmed (approved small molecules among them entered the pool); {len(c['completion']['not_confirmed'])} were not confirmed.",
         f"- Engine: frozen B0 score and rank stored for every drug (consistent with the stored top 20: {c['completion']['engine_consistent_with_stored_top20']}); 'top 50 of the engine' now uses the full ranking.",
         f"- Line C at drug level: {c['completion']['line_c_n_in_prism']} pool drugs in PRISM, {c['completion']['line_c_n_tested']} tested, {c['completion']['line_c_n_significant']} with BH q < 0.10. Correlational, from cancer cell lines, "
         f"not adjusted for proliferation (the score correlates {c['completion']['line_c_spearman_score_vs_proliferation']} with a proliferation score) or lineage.", "",
         "## Side by side", "", "| | Original | Completed |", "|---|---|---|",
         f"| Pool size (after exclusions) | {original['pool']['size']} ({original['pool']['n_after_exclusion']}) | {c['pool']['size']} ({c['pool']['n_after_exclusion']}) |",
         f"| Proposed candidates | {t(o['proposed'])} | {t(n['proposed'])} |", f"| Primary | {t(o['primary'])} | {t(n['primary'])} |",
         "| Best by P(top 5) | " + "; ".join(f"{d.title()} {p:.2f}" for d, p in o["top_by_p_top5"][:5]) + " | " + "; ".join(f"{d.title()} {p:.2f}" for d, p in n["top_by_p_top5"][:5]) + " |",
         f"| Equal weights, top 5 | {t(o['equal_weights_top5'])} | {t(n['equal_weights_top5'])} |",
         f"| Literature-only weights, proposed | {t(o['literature_only_proposed'])} | {t(n['literature_only_proposed'])} |",
         f"| Track 1 non-BUB1B genes off, proposed | {t(o['track1_off_proposed'])} | {t(n['track1_off_proposed'])} |"]
    for k in ("10", "20", "30"):
        a, b = o["mc_retained"][k], n["mc_retained"][k]
        L.append(f"| Label errors {k} %: baseline candidates retained | {a[0]:.2f} of {a[1]} | {b[0]:.2f} of {b[1]} |")
    L += ["", "Leave-one-criterion-out, proposed set:", "", "| Left out | Original | Completed |", "|---|---|---|"]
    for k in CRITERIA:
        L.append(f"| {k} | {t(o['leave_one_out_proposed'][k])} | {t(n['leave_one_out_proposed'][k])} |")
    L += ["", "Age bands (proposed in the full pool and passing the label filter):", "", "| Band | Original | Completed |", "|---|---|---|"]
    for b in o["age_bands"]:
        L.append(f"| {b} | {t(o['age_bands'][b])} | {t(n['age_bands'][b])} |")
    L += ["", "The nine drugs of the report (status, rank by P(top 5), P(top 5)):", "", "| Drug | Original | Completed |", "|---|---|---|"]

    def f(v: list[Any]) -> str:
        return f"{v[0]}" + (f", rank {v[1]}, P {v[2]:.3f}" if v[1] else "")

    for d in o["nine"]:
        L.append(f"| {d.title()} | {f(o['nine'][d])} | {f(n['nine'][d])} |")
    L += ["", "---", "", "# Completed run in detail", ""]
    return "\n".join(L) + "\n" + "\n".join(to_markdown(c).split("\n")[1:])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("out"))
    ap.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    ap.add_argument("--rationale", type=Path, default=Path("config/repositioning_rationale.json"))
    ap.add_argument("--gene-block", type=Path, default=Path("config/multi_evidence_track1_gene_block.json"))
    ap.add_argument("--cache", type=Path, default=Path("data/cache"))
    ap.add_argument("--draws", type=int, default=N_DRAWS)
    ap.add_argument("--reps", type=int, default=MC_REPS)
    ap.add_argument("--completed", action="store_true", help="second run with the completed data gaps; writes *_completed outputs, the first outputs are untouched")
    args = ap.parse_args()
    if args.completed:
        res = run_completed(args.out, args.universe, args.rationale, args.gene_block if args.gene_block.exists() else None, Path("."), n_draws=args.draws, mc_reps=args.reps)
        orig = load_json(args.out / "multi_evidence_mva.json")
        (args.out / "multi_evidence_mva_completed.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        (args.out / "multi_evidence_mva_completed.md").write_text(to_markdown_completed(res, orig), encoding="utf-8", newline="\n")
        print({"pool": res["pool"]["size"], "proposed": res["proposed"], "primary": res["primary"], "trametinib_survives": res["trametinib_survives"], "new": res["new_proposed_vs_original"]})
        return
    res = run(args.out, args.universe, args.rationale, args.gene_block if args.gene_block.exists() else None, args.cache if args.cache.exists() else None, n_draws=args.draws, mc_reps=args.reps)
    (args.out / "multi_evidence_mva.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    (args.out / "multi_evidence_mva.md").write_text(to_markdown(res), encoding="utf-8", newline="\n")
    print({"pool": res["pool"]["size"], "proposed": res["proposed"], "primary": res["primary"]})


if __name__ == "__main__":
    main()
