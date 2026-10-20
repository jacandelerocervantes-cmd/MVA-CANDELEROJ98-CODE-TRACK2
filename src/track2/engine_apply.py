r"""Application of the frozen B0 engine to mosaic variegated aneuploidy (MVA) and Fanconi anaemia.

LABEL: application of the frozen engine; exploratory; not a candidate list.

Run:  PYTHONPATH=src python -m track2.engine_apply [mva|fanconi|both]   (from the repository root)
Writes out/engine/apply_mva.{json,md} and out/engine/apply_fanconi.{json,md} (sorted keys).

Protocol: protocol/addendum_engine_apply.md (written before the run). This is NOT a validation and makes no claim of efficacy.
No existing module is modified: engine_case_study and engine_v4 functions are imported and reused. No label of either disease is
read (neither disease is in the v4 benchmark). No patient data: gene symbols and disease names only.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from track2 import empirical_literature_control as elc
from track2 import engine_benchmark as eb
from track2 import engine_case_study as cs
from track2 import engine_score as es
from track2 import engine_score_v2 as v2
from track2 import engine_v4 as v4

ROOT = Path(__file__).resolve().parents[2]
LABEL = "application of the frozen engine; exploratory; not a candidate list"
TOP_N = 20
MAX_YEAR = es.YEAR_PRIMARY
HIT_STRATUM_L = 50           # engine_v4.stratum_of: L < 50
DISEASES = {
    "mva": {"title": "mosaic variegated aneuploidy (MVA)", "landscape": "out/mva_landscape.json", "orpha": ["ORPHA:1052"],
            "known_cases": "the proband is one known case (stated by the project; not computed from any record)"},
    "fanconi": {"title": "Fanconi anaemia", "landscape": "out/landscape_fanconi_anaemia.json", "orpha": ["ORPHA:84"], "known_cases": None},
}
# landscape mechanism class -> general class of the borrowing pools of validation_v4 (None = no pool)
POOL_OF_CLASS = {"spindle-assembly checkpoint": "cell_cycle_mitotic_checkpoint", "DNA repair/replication": "dna_repair_replication",
                 "DNA interstrand-crosslink repair": "dna_repair_replication"}
NINE = [("sirolimus", ("sirolimus",)), ("everolimus", ("everolimus",)), ("metformin", ("metformin",)),
        ("sodium phenylbutyrate", ("sodium phenylbutyrate", "phenylbutyrate", "phenylbutanoic acid", "sodium phenylbutanoate")),
        ("trametinib", ("trametinib",)), ("selumetinib", ("selumetinib",)), ("bortezomib", ("bortezomib",)),
        ("chloroquine", ("chloroquine",)), ("bosutinib", ("bosutinib",))]
LIMITS = [
    "Application of the frozen B0 engine; exploratory; not a candidate list. It is not a validation and makes no claim of efficacy.",
    "The engine was validated on diseases with known drugs (benchmark gold from approved indications and trials). Where no drug is known "
    "the ranking cannot be checked.",
    "On Wilms tumour the engine failed to rank the known drugs of the disease in its top 20; a high rank is not evidence that a drug works "
    "and a low rank is not evidence against it.",
    "Literature counts are co-mentions among at most 300 abstracts; they are not effect evidence. The knowledge-graph prior is target based.",
    "ChEMBL molecule type and openFDA label presence are used as recorded; Drugs@FDA is a confirmation of marketing status only.",
    "Label age support is label text read by a heuristic classifier unless config/pediatric_label_review.json holds a manual entry "
    "(manual review required otherwise); the age of a patient is never assumed.",
    "Oncology-cytotoxic drugs carry serious label-documented toxicity (marrow suppression, second malignancy, embryo-fetal harm). In "
    "Fanconi anaemia DNA-crosslinking agents are of particular concern because of cellular hypersensitivity. This is a flag, not advice.",
    "No dosing and no treatment recommendation. Hypotheses only.",
]


# ============================================================================= pure helpers (tested offline)
def _r(x: Any, nd: int = 6) -> Optional[float]:
    return None if x is None else es._r(x, nd)


def curated_genes(landscape: dict) -> Tuple[List[dict], List[dict]]:
    """(curated genes, proposed-not-curated genes) of a landscape file, sorted by symbol."""
    g = sorted(landscape["genes"], key=lambda x: x["symbol"])
    cur = [x for x in g if x.get("tier") == "curated"]
    return cur, [x for x in g if x.get("tier") != "curated"]


def genes_by_class(genes: Sequence[dict]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for g in genes:
        out.setdefault(g["mechanism_class"], []).append(g["symbol"])
    return {k: sorted(set(v)) for k, v in sorted(out.items())}


def low_literature(n_papers_with_gene: int, hit_count: Optional[int], threshold: int) -> Dict[str, Any]:
    """Flag fixed in protocol/addendum_engine_apply.md section 5.8: few gene papers (<= threshold) or engine stratum L."""
    few = n_papers_with_gene <= threshold
    stratum_l = hit_count is not None and hit_count < HIT_STRATUM_L
    return {"flag": bool(few or stratum_l), "papers_with_gene": int(n_papers_with_gene), "threshold_papers": int(threshold),
            "europe_pmc_hit_count": hit_count, "hit_count_stratum": None if hit_count is None else v4.stratum_of(hit_count),
            "text": ("low literature: ranking driven by the prior" if (few or stratum_l)
                     else "literature volume above the dilution threshold: the literature features carry part of the ranking")}


def dilution_threshold(dilution: dict) -> Tuple[int, Dict[str, Any]]:
    """Papers at the 10 % thinning (mean over benchmark diseases) from validation_v4_dilution.json, and the gains it rests on."""
    fr = dilution["heldout"]["table_disease_level_bootstrap_sensitivity"]["fractions"]
    thr = int(math.ceil(fr["0.1"]["n_papers_mean"]))
    gains = {f: {"n_papers_mean": _r(fr[f]["n_papers_mean"]), "B0_auc_gap_vs_popularity": _r(fr[f]["B0"]["auc_gap_vs_popularity"]["estimate"])}
             for f in sorted(fr)}
    return thr, gains


def oncology_flag(atc: Sequence[str], in_caution_class: bool) -> bool:
    """Project rule (engine_case_study.evidence_gap): ATC L01; plus members of the caution-list classes."""
    return bool(any(a.startswith("L01") for a in atc) or in_caution_class)


def resolve_nine(names: Sequence[str]) -> Dict[str, Optional[str]]:
    """Report drug -> universe name (lower case) or None."""
    s = set(names)
    out: Dict[str, Optional[str]] = {}
    for label, alts in NINE:
        out[label] = next((a for a in alts if a in s), None)
    return out


def positions(score: np.ndarray, names: Sequence[str], masks: Dict[str, np.ndarray], drugs: Dict[str, Optional[str]]) -> Dict[str, Any]:
    """Rank_min of each drug among the candidates of each mask (None when the drug is not in the mask) and the mask size."""
    ranks = {k: cs.rank_in_mask(score, m, names) for k, m in masks.items()}
    idx = {n: i for i, n in enumerate(names)}
    out = {}
    for label, n in drugs.items():
        row: Dict[str, Any] = {"universe_name": n}
        if n is not None:
            row["score"] = _r(score[idx[n]])
            for k in masks:
                row[f"rank_{k}"] = ranks[k].get(n)
                row[f"of_{k}"] = int(masks[k].sum())
        out[label] = row
    return out


def pool_for(classes: Dict[str, List[str]], borrowing: dict) -> Dict[str, Any]:
    """Borrowing pool of the plurality mappable mechanism class (validation_v4 evidence_gap block); own AUC does not exist."""
    best = None
    for c, g in classes.items():
        p = POOL_OF_CLASS.get(c)
        if p and (best is None or (-len(g), c) < (-len(classes[best[0]]), best[0])):
            best = (c, p)
    if best is None or best[1] not in borrowing["classes"]:
        return {"landscape_class": best[0] if best else None, "pool": None, "note": "no mappable class pool: nothing is borrowed"}
    b = borrowing["classes"][best[1]]
    return {"landscape_class": best[0], "pool": best[1], "pool_k_diseases": b["k"], "pool_mean_auc": _r(b["mu"]),
            "pool_between_disease_variance_tau2": _r(b["tau2"]),
            "note": "class-level mean and variance of the engine's AUC over benchmark diseases of the class; the disease's own AUC does not "
                    "exist (no label), so the weight on it is not computed"}


def evidence_gap_block(key: str, tier: str, tier_codes: List[str], pool: dict, band_gap_counts: Dict[str, int], n_drugs: int,
                       known_cases: Optional[str]) -> Dict[str, Any]:
    """Gaps and the project's designs that could close them. Computed from tier, pool and label states; no dose, no recommendation."""
    gaps = []
    eff_tier = "1 known case" if known_cases else tier
    gaps.append({"gap": "number of patients available for any study",
                 "state": known_cases or f"Orphanet tier {tier}", "orphanet_tier": tier, "orpha_codes": tier_codes,
                 "design_from_borrowing_logic": v4.DESIGN_BY_TIER[eff_tier]})
    if pool.get("pool"):
        gaps.append({"gap": "no label to measure the engine on this disease; only the class-level behaviour is known",
                     "state": f"class pool {pool['pool']}: mean AUC {pool['pool_mean_auc']}, tau2 {pool['pool_between_disease_variance_tau2']}, "
                              f"k={pool['pool_k_diseases']}",
                     "design_from_borrowing_logic": "Bayesian hierarchical model borrowing from the mechanism class pool; any study result of one "
                                                    "disease of the class updates the class level (hypothesis for study design only)."})
    else:
        gaps.append({"gap": "no label and no class pool", "state": pool.get("note"),
                     "design_from_borrowing_logic": "Ascertain patients through registries before choosing a design."})
    for b, n in band_gap_counts.items():
        gaps.append({"gap": f"age band {cs.BAND_TEXT[b]} without labelled support for {n} of the {n_drugs} top-ranked drugs",
                     "state": "states verify, not_labelled or age_not_stated",
                     "design_from_borrowing_logic": "paediatric extrapolation or an age-stratified label supplement is a sponsor/regulatory step; "
                                                    "within a study, the N-of-1 or series design of the first gap applies with an independent safety committee."})
    return {"gaps": gaps, "note": "Hypotheses for study design only: no dose, no drug recommendation."}


def band_gap_counts(supports: Sequence[dict]) -> Dict[str, int]:
    out = {}
    for b in cs.BANDS:
        out[b] = sum(1 for s in supports if cs.status_class(s["bands"][b]) in ("verify", "not_labelled", "age_not_stated"))
    return out


def approval_state(label_match: bool, fda: Optional[dict]) -> str:
    if fda is None:
        return "openFDA label match" if label_match else "not approved"
    el = fda.get("eligible")
    if el is None:
        return ("openFDA label match; " if label_match else "") + "Drugs@FDA check failed (not read as approved)"
    return ("openFDA label match + " if label_match else "") + ("Drugs@FDA Prescription" if el else "Drugs@FDA: no Prescription product")


def orphanet_tier(codes: Sequence[str], root: Path) -> Tuple[str, str]:
    """Case-count tier from the cached Orphanet prevalence file only (no download); (tier, source note)."""
    p = root / "data" / "cache_engine" / "landscape" / (hashlib.sha256(v4.ORPHA_PREV_URL.encode()).hexdigest()[:24] + ".bin")
    if not p.exists():
        return "unknown", "Orphanet prevalence file not in the local cache"
    prev = v4.parse_orphanet_prevalence(p.read_bytes())
    tiers = [v4.tier_from_records(prev[c]) for c in codes if c in prev]
    if not tiers:
        return "unknown", "no prevalence record for " + ", ".join(codes)
    return max(tiers, key=lambda t: v4.TIER_ORDER[t]), "Orphanet prevalence records of " + ", ".join(c for c in codes if c in prev)


# ============================================================================= run
def _no_network(url: str, params: dict) -> dict:
    raise RuntimeError("openFDA label cache miss")


def run(which: str, root: Path = ROOT, getter: Optional[Callable] = None, fetch: Optional[Callable] = None,
        lookup: Optional[Callable] = None, out_dir: Optional[Path] = None) -> Dict[str, Any]:
    from track2 import age_label_table as alt
    from track2 import caution_list as cl
    from track2 import engine_benchmark_v2 as b2mod
    from track2 import fda_approved
    from track2.pediatric_review import load_review

    cfg = DISEASES[which]
    landscape = json.loads((root / cfg["landscape"]).read_text(encoding="utf-8"))
    cur, proposed = curated_genes(landscape)
    symbols = [g["symbol"] for g in cur]
    classes = genes_by_class(cur)
    frozen_path = root / "out" / "engine" / "engine_v4_frozen.json"
    model = v4.thaw_v4(json.loads(frozen_path.read_text(encoding="utf-8")))      # raises if the model hash does not match
    validation = json.loads((root / "out" / "engine" / "validation_v4.json").read_text(encoding="utf-8"))
    dilution = json.loads((root / "out" / "engine" / "validation_v4_dilution.json").read_text(encoding="utf-8"))
    thr, gains = dilution_threshold(dilution)

    drugs = eb.load_universe(root / "data" / "universe.json")
    universe = eb.universe_index(drugs)
    names = sorted(universe)
    index = {n: i for i, n in enumerate(names)}
    dre = elc.drug_regex(names)
    gate = es.name_gate(names)
    guard = v2.RequestGuard(getter or v2.epmc_getter, None)
    cdir = root / "data" / "cache_engine"
    targets, ppi = es.load_prior_graph(root, cdir / "primekg_prior_edges.json")
    pop, tot = es.fetch_popularity(names, None, cdir / "epmc_popularity_all.json", guard)

    # ---- features (no label) and frozen B0 score
    papers = es.fetch_papers(symbols, MAX_YEAR, cdir / "epmc", guard)
    lit = es.literature_counts(papers, symbols, dre, gate)
    prior = es.prior_features(symbols, names, targets, ppi)
    X = es.disease_matrix(lit, names, pop, tot, prior)
    pred = v4.predict_frozen_v4(model, "B0", X)
    sc, pl = pred["score"], pred["platt"]
    hits = v4.hit_count(symbols, guard, {})
    counts = lit["counts"]
    litvol = {"n_papers_retrieved": len(papers), "n_papers_with_gene": lit["n_papers_with_gene"], "europe_pmc_hit_count": hits,
              "n_drugs_with_at_least_one_comention": sum(1 for v in counts.values() if v["n_papers"] >= 1),
              "total_drug_mentions_papers_sum": sum(v["n_papers"] for v in counts.values()),
              "dilution_threshold_papers": thr, "dilution_b0_gain_over_popularity_by_papers": gains,
              "dilution_prior_only_gain_over_popularity_at_5_to_9_papers": None,
              "flag": low_literature(lit["n_papers_with_gene"], hits, thr)}
    abl = json.loads((root / "out" / "engine" / "validation_v4_ablation.json").read_text(encoding="utf-8"))
    po = abl["heldout"]["table_disease_level_bootstrap_sensitivity"]["fractions"]
    litvol["dilution_prior_only_gain_over_popularity_by_papers"] = {
        f: {"n_papers_mean": _r(po[f]["n_papers_mean"]), "PRIOR_ONLY_auc_gap_vs_popularity": _r(po[f]["PRIOR_ONLY"]["auc_gap_vs_popularity"]["estimate"])}
        for f in sorted(po)}
    litvol.pop("dilution_prior_only_gain_over_popularity_at_5_to_9_papers")

    # ---- filter chain (openFDA label cache + Drugs@FDA rescue for the top 100 without a label match)
    label_only = {n: cs.has_fda_label(n, b2mod.fetch_labels(n, getter=_no_network)) for n in names}
    lookup_fn = lookup or fda_approved.lookup
    order_all = sorted(range(len(names)), key=lambda i: (-sc[i], names[i]))
    check = sorted({names[i] for i in order_all[:cs.RESCUE_TOP] if not label_only[names[i]]})
    fda: Dict[str, dict] = {}

    def ask(n: str) -> dict:
        if n not in fda:
            try:
                fda[n] = lookup_fn(universe[n]["name"])
            except Exception as e:  # noqa: BLE001 - a service error is never read as approved
                fda[n] = {"eligible": None, "error": repr(e)[:200]}
        return fda[n]

    for n in check:
        ask(n)
    has_label = {n: label_only[n] or bool((fda.get(n) or {}).get("eligible")) for n in names}
    chain = cs.filter_chain(names, universe, has_label)
    m_all = np.array([has_label[n] for n in names])
    sm = set(chain["small_molecule_names"])
    m_sm = np.array([n in sm for n in names])
    m_univ = np.ones(len(names), dtype=bool)
    chain_out = {k: v for k, v in chain.items() if k not in ("fda_label_names", "small_molecule_names")}
    chain_out.update({"n_fda_label_only_rule": sum(label_only.values()), "n_checked_at_drugs_at_fda": len(check),
                      "n_added_by_drugs_at_fda_check": sum(1 for n in names if has_label[n] and not label_only[n])})

    caution_members = {m["name"].lower() for k in ("a_tubulina", "b_checkpoint_mitotico", "c_quimioterapia_genotoxica", "d_proteasoma")
                       for m in cl.classify_universe(drugs)[k]}

    top = cs.top_table(sc, names, m_sm, TOP_N)
    for r in top:
        i = r.pop("_i")
        n = r["name"]
        c = counts.get(n, {})
        r.update({"platt_value": _r(pl[i]), "molecule_type": universe[n].get("molecule_type"), "atc": universe[n].get("atc") or [],
                  "first_approval": universe[n].get("first_approval"), "n_papers_with_gene": c.get("n_papers", 0), "n_benefit": c.get("n_benefit", 0),
                  "n_human": c.get("n_human", 0), "openfda_label_match": label_only[n],
                  "oncology_cytotoxic_chemotherapy_caution": oncology_flag(universe[n].get("atc") or [], n in caution_members)})
        r["approval"] = approval_state(label_only[n], ask(n))
        r["drugs_at_fda"] = fda[n]

    # ---- the nine drugs of the report (MVA comparison only)
    nine = None
    if which == "mva":
        nine = positions(sc, names, {"universe": m_univ, "fda_label": m_all, "fda_small_molecule": m_sm}, resolve_nine(names))

    # ---- mechanism class x age band
    class_rows: Dict[str, List[dict]] = {}
    class_blocks: Dict[str, Any] = {}
    sm_idx = [i for i in range(len(names)) if m_sm[i]]
    rmin_sm, _ = cs.rank_arrays(sc[sm_idx])
    pos_sm = {names[i]: k for k, i in enumerate(sm_idx)}
    for c, gl in classes.items():
        cnt = cs.class_comention(papers, gl, dre, gate)
        rows = cs.class_top_drugs(cnt, sm)
        for r in rows:
            r["engine_rank_small_molecule_min"] = int(rmin_sm[pos_sm[r["name"]]])
        class_rows[c] = rows
        class_blocks[c] = {"genes": gl, "n_retrieved_papers_mentioning_class_gene": es.literature_counts(papers, gl, dre, gate)["n_papers_with_gene"],
                           "top_drugs_by_comention_small_molecule_fda": rows,
                           "scope_note": "independent papers among the retrieved abstracts that name the drug and a gene of this class"}

    # ---- label support
    need = sorted({r["name"] for r in top} | {r["name"] for rows in class_rows.values() for r in rows} |
                  ({n["universe_name"] for n in nine.values() if n["universe_name"] and n["universe_name"] in sm} if nine else set()))
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
    supports = {n: cs.label_support(universe[n], fetch_fn, review, alt.row_for) for n in need}
    for r in top:
        r["bands"] = supports[r["name"]]["bands"]
        r["label_status"] = supports[r["name"]]["label_status"]
    if nine:
        for v in nine.values():
            n = v["universe_name"]
            v["bands"] = supports[n]["bands"] if n in supports else None
            v["molecule_type"] = universe[n].get("molecule_type") if n else None
            v["fda_label"] = has_label[n] if n else None
            v["in_top20_filtered"] = bool(n and n in {r["name"] for r in top})

    # ---- evidence gap
    tier, tier_note = orphanet_tier(cfg["orpha"], root)
    pool = pool_for(classes, validation["evidence_gap_hierarchical_borrowing"])
    gap = evidence_gap_block(which, tier, cfg["orpha"], pool, band_gap_counts([supports[r["name"]] for r in top]), len(top), cfg["known_cases"])
    gap["orphanet_tier_source"] = tier_note

    result = {
        "label": LABEL, "disease": cfg["title"], "protocol": "protocol/addendum_engine_apply.md", "model": "B0 (frozen; no refit; no label of the disease)",
        "inputs": {"frozen_model_sha256": model_hash(frozen_path), "universe_sha256": cs.sha256_file(root / "data" / "universe.json"),
                   "landscape_file": cfg["landscape"], "landscape_sha256": cs.sha256_file(root / cfg["landscape"]),
                   "genes_scored_curated": symbols, "genes_proposed_not_scored": [g["symbol"] for g in proposed],
                   "mechanism_classes": classes, "n_universe_drugs": len(names), "max_publication_year": MAX_YEAR},
        "ranking_top20_fda_small_molecule": top, "filter_chain": chain_out, "literature_volume": litvol,
        "mechanism_x_age_band": cs.age_table(class_rows, supports), "mechanism_classes": class_blocks,
        "label_support": {n: supports[n] for n in need}, "drugs_at_fda_checked": {n: fda[n] for n in sorted(fda)},
        "nine_report_drugs_positions": nine, "evidence_gap": gap, "limits": LIMITS,
        "age_rule": "bands 0-4, 5-11, 12-17 y from label text; unknown age is never assumed (a patient of unknown age reads the union of the bands)",
    }
    out_dir = out_dir or root / "out" / "engine"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"apply_{which}.json").write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False, default=v2._default) + "\n", encoding="utf-8")
    (out_dir / f"apply_{which}.md").write_text(render_md(result), encoding="utf-8")
    return result


def model_hash(p: Path) -> str:
    return cs.sha256_file(p)


# ============================================================================= markdown
def render_md(r: Dict[str, Any]) -> str:
    q = lambda t, n=420: cs._q(t, n)  # noqa: E731
    L: List[str] = []
    a = L.append
    a(f"# {r['disease']}: application of the frozen engine (B0)\n")
    a(f"**{r['label']}.** Not a validation, no claim of efficacy, no dosing, no treatment recommendation. Protocol: `{r['protocol']}`.\n")
    lv = r["literature_volume"]
    a(f"**Literature flag: {lv['flag']['text']}** ({lv['n_papers_with_gene']} retrieved papers carry a gene symbol, threshold <= {lv['dilution_threshold_papers']}; "
      f"Europe PMC hit count {lv['europe_pmc_hit_count']}, stratum {lv['flag']['hit_count_stratum']}).\n")
    inp = r["inputs"]
    a("## 1. Inputs\n")
    a(f"* Genes scored (tier curated): {', '.join(inp['genes_scored_curated'])}")
    a(f"* Genes proposed, not curated, not scored: {', '.join(inp['genes_proposed_not_scored']) or 'none'}")
    a("* Mechanism classes: " + "; ".join(f"{c}: {', '.join(g)}" for c, g in inp["mechanism_classes"].items()))
    a(f"* Frozen model sha256 `{inp['frozen_model_sha256'][:16]}`, universe {inp['n_universe_drugs']} drugs, literature <= {inp['max_publication_year']}\n")
    fc = r["filter_chain"]
    a("## 2. Filter chain\n")
    a("| step | drugs |\n|---|---|")
    a(f"| universe | {fc['n_universe']} |")
    a(f"| with FDA label (openFDA match {fc['n_fda_label_only_rule']}, +{fc['n_added_by_drugs_at_fda_check']} by Drugs@FDA check of {fc['n_checked_at_drugs_at_fda']} asked) | {fc['n_fda_label']} |")
    a(f"| small molecule (ChEMBL type) | {fc['n_small_molecule']} |")
    a(f"\nDropped by type after the label step: {json.dumps(fc['dropped_by_type_after_fda'], sort_keys=True)}\n")
    a("## 3. Top 20 (FDA-approved small molecules), frozen B0 score\n")
    a("| # | drug | score | Platt | approval | chemo caution | 0-4 y | 5-11 y | 12-17 y | gene papers (benefit/human) |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for t in r["ranking_top20_fda_small_molecule"]:
        b = t["bands"]
        a(f"| {t['position']} | {t['name']} | {t['score']} | {t['platt_value']} | {t['approval']} | {'YES' if t['oncology_cytotoxic_chemotherapy_caution'] else 'no'} | "
          f"{b['0_4']} | {b['5_11']} | {b['12_17']} | {t['n_papers_with_gene']} ({t['n_benefit']}/{t['n_human']}) |")
    a("\nTies share a score; `rank_min` is in the JSON. 'verify' = needs manual verification; 'not_labelled' = label does not establish use in the band; "
      "'age_not_stated' = established use without a minimum age. Unknown age is never assumed.\n")
    a("## 4. Label quotes (verbatim, truncated where marked)\n")
    for t in r["ranking_top20_fda_small_molecule"]:
        s = r["label_support"][t["name"]]
        a(f"### {t['position']}. {t['name']} (status {s['label_status']}; {s['basis']}{'; manual review required' if s['manual_review_required'] else ''})")
        pu = s["verbatim_pediatric_use"]
        if pu:
            a(f"* Paediatric use ({pu['fuente']}, set id {pu['setid']}, section {pu['seccion']}): \"{q(pu['texto_literal'])}\"")
        else:
            a("* Paediatric use: no readable paediatric-use sentence (not labelled or needs manual verification)")
        for w in (s["verbatim_pediatric_warnings"] or [])[:2]:
            a("* Paediatric warning: " + q(w.get("texto_literal") if isinstance(w, dict) else str(w)))
        if t["oncology_cytotoxic_chemotherapy_caution"]:
            hz = s["verbatim_hazard_quotes"]
            a("* Chemotherapy caution (label-documented hazards): " + ("" if hz else "no hazard sentence retrieved; not read as safe"))
            for h in hz[:4]:
                a(f"  * [{h['category']}, {h['section']}] {q(h['quote'], 320)}")
        a("")
    a("## 5. Mechanism class x age band\n")
    a("Per class: the five small-molecule FDA drugs with most papers co-mentioning the drug and a gene of the class.\n")
    for c, blk in r["mechanism_classes"].items():
        a(f"### {c} (genes {', '.join(blk['genes'])}; {blk['n_retrieved_papers_mentioning_class_gene']} retrieved papers mention a class gene)\n")
        a("| drug | papers | benefit | human | engine rank (small molecule) | 0-4 y | 5-11 y | 12-17 y |\n|---|---|---|---|---|---|---|---|")
        for d in blk["top_drugs_by_comention_small_molecule_fda"]:
            b = r["label_support"][d["name"]]["bands"]
            a(f"| {d['name']} | {d['n_papers']} | {d['n_benefit']} | {d['n_human']} | {d['engine_rank_small_molecule_min']} | {b['0_4']} | {b['5_11']} | {b['12_17']} |")
        if not blk["top_drugs_by_comention_small_molecule_fda"]:
            a("| (no small-molecule FDA drug co-mentioned with a class gene) | | | | | | | |")
        a("\nBand cells by state (drug names):\n")
        for bnd, cell in r["mechanism_x_age_band"][c].items():
            a(f"* {cs.BAND_TEXT[bnd]}: " + ("; ".join(f"{k}: {', '.join(v)}" for k, v in cell.items()) or "none"))
        a("")
    if r.get("nine_report_drugs_positions"):
        a("## 6. The nine drugs of the report in this ranking (comparison only, not a validation)\n")
        a("| drug | universe name | score | rank in universe | rank among FDA-label | rank among FDA small molecules | in top 20 | 0-4 y | 5-11 y | 12-17 y |\n|---|---|---|---|---|---|---|---|---|---|")
        for k, v in r["nine_report_drugs_positions"].items():
            if v["universe_name"] is None:
                a(f"| {k} | not in universe | | | | | | | | |")
                continue
            b = v.get("bands") or {}
            a(f"| {k} | {v['universe_name']} | {v['score']} | {v['rank_universe']} of {v['of_universe']} | {v['rank_fda_label']} of {v['of_fda_label']} | "
              f"{v['rank_fda_small_molecule'] if v['rank_fda_small_molecule'] is not None else 'filtered out'} of {v['of_fda_small_molecule']} | "
              f"{'yes' if v['in_top20_filtered'] else 'no'} | {b.get('0_4', 'not read')} | {b.get('5_11', 'not read')} | {b.get('12_17', 'not read')} |")
        a("")
    a("## 7. Evidence gap (computed; no dosing; no recommendation)\n")
    g = r["evidence_gap"]
    a(f"Orphanet case-count tier: {g['gaps'][0]['orphanet_tier']} ({g['orphanet_tier_source']}).\n")
    for x in g["gaps"]:
        a(f"* **{x['gap']}** ({x['state']}). Design: {x['design_from_borrowing_logic']}")
    a("\n## 8. Literature volume\n")
    a(f"* Papers retrieved: {lv['n_papers_retrieved']}; with a gene symbol: {lv['n_papers_with_gene']}; Europe PMC hit count: {lv['europe_pmc_hit_count']}")
    a(f"* Drugs co-mentioned at least once: {lv['n_drugs_with_at_least_one_comention']}; sum of paper counts over drugs: {lv['total_drug_mentions_papers_sum']}")
    a(f"* **{lv['flag']['text']}**")
    a("* Dilution result (benchmark, held-out): B0 gain over popularity by mean papers per disease: " +
      "; ".join(f"{v['n_papers_mean']:.0f} papers: {v['B0_auc_gap_vs_popularity']}" for v in lv["dilution_b0_gain_over_popularity_by_papers"].values()) +
      ". Prior-only gain (flat over volume): " +
      "; ".join(f"{v['n_papers_mean']:.0f}: {v['PRIOR_ONLY_auc_gap_vs_popularity']}" for v in lv["dilution_prior_only_gain_over_popularity_by_papers"].values()) + "\n")
    a("## 9. Limits\n")
    for x in r["limits"]:
        a(f"* {x}")
    return "\n".join(L) + "\n"


def main() -> None:
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    for k in (["mva", "fanconi"] if which == "both" else [which]):
        r = run(k)
        print(json.dumps({"disease": k, "top10": [(t["name"], t["score"], t["approval"], t["bands"]) for t in r["ranking_top20_fda_small_molecule"][:10]],
                          "flag": r["literature_volume"]["flag"]["flag"]}, default=str))


if __name__ == "__main__":
    main()
