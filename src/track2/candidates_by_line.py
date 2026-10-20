"""Candidates by line: each line proposes its own hypotheses (protocol/addendum_multi_evidence.md, section 12).

    PYTHONPATH=src python -m track2.candidates_by_line --out out

Each line (literature, engine, A, B1, B2, C, D, E, E-P, E-Track1) contributes the drugs that meet ITS OWN native criterion after common filters (approved small
molecule of the completed pool; hard safety exclusions including label warnings about conditions of the child; cytotoxic agents listed as context).
Convergence between lines is information, never a requirement. The multi-evidence SMAA probability is a cross-check column only. Hypotheses for follow-up;
no claim of efficacy, no recommendation, no dose. Gene symbols, drug names and categorical outputs only.
"""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Any, Callable

from scipy.stats import hypergeom

from track2 import multi_evidence as ME

LINES = ["LIT", "ENG", "A", "B1", "B2", "C", "D", "E", "EP", "ET1"]
LINE_NAMES = {"LIT": "Literature (research line)", "ENG": "Engine (frozen B0)", "A": "A knowledge graph", "B1": "B1 phenotype module", "B2": "B2 network proximity (mitosis module)",
              "C": "C PRISM drug-level vulnerability", "D": "D signed direction", "E": "E signature reversal, BUB1B knockout signature", "EP": "E-P signature reversal, patient signature (GSE22206)",
              "ET1": "E-T1 signature reversal, Track 1 combined signature"}
LIMIT = {"LIT": "AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read",
         "ENG": "prior from literature co-mention and graph targets; cannot be checked where no drug is known",
         "A": "weak discrimination (excess over a degree baseline 0.022); generic-drug filter applied",
         "B1": "did not beat a degree baseline; validation not met", "B2": "exploratory; support set of a few drugs, from the pathogenic gene only",
         "C": "correlational, from cancer cell lines, not adjusted for proliferation or lineage",
         "D": "exploratory; 4 of 6 controls recovered, two by construction",
         "E": "knockout signature may not represent partial loss of function; consensus library, no cell-line-specific signatures",
         "EP": "one adult patient lymphoblastoid cell line, 2 arrays against 2 controls", "ET1": "weighted combination of gene signatures; consensus library"}
EXPERIMENT = {"LIT": "dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress)",
              "ENG": "test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts",
              "A": "dose-response in BUB1B-deficient versus isogenic control cells for the phenotype the graph links to the drug",
              "B1": "target-engagement and rescue of the phenotype module in patient-derived cells",
              "B2": "target engagement and viability or mitotic-fidelity readout in cells with BUB1B loss",
              "C": "isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated",
              "D": "spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells",
              "E": "RNA-seq before and after the drug in BUB1B-depleted cells, testing whether the knockout signature reverses",
              "EP": "RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses",
              "ET1": "RNA-seq before and after the drug in patient-derived cells for the combined gene signature"}
# keyword map by documented term of the child (addendum 12.1, corrected in 13.1): broad keywords only FLAG a label for manual verification
PHENOTYPE_KEYWORDS = {
    "HP:0000121 nephrocalcinosis": ["nephrocalcinosis", "hypercalcemia", "hypercalcaemia", "hypercalciuria", "renal calcification"],
    "HP:0004322/HP:0001508/HP:0001518 growth": ["growth retardation", "growth suppression", "growth inhibition", "decreased growth", "reduced growth", "impaired growth", "stunted growth", "failure to thrive"],
    "HP:0003202 muscle atrophy": ["muscle atrophy", "muscle wasting", "myopathy", "rhabdomyolysis"],
    "HP:0002859 rhabdomyosarcoma": ["rhabdomyosarcoma"]}
CONTRA_ONLY = {"HP:0002859 rhabdomyosarcoma": ["active malignancy", "malignant neoplasm"]}
# documented child-specific contraindications = a SPECIFIC keyword in the contraindications or warnings section (13.1(c)); everything else is a flag
EXCLUSION_TERMS = ("HP:0000121 nephrocalcinosis", "HP:0002859 rhabdomyosarcoma")
EXCLUSION_SECTIONS = ("contraindications", "warnings")
LABEL_NAME = {"PHENYLBUTANOIC ACID": "sodium phenylbutyrate"}
SEARCHED = ("boxed_warning", "contraindications", "warnings")
PAIR_NOTE = "hypergeometric test over the completed pool; no multiplicity claim"


# ---------------------------------------------------------------- label keywords

def sentence_around(text: str, i: int) -> str:
    start = max(text.rfind(". ", 0, i) + 2, 0) if text.rfind(". ", 0, i) >= 0 else 0
    end = text.find(". ", i)
    end = len(text) if end < 0 else end + 1
    return text[start:end].strip()[:400]


def label_flags(sections: dict[str, str]) -> list[dict[str, Any]]:
    """Keyword matches (whole words) in boxed warning, contraindications and warnings; each with its sentence. `exclusion` is true only for the specific
    keywords of the child-specific terms in the contraindications or warnings section (13.1(c)); every other match is a flag for manual verification."""
    out = []
    for term, kws in PHENOTYPE_KEYWORDS.items():
        for sec in SEARCHED:
            text = sections.get(sec) or ""
            low = text.lower()
            for kw in kws + (CONTRA_ONLY.get(term, []) if sec == "contraindications" else []):
                m = re.search(r"\b" + re.escape(kw) + r"\b", low)  # whole words: 'cardiomyopathy' is not 'myopathy'
                if m:
                    out.append({"term": term, "keyword": kw, "section": sec, "sentence": sentence_around(text, m.start()),
                                "exclusion": term in EXCLUSION_TERMS and sec in EXCLUSION_SECTIONS})
                    break
    return out


def default_label_fn() -> Callable[[str], tuple[dict[str, Any] | None, list[str]]]:
    from track2 import caution_list as cl

    def fn(name: str) -> tuple[dict[str, Any] | None, list[str]]:
        lab, failures = cl.fetch_label(name)
        time.sleep(0.3)  # at most about 3 requests per second to DailyMed / openFDA; answers are cached by the project
        return lab, failures
    return fn


# ---------------------------------------------------------------- native sets of the lines: ordered lists of {drug, score ...}

def native_literature(entries: dict[str, Any]) -> list[dict[str, Any]]:
    from track2.empirical_label_robustness import decide
    proposed, _ = decide(entries, set(entries))
    rows = []
    for d, grp in proposed.items():
        g, info = ME.lit_grade(entries[d]["entries"])
        rows.append({"drug": d, "rule_group": grp, "best_help_level": info["best_help_level"], "n_help_papers": info["n_help_papers"], "net": info["net"],
                     "key": (-(info["best_help_level"] if info["best_help_level"] is not None else -1), -info["n_help_papers"], -info["net"], d)})
    rows.sort(key=lambda r: r["key"])
    for i, r in enumerate(rows):
        r["rank_in_line"] = i + 1
        r["native"] = f"rule group {r['rule_group']}; best help level {r['best_help_level']}; {r['n_help_papers']} help papers; net {r['net']}"
        del r["key"]
    return rows


def native_engine(full: dict[str, Any] | None) -> list[dict[str, Any]]:
    if not full:
        return []
    rows = [(v["rank_fda_small_molecule"], d, v["score"]) for d, v in full["drugs"].items() if v.get("rank_fda_small_molecule")]
    rows.sort()
    return [{"drug": d, "rank_in_line": r, "engine_score": s, "native": f"engine rank {r} of {full.get('n_fda_small_molecule')}, score {s:.3f}"} for r, d, s in rows]


def native_A(out: Path, names: dict[str, str]) -> list[dict[str, Any]]:
    q, nq = ME.load_json(out / "kg_query_e5.json"), ME.load_json(out / "kg_nullquery_e5.json")
    if not q:
        return []
    best: dict[str, int] = {}
    for did, r in q["rank_combined"].items():
        if did in names:
            best[names[did]] = min(r, best.get(names[did], r))
    frac = {ME.up(k): v for k, v in (nq or {}).get("fraction_of_random_queries_with_drug_in_top50", {}).items()}
    top = sorted(best, key=lambda d: (best[d], d))[:ME.TOP_N]
    sel = [d for d in top if d in frac and frac[d] < 0.05]
    return [{"drug": d, "rank_in_line": i + 1, "graph_rank": best[d], "random_query_fraction": frac[d], "native": f"graph rank {best[d]}; fraction of random queries with the drug in their top 50 {frac[d]:.3f}"} for i, d in enumerate(sel)]


def native_B1(out: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    p = ME.load_json(out / "phenotype_proximity_e5.json")
    if not p:
        return [], []
    rule = p.get("rule", {})
    rows = sorted(p["drugs"], key=lambda d: (-d["robustness_count"], d["name"]))
    mk = lambda i, d: {"drug": ME.up(d["name"]), "rank_in_line": i + 1, "robustness_count": d["robustness_count"], "native": f"robust in {d['robustness_count']} of {rule.get('configurations', 9)} configurations (needs {rule.get('robust_min', 7)})"}  # noqa: E731
    robust = [mk(i, d) for i, d in enumerate([d for d in rows if d["robust"]])]
    near = [mk(i, d) for i, d in enumerate([d for d in rows if not d["robust"]][:5])]
    return robust, near


def native_B2(out: Path) -> list[dict[str, Any]]:
    s = ME.load_json(out / "track1_scenarios.json")
    if not s:
        return []
    w3 = s["variants"]["two_record_links_no_predicted"]["scheme_robust"]["W3_pathogenic_gene_only"]
    return [{"drug": ME.up(n), "rank_in_line": i + 1, "native": "robust in scheme W3 (pathogenic gene only), two-record links"} for i, n in enumerate(sorted(w3))]


def native_C(out: Path) -> list[dict[str, Any]]:
    c = ME.load_json(out / "line_c_drug_level.json")
    if not c:
        return []
    rows = [(v["q_bh"], -v["rho"], d, v) for d, v in c["drugs"].items() if v.get("q_bh") is not None and v["q_bh"] < 0.10 and v["rho"] > 0]
    rows.sort()
    return [{"drug": d, "rank_in_line": i + 1, "rho": v["rho"], "q_bh": v["q_bh"], "n_lines": v["n_lines"], "native": f"rho {v['rho']:.3f}, BH q {v['q_bh']:.3f}, {v['n_lines']} cell lines"} for i, (_, _, d, v) in enumerate(rows)]


def native_D(out: Path) -> list[dict[str, Any]]:
    s = ME.load_json(out / "signed_effect.json")
    if not s:
        return []
    rows = sorted(((-(r.get("net") or 0), ME.up(r["name"])) for r in s["drugs"].values() if r["call"] == "strengthens_checkpoint"))
    return [{"drug": d, "rank_in_line": i + 1, "net": -n, "native": f"strengthens checkpoint, net signed effect {-n:g}"} for i, (n, d) in enumerate(rows)]


def native_E(out: Path, key: str) -> list[dict[str, Any]]:
    e = ME.load_json(out / "line_e_signature_reversal.json")
    if not e:
        return []
    rows = []
    for r in e["table_all_drugs_n150"]:
        s = r.get(key)
        if s and s.get("usable") and not s.get("hazard"):
            rows.append((s["p_rev"], ME.up(r["name"]), s))
    rows.sort(key=lambda x: (x[0], x[1]))
    return [{"drug": d, "rank_in_line": i + 1, "p_rev": s["p_rev"], "q_rev": s["q_rev"], "R": s["R"], "native": f"reversal p {s['p_rev']:.4f}, FDR q {s['q_rev']:.3f}, R {s['R']}"} for i, (_, d, s) in enumerate(rows)]


# ---------------------------------------------------------------- filters

def classify_drug(d: str, ctx: dict[str, Any]) -> tuple[str, str]:
    """('candidate', '') or (status, reason). Order: eligibility, hard exclusions, cytotoxic context, label keywords."""
    by_name = ctx["by_name"]
    if d not in ctx["pool"]:
        rec = by_name.get(d)
        if rec is None:
            return "not_eligible", "not in the drug universe"
        if rec.get("molecule_type") != "Small molecule":
            return "not_eligible", f"not a small molecule ({rec.get('molecule_type')})"
        if d in ctx["unconfirmed"]:
            return "approval_unconfirmed", "the approval check could not confirm it (may be approved under another name)"
        return "not_eligible", "not in the completed pool (approval not confirmed in the project check)"
    if d in ctx["excluded"]:
        return "hard_excluded", ctx["excluded"][d]
    if d in ctx["caution_names"]:
        return "context_hazard", "cytotoxic or caution-list class: not proposed for constitutional use (project caution)"
    if d in ctx.get("level6", set()):
        return "hard_excluded", "project label level 6 (not indicated or contraindicated in children, src/track2/label_levels.py)"
    if d in ctx["safety_entries"]:
        return "candidate", ""  # covered by the hand-read layer config/proposal_safety.json: its label warnings are shown; no keyword screen
    lab = ctx["labels"](d)
    hits = [f for f in lab["flags"] if f.get("exclusion")]
    if hits:
        f = hits[0]
        return "hard_excluded", f"documented child-specific contraindication tied to {f['term']} ('{f['keyword']}' in the {f['section']} section), for human review. Quote: \"{f['sentence']}\""
    return "candidate", ""


def make_labels_fn(label_fn: Callable[[str], tuple[dict[str, Any] | None, list[str]]]) -> Callable[[str], dict[str, Any]]:
    cache: dict[str, dict[str, Any]] = {}

    def fn(d: str) -> dict[str, Any]:
        if d not in cache:
            lab, failures = label_fn(d)
            if lab is None and d in LABEL_NAME:  # the universe name is the parent acid; the marketed label is the sodium salt
                lab, failures = label_fn(LABEL_NAME[d])
            if lab is None:
                cache[d] = {"read": False, "flags": [], "failures": failures[:2], "source": None}
            else:
                cache[d] = {"read": True, "flags": label_flags(lab["sections"]), "failures": [], "source": f"{lab.get('source')} {lab.get('setid')}"}
        return cache[d]
    return fn


def tier_for(line: str, d: str, rows: dict[str, dict[str, Any]], ctx: dict[str, Any]) -> str:
    r = rows.get(d)
    if r is None or r["grades"]["LIT"] is None:
        head = "no verified literature entry"
    else:
        lit = r["lit"]
        head = (f"level {lit['best_help_level']} literature ({lit['n_help_papers']} help, {lit['n_harm_papers']} harm papers)" if lit["best_help_level"] is not None
                else f"literature without a supportive finding ({lit['n_harm_papers']} harm papers)")
    contrary = []
    if r:
        g = r["grades"]
        contrary = [c for c in ("LIT", "D", "E", "EP") if g[c] is not None and g[c] < 0.5 and (c != "LIT" or r["lit"]["n_harm_papers"] > 0)]
    return f"{head}; contrary: {', '.join(contrary) if contrary else 'none found'}; line limit: {LIMIT[line]}"


def make_age_fn(universe_by_name: dict[str, dict[str, Any]], label_fn: Callable[[str], tuple[dict[str, Any] | None, list[str]]], review: dict[str, Any]) -> Callable[[str], dict[str, Any]]:
    """Age-band reading with the project's functions: engine_case_study.label_support -> age_label_table.row_for + manual review + labels.classify_pediatric."""
    from track2 import age_label_table as alt
    from track2 import engine_case_study as cs

    setids: dict[str, str] = {}
    for e in review.values():
        m = re.search(r"setid=([0-9a-f-]{36})", e.get("source") or "")
        if m:
            setids[e["name"].upper()] = m.group(1)
    cache: dict[str, dict[str, Any]] = {}

    def fetch(name: str) -> tuple[dict[str, Any] | None, list[str]]:
        if name.upper() in setids:
            from track2 import caution_list as cl
            lab, failures = cl.fetch_by_setid(setids[name.upper()])
            time.sleep(0.3)
            if lab:
                return lab, failures
        lab, failures = label_fn(name)
        if lab is None and name in LABEL_NAME:
            lab, failures = label_fn(LABEL_NAME[name])
        return lab, failures

    def fn(d: str) -> dict[str, Any]:
        if d not in cache:
            rec = universe_by_name[d]
            res = cs.label_support(rec, fetch, review, alt.row_for)
            raw = {"0-4": res["bands"].get("0_4"), "5-11": res["bands"].get("5_11"), "12-17": res["bands"].get("12_17")}
            retrievable = res.get("label") is not None
            bands = {b: (ME.engine_band(v) if retrievable else "not_assessed (label not retrievable)") for b, v in raw.items()}
            q = res.get("verbatim_pediatric_use") or {}
            cache[d] = {"bands": bands, "bands_raw": raw, "label_status": res["label_status"], "minimum_age_years": res.get("minimum_age_years"), "basis": res["basis"],
                        "manual_review": res.get("scope_of_manual_review"), "retrievable": retrievable, "label": res.get("label"),
                        "quote_section": q.get("seccion"), "quote": (q.get("texto_literal") or "")[:500] if isinstance(q, dict) else None,
                        "failures": res.get("reading_failures")}
        return cache[d]
    return fn


def evaluate_line(line: str, native: list[dict[str, Any]], ctx: dict[str, Any], rows: dict[str, dict[str, Any]], age: dict[str, dict[str, str]], limit_candidates: int = 5, qualify_cap: int | None = None) -> dict[str, Any]:
    cands, others = [], []
    cap = qualify_cap if qualify_cap is not None else 10 ** 6
    for n in native:
        if len(cands) >= cap:
            break
        status, why = classify_drug(n["drug"], ctx)
        if status == "candidate":
            cands.append(n)
        else:
            others.append({"drug": n["drug"], "status": status, "reason": why, "native": n["native"]})
    out_c = []
    for n in cands:
        d = n["drug"]
        r = rows.get(d)
        rec = ctx["by_name"].get(d)
        lab = ctx["labels"](d)
        sm = ctx["smaa"].get(d)
        flags = []
        if (rec or {}).get("black_box_warning"):
            flags.append("boxed warning in the label")
        if d in ctx["safety_entries"]:
            ent = ctx.get("safety_text", {}).get(d, [])
            flags.append("project hand-read label warnings (config/proposal_safety.json)" + ("".join(f"; [{e}]" for e in ent[:2]) if ent else ""))
        else:
            kw = [f for f in lab["flags"] if not f.get("exclusion")]
            for f in kw[:3]:
                flags.append(f"needs manual verification: {f['term']} '{f['keyword']}' in the {f['section']} section: \"{f['sentence'][:260]}\"")
            if not lab["read"]:
                flags.append("label not read")
            elif not kw:
                flags.append("label read, no keyword match for the phenotype terms")
        reading = ctx["age_fn"](d) if ctx.get("age_fn") else None
        if reading is not None:
            bands = reading["bands"]
        else:
            bands = ME.bands_for(d, age)
        out_c.append({**{k: v for k, v in n.items()}, "line": line, "evidence_tier": tier_for(line, d, rows, ctx), "age_bands": bands,
                      "age_union": "supported in " + ", ".join(b for b in ME.BANDS if bands[b] in ME.BAND_PASS) if any(bands[b] in ME.BAND_PASS for b in ME.BANDS) else "no band with label support (" + "/".join(bands[b] for b in ME.BANDS) + ")",
                      "age_reading": reading,
                      "safety_flags": flags, "label_source": lab["source"],
                      "consequences_addressed": (r or {}).get("consequences", {}), "track1_genes_bearing": (r or {}).get("track1_genes", []),
                      "refuting_experiment": EXPERIMENT[line], "smaa_p_top5_cross_check": None if sm is None else sm["p_top5"], "smaa_rank_by_p_top5": None if sm is None else sm["rank"]})
    return {"line": line, "name": LINE_NAMES[line], "n_native_set": len(native), "n_qualifying": len(out_c), "shown": out_c[:limit_candidates], "also_qualifying": [c["drug"] for c in out_c[limit_candidates:]],
            "found_but_not_proposed": others, "limit": LIMIT[line]}


# ---------------------------------------------------------------- convergence

def convergence(sets: dict[str, set[str]], universe_n: int) -> dict[str, Any]:
    drugs = sorted(set().union(*sets.values())) if sets else []
    matrix = {d: {ln: d in sets[ln] for ln in sets} for d in drugs}
    pairs = []
    lines = [ln for ln in sets if sets[ln]]
    for i, a in enumerate(lines):
        for b in lines[i + 1:]:
            k = len(sets[a] & sets[b])
            p = float(hypergeom.sf(k - 1, universe_n, len(sets[a]), len(sets[b]))) if k > 0 else 1.0
            pairs.append({"lines": [a, b], "sizes": [len(sets[a]), len(sets[b])], "overlap": k, "shared": sorted(sets[a] & sets[b]), "p_hypergeometric": p, "note": PAIR_NOTE})
    n_lines = {d: sum(matrix[d].values()) for d in drugs}
    return {"matrix": matrix, "pairs": pairs, "universe_n": universe_n, "n_candidates_distinct": len(drugs), "n_unique_to_one_line": sum(1 for d in drugs if n_lines[d] == 1),
            "in_two_or_more_lines": sorted(d for d in drugs if n_lines[d] >= 2)}


# ---------------------------------------------------------------- assembly

def changes_vs_first(first: dict[str, Any], results: dict[str, Any], ctx: dict[str, Any] | None = None) -> dict[str, Any]:
    """Per line: candidates added and removed relative to the first output, and age readings that changed."""
    out: dict[str, Any] = {}
    for ln in LINES:
        old = {c["drug"]: c for c in first["lines"][ln]["shown"]}
        new = {c["drug"]: c for c in results[ln]["shown"]}
        changed = {d: {"first": [old[d]["age_bands"][b] for b in ME.BANDS], "v2": [new[d]["age_bands"][b] for b in ME.BANDS]} for d in new if d in old and old[d]["age_bands"] != new[d]["age_bands"]}
        out[ln] = {"first": list(old), "v2": list(new), "added": [d for d in new if d not in old], "removed": [d for d in old if d not in new], "age_reading_changed": changed}
    status = {}
    for d in ("SIROLIMUS", "SELUMETINIB", "CALCITRIOL"):
        st, why = classify_drug(d, ctx) if ctx else ("unknown", "")
        status[d] = {"screen_result": st, "reason": why, "candidate_in_lines": [ln for ln in LINES for c in results[ln]["shown"] if c["drug"] == d],
                     "also_qualifying_in_lines": [ln for ln in LINES if d in results[ln]["also_qualifying"]]}
    out["named_drugs"] = status
    return out


def run(out: Path, universe_path: Path, rationale_path: Path, gene_block_path: Path | None, label_fn: Callable[[str], tuple[dict[str, Any] | None, list[str]]] | None = None,
        completed_name: str = "multi_evidence_mva_completed.json") -> dict[str, Any]:
    comp = ME.load_json(out / completed_name) or {}
    appr = ME.load_json(out / "fda_approval_check_pool_completion.json") or {}
    eng_full = ME.load_json(out / "engine" / "apply_mva_full_ranking.json")
    universe = ME.load_json(universe_path)
    case = ME.build_case(out, universe, ME.load_json(rationale_path), ME.load_json(gene_block_path) if gene_block_path else None, extra_approved=set(appr.get("confirmed", [])), engine_full=eng_full)
    rows = ME.drug_rows(case)
    age, _ = ME.age_table(out)
    pool = set(comp.get("pool", {}).get("drugs", case["pool"]))
    ps = case["proposal_safety"]
    drugs = comp.get("drugs", {})
    ranked = sorted(drugs, key=lambda d: (-drugs[d]["p_top5"], d))
    smaa = {d: {"p_top5": drugs[d]["p_top5"], "rank": ranked.index(d) + 1} for d in drugs}
    ctx = {"by_name": case["by_name"], "pool": pool, "unconfirmed": set(comp.get("completion", {}).get("not_confirmed", [])), "excluded": dict(comp.get("excluded", {})),
           "caution_names": case["caution_names"], "labels": make_labels_fn(label_fn or default_label_fn()), "smaa": smaa, "safety_entries": {ME.up(k) for k in (ps or {}).get("entries", {})}}
    from track2.label_levels import LEVEL
    from track2.pediatric_review import load_review
    ctx["level6"] = {n.upper() for n, lv in LEVEL.items() if lv == 6}
    ctx["safety_text"] = {ME.up(k): [f"{e.get('section')}: {e.get('flag')} - \"{(e.get('quote') or '')[:160]}\"" for e in v] for k, v in (ps or {}).get("entries", {}).items()}
    ctx["age_fn"] = make_age_fn(case["by_name"], label_fn or default_label_fn(), load_review(Path("config/pediatric_label_review.json")))
    for d in case["hard_class"]:
        ctx["excluded"].setdefault(d, "growth-hormone/IGF-1 class: label contraindicates active malignancy")
    names = case["names"]
    b1_robust, b1_near = native_B1(out)
    natives = {"LIT": native_literature(case["entries"]), "ENG": native_engine(eng_full), "A": native_A(out, names), "B1": b1_robust, "B2": native_B2(out), "C": native_C(out),
              "D": native_D(out), "E": native_E(out, "bub1b"), "EP": native_E(out, "P"), "ET1": native_E(out, "track1_combined")}
    results = {ln: evaluate_line(ln, natives[ln], ctx, rows, age, qualify_cap=10 if ln == "ENG" else (None if ln == "LIT" else 5)) for ln in LINES}
    near = []
    for n in b1_near:
        st, why = classify_drug(n["drug"], ctx)
        near.append({**n, "status": st, "reason": why})
    sets = {ln: {c["drug"] for c in results[ln]["shown"]} for ln in LINES}
    qual_sets = {ln: {c["drug"] for c in results[ln]["shown"]} | set(results[ln]["also_qualifying"]) for ln in LINES}
    conv = convergence(sets, len(pool) or 224)
    conv_all = convergence(qual_sets, len(pool) or 224)
    stab = comp.get("side_by_side", {}).get("completed", {})
    candidate_rows = {}
    for ln in LINES:
        for c in results[ln]["shown"]:
            d = c["drug"]
            candidate_rows.setdefault(d, {"drug": d, "lines": [], "age_bands": c["age_bands"], "age_union": c["age_union"], "smaa_p_top5": c["smaa_p_top5_cross_check"], "smaa_rank": c["smaa_rank_by_p_top5"], "safety_flags": c["safety_flags"]})["lines"].append(ln)
    first = ME.load_json(out / "candidates_by_line_mva.json")
    changes = changes_vs_first(first, results, ctx) if first else None
    return {
        "label": "hypotheses for follow-up; each line proposes its own candidates; convergence is information, not a requirement; no claim of efficacy, no recommendation, no dose",
        "protocol": "protocol/addendum_multi_evidence.md sections 12 and 13 (version 2: corrected safety screen and age readings)", "changes_vs_first_output": changes, "filters": {"pool_size": len(pool), "approval_unconfirmed_28": sorted(ctx["unconfirmed"]), "hard_excluded": ctx["excluded"],
                                                                              "phenotype_keywords": PHENOTYPE_KEYWORDS, "contra_only": CONTRA_ONLY, "searched_sections": list(SEARCHED), "cytotoxic_context_n": len(case["caution_names"])},
        "lines": results, "b1_near_threshold_information_only": near, "convergence": conv, "convergence_including_also_qualifying": conv_all, "candidates": candidate_rows,
        "cross_check": {"note": "multi-evidence SMAA (completed run) as a cross-check column only", "proposed_by_smaa": comp.get("proposed", []), "trametinib_survives": comp.get("trametinib_survives"),
                        "stability": {k: stab.get(k) for k in ("leave_one_out_proposed", "mc_retained", "equal_weights_top5", "literature_only_proposed", "track1_off_proposed")}},
        "track1_gene_block_statement": "Track 1 genes other than BUB1B are hypotheses to test (karyotype or single-cell sequencing for the megabase events; function tests for IGBP1), not findings.",
        "data_gaps": [], "n_native": {ln: len(natives[ln]) for ln in LINES}}


# ---------------------------------------------------------------- markdown

def to_markdown(r: dict[str, Any]) -> str:
    L = ["# Candidates by line (hypotheses for follow-up), version 2", "", f"*{r['label']}.*", "", f"Protocol: `{r['protocol']}`. Common filters: approved small molecules of the completed pool ({r['filters']['pool_size']} drugs); hard exclusions (the three hazard drugs, "
         "the growth-hormone class, drugs the project's layers exclude, and documented child-specific contraindications quoted from the label; broad label keywords only flag for manual verification); cytotoxic agents are listed as context/hazard under the line that found them. "
         f"{len(r['filters']['approval_unconfirmed_28'])} candidates the approval check could not confirm are listed separately: {', '.join(x.title() for x in r['filters']['approval_unconfirmed_28'])}.", ""]
    ch = r.get("changes_vs_first_output")
    if ch:
        L += ["## What changed relative to the first output (`candidates_by_line_mva.md`, kept)", "",
              "Corrections (addendum section 13): the keyword exclusion of the label screen is replaced by the project's hand-read layers as authority (keywords only flag, with the matched sentence), and every candidate has an age-band reading from the project's label functions.", "",
              "| Line | First output | Version 2 | Added | Removed | Age readings changed |", "|---|---|---|---|---|---|"]
        for ln in LINES:
            c = ch[ln]
            tt = lambda xs: ", ".join(x.title() for x in xs) or "none"  # noqa: E731
            L.append(f"| {ln} | {tt(c['first'])} | {tt(c['v2'])} | {tt(c['added'])} | {tt(c['removed'])} | {', '.join(c['age_reading_changed']) or 'none'} |")
        nd = ch["named_drugs"]
        L += ["", "Named drugs under the corrected screen: " + "; ".join(f"{d.title()}: screen result {v['screen_result']}" + (f" ({v['reason'][:230]})" if v['reason'] else "") + f"; candidate in {', '.join(v['candidate_in_lines']) or 'no line'}; also qualifying in {', '.join(v['also_qualifying_in_lines']) or 'no line'}" for d, v in nd.items()) + ".", ""]
    for ln in LINES:
        x = r["lines"][ln]
        L += [f"## {x['name']}", "", f"Native set: {x['n_native_set']} drugs; qualifying after the filters: {x['n_qualifying']} (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: {x['limit']}.", ""]
        if x["shown"]:
            L += ["| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |", "|---|---|---|---|---|---|---|--:|"]
            for c in x["shown"]:
                cons = ", ".join(c["consequences_addressed"]) or "none identified"
                if c["track1_genes_bearing"]:
                    cons += " (Track 1 genes bearing: " + "/".join(c["track1_genes_bearing"]) + ", hypothesis)"
                p = "-" if c["smaa_p_top5_cross_check"] is None else f"{c['smaa_p_top5_cross_check']:.3f} (rank {c['smaa_rank_by_p_top5']})"
                L.append(f"| {c['drug'].title()} | {c['native']} | {' / '.join(c['age_bands'][b] for b in ME.BANDS)} ({c['age_union']}) | {'; '.join(c['safety_flags'])} | {cons} | {c['evidence_tier']} | {c['refuting_experiment']} | {p} |")
            if x["also_qualifying"]:
                L += ["", "Also qualifying (not shown): " + ", ".join(d.title() for d in x["also_qualifying"]) + "."]
        else:
            L += ["No drug meets the line's native criterion after the filters. This is a result, not a failure of the table."]
        if ln == "B1" and r["b1_near_threshold_information_only"]:
            L += ["", "Below the native threshold, information only (not candidates): " + "; ".join(f"{n['drug'].title()} ({n['native']}; {n['status']})" for n in r["b1_near_threshold_information_only"]) + "."]
        if x["found_but_not_proposed"]:
            show = x["found_but_not_proposed"][:12]
            L += ["", "Found by the line but not proposed (context/hazard or filtered): " + "; ".join(f"{o['drug'].title()} [{o['status']}: {o['reason']}]" for o in show) + (f"; and {len(x['found_but_not_proposed']) - len(show)} more (see the JSON)." if len(x['found_but_not_proposed']) > len(show) else ".")]
        L.append("")
    cv = r["convergence"]
    L += ["## Convergence matrix (candidate shown by the line: x)", "", "Convergence across lines is information and is never required.", ""]
    if cv["matrix"]:
        L += ["| Drug | " + " | ".join(LINES) + " | Lines |", "|---|" + "---|" * len(LINES) + "--:|"]
        for d, row in sorted(cv["matrix"].items(), key=lambda kv: (-sum(kv[1].values()), kv[0])):
            L.append(f"| {d.title()} | " + " | ".join("x" if row.get(ln) else "" for ln in LINES) + f" | {sum(row.values())} |")
    else:
        L.append("No candidates.")
    L += ["", f"Distinct candidates: {cv['n_candidates_distinct']}; unique to one line: {cv['n_unique_to_one_line']}; in two or more lines: {', '.join(x.title() for x in cv['in_two_or_more_lines']) or 'none'}.", "",
          "Overlap between lines against chance (hypergeometric, universe = completed pool of " + str(cv["universe_n"]) + " drugs; pairs with candidates in both lines; no multiplicity claim):", ""]
    if cv["pairs"]:
        L += ["| Lines | Sizes | Overlap | Shared | p |", "|---|---|--:|---|--:|"]
        for p in cv["pairs"]:
            L.append(f"| {p['lines'][0]} - {p['lines'][1]} | {p['sizes'][0]}, {p['sizes'][1]} | {p['overlap']} | {', '.join(x.title() for x in p['shared']) or '-'} | {p['p_hypergeometric']:.3g} |")
    else:
        L.append("Fewer than two lines have candidates.")
    L += ["", "## Age-band view (label support, a filter and a column, not evidence of efficacy)", "", "| Drug | Lines | 0-4 | 5-11 | 12-17 | Unknown age (union) |", "|---|---|---|---|---|---|"]
    for d, c in sorted(r["candidates"].items()):
        L.append(f"| {d.title()} | {', '.join(c['lines'])} | " + " | ".join(c["age_bands"][b] for b in ME.BANDS) + f" | {c['age_union']} |")
    L += ["", "Age readings come from `engine_case_study.label_support` (age_label_table.row_for, manual review where it exists, otherwise the heuristic classifier). Section used and quote:", ""]
    seen = set()
    for ln in LINES:
        for c in r["lines"][ln]["shown"]:
            d = c["drug"]
            if d in seen:
                continue
            seen.add(d)
            a = c.get("age_reading")
            if not a:
                continue
            how = "manual review" if "manual review" in (a.get("basis") or "") else "heuristic classifier, not hand-verified"
            if not a.get("retrievable"):
                L.append(f"- {d.title()}: label not retrievable ({'; '.join(a.get('failures') or [])[:120]}).")
            else:
                L.append(f"- {d.title()}: status {a['label_status']} ({how}); section {a.get('quote_section') or 'none found'}: \"{(a.get('quote') or 'no paediatric-use sentence found').strip()[:260]}\"")
    cc = r["cross_check"]
    L += ["", "## Cross-check with the multi-evidence SMAA (completed run)", "", f"The SMAA rule (P(top 5) >= 0.5) proposed: {', '.join(x.title() for x in cc['proposed_by_smaa']) or 'none'}; trametinib survives: {cc['trametinib_survives']}. The probability of each candidate is in its table row and is not a criterion here. "
          "Stability of the SMAA run (completed): "]
    st = cc["stability"]
    if st.get("mc_retained"):
        L.append("- label errors 10/20/30 %: baseline candidates retained " + ", ".join(f"{k} % {v[0]:.2f} of {v[1]}" for k, v in st["mc_retained"].items()) + " (trivial when the SMAA candidate has no literature labels).")
    L.append("- equal weights, top 5: " + ", ".join(x.title() for x in st.get("equal_weights_top5") or []) + "; literature-only weights, proposed: " + ", ".join(x.title() for x in st.get("literature_only_proposed") or []) + "; gene block off, proposed: " + ", ".join(x.title() for x in st.get("track1_off_proposed") or []) + ".")
    L += ["", "## Limits", "", "- Each line has its own validity and limits (shown with each table). The lines differ in what they measure, so different candidates are expected; their overlap is reported against chance, not required.",
          "- Candidates are hypotheses for follow-up: no claim of efficacy, no recommendation, no dose. Label keyword exclusions are a screen on label text, not a clinical judgement; a label that could not be read is flagged and not excluded.",
          "- " + r["track1_gene_block_statement"], ""]
    for g in r.get("data_gaps", []):
        L.append(f"- Data gap: {g}")
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("out"))
    ap.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    ap.add_argument("--rationale", type=Path, default=Path("config/repositioning_rationale.json"))
    ap.add_argument("--gene-block", type=Path, default=Path("config/multi_evidence_track1_gene_block.json"))
    args = ap.parse_args()
    res = run(args.out, args.universe, args.rationale, args.gene_block if args.gene_block.exists() else None)
    unread = sorted(d for d in res["candidates"] if "label not read" in " ".join(res["candidates"][d]["safety_flags"]))
    if unread:
        res["data_gaps"].append("label not read for: " + ", ".join(unread))
    res["data_gaps"].append("lines with an empty native set are reported as empty; the B1 near-threshold drugs are information only")
    (args.out / "candidates_by_line_mva_v2.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    (args.out / "candidates_by_line_mva_v2.md").write_text(to_markdown(res), encoding="utf-8", newline="\n")
    print({ln: [c["drug"] for c in res["lines"][ln]["shown"]] for ln in LINES})


if __name__ == "__main__":
    main()
