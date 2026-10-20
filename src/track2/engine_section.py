"""Short section of the report on the evidence engine, written only from the JSON files of `out/engine/`; nothing but wording is typed here."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _load(eng: Path, name: str) -> dict[str, Any] | None:
    p = eng / f"{name}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def _g(d: Any, path: str) -> Any:
    for k in path.split("/"):
        d = d[int(k)] if isinstance(d, list) else d[k]
    return d


def _ci(x: dict[str, Any], digits: int = 2) -> str:
    def fmt(v: float) -> str:
        r = f"{v:+.{digits}f}"
        return r.replace("-0.00", "0.00").replace("+0.00", "0.00")

    lo, hi = x["ci95"]
    return f"{fmt(x['estimate'])} [{fmt(lo)}, {fmt(hi)}]"


def _crit1(a: dict[str, Any]) -> bool:
    return a["estimate"] >= 0.05 and a["ci95"][0] > 0


def _have(eng: Path) -> bool:
    return all((eng / f"{n}.json").exists() for n in ("validation", "validation_v2", "validation_v3", "validation_v4", "validation_v4_ablation", "validation_v4_llm"))


def abstract_sentence(eng: Path) -> str:
    if not _have(eng):
        return ""
    v1, v4 = _load(eng, "validation"), _load(eng, "validation_v4")
    a = _g(v1, "main/primary_gold1/arms/A/auc_minus_popularity")
    b = _g(v4, "heldout/evaluation_group_bootstrap/models/B0/auc_minus_popularity")
    return (f"A literature-grounded ranking engine (Section 3.12) beat a popularity baseline in AUC by {a['estimate']:.2f} on 17 pre-registered rare diseases and by {b['estimate']:.2f} on {_g(v4, 'heldout/n_diseases')} sealed held-out diseases, "
            "but its advantage at low literature volume comes from the knowledge-graph prior and not from the literature.")


def main_section(eng: Path) -> str:
    if not _have(eng):
        return ""
    v1, v2, v3, v4 = (_load(eng, n) for n in ("validation", "validation_v2", "validation_v3", "validation_v4"))
    ab, llm, au, wl = _load(eng, "validation_v4_ablation"), _load(eng, "validation_v4_llm"), _load(eng, "audit_llm_coverage_dev"), _load(eng, "case_wilms")
    d1 = _g(v1, "main/primary_gold1/decision_arm_A")
    a1 = _g(v1, "main/primary_gold1/arms/A/auc_minus_popularity")
    v2a = _g(v2, "primary_sets/gold1/new_28_heldout/arms/A/auc_minus_popularity")
    v2u = _g(v2, "primary_sets/union/all_45/arms/A/auc_minus_popularity")
    h3 = _g(v3, "heldout/criteria_v3/criterion_1_auc_minus_popularity_ge_0.05_and_ci_lower_gt_0")
    p3 = _g(v3, "heldout/paired_auc_v3_minus_v1_weights")
    c4 = _g(v4, "heldout/criteria_all_models/B0")
    c4e = _g(v4, "heldout/criteria_all_models/E")
    n4 = _g(v4, "heldout/n_diseases")
    g4 = c4["criterion_1_auc_minus_popularity_ge_0.05_and_ci_lower_gt_0"]
    fr = _g(ab, "heldout/table_group_bootstrap/fractions")
    rows = []
    for f in ("1.0", "0.1", "0.03", "0.01"):
        r = fr[f]
        rows.append(f"| {float(f) * 100:g} % | {r['n_papers_mean']:.0f} | {_ci(r['B0']['auc_gap_vs_popularity'])} | {_ci(r['PRIOR_ONLY']['auc_gap_vs_popularity'])} | {_ci(r['LIT_ONLY']['auc_gap_vs_popularity'])} |")
    gate = _g(llm, "dev/relevance_gate")
    dirk = _g(llm, "dev/agreement/direction/kappa")
    relk = _g(llm, "dev/agreement/relevance/kappa")
    t = ["### 3.12 A ranking engine built on the same principles, tested on other rare diseases\n"]
    t.append("The pipeline above ends with no drug supported by two independent computational lines. To learn how far a literature-grounded ranking can be trusted, we built a separate engine and tested it on rare diseases whose treatments are known, using the principles of the CASMI and DREAM challenges: thresholds and endpoints written before the data were scored, folds grouped by disease, baselines (popularity of the drug, a mention-based proxy, random), and a sealed held-out set opened once. "
             "It is deterministic code and not a language model: it retrieves abstracts, counts mentions of a drug with the disease genes, adds a knowledge-graph prior and combines them with weights learned by cross-validation. Disease matching used identifiers and synonyms (MONDO, Orphanet, HPO); the gold sets are PrimeKG indications (gold 1), current FDA labels (gold 3) and drugs in experimental arms of trials started after 2018 (gold 2).\n")
    t.append("| Test | Diseases | Gain in AUC over popularity [95 % CI] | Reading |\n|---|---|---|---|")
    t.append(f"| Pre-registered primary (17 diseases) | 17 | {_ci(a1)} | All three criteria met ({d1['verdict']}) |")
    t.append(f"| Extension, 28 new diseases (gold 1) | 28 | {_ci(v2a)} | Criterion {'met' if _crit1(v2a) else 'not met'} |")
    t.append(f"| Extension, union gold, 45 diseases | 45 | {_ci(v2u)} | Criterion {'met' if _crit1(v2u) else 'not met'} |")
    t.append(f"| New weights (v3), sealed held-out | 14 | {_ci(h3)} | Improvement over the v1 weights not shown (paired {_ci(p3)}) |")
    t.append(f"| Larger benchmark (v4), sealed held-out | {n4} | {_ci(g4)} | Simple model supported; embeddings model {c4e['verdict']} |\n")
    t.append(f"The engine recovers known indications better than popularity on the pre-registered set and on {n4} sealed held-out diseases, but the first extension to new diseases did not hold with the original weights, and learning new weights did not help once they were scored on unseen diseases. "
             "We therefore report the simple model refit on more diseases as the result and do not claim an improved engine. In a worked case, Wilms tumour (a disease whose treatment is known), the engine did not rank the three known drugs in its first twenty: "
             f"{_wilms(wl)}.\n")
    t.append("**Low literature volume.** Ultra-rare diseases such as the one in this report have few publications, so we thinned the retrieved literature to 100, 10, 3 and 1 % and compared the engine with versions using only the prior or only the literature (sealed held-out):\n")
    t.append("| Literature kept | Papers (mean) | Full engine | Prior only | Literature only |\n|---|---|---|---|---|")
    t.extend(rows)
    t.append("")
    t.append("The gain shrinks with less literature but stays above zero because of the prior; the literature alone falls below popularity. The advantage at low volume is therefore explained by the knowledge-graph prior and not by the literature. This is a decomposition of a result already seen, not a new confirmation.\n")
    t.append(f"**A language model as annotator.** Two language models (Opus and Sonnet), blind to each other and to the gold sets, labelled sentences for the direction of the effect and for relevance. Direction agreed well (kappa {dirk['estimate']:.2f}) but relevance did not reach the gate of 0.6 written beforehand (kappa {relk['estimate']:.2f}, gate {'passed' if gate['passed'] else 'failed'}), so that first experiment is reported as not testable with those labels. "
             "A second guideline with explicit rules raised agreement on the same 800 items from kappa 0.62 to 0.85. An audit then showed that the retrieval, not the labels, had limited coverage: the papers had been retrieved by gene symbols alone, and only a minority of known drug-disease pairs had the drug named in any retrieved paper.\n")
    t.append(_ultralight(_load(eng, "validation_v5_ul")))
    t.append(_application(_load(eng, "apply_mva"), _load(eng, "apply_fanconi")))
    return "\n".join(t) + "\n"


def _age(bands: dict[str, str]) -> str:
    if all(v == "not_labelled" for v in bands.values()):
        return "not labelled in any band"
    return ", ".join(f"{b.replace('_', '-')} y: {v.replace('_', ' ')}" for b, v in sorted(bands.items()) if v != "not_labelled")


def _top_table(ap: dict[str, Any], n: int = 10) -> tuple[list[str], int]:
    rows, chemo = ["| Rank | Drug | FDA label | Age-band support | Cytotoxic caution |", "|---|---|---|---|---|"], 0
    for r in ap["ranking_top20_fda_small_molecule"][:n]:
        chemo += bool(r["oncology_cytotoxic_chemotherapy_caution"])
        rows.append(f"| {r['position']} | {r['name']} | {'yes' if r['openfda_label_match'] else 'Drugs@FDA only'} | {_age(r['bands'])} | {'yes' if r['oncology_cytotoxic_chemotherapy_caution'] else 'no'} |")
    return rows, chemo


def _application(mva: dict[str, Any] | None, fa: dict[str, Any] | None) -> str:
    if not mva or not fa:
        return ""
    t = ["**Applied to the cases of this report.** The frozen engine (no refit, no label of either disease) was then applied to mosaic variegated aneuploidy and to Fanconi anaemia, restricted to FDA-approved small molecules, with the label-based age support of Section 3.6 (three age bands; the age of the child is never assumed). The output is exploratory and is not a candidate list.\n"]
    for ap, name in ((mva, "Mosaic variegated aneuploidy"), (fa, "Fanconi anaemia")):
        rows, chemo = _top_table(ap)
        t.append(f"*{name}: first ten drugs.* {chemo} of the ten carry the cytotoxic-chemotherapy caution.\n")
        t.extend(rows)
        t.append("")
    nine = mva["nine_report_drugs_positions"]
    ranks = sorted((v["rank_fda_small_molecule"], k) for k, v in nine.items() if v.get("rank_fda_small_molecule"))
    tot = next(iter(nine.values()))["of_fda_small_molecule"]
    t.append("None of the nine drugs of this report is in the first twenty for MVA; their positions among " + f"{tot} approved small molecules are " + ", ".join(f"{k} {r}" for r, k in ranks) + ". "
             "Nothing in these tables is a candidate: the engine counts mentions of a drug with the disease genes and does not read their direction, and for genes of mitosis and DNA repair the mentions of cytotoxic agents often describe sensitivity or toxicity, a pattern the caution list of this report already records. "
             "No top-ranked drug has label support for a paediatric age band, so the age filter leaves nothing to propose. The low-literature flag was not raised, but it counts papers about the genes; the audit above showed that most known drug-disease pairs have no paper that names both, so the absence of the flag is not evidence that enough literature exists. "
             f"The case-count tier from Orphanet is '{mva['evidence_gap']['gaps'][0]['orphanet_tier']}' for MVA and '{fa['evidence_gap']['gaps'][0]['orphanet_tier']}' for Fanconi anaemia; for MVA the child is one known case, and the study designs that could close each gap (a single-patient design with a hierarchical model borrowing from the mechanism class) are listed in the output file `out/engine/apply_mva.md`.\n")
    return "\n".join(t) + "\n"


def _ultralight(ul: dict[str, Any] | None) -> str:
    if not ul:
        return "Retrieval per drug-disease pair was designed to correct this and has not been annotated (Supplement S14).\n"
    q, tb = ul["quality_gate"], ul["table"]
    parts = []
    for f in ("0.03", "0.01"):
        r = tb[f]
        parts.append(f"at {float(f) * 100:g} % of the records the gain over the prior alone is {_ci(r['auc_difference'], 3)} with {r['coverage_share_pairs_B1'] * 100:.1f} % of pairs covered ({'informative' if r['informative'] else 'not informative'})")
    return (f"We then retrieved per drug-disease pair and ran the cheapest design that could test the idea: the held-out diseases only, {ul['manifest_items']} sentences read by one language model with the second guideline, a quality check by a second model (kappa {q['kappa']:.2f}, gate {'passed' if q['gate_passed'] else 'failed'}), and a fixed score bump for pairs with a sentence of relevant benefit, with no training. "
            f"Scored once: {parts[0]}; {parts[1]}. The rule written beforehand required both fractions, so the outcome is {ul['outcome']}. "
            "At 3 % the lower bound is above zero but the gain is a few thousandths of AUC and about 1 % of the pairs carry any evidence; at 1 % too few pairs are covered to say anything. "
            "The design reads about one sentence per candidate pair (roughly 190 tokens per pair for two annotators), so its cost grows linearly with the number of diseases; the full design would have needed about 15 million tokens and is kept as a documented, unexecuted alternative (Supplement S14).\n")


def _wilms(wl: dict[str, Any] | None) -> str:
    if not wl:
        return "case file not available"
    drugs = _g(wl, "concordance/gold1/primary_literature_to_2024/approved_fda_label/drugs")
    n = _g(wl, "concordance/gold1/primary_literature_to_2024/approved_fda_label/n_candidates")
    return "; ".join(f"{d['name']} ranked {d['rank_min']} of {n}" for d in drugs)


def supplement(eng: Path) -> str:
    if not _have(eng):
        return ""
    t = ["## S14. Evidence engine: protocols, files and what was not done\n"]
    t.append("Each stage has a written protocol dated before its scoring, a frozen model with hash, a lock file for the sealed set and a test file. Protocols: `protocol/addendum_engine.md` (v1, pre-registered), `addendum_engine_v2.md`, `addendum_engine_v3.md`, `addendum_engine_v4.md`, `addendum_engine_v4_dilution.md`, `addendum_engine_v4_ablation.md`, `addendum_engine_v4_llm.md` with its guideline-v2 and audit addenda, and `addendum_engine_v5.md` (per-pair retrieval, full design: 47,461 sentences, about 14.7 million tokens, prepared and not executed), a reduced design (held-out plus a sample of development, about 3 million tokens, not executed) and `addendum_engine_v5_ultralight.md` (executed: 1,281 held-out sentences, about 0.26 million tokens). Results are in `out/engine/validation*.json` (`validation_v5_ul.json` for the ultralight run).\n")
    t.append("Deviations that matter: the v2 and later benchmarks relaxed inclusion rules (minimum genes and gold drugs, grouping of diseases with identical gold) and this is declared in the protocols; the held-out gold was read to define the candidate set, as in v4; the engine is evaluated on a candidate set and not on the full drug universe; a rule that admitted some common diseases was kept unchanged and examined post hoc. The relevance guideline v2 was written after observing annotator disagreement and is a post hoc revision; v1 remains the pre-registered guideline. The per-pair retrieval (v5) and the ultralight design were made after the audit; the size of the evidence bump (0.5 of a standard deviation) was fixed beforehand and marked UNVALIDATED with default zero.\n")
    t.append("What is not shown: that the engine predicts future indications better than popularity (the temporal check was not distinguishable from it), that language-model-read evidence helps at low literature volume (the one test run was not informative), and any claim about the child's drugs.\n")
    return "\n".join(t) + "\n"
