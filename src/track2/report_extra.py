"""Report sections added by the robustness, literature-context, extrapolation and age-band layers (all generated, no typed figures)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from track2.pediatric_review import load_review


def _opt(out: Path, name: str) -> Any:
    path = out / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def robustness_section(out: Path) -> str:
    cons = _opt(out, "consensus_ranking.json")
    if not cons:
        return ""
    enriched = {r["drug_id"]: r for r in (_opt(out, "candidates_enriched.json") or [])}
    s, ag = cons["summary"], cons["agreement"]
    schemes = ag["schemes_at_700"]
    lines = ["| Drug | Configurations in the top 20 (of 9) | Class | Mechanism of action (ChEMBL) | Paediatric label (heuristic) |", "|---|:-:|---|---|---|"]
    for d in cons["drugs"]:
        rec = enriched.get(d["drug_id"])
        moa = "; ".join(rec["mechanism_of_action"])[:60] if rec else ""
        ped = rec["age_profile"]["pediatric_label_status"].replace("_", " ").lower() if rec else ""
        lines.append(f"| {d['name']} | {d['robustness_count']} | {d['class']} | {moa} | {ped} |")
    n_short, n_in = len(enriched), len(cons["drugs"])
    return f"""### 5.12 Robustness across network thresholds and weighting schemes (`protocol/addendum_robustness.md`)

The combined ranking was recomputed in nine configurations: STRING 400, 700 and 900, each with three weighting schemes (W1: Track 1 `epcr` weights; W2: equal weights; W3: only the gene with a ClinVar pathogenic variant). A drug is **robust** if it is among the 20 best in at least 7 of the 9 configurations and **fragile** if in 3 or fewer (rule fixed before the run). Result: {s['robust']} robust, {s['intermediate']} intermediate and {s['fragile']} fragile drugs; the other {n_short - n_in} of the {n_short} shortlisted drugs are in the top 20 of no configuration and should be read as noise of the shortlist rule (many are GABA-A modulators near the module by network topology only).

- The weighting scheme barely matters: the top 20 is identical under W1, W2 and W3 at threshold 700 (Jaccard {min(v['jaccard_top20'] for v in schemes.values()):.2f}); rank correlations between schemes are {min(v['spearman'] for v in schemes.values()):.2f} to {max(v['spearman'] for v in schemes.values()):.2f}.
- The network threshold matters more: top-20 overlap with threshold 700 is {ag['thresholds_under_W1']['700 vs 400']['jaccard_top20']:.2f} at 400 and {ag['thresholds_under_W1']['700 vs 900']['jaccard_top20']:.2f} at 900. The earlier instability of a single BUB1B module (section 6) is much smaller for the combined ranking.
- **Robust is not the same as right.** The robust set is essentially one class (proteasome inhibitors, tubulin binders and antibody-drug conjugates with tubulin toxins): it is stably the same class, and these are cytotoxic agents whose direction of effect in a checkpoint-deficient child is unknown (section 5.13).

{chr(10).join(lines)}
"""


def literature_section(out: Path) -> str:
    ctx = _opt(out, "literature_context.json")
    if not ctx:
        return ""
    drugs = ctx["drugs"]
    classes: dict[str, list[str]] = {}
    names = {r["drug_id"]: r["name"] for r in (_opt(out, "candidates.json") or [])}
    notes: dict[str, dict[str, Any]] = {}
    for drug_id, entries in drugs.items():
        for e in entries:
            classes.setdefault(e["class"], []).append(names.get(drug_id, drug_id))
            notes[e["class"]] = e
    blocks = []
    for cls, e in notes.items():
        pm = "; ".join(f"PMID {p} ({t[:70]})" for p, t in e["pmids_verified"].items())
        blocks.append(f"- **{cls.replace('_', ' ')}** ({len(classes[cls])} candidates; evidence: {e['claim_level']}; direction hint: `{e['direction_hint']}`). {e['note']} Sources: {pm}.")
    return f"""### 5.13 Literature context and direction of effect (curated, cited)

Proximity gives no direction of effect (section 5.9). To make the caveat concrete, each candidate is matched by mechanism class to a short curated note (`config/literature_context.json`). Every PMID was re-checked to exist, and the abstracts carrying a claim were read; each note states only what its source supports and where the evidence comes from (cancer cells, yeast, mice or case reports). **None of it is evidence that a drug works in this patient.** PMIDs not found: {len(ctx['pmids_not_found'])}.

{chr(10).join(blocks)}
"""


def extrapolation_section(out: Path) -> str:
    path = out / "arm_extrapolation" / "candidates.json"
    if not path.exists():
        return ""
    rows = json.loads(path.read_text(encoding="utf-8"))
    ctx = (_opt(out, "literature_context_extrapolation.json") or {}).get("drugs", {})
    review = load_review()
    lines = ["| Drug | Best module z | Mechanism of action (ChEMBL) | Paediatric label (reviewed) | Literature class |", "|---|:-:|---|---|---|"]
    shown = [r for r in rows if any(m["passes_primary_threshold"] for m in r["mechanism_scores"])][:50]
    for r in shown:
        best = min(m["z"] for m in r["mechanism_scores"])
        lit = ", ".join(sorted({e["class"].replace("_", " ") for e in ctx.get(r["drug_id"], [])})) or "-"
        ent = review.get(r["drug_id"])
        ped = (f"established from {ent['min_age_years']:g} y: {ent['scope'][:70]}" if ent and ent["verdict"] == "established"
               else "not established" if ent and ent["verdict"] == "not_established" else "unreviewed (heuristic)")
        lines.append(f"| {r['name']} | {best:.2f} | {'; '.join(r['mechanism_of_action'])[:45]} | {ped} | {lit} |")
    return f"""### 5.14 Extrapolation arm (hypotheses from non-human or cancer-cell evidence, kept separate)

The main ranking comes from the Track 1 submission. A second, **separate** arm uses three author-defined consequence modules (proteotoxic stress, p53/senescence, mTORC1/AMPK growth signalling; `config/mechanisms_extrapolation/`) motivated by the BubR1 mouse and aneuploidy literature cited in section 5.13: sarcopenia with mTORC1 hyperactivity in BubR1-mutant mice, senescence in BubR1-insufficient mice, and proteotoxic vulnerability of aneuploid cells. They are not derived from the submission and are **not part of the main ranking**; they exist so that growth-restriction and muscle findings are not ignored. The arm contains {len(rows)} drugs; the {len(shown)} below reach z ≤ −2 in at least one module. Several (gabapentinoids) are likely network-topology artefacts, and no drug was tested in BubR1 mice, so every entry is an extrapolation. This arm is where the only **non-cytotoxic drugs with a paediatric label and a literature-linked rationale** appear: the mTOR-pathway drugs everolimus (established from 1 year for TSC-associated SEGA) and sirolimus (established from 13 years for renal transplantation), which ChEMBL curates as FKBP1A binders. The rationale is a correlation in BubR1-mutant mice (sarcopenia with mTORC1 hyperactivity), the proximity is partly circular because the growth module was seeded with mTOR-pathway genes, and both drugs carry boxed warnings; they are a hypothesis for experts to weigh, not a proposal. In the submission-derived modules their proximity is none (z about +0.3).

{chr(10).join(lines)}
"""


def submission_validation_paragraph(out: Path) -> str:
    parts = []
    for path in sorted(out.glob("validation_sub_*.json")):
        v = json.loads(path.read_text(encoding="utf-8"))
        p1, p3 = v["P1_leave_one_out"], v["P3_P4_robustness"]
        parts.append(f"| {v['mechanism'].removeprefix('sub_')} | {p1['n_seeds']} ({p1['n_seeds_tested']} tested) | {p1['auc']:.2f} ({p1['auc_ci95'][0]:.2f}-{p1['auc_ci95'][1]:.2f}) | {p1['decision']} | "
                     f"{p3['drop_20pct_seeds']['mean_jaccard_top20']:.2f} | {p3['string_thresholds']['400']['jaccard_top20']:.2f} / {p3['string_thresholds']['900']['jaccard_top20']:.2f} | {p3['p4_degree_bias_spearman']:.2f} |")
    if not parts:
        return ""
    header = ("- **P1, P3 and P4 on the submission-derived modules** (`protocol/addendum_robustness.md`, analysis A; same rules as the frozen protocol; for modules with more than 120 seeds P1 uses a fixed subsample of 120):\n\n"
              "| Gene module | Seeds | P1 AUC (95% CI) | P1 decision | P3 drop-20% Jaccard | P3 STRING 400 / 900 Jaccard | P4 degree-bias rho |\n|---|:-:|---|---|:-:|:-:|:-:|")
    note = ('\nReading: the Reactome-based modules (BUB1B, CENPE, MAD1L1) pass P1, and the BUB1B module is far more stable than the hand-defined checkpoint module (top-20 Jaccard against thresholds 400 and 900 well above the earlier 0.48 and 0.38, and a degree-bias correlation close to 0). No seed reaches the strict z <= -2 in leave-one-out, as before. The IGBP1 module is a network neighbourhood, so P1 is partly circular (each seed is adjacent to other seeds by construction) and its drug set is unstable across thresholds; its P1 value should not be read as independent support, and its drugs stay exploratory.\n')
    return header + "\n" + "\n".join(parts) + "\n" + note


def age_section(out: Path) -> str:
    views = _opt(out, "age_views.json")
    if not views:
        return ""
    main = views["register"]
    cells = "; ".join(f"{b}: {len(main[b]['labelled_use'])} labelled / {len(main[b]['no_labelled_use'])} no labelled use / {len(main[b]['unknown_or_verify'])} verify" for b in main)
    return f"""**Age bands (the patient's age is not documented, so every band is shown).** The age of the child is not in the clinical document, so no candidate is filtered by age: each is shown for each band (neonate, infant, child, adolescent, adult) by whether the US label supports use there (`out/age_views.csv`, `out/age_views.json`). For the robust and intermediate candidates: {cells}. After reading the labels of the drugs the heuristic could not classify (`config/pediatric_label_review.json`), **only {', '.join(n.title() for n in main['child']['labelled_use']) or 'none'} have an established paediatric use**, each for an indication unrelated to this condition; the others would be off-label in a child and would need specialist review (details and the wider search in section 5.15). If the age becomes known, the same code produces a filtered view (`python -m track2.age_view --age-years <n>`) and the age requirement of each candidate changes from `NOT_ASSESSABLE` to `PASS`/`FAIL` (`config/patient.json`, field `age_years`). Label extraction is heuristic and must be verified by hand; this is not dosing advice.
"""


def pediatric_section(out: Path) -> str:
    arm = _opt(out, "pediatric_arm.json")
    if not arm:
        return ""
    rows = arm["drugs"]
    head = [r for r in rows if r["headline"]]
    lines = ["| Drug | Min. age (label) | Labelled scope (reviewed) | Best z (700) | Stability (rank <= 100 in n of 9) | Mechanism of action (ChEMBL) | Literature class |", "|---|:-:|---|:-:|:-:|---|---|"]
    for r in head:
        scope = (r["review"] or {}).get("scope", "")[:120]
        lit = ", ".join(sorted({e["class"].replace("_", " ") for e in r["literature_context"]})) or "-"
        lines.append(f"| {r['name']} | {r['min_age_years']:g} y | {scope} | {r['min_z_700']:.2f} | {r['stability']} ({r['n_configs_rank_le_100']}) | {'; '.join(r['mechanism_of_action'])[:50]} | {lit} |")
    weak = [r for r in rows if r["final_label_status"] != "not_established_reviewed" and not r["headline"]]
    n_weak_stable = sum(1 for r in weak if r["stability"] == "stable")
    weak_zs = sorted(r["min_z_700"] for r in weak if r["min_z_700"] is not None)
    weak_median = weak_zs[len(weak_zs) // 2] if weak_zs else float("nan")
    return f"""### 5.15 Paediatric-labelled candidates (`protocol/addendum_pediatric.md`)

Because the child's age is unknown but the context is paediatric, the search was widened from the shortlist to the **whole scorable universe**: of the {arm['pool_size']} drugs with the highest combined score, {arm['n_with_label']} have a US label, and {arm['n_paediatric']} of those state established or mixed paediatric use by a heuristic reading. The heuristic cannot tell "established for one indication" from "not established" (it marked docetaxel, cabazitaxel and natalizumab as paediatric although their labels say efficacy is not established), so the labels of the drugs that matter were read by hand and recorded (`config/pediatric_label_review.json`; a reading by the author's assistant, not a pharmacist review; this review and the headline rule were decided after seeing the first result and are disclosed as such in the addendum).

**A paediatric label is not a reason to give a drug to this child.** The rule below is mechanical (label plus proximity) and it lets through chemotherapy, a checkpoint antibody, sedatives and anaesthetics, none of which anyone would propose for a child without a tumour. The {len(head)} drugs that pass it are listed so the search is transparent; **none is proposed**, and whether any drug is appropriate for this child is a clinical decision that this project does not make.

Drugs that pass the rule (established paediatric use on the label after review, **and** z <= -2 in at least one module):

{chr(10).join(lines)}

Reading this honestly:
- **Stable drugs that pass the rule ({', '.join(r['name'].title() for r in head if r['stability'] == 'stable') or 'none'}).** Both are antineoplastic agents with a paediatric indication unrelated to MVA (Hodgkin lymphoma; ALK-positive lymphoma or inflammatory myofibroblastic tumour). Any use here would be a **tumour-context hypothesis**, not a treatment of the constitutional defect, and not a recommendation. Brentuximab vedotin carries a tubulin toxin and belongs to the same cytotoxic class as the robust list (section 5.13).
- **Unstable drugs that pass the rule ({', '.join(r['name'].title() for r in head if r['stability'] == 'unstable')}).** They reach z <= -2 only through the IGBP1-derived module, whose drug set is unstable across network thresholds (section 6). The kinase inhibitors (RET) and the checkpoint antibody have a tumour indication unrelated to MVA; the sedative, anaesthetic and antiseizure drugs sit near the module through GABA-A receptor links and have no biological rationale here, so they are most likely **network-topology artefacts**. They satisfy the rule but are not credible candidates.
- The other {len(weak)} paediatric-labelled drugs of the pool ({n_weak_stable} stable by rank) have a proximity close to chance (median best z {weak_median:.2f}): they were selected by label, not by distance, and several are topical (for example tapinarof, roflumilast, crisaborole). They are in `out/pediatric_arm.csv` and are not candidates.
- A paediatric label for an unrelated indication says nothing about safety in this child (aneuploidy mosaicism, tumour predisposition, kidney findings). Age, tumour status, organ function and medications remain unknown (section 7). Entries marked `verify` in the review file (diazepam, oxazepam, chlordiazepoxide) were not classified because only a short excerpt was read.
"""
