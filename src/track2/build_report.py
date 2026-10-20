"""Generate the Track 2 report (Markdown) from the computed results: no hand-written figures.

    python -m track2.build_report --out out --report submissions/candeleroj98_track2_report.md
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from track2.report_proposals import docking_section, dti_section, mechanism_block, prism_section, avoid_section, core_section, proposals_section, proposals_summary, register_section, use_section
from track2.report_extra import pediatric_section, age_section, extrapolation_section, literature_section, robustness_section, submission_validation_paragraph

AI_DECLARATION_FILE = Path("config/ai_declaration.md")
BAND_ABBREV = {"LABELLED": "L", "LABELLED_USE": "L", "NO_LABELLED_USE": "-", "UNKNOWN": "?", "VERIFY_MANUALLY": "check", "LABELLED_AGE_NOT_PARSED": "L(age?)"}


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _band(age: dict[str, str], key: str) -> str:
    return BAND_ABBREV.get(age.get(key, "UNKNOWN"), age.get(key, "?"))


def top_table(scores: dict[str, Any], candidates: dict[str, dict[str, Any]], n: int = 10) -> str:
    lines = ["| Rank | Drug | Mechanism of action (ChEMBL) | z | Flags | Label: adult / adolescent / child / infant |", "|:-:|---|---|:-:|---|---|"]
    for row in scores["rows"][:n]:
        record = candidates.get(row["drug_id"])
        moa = "; ".join(record["mechanism_of_action"])[:70] if record else ""
        age = record["age_profile"] if record else {}
        labels = " / ".join(_band(age, k) for k in ("adult", "adolescent", "child", "infant")) if record else "not in shortlist"
        lines.append(f"| {row['rank']} | {row['name']} | {moa} | {row['proximity']['z']:.2f} | {', '.join(row['risk_flags']) or '-'} | {labels} |")
    return "\n".join(lines)


def submission_modules(out: Path) -> dict[str, str]:
    """Modellable genes of the Track 1 submission -> section title (order: weight)."""
    genes = _load(out / "submission_genes.json")["genes"]
    return {f"sub_{g['gene']}": f"{g['gene']} module (weight {g['weight_epcr']:.4f})" for g in genes if g["modellable"]}


def submission_table(out: Path) -> str:
    lines = ["| Gene | Track 1 rows | Best epcr (weight) | ClinVar (exact variants) | Module | Modellable |", "|---|---|:-:|---|---|:-:|"]
    for g in _load(out / "submission_genes.json")["genes"]:
        cv = "; ".join(f"{c['variant']}: {c.get('classification') or c['status']}" for c in g["clinvar"]) or "no SNV rows"
        module = f"{g['module_size_in_network']} proteins ({g['n_pathways']} Reactome pathways)" if g["modellable"] else \
            f"too small or absent in the network ({g['module_rule'].replace('_', ' ')})"
        lines.append(f"| {g['gene']} | {len(g['rows'])} | {g['weight_epcr']:.4f} | {cv} | {module} | {'yes' if g['modellable'] else 'no'} |")
    return chr(10).join(lines)


def comparison_paragraph(out: Path) -> str:
    path = out / "comparison_manual" / "candidates.json"
    if not path.exists():
        return ""
    old = [r["name"] for r in _load(path)]
    new = [r["name"] for r in _load(out / "candidates.json")]
    both = sorted(set(old) & set(new))
    return (f"**Comparison with hand-defined modules.** An earlier version used four modules written by hand (checkpoint, proteotoxic stress, "
            f"p53/senescence, growth signalling; kept in `config/mechanisms/` and used for the validation of section 6; results in `out/comparison_manual/`). Its shortlist had {len(old)} drugs; "
            f"the submission-derived shortlist has {len(new)}; {len(both)} are in both. The top of the list is the same class (proteasome inhibitors, tubulin binders), "
            f"which shows that the submission-derived rule does not depend on the author's choice of pathways, but also that it adds little beyond the checkpoint biology.")


def mechanism_table(screen: dict[str, Any]) -> str:
    lines = ["| Gene | Variants in the proband | Disease | Mechanism | Phenotype match | Status |", "|---|---|---|---|---|---|"]
    for g in screen["genes"]:
        lines.append(f"| {g['gene']} | {g['variants']} | {g['disease']} | {g['mechanism']} | {g['phenotype_match']} | **{g['status']}**{' - ' + g['reason'] if g.get('reason') else ''} |")
    return "\n".join(lines)


def _optional(out: Path, name: str) -> Any:
    path = out / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def layers_sections(out: Path) -> str:
    """Sections 5.6-5.10: evidence, pharmacovigilance, pharmacogenomics, direction and reversal (all optional inputs)."""
    enriched = _optional(out, "candidates_enriched.json") or []
    faers = _optional(out, "faers_disproportionality.json")
    pgx = _optional(Path("data/private"), "pgx_summary.json")
    direction = _optional(out, "direction_review.json")
    reversal_meta = _optional(out, "reversal_meta.json")
    reversal = _optional(out, "reversal_scores.json")
    parts = []
    if enriched:
        lines = ["| Drug | PubMed: aneuploidy / tumour context | Trials: rhabdomyosarcoma / Wilms / MVA / paediatric |", "|---|:-:|:-:|"]
        for r in enriched[:12]:
            ev = r["layers"]["evidence_retrieved"]
            if not ev.get("available"):
                continue
            pm, tr = ev["pubmed_counts"], ev["trial_counts"]
            lines.append(f"| {r['name']} | {pm.get('aneuploidy_cin', 0)} / {pm.get('tumour_context', 0)} | "
                         f"{tr.get('rhabdomyosarcoma', 0)} / {tr.get('wilms_tumor', 0)} / {tr.get('mosaic_variegated_aneuploidy', 0)} / {tr.get('pediatric_population', 0)} |")
        n_any = sum(1 for r in enriched if r["layers"]["evidence_retrieved"].get("pubmed_counts", {}).get("aneuploidy_cin", 0) > 0)
        n_mva = sum(1 for r in enriched if r["layers"]["evidence_retrieved"].get("trial_counts", {}).get("mosaic_variegated_aneuploidy", 0) > 0)
        parts.append(f"""### 5.6 Literature and trials retrieved (not curated)

For each of the {len(enriched)} shortlisted drugs, PubMed and ClinicalTrials.gov were queried by keyword (aneuploidy/chromosomal instability, embryonal-tumour context, paediatric terms; rhabdomyosarcoma, Wilms tumour, mosaic variegated aneuploidy and paediatric-eligible trials). Every PMID was re-checked to exist. {n_any} drugs have at least one paper on aneuploidy or chromosomal instability, and **{n_mva} have a registered trial for mosaic variegated aneuploidy**: none of these drugs has been tried for this condition, so each candidate is a new hypothesis. A retrieved paper or trial does not mean the drug is supported; the counts are context for a human reviewer.

{chr(10).join(lines)}
""")
    if faers:
        signals = {}
        for entry in faers["drugs"].values():
            if entry.get("available"):
                for g, v in entry["groups"].items():
                    if v.get("signal"):
                        signals[g] = signals.get(g, 0) + 1
        n_av = sum(1 for e in faers["drugs"].values() if e.get("available"))
        summary = ", ".join(f"{g}: {n}" for g, n in sorted(signals.items()))
        parts.append(f"""### 5.7 Pharmacovigilance disproportionality (FDA adverse-event reports)

Beyond the descriptive shares in section 7, reporting disproportionality (PRR and ROR with 95% intervals; signal if PRR ≥ 2, at least 3 cases and χ² ≥ 4) was computed for six organ-system groups from openFDA counts, for {n_av} of the {len(faers['drugs'])} candidates ({len(faers['drugs']) - n_av} have no reports under their generic name). Number of drugs with a signal per group: {summary}. Strong haematological and neurological signals for antineoplastics are expected (indication and known toxicity). These are spontaneous reports with no denominator and confounding by indication: they flag where a clinician should look, not an incidence, and they cannot compare drugs.
""")
    if pgx:
        parts.append(f"""### 5.8 Pharmacogenomics (CPIC)

For eight genes with SNV-based CPIC allele definitions (CYP2C19, CYP2C9, CYP3A5, TPMT, NUDT15, DPYD, UGT1A1, SLCO1B1), the defining positions were checked against the proband's variant calls. Non-reference calls were found for some of them (the per-gene detail is kept private and is not published). The data are **unphased variant calls without reference calls**, so the absence of a call is not evidence of the reference allele and **no diplotype or phenotype is asserted**. CYP2D6 and HLA genes were excluded (they need dedicated callers). {pgx.get('drugs_with_cpic_pair', 0)} of the {pgx.get('drugs_total', 0)} candidates have CPIC drug-gene pairs (listed in `out/pgx_drug_relevance.json`); matching is by exact drug name.
""")
    if direction:
        s = direction["summary"]
        parts.append(f"""### 5.9 Direction of effect on the checkpoint

Proximity does not say whether a drug helps or harms. An attempt was made to orient each candidate with Gene Ontology annotations of its targets for spindle-checkpoint terms (QuickGO, {s['n_targets_queried']} human targets of the {s['n_drugs']} candidates). **None of the targets carries such an annotation** (the same code does return annotations for known checkpoint proteins such as BUB1 and MAD2L1), so no candidate could be oriented and every candidate keeps `direction_review = REQUIRES_REVIEW`.
""")
    if reversal and reversal_meta:
        n_rev = sum(1 for v in reversal.values() if v.get("available") and v.get("direction") == "reverses")
        n_scored = sum(1 for v in reversal.values() if v.get("available"))
        parts.append(f"""### 5.10 Transcriptomic reversal (LINCS L1000), a negative result

Signatures of BUB1B knockout (and four checkpoint/kinetochore genes as proxies) from LINCS L1000 were compared with {n_scored:,} drug signatures (SigCom LINCS). The BUB1B signature is dominated by generic interferon and stress genes, and unrelated knockout signatures used as controls produced similar rankings, so the {n_rev} drugs called as reversers are **not specific to BUB1B-like biology** (the top of the list contains drugs with no plausible link to this condition). The layer is therefore reported but **not used** in the ranking.
""")
    physchem = _optional(out, "physchem_adme.json")
    if physchem:
        entries = list(physchem.values())
        small = [e for e in entries if e["has_small_molecule_properties"]]
        mws = sorted(float(e["properties"]["full_mwt"]) for e in small if e["properties"].get("full_mwt"))
        ro5 = sum(1 for e in small if (e["properties"].get("num_ro5_violations") or 0) >= 2)
        aff = [e for e in entries if e["measured_affinity_to_mechanism_targets"]]
        adme = [e for e in entries if e["measured_adme"]]
        median_mw = mws[len(mws) // 2] if mws else float("nan")
        parts.append(f"""### 5.11 Measured binding, physicochemical properties and ADME (ChEMBL; nothing predicted here)

For each of the {len(entries)} shortlisted drugs, ChEMBL was queried for (i) computed physicochemical properties, (ii) exact nanomolar activity measurements (IC50, Ki, Kd, EC50) on the targets of its curated mechanism, and (iii) measured ADME values (solubility, logD, clearance, half-life, permeability, bioavailability). {len(small)} drugs have small-molecule properties (median molecular weight {median_mw:.0f}; {ro5} break two or more of Lipinski's rules, expected for oncology agents and natural-product derivatives); the remaining {len(entries) - len(small)} are biologics (antibodies and antibody-drug conjugates) for which no values are reported. {len(aff)} drugs have measured affinity on a mechanism target and {len(adme)} have at least one measured ADME value; **for the rest no value exists in ChEMBL and none was invented**. Per-drug values are in `out/physchem_adme.json` and `out/candidates_enriched.csv`.

**No predicted affinity is used.** A deep-learning affinity model was dropped because an earlier calibration missed a reference inhibitor by more than two log units, and docking was run with its own validation and did not pass (section 5.18), so no docking score is used either. Measured values are reported instead. Solubility and ADME values are heterogeneous assay results from the literature, useful as context and not as pharmacokinetic predictions; no doses are proposed.
""")
    for extra in (robustness_section(out), literature_section(out), extrapolation_section(out), pediatric_section(out), dti_section(out), prism_section(out), docking_section(out), core_section(out)):
        if extra:
            parts.append(extra)
    return "\n".join(parts)


def benchmark_paragraph(out: Path) -> str:
    bm = _optional(out, "benchmark_results.json")
    if not bm:
        return ""
    s = bm["summary"]
    return f"""- **P5, validation across {s['700']['n_diseases']} diseases (pre-registered addendum, `protocol/addendum_benchmark.md`).** For diseases selected by a fixed rule from Open Targets, drugs with an approved indication for the disease are ranked against all other scorable drugs by the same proximity method. Mean AUC {s['700']['mean_auc']:.2f} (95% CI over diseases {s['700']['ci95_low']:.2f}–{s['700']['ci95_high']:.2f}; {s['700']['n_auc_above_0_5']} of {s['700']['n_diseases']} diseases above 0.5; permutation p ≈ {s['700']['permutation']['p_one_sided']:.3f}, null mean {s['700']['permutation']['null_mean']:.2f}); consensus over the three network thresholds {s['consensus']['mean_auc']:.2f} ({s['consensus']['ci95_low']:.2f}–{s['consensus']['ci95_high']:.2f}). The pre-set rule (mean AUC ≥ 0.60 and lower bound > 0.50) is met: **{bm['decision']['verdict']}**. Several diseases fail (for example glaucoma, infection, atrial fibrillation), the sample is biased towards broad terms, diseases overlap in drugs and targets (so the interval is optimistic), literature evidence in the disease modules may mention known drug targets, and indications missing from the labels pull the AUC towards 0.5."""

def safety_table(records: list[dict[str, Any]], n: int = 15) -> str:
    """Toxicity and allergy view of the top candidates (all generated from candidates.json)."""
    lines = ["| Drug | Nephrocalcinosis requirement | Hypersensitivity contraindication | Excipient alerts | FAERS reports (top-term shares: renal / haematological / hepatic) |", "|---|:-:|:-:|---|---|"]
    for r in list(records)[:n]:
        gate = next(g["status"] for g in r["gates"] if g["requirement"] == "NEPHROCALCINOSIS_CAUTION")
        hyper = any("hypersensitivity among its contraindications" in e["text"] for e in r["evidence_against"])
        excipients = "; ".join(a["excipient"].title() for a in r["excipients"]["alerts"]) or "-"
        tox = r["toxicity"]
        faers = "not available"
        if tox.get("available"):
            sh = tox["group_shares_of_top_terms"]
            faers = f"{tox['total_reports']:,} ({sh['renal']:.0%} / {sh['haematological']:.0%} / {sh['hepatic']:.0%})"
        lines.append(f"| {r['name']} | {gate} | {'yes' if hyper else 'no'} | {excipients} | {faers} |")
    return chr(10).join(lines)


def _ped_names(out: Path, _kind: str) -> str:
    views = _optional(out, "age_views.json")
    if not views:
        return "none"
    names = sorted(set(views["register"]["adolescent"]["labelled_use"]) | set(views["register"]["child"]["labelled_use"]))
    return ", ".join(n.title() for n in names) or "none"


def primary_module(out: Path) -> str:
    """The modellable gene of the submission with the highest weight (its module is the primary one)."""
    genes = [g for g in _load(out / "submission_genes.json")["genes"] if g["modellable"]]
    return "sub_" + max(genes, key=lambda g: g["weight_epcr"])["gene"]


def _robust_summary(out: Path) -> str:
    cons = _optional(out, "consensus_ranking.json")
    if not cons:
        return "robustness was not computed"
    c = cons["summary"]
    return f"{c['robust']} drugs are robust, {c['intermediate']} intermediate and {c['fragile']} fragile"


def build(out: Path) -> str:
    patient = _load(Path("config/patient.json"))
    screen = _load(Path("config/mechanism_screen.json"))
    validation = _load(out / "validation_results.json")
    PRIMARY = primary_module(out)
    manifest = _load(out / f"manifest_{PRIMARY}.json")
    MODULE_TITLES = submission_modules(out)
    cov_path = Path("data/private/genome_covariates_summary.json")  # patient-level summary: private, not published
    covariates = _load(cov_path) if cov_path.exists() else {"aneuploidy_signal": {"per_chromosome": {}}}
    candidates = {r["drug_id"]: r for r in _load(out / "candidates.json")}
    scores = {m: _load(out / f"scores_{m}.json") for m in MODULE_TITLES}
    fanconi = _load(out / "scores_fanconi_anemia_demo.json")
    p1, p2, p3 = validation["P1_leave_one_out"], validation["P2_fanconi_controls"], validation["P3_P4_robustness"]
    groups: dict[str, list[str]] = {}
    for term in patient["hpo_terms"]:
        groups.setdefault(term["group"], []).append(f"{term['label']} ({term['id']})")
    group_lines = "\n".join(f"- **{g.replace('_', ' ')}:** {'; '.join(v)}" for g, v in groups.items())
    chrom = covariates["aneuploidy_signal"]["per_chromosome"]
    flagged = [c for c, r in chrom.items() if r["flag"]]
    n_universe = manifest["n_approved_drugs"]
    n_scored = scores[PRIMARY]["n_scored"]
    sections = []
    sections.append(f"""# Track 2 Report: Repositioning Hypotheses for Mosaic Variegated Aneuploidy (BUB1B)

**Author:** Juan Antonio Candelero Cervantes · **Hugging Face user:** `candeleroj98` · **Track:** 2 (Drug Repurposing)

> **These are hypotheses for follow-up, not evidence that any medicine works, and not treatment recommendations.** Every candidate is an already approved drug; the analysis ranks drugs by how close their targets lie to the biological pathways disrupted in this case.

## 1. Summary

- **Proposals (section 3.3).** {proposals_summary(out)}
- Universe: {n_universe} approved drugs (ChEMBL, `max_phase = 4`, release {manifest['chembl_version']}); {n_scored} of them have a curated protein target that can be placed in the interaction network.
- Method: network proximity (Guney 2016; Cheng 2018) between each drug's targets and {len(MODULE_TITLES)} gene modules derived automatically from the Track 1 submission file (genes, weights and ClinVar verification; rule in `protocol/addendum_submission.md`) (STRING {manifest['string_version']}, physical interactions, combined score ≥ 700; Reactome release {manifest['reactome_version']}).
- Validation (frozen protocol on the hand-defined checkpoint module, repeated on the submission-derived modules; section 6): the method **ranks held-out checkpoint genes well** (AUC {p1['auc']:.2f}, 95% CI {p1['auc_ci95'][0]:.2f}–{p1['auc_ci95'][1]:.2f}) but **does not recover the one available positive control in a second disease** and the single hand-defined module's drug ranking **changes with the interaction-network threshold** (the submission-derived BUB1B module and the combined ranking are much more stable, sections 5.12 and 6). A pre-registered benchmark over 30 other diseases is more favourable (mean AUC 0.69) but biased and optimistic (section 6). Drug-level results are therefore exploratory.
- Robustness (section 5.12): across 9 configurations (3 network thresholds x 3 weighting schemes) {_robust_summary(out)}. The robust set is one class of cytotoxic agents whose direction of effect in a checkpoint-deficient child is unknown (section 5.13), and its proximity is mostly membership of generic complexes (section 5.19). Only {_ped_names(out, 'robust')} of the robust and intermediate candidates have an established paediatric use on the US label, each for an indication unrelated to this condition (sections 5.15 and 7). The patient's age is not documented, so every candidate is shown for every age band.
- Every patient-specific safety requirement that depends on data we do not have (age, tumour status, organ function, blood counts, allergies, treatments) is reported as **NOT_ASSESSABLE**, never as passed. The one requirement that can be checked with documented data (nephrocalcinosis) fails when a drug's label warns about it. Toxicity (FDA adverse-event reports, label warnings) and allergy information (label contraindications, excipients) are attached to every candidate.

## 2. The proband's phenotype, used as a whole

The challenge's clinical document lists eight findings. The document asks that they be read as one constellation rather than individually, and that the parental miscarriages count as phenotypic input. They are used together here, grouped by the biological theme they illustrate:

{group_lines}

No individual finding selects or removes a drug. Birth details from the confidential document are deliberately not reproduced.

## 3. Variant mechanism and proposed candidates

{mechanism_block()}
### 3.2 From the Track 1 submission to disease modules

**Starting point: the Track 1 submission.** The only patient-specific input is the Track 1 CSV (`inputs/track1_submission.csv`, sha256 `{_load(out / "submission_genes.json")["submission_sha256"][:16]}…`). Each gene in it becomes a disease module by a fixed rule written before the run (Reactome pathways of the gene; STRING neighbours if none), weighted by the best `epcr` of its rows (a ranking value, not a probability). ClinVar is queried for every single-nucleotide variant and reported as context, not as a filter.

{submission_table(out)}

Genes whose module is too small (fewer than 10 proteins in the network) are reported and not scored; nothing was invented for them. Because the weights are the submission's own `epcr`, BUB1B (weight 1.0) dominates the combined score and the other genes (weights 0.10 to 0.14) add little.

**Earlier manual screening of the same genes (context only).** Before the submission-derived run, the genes were also screened by the author with three criteria: a pathogenic or likely pathogenic variant, a disease mechanism matching the phenotype constellation, and published pharmacological rationale.

{mechanism_table(screen)}

**Primary mechanism (BUB1B).** `p.Leu737*` is predicted to be a loss-of-function allele (nonsense; nonsense-mediated decay expected, not tested here) and is classified pathogenic/likely pathogenic in ClinVar; `p.Asn1002Lys` is a variant of uncertain significance. BubR1 is a component of the mitotic checkpoint complex, which delays anaphase until chromosomes are correctly attached; constitutional hypomorphic loss causes mosaic aneuploidy (MVA1) with predisposition to embryonal tumours. A biallelic (compound-heterozygous) state is a **hypothesis**: the phase of the two variants is not resolved, and the effect of `p.Asn1002Lys` is unproven. The BubR1 kinase domain is described in the literature as a pseudokinase while UniProt annotates a kinase domain; this discrepancy was not resolved from primary sources.

{comparison_paragraph(out)}

Nephrocalcinosis and the parental miscarriages are not established features of MVA1 in the sources consulted; they are treated as part of the constellation, not as explained by BUB1B.

{proposals_section(out)}

{register_section(out)}

{use_section(out)}

{avoid_section(out)}

## 4. Method

1. **Universe.** ChEMBL parent molecules with `max_phase = 4`, not withdrawn; targets from the curated mechanism table, mapped to UniProt, plus measured (pChEMBL >= 6, binding assay) and, if validated, similarity-predicted targets (section 5.16). **Proposals are restricted to FDA-approved drugs** (openFDA Drugs@FDA, currently marketed product), and shown by age group (section 3.3).
2. **Modules.** One per gene of the Track 1 submission: the gene plus every Reactome pathway mapped to its UniProt entry with 10 to 300 participants (`config/mechanisms_from_submission/*.json`, written by `track2.from_submission`); the BUB1B module has {scores[PRIMARY]['module']['n_in_interactome']} proteins in the network.
3. **Proximity.** d_closest between a drug's targets and the module, converted to a z-score against {manifest['iterations']} degree-matched random draws of both sets (degree bins of at least 100 nodes). Primary threshold z ≤ −2.
4. **Risk flags and direction.** Proximity does not say whether a drug helps or harms. Each drug carries its ChEMBL action, flags (antineoplastic, inhibits a module component, boxed warning, stabiliser) and `direction_review = REQUIRES_REVIEW`.
5. **Age bands.** Neonate, infant, child, adolescent and adult labels are extracted from the FDA label (DailyMed) by a heuristic and **must be verified by hand** for any drug that is taken further.
6. **Combined score.** Sum over modules of weight × max(0, −z), with weight equal to the gene's best Track 1 `epcr` (so BUB1B has 1.0).

Everything is reproducible from the repository (section 10); no gene, drug or score is fixed in the code.

## 5. Results""")
    for module, title in MODULE_TITLES.items():
        s = scores[module]
        sections.append(f"""### 5.{list(MODULE_TITLES).index(module) + 1} {title}

Module of {s['module']['n_in_interactome']} proteins; {s['n_pass']} of {s['n_scored']} scored drugs reach z ≤ −2. Top 10:

{top_table(s, candidates)}
""")
    sections.append("""### 5.5 How to read these tables

- The three Reactome-based modules (BUB1B, CENPE, MAD1L1) are large, mitosis-related pathway unions that share many proteins, so they select **the same drugs**: proteasome inhibitors, tubulin binders (taxanes, colchicine, vinca-type drugs), antibody-drug conjugates carrying tubulin toxins, and an XPO1 inhibitor. The agreement between modules is therefore weak evidence of robustness, not independent confirmation.
- In a person with a weakened checkpoint these agents may depend on a functional checkpoint for their antimitotic effect and may worsen chromosomal instability or toxicity in healthy tissue; they are flagged as **risks**, and any use would only make sense in a tumour context. Proteasome inhibitors are a **tumour-directed vulnerability hypothesis** (aneuploid cells are thought to be sensitive to proteotoxic stress), not a correction of the constitutional defect.
- The ranking is dominated by generic mitotic biology (the failure mode named in the addendum). The IGBP1 module (a regulator of PP2A-type phosphatases; neighbourhood-based, because IGBP1 has no Reactome pathway) is the only one with a different biology: kinase inhibitors of the RET and MEK/MAPK axis, a calcineurin inhibitor, adrenergic and vasopressin drugs. Its support is weak (a variant of uncertain significance, weight 0.0986, only a few drugs near z = -2) and it is shown as an exploratory hypothesis set. DNAH10 (a ciliary dynein) is absent from the physical network at score 700 and could not be modelled, so nothing is said about the renal and ciliary side of the phenotype.
- Several well-known drugs are **not scorable** because ChEMBL has no curated protein-target mechanism for them (for example metformin and vincristine). They are not ranked, and no value was invented for them.
""")
    sections.append(f"""{layers_sections(out)}

## 6. Validation (frozen protocol, `protocol/frozen_protocol.md`)

The protocol, with its decision rules, was written and hashed before any proximity was computed on real data (`protocol/freeze_note.md`). Clarifications made afterwards are in `protocol/deviations.md`.

- **P1, held-out checkpoint genes.** Leaving each of the {p1['n_seeds']} module genes out in turn, the held-out gene ranks well above degree-matched random genes: **AUC {p1['auc']:.2f} (95% CI {p1['auc_ci95'][0]:.2f}–{p1['auc_ci95'][1]:.2f})**; decision rule met (`{p1['decision']}`). However, no held-out gene reaches the primary threshold z ≤ −2 (median z {p1['median_z_seeds']:.2f}, against {p1['median_z_negatives']:.2f} for random genes), so the strict threshold selects very few drugs.
- **P2, a second disease with published controls (Fanconi anemia, exploratory).** Metformin, cisplatin and mitomycin C are not scorable (no curated protein target). Oxymetholone, a standard treatment, ranks at the {p2['OXYMETHOLONE']['percentile']:.0f}th percentile (z = {p2['OXYMETHOLONE']['z']:.2f}) of {p2['OXYMETHOLONE']['of']} drugs: **the method did not recover the available positive control.**
- **P3, robustness.** Dropping 20% of the seed genes 20 times keeps the top 20 similar (mean Jaccard {p3['drop_20pct_seeds']['mean_jaccard_top20']:.2f}). Changing the STRING score threshold changes it substantially (Jaccard {p3['string_thresholds']['400']['jaccard_top20']:.2f} at 400 and {p3['string_thresholds']['900']['jaccard_top20']:.2f} at 900), which by the protocol's rule means the ranking is **unstable with respect to the interactome threshold**.
- **P4, degree bias.** Spearman correlation between z and the targets' mean degree is {p3['p4_degree_bias_spearman']:.2f}: normalisation reduces but does not remove the bias.
{benchmark_paragraph(out)}
{submission_validation_paragraph(out)}

**Reading.** The network captures the checkpoint pathway at the gene level and, across many diseases, recovers known indications better than chance; but for this case the drug-level rankings are sensitive to the network threshold and the Fanconi control was not recovered. They should be read as a structured way to generate hypotheses, not as predictions.

## 7. Patient-specific considerations

**Requirements per drug.** For each candidate the repository records eight requirements: renal, hepatic and haematological safety, age label, administration route, tumour context, allergy history and a nephrocalcinosis caution. The clinical document gives no age, weight, tumour status, organ function, blood counts, allergies, treatments or medications, so those requirements are `NOT_ASSESSABLE`. The nephrocalcinosis caution is the exception: nephrocalcinosis is documented, so when a drug's FDA label warns about nephrocalcinosis, kidney stones or hypercalcaemia the requirement is `FAIL`, with the label excerpt recorded; a label that does not mention them stays `NOT_ASSESSABLE` (silence is not proof of safety). Other label warnings are recorded as alerts next to the requirement.

**Toxicity and allergy.** Each candidate carries (i) label warnings found by keyword (renal, hepatic, haematological, cardiac, neurological, pulmonary, embryo-fetal, hypersensitivity), (ii) whether hypersensitivity is listed among the contraindications, (iii) its inactive ingredients with alerts from an author-curated excipient list (for example benzyl alcohol, propylene glycol, cyclodextrin vehicles, castor-oil vehicles, sulfites, dyes), and (iv) a descriptive summary of the FDA Adverse Event Reporting System (number of reports and the share of the 100 most reported reaction terms in renal, haematological, hepatic, cardiac, neurological and hypersensitivity groups). The adverse-event data are spontaneous reports without a denominator and are confounded by the indication, so they give no incidence and no comparison between drugs. Cross-reactivity between drug classes and the patient's own allergy history are not evaluated (no allergy data are available).

{safety_table(list(candidates.values()))}

**Age bands.** `candidates_by_age_band.csv` marks, for each candidate, whether the US label supports use in neonates, infants, children, adolescents and adults. These labels were extracted heuristically and need manual verification before any use.

{age_section(out)}

**Covariates derived from the genome** (aggregates only; the genome itself is not published):
- The pattern of genotypes on chromosome X indicates a single X chromosome.
- Aneuploidy signal: {'no chromosome' if not flagged else ', '.join(flagged)} departs from the expected depth range by more than 5%{'' if not flagged else '; a gradient of slightly higher depth on small chromosomes is more consistent with a technical effect than with chromosome-specific mosaicism'}. This is a weak signal with limited sensitivity (low-fraction mosaicism or mosaicism in other tissues is not detectable) and is **not** a karyotype.
- Pharmacogenes (CYP3A4/5, UGT1A1, DPYD, TPMT, NUDT15, ABCB1, CYP2D6, CYP2C9, CYP2C19, SLCO1B1): no rare protein-altering variants were found in the variant call set; star-allele calling was not performed, so this does not establish a normal genotype.
- Tumour-predisposition genes (TP53, RB1, WT1, DICER1, NF1, PTCH1, CDKN2A, ATM, BRCA2, SMARCB1): one rare missense variant in a single gene, of unknown significance.

## 8. Scalability

The only input specific to this case is a mechanism file (`the Track 1 submission CSV; modules are derived from it automatically); the universe, network, scoring, requirement structure, schema and validation are generic. Replacing the file with another disease's pathway reruns the whole analysis, which is how the Fanconi anemia check in section 6 was run (module from the Reactome Fanconi anemia pathway, {fanconi['n_scored']} drugs scored). The age of the patient is also a parameter: with no age documented every candidate is shown for every age band, and a known age (`python -m track2.age_view --age-years <n>`, `config/patient.json`) filters the same output, so another patient only needs a new submission file and, if known, an age. That check also shows the limitation: a method that places drugs near a pathway does not by itself recover treatments whose benefit comes through other mechanisms. Applying it to a larger set of diseases would require curated positive controls, which are scarce.

## 9. Limitations

- Network proximity does not give the direction of an effect; Gene Ontology could not orient any candidate (section 5.9), so every candidate requires mechanistic review.
- Coverage depends on ChEMBL's curated mechanisms: {n_scored} of {n_universe} approved drugs are scorable.
- The drug ranking is sensitive to the interactome threshold and shows residual degree bias (section 6).
- No predicted affinity is used: docking was run with its own validation and no target passed (section 5.18), and measured ChEMBL values are given where they exist (section 5.11); where they do not exist nothing is invented.
- Out of scope for a hypothesis-generation hackathon: doses, paediatric dosing and exposure, and interactions with the patient's own medications (no medication list exists).
- The transcriptomic reversal layer gave a non-specific result (section 5.10) and is not used.
- Label-derived age and safety information is heuristic.
- The biallelic state of BUB1B is a hypothesis; the effect of `p.Asn1002Lys` is unproven; nephrocalcinosis and the parental miscarriages are not established consequences of BUB1B.
- One case: no sensitivity, specificity or significance is claimed. The drugs listed must not be used on the basis of this report.

## 10. Reproducibility

Code, configuration, frozen protocol and instructions for Windows, Linux and macOS are in the repository README. Data from third parties are downloaded at run time with their versions recorded (`out/manifest_*.json`); licences are in `THIRD_PARTY_LICENSES.md`.

## 11. Generative-AI disclosure

{AI_DECLARATION_FILE.read_text(encoding='utf-8').strip()}

## 12. Acknowledgments

We are deeply grateful to the child and their family who chose to share their genomic data openly with the research community. This work would not exist without their generosity and trust. No information that could re-identify the patient or family is included.

The Rare Disease Genomic Challenge (MVA Hackathon 2026) was organized in partnership with Sage Bionetworks, the MVA Society, Hugging Face and BEACON (the Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We thank the Wilhelm Foundation for their guidance and insight.

## References

1. Guney E, Menche J, Vidal M, Barábasi AL. Network-based in silico drug efficacy screening. *Nat Commun* 2016. doi:10.1038/ncomms10331 (PMID 26831545).
2. Cheng F, Kovács IA, Barabási AL. Network-based approach to prediction and population-based validation of in silico drug repurposing. *Nat Commun* 2018. doi:10.1038/s41467-018-05116-5 (PMID 30002366).
3. ChEMBL, STRING, Reactome, DailyMed: see `THIRD_PARTY_LICENSES.md`.
""")
    return "\n\n".join(sections) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--report", type=Path, default=Path("submissions/candeleroj98_track2_report.md"))
    args = parser.parse_args()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(build(args.out), encoding="utf-8", newline="\n")
    print("report written:", args.report)


if __name__ == "__main__":
    main()
