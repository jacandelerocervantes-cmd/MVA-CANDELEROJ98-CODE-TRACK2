"""Report sections: the variant's mechanism as the basis of the rationale, and the proposed candidates (all generated from the outputs)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from track2.age_groups import GROUPS

GROUP_TITLES = {"0_to_under_3": "0-<3 y", "3_to_under_12": "3-<12 y", "12_and_older": "12 y and older"}
STATUS_SHORT = {"labelled": "labelled", "not_labelled": "not labelled", "adult_only": "adult label only", "verify": "verify", "age_not_stated": "labelled (age not stated)"}


def mechanism_block() -> str:
    return """### 3.1 Characterisation of the variant's mechanism (basis of the repurposing rationale)

| Variant (GRCh38, MANE transcript NM_001211.6) | Molecular effect | Evidence | Status |
|---|---|---|---|
| *BUB1B* p.Leu737* (nonsense) | **Loss of function**: premature stop codon; nonsense-mediated decay expected (not tested here) | ClinVar pathogenic / likely pathogenic, multiple submitters, no conflicts | LoF allele |
| *BUB1B* p.Asn1002Lys (missense) | Working hypothesis: **partial loss of function (hypomorphic)** as the second allele; no evidence of gain of function | Variant of uncertain significance; REVEL 0.47 | Effect unproven; phase of the two alleles unresolved |

**Pathway disrupted.** BubR1 (encoded by *BUB1B*) is a component of the mitotic checkpoint complex (with BUB3, MAD2 and CDC20) of the **spindle assembly checkpoint**, which keeps the anaphase-promoting complex inactive until every kinetochore is correctly attached. With less functional BubR1, anaphase can start prematurely.

**Downstream biological consequences** (what a repurposed drug would have to act on), with the strength of the evidence for each:

1. **Chromosome missegregation and mosaic aneuploidy** (the defining feature of MVA; the proband's genome shows only a weak aneuploidy signal, which is not a karyotype).
2. **Proteotoxic and energy stress in aneuploid cells**: aneuploidy causes proteotoxic stress in yeast (PMID 23222101), and compounds acting on protein folding, autophagy or energy stress were selectively toxic to aneuploid cells (PMID 21315436). Cancer cells and yeast, not MVA patients.
3. **p53 / senescence response and progeroid features**: senescence and ageing in BubR1-insufficient mice depend on p16Ink4a and p19Arf (PMID 18516091), and genetic clearance of p16-positive senescent cells delayed age-related disorders (PMID 22048312). Mouse genetics.
4. **Growth restriction and muscle atrophy**: in mice with BubR1 mutations, a predisposition to sarcopenia correlated with mTORC1 hyperactivity (PMID 31738183). The proband has short stature, failure to thrive and skeletal muscle atrophy. Mouse correlation.
5. **Tumour predisposition**: rhabdomyosarcoma and Wilms tumours have been reported in children with this syndrome (PMIDs 31184400, 31081598). Case reports.

**How the mechanism shapes the rationale.** No approved drug restores BubR1. The pipeline therefore asks three questions, each tied to a consequence above: (a) which approved drugs act on proteins of the checkpoint and mitotic pathways (network proximity to modules built from the genes of the submission), (b) which exploit the vulnerability of aneuploid cells in a **tumour context** (consequences 2 and 5), and (c) which act on downstream **consequences** such as the mTOR and senescence axes (consequences 3 and 4), kept as a separate extrapolation arm because the evidence is from mice. Direction of effect is not established for any candidate.
"""


def _row(r: dict[str, Any]) -> str:
    ev = r["evidence"]
    basis = (f"proximity {ev['robustness_count_of_9']}/9 configs, z {ev['expanded_best_z']}; binding {'yes' if ev['L2_binding'] else 'no'}; rationale {'yes' if ev['L3_rationale'] else 'no'}"
             if r["arm"] == "main" else "extrapolation (mouse evidence): " + ", ".join(c.replace("_", " ") for c in ev["L3_classes"]))
    l4 = ev.get("L4") or {}
    l4_text = {"supports_resistance_when_SAC_weak": "weak-SAC lines MORE resistant", "supports_sensitivity_when_SAC_weak": "weak-SAC lines more sensitive",
               "no_association": "no association", "lineage_confounded": "lineage-confounded", "mixed": "mixed", "not_in_PRISM": "not in PRISM", "not_assessed": "not assessed"}.get(l4.get("call"), "-")
    if l4.get("diff_weak_minus_strong") is not None and l4.get("call") not in ("no_association", "not_in_PRISM"):
        l4_text += f" ({l4['diff_weak_minus_strong']:+.2f})"
    core_label = (ev.get("core") or {}).get("label")
    if r["arm"] == "main" and core_label == "proximity_by_generic_complex_only":
        basis += "; **proximity only through generic complexes**"
    elif r["arm"] == "main" and core_label == "proximity_also_without_generic_complexes":
        basis += "; proximity also without generic complexes"
    basis = basis + f"; L4 PRISM: {l4_text}"
    acts = "; ".join(x["role"] for x in r["acts_on"][:2])
    age = " / ".join(STATUS_SHORT.get(r["age_groups"][g], r["age_groups"][g]) for g in GROUPS)
    c = r["cautions"]
    cautions = ", ".join(x for x in ("boxed warning" if c["boxed_warning"] else "", c.get("label_caution") or "", "broad target profile" if c["broad_target_profile"] else "",
                                     ("excipients: " + "/".join(c["excipient_alerts"])) if c["excipient_alerts"] else "") if x) or "-"
    tier = f"T{r['tier']}" if r["tier"] else "ext."
    return f"| {tier} | {r['name'].title()} | {'; '.join(r['approved_for'][:3]).lower() or '-'} | {acts[:150]} | {basis} | {age} | {cautions} |"


def _table(rows: list[dict[str, Any]]) -> str:
    head = ("| Tier | Drug | Approved today for | Acts on | Basis | Age " + " / ".join(GROUP_TITLES[g] for g in GROUPS)
            + " (label) | Cautions shown, not used to exclude |\n|:-:|---|---|---|---|---|---|")
    return head + "\n" + "\n".join(_row(r) for r in rows)


def proposals_section(out: Path) -> str:
    path = out / "proposals.json"
    if not path.exists():
        return ""
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = data["proposals"]
    plaus = json.loads((out / f"plausibility_{data['variant']}.json").read_text(encoding="utf-8"))["summary"]
    tumour = [r for r in rows if r["arm"] == "main" and r["list"] == "tumour_context"]
    constit = [r for r in rows if r["arm"] == "main" and r["list"] == "constitutional_defect"]
    extra = [r for r in rows if r["arm"] == "extrapolation"]
    t1 = [r["name"].title() for r in rows if r["arm"] == "main" and r["tier"] == 1]
    lab = {g: sorted({r["name"].title() for r in rows if r["age_groups"][g] in ("labelled",) or r["age_groups"][g].startswith("partly")}) for g in GROUPS}
    verify_names = ", ".join(sorted({r["name"].title() for r in rows if "verify" in r["age_groups"].values()})) or "none"
    return f"""### 3.3 Proposed repositioning candidates (FDA-approved drugs only; hypotheses for follow-up, not medical advice)

**Scope and selection rule** (`protocol/addendum_plausibility.md`, fixed before the run). Only drugs with an openFDA Drugs@FDA record and a currently marketed (Prescription) product are eligible ({plaus['n_not_fda_eligible']} drugs that met the proximity criteria were removed because they are discontinued or not FDA-approved). "Possible to reuse" is not one score: it is the number of independent evidence lines that agree, plus feasibility gates. **L1** the drug is in the best 20 in at least 7 of 9 configurations (3 network thresholds x 3 weighting schemes); **L2** it is linked to a module protein by its curated mechanism or by a measured activity (pChEMBL >= 6); **L3** a PMID-verified literature rationale applies; **L4** independent functional data (cell-line dependency or drug-sensitivity screens) are **not assessed in this version**. **Tier 1** = L1+L2+L3, **Tier 2** = L1+L2 (gates: FDA-eligible, not topical-only, not an acute-use anaesthetic or sedative). Drugs for tumours are listed apart from drugs for the constitutional defect, and each row shows what the drug is approved for today and which gene products it acts on, taken from ChEMBL.

**Age.** The patient's age is not documented, so each candidate is shown for three age groups (0-<3, 3-<12, 12 and older), grouped from the five label bands used in the methods (neonate, infant, child, adolescent, adult; `out/age_views.csv`). The label status is `labelled` (use supported across the group), `partly (from n y)`, `not labelled`, `adult label only` or `verify` (heuristic extraction, reviewed by hand where marked; verify against the current label). With a known age the same tables filter (`python -m track2.age_view --age-years N`). This is label support, not a dose.

**A. Tumour-context candidates** (relevant only if a tumour is present; the proband has a tumour history). {len(tumour)} drugs; Tier 1: {', '.join(t1) or 'none'}.

{_table(tumour)}

**B. Candidates for the constitutional defect and its consequences** (non-antineoplastic). {len(constit)} drugs reach Tier 1 or 2 and none reaches Tier 1: the network and binding evidence does not support a robust non-cytotoxic candidate for the constitutional defect.

{_table(constit)}

**C. Extrapolation arm (growth, muscle and senescence consequences; not part of the ranking).** {len(extra)} FDA-approved drugs whose mechanism matches the mTOR or senescence rationale from the mouse studies cited in section 3.1; no drug was tested in those mice, and the proximity of this arm is partly circular because its module was seeded with pathway genes.

{_table(extra)}

**What the table does and does not say.** Label support by age group for the candidates: {'; '.join(f"{GROUP_TITLES[g]}: {', '.join(lab[g]) or 'none'}" for g in GROUPS)}. Every other candidate is not labelled for children (adult label only) or needs the label verified by hand ({verify_names}; their paediatric section could not be extracted). Direction of effect (whether a drug would help or worsen chromosomal instability) is not established for any row; the curated literature notes in section 5.13 show where evidence exists and that it is mostly from cancer cells or mice. A candidate is a hypothesis for specialist review, never a recommendation.
"""


def dti_section(out: Path) -> str:
    e1 = json.loads((out / "dti_e1_module_edges.json").read_text(encoding="utf-8"))["summary"] if (out / "dti_e1_module_edges.json").exists() else None
    val = json.loads((out / "dti_e2_validation.json").read_text(encoding="utf-8")) if (out / "dti_e2_validation.json").exists() else None
    path = out / "plausibility_e0_e1_e2.json"
    if not (e1 and val and path.exists()):
        return ""
    plaus = json.loads(path.read_text(encoding="utf-8"))
    row = next(r for r in val["by_tau"] if r["tau"] == val["tau"]) if val["valid"] else None
    new = [r for r in plaus["drugs"] if r["novelty"].startswith("new")]
    lines = ["| Drug | Scorable before | Curated best z | Expanded best z | Closest targets (gene, source, steps from a module) | FDA | Tier |", "|---|:-:|:-:|:-:|---|:-:|:-:|"]
    for r in sorted(new, key=lambda x: (x["tier"] == 0, x["tier"], x["expanded_best_z"])):
        prox = "; ".join(f"{p['gene']} ({p['source']}{', pChEMBL ' + str(p['median_pchembl']) if p.get('median_pchembl') else ''}; {p['distance']})" for p in r["proximal_targets"][:3])
        lines.append(f"| {r['name'].title()} | {'yes' if r['curated_best_z'] is not None else 'no'} | {r['curated_best_z'] if r['curated_best_z'] is not None else '-'} | {r['expanded_best_z']} | {prox or '-'} | {'yes' if r['gates']['fda_eligible'] else 'no'} | {r['tier'] or '-'} |")
    ho = val.get('heldout_by_target') or {'heldout_precision_mean': float('nan'), 'heldout_precision_p10_p90': [float('nan'), float('nan')]}
    e2_text = (f"E2 passed: at tau = {val['tau']} the precision was {row['precision']:.2f} and the recall {row['recall']:.2f} ({row['true_positive']} true links of {row['predicted']} predicted), against a base rate of {val['base_rate']:.4f}; "
               f"it is a weak signal that met the pre-set rule in-sample. **Out of sample it does not hold:** when the targets are split in two halves and the cut-off chosen on one half is applied to the other (20 repeats), the precision is {ho['heldout_precision_mean']:.2f} (10th to 90th percentile {ho['heldout_precision_p10_p90'][0]:.2f} to {ho['heldout_precision_p10_p90'][1]:.2f}), below the 0.5 criterion. No drug of the register depends on E2, and its links are labelled as predicted and never count as evidence." if row else "E2 did not reach the pre-set precision and was not used.")
    return f"""### 5.16 Expanded drug-target links: measured and predicted binding (`protocol/addendum_dti.md`)

The ranking so far used only each drug's curated mechanism, so it could not find a **new use**. Two kinds of links were added beyond the curated table. **E1, measured:** ChEMBL binding activities of approved small-molecule drugs on human single-protein targets with pChEMBL >= 6 (<= 1 uM): {e1['n_drugs_any_target']:,} of {e1['n_drugs_queried']:,} drugs have at least one measured target, and {e1['n_drugs_with_module_link']} are linked directly to a protein of the four modules ({e1['n_links']} links; the modules have {e1['n_module_nodes']} proteins, mostly kinetochore, APC/C and motor proteins that few drugs bind). **E2, predicted by ligand similarity** (Morgan fingerprints; each drug scored against reference ligands of {val['n_targets_with_ligands']} module proteins that have at least 10; the drug itself removed): the retrospective test over {val['n_pairs']:,} drug-protein pairs with {val['n_positive']} known links set the threshold. {e2_text} Docking was run with its own validation and did not pass (section 5.18).

With the expanded links {plaus['summary']['n_pool']} drugs are in the pool (new or in the top 20 of at least 5 configurations), and the number of drugs with a computable proximity grows from 1,040 to about 1,440. **{plaus['summary']['n_new_hypotheses']} drugs are new hypotheses** (close to a module only thanks to the expanded links, or not scorable before); {plaus['summary']['n_new_hypotheses'] - sum(1 for r in new if r['gates']['fda_eligible'])} of them fail the FDA scope (discontinued or not FDA-approved), and **{plaus['summary']['n_robust_new']} is robust** across the nine configurations. A first run of this step used an incomplete protein mapping (one accession per network node) that hid most links; it was corrected before the results were read (`protocol/deviations.md`). After an independent review, a measured link counts as evidence only if it rests on at least two activity records: every link that carried the new hypotheses of fostamatinib (NEK9, AURKA, PLK4), selumetinib (SMC1A, SMC2) and daunorubicin (HDAC1, HDAC6, HDAC8) was a single measurement, so those drugs no longer reach L2.

{chr(10).join(lines)}

Reading: the expanded links add almost nothing. The new items rest on single measurements (fostamatinib binding NEK9, PLK4 and AURKA; selumetinib binding SMC1A and SMC2; daunorubicin binding HDAC1, HDAC6 and HDAC8), so they are shown as exploratory and do not count as evidence. Fostamatinib is also a broad kinase binder whose label does not recommend its use under 18 years because of growth-plate effects in juvenile animals. The honest summary is that **measured and predicted binding beyond the curated mechanism does not reveal a strong non-cytotoxic repositioning candidate for this case**.
"""


def proposals_summary(out: Path) -> str:
    path = out / "proposals.json"
    if not path.exists():
        return "not computed"
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = data["proposals"]
    main = [r for r in rows if r["arm"] == "main"]
    t1 = [r["name"].title() for r in main if r["tier"] == 1]
    constit = [r["name"].title() for r in main if r["list"] == "constitutional_defect"]
    extra = [r["name"].title() for r in rows if r["arm"] == "extrapolation"]
    labelled = sorted({r["name"].title() for r in rows if any(v == "labelled" or v.startswith("partly") for k, v in r["age_groups"].items() if k != "12_and_older")})
    return (f"{len(main)} FDA-approved drugs reach Tier 1 or 2 of the plausibility framework (independent evidence lines plus feasibility gates fixed in advance): "
            f"Tier 1 = {', '.join(t1) or 'none'} (all for a tumour context); for the constitutional defect and its consequences only {', '.join(constit) or 'none'} reach Tier 2 and none reaches Tier 1, "
            f"so **no robust non-cytotoxic candidate is supported**; an extrapolation arm from mouse evidence adds {', '.join(extra) or 'none'} (mTOR pathway). "
            f"Candidates are shown for three age groups (0-<3, 3-<12, 12 and older) because the child's age is unknown; with a paediatric label for at least one group under 12: {', '.join(labelled) or 'none'}. "
            "**The network proximity of the Tier 1 and Tier 2 drugs is mostly membership of generic complexes** (proteasome, tubulin; section 5.19): without them almost none stays close, so their rationale rests on the literature and on the cell-line direction, not on the network. A complete justified register (section 3.4) lists every FDA-approved drug that passed the gates and reached at least one evidence line, with how each could be validated or refuted (section 3.5). Independent cell-line data (section 5.17, adjusted for proliferation) show that the microtubule agents are less effective where the checkpoint is weak, and docking (section 5.18) did not pass its validation. The mechanism of the variant (loss of function of BUB1B, spindle assembly checkpoint) is characterised in section 3.1 as the basis of the rationale. Hypotheses for follow-up, not evidence that any medicine works.")


def prism_section(out: Path) -> str:
    path = out / "prism_l4.json"
    if not path.exists():
        return ""
    data = json.loads(path.read_text(encoding="utf-8"))
    main, alt = data["mean_all"], data["BUB1B"]
    label = {"supports_resistance_when_SAC_weak": "weak-checkpoint lines more resistant", "supports_sensitivity_when_SAC_weak": "weak-checkpoint lines more sensitive", "no_association": "no association",
             "confounded_by_proliferation_or_lineage": "explained by proliferation or lineage", "mixed": "mixed", "not_in_PRISM": "not in PRISM", "not_testable": "not testable"}
    lines = ["| Drug | Call (mean score) | Viability difference weak minus strong | FDR unadjusted | FDR adjusted for lineage and proliferation | Adjusted effect per SD of the score | Call (BUB1B alone) |", "|---|---|:-:|:-:|:-:|:-:|---|"]
    resistant, confounded, trend = [], [], []
    for did, d in sorted(main["drugs"].items(), key=lambda kv: (kv[1]["call"] == "not_in_PRISM", kv[1]["name"])):
        if d["call"] == "not_in_PRISM":
            continue
        best = min(d["rows"], key=lambda r: r["fdr"]) if d["rows"] else None
        if d["call"] == "supports_resistance_when_SAC_weak":
            resistant.append(d["name"].title())
        elif d["call"] == "confounded_by_proliferation_or_lineage":
            confounded.append(d["name"].title())
        elif best and best["fdr_lineage_prolif"] < 0.10 and best["coef_lineage_prolif_per_sd"] < 0:
            trend.append(f"{d['name'].title()} (adjusted FDR {best['fdr_lineage_prolif']:.3g})")
        if best:
            lines.append(f"| {d['name'].title()} | {label.get(d['call'], d['call'])} | {best['diff_weak_minus_strong']:+.2f} | {best['fdr']:.3g} | {best['fdr_lineage_prolif']:.3g} | {best['coef_lineage_prolif_per_sd']:+.2f} | "
                         f"{label.get(alt['drugs'][did]['call'], alt['drugs'][did]['call'])} |")
    return f"""### 5.17 Independent functional evidence (L4): drug sensitivity of cancer cell lines with a weaker spindle checkpoint (`protocol/addendum_prism.md`)

Do the candidates work better or worse when the spindle assembly checkpoint is weak? PRISM Repurposing Public 24Q2 (single-dose pooled screen of up to {main['n_cell_lines']} cancer cell lines) was joined to DepMap Public 24Q2 expression (both CC BY 4.0, downloaded from figshare with MD5 checks). Each cell line gets a SAC expression score (mean z-score of BUB1B, BUB1, BUB3, MAD1L1, MAD2L1 and TTK) and the {main['n_weak']} lowest-scoring lines are compared with the {main['n_strong']} highest-scoring ones. {main['n_in_prism']} of the {main['n_candidates']} candidates are in PRISM.

**A confounder that was found and corrected.** An independent review found, and it was verified, that the SAC score correlates with proliferation (Spearman {main['spearman_score_vs_proliferation']:.2f} with a score built from {', '.join(main['proliferation_genes'])}): checkpoint genes are proliferation genes, and antimitotic drugs act on proliferating cells, so the plain comparison cannot separate a weak checkpoint from slow growth. The analysis was therefore amended (`protocol/addendum_prism.md`, amendment 1): the call needs significance both unadjusted and in a model adjusted for lineage and proliferation, with the same sign.

**Result.** No candidate is more sensitive in weak-checkpoint cells. **{len(resistant)} are more resistant and the result survives the adjustment:** {', '.join(resistant) or 'none'}; they are all microtubule-targeting agents, in line with taxane-resistance reports (PMIDs 17234259 and 19458043). {('Explained by proliferation or lineage once adjusted: ' + ', '.join(confounded) + '. ') if confounded else ''}{('Proteasome inhibitors and others show no association unadjusted, but after adjusting the direction is the same (less effective with a weak checkpoint): ' + '; '.join(trend) + '. This is shown and not counted as support, because it is not confirmed unadjusted.') if trend else ''} The microtubule agents are therefore expected to be less, not more, effective where the checkpoint is weak. Limits: tumour cell lines, expression as a proxy for checkpoint function, a single concentration, and a correlation across lines; this is a hypothesis about direction and not a clinical statement.

{chr(10).join(lines)}
"""


def docking_section(out: Path) -> str:
    path = out / "docking_validation.json"
    struct = Path("data/docking/structures.json")
    if not path.exists():
        return ""
    val = json.loads(path.read_text(encoding="utf-8"))
    rec = json.loads(struct.read_text(encoding="utf-8")) if struct.exists() else {}
    rules = val["rules"]
    lines = ["| Target (UniProt) | Structure, ligand | Redocking RMSD (limit " + f"{rules['redock_rmsd_max']:g}" + " A) | Known approved binders | Enrichment AUC (min " + f"{rules['auc_min']:g}" + ") | Cut-off with precision >= 0.5 | Validated |",
             "|---|---|:-:|:-:|:-:|:-:|:-:|"]
    for name, e in val["targets"].items():
        sel = (rec.get(name) or {}).get("selected") or {}
        r = e["redock"]
        auc_text = f"{e['auc']:.2f}" if e.get("auc") is not None else "not computed (fewer than " + str(rules["min_actives"]) + " binders)"
        thr = f"{e['threshold_affinity']:.1f}" if e.get("threshold_affinity") is not None else "none"
        lines.append(f"| {name} ({(rec.get(name) or {}).get('uniprot', '')}) | {e['structure']}, {sel.get('ligand', '')} | {r['rmsd_best_pose']:.2f} {'(pass)' if e['passes_redock'] else '(fail)'} | {e['n_actives']} | {auc_text} | {thr} | {'yes' if e['validated'] else 'no'} |")
    n_valid = sum(1 for e in val["targets"].values() if e["validated"])
    return f"""### 5.18 Docking (E3), with its own validation (`protocol/addendum_docking.md`): a negative result

Docking was run with AutoDock Vina 1.2.7 (official Scripps build, SHA-256 recorded) on the best-resolution human X-ray structure with a drug-like ligand for seven kinases of the checkpoint and mitotic network (CDK1, AURKA, AURKB, PLK1, PLK4, TTK, NEK7), chosen by a rule fixed beforehand. Before any score could be used each target had to pass two checks: the co-crystal ligand must be re-docked within {rules['redock_rmsd_max']:g} angstrom, and the Vina affinity must rank known approved binders above {rules['n_decoys']} other approved drugs (ROC AUC >= {rules['auc_min']:g}, at least {rules['min_actives']} binders) with a usable cut-off.

{chr(10).join(lines)}

**{n_valid} of {len(val['targets'])} targets passed both checks, so no docking score was used and no docking-predicted link (E3) was created.** A symmetry-aware, best-of-nine-poses check of the RMSD did not rescue AURKA (7.4 angstrom) or PLK1 (4.3), so the failures are not an artefact of how the RMSD was measured; but five failures above 5 angstrom also point to the limits of this preparation (rigid receptor, ligand protonation, cofactors, search box), so the result says that **this protocol** cannot be used, not that docking cannot work. Where the enrichment looked good (AURKA, PLK4) the poses were wrong, so the enrichment cannot be trusted as evidence of binding; AURKB reproduces its crystal pose and ranks binders reasonably (AUC 0.81) but no cut-off reaches the precision needed to call hits; CDK1 reproduces its pose but has too few known approved binders to validate. Limits of this implementation: rigid receptor, ligands protonated by default rules, Vina scoring that favours large molecules, and no covalent docking (proteasome inhibitors). The conclusion is that **this docking protocol cannot be trusted to propose new binders**, and the repurposing hypotheses rest on measured binding, ligand-similarity prediction and network proximity instead. The code, the structures and the validation are in the repository so that a better-validated protocol (flexible receptor, other targets, consensus scoring) can be tried.
"""


def register_section(out: Path) -> str:
    path = out / "candidate_register.json"
    if not path.exists():
        return ""
    reg = json.loads(path.read_text(encoding="utf-8"))
    l4 = {"supports_resistance_when_SAC_weak": "more resistant", "supports_sensitivity_when_SAC_weak": "more sensitive", "no_association": "none", "lineage_confounded": "confounded",
          "mixed": "mixed", "not_in_PRISM": "not in PRISM", "not_assessed": "-", "not_testable": "-"}
    lines = ["| Drug | List | Tier strict / wide | L1 / L1 core / L2 / L3 | PRISM (weak checkpoint) | Approved today for | Age 0-<3 / 3-<12 / 12+ | Justification |", "|---|---|:-:|:-:|---|---|---|---|"]
    for r in reg["drugs"]:
        age = " / ".join(STATUS_SHORT.get(r["age_groups"][g], r["age_groups"][g]) for g in GROUPS)
        lines.append(f"| {r['name'].title()} | {r['list'].replace('_', ' ')} | {r['tier_strict'] or '-'} / {r['tier_wide'] or '-'} | {'Y' if r['L1'] else 'n'} / {'Y' if (r.get('core') or {}).get('L1_core') else 'n'} / {'Y' if r['L2'] else 'n'} / {'Y' if r['L3'] else 'n'} | "
                     f"{l4.get(r['L4_call'], r['L4_call'])} | {'; '.join((r.get('approved_indications') or [])[:2]).lower() or '-'} | {age} | {r['justification']} |")
    return f"""### 3.4 Complete justified register ({reg['n_drugs']} FDA-approved drugs)

Every FDA-approved drug that passed all the feasibility gates and reached at least one evidence line in any arm of the analysis (strict tiers: {reg['by_source']['strict']}; wider ring of addendum 6: {reg['by_source']['wide_ring']}; extrapolation arm: {reg['by_source']['extrapolation_arm']}) is listed here with the evidence it satisfies and the evidence it lacks. The tiers matter less than the justification: a reader can see for each drug why it is on the list and what is missing. The same table is `out/candidate_register.csv`. Drugs that met the proximity criteria but are not FDA-eligible, are topical-only or are acute-use anaesthetics and sedatives are excluded and counted in section 3.3.

{chr(10).join(lines)}
"""


def use_section(out: Path) -> str:
    path = out / "candidate_register.json"
    if not path.exists():
        return ""
    reg = json.loads(path.read_text(encoding="utf-8"))
    cfg = json.loads(Path("config/repositioning_rationale.json").read_text(encoding="utf-8"))
    cons = cfg["consequences"]
    lines = ["| Drug | Approved today for | Proposed repositioned use in this case | Consequence of the variant it addresses | Caveat on direction and evidence | How it could be validated or refuted |", "|---|---|---|---|---|---|"]
    groups: dict[str, list[str]] = {}
    for r in reg["drugs"]:
        u = r["repositioning"]
        c = "; ".join(f"{k} ({cons[k]})" for k in u["consequences"]) or "-"
        lines.append(f"| {r['name'].title()} | {'; '.join((r.get('approved_indications') or [])[:2]).lower() or '-'} | {u['new_use']} | {c} | {u['caveat'] or '-'} | {u.get('validation') or '-'} |")
        groups.setdefault(u["rule"], []).append(r["name"].title())
    counts = "; ".join(f"{k.replace('_', ' ')}: {', '.join(v)}" for k, v in groups.items())
    return f"""### 3.5 What each drug is proposed for

The tables above say what each drug is approved for today and which gene products it acts on. This section says **what each one is proposed for in this case**, tied to the consequences of the variant of section 3.1 (1 missegregation and mosaic aneuploidy; 2 proteotoxic and energy stress; 3 p53 and senescence response; 4 growth restriction and muscle atrophy; 5 tumour predisposition). The wording comes from fixed rules over each drug's own mechanism and measured links (`config/repositioning_rationale.json`, `src/track2/repositioning.py`), not from text written drug by drug, and it only says what the evidence in this repository supports. Grouped by rule: {counts}.

{chr(10).join(lines)}

In short: the proposals are **tumour-directed** (proteasome inhibitors as the best supported, antimitotic agents and antibody-drug conjugates as listed but not favoured), plus **one extrapolation for the growth, muscle and senescence consequences** (mTOR-pathway drugs), plus **two hypotheses of modulating the mitotic kinase network** by measured binding (fostamatinib and, more weakly, the new drugs of section 5.16). No drug is proposed to repair the checkpoint itself, because none is known to do so.
"""


def core_section(out: Path) -> str:
    path = out / "core_sensitivity.json"
    reg_path = out / "candidate_register.json"
    if not path.exists() or not reg_path.exists():
        return ""
    core = json.loads(path.read_text(encoding="utf-8"))
    reg = json.loads(reg_path.read_text(encoding="utf-8"))["drugs"]
    only = [r for r in reg if (r.get("core") or {}).get("label") == "proximity_by_generic_complex_only"]
    both = [r for r in reg if (r.get("core") or {}).get("label") == "proximity_also_without_generic_complexes"]
    others = [d for d in core["drugs"] if d["L1_core"] and d["drug_id"] not in {r["drug_id"] for r in reg}]
    return f"""### 5.19 Sensitivity without generic complexes: how much of the proximity is pathway membership

An independent review found, and it was verified, that the modules built from Reactome pathways contain **30 proteasome subunits plus ADRM1 and 24 tubulins** (in the BUB1B, CENPE and MAD1L1 modules), so a drug whose targets are those proteins has distance 0 to the module **by construction**. The strong proximity of proteasome inhibitors, taxanes, colchicine, vinca-type agents and the antibody-drug conjugates with tubulin toxins (z of about -7 to -10) is therefore pathway membership and not independent evidence that they act on the disease mechanism. The modules were rebuilt without those proteins (`protocol/addendum_core.md`, written after the review) and the nine-configuration ranking recomputed.

- Drugs of the register whose robust proximity **disappears** without the generic complexes ({len(only)}): {', '.join(r['name'].title() for r in only) or 'none'}.
- Drugs of the register that remain robust ({len(both)}): {', '.join(r['name'].title() for r in both) or 'none'}.
- Other drugs that are robust only with the core modules and are not in the register ({len(others)}): {', '.join(d['name'].title() for d in others) or 'none'}; these are broad kinase inhibitors or topical agents, so their closeness may reflect breadth of targets and they are not proposed.

**Consequence for the proposals.** The rationale for the proteasome inhibitors (the Tier 1 drugs) is the literature on the proteotoxic stress of aneuploid cells (section 3.1), not the network; for the microtubule agents and conjugates the network adds nothing beyond membership of tubulin. Network proximity, as built here, finds almost nothing specific to this case beyond the generic machinery that every pathway of the checkpoint shares.
"""


def avoid_section(out: Path) -> str:
    """Drugs to be wary of, generated from the data (cell-line direction, generic-complex proximity, label cautions)."""
    path = out / "candidate_register.json"
    if not path.exists():
        return ""
    reg = json.loads(path.read_text(encoding="utf-8"))["drugs"]
    less = [r["name"].title() for r in reg if r["L4_call"] == "supports_resistance_when_SAC_weak"]
    growth = [r["name"].title() for r in reg if "growth" in (r.get("label_caution") or "").lower() or "bones" in (r.get("label_caution") or "").lower()]
    return f"""### 3.6 What to avoid or treat with caution (generated from the data)

- **Microtubule agents are expected to be less effective, not more, where the checkpoint is weak:** {', '.join(less) or 'none'}. In cancer cell lines with a weaker checkpoint they were more resistant even after adjusting for lineage and proliferation (section 5.17), in line with taxane-resistance reports (PMIDs 17234259, 19458043). They should not be read as candidates for the constitutional defect.
- **Proximity that is only membership of generic complexes** (proteasome, tubulin; section 5.19) is not evidence of relevance to this case.
- **Single-measurement binding** (one activity record) is not counted as evidence (section 5.16).
- **Label cautions that matter for a child with growth restriction:** {', '.join(growth) or 'none in the register'} (for example the fostamatinib label does not recommend use under 18 years because of growth-plate effects in juvenile animals).
- Sedatives, anaesthetics and other acute-use drugs that reached the proximity criteria (GABA-A modulators) are network-topology artefacts and are excluded.
"""
