# Addendum 11: amendments after the second independent review (APPROVED by the owner on 2026-10-07 with two conditions, E14 and E15; nothing has been run under it yet)

Status: approved 2026-10-07. When approved it moves to `protocol/addendum_review2.md` **before** any analysis below is run. There is no git repository (owner's decision), so the order of events is attested by the author only; the report must say so.

## Principle

No threshold of addendum 10 is loosened. Every amendment below is stricter than or equal to the rule it replaces. The review that motivated them read results already produced, so all of them are **post hoc with respect to those results**; they are pre-specified with respect to the new runs. All earlier results stay in `deviations.md` and are reported as sensitivity analyses.

## Amendments to the lines

**E1. Line B is split.** B1 = proximity to the phenotype module (validated, usable). B2 = proximity to the mitosis module without generic complexes (not validated). B2 supports a candidate **only if** `module_rule_validation` passes (mean AUC >= 0.60, lower bound > 0.50, >= 15 evaluable diseases); otherwise it is reported as exploratory. B1 and B2 share the proximity method and the drug-target table, so together they count as **one** line "B" in the two-line rule, never two.

**E2. Validation of B1 without the causal gene.** The validation of amendment 1 is repeated removing from each disease module its causal genes (those annotated to its OMIM entries in `genes_to_phenotype.txt`) and removing from the ground truth every known drug whose curated or measured targets include a causal gene. Diseases left with fewer than 3 known drugs are dropped. Thresholds unchanged (mean AUC >= 0.60, lower bound > 0.50, >= 15 evaluable diseases). **B1 is usable only if this leave-causal-gene-out version passes**; the original version is reported next to it.

**E3. Excess over the degree null is a second condition for A and B1.** In addition to the existing thresholds, the bootstrap lower bound of (AUC - degree-matched null AUC) must be > 0. Line A already satisfies this (0.014 to 0.030); B1 is evaluated under E2.

**E4. Specificity of A for the patient.** A drug counts as supported by A only if it is in the real top 50 **and** in the top 50 of at most 5 % of the 500 degree-matched random queries of `kg_nullquery`. Drugs failing the second condition are listed as "generic" and do not count.

**E5. Patient phenotype set.** HP:0200067 (parental history of recurrent miscarriage) is removed from the queries of A and B1: it describes the parents, not the child. The queries use the remaining seven terms (six in PrimeKG, which has no node for rhabdomyosarcoma). Reported as sensitivity: (a) each term left out in turn, with the overlap (Jaccard) of the top 50 of A and of the robust set of B1; (b) only the core terms (rhabdomyosarcoma, nephrocalcinosis, skeletal muscle atrophy) against only the non-specific constitutional terms (short stature, failure to thrive, small for gestational age, premature birth). No claim of a phenotype signal is made for a class that does not survive (a).

**E6. Line C with TP53 and a null of 1,000 sets.** (i) The model adds TP53 status as a covariate (damaging somatic mutation in DepMap 24Q2 `OmicsSomaticMutations`). Where the files exist in the same release, a chromosomal-instability score and MSI status are added; if they do not exist they are reported as not available, not approximated. (ii) The variance inflation factor of the SAC score against proliferation is reported, together with a sensitivity model with the SAC score orthogonalised to proliferation. (iii) The null uses 1,000 matched sets and the same procedure; **line C is used only if the observed total in the TP53-adjusted model exceeds the 95th percentile of that null**, and the 1,000-set estimate replaces the 200-set one whichever way it falls. (iv) PRISM and CRISPR totals are also reported separately and collapsed by MOA class; the decision stays on the uncollapsed adjusted total.

**E7. Track 1 input.** The primary analysis uses **BUB1B alone** as the causal module (scheme W3). `config/mechanism_screen.json` classifies MAD1L1 and CENPE as "pathway context only" (non-coding structural calls) and IGBP1 as "insufficient evidence" (VUS), so they cannot be treated as causal mechanisms. W1 (weights from Track 1) and W2 are reported as sensitivity labelled "pathway context". The two BUB1B alleles (p.Leu737*, p.Asn1002Lys VUS; phase unresolved) are stated as two direction scenarios, complete loss and hypomorph. Line D stays **exploratory** and does not assign tiers A/B/avoid; it is shown as a direction note with its stated assumption.

**E8. Two-record rule as the primary one.** Measured drug-target links need >= 2 activity records and pChEMBL >= 6, with no similarity-predicted links, in **both** B1 and B2 and in the register. The single-record version is the sensitivity analysis. Declared as post hoc: the rule was found after fostamatinib's single-record support was read.

**E9. Candidate rule.** Candidate = FDA-eligible, systemic, not N01/N05C, supported by >= 2 usable lines among A, B (E1), C. Amendment 3 (one line plus literature) is kept only as a labelled list "one line plus literature", never as candidates. The hypergeometric overlap is computed on the **intersection** of the universes that each line can score, not on the universe of A.

## New pre-specified analysis

**E10. Class-level analysis.** Unit = mechanism class (ChEMBL mechanism-of-action group, with ATC level 4 as a second definition). For every class with >= 3 scorable drugs: the mean score in A (PageRank rank), the proximity of B1, and the sensitivity of C are compared with 10,000 random classes matched in size and mean drug degree; BH-FDR < 0.10 across classes within each line. A class is *supported* if it passes in >= 2 lines. Output: hypotheses by class, each with an evidence level 0 to 3 (0 = proximity or graph score only; 1 = plus an independent line; 2 = plus a measured effect in a human cell system; 3 = plus functional evidence with a verified PMID and the exact sentence that supports it). This replaces drug-by-drug claims; a class with no verified sentence cannot exceed level 1.

**E11. Caution list.** Pre-specified classes whose use in this patient needs explicit caution: tubulin binders and vinca alkaloids, taxanes, SAC/Aurora/PLK/TTK inhibitors, genotoxic chemotherapy, proteasome inhibitors in the constitutional context, growth hormone with a neoplasm term. A class enters only with the exact label or paper sentence quoted. Every proposed or cautioned drug needs a manual label review; the automatic classifier output `MIXED` is recorded as NOT_ASSESSABLE, and no age band is assumed.

## Scope and age (conditions set by the owner on approval)

**E14. Scope: repurposing of approved drugs only.** Every drug, class and caution entry at every stage (universe, lines A to D, class analysis E10, caution list E11) must be an approved drug: FDA-approved (openFDA) and ChEMBL `max_phase = 4`. Investigational, withdrawn and unapproved compounds are excluded, including PRISM compounds that are not approved. *Repurposing* means that the hypothesised use differs from the labelled indication. For each proposed item the report states the labelled indication and the hypothesised use side by side. A drug whose label already covers a feature of the child (for example somatropin for short stature in children born small for gestational age) is **on-label supportive care, not repurposing**, and goes to a separate appendix, never to the candidate list; its label contraindications (for somatropin, active malignancy, with a rhabdomyosarcoma term in the child) are quoted. Output is hypotheses for investigation and for a clinician, never advice, dosing or a recommendation to treat.

**E15. Age is unknown: fixed thresholds and rules.** The child's age is not known and the HPO annotations of MVA have no onset, so no age scenario is built and no drug is described as suitable for the patient. The label support is reported for fixed bands (ICH E11 / FDA pediatric subpopulations) in `labels.BANDS`: neonate 0 to < 1 month (28 days), infant 1 month to < 2 years, child 2 to < 12 years, adolescent 12 to < 18 years, adult 18 years and older; and for the three reporting groups already used by `age_groups.py`: **0 to < 3 years, 3 to < 12 years, 12 years and older**. Rules, fixed before any run:
- The status of a drug in a group comes only from the label text: `labelled` (label supports use across the whole group), `partly` (from a stated minimum age inside the group, which is reported), `not_labelled` (paediatric use not established), `adult_only`, `age_not_stated` or `verify` (the label is mixed or could not be read, or its paediatric section only refers to paediatric information in another section without stating that paediatric use is established). The age band is never assumed.
- A drug is *proposed* for a group only if its status there is `labelled`. With `partly` it is listed with the minimum age and the group it covers. With `not_labelled`, `adult_only`, `age_not_stated` or `verify` it is listed as a hypothesis for a clinician and flagged "no paediatric label for this group".
- `MIXED` from the automatic classifier is recorded as `verify`; each proposed or cautioned drug needs a manual label review in `config/pediatric_label_review.json`, with the label sentence and its source (DailyMed or openFDA) and the date read.
- The report shows the table for all three groups for every item; the effect of age is on the label support only, and the paediatric warnings (growth-plate findings, cardiomyopathy, bone-marrow suppression, etc.) are quoted from the label.
- A drug with a paediatric label for none of the groups can still appear as a class hypothesis (E10) but not as a proposal.

## Reporting and reproducibility

**E12. Forking-path ledger.** `out/forking_ledger.json` lists every variant run (lines, STRING thresholds, weighting schemes, top-N, k, link rules). One primary analysis is declared per line (E1 to E9); everything else is sensitivity and is never presented as a finding.

**E13. Checksums.** `out/checksums.sha256` is regenerated after all runs and `verify_outputs` must pass.

## Runs, in order (after approval)

1. Restore sources by `fetch_sources` (SHA-256 against `data/sources/record.json`): PrimeKG, MONDO, HPO (`hp.obo`, `genes_to_phenotype.txt`, `phenotype.hpoa`), DepMap CRISPR 24Q2, plus the DepMap 24Q2 somatic mutation file (new, for E6); OmniPath and the DepMap expression/PRISM files are already on disk.
2. `module_rule_validation` (E1), `kg_nullquery` (E4), line C with 1,000 sets and TP53 (E6).
3. Re-run B1 validation (E2) and the lines with E5, E7, E8.
4. `four_lines` with E9, class analysis (E10), caution list (E11), ledger (E12).
5. Rewrite report, methods, abstract (<= 500 words) and README on the actual result; regenerate checksums.

## What would count as failure

If fewer than two lines are usable under these rules, the result is reported as a negative result for candidate drugs, and the deliverable is the class-level hypotheses, the caution list and the reproducible framework. No threshold is changed to avoid that outcome.
