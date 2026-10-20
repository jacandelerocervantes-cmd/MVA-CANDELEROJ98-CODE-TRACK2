# Deviations and clarifications after the freeze

1. **Stability verdict (P3).** The `verdict` field of `out/validation_results.json` is computed only from the test that drops 20 % of the seeds (mean Jaccard 0.81, "STABLE"). The protocol rule ("a top 20 with mean Jaccard < 0.5 is declared unstable") applied to the STRING threshold sensitivity gives Jaccard 0.48 (threshold 400) and 0.38 (threshold 900), so the ranking is **declared unstable with respect to the interactome threshold**. The code was not changed; the report applies the rule to both families of perturbation.
2. **P1.** The decision rule is met (AUC 0.92, interval 0.89–0.95: "supported at pathway level"), but no left-out seed reaches z ≤ −2 (median −0.84): the primary threshold is very strict for a single gene. The threshold was not modified; this is reported.
3. **P2.** Metformin, cisplatin and mitomycin cannot be scored (ChEMBL has no curated mechanism with a protein target); oxymetholone falls at the 59th percentile (z = +0.49). The method **does not recover** the available Fanconi positive control.
4. **P4.** The Spearman correlation between z and the mean degree of the targets is 0.28 (degree normalisation does not completely remove the bias).
5. **Freeze date.** See `freeze_note.md`.

## Submission-derived modules (2026-10-06)

- `addendum_submission.md` replaces the hand-defined modules as the source of the final ranking (see its rule). The frozen protocol's validation (P1 to P4) still refers to the hand-defined checkpoint module in `config/mechanisms/`; it was not rerun for the submission-derived modules.
- The addendum first stated the shortlist rule as z <= -1.5; the code (unchanged since the first run) uses z <= -2 (the primary threshold) or top 15 per module. The addendum text was corrected to match the code; no result depends on the wrong number.
- Reactome answers 404 for proteins that belong to no pathway (IGBP1); this was handled as "no pathways", which triggers the pre-registered STRING-neighbour fallback.
- IGBP1 (5 proteins with its STRING neighbours) and DNAH10 (not in the network) were not modellable and are reported, not scored.
- The transcriptomic reversal layer was computed for the earlier shortlist and was not rerun (it is not used in the ranking); CPIC drug relevance was recomputed for the new shortlist without rerunning the VCF genotyping.
- Amendment 1 to `addendum_submission.md` (second-order STRING neighbours as a fallback) was introduced after the first submission-derived run, because two of five genes had contributed nothing. It is disclosed as post hoc.

## Expanded links, FDA scope and plausibility (2026-10-06)

- `addendum_dti.md`, `addendum_plausibility.md` and `addendum_pediatric.md` were written before the analyses they describe. Changes made afterwards, all disclosed:
  - A first run of E2 mapped each STRING node to a single UniProt accession (`inv = {node: accession}`), so most module proteins were missed (3 targets with ligands, E2 not valid). The mapping was corrected to use every accession that maps to a module node (51 targets with ligands, 149 known links, E2 valid at tau 0.65). The same bug affected the E1 module-link count (9 drugs before, 57 after) and the module-protein matching in `novel_hypotheses.py`. It was fixed before any result was interpreted.
  - The "systemic route" gate was reworded to "not topical-only" because ChEMBL's route flags are incomplete (bortezomib has none recorded).
  - Drug class (antineoplastic or not) is also decided from the approved indications, because ChEMBL lacks an ATC code for several antineoplastics (ixazomib, docetaxel).
  - Paediatric label review (`config/pediatric_label_review.json`) and the age groups (0-<3, 3-<12, 12 and older) were added after the first paediatric result; see `addendum_pediatric.md`.

## PRISM (L4) and docking (2026-10-07)

- `addendum_prism.md` was fixed before the data were analysed; the only change made before the analysis was to replace a within-lineage comparison (too few cell lines per group) by a lineage-adjusted regression, recorded in the addendum.
- `addendum_docking.md` was fixed before any structure was downloaded. All seven targets failed its validation and no docking score was used. A check of the redocking RMSD with symmetry-aware matching and the best of nine poses was made to rule out a measurement artefact; it changed no verdict.
- `addendum_wide.md` (a wider ring, top 50 instead of top 20) was run: it adds a single drug to the strict list, which says the short list reflects the structure of the signal and not the cut-off.

## Changes made after the first independent review (2026-10-07)

Each item was first verified by the author; the review's claims were not taken as true.

- **Generic complexes (verified).** The modules contain 30 proteasome subunits plus ADRM1 and 24 tubulins (BUB1B, CENPE and MAD1L1 modules), so proteasome and tubulin drugs have distance 0 by construction. `addendum_core.md` adds a sensitivity analysis without them and `L1 core` is reported next to L1; without them only 5 drugs stay robust, and none of the Tier 1 drugs.
- **L3 rule (verified).** The code excludes the class `mva_tumour_case_reports` from L3 (case reports are not a mechanistic rationale), which `addendum_plausibility.md` did not state. Read literally, the addendum would give 13 Tier 1 drugs instead of 3. The exclusion is kept and is now stated in the addendum.
- **Single-measurement links (verified).** Every measured link that carried the new hypotheses (fostamatinib, selumetinib, daunorubicin) rested on one activity record. A measured link now needs at least two records to count as L2 (`MIN_MEASUREMENTS`); the change removes those drugs from Tier 2.
- **PRISM and proliferation (verified).** The SAC score correlates with a proliferation score (Spearman 0.81 to 0.84). `addendum_prism.md`, amendment 1: the call is made on a model adjusted for lineage and proliferation. Sirolimus is explained by proliferation; the microtubule agents are not.
- **Paediatric labels (verified against the openFDA label text).** Colchicine (familial Mediterranean fever from 4 years), selumetinib (from 1 year) and daunorubicin (acute lymphocytic leukaemia in children, age not stated) are labelled; paclitaxel and natalizumab are not; fostamatinib is not recommended under 18 years (growth-plate effects in juvenile animals). `config/pediatric_label_review.json` was corrected and the summary sentence is now generated from the data.
- **E2 out of sample.** The cut-off was chosen on the same pairs it was tested on. A split-by-target check gives a precision of 0.44 to 0.48, below the criterion. No drug of the register depends on E2.
- **Extrapolation arm.** L2 and L3 were set to true in code; they are now computed (L2 from the direct module targets in the comparison run, L3 from the literature class).
- **Privacy.** The patient-level summaries `pgx_summary.json` and `genome_covariates_summary.json` (possible star alleles, depth per chromosome) were moved out of `out/` to `data/private/`, which is not published.
- **Generality.** Gene names were removed from `run_all.py` and `build_report.py`; the genes of the expression scores are in `config/prism_genes.json`; a malformed submission file or a submission with no modellable gene now stops with a clear message.
- **Not changed.** The hand-defined modules of the extrapolation arm remain circular (seeded with pathway genes) and are declared so; docking is reported as not validated, not as proof that docking cannot work.

## Line B implementation (2026-10-07)

- **MGI not used in the phenotype module** (no MP-to-HPO mapping downloaded); the pre-registered text said "if the licence allows", and the module is built from HPO only.
- **Nine configurations for the single phenotype module** are defined as 3 STRING thresholds x 3 module definitions (see `addendum_four_lines.md`, implementation notes), because the weighting schemes need several modules. Written before the results were seen.
- **Placeholder gene symbol.** The HPO file marks some annotations with `-` instead of a gene symbol; the first run counted `-` as a gene (1,042 genes instead of 1,041). It was found by reading the top of the ranking, fixed and covered by a test; no result was used before the fix.

## Line A (2026-10-07)

- **Side-effect edges removed after the first patient ranking** (see `addendum_four_lines.md`, amendment 2). First validation, graph with `drug_effect` edges: mean AUC 0.821 (95 % interval 0.801-0.839), 200 diseases. The patient ranking of that graph was not used.
- **Rhabdomyosarcoma (HP:0002859) is not in PrimeKG**; the patient query uses 7 of the 8 terms.

## Candidate rule widened (2026-10-07)

- **Tier 2 added: one validated line plus literature** (`addendum_four_lines.md`, amendment 3), after the original rule (two usable lines) gave a single candidate. The literature search is scripted and its rule was fixed before it was run; the abstracts are read by the author and the claims are written down.
- **Age scenarios not built.** The HPO annotations of mosaic variegated aneuploidy (216 rows, OMIM and ORPHA) have no onset value, so "phenotypes with an onset compatible with the age" cannot be derived. Label support by age group is kept as the age information.
- **Line C not used** (22 hits against a null 95th percentile of 23.0) and **line D not used** (4 of 6 controls). Six of the 19 PRISM hits are MDM2 inhibitors and MDM2 is a CRISPR hit: a TP53-status confounder that lineage and proliferation do not capture; it was not re-analysed.

## Validation of line B replaced (2026-10-07; amendment 1 of `addendum_four_lines.md`, written before the new validation was run)

- The pre-registered validation of line B (the 30 diseases of the earlier benchmark, each with a module from "its own HPO terms") could not be run: those diseases are broad and have no HPO annotation; the direct match through Orphanet gave 0 diseases with at least three scorable approved drugs. It was replaced by a validation on diseases that have OMIM cross-references in MONDO, at least 5 HPO terms and at least 3 known approved drugs (19 diseases). **The thresholds were not changed** (mean AUC >= 0.60, lower bound of the interval > 0.50, at least 15 evaluable diseases).
- The replacement favours diseases whose causal genes sit in their own phenotype module. The independent review of 2026-10-07 also notes that the validation tests the phenotype-module method only: the mitosis module without generic complexes (which supplies 5 of the 6 drugs supported by line B, among them fostamatinib) was **not** validated by it, and the measured-link rule of two records (`MIN_MEASUREMENTS`) was applied in the plausibility tiers but not in the proximity analyses. Both points are open limitations, not corrections.

## Independent review of thresholds (2026-10-07)

Two read-only reviews in a fresh context (network-pharmacology and validation-audit roles in `.claude/agents/`) found, and the author verified in `out/`: the degree-matched permutation null of line A has mean AUC 0.773 against 0.802 observed (the 0.70 threshold is passed by drug degree alone; the excess is about 0.03), while that of line B is 0.539 against 0.717; line C failed by 22 hits against a null 95th percentile of 23.05 with only 200 null sets; line D recovered 4 of 6 controls, two of which are direct readout nodes and pass by construction. The overlap of lines A and B (one drug, fostamatinib) has p = 0.099 and does not change under any setting of the top-N of A or of k for B that was tried. No threshold was changed as a result.

## Addendum 11 approved (2026-10-07)

The owner approved protocol/addendum_review2.md (amendments E1 to E15) after the second independent review, with two conditions that became E14 (repurposing of approved drugs only) and E15 (age unknown: fixed bands and label-only rules). No threshold was loosened. The amendments are post hoc with respect to the results already produced and pre-specified with respect to the new runs. Earlier results (including the 200-set estimate of line C, the single-record links and the candidate fostamatinib of out/four_lines.json) are kept as sensitivity analyses. There is no git repository by the owner's decision; the order of events is attested by the author only.

