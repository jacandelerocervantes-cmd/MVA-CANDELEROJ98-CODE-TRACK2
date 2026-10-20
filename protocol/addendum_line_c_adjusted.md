# Addendum: line C at drug level, adjusted for lineage and proliferation

Fixed on 2026-10-18 before any adjusted statistic is computed. Offline: only data already on disk (`data/depmap`: PRISM Repurposing 24Q2 primary screen, `Model.csv`, expression). No patient data. Nothing is downloaded. Existing files are not edited.

## 1. Question

The original line C statistic (`track2.multi_evidence_completion.line_c_drug_level`, protocol section 11.3 of `addendum_multi_evidence.md`) is the Spearman correlation between a drug's PRISM viability (log fold change at 2.5 uM) and the checkpoint score of the cell line (mean z of BUB1B, BUB1, BUB3, MAD1L1, MAD2L1, TTK; low = weak checkpoint). It was not adjusted for cell lineage or for proliferation, and the score correlates 0.81 with a proliferation score. Question: do the line C candidates (lazertinib, neratinib, prednisolone, trametinib, cobimetinib) and the other original hits (23 drugs with BH q < 0.10) survive adjustment for these two confounders?

## 2. Data and variables (fixed)

- Cell lines, score, drug matrix, drug-to-PRISM-row matching and pool of 224 drugs: exactly as in the original (895 lines with expression; several PRISM rows of one drug averaged; drugs with fewer than 100 lines with both values are not estimable).
- Lineage: `OncotreeLineage` of `Model.csv`. Within a drug's analysed lines, lineages with fewer than 10 lines are pooled into one level "small lineages"; lines with no lineage are placed in that level.
- Proliferation covariate: `data/depmap/Model.csv` has no doubling-time or growth-rate column (only `GrowthPattern`, a categorical culture description), so no measured doubling time is available. The covariate is the mean z-score (over cell lines, log2(TPM+1) from `OmicsExpressionProteinCodingGenesTPMLogp1.csv`) of the fixed gene set MKI67, PCNA, TOP2A, MCM2, CCNB1. None of these is in the checkpoint score. (This differs slightly from the six-gene set of `prism_l4`; it is the set stated here, fixed before computing.)

## 3. Adjusted statistic (primary)

For each drug, on the lines with a PRISM value: rank-transform the viability, the score and the proliferation covariate; regress the ranked viability and the ranked score each on lineage dummies (fixed effects) plus the ranked proliferation covariate; the adjusted effect is the Pearson correlation of the two residual vectors (a partial rank correlation, rho_adj). Sign convention unchanged: rho_adj > 0 means the drug is more lethal in weak-checkpoint lines.

Null: 1,000 permutations of the score WITHIN lineage level (seed 20261018), the permuted score re-ranked and re-residualised each time; two-sided p = (1 + number of permuted |rho_adj| >= observed |rho_adj|) / 1001. Benjamini-Hochberg across the same drug pool as the original (all drugs estimable, 190 in the original).

**Survival rule:** a drug survives if the adjusted BH q < 0.10 and sign(rho_adj) = sign(original rho). Verdicts: survives, does not survive, not estimable (not in PRISM or fewer than 100 lines).

## 4. Secondary and sensitivity analyses

1. **Within-lineage analysis:** lineages with at least 30 lines having the drug's value; in each, the partial rank correlation (proliferation partialled out); Fisher z with variance 1/(n-3), combined by DerSimonian-Laird random effects (pooled rho, p, I-squared) and, as a check, by Stouffer's weighted Z. Reported; it does not change the survival rule.
2. **Leave-one-lineage-out (LOLO)** for lineages with at least 30 lines: the primary analysis repeated without that lineage (same 1,000-permutation null, seed). Reported: number of lineages omitted in which the sign is kept and the omitted-lineage p < 0.05; minimum and maximum rho_adj. Run for the original hits.
3. **Excluding haematopoietic and lymphoid lines** (`OncotreeLineage` in Myeloid, Lymphoid): the whole primary analysis repeated on the remaining lines, with BH across the pool. Suspicion to test: prednisolone reflects the glucocorticoid sensitivity of lymphoid lines. The suspicion is called supported if prednisolone's original association disappears in this analysis (q >= 0.10 or sign reversed) while the original q < 0.10, or if its adjusted rho falls below half of the full-sample rho_adj.

## 5. Reported

How many of the original hits (q < 0.10 in `out/line_c_drug_level.json`) survive; Spearman correlation between original rho and adjusted rho over all estimable drugs and over the original hits; survival table for the five line C candidates and for all original hits; LOLO and lineage-excluded results; a plain verdict per candidate.

## 6. Declared limits

PRISM association from cancer cell lines is correlational. It does not model partial BUB1B loss (a constitutional mosaic) and the checkpoint score is an expression proxy. The gene-set proliferation covariate is an imperfect proxy for growth rate, so residual confounding by proliferation is possible; lineage fixed effects remove between-lineage differences but not within-lineage heterogeneity. A surviving drug is still only a hypothesis signal.

## 7. Clarification added after the first run (does not change any threshold)

The first run showed that the screen rows of prednisolone (and of the drugs measured in the same PRISM screen) contain no haematopoietic or lymphoid lines, so excluding them changes nothing for these drugs. The lymphoid suspicion (section 4.3) is therefore **not testable** for prednisolone from these data and is reported as not testable (the criterion "falls below half" and the q rule are evaluated only when the drug has such lines in its data). Its non-survival in the primary analysis is attributable to lineage and proliferation adjustment among solid-tumour lineages, not to lymphoid lines.
