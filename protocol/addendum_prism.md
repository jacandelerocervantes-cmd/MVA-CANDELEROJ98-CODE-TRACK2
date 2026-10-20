# Addendum 7: independent functional evidence (L4) from DepMap and PRISM

Written on 2026-10-07 **before** the data are downloaded or analysed. It fills the evidence line L4 of `addendum_plausibility.md`, which was "not assessed".

## Question

The proximity method does not say whether a drug helps or harms a cell with a weakened spindle assembly checkpoint (SAC). Public screens of drug sensitivity
in hundreds of cancer cell lines can say whether cells with a **weaker SAC** are more or less sensitive to a candidate drug. This is an independent, empirical
check of the direction of effect, with the limits stated below.

## Data (public, CC BY 4.0, downloaded from figshare with MD5 checks)

- PRISM Repurposing Public 24Q2 (figshare article 25917643): `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` (cell-line viability, log fold change,
  one value per compound and cell line), `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv`, `Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv`,
  `Repurposing_Public_24Q2_Treatment_Meta_Data.csv`, README.
- DepMap Public 24Q2 (figshare article 25880521): `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (log2(TPM+1)) and `Model.csv`.
The two releases are the same quarter so that cell-line identifiers match.

## Rule (fixed before running)

1. **SAC expression score** of each cell line: the mean over the genes BUB1B, BUB1, BUB3, MAD1L1, MAD2L1 and TTK of the z-score (across cell lines) of log2(TPM+1).
   A **weak-SAC group** is the lowest quartile of the score and a **strong-SAC group** the highest quartile. A second score uses BUB1B alone.
2. **Drugs tested:** the candidates of `out/proposals.json` and of the wider ring that are present in the PRISM compound list (match by name, case-insensitive,
   salt words removed, and by InChIKey first block when the list has it). Drugs not in PRISM are reported as `not_in_PRISM`; nothing is imputed.
3. **Statistic:** for each tested drug, the difference in mean viability log fold change between the weak-SAC and strong-SAC groups (a negative difference means the drug
   kills weak-SAC cells more), a Welch t-test, and Benjamini-Hochberg correction across the tested drugs. A Spearman correlation between the SAC score and the viability over all lines
   is reported as a sensitivity analysis.
4. **Rows and calls:** each compound-dose row of PRISM is tested separately (most compounds have one row at 2.5 uM); the drug-level call uses all its rows and is `mixed` if rows disagree.
   **L4 call:** `supports_sensitivity_when_SAC_weak` if FDR < 0.10 and the difference is negative; `supports_resistance_when_SAC_weak` if FDR < 0.10 and the
   difference is positive; `no_association` otherwise; `not_in_PRISM` if untested. L4 does not change the tiers; it is reported next to them.
5. **Reported with every call:** the numbers of cell lines per group, the effect size and the FDR.

## Limits (stated now)

- These are tumour cell lines, not MVA cells; expression is a proxy for checkpoint strength, not a measurement of function; aneuploidy, proliferation rate and tissue of origin
  confound the comparison. A **lineage-adjusted** model is run as the sensitivity analysis (written before the data were analysed; it replaces a within-lineage comparison that would
  have too few cell lines per group): viability ~ SAC score + lineage indicators (OncotreeLineage), the t-test of the SAC coefficient, Benjamini-Hochberg across the tested rows. A primary
  call that is not confirmed by the adjusted model is reported as `lineage_confounded`.
- The primary PRISM screen uses one concentration and pooled cells; a drug that is toxic in all lines shows no difference even if it would act in a patient.
- A result here is a hypothesis about direction, never a clinical statement.

## What would count as failure

- Fewer than 5 candidates present in PRISM: reported; L4 stays largely "not tested".
- A strong association that disappears within lineages: reported as lineage-confounded.

## Amendment 1 (2026-10-07, after the first independent review and before the re-run; disclosed)

The review found, and it was verified, that the SAC expression score has a Spearman correlation of 0.84 (BUB1B alone: 0.81) with a **proliferation score** (mean z-score of
MKI67, TOP2A, PCNA, MCM2, MCM6 and RRM2, genes that are not in the SAC score) across 1,517 cell lines. Checkpoint genes are proliferation genes, so the lowest-quartile
group is largely the slowest-growing lines, and antimitotic drugs act on proliferating cells. The unadjusted comparison therefore cannot separate a weak checkpoint from slow growth.

Changes (decided before the re-run):
1. The proliferation score is added as a covariate: the sensitivity model becomes viability ~ SAC score + proliferation score + lineage indicators.
2. **The call is made on the adjusted model.** `supports_sensitivity_when_SAC_weak` or `supports_resistance_when_SAC_weak` needs FDR < 0.10 in **both** the unadjusted
   comparison and the model adjusted for lineage and proliferation, with the same sign. A drug that is significant only unadjusted is `confounded_by_proliferation_or_lineage`.
3. The unadjusted result, the adjusted coefficient and the proliferation coefficient are all reported, so that the effect of the adjustment is visible.
