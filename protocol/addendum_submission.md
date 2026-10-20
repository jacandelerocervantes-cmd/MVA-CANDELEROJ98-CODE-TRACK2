# Addendum: candidates generated from the Track 1 submission

Written on 2026-10-06, **before** the submission-derived run. It does not change the frozen protocol
(`frozen_protocol.md`); it changes only where the disease modules come from.

## Why

In the first version of Track 2 the disease modules were written by hand from the author's reading of the
condition. That makes the drug list depend on a human choice. This addendum replaces it: the only input about the
patient is the Track 1 submission file (`submissions/track1/candeleroj98_gatk-manta-fusion.csv`), and every module
is derived from it by a fixed rule.

## Rule (fixed before running)

1. **Genes and weights.** Each row of the Track 1 CSV names a gene (first word of `notes`). The gene weight is the
   best `epcr` among its rows. `epcr` is the submission's own ranking value; it is not a probability.
2. **Verification.** Every single-nucleotide variant of the file is looked up in ClinVar by GRCh38 position (NCBI
   E-utilities). The result (classification, review status) is reported next to the gene. It is **context, not a
   filter**: a gene is not dropped because ClinVar has no entry, and not promoted because it has one.
3. **Module per gene.** The UniProt accession is resolved from UniProt (reviewed, human). Every Reactome pathway
   mapped to that accession with 10 to 300 UniProt participants is merged, together with the gene itself, into one
   module. If no pathway qualifies, the module is the gene plus its STRING physical neighbours (score >= 700). A
   module with fewer than 10 proteins in the network is reported as "not modellable" and gets no score.
4. **Scoring.** Unchanged: closest-distance proximity to the module on STRING v12.0 physical (700, sensitivity 400
   and 900), z-score against 1000 degree-matched random draws (seed 20261006), approved ChEMBL_37 drugs only.
5. **Combination.** Combined score = sum over genes of `epcr_weight * max(0, -z)`. Shortlist = z <= -2 in any
   module or top 15 in any module (as before).
6. **Hand-defined modules.** The earlier modules (spindle checkpoint, proteotoxic stress, p53/senescence, growth
   signalling) are kept only in `config/mechanisms/` and `out/comparison_manual/` as a comparison. They do not
   feed the final ranking. The comparison reports the overlap of the two shortlists.

## What would count as a failure of this approach

- Fewer than 3 genes of the submission yield a modellable module (the file would be too thin to drive the analysis).
- Rankings dominated by a single large generic pathway (for example "Cell Cycle"): reported, not corrected silently.
- Any deviation from this rule is added to `deviations.md` with its reason.

## Amendment 1 (2026-10-06, made AFTER the first run; disclosed)

The first run scored only 3 of the 5 genes in the Track 1 file. IGBP1 (no Reactome pathway; 5 proteins with its first-order STRING neighbours) and DNAH10 (absent from the physical network at score >= 700) produced no module, so the final ranking reflected only mitotic biology and ignored the other genes of the submission. Because that was a limitation of the rule rather than evidence about the genes, the fallback in step 3 was extended: if the gene plus its first-order neighbours has fewer than 10 proteins, the neighbours of those neighbours (second order) are added. The change was decided after seeing the first run; no drug ranking was used to choose it. Genes that still have no module (absent from the network at the primary threshold) remain reported as "not modellable"; no value is invented for them.
