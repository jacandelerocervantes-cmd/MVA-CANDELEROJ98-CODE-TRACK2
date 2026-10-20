# Addendum: deterministic differential diagnosis of the Track 1 genes

Written on 2026-10-15 **before** the code (`src/track2/differential.py`) was run. It adds one analysis; no earlier rule changes.

## Why

The report should state by rule which gene is the causal candidate, which are contextual, and which are insufficient or excluded, instead of
using a hand-written table. The rules below are fixed first and are justified from the criteria, not from the expected answer. They are not
tuned against `config/mechanism_screen.json`, which is read only at the end, for comparison.

## Inputs (gene-level fields only)

* `out/mva_landscape.json`: tier (`curated` / `proposed, not curated`), inheritance as annotated, sources.
* `out/submission_genes.json`: ClinVar classification and finding types per gene.
* Track 1, read-only: `candeleroj98_genomewide_ranking.csv` (columns gene, type, n_high_impact, log10_joint_rarity, phenotype_score only),
  `candeleroj98_phenotype_summary.json` (candidate scores, phenotype rank among annotated genes, number of annotated genes), and from
  `candeleroj98_gatk-manta-fusion.csv` only `finding_type` and the first word of `notes` (a gene symbol). Genomic positions, alleles,
  genotypes and read metrics are never read, printed or stored.
* Genes evaluated: the union of the genes of the Track 1 ranking and of the structural-call file.

## Criteria

* **C1, MVA status.** `curated` (listed by HPO and/or Orphanet for an MVA / premature chromatid separation disease), `proposed`
  (literature only, "proposed, not curated"), or `none` (absent from the landscape).
* **C2, genotype pattern versus inheritance.** `consistent` = the Track 1 type is `compound_het` or `homozygous` (two rare alleles in one gene)
  AND at least one annotated inheritance string contains "recessive". `biallelic_inheritance_not_recessive` = two alleles but no recessive
  annotation (including "not annotated"). `hemizygous` = X-linked single copy; it is consistent only if the annotation contains "X-linked",
  otherwise `inconsistent`. `unpaired` = the gene appears only as a structural call (one call, no second allele): never consistent.
  Phase of a compound heterozygote is not resolved at this level and is stated as a caveat, not tested.
* **C3, consequence class** (gene level only): `high_impact` if `n_high_impact` >= 1 (loss-of-function class of Track 1; predicted decay is not
  available at gene level, so this is a proxy and is stated as such); `non_high_impact` if the gene is in the ranking with `n_high_impact` = 0
  (missense, in-frame or other; variants of uncertain consequence); `structural_unpaired` if the gene appears only as a structural call. A
  structural call is not assumed to overlap coding sequence, because no coordinate is read.
* **C4, phenotype fit.** Track 1 phenotype score; rank among annotated genes (null if the gene has no HPO annotation); `has_hpo` is true if the
  rank is not null, or, for a gene outside the candidate score table, if HPO is among its landscape sources.
* **C5, ClinVar.** `P/LP` if any record is classified pathogenic or likely pathogenic; `uncertain` if any is uncertain/conflicting and none is
  P/LP; `none` otherwise (including no record or an unclassified record).

## Verdicts (applied in this order; first match wins)

1. **excluded**: (a) the gene has a Track 1 finding type of `artefact`, `common_variant` or `filtered`; or (b) `log10_joint_rarity` < 2,
   i.e. the joint frequency of the alleles is above 1e-2, which is the usual common-variant convention (a recessive disease allele pair cannot
   be that frequent). Justification: such a gene cannot explain a rare recessive disease whatever its annotation.
2. **causal candidate**: C1 = `curated` AND C2 = `consistent` AND (C3 = `high_impact` OR C5 = `P/LP`). Justification: each part answers one
   requirement of a recessive monogenic diagnosis (known disease gene, two alleles, at least one allele with a damaging effect or established
   pathogenicity). A `proposed` gene never reaches this verdict because its gene-disease link is not curated.
3. **insufficient evidence**: C1 = `none`; or `has_hpo` is false; or the gene has a consistent biallelic genotype (C2 = `consistent`) with
   neither a high-impact allele nor P/LP ClinVar (variants of uncertain significance only).
4. **contextual**: any other MVA gene (C1 curated or proposed), i.e. an MVA gene without a qualifying coding genotype (C2 not consistent, or
   a single/unpaired structural call). It supports the pathway, it does not explain the case.

The order matters: a non-MVA gene is insufficient even if its genotype is biallelic; an MVA gene without a biallelic genotype is contextual,
not insufficient, because its relevance is the pathway, not the variant.

## Limits

Gene-level fields only; C3 is a proxy; phase and protein-level effect are not assessed; the verdict is a rule-based triage, not a clinical
classification. The comparison with the hand-made screen is reported afterwards and neither file is edited.
