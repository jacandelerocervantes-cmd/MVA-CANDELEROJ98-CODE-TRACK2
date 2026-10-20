# Addendum 10: four independent lines of evidence (design, written before any new data are downloaded or analysed)

Written on 2026-10-07, after the first independent review. It does not change any earlier result; it adds a new analysis that replaces literature as the basis of a candidate.

## Why

The review showed, and it was verified, that the proximity of the Tier 1 and Tier 2 drugs is membership of generic complexes (proteasome, tubulin) in modules built from Reactome
pathways, and that without them almost nothing stays close (`addendum_core.md`). What then supports those candidates is the literature, which is context and not a technical
basis. The clinical picture of the child (eight HPO terms: rhabdomyosarcoma, nephrocalcinosis, short stature, failure to thrive, skeletal muscle atrophy, premature birth, small for
gestational age, and the parental history of recurrent miscarriage) was also not used to generate candidates. This addendum defines four lines of evidence that use different information, a
rule for combining them, and the validation each line must pass before it can be used.

## The four lines

**A. Phenotype-guided knowledge graph.** A graph that links diseases, HPO phenotypes, genes, pathways and approved drugs (PrimeKG, or its successor if the data licence allows). Query nodes: the
disease node of mosaic variegated aneuploidy and the patient's HPO terms, weighted by the age scenario (below). Method: the pretrained TxGNN model if it runs on CPU in the VM; otherwise personalised
PageRank (restart on the query nodes) over the graph. A drug is *supported by A* if it is among the 50 best-scored FDA-approved drugs.
*Validation required:* leave-out of known indications: for 200 diseases with at least three approved drugs, hide those drug-disease links, rank all drugs for the disease and compute the ROC AUC of the
hidden drugs. The line is used only if the mean AUC is >= 0.70 with a bootstrap lower bound > 0.55.

**B. Network proximity to a phenotype module and to the mitosis module.** Same proximity method as before (STRING physical, degree-matched null, 1,000 draws, seed 20261006, nine configurations),
with two modules: the mitosis modules of the earlier analysis **without the generic complexes** (`addendum_core.md`), and a **phenotype module** made of the genes annotated to the patient's HPO terms (HPO
`genes_to_phenotype`, with the descendants of each term) plus, if the licence allows, the genes whose mouse orthologue shows the matching phenotype (MGI). To keep the module specific, a gene enters only if it is
annotated to at least two of the eight terms; if that gives fewer than 10 genes, the 100 genes with the highest information-content-weighted score are used. A drug is *supported by B* if it is among the 20 best in at least 7 of 9
configurations for either module.
*Validation required:* for the 30 diseases of the earlier benchmark, build a phenotype module by the same rule from the disease's own HPO terms and compute the AUC of its known approved drugs against the other
drugs. The line is used only if the mean AUC is >= 0.60 with a lower bound > 0.50.

**C. Measured vulnerability of cells with a weak checkpoint.** DepMap CRISPR gene effect and the full PRISM Repurposing screen (all compounds, not only the candidates), against the SAC expression score and
the proliferation covariate of `addendum_prism.md` (amendment 1). Drug level: compounds whose sensitivity is higher in weak-checkpoint lines, with FDR < 0.10 in the model adjusted for lineage and proliferation, the same sign
unadjusted, and a viability difference of at least 0.2. Gene level: genes whose knockout is more lethal when the checkpoint is weak (same adjusted test); drugs whose curated target is such a gene are *supported by C*.
*Calibration required:* repeat the whole analysis with 200 random gene sets of six genes matched to the SAC score in correlation with proliferation (+-0.05); the line is used only if the number of hits observed is above the 95th percentile of the null.
*Limits stated now:* tumour cell lines are a proxy; this line speaks to the tumour context.

**D. Direction of effect from a signed network.** Directed, signed interactions (OmniPath, SIGNOR). For each candidate, the sign of the net effect of its action (inhibitor, agonist...; ChEMBL action type) on a fixed set of checkpoint-strength nodes (BUB1B, BUB3, MAD2L1 and CDC20 as
the mitotic checkpoint complex, with APC/C activity counted with the opposite sign), along signed paths of at most three steps; the call is `increases checkpoint strength` (favourable), `weakens it` (unfavourable) or `unclear`.
*Validation required:* six textbook perturbations with a known sign, fixed here: TTK inhibition weakens the checkpoint; BUB1 kinase inhibition weakens it; AURKB inhibition weakens it; CDC20 or APC/C inhibition prolongs the arrest (strengthens it); MAD2L1 loss weakens it; TRIP13 inhibition strengthens it.
The line is used only if at least five of the six are recovered; otherwise it is reported as exploratory and cannot move a candidate to a higher tier.

## Combining the lines (fixed now)

- A drug is a **candidate** if it is FDA-approved (openFDA), is not topical-only and not an acute-use anaesthetic or sedative, and is *supported by at least two of A, B, C*. Literature is context and does not count as a line.
- **Tier A:** candidate and D is favourable. **Tier B:** candidate and D is unclear. **Tier C:** supported by one line only (listed, not proposed). A candidate with D unfavourable is listed as **avoid**.
- A drug whose support in B comes only from generic complexes is not supported by B. The overlap between lines is reported against what chance would give (hypergeometric test).
- Antineoplastics are listed in the tumour-context list; the others in the list for the constitutional picture.

## Age

The patient's age is unknown and the picture can change. Lines A and B are run for three **age scenarios** (0 to under 3, 3 to under 12, 12 and older): the phenotype set of each scenario is the patient's eight terms plus the phenotypes that the HPO annotations of mosaic variegated aneuploidy list with an onset
compatible with that age, weighted by their frequency, and the label support of every candidate is given for the same three groups. Lines C and D do not depend on age. The scenarios are not predictions of how the child will evolve; HPO onset terms are coarse and the annotations of this disease are few.

## Compute and reproducibility

CPU only, on a Linux x86-64 virtual machine in Google Cloud; no GPU and no training of large models. Each line is a module of `track2` with a stage in `run_all`; downloads are recorded with their source, size and checksum; all randomness is seeded. Data are used for non-commercial research; derived tables of any source whose licence restricts redistribution are not published in the repository.

## What would count as failure

- A line that fails its validation is reported with its numbers and is not used; if fewer than two lines are usable, the result is reported as a negative result and the list stays exploratory.
- Candidates whose support comes from a single data resource that all lines share (for example the same drug-target table) are flagged: the lines would not be independent.

## Implementation notes for line B (written 2026-10-07, before the proximity to the phenotype module was computed)

- **Phenotype module from HPO only.** The module uses the HPO annotations of genes (`genes_to_phenotype.txt`, release v2026-09-01) with the descendants of each patient term. The addition of genes by mouse phenotype (MGI) is **not implemented**: it needs a table that maps mouse phenotype terms (MP) to human ones (HPO), which was not downloaded. The MGI files were downloaded but are not used by the module.
- **Nine configurations for one module.** The weighting schemes of the mitosis modules do not apply to a single module. The nine configurations of the phenotype module are 3 STRING thresholds (400, 700, 900) x 3 definitions of the module: genes annotated to at least 2 terms (the pre-registered rule), to at least 3 terms, and the 100 genes with the highest information-content score. Robust = among the 20 best in at least 7 of the 9. The two additional module definitions are sensitivity analyses, not changes of the rule.
- **Drug-target links** for the phenotype proximity are the curated plus measured ones (E1), the same as in the register.
- **Mitosis module without generic complexes** is the analysis of `addendum_core.md` (`L1_core`), reused as it is.

## Amendment 1 to the validation of line B (written 2026-10-07, before it was run)

The validation as first written (the 30 diseases of the earlier benchmark, each with a module built from "the disease's own HPO terms") cannot be run: those diseases are broad (cancer, asthma, hypertension, infection...), HPO annotates mostly rare diseases under OMIM and ORPHA identifiers, and the direct match through Orphanet gives 0 diseases with at least three scorable approved drugs. The validation is replaced by the following; the thresholds are **not** changed.

- **Diseases.** Every disease of the approved-indication table (Open Targets, the same table as the earlier benchmark) whose MONDO identifier has at least one OMIM cross-reference in MONDO (`mondo.obo`, the file's own `xref: OMIM:` lines, the disease itself and not its descendants), whose OMIM entries have at least 5 distinct HPO terms in `phenotype.hpoa` (non-negated, phenotypic abnormality aspect), with at least **3 known approved drugs** that are scorable in the network, and whose module has at least 10 nodes in the network. Ascending MONDO id, at most 60 diseases.
- **Module of a disease.** The union of the HPO terms of its OMIM entries, and then exactly the rule of the patient's module (`phenotype_module.py`): genes annotated to at least 2 of the terms; fewer than 10 genes gives the 100 best by information content. The variant with the 100 best by information content is also run and reported.
- **Ground truth and metric.** Known approved drugs of the disease against all other scorable drugs; proximity at STRING 700 with the same null model and seed; ROC AUC by disease; mean, bootstrap interval over diseases and degree-matched permutation test as in `benchmark_stats.py`.
- **Decision.** The test is valid only if **at least 15 diseases** are evaluable. If it is valid, line B is *validated* if the mean AUC is >= 0.60 with a lower bound > 0.50 for the patient's rule (genes annotated to at least 2 terms). If fewer than 15 diseases are evaluable, the result is "not testable": line B is reported as exploratory and cannot support a candidate. The top-100 variant is reported but does not decide.
- **Resource.** MONDO `mondo.obo` (monarch-initiative/mondo, latest release at the time of download; CC BY 4.0). Because "latest" changes, the SHA-256 of the file used is recorded.

## Implementation notes for line D (written 2026-10-07, before it was run)

- **Network.** OmniPath interactions (`datasets=omnipath,pathwayextra,kinaseextra`, which include SIGNOR), directed rows only, nodes are UniProt accessions (complexes are dropped). An edge has sign +1 if the row is annotated as stimulation only and -1 if inhibition only; rows annotated with both are dropped.
- **Checkpoint-strength readout.** Mitotic checkpoint complex nodes BUB1B, BUB3 and MAD2L1 count +1. APC/C activity counts -1: CDC20 (as the APC/C activator, so inhibiting CDC20 strengthens the arrest, as the validation requires), ANAPC1, ANAPC2, ANAPC4, ANAPC5, ANAPC10, ANAPC11, CDC16, CDC23, CDC26 and CDC27.
- **Effect of one perturbation.** For a perturbed node and the sign of the perturbation (-1 inhibition or loss, +1 activation): if the node is a readout node the effect is direct. Otherwise, for every readout node, the signed walks of length 1, 2 and 3 from the perturbed node are counted; the effect on that readout node is the sign of the (positive minus negative) walk count at the **shortest length at which a walk exists**; zero (balanced) counts as no effect. The net effect is the sum over readout nodes of (readout sign x perturbation sign x effect). Net > 0: strengthens the checkpoint; net < 0: weakens it; 0 or no path: unclear.
- **Drug call.** For every curated mechanism of the drug whose ChEMBL action type is inhibitory (-1) or activating (+1) (the sets of `direction.py`) and whose target is a node of the network, the net effects of its targets are added. Drugs with only measured links (no action type) are *unclear*. This is a signed-path heuristic on a literature-curated network, not a pharmacodynamic model.
- **Validation.** The six perturbations of the main text, evaluated once. CDC20 inhibition and MAD2L1 loss act on readout nodes directly and are therefore recovered by construction; TTK, BUB1, AURKB and TRIP13 depend on the network. The line is used if at least five of the six are recovered.

## Implementation notes for line C (written 2026-10-07, before it was run)

- **Score and covariate.** The checkpoint-expression score is the mean z-score of the six genes of `config/prism_genes.json` (log2(TPM+1), DepMap 24Q2); the proliferation covariate is the mean z-score of its six proliferation genes. Weak = lowest quartile of the score, strong = highest quartile, among the cell lines of each data set.
- **A hit** (compound in PRISM, or gene in CRISPR) needs: BH-FDR < 0.10 over all tested rows of its data set for the score coefficient in the model adjusted for lineage and proliferation; a positive adjusted coefficient (the row is lower when the score is lower); a negative unadjusted difference weak minus strong (same sign); and |difference| >= 0.2. A row needs at least 10 measurements in each quartile.
- **Support for a drug.** A drug is supported by C if one of its PRISM rows (matched by name or synonym) is a hit, or if one of its curated targets (UniProt, mapped through the STRING nodes) is a hit gene.
- **Calibration.** 200 random six-gene scores. To make a score with the same confounding as the real one possible, the genes are drawn at random from the 1,500 genes most correlated with the proliferation score (the genes of the real score and the proliferation genes are excluded), and a set is kept if the Spearman correlation of its score with the proliferation score is within +-0.05 of the real score's. The proliferation covariate is the real one in every repetition. The statistic is the total number of hits (PRISM rows plus CRISPR genes); the line is used only if the observed total is above the 95th percentile of the null. Seed 20261007.
- **Cell lines are a proxy** and this line speaks to the tumour context only.

## Implementation notes for line A (written 2026-10-07, before its validation was run)

- **Graph.** PrimeKG (Harvard Dataverse doi:10.7910/DVN/IXA7BM, `nodes.tab` and `edges.csv`; despite its name `nodes.tab` is comma-separated; MD5 checked against the original file). Edges are made undirected and binary. **All drug-disease edges (indication, contraindication, off-label use) are removed from the graph used for ranking**; `indication` edges are kept apart as ground truth. This is stricter than hiding the indications of one disease at a time.
- **Method.** TxGNN needs DGL with CUDA and is not used. Personalised PageRank with restart probability 0.15, convergence 1e-9 (L1) or 100 iterations, uniform edge weights. The restart distribution of a disease is its own node.
- **Validation.** 200 diseases drawn at random (seed 20261007) among the disease nodes with at least 3 drugs linked by `indication`. Candidates are all drug nodes of the graph (7,957); the AUC compares the indicated drugs with the other drug nodes. Mean AUC >= 0.70 with the lower bound of the 95 % bootstrap interval over diseases > 0.55 (`benchmark_stats.summarize`, which also reports a permutation test within drug-degree bins). The validation does not use the patient's phenotypes; it tests the graph-ranking method on diseases whose indications are known.
- **Use for the patient (only if validated).** Restart nodes: the node `mosaic variegated aneuploidy syndrome` plus the PrimeKG phenotype nodes of the patient's HPO terms (PrimeKG stores the number of the HPO identifier as node id). A drug is *supported by A* if it matches (by name) a FDA-eligible drug of the register universe and is among the 50 best of those.

## Amendment 2 to line A: side-effect edges removed (written 2026-10-07, after the first patient ranking was read)

The first patient ranking (graph without drug-disease edges but with `drug_effect` edges) was dominated by drugs that **cause** the patient's phenotypes: amphotericin B, furosemide, calcitriol and ergocalciferol (nephrocalcinosis), and antidepressants and antipsychotics (weight loss, preterm birth). PrimeKG links a drug to a phenotype with the relation `drug_effect`, which is a side effect; a walk from a phenotype to the drug through that edge ranks the cause, not a treatment. This was found by reading the output, not before, and it is a post hoc change:

- The relation `drug_effect` is removed from the graph used for ranking (together with the three drug-disease relations). Nothing else changes.
- The validation (200 diseases, same seed and thresholds) is **run again on the new graph and must pass again**; if it does not, line A is not used. The result of the first run (mean AUC 0.821, 95 % interval 0.801-0.839, graph with side-effect edges) is kept in `deviations.md`. The post hoc split by degree of the disease node (AUC 0.805 for degrees 1-8, 0.819 for 9-52, 0.838 for 52-1,168) is reported as a check that the method works for poorly connected diseases such as an ultrarare one; it is not a decision rule.
- The rhabdomyosarcoma term (HP:0002859) has no node in PrimeKG (only subtypes), so the patient query uses seven of the eight terms.

## Amendment 3: candidate rule widened to "one validated line plus literature" (written 2026-10-07, after the combination of the lines was read)

The combination with the original rule (at least two usable lines among A, B, C) gave one candidate (fostamatinib, lines A and B; overlap of A and B p = 0.10, not distinguishable from chance). On the author's instruction the rule is widened; this is a **post hoc change** and is declared in `deviations.md`:

- **Tier 1 (two lines):** supported by at least two usable lines (A, B, C). Unchanged.
- **Tier 2 (one line plus literature):** supported by **exactly one usable line (A or B)** and having a **literature rationale** (below), and passing the same gates (FDA-eligible, systemic, not acute-use anaesthetic or sedative). The report states for each of them which line supports it and that the other half of the support is literature, which is context and not a measured result.
- **Literature rationale, fixed before any search.** For every drug supported by one usable line, and for every curated human target gene of the drug (`universe.json`), the PubMed query `GENE[Title/Abstract] AND (spindle assembly checkpoint OR BUBR1 OR BUB1B OR mosaic variegated aneuploidy OR chromosome missegregation OR aneuploidy)` is run with NCBI E-utilities and the count and the first 5 PMIDs are stored (`out/literature_search.json`). A drug has a rationale only if the author reads the abstract of at least one of these papers and it supports, for the drug's target, an effect on the spindle checkpoint, BubR1 function, chromosome segregation or the aneuploidy response (or on the growth and muscle phenotypes of BubR1-deficient models). The papers read and the claim each one supports are written in `config/literature_rationale.json`; the same rule applies to all drugs, and a count of zero or no supporting abstract means no rationale. The earlier curated classes (`config/literature_context.json`) are kept as they are and do not match any of these drugs.
- The age question: the HPO annotations of MVA have no onset, so the age scenarios of lines A and B cannot be built from data (see below); age enters only through the label support by age group.

## Sensitivity analyses after the independent review (written 2026-10-07, before they were run)

The review of the thresholds (see `deviations.md`) did not change any threshold. Four analyses are added; none of them is a new pass/fail threshold and each is read as declared here.

1. **Mitosis module with the two-record rule.** `core_sensitivity` is run again (same nine configurations, same seed) with the measured drug-target links restricted to those with at least two activity records (`MIN_MEASUREMENTS`, the rule already applied in the plausibility tiers) and without the similarity-predicted links (they failed their out-of-sample check). Output: `out/core_n2/`. Question: is fostamatinib still among the 20 best in at least 7 of the 9 configurations? The answer is reported whichever way it falls; if it is not, the support of line B for fostamatinib does not stand and the candidate list is empty under the original rule.
2. **Line C with 1,000 null sets.** The calibration is repeated with 1,000 matched random scores (the same procedure and tolerance; seeds continue the same generator) and the per-set counts of PRISM and CRISPR are stored. **The estimate of the 95th percentile from 1,000 sets replaces the one from 200 sets, whichever way the decision goes** (observed total > 95th percentile means used). As a sensitivity analysis only, the observed count with the six MDM2 inhibitors collapsed to one hit is also reported.
3. **Line A, what is actually used.** On the same 200 diseases: (a) recall of the indicated drugs in the top 50 (all drug nodes, and only the drugs matched to the register universe); (b) the AUC of a degree-only baseline (score = degree of the drug) and the per-disease excess of the PageRank AUC over it, with a bootstrap interval over diseases; (c) the AUC and the recall when the restart is the disease node plus its own phenotype nodes (PrimeKG `disease_phenotype_positive` neighbours), which is the nearest to the patient query. They are limitations to be reported; the pre-registered pass rule of line A stays as it was.
4. **Report.** All of it is stated in the report: A and B valid with a modest excess over degree, at most one candidate, overlap compatible with chance.

Amendment to item 2 above (still before the run): the collapse is generalised from MDM2 inhibitors to **every mechanism class of the PRISM annotation (`MOA` column)**: the PRISM hits of a class count as one, and a hit without annotation counts as its own class. The same collapse is applied to every null set, and the observed class-level total is compared with the 95th percentile of the class-level null. This is a sensitivity analysis; the decision of line C stays on the uncollapsed total.

## Result of the extended check of line A (2026-10-07)

On the 200 validation diseases: mean AUC 0.802; degree-only baseline 0.780; mean excess 0.022 (bootstrap interval 0.014 to 0.030; PageRank above the baseline in 137 of 200 diseases); recall of the indicated drugs in the top 50 among the drugs of the register universe 0.099 against 0.026 for the degree baseline. Adding the disease's own phenotype nodes to the restart gives AUC 0.801 (excess 0.020), i.e. no gain. Reported in the report as the headline numbers of line A.

## Three further analyses, fixed before they are run (2026-10-07)

None of them is a new pass/fail threshold; each is a limitation analysis or a test of specificity.

1. **Track 1 scenarios (`track1_scenarios.py`).** Using the stored threshold scores of the mitosis modules (original links and the two-record variant without predicted links), for each weighting scheme separately (W1 weights from Track 1, W2 equal weights, W3 the pathogenic gene only, which is the single-hit-like scenario) a drug is *scheme-robust* if it is among the 20 best at all three STRING thresholds. Reported: the scheme-robust set of each scheme, and the drugs that are scheme-robust in all three (those do not depend on how much weight the secondary genes CENPE, MAD1L1 and IGBP1 get).
2. **Specificity of the patient query in line A (`kg_nullquery.py`).** 500 random queries, each a disease node whose degree (after edge removal) is within 30 % of the real one (115) plus seven phenotype nodes, one per real phenotype node, with degree within 30 % of it; same restart weights and graph. For fostamatinib and for the 50 drugs of the real top 50, the rank in each random query is stored. Reported: how often a random query puts fostamatinib in its top 5 and top 50 (empirical p), and how many of the real top-50 drugs are in the top 50 of at least 5 % of the random queries (a measure of how much of the real list is generic).
3. **Validation of the module rule of the mitosis lines (`module_rule_validation.py`).** For the diseases of the validation of line B whose causal genes are known (genes annotated in `genes_to_phenotype.txt` to the disease's OMIM entries), the module is built with the same rule used for the Track 1 genes (`auto_modules`: Reactome pathways of the gene, 10 to 300 participants; STRING neighbours as fallback), the proximity of all scorable drugs is computed at STRING 700 with the two-record measured links and no predicted links, and the AUC of the known approved drugs is computed. Same reading as the validation of line B (mean AUC >= 0.60 with the lower bound > 0.50 and at least 15 evaluable diseases); if fewer than 15 diseases have usable modules the result is "not testable". It decides only whether the mitosis lines are described as validated or as unvalidated in the report.
