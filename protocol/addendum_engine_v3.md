# Addendum - evidence engine v3: error-correcting development with a sealed held-out set

Written on 2026-10-15, **before any v3 model, feature, fit or evaluation exists** (no v3 code was run; `src/track2/engine_v3.py` is written after this file). It extends, and does not replace, `addendum_engine.md` (v1, pre-registered primary result: supported on 17 diseases) and `addendum_engine_v2.md` (v2: criterion 1 fails on the 28 new diseases). Nothing here modifies `frozen_protocol.md`, `engine_score.py`, `engine_score_v2.py`, `engine_case_study.py`, `benchmark_v2.json`, `validation.json`, `validation_v2.json` or `case_wilms.*`.

## 0. Status and honesty clauses

- This is **post hoc development, not a pre-registered primary result**. The v1 decision (supported, 17 diseases) and the v2 result (criterion 1 not met on the new diseases) stand unchanged. Engine v3 is built because v2 and the Wilms worked case showed errors in features and weights; whether the corrections generalise is decided only by the held-out set below.
- **Kidney Wilms tumor was inspected before this protocol** (`case_wilms.md`: the known drugs are missed, the top ranks are hormonal noise such as tamoxifen, progesterone, danazol). It therefore belongs to DEV, and its improvement is a **diagnostic, not evidence**.
- Failure to improve is reported as is. If the held-out set shows no improvement, the report says "no demonstrated improvement" and v3 is not recommended over v1.
- The v1 / v2 per-disease and aggregate results were seen when v3 was designed. The candidate families below are motivated by those observations; this is why only the held-out set, never touched during design, tuning or selection, can support a claim.

## 1. Split (fixed now, seed 20261015)

Source: the 45 diseases of `out/engine/benchmark_v2.json`. Rule:
1. DEV contains all 17 v1 diseases (`origin = v1`) and kidney Wilms tumor.
2. The remaining 27 `v2_new` diseases are sorted alphabetically (Python `sorted`), shuffled with `random.Random(20261015).shuffle`, and the first `floor(27/2) = 13` go to DEV; the other 14 are HELD-OUT.
3. Resulting sizes: **DEV 31, HELD-OUT 14** (at least 12 as required). The code recomputes the split with this rule and aborts if the lists below differ.

DEV additions (13): X-linked adrenoleukodystrophy; absence epilepsy; autoimmune thrombocytopenic; diffuse gastric adenocarcinoma; hereditary angioedema; hereditary breast ovarian cancer syndrome; hypogonadotropic hypogonadism with or without anosmia; ovarian cancer; patent ductus arteriosus; small cell lung carcinoma; squamous cell carcinoma; temporal arteritis; testicular germ cell tumor.

HELD-OUT (14): Alzheimer disease; Parkinson disease; familial thyroid carcinoma; hyperparathyroidism; hypoalphalipoproteinemia; hypoparathyroidism; multiple endocrine neoplasia; myelomonocytic leukemia; narcolepsy; orthostatic hypotension; prolactin producing pituitary gland tumor; systemic lupus erythematosus (disease); thyroid gland carcinoma; von Willebrand disease.

HELD-OUT diseases are **never used** for feature design, tuning, model selection, hyper-parameters, the Platt map or any decision. The development command does not compute their scores. The held-out command refuses to run unless the frozen file exists and its hash verifies, and refuses a second run (a lock file records the first).

## 2. Unchanged endpoints

The three criteria of `addendum_engine.md` section 4, evaluated on arm A-type scoring (prior and popularity features available), gold 1, candidate set = all drugs of `data/universe.json`:
1. macro AUC minus the AUC of the literature-popularity baseline >= 0.05, with the paired 95 % disease-level bootstrap interval (2000 resamples, seed 20261015, percentile) lower bound > 0;
2. Platt calibration slope in [0.8, 1.2] (pooled over pairs, logit of the calibrated probability);
3. macro precision at 10 above the best of the three baselines (popularity, proxy, random; best by macro value).
Baselines (popularity, proxy, random), tie handling, metrics, bootstrap and the `decide` function are the v1 ones (imported, not modified). Secondary: recall at 20 and Brier score. No threshold is changed.

## 3. Candidate improvements (closed list; nothing else will be tried)

Base = the 8 engine-v1 features of arm A (`engine_score.FEATURE_NAMES`). Notation per disease and drug: c = papers with drug and a disease gene (v1 `n_papers`), nb = benefit papers, nh = human papers, pop = Europe PMC record count of the drug, N = total Europe PMC records, n_g = papers with a gene symbol, E = n_g * pop / N (v1 expectation).

**Family F1 - specificity normalisation against global drug frequency (4 features).** idf = ln(N / (pop + 1)).
- `f1_tfidf_papers` = ln(1 + c) * idf; `f1_tfidf_benefit` = ln(1 + nb) * idf; `f1_tfidf_human` = ln(1 + nh) * idf;
- `f1_pmi_benefit` = ln((nb + 0.5) / (E * rho + 0.5)), rho = sum over drugs of nb / sum over drugs of c for that disease (0 if the sum of c is 0).

**Family F2 - mechanism-class-aware features (3 features).** The class of a disease is the plurality `mechanism_class_general` among its genes, ignoring `other`, ties by alphabetical order; with no non-`other` gene all three features are 0. Class genes G = the disease genes of that class.
- `f2_class_gene_papers` = ln(1 + cg), cg = papers (of the same retrieved set, v1 sentence-level drug match) that mention the drug and, case-sensitively as whole words, a symbol of G;
- `f2_class_gene_pmi` = ln((cg + 0.5) / (E_cg + 0.5)), E_cg = n_cg * pop / N, n_cg = papers mentioning a symbol of G;
- `f2_class_keyword_papers` = ln(1 + kw), kw = papers with a sentence that contains the drug and a keyword of the class (case-insensitive, word-start match). Keywords (fixed here, no gene or disease names):
  lysosomal_vesicular_trafficking: lysosom, vesicl, endosom, autophag, traffick, storage, exocyt, endocyt; ion_channel_transporter: channel, transporter, sodium, potassium, calcium, chloride, electrolyte; metabolic_enzyme: enzyme, metabol, substrate, synthase, dehydrogenase, hydroxylase, cofactor; receptor_signalling: receptor, signalling, signaling, kinase, agonist, antagonist, inhibitor, pathway; transcription_chromatin_regulation: transcription, chromatin, histone, methylat, epigenet, gene expression; structural_cytoskeletal: collagen, cytoskelet, actin, myosin, muscle, fibrosis, cardiomyopath; mitochondrial: mitochond, oxidative phosphorylation, respiratory chain, oxidative stress; dna_repair_replication: DNA repair, DNA damage, replication, genotoxic, crosslink, PARP; cell_cycle_mitotic_checkpoint: cell cycle, mitosis, mitotic, checkpoint, proliferat, apoptosis; rna_processing_translation: RNA, splic, translation, ribosom; extracellular_secreted: secreted, extracellular, matrix, cytokine, complement, growth factor, hormone; other: none.

**Family F3 - drug-class plausibility (1 feature).** `f3_atc_class_rate` = ln(rate_c / p0), the ATC level-3 class c of the drug (the v1 class definition; drugs without ATC are their own class), rate_c = (G_c + 20 * p0) / (M_c + 20), where G_c = number of gold 1 pairs in class c and M_c = number of (disease, drug) pairs in class c over the **reference diseases**, p0 = overall gold 1 pair fraction in the reference diseases. Reference diseases are **only the training diseases of the fold being fitted**; for a training row the row's own disease is left out of the reference (leave-one-disease-out, so a label never encodes itself); for a test row all training diseases are the reference. No held-out or test-fold label is used.

**Family F4 - L2 regularisation strength (not a feature).** Logistic regression with standardised features (standardisation fitted on training rows) and L2 penalty C in {0.01, 0.03, 0.1, 0.3, 1}, chosen **inside every training set by nested grouped CV** (disease-grouped inner folds, seed 20261016, v1 fold maker, same ATC class masking as v1): highest inner macro AUC, then the **1-SE rule toward stronger regularisation** (smallest C whose inner mean AUC >= best - SE of the best, SE = sd / sqrt(number of diseases) of its per-disease AUCs). v1 used the fixed C = 1.

**Configurations.** Base features plus every subset of {F1, F2, F3}: 8 configurations (`base`, `F1`, `F2`, `F3`, `F1+F2`, `F1+F3`, `F2+F3`, `F1+F2+F3`), each with nested C selection (F4). Complexity order = number of families, then number of features.

## 4. Selection rule

Outer cross-validation on DEV only: 5 folds by disease (v1 fold maker, seed 20261015), ATC-class masking of training rows as in `addendum_engine.md` section 6, nested C selection inside each outer training set, F3 reference = outer training diseases. Statistic = macro mean of per-disease AUC of the out-of-fold scores vs gold 1. **Selected configuration = the simplest configuration whose DEV macro AUC is >= (best DEV macro AUC - its standard error)**, SE = sd(per-disease AUC, ddof 1) / sqrt(31) of the best configuration; ties in complexity are broken by higher AUC. The `base` configuration with nested C is itself a candidate, so "no feature family helps" is a possible outcome. The v1 engine (base features, C = 1) is also scored on DEV as a reference only; it is not a candidate.

## 5. Freeze and the single held-out scoring

1. After selection, the final C is chosen by the same nested procedure on all 31 DEV diseases (inner folds seed 20261017), the model is refit on all DEV diseases (no ATC masking; F3 reference = all DEV diseases), and a Platt map is fitted on the pooled outer-CV out-of-fold raw scores of the selected configuration on DEV.
2. Everything needed to score a new disease (configuration, feature list, C, standardisation, coefficients, intercept, F3 table, Platt map, DEV and HELD-OUT name lists, split rule) is written to `out/engine/engine_v3_frozen.json` and its SHA-256 is printed and stored. The same is done for the v1 reference weights (base features, C = 1, fitted on all DEV, Platt from DEV out-of-fold scores of the v1 pipeline). Both hashes are appended to this file in section 8 before HELD-OUT is scored.
3. HELD-OUT is scored **exactly once** with the frozen files: engine v1-weights, engine v3, popularity, proxy, random; the three criteria for v3 (and, for information, for v1-weights); the paired AUC difference v3 minus v1-weights with 95 % disease-level bootstrap. No ATC masking is applied because the weights never read a HELD-OUT label (stricter than the v1 class-split rule).
4. **Rule for "improvement"**: v3 is reported as improved over v1-weights only if the 95 % interval of the paired held-out AUC difference (v3 minus v1-weights) lies above 0; otherwise the report states "no demonstrated improvement". Meeting the three criteria is reported separately and verbatim.

## 6. Ablation on DEV

(a) Add-one: the DEV results of `base`, `F1`, `F2`, `F3` from the configuration grid. (b) Drop-one: for the selected configuration, each included family removed in turn, with the identical nested-CV procedure. (c) Effect of F4: the selected configuration with fixed C = 1 instead of the nested choice. Reported with DEV macro AUC and the paired difference to the selected configuration; none is used to change the selection.

## 7. Wilms diagnostic (diagnostic, not evidence)

Using the selected configuration with its frozen C, a model is fitted on the DEV diseases **excluding Wilms** (no Wilms label read; strict refit), applied to Wilms, and the ranks (1 = best, among all universe drugs, ties given the lowest rank of the tie group) of doxorubicin, dactinomycin, vincristine, tamoxifen, progesterone and danazol are reported next to the ranks under the v1 features and C = 1 fitted on the same training diseases.

## 8. Reporting and deviations

Results go to `out/engine/validation_v3.json` (sorted keys). All numbers are computed by `src/track2/engine_v3.py`. Any deviation from this text is listed in the output and in the final report. Limits fixed in advance: DEV has 31 diseases, HELD-OUT 14, so intervals are wide; gold 1 is incomplete and partly circular; the held-out diseases include several common (non-rare) diseases; literature features are abstract co-mentions.

### Freeze record
(appended after the freeze, before HELD-OUT scoring)

Freeze record written on 2026-10-15 after the DEV run and **before any HELD-OUT score exists**:
- Selected configuration: `F3` (base features plus the ATC-class plausibility feature), nested-CV C = 0.01. DEV selection: best mean `F1+F3`, threshold = best - SE, eligible `F1+F2+F3`, `F1+F3`, `F2+F3`, `F3`; the simplest eligible is `F3`.
- `out/engine/engine_v3_frozen.json`: SHA-256 of the model body `cc0a54e6afb8c3fbf540ac93963d4cc821cea68cfc96d7276a8ed9568adb579a`.
- `out/engine/engine_v1_weights_frozen.json` (v1 reference weights, base features, C = 1, fitted on DEV): `b4c63900b0f1e57d8df7f07607199b72eef4a0e4b9035868ee15a8ce40a7d65e`.
- The held-out command verifies both hashes, writes a lock file `out/engine/heldout_v3_scored.lock` and refuses a second run.
