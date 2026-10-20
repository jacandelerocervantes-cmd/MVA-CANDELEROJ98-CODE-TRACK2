# Addendum - engine v4, prior / literature ablation of the dilution result

Written on 2026-10-15, **before any ablation score, feature or number exists**. Nothing from `validation_v4_dilution.json` is used to set any choice below. This file extends `addendum_engine_v4.md` and `addendum_engine_v4_dilution.md` and modifies none of them, nor `engine_v4.py`, `engine_v4_dilution.py`, `benchmark_v4.json`, `validation_v4*.json`, `engine_v4_frozen.json` or any earlier file.

## 0. Status and honesty clauses

- **Post hoc development experiment, not a pre-registered primary result.** The v4 decision (sealed HELD-OUT) and the dilution result stand as they are.
- **Why.** In the dilution experiment only the paper-derived features were thinned; `log_popularity` and the two PrimeKG prior columns stayed fixed. The residual gap of B0 over popularity at 1 % of the literature (f = 0.01, 5 papers) may therefore come from the prior and not from the literature. This experiment separates the two sources.
- **Any failure is reported as is** (a gap whose lower bound never exceeds 0, a non-monotone curve, a decision-rule outcome that is neither of the two hypotheses).
- Every number is produced by `src/track2/engine_v4_ablation.py`. Seed everywhere: **20261015**.
- The HELD-OUT set is scored **once** for all ablations, after everything in this file is fixed, with a lock file (section 6). The HELD-OUT labels of benchmark_v4 were already used for B0 / E (v4) and for the dilution run; this experiment re-uses the same sealed diseases, so it is a further look at an already opened set, not a fresh confirmation. Its role is to decompose a known result, not to select a model.
- Limitations stated now: random thinning of papers is a stress test of volume dependence (as in the dilution protocol); the frozen ablations replace features by constants, which is a counterfactual evaluation of fixed weights, not a model fitted for that feature set.

## 1. Feature groups (B0 columns, `engine_score.FEATURE_NAMES`)

- **Paper-derived (literature) group, columns 0, 1, 2, 3, 5**: `log1p_n_papers`, `log1p_n_benefit`, `log1p_n_human`, `log1p_n_harm`, `excess_over_popularity` (the last is computed from the number of papers that contain a disease gene symbol and from the paper count of the drug, so it is paper-derived and is thinned in the dilution experiment).
- **Fixed (prior) group, columns 4, 6, 7**: `log_popularity` (global drug popularity, independent of the disease's papers), `prior_direct_targets`, `prior_ppi_targets` (PrimeKG).

The "DEV mean" of a column is its mean over all (DEV disease, drug) rows of the v4 feature cache at f = 1.0 (no class masking); this equals the standardisation mean `mu` stored in the frozen B0 head, which the code checks (maximum absolute difference reported).

## 2. The ablations

1. **PRIOR-ONLY (frozen).** The frozen B0 weights (`engine_v4_frozen.json`, hashes verified exactly as in the dilution protocol); columns 0, 1, 2, 3, 5 set to their DEV mean for every drug of every disease at every fraction. The score then depends only on columns 4, 6, 7. **No refit.** Because the literature columns are constant, the result is by construction the same at every fraction and every replicate; it is still computed through the same path and the invariance is checked and reported.
2. **LITERATURE-ONLY (frozen).** The frozen B0 weights; columns 4, 6, 7 set to their DEV mean. Literature columns come from the diluted papers of each fraction as in the dilution protocol. **No refit.**
3. **Refit-on-DEV reduced models (fair comparison).** L2 logistic regression, C = 1, standardised features, ATC level-3 class masking of training rows, exactly as B0 (`engine_v4.fit_rows`), on only the columns of the reduced set: **PRIOR-ONLY-refit** (columns 4, 6, 7) and **LITERATURE-ONLY-refit** (columns 0, 1, 2, 3, 5). Training uses the full-literature (f = 1.0) features of the training diseases; the fitted model is then applied to the thinned features of other diseases. Estimates:
   - DEV (honest): **grouped 5-fold cross-validation** with the v4 folds (`engine_v4.group_folds`; a group is never split), each fold's model fitted on the other folds only, applied to the held-out fold's features at every fraction. Platt calibration is not needed for any endpoint here (AUC, precision at 10 and recall at 20 are rank-based and invariant to a monotone map), so the nested Platt layer of `engine_v4.oof_head` is not repeated; the code checks that its f = 1.0 out-of-fold scores equal those of `engine_v4.oof_head(..., nested_platt=True)` (maximum absolute difference reported). This is the only deviation from "nested" CV and it cannot change any reported number.
   - HELD-OUT: the model refitted on **all DEV** diseases (no masking, as the frozen B0), applied to HELD-OUT.
   - The frozen ablations 1 and 2 and the frozen B0 are in-sample for the weights on DEV (as stated in the dilution protocol); the refit models are out-of-fold on DEV. The primary reading is HELD-OUT.

B0 itself (frozen, all 8 columns) is re-scored with the same code as the reference row, and its numbers are checked against `validation_v4_dilution.json` (maximum absolute difference of the macro AUC per fraction reported; a non-zero difference is reported, not hidden).

## 3. Design (identical to the dilution protocol)

Fractions f in {1.0, 0.3, 0.1, 0.03, 0.01}; floor of 5 papers; 5 replicates for f < 1 (f = 1.0 single deterministic replicate); the same subsamples, i.e. the thinned feature cache `data/cache_engine/v4_dilution/dilution_features.pkl` produced under `addendum_engine_v4_dilution.md` (seed 20261015, `default_rng([20261015, rank, frac_index, rep])`); 145 diseases, 103 DEV and 42 HELD-OUT; cache only, no network, no embeddings (B0 feature set only). Endpoints per fraction: macro AUC; **AUC gap vs popularity** with a 95 % group-level bootstrap interval (2000 resamples, percentile, seed 20261015, resampling the disease groups, same index matrix as the dilution run); precision at 10 and recall at 20 (macro, same bootstrap); the plain disease-level bootstrap as sensitivity for the AUC gap; mean `n_sub`.

**Paired differences** per fraction, same bootstrap indices: B0 minus PRIOR-ONLY (primary, requested), B0 minus LITERATURE-ONLY, and LITERATURE-ONLY minus PRIOR-ONLY, in macro AUC (also P@10), each with 95 % CI; and the same for the refit pair.

The OLS gap-versus-log10(papers) curve and break-even volume of the dilution protocol are also reported for LITERATURE-ONLY and the refit models (descriptive; `engine_v4_dilution.fit_gap_curve`).

## 4. Decision rules (fixed now, evaluated on HELD-OUT only)

Primary models are the **frozen** ablations; the refit models are the sensitivity and are read only to say whether they agree.

- **H-LIT: "the literature features contribute at low volume"** holds if and only if, for LITERATURE-ONLY (frozen) at **both** f = 0.03 and f = 0.01 (f <= 0.03), the HELD-OUT AUC gap vs popularity has a 95 % group-bootstrap **lower bound above 0**. If only one of the two fractions satisfies it, H-LIT does not hold (reported as "partial").
- **H-PRIOR: "the low-volume gap is explained by the prior"** holds if and only if, at f = 0.01 on HELD-OUT, the PRIOR-ONLY (frozen) AUC gap vs popularity **lies inside the 95 % CI of the B0 AUC gap** at f = 0.01 (point estimate of PRIOR-ONLY within [lower, upper] of B0). Corroboration reported next to it (not part of the rule): the paired difference B0 minus PRIOR-ONLY at f = 0.01 with its CI, and whether that CI contains 0.
- The two are **not mutually exclusive**; the four combinations are reported as such (both: the prior carries the gap but the literature still adds something; H-PRIOR only: no evidence that the literature contributes; H-LIT only: the gap needs the literature; neither: reported as inconclusive). Agreement or disagreement of the refit models is stated. No threshold on the gap size is used to accept or reject the engine.
- Directions anticipated before running (a result opposite to them is reported as such): (A1) PRIOR-ONLY retains a positive gap over popularity (the prior adds information beyond popularity); (A2) the B0 minus PRIOR-ONLY difference shrinks as f decreases; (A3) LITERATURE-ONLY (frozen) at f = 1.0 is below B0 because its weights were fitted jointly with the prior.

## 5. UNVALIDATED values (status and conservative default; none promoted by this run)

Which columns form the literature group (0, 1, 2, 3, 5) and the prior group (4, 6, 7); mean replacement as the removal operation; the interval-inclusion form of H-PRIOR and the both-fractions form of H-LIT; the use of the DEV mean (not zero or the HELD-OUT mean); the refit models as sensitivity only. A test fails if any lacks a status or default.

## 6. Single HELD-OUT use and lock

1. Nothing in this file needs the HELD-OUT labels to be built; the thinned features exist and are label-free.
2. The DEV command scores everything on DEV (frozen ablations; refit models by grouped CV), writes the DEV part of `out/engine/validation_v4_ablation.json`, and reads no HELD-OUT label. After it, no design element of this file changes.
3. The HELD-OUT command verifies the hashes (frozen file SHA-256, DEV and HELD-OUT id hashes), **creates `out/engine/heldout_v4_ablation_scored.lock` before reading any HELD-OUT label**, scores B0 and the four ablations once, and refuses to run if the lock exists. The decision rules of section 4 are evaluated by code and written to the JSON.
4. Output `out/engine/validation_v4_ablation.json` is written with sorted keys.

## 7. Run log and deviations
(appended after the runs)

- DEV run, then HELD-OUT run once (lock `out/engine/heldout_v4_ablation_scored.lock`, hashes verified). Checks: DEV mean vs frozen `mu` max abs diff 8.8e-13; B0 re-scored here reproduces `validation_v4_dilution.json` exactly (max abs diff 0.0, DEV and HELD-OUT); our f = 1.0 out-of-fold refit scores equal `engine_v4.oof_head` (diff 0.0); PRIOR-ONLY AUC is identical across fractions (range 0.0).
- Deviations from sections 1-6: none. Observation, not a design change: column 5 (`excess_over_popularity`) is paper-derived by the section 1 definition but is a decreasing function of drug popularity when the drug has few papers, so with thinned papers the LITERATURE-ONLY score ranks roughly against popularity and its gap is strongly negative; LITERATURE-ONLY (frozen) and its refit give almost identical rankings for that reason. H-LIT is evaluated as written. The H-PRIOR rule holds, while the corroborating paired difference B0 minus PRIOR-ONLY at f = 0.01 is small but its CI excludes 0.
