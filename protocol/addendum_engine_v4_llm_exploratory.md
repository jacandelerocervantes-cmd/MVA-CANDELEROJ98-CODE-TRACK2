# Addendum - engine v4, LLM-read evidence: EXPLORATORY, DEV-ONLY, POST HOC analysis

Written on 2026-10-16, **before** `src/track2/engine_v4_llm_exploratory.py` was run. It extends and does not modify `addendum_engine_v4_llm.md`, `engine_v4_llm.py`, any label file or any earlier output.

## 0. Status (read first)

- **EXPLORATORY, DEV-ONLY, POST HOC.** It cannot confirm anything. It is a hypothesis-generating look, run after the pre-registered experiment closed.
- **Nothing here changes the pre-registered outcome.** With guideline v1 the DEV relevance kappa between the two annotators is 0.565, below the pre-registered gate of 0.6. Therefore the pre-registered experiment is closed: the relevance labels are unreliable, no HELD-OUT scoring is run, and the hypothesis is reported as "not testable with these labels". The direction kappa (0.841) is reliable. This addendum does not reopen the gate, does not re-score anything pre-registered, and does not touch the decision rule of section 7 of the parent addendum.
- **The sealed HELD-OUT is not touched.** No HELD-OUT label, HELD-OUT gold or HELD-OUT annotation file is read or computed on; `engine_v4_llm heldout` is not run; no lock file is created. The code reads the candidate file, the pair index and the benchmark file (which contain both splits) only to select DEV diseases, and discards the rest at once.
- Question: does LLM-read evidence carry signal **beyond the PrimeKG prior** at low literature volume (thinned fractions)?
- Multiplicity: three features x four fractions are compared and no multiplicity correction is applied. No threshold of success is set. Intervals are descriptive.

## 1. Features (all declared here; no other feature is computed)

Unit: a (DEV disease, candidate drug) pair of the candidate sets P_d of the parent addendum. A "counted sentence" is a verified sentence (annotated item) of one of the first 3 kept papers of the pair (kept = papers retained by the sample), at most 2 sentences per paper (the parent design, so every count is capped at 6). The model input is `log1p(count)`; a pair with no annotated sentence has value 0.

Direction used by all three features: the consensus of the parent `consensus()`: the label shared by at least 2 of (Opus, Sonnet, rule extractor); this includes the case "both annotators agree". `unresolved` sentences are never counted. Only the direction label is used for F_dir.

- **F_dir**: number of counted sentences whose consensus direction is `benefit`. Relevance labels are not used.
- **F_rel_strict**: number of counted sentences whose consensus direction is `benefit` AND both annotators say `relevance = yes`.
- **F_rel_any**: number of counted sentences whose consensus direction is `benefit` AND at least one annotator says `relevance = yes`.

F_dir is the only feature built on a reliable label (direction kappa 0.841); the two relevance features rest on labels that failed the gate (kappa 0.565) and are reported as such.

## 2. Models

- **PO** (comparator): the prior group only (`log_popularity` and the two PrimeKG columns, columns 4, 6, 7 of the B0 matrix), L2 logistic (`engine_v4.fit_head`), refit by the parent nested procedure.
- **PO + feature**: the same three columns plus one of the features above, one model per feature (three models).
- Fit and selection as the parent: training on the f = 1.0 (capped) features, ATC-class masking of training rows, nested grouped CV on DEV only (outer: `engine_v4.group_folds`, a group of identical gold sets is never split; inner grouped 3-fold chooses C in {0.1, 1, 10} by mean macro AUC). The scores used are out-of-fold: each DEV disease is scored by the head of its outer fold. Group = the engine_v4 group.
- Evaluation universe: the candidate set P_d (parent, section 5), so that gold membership is not encoded by "has an annotated sentence". Only paired differences are interpreted; levels are not comparable to full-universe AUCs.

## 3. Fractions and seeds

f in {1.0, 0.1, 0.03, 0.01}. The thinning is exactly the dilution/ablation subsample: `default_rng([20261015, rank, frac_index, rep])`, 5 replicates for f < 1 and one for f = 1.0, `engine_v4_dilution.subsample_indices`. The per-disease metric is the mean over replicates. (The grid of the dilution module also contains 0.3; it is computed as a by-product of the shared table code and is not reported.)

## 4. Metric

Paired macro AUC difference (PO + feature) minus PO, per feature and fraction, with a **95 % percentile bootstrap resampling the engine_v4 groups** (2000 resamples, seed 20261015, same indices for both models; `engine_v4_ablation.table`). The disease-level bootstrap is stored as a sensitivity. Precision at 10 differences are stored but are secondary.

## 5. Coverage

For each feature and fraction: the number (and share) of DEV candidate pairs with a non-zero feature value, as the mean over replicates (with minimum and maximum over replicates), overall and split into gold-1 positives and negatives. This shows whether a feature has any coverage at low volume; if coverage is near zero, a null difference carries no information about signal.

## 6. Reading rule (declared before the run)

A positive interval is a hint worth a future pre-registered test with reliable relevance labels, not evidence. A null interval at low coverage is uninformative. Nothing is adopted into the engine. All constants are inherited from the parent and stay UNVALIDATED with their conservative defaults (parent registry, section 9); the new choices (feature definitions, "either annotator" rule for F_rel_any, one extra column per model) are UNVALIDATED and registered in the output.

## 7. Leakage audit

Labels of the annotators are independent of the gold. Features have no gold argument (signature tested). Candidate sets P_d include DEV gold and the top-30 B0 negatives (parent S1); the same set is used for both models. Annotators know what drugs treat which disease (parent S4); this affects F_dir too, so a positive result could reflect parametric knowledge rather than literature reading.
