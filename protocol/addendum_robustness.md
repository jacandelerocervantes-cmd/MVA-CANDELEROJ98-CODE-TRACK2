# Addendum 2: robustness of the submission-derived ranking

Written on 2026-10-06 **before** the analyses it describes are run. It adds three analyses; it does not change the frozen
protocol (`frozen_protocol.md`) or `addendum_submission.md`.

## Why

A review of the first submission-derived result found three weaknesses: (1) the validation tests P1 to P4 were run on the hand-defined
checkpoint module, not on the modules actually used for the ranking; (2) the gene weights are the Track 1 `epcr` values, which are
ranking values and not calibrated probabilities, and BUB1B dominates the sum; (3) the ranking changes with the STRING threshold
(top-20 Jaccard 0.48 and 0.38 against the 700 threshold).

## Analysis A: validation on the submission-derived modules

P1 (leave-one-out recovery of seed genes), P3 (STRING thresholds 400/900, dropping 20 % of the seeds 20 times) and P4 (degree bias)
are run, with the decision rules of the frozen protocol unchanged (AUC >= 0.70 and lower bound > 0.50; unstable if mean Jaccard < 0.5),
for each modellable gene module of `config/mechanisms_from_submission/`. P2 (Fanconi anemia) does not depend on the module and is not rerun.
For a module with more than 120 seeds, P1 uses a fixed random subsample of 120 seeds (seed 101) to keep the run time practical; this is the only change.

## Analysis B: weighting schemes

The combined score is computed under three schemes: (W1) `epcr` weights (the main scheme), (W2) equal weights for every modellable gene,
(W3) BUB1B only (the one gene with a ClinVar pathogenic variant). Rank agreement between schemes is reported (Spearman and top-20 Jaccard).

## Analysis C: consensus across thresholds and weights

Nine configurations: STRING 400, 700, 900 x W1, W2, W3. For each drug the **robustness count** is the number of configurations in which it is
among the 20 best of the combined ranking (drugs not scorable at a threshold count as absent there).

- **Robust candidate:** count >= 7 of 9.
- **Fragile:** count <= 3 of 9.
- Everything else: **intermediate**.

The report's main list is the robust set; fragile drugs are reported but labelled as sensitive to modelling choices. No threshold in this
rule is changed after seeing the counts. The number of robust candidates is reported whatever it is (including zero).

## What would count as failure

- No drug is robust: the ranking is then reported as not stable and the drug list is presented as exploratory only.
- P1 not supported for the BUB1B module: reported as such, with the same wording as for the manual module.
