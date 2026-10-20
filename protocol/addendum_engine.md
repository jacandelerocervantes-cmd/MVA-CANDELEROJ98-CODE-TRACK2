# Addendum - validation of the literature-grounded evidence engine (pre-registered)

Written on 2026-10-15, **before any engine scorer, popularity count or baseline score exists**. At this point only the benchmark has been built (`src/track2/engine_benchmark.py`, output `out/engine/benchmark.json`, SHA-256 `48c53a059a3f35218a3fa57efaf26d8d28d4002d86dde536ab4432211dde9a74`). Nothing here modifies `frozen_protocol.md`. All thresholds below are fixed now and **must not be changed after seeing results**.

## 1. Benchmark (already built)

- Target 30 rare genetic diseases; **17 qualified** (`n_benchmark` in the file). The shortfall is reported, not corrected: the file's `selection_log` lists every evaluated candidate (840) with the reason for exclusion. No rule was relaxed to reach 30.
- Eligibility, seed 20261015: at least 2 genes (HPO and/or Orphanet, via the same sources and mechanism rules as `disease_landscape`); at least 2 approved drugs of `data/universe.json` with a PrimeKG `indication` edge (gold 1); rare genetic = an Orphanet-catalogued disorder of exactly the disease name with Mendelian/mitochondrial/sex-linked inheritance (curated candidates of `empirical_literature_control.py` may pass with any matched Orphanet-catalogued Mendelian disorder); MVA and Fanconi anaemia excluded (application cases); a disease whose gold 1 contains an enzyme, gene or cell therapy is excluded; diseases with an identical gold 1 set are collapsed to one.
- Gold 1 (primary): universe drugs with a PrimeKG `indication` edge to the disease.
- Gold 2 (temporal): universe drugs named as an intervention of an interventional ClinicalTrials.gov study of the disease, phase 2 or later, start date after 2018-12-31. Evidence cutoff year **Y = 2018**.
- Candidate set of every ranking: all drugs of `data/universe.json` (approved drugs), for every disease.

## 2. Primary endpoint

Macro AUC over benchmark diseases of the engine score versus gold 1 (Mann-Whitney, ties = 0.5; positives = gold 1 drugs, negatives = all other universe drugs; a disease enters the macro mean only if it has at least one positive and one negative). Uncertainty: **disease-level bootstrap, 2000 resamples of diseases with replacement, percentile 95 % interval**, seed 20261015. Differences against a baseline use the same resamples (paired).

## 3. Baselines (identical candidate set and gold)

1. **Literature popularity**: total number of Europe PMC records of the drug (hit count of the exact-phrase query on the drug name, all years; disease-independent score).
2. **Proxy score** of `src/track2/empirical_literature_control.py` (sentence-level co-mention count of drug, benefit verb and mechanism term, using the same mechanism terms as the engine).
3. **Random**: uniformly random scores (expected AUC 0.5; precision at k by exact tie-averaging, as in the proxy module).

## 4. Success thresholds (all three required to declare "supported")

1. Macro AUC of the engine exceeds that of the popularity baseline by **at least 0.05**, and the **95 % bootstrap interval of the paired difference lies above 0**.
2. **Calibration slope between 0.8 and 1.2** after out-of-fold calibration: calibrated probabilities come from a logistic (Platt) map fitted only on training folds of the grouped cross-validation (section 6); the slope is the coefficient of a logistic regression of the gold 1 label on the logit of the out-of-fold calibrated probability, pooled over all (disease, drug) pairs.
3. **Precision at 10** (macro over diseases) of the engine is **above that of the best baseline** (best of the three, by macro value).

If any criterion fails the engine is declared "not supported" and its ranking stays exploratory.

## 5. Secondary endpoints (descriptive; they do not change the decision)

- **Temporal AUC** against gold 2, using only literature published up to Y = 2018 (Europe PMC `PUB_YEAR` <= 2018, for the engine, the proxy and the popularity baseline alike). Diseases enter only with at least one gold 2 drug. A sensitivity analysis removes drugs that are also in gold 1. Because PrimeKG is a present-day snapshot, the knowledge-graph prior is excluded from the temporal analysis (only the no-prior arm is valid there).
- **Recall at 20** (macro over diseases) against gold 1 and gold 2.
- **Brier score** of the out-of-fold calibrated probabilities (pooled pairs), reported next to the Brier score of a constant prevalence predictor.

## 6. Grouped cross-validation

- Folds are by **disease**: 5 folds (fewer if fewer than 10 diseases), diseases assigned by a seeded shuffle (seed 20261015). A disease is never in training and test of the same fold.
- Pharmacological class = first four characters (ATC level 3) of the lexicographically first ATC code in `data/universe.json`; drugs without ATC are their own class. In each fold, training pairs exclude every drug whose class contains a gold 1 drug of any test-fold disease, so a class is never split between training and test of the same fold. Test rankings use all universe drugs for the held-out disease.
- Engine components are fixed before running; the cross-validation fits only learned parts (the calibration map and any learned combination weights).

## 7. Ablations (all reported with every endpoint, same thresholds, none selected post hoc)

| Arm | Knowledge-graph prior | Background correction for popularity |
|---|---|---|
| A (full) | yes | yes |
| B | no | yes |
| C | yes | no |
| D | no | no |

The success decision of section 4 is made on arm A only.

## 8. Reporting rules

- Results are reported whatever they are, including failure of any threshold, with the benchmark size, the exclusion log, and the number of diseases with empty gold 2.
- The thresholds, baselines, bootstrap settings, folds and ablations above must not be changed after seeing results. Any later deviation is listed in `protocol/deviations.md` with its reason and does not replace this pre-registration.
- Known limits fixed in advance: with this benchmark size the bootstrap interval is wide; gold 1 is incomplete and partly circular; gold 2 is a name match in intervention fields and may include comparators; PrimeKG-based gold may reflect post-2018 knowledge.
