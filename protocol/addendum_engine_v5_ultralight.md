# Addendum: engine v5 ULTRALIGHT (written before any v5 label exists)

Fixed on 2026-10-17, before the select stage ran. Seed everywhere: 20261015. Every number below is produced by `src/track2/engine_v5_ultralight.py`.

## 1. Status
- Post hoc development extending the pre-registered v4 programme. The v4 outcome stands unchanged whatever this addendum yields; nothing here is a pre-registered primary result of the original programme.
- The full design `protocol/addendum_engine_v5.md` (SHA-256 `1e80678f4f7f54d821496e5db2f0e9897ad308e448ec3820a572fa4ad446e35b`) is left untouched. It is a documented but UNEXECUTED alternative.
- Executed variant: ULTRALIGHT (this file). Documented, unexecuted alternatives, to show that the procedure scales with budget (tokens use the manifest method of `out/engine/llm_input_v5/manifest.json`, both annotators, one pass over both tasks):
  - Full v5 (DEV + HELD-OUT, 47,461 items, 476 shards): 14,741,406 tokens.
  - Reduced (my reading: HELD-OUT only, all 7,027 items, both annotators, no DEV): 2,191,890 tokens.
  - Ultralight (this file): see section 8; about 0.26 M tokens at 155 tokens per item.

## 2. Data
- Only HELD-OUT pairs (42 diseases, 1,350 candidate pairs of the v4 candidate file). No DEV annotation, no fitting, no training.
- Records come from the frozen v5 retrieval (K = 15 records per pair). The v5 thinning is reproduced exactly: `engine_v5.thin_records`, each retrieved record kept independently with probability f, `numpy.random.default_rng([20261015, disease_rank, drug_index, fraction_index, replicate])`, fraction_index 3 (f = 0.03) and 4 (f = 0.01), replicates 0-4 (the replicate numbering of `engine_v4_dilution.sample_keys`).
- Selected items: the verified sentences (item ids exactly as in `out/engine/llm_input_v5/heldout_shard_*.json`) of records that survive the thinning at f = 0.03 or f = 0.01 in at least one of the 5 replicates. Count (computed): 1,281 of 7,027 HELD-OUT items (18.2 %, 332 pairs have at least one item).
- Items carry only id, drug, disease, sentence, pmid. No gold flag, no score.

## 3. Annotation
- Sonnet labels all selected items (shards of at most 100) with relevance guideline v2 (`protocol/relevance_guideline_v2.md`, SHA-256 `5f55313ea0046b795f32e5fd28baa7330b249c8c9b16c8ce3e51d452808be5d9`) and direction guideline v1 (`protocol/annotation_guideline.md`, SHA-256 `f56abc1ff93ec8e2a06d8eabb062c318e51bba53cd1c0e738a6a6092fea2e45d`).
- Opus labels a seeded random 30 % of the Sonnet shards (`random.Random(20261015 + 7)`, round(0.3 * 13) = 4 of 13 shards) for the quality gate. Labels placed at `out/engine/llm_labels_v5_ul/sonnet_NNN.json`, `opus_NNN.json` = `{"annotator", "labels": [{"id", "direction", "relevance"}]}`.

## 4. Feature (fixed, a-priori, no training)
- For each pair, B = 1 if at least one surviving verified sentence has relevance = yes AND direction = benefit by Sonnet, and Opus does not label the same item relevance != yes (if Opus did not label the item, Sonnet alone decides). Surviving = in a record kept by the thinning of that fraction and replicate.
- score = PO_score + delta * B, where PO_score is the frozen PRIOR_ONLY score (`engine_v4_ablation.frozen_scorer`, frozen B0 weights, literature columns set to the DEV mean) and delta = 0.5 * SD (ddof 1) of the PO scores among the candidate pairs of the same disease. B = 0 for non-candidate drugs.
- delta = 0.5 SD is UNVALIDATED, chosen now without any label; conservative default delta = 0 (no effect, score = PRIOR_ONLY). It is never changed after a result.

## 5. Primary endpoints (fixed)
- Paired macro AUC difference (score minus PRIOR_ONLY) on HELD-OUT at f = 0.03 AND f = 0.01, mean over the 5 replicates per disease, groups of identical gold sets as units, group bootstrap 95 % CI (`engine_v4_dilution.boot_idx`, `engine_score.diff_with_ci`).
- Supported only if the lower bound is above 0 at both fractions. Otherwise: not supported.
- Coverage = share of candidate pairs with B = 1, reported at each fraction (mean over replicates). If coverage is below 1 % of pairs at a fraction, the result is "not informative" whatever the interval.

## 6. Quality gate (written now)
Opus-vs-Sonnet relevance kappa (`engine_extractors.kappa`), point estimate on the overlap items, must be >= 0.6. Otherwise the labels are unreliable, B is not used, no scoring is done and the hypothesis is "not testable".

## 7. Governance
- Single scoring, protected by a lock file (`out/engine/heldout_v5_ul_scored.lock`); HELD-OUT gold is read only after the lock exists. Failure is reported as is. The result cannot change the pre-registered v4 outcome.
- Limitations: random thinning is a stress test of volume dependence, not rare diseases; annotators share biases and know the drugs; the selection keeps only pairs with surviving records, so coverage is small by design.

## 8. Cost (scalability claim)
Tokens = 155 per item per annotator. Computed by the select stage (`out/engine/llm_input_v5_ul/manifest.json`): items per HELD-OUT pair, Opus items and the total token estimate. Cost grows linearly with the number of diseases (items per pair is a constant of the retrieval and the thinning, not of the disease count).
