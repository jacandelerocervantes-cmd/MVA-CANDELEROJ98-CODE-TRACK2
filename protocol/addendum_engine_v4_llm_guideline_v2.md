# Addendum - engine v4 LLM experiment, relevance guideline v2 (post hoc)

Written on 2026-10-16. This file extends `addendum_engine_v4_llm.md` and modifies nothing in it.

## 1. What this is

`protocol/relevance_guideline_v2.md` is a **post hoc revision** of `protocol/relevance_guideline.md`. It was written **after** observing annotator disagreement on the **HELD-OUT** sentences (Opus vs Sonnet relevance kappa 0.591 on 3396 items; Sonnet `yes` 642, Opus `yes` 239; Opus `unclear` 573, Sonnet 157). The disagreement analysis is `out/engine/llm_labels/disagreement_analysis_heldout.json` (aggregate counts and held-out sentence features only).

## 2. What the revision used and did not use

- Used: the two annotators' labels and the sentence text of HELD-OUT items, and nothing else. No engine score, gold set, validation result or benchmark membership was read for the revision. Examples in v2 are synthetic.
- **No HELD-OUT engine score exists yet.** The `heldout` stage of `engine_v4_llm` has not run and its lock file does not exist.
- The revision is, however, informed by HELD-OUT sentences. The HELD-OUT labels therefore cannot be treated as a clean test of v2: **the HELD-OUT would be re-annotated under v2 by new blind annotators** (fresh sessions that have not seen v1 labels or this analysis), and the v1 HELD-OUT labels are not used for any feature, fit or score under v2.

## 3. Gate and order of use

- The gate is unchanged in form: DEV relevance kappa (Opus vs Sonnet-type blind annotators) **point estimate >= 0.6**.
- Under v2 it is evaluated first on the **DEV pilot shards** (a small fixed subset chosen by shard number before annotation, not by looking at labels), then on **all DEV**. If the point estimate is below 0.6 on all DEV, the relevance labels are reported unreliable, the LLM features are not used, and no HELD-OUT scoring is run (as in section 6 of the base addendum). Passing the pilot alone does not release the HELD-OUT.
- DEV has not been annotated under v2 at the time of writing; nothing here depends on its outcome.

## 4. Status of v1

**v1 remains the pre-registered guideline.** Results under v1 (the DEV relevance kappa and the gate) stay reported as the pre-registered path. Any use of v2 is labelled a post hoc development branch, reported separately, not merged into the pre-registered decision, and cannot turn a failed v1 gate into a pass; a v2 outcome supports only a post hoc statement ("under a guideline revised after seeing annotator disagreement, ...").

## 5. Limitations stated now

- v2 was written by looking at HELD-OUT disagreement, so any v2 gain in agreement on HELD-OUT is optimistically biased; DEV agreement under v2 is the clean estimate.
- New annotators are LLMs of the same families; they share biases and parametric drug knowledge.
- The new rules (notably C3, the sentence must supply the condition link) change what counts as `yes`; the downstream features count fewer `yes` sentences, which changes coverage. This is reported, not tuned.
- Every threshold value in v2 is UNVALIDATED (status and conservative default: the v1 wording).

## Pilot decision criterion (written 2026-10-16, before any v2 annotation)

Facts known when this was written: with guideline v1 the DEV relevance kappa over all 7637 items is 0.565 (95 % CI 0.546 to 0.584), below the 0.6 gate, so the v1 LLM-feature experiment is closed as "not testable with these labels"; the 8-shard pilot under v1 gave 0.625, which over-estimated the full-DEV value.

Pilot: DEV shards 002, 004, 006, 025, 026, 028, 057, 062 (seed 20261015, the same as the v1 pilot), annotated blind by Opus and Sonnet under `protocol/relevance_guideline_v2.md`. The pilot advances to full DEV re-annotation only if BOTH: (a) the relevance kappa point estimate is at least 0.75 and (b) the lower end of its paper-bootstrap 95 % interval is at least 0.70. Otherwise the v2 route is closed and the result is reported as is. The 0.6 gate on full DEV, the DEV-only rule, and the single opening of the sealed HELD-OUT are unchanged. No threshold is lowered after seeing results. HELD-OUT items are re-annotated under v2 only after full DEV passes the 0.6 gate.
