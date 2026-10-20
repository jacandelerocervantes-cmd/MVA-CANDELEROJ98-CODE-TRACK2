# Addendum - engine v5, per-pair retrieval of drug-and-disease evidence, read by LLMs (post hoc)

Written on 2026-10-16, **before** `src/track2/engine_v5.py` exists or ran, and before any v5 query, item, annotation, label, feature, refit or score. Facts looked up before writing (all from existing files, DEV-side or label-free): `addendum_engine_v4*.md`, `protocol/relevance_guideline_v2.md`, `protocol/annotation_guideline.md`, `out/engine/audit_llm_coverage_dev.json/.md`, `out/engine/validation_v4_llm.json`, `exploratory_v4_llm_dev.json`, `validation_v4_dilution.json`, `validation_v4_ablation.json`. This file extends and does not modify any of them, nor `engine_v4*.py`, `benchmark_v4.json` or `engine_v4_frozen.json`.

## 0. Status and honesty clauses

- **Post hoc development extending the pre-registered v4 programme.** The v4 decision stands unchanged. Nothing here is a pre-registered primary result of the original programme; the endpoints of section 6 are fixed *now*, before any v5 data exist, and are confirmatory only in that narrow sense.
- **Motivation (observed in the v4 audits, DEV only).** With v1/v4 retrieval (300 records per disease selected by pooled gene symbols, no disease or drug term), only 3.8 % of the papers that mention both the drug and the disease were retrieved; 35 % of gold-positive DEV pairs had the drug named in any retrieved paper; a direct "drug AND disease" query gives about 3.9 sentences per pair and covers 53 % of gold pairs (about 38 % have zero direct hits: the genuine low-literature regime); language and drug-synonym loss are immaterial. Relevance guideline v1 gave a DEV kappa of 0.565 (failed the 0.6 gate); guideline v2 gave 0.851 in the 800-item pilot.
- **Hypothesis (written beforehand).** Evidence retrieved *per pair* (drug AND disease) and read by LLMs, counted as relevant-benefit and relevant-harm sentences, adds to a prior-only model (PO) at low literature volume (dilution f = 0.03 and 0.01 of the K records of each pair) where the v4 gene-pooled retrieval did not.
- **Failure is reported as is.** If the lower bound of the paired difference is not above 0 at both fractions, the hypothesis is reported as not supported (partial if only one) and the features are not adopted. If the gate of section 7 fails, the hypothesis is reported as not testable with these labels.
- **Sealed HELD-OUT.** The HELD-OUT of `benchmark_v4.json` has not been scored for any LLM feature, and **this phase scores nothing and creates no lock**. The HELD-OUT has been opened by v4, dilution and ablation scoring (non-LLM); the v1-guideline LLM disagreement analysis used HELD-OUT sentences (disclosed in `addendum_engine_v4_llm_guideline_v2.md`). For v5 the HELD-OUT items are new (different retrieval, different ids `ho5-`), to be annotated by fresh blind annotators.
- Seed everywhere: **20261015**. Every number is produced by `src/track2/engine_v5.py`. Limitations stated now: random thinning of records is a stress test of volume dependence, not a measurement on rare diseases; the two LLM annotators are of the same families, share biases and know the drugs; the drug AND disease query itself favours co-mentioned (often indicated) pairs.

## 1. Per-pair retrieval

- **Candidate set P_d: identical to v4** (`out/engine/llm_build/candidate_pairs.json`, SHA-256 recorded in the manifest): the 30 drugs with the highest frozen B0 score at f = 1.0 union the gold-1, gold-3 and (where it exists, 45 of 145 diseases) gold-2 drugs. Disclosed as in v4 (section S2 there): **gold sets define which pairs enter**, for DEV and HELD-OUT, and are never used in a query, an item or a feature.
- **One Europe PMC query per (disease, drug).** `(drug clause) AND (disease clause) AND HAS_ABSTRACT:y`, each clause an OR of `TITLE_ABS:"term"`; default relevance order; `pageSize` = K; `cursorMark` = `*`; **no year limit**; at most 4 requests per second; each response cached under `data/cache_engine/v5/`.
- **Drug terms.** The project (ChEMBL preferred) name first, then, in this priority: salt-stripped form, hyphen/space variants, ChEMBL molecule synonyms available locally (cached ChEMBL molecule records), ClinicalTrials.gov intervention names and other names (cached ClinicalTrials.gov API v2 studies); terms of at least 4 characters; within a source, terms without digits first, then shorter, then alphabetical. The **query uses at most 8 drug terms** (UNVALIDATED); sentence matching uses **all** terms (so a sentence is never lost because its synonym fell outside the first 8). Synonyms are built by the existing `engine_v4_llm_audit.build_synonyms`.
- **Disease terms.** The benchmark name plus the shortest other names of the benchmark `name_set` (MONDO, Orphanet and HPO names and synonyms already stored in `benchmark_v4.json`), at most 6 terms in the query, via the existing `engine_v4_llm_audit.query_disease_terms` (UNVALIDATED).
- **K = 15 records per pair** (UNVALIDATED; conservative default 10, the size of the audit's direct-query probe; 15 is used because cost is bounded by sentence caps and K is not tuned). K is never changed after seeing a result.
- Recorded per pair: the drug and disease terms, the query string, `hitCount`, the number and ids of the retrieved records in order, and the retrieval date. The "direct record count" of a pair is the number of retrieved records (at most K).

## 2. Sentence selection and verification

- Text of a record = `title + ". " + abstract` (as in v4). Sentences by `empirical_literature_control.split_sentences`. A sentence is selected if it names the drug by **any** of its terms (whole word, case-insensitive). **At most 3 sentences per record**, the first three in sentence order. No gene-symbol gate and no disease-term gate on the sentence (the query already required the disease in the record; the annotators judge relevance).
- **Every sentence is verified verbatim** against Europe PMC text re-fetched by external id, exactly as in v4: `track2.verify_literature.squash` rule through `engine_v4_llm.fetch_source_texts`, `verify_sentences` and the declared join-period repair `repair_join_period`. Sentences `not_found` or `not_checkable` are dropped before annotation and counted.
- Items are de-duplicated on (disease, drug, pmid, sentence), shuffled (seed 20261015, different from v4), given opaque sequential ids prefixed **`dev5-`** and **`ho5-`**, and cut into **shards of at most 100 items** in `out/engine/llm_input_v5/`. An item carries exactly `id, drug, disease, sentence, pmid`: no gold flag, no score, no candidate-membership flag (checked by test).

## 3. Annotation (done by the author after this phase)

Two blind LLM annotators (Opus, Sonnet), fresh sessions, each seeing only the items and the guidelines: **direction guideline v1 (`protocol/annotation_guideline.md`, unchanged)** and **relevance guideline v2 (`protocol/relevance_guideline_v2.md`)**. Both files are frozen by SHA-256 now:

| file | SHA-256 |
|---|---|
| `protocol/relevance_guideline_v2.md` | `5f55313ea0046b795f32e5fd28baa7330b249c8c9b16c8ce3e51d452808be5d9` |
| `protocol/annotation_guideline.md` | `f56abc1ff93ec8e2a06d8eabb062c318e51bba53cd1c0e738a6a6092fea2e45d` |

The manifest recomputes both; any difference invalidates the annotation. The third direction reader is the rule extractor `engine_extractors.rule_direction` (as in v4). HELD-OUT items may be annotated before the gate result, because labels do not depend on outcomes.

## 4. Features and models

- A verified sentence is **relevant** if both annotators say `relevance = yes`; its **direction** is the label shared by at least two of (Opus, Sonnet, rule); `unresolved` sentences are excluded (the v4 rule, unchanged).
- Per (disease, drug) and per sample: `n_benefit` = relevant sentences with direction benefit; `n_harm` = relevant with harm; `n_records` = retrieved direct records. Model inputs `log1p(n_benefit)`, `log1p(n_harm)`, `log1p(n_records)`. No gold, label or score enters the construction (no gold argument; tested).
- **PO** (comparator): the prior group of the frozen B0 feature matrix (columns 4, 6, 7: log popularity and the two PrimeKG prior columns), L2 logistic, refit. **PO+F** (primary): PO plus the three features. Secondary diagnostics: **PO+volume** (PO plus `log1p(n_records)` only, to separate literature volume from LLM-read content) and PO+F restricted to the two labelled features.
- **Training and tuning, DEV only.** Heads by `engine_v4.fit_head` (L2, standardised features); training rows use the f = 1.0 features of the training diseases with ATC level-3 class masking as in v1; the fitted head is applied to the thinned features at each fraction. **Nested grouped CV on DEV**: outer = the v4 grouped 5-fold (groups of identical gold sets never split), inner grouped 3-fold choosing C in {0.1, 1, 10} (UNVALIDATED) by mean macro AUC. The final DEV fit chooses C by the inner 3-fold over all DEV and is frozen with its SHA-256 before any HELD-OUT label of v5 is read.

## 5. Dilution

Thinning is applied to the per-pair records: for fraction f in {1.0, 0.3, 0.1, 0.03, 0.01} and replicates 1-5 (f = 1.0 has one deterministic replicate), **each retrieved record of a pair is kept independently with probability f** (Bernoulli, `numpy.random.default_rng([20261015, disease_rank, drug_index, fraction_index, replicate])`), so the expected number of kept records is 15, 4.5, 1.5, 0.45, 0.15 and a pair can lose all its records. (A rounded count would make f = 0.03 and f = 0.01 identical; Bernoulli keeps the fractions distinct. Stated, UNVALIDATED.) Annotation is done once per (disease, drug, pmid, sentence) and reused by every sample. The B0-prior columns are not thinned (they do not use v5 records).
**Natural-scarcity stratum:** pairs whose number of retrieved direct records at f = 1.0 is 0, 1 or 2 (UNVALIDATED threshold), reported separately and never pooled into the primary endpoint; the stratum metric is the macro AUC over the diseases that have at least one positive and one negative pair in the stratum (count reported), at every f.

## 6. Endpoints (fixed beforehand)

- Evaluation universe = P_d (as in v4 LLM): all models and the popularity baseline are scored on the same pairs; only paired differences are interpreted.
- **Primary (HELD-OUT, confirmatory):** the paired macro AUC difference **PO+F minus PO at f = 0.03 and at f = 0.01**, 95 % group bootstrap (2000 resamples, percentile, seed 20261015, resampling the v4 disease groups, same indices for both models; the plain disease-level bootstrap as a sensitivity). **Decision: the hypothesis is supported only if the lower bound is above 0 at both fractions**; one fraction: partial, not supported; none: not supported.
- **Secondary (descriptive, no multiplicity correction):** the same difference at f = 1.0, 0.3, 0.1; the three criteria of the v1 protocol (`addendum_engine.md` section 4, via `engine_score.decide`) for the full model, with the caveat that they are evaluated on P_d; precision at 10 and recall at 20 within P_d; PO+volume and the labelled-only variant; the natural-scarcity stratum; DEV out-of-fold values of all of these; coverage (share of pairs with at least one counted relevant sentence) by split and by gold membership on DEV.
- **Multiplicity:** the primary endpoint alone is confirmatory; everything else is exploratory and labelled so.

## 7. Quality gate (unchangeable)

On **all DEV items** the Cohen's kappa of the relevance labels (Opus vs Sonnet, paper bootstrap 95 % interval, 2000 resamples, seed 20261015) must have a **point estimate >= 0.6**. Otherwise the relevance labels are unreliable, the features are not used, the HELD-OUT is not scored, and the hypothesis is **not testable**. The direction kappa (and interval) is reported, not a gate. The same agreement on HELD-OUT items is reported after scoring and never used for a choice. The gate is not re-run with another guideline, subset or annotator after seeing it.

## 8. Single HELD-OUT use and lock

1. `engine_v5 retrieve` and `engine_v5 build` (this phase) are label-free and call no LLM.
2. The author dispatches the annotators; labels go in `out/engine/llm_labels_v5/`.
3. A later stage computes agreement, the gate, the nested CV on DEV and freezes the PO and PO+F heads (coefficients, C, SHA-256 of the model, of the candidate file, of the retrieval file, of this protocol).
4. The HELD-OUT is scored **once**, after the freeze and after a lock file is created, which refuses a second run. If any frozen hash differs, the run is refused.

## 9. Leakage audit items

(L-a) Gold sets decide which pairs are retrieved and annotated (disclosed), never the query, item or feature. (L-b) Gold pairs are always in P_d, negatives are the 30 hardest B0 negatives: evaluation inside P_d with the same set for both models; residual selection effect stated. (L-c) The direct record count is itself informative of indication (co-mention); hence the PO+volume diagnostic. (L-d) No year limit: the records can postdate approval or label changes, and openFDA gold 3 derives from labels informed by that literature; stated as inherent. (L-e) Drug synonyms come from ClinicalTrials.gov interventions cached for the gold-2 build; they are keyed by drug, used only to name the drug, and no condition field is read. (L-f) Annotator parametric knowledge: DEV audit of the share of `yes` and `unclear` among gold-1 versus other pairs, reported, not corrected. (L-g) Groups never straddle DEV and HELD-OUT (inherited split, asserted); papers shared between the splits are counted and reported, not removed. (L-h) Europe PMC changes over time: the retrieval date and the cached responses are stored; items and shards are hashed in the manifest. (L-i) HELD-OUT labels do not exist for v5; the HELD-OUT gold flags are read only through the v4 candidate file (disclosed above) and, for scoring, after the lock.

## 10. UNVALIDATED values (each stored with `status = UNVALIDATED` and a `conservative_default`; a test fails if any lacks them)

K = 15 (default 10); drug terms in the query = 8; disease terms in the query = 6; sentences per record = 3; minimum term length 4; Bernoulli thinning rule; natural-scarcity threshold of 2 records; top-30 candidate size (inherited); relevance conjunction and the 2-of-3 direction vote (inherited); C grid {0.1, 1, 10} and the inner 3-fold; kappa gate 0.6; evaluation on P_d; 4 requests per second.

## 11. Run log and deviations
(appended after the label-free stages)
