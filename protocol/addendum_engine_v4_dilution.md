# Addendum - engine v4, literature-dilution experiment (low-report regime)

Written on 2026-10-15, **before any diluted feature, embedding or score exists**. At the time of writing only the following non-result facts were looked up: the number of retrieved papers per benchmark disease already stored in `data/cache_engine/v4_features.pkl` (128 of 145 diseases sit at the retrieval cap of 300; the smallest has 11), and that the Fanconi anaemia and MVA gene queries are not in the retrieval cache. This file extends `addendum_engine_v4.md` and does not modify it, nor `engine_v4.py`, `benchmark_v4.json`, `validation_v4.json`, `engine_v4_frozen.json` or any earlier file.

## 0. Status and honesty clauses

- **Post hoc development experiment, not a pre-registered primary result.** The v4 decision (sealed HELD-OUT, three criteria) stands unchanged.
- **Why.** The low-literature stratum of v4 (L, fewer than 50 Europe PMC hits) has 1 disease, so performance for diseases with few publications (the regime of the ultra-rare disease at the centre of this project) cannot be tested directly: diseases enter the benchmark because they have known drugs, which implies much literature. The experiment **simulates** low-report diseases by thinning the retrieved papers of each benchmark disease. A simulation of this kind removes reports at random; real rare-disease literature also differs in kind (other gene-drug co-mention patterns, fewer trials), which random thinning does not reproduce. The result is therefore a stress test of volume dependence, not a measurement on rare diseases.
- **Any failure is reported as is** (for example a gap whose lower bound never exceeds 0, a non-monotone curve, or a collapse of precision at 10).
- Every number is produced by `src/track2/engine_v4_dilution.py`. Seed everywhere: **20261015**.
- The HELD-OUT set is used **once, only for the frozen B0 and E, after all design in this file is fixed** (section 7).

## 1. Design

**Unit and data.** Every disease of `benchmark_v4.json`: 103 DEV and 42 HELD-OUT diseases. For each disease the papers are the cached Europe PMC retrieval already used by v4 (`engine_score.fetch_papers`, cache only, at most 300 papers, query `engine_score.build_gene_query(symbols, 2024)`). No new retrieval is made for benchmark diseases and nothing is downloaded.

**Fractions** f in {1.0, 0.3, 0.1, 0.03, 0.01}. For a disease with n retrieved papers the subsample size is

`n_sub = min(n, max(5, ceil(f * n)))`

so a **floor of 5 papers** applies at every fraction below 1 (it binds when f * n < 5, for example f = 0.01 and f = 0.03 for every disease with n <= 300, and f = 0.1 when n < 50), and a disease with fewer than 5 papers keeps all of them (floor "where applicable": never more than n). At f = 1.0 the subsample is the full retrieval, so it is the v4 setting and has a single deterministic replicate (the 5 replicates would be identical; this is stated, not hidden). **Number of papers** on the volume axis is `n_sub` (papers retrieved after thinning). Because the engine caps retrieval at 300, `n_sub <= 300`; the Europe PMC hit count of each disease is also recorded.

**Replicates.** 5 independent subsample replicates per fraction below 1. Papers are drawn without replacement with `numpy.random.default_rng([20261015, rank, frac_index, rep])`, where `rank` is the position of the disease in the sorted list of all benchmark disease names, `frac_index` the index of the fraction in the list above and `rep` in 0..4. The selected papers keep their original retrieval order (the v4 drug-sentence rule depends on paper order).

**Recomputed (paper-derived) features, per diluted disease.** With the unchanged v4/v1 code (`engine_score.literature_counts`, `engine_score.disease_matrix`, `engine_v4.drug_sentences`, `engine_v4.disease_text`): the four literature columns (log1p of paper, benefit, human and harm counts), `log_popularity` (unchanged: global drug popularity), `excess_over_popularity` (recomputed because it uses the number of papers with a gene symbol), and the two PrimeKG prior columns (not paper-derived: unchanged, copied from the v4 feature cache). For E: the three embedding features (`emb_cos_biomedbert`, `emb_cos_sapbert`, `emb_has_text`) recomputed from the drug sentences of the diluted papers, embedded with the same frozen BiomedBERT and SapBERT (local `data/models/`, isolated `python -I` worker, 128 tokens, same pooling). Texts already embedded by v4 are reused from `texts.json`; only new texts are embedded.

**Models and baselines (FROZEN, no refitting).** Weights read from `out/engine/engine_v4_frozen.json` after verifying its model hash (`0543eb2c...afedb9`), its file SHA-256 (`09e3a8e8f5a32ea79d455d69439a2edd86a96058fc85cd5a62be1c13b398e673`) and the DEV and HELD-OUT id hashes (`3623a3a6...`, `2719907d...`). **B0 = primary model**, **E = secondary**; M, conformal and every other model are out of scope. The popularity baseline is the v4 one, unchanged by dilution (it is global drug popularity and does not depend on the disease's papers); the random baseline is the analytic one of v1. The proxy baseline is not used here.

**Important caveat on DEV.** The frozen weights were fitted on the DEV diseases, so DEV scores here are **in-sample for the weights** (the diluted features are new, the weights are not). DEV is reported for completeness and for a larger n, but it is optimistic; the HELD-OUT set (never used for fitting) is the only out-of-sample reading, and the primary break-even volume is the HELD-OUT one.

## 2. Endpoints (per fraction, separately for DEV and HELD-OUT)

For each disease and fraction the per-disease metric is the mean over the replicates of the per-replicate metric. Then:
1. **Macro AUC** of B0 and of E, and of popularity.
2. **AUC gap vs popularity** = macro AUC(model) minus macro AUC(popularity), with a **95 % disease-level bootstrap interval**: 2000 resamples, percentile, seed 20261015, resampling the disease **groups** of v4 section 4.1 (a group = one observation = the mean over its member diseases; identical to a disease-level bootstrap when every group has one disease). The plain disease-level bootstrap (every disease its own unit) is also stored as a sensitivity.
3. **Precision at 10** and **recall at 20** (macro, same bootstrap), for B0, E, popularity and random.
4. Mean and median `n_sub` over diseases (and replicates), and the SD across the 5 replicates of the macro gap (replicate noise).

## 3. Gap versus log10(papers) and break-even literature volume

- Observations: one per (disease, fraction, replicate), i.e. gap_i = AUC(model, diluted features) - AUC(popularity) of that disease, and x_i = log10(n_sub_i).
- Fit: ordinary least squares `gap = a + b * log10(n_sub)`.
- Interval: **cluster bootstrap** (2000 resamples, seed 20261015, resampling disease groups with all their observations), refit each time, percentile 95 % band of the fitted line on a grid of 200 values of n_sub equally spaced in log10 between 5 and 300; the slope b is also reported with its percentile interval.
- **Break-even literature volume** = the smallest grid value N such that the lower bound of the fitted gap is greater than 0 for every grid value >= N (so the model beats popularity with the interval excluding 0 from N papers upward). If the lower bound exceeds 0 over the whole grid, it is reported as "at or below 5 papers (smallest simulated volume)"; if it never does at 300, as "not reached up to 300 papers". The point estimate crossing (fitted gap = 0) is reported next to it. Computed for B0 (primary) and E (secondary), for DEV and for HELD-OUT. The decision-relevant value is the **HELD-OUT B0** one; DEV is a sensitivity (in-sample weights).
- Limitation stated now: disease identity and natural n are partly confounded in the pooled fit (128 of 145 diseases have n = 300); the dilution gives each disease its own range of n_sub, but the line is a pooled summary, not a per-disease law.

## 4. Mapping the project's diseases on the volume axis

Reported with the same retrieval (`engine_score.fetch_papers`, at most 300 papers, query for the gene set, year <= 2024) and hit count:
- **Fanconi anaemia**: the gene set of the project's Fanconi case (`out_case_fanconi_anemia_run1/submission_genes.json`: ATR, FANCA).
- **MVA gene set**: BUB1B, CEP57, TRIP13.
- **The 6 benchmark diseases with the fewest retrieved papers** (from the v4 feature cache).

These queries are not in the cache; if absent they are retrieved from Europe PMC with cursorMark pagination at no more than 4 requests per second (the existing limiter), a handful of requests in total. For each: number of papers retrieved, hit count, and its position relative to the break-even volume (above, below, or at the retrieval cap). Where the retrieval is at the cap of 300, the statement is that the engine's gene-level retrieval is saturated and the disease sits in the abundant-literature regime by this axis; the report does not extrapolate beyond 300.

## 5. Directions stated before running

(D1) The gap and precision at 10 of B0 decrease as f decreases (fewer papers, fewer co-mentions); (D2) at the smallest volumes (5 papers) B0 approaches the popularity-only performance because only `log_popularity` and the two prior columns carry information; (D3) E adds nothing over B0 at low volume (few sentences to embed, `emb_has_text` mostly 0). A result opposite to D1-D3 is reported as such. No threshold on the gap is used to accept or reject the engine here.

## 6. UNVALIDATED values (each with status and conservative default; none promoted by this run)

Fractions {1.0, 0.3, 0.1, 0.03, 0.01}; floor of 5 papers; 5 replicates; the log-linear form of the gap curve; the 200-point grid and the "lower bound above 0 for all larger volumes" definition of break-even; the use of in-sample weights on DEV. A test fails if any lacks a status or default.

## 7. Single HELD-OUT use and lock

1. Feature building and embedding of diluted texts are label-free and run for all 145 diseases.
2. The DEV command scores the frozen models on diluted DEV diseases and writes the DEV part of `out/engine/validation_v4_dilution.json`. After it, no design element of this file changes.
3. The HELD-OUT command verifies the hashes of section 1, **creates `out/engine/heldout_v4_dilution_scored.lock` before reading any HELD-OUT label**, scores only the frozen B0 and E, once, and refuses to run if the lock exists.
4. Output `out/engine/validation_v4_dilution.json` is written with sorted keys.

## 8. Run log and deviations
(appended after the runs; any deviation from this file is listed here with the reason)

- Features: 145 diseases x 21 samples (1 at f = 1.0 + 4 x 5), cache only. f = 1.0 reproduced the v4 feature matrix and drug texts exactly (maximum absolute difference 0.0 in DEV and HELD-OUT). 5785 new drug texts were embedded with the frozen models (local, no download).
- DEV run, then HELD-OUT run once (lock `out/engine/heldout_v4_dilution_scored.lock`, frozen file SHA-256 verified). The Europe PMC mapping used 8 requests (the two project gene sets were not cached).
- Deviations: none from sections 1-7. Clarifications: (i) the Fanconi gene set is ATR + FANCA, the genes of the project's Fanconi case files; (ii) the abundance axis is capped at 300 retrieved papers by the engine; (iii) the PrimeKG prior columns and `log_popularity` of B0 are not paper-derived and stay fixed under dilution, so the model keeps non-literature signal at 5 papers (protocol direction D2 anticipated approach to popularity-only performance; the result is reported as it came out in `out/engine/validation_v4_dilution.json`).
