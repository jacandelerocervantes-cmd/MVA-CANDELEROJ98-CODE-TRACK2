# Addendum - evidence engine v4 (phase 1): scale-up, literature-volume strata, embeddings, abstention

Written on 2026-10-15, **before any v4 benchmark, feature, embedding, fit or score exists**. At the time of writing, only eligibility counts (no literature, no score) were computed, to learn whether the requested scale is reachable; they are quoted in section 1.2. This file extends, and does not replace, `addendum_engine.md` (v1, pre-registered primary result: supported on 17 diseases), `addendum_engine_v2.md` and `addendum_engine_v3.md`. Nothing here modifies `frozen_protocol.md`, `engine_score.py`, `engine_score_v2.py`, `engine_v3.py`, `benchmark_v2.json`, `validation*.json`, `case_wilms.*` or any earlier protocol file. Phase 2 (LLM annotation uses) is out of scope.

## 0. Status and honesty clauses

- **Post hoc development extension, not a pre-registered primary result.** The v1 decision (supported on 17 diseases) stands unchanged. v4 asks whether the engine, with its v1 features and refit weights, still beats literature popularity when the benchmark is scaled up and stratified by literature volume (the setting of low-report diseases such as MVA), and whether three declared additions help.
- **Failure to improve is reported as is.** If an addition does not pass its inclusion rule, the final model is the simpler one and the report says so. If the sealed HELD-OUT does not meet a criterion, it is reported as not met.
- Everything seen when v4 was designed: v1, v2 and v3 results (including the opened v3 HELD-OUT, whose 14 diseases are now part of DEV), and the Wilms case. Only the new v4 HELD-OUT diseases (never seen in any benchmark) can support a claim.
- Every number is produced by `src/track2/engine_v4.py`. Seed everywhere: **20261015** (derived seeds are +1, +2 as in v3).

## 1. Benchmark v4 (scale-up) and split

### 1.1 Rules (same code path as v2: identifier-and-synonym matching, unchanged)
Pool = every PrimeKG disease node that has at least **one** approved universe drug (`data/universe.json`, 2705 drugs) linked by an `indication` edge (gold 1). Matching, name sets, umbrella cap (more than 30 matched disorders excluded), the "rare genetic" rule (Orphanet-catalogued identity disorder with Mendelian / mitochondrial / sex-linked inheritance; curated candidates keep the looser alternative), the exclusion of MVA / premature chromatid separation / Fanconi anaemia, the exclusion of gold 1 containing an Enzyme / Gene / Cell product, and gold 1 as the endpoint gold are **those of `addendum_engine_v2.md`, imported from `engine_benchmark_v2.py` without modification**. The 45 diseases of `benchmark_v2.json` are all retained (they become DEV).

### 1.2 Three declared changes, and why
Eligibility counts computed before writing this file (diseases not in benchmark_v2; pool of 1288 PrimeKG diseases with at least 1 universe drug, 820 with at least 2):

| minimum genes / minimum gold 1 drugs | eligible new diseases | with distinct gold 1 sets |
|---|---|---|
| 2 / 2 (v1, v2 rules) | 8 | 7 |
| 1 / 2 | 26 | 19 |
| 2 / 1 | 50 | 37 |
| 1 / 1 | 100 | 61 |

Under the unchanged rules only 7 further distinct diseases exist, so neither the target of at least 150 diseases nor a sealed HELD-OUT of at least 40 new diseases is reachable. To reach the largest possible benchmark from the same sources, **three changes are declared now, all applied together and never tuned**:
- **D1** minimum number of genes 2 -> 1 (the engine query uses gene symbols; one gene is enough to retrieve literature, and most low-report diseases are single-gene);
- **D2** minimum number of gold 1 drugs 2 -> 1 (a disease with one positive still has a defined AUC; it enters the macro mean as in v1 if it has one positive and one negative);
- **D3** diseases with an identical gold 1 set are **not dropped** but form one **disease group**. The group is the unit of splitting, of cross-validation folds and of the bootstrap (a group mean of the member diseases' metrics is one observation); identical gold sets are the dependence that the v1 collapse rule tried to remove, and grouping removes it without discarding diseases.

No other rule is relaxed (umbrella cap, rare-genetic rule, biologic exclusion, application-case exclusion stay). Expected size, from the table: 45 + 100 = **145 diseases, below the target of 150**; the shortfall of at least 5 is reported and no further rule is relaxed. The code writes the retained and dropped counts by reason (selection log of every evaluated PrimeKG disease with its first failing rule) to `out/engine/benchmark_v4.json`.

### 1.3 Split (seed 20261015)
1. Groups = identical gold 1 set (sorted tuple of drug names). A new disease whose group contains a benchmark_v2 disease goes to DEV with it.
2. The remaining ("free") groups are sorted by their key, shuffled with `random.Random(20261015).shuffle`, and added to DEV in that order while the number of new DEV diseases is below `floor(n_free_diseases / 2)` (the last group may overshoot); all other free groups form the **sealed HELD-OUT**.
3. DEV = benchmark_v2 diseases + group-mates + the DEV share of free groups. HELD-OUT = eligible diseases not in benchmark_v2 and not in a v2 group. Required: HELD-OUT of at least 40 diseases; if fewer, the shortfall is reported and no rule is changed.
4. SHA-256 of both id lists (`hash_ids`: sorted names joined with newline) are computed by the benchmark command and appended in the freeze record (section 9) **before** any feature or score exists. HELD-OUT diseases are never used for design, tuning, selection, Platt or conformal calibration; the DEV command does not read their gold labels (it reads only their features if needed to cache literature; no score).

## 2. Literature-volume strata

- Stratifier = **Europe PMC hit count of the primary retrieval query** of the disease (`engine_score.build_gene_query(symbols, 2024)`, the same query the engine uses; one extra `pageSize=1` request, cached). Clarification: the engine keeps at most 300 abstracts per disease (`MAX_PAPERS`), so the number of *retrieved* papers cannot exceed 300 and a ">500" stratum would be empty; the hit count of the same query is the literature volume that the retrieval faces. The number actually retrieved is also reported per disease.
- Strata: **L < 50**, **M 50-500**, **H > 500** hits (cut points UNVALIDATED, conservative default: same values).
- Every result (AUC, difference to popularity, precision at 10, abstention, coverage) is stratified, for DEV (out-of-fold) and HELD-OUT, with group-level bootstrap intervals and the number of diseases / groups; strata with fewer than 5 groups are labelled "descriptive only, n too small".
- **Direction predicted before scoring**: (P1) the engine advantage over popularity (paired AUC difference) is **smaller in L than in H**, with the L interval likely to include 0; (P2) popularity AUC is not higher in L than in H for lack of rare-disease signal, i.e. the engine-minus-popularity gap is ordered L < M <= H; (P3) the embeddings feature adds least in L (few sentences to embed); (P4) abstention is higher in L than in H. A result opposite to P1-P4 is reported as such.
- Per-stratum results never change the overall decision; the three criteria are taken overall.

## 3. Candidate additions (closed list; nothing else is tried)

Base model **B0** = the 8 engine-v1 features (`engine_score.FEATURE_NAMES`, arm A), L2 logistic regression C = 1 on standardised features, ATC level-3 class masking of training rows as in v1, refit on the DEV diseases. B0 is the simpler model.

**E - embeddings feature (small logistic head).** Models (downloaded once, git-ignored `data/models/`, revisions recorded in section 9): `microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext` (mean pooling of the last hidden layer over non-padding tokens) and `cambridgeltl/SapBERT-from-PubMedBERT-fulltext` ([CLS] vector of the last layer, as recommended by its authors). Frozen, inference only, maximum 128 tokens, CPU. Disease text = cleaned disease name + ". Genes: " + sorted gene symbols. Drug text of a (disease, drug) pair = the first 3 sentences (retrieval order of papers, then sentence order) of the disease's retrieved papers that contain a disease gene symbol and mention the drug (the v1 sentence-level drug match). Features: `emb_cos_biomedbert`, `emb_cos_sapbert` (cosine of disease text and drug text), `emb_has_text` (1 if the pair has at least one sentence, else 0 and both cosines 0). **E = B0 features + these 3, same logistic head (C = 1).** Texts come only from retrieved abstracts: no PrimeKG, openFDA, gold or label text enters. (Sentence count 3 and 128 tokens: UNVALIDATED, conservative default same.)

**M - MLP head (ablation only).** `sklearn.neural_network.MLPClassifier`, one hidden layer of 16 ReLU units, L2 alpha 0.01, adam, 100 iterations maximum, seed 20261015, on the E feature set (standardised), same masking. M is **never eligible** for the final model; its DEV and HELD-OUT results are reported next to E and B0 as an ablation of head capacity. (All MLP hyper-parameters: UNVALIDATED, not tuned.)

**K - split-conformal abstention, alpha = 0.1 (UNVALIDATED, conservative default 0.1).** Class-conditional (Mondrian, by gold label) split conformal on the Platt-calibrated probability p of the final model: nonconformity a(x,1) = 1 - p, a(x,0) = p; for each label y the threshold q_y is the `ceil((n_y + 1)(1 - alpha))`-th smallest nonconformity among calibration pairs of label y (infinite if that rank exceeds n_y); the set of a pair is {y : a(x,y) <= q_y}; **abstain when the set does not have exactly one label** (size 0 or 2; both rates reported). Calibration pairs = pooled DEV out-of-fold predictions of the final model (the held-out never calibrates). Reported per stratum and overall, DEV and HELD-OUT: abstention rate (pairs), share of empty and of two-label sets, marginal coverage, positive-class coverage, negative-class coverage, abstention among the model's top 10 pairs of each disease. Exchangeability across diseases is not guaranteed (distribution shift between strata); coverage is measured, not assumed.

## 4. Endpoints and inclusion rules

1. **Endpoints unchanged**: the three criteria of `addendum_engine.md` section 4 evaluated by `engine_score.decide` (macro AUC minus popularity AUC >= 0.05 with the paired 95 % interval lower bound > 0; Platt calibration slope in [0.8, 1.2]; macro precision at 10 above the best of popularity, proxy, random). Bootstrap: 2000 resamples, percentile 95 %, seed 20261015, **resampling disease groups** (group = one observation = the mean over its member diseases; equal to a disease-level bootstrap when every group has one disease). Secondary: recall at 20, Brier score. A disease-level bootstrap is also reported for the HELD-OUT as a sensitivity.
2. **Inclusion rule (1-SE rule).** A component (E; K as a validity add-on) enters the final model only if, on the DEV out-of-fold results (same folds), the mean paired group-level AUC difference (candidate minus the simpler model) is **greater than 1 SE of the paired differences** (SE = sd(ddof 1) / sqrt(number of groups)) and greater than 0. K has no AUC; K is included if its DEV positive-class coverage is at least 0.90 minus the bootstrap SE of that coverage (validity), otherwise it is reported as failing. M is never included.
3. Complementary information (not a rule): noise threshold `2 * SE_diff` and minimum detectable difference `2.8 * SE_diff` from the DEV paired differences, as in the holdout-evaluation protocol.
4. **Single HELD-OUT scoring.** After freezing (section 9) the HELD-OUT is scored exactly once, guarded by the lock file `out/engine/heldout_v4_scored.lock`. In that one call all frozen models are scored: the final model, B0, E, M (ablation), the frozen v1 weights of the v3 freeze (`engine_v1_weights_frozen.json`, 31 diseases) as a reference, and the baselines popularity, proxy, random. The selection of the final model is made only from DEV.

## 5. Cost per disease (scalability)

For every disease the code records seconds and Europe PMC requests for retrieval (cold cache; for cached v1/v2 diseases the cold cost is estimated as 1 hit-count request + ceil(papers / 100) page requests and labelled "estimated"), seconds for feature building, for embedding texts and for scoring, and the total; the report gives mean, median, maximum per disease and the linear extrapolation to 1000 diseases. Zero manual edits: the pipeline from disease name to ranking is one command; the test file asserts that no disease-specific constant exists in the code.

## 6. Evidence-gap block with hierarchical borrowing (computed, no dosing)

For each benchmark disease: (i) **case-count tier** `1 known case`, `2-5`, `>5` or `unknown`, from Orphanet prevalence (`en_product9_prev.xml`): if a "Cases/families" record has a mean value, tier = 1 / 2-5 / >5 by that value; otherwise a documented point or birth prevalence class of at least 1-9 per 1 000 000 implies more than 5 cases worldwide (`>5`); classes "<1 / 1 000 000", "Unknown", "Not yet documented" or no record give `unknown`. OMIM clinical synopses are not in the pinned sources and are not used. (ii) **Borrowing pool** = other benchmark diseases with the same plurality `mechanism_class_general` (genes with a cached annotation; none -> no pool). (iii) **Empirical-Bayes weight**: DerSimonian-Laird between-disease variance tau^2 of the out-of-fold (DEV) or frozen-model (HELD-OUT) AUCs within the pool class, sampling variance of each AUC from the Hanley-McNeil formula, weight w = tau^2 / (tau^2 + v_i) of the disease's own AUC versus the class mean (reported, never fed back into any model). (iv) **Design that would close the gap**, a fixed text by tier: `1 known case` -> single-patient (N-of-1) crossover design with pre-specified outcome, borrowing from the mechanism class through a hierarchical model; `2-5` -> series of N-of-1 trials analysed jointly with a Bayesian hierarchical model borrowing across the class; `>5` -> multi-centre single-arm study against a natural-history control; `unknown` -> first ascertain the number of patients through registries. No dose, drug choice or recommendation is produced.

## 7. Leakage audit (reported in `validation_v4.json`)

Computed checks, each with its outcome: (L1) the PrimeKG prior cache (`primekg_prior_edges.json`) was built from the relations `drug_protein` and `protein_protein` only (the filter is asserted in `engine_score.load_prior_graph`; the check counts how many of the benchmark's PrimeKG disease names appear as a key or value in the cache - must be 0 - and the cache's keys are drugs and genes); no `indication`, `contraindication` or `off-label use` edge is read by any feature; (L2) gold 3 (openFDA labels) is never a feature nor a training label; the embedding texts contain no label or PrimeKG text; the code path that builds features has no gold argument (checked by signature in the tests); (L3) no group is in both DEV and HELD-OUT, and no group straddles training and test of any cross-validation fold; (L4) HELD-OUT labels are never read before the freeze; (L5) overlap between HELD-OUT and DEV diseases in genes and in retrieved papers (counts and fractions: same abstract in both sets is possible when genes are shared; reported, not removed); (L6) ATC-class masking in CV as in v1.

## 8. UNVALIDATED values (registry, written into the output with status and conservative default)

alpha 0.1; stratum cut points 50 / 500; embedding sentence count 3 and 128 tokens; MLP hyper-parameters; the tier rule for prevalence classes; the 1-SE margin as the inclusion margin (fixed by the task, not validated on held-out); coverage validity margin. Each is stored with `status = UNVALIDATED` and the `conservative_default` equal to the value used, and a test fails if any lacks a status or default. None is promoted by this run.

## 9. Freeze records
(appended while running; each entry is written before the step that depends on it)

### Freeze record 1 - benchmark and split (written by `engine_v4 benchmark`, before any feature or score exists)

- `out/engine/benchmark_v4.json` SHA-256 `3e330ffe0c94f2256c16f366f21d4a5115fec364a4459a16fdc18d42183ad419`; 145 diseases (45 from benchmark_v2 + 100 new) in 96 groups; target of 150 met: False.
- DEV 103 diseases, SHA-256 of the sorted id list `3623a3a6db60bc727e9942b45496a2347046e737c82f0cd859e86935aa599883`.
- HELD-OUT 42 diseases, SHA-256 of the sorted id list `2719907de5ab927a17ecb0d3714fdd4d922e60220f81d64ca86e3348db5f0bb9` (at least 40: True).
- Retained and dropped by first failing rule: `{"dropped": {"fewer than 1 genes": 9, "gold drug is an enzyme replacement / gene / cell therapy": 12, "no HPO or Orphanet disorder matched by identifier, name or synonym": 631, "not a rare genetic disease": 464, "umbrella term": 27}, "eligible_new": 100, "retained_from_benchmark_v2": 45}`.

### Freeze record 2 - model (written by `engine_v4 dev`, before any HELD-OUT score exists)

- Final model by the 1-SE rule on DEV: `E` (E minus B0 mean paired group AUC 0.0085, SE 0.0056, included: True); MLP is ablation only.
- Conformal K (alpha 0.1): DEV cross-conformal positive coverage 0.735 (0.9 minus SE = 0.787), valid: False.
- `out/engine/engine_v4_frozen.json`: SHA-256 of the model body `0543eb2cf094aa2ab213cc0cfec4345bc7637f9d0e55c31fba647fe4b8afedb9`; file SHA-256 `09e3a8e8f5a32ea79d455d69439a2edd86a96058fc85cd5a62be1c13b398e673`.
- Embedding model revisions: `{"biomedbert": "e1354b7a3a09615f6aba48dfad4b7a613eef7062", "sapbert": "090663c3ae57bf35ffe4d0d468a2a88d03051a4d"}`.
- Clarification applied in code: the drug text takes sentences that mention the drug from retrieved papers that contain a disease gene symbol (paper level, as the v1 paper count).
- The held-out command verifies this hash and the split hashes, writes `out/engine/heldout_v4_scored.lock` and refuses a second run.
