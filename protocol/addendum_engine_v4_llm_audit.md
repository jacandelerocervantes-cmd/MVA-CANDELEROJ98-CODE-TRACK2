# Addendum - engine v4, LLM-read evidence: AUDIT of the low-coverage claim (DEV only, read-only)

Written on 2026-10-16, **before** `src/track2/engine_v4_llm_audit.py` was run. It extends and does not modify `addendum_engine_v4_llm.md`, `addendum_engine_v4_llm_exploratory.md`, `addendum_engine_v4_llm_guideline_v2.md`, any module, label file or output. Facts known when this was written: the parent protocol, the code of `engine_v4_llm.py`, the manifest counts (DEV: 4442 candidate pairs, 2427 with a sentence, 7637 verified items) and the exploratory coverage numbers (non-zero relevance features in under 1 % of pairs at 3 % of the literature). No audit number existed.

## 0. Status and rules

- **Claim audited** (stated by the assistant, doubted by the user): *at low literature volume there are almost no relevant sentences, so LLM-read evidence cannot be tested.* The user's suspicion: the low coverage is an artefact of how sentences are retrieved and selected, and of the labels, not a property of the literature.
- **Read-only, DEV only.** No HELD-OUT label, HELD-OUT gold, HELD-OUT shard or HELD-OUT annotation file is read; no stage of any pipeline and no `run_all` is run; no lock is touched. Multi-split files (candidate file, pair index, benchmark) are filtered to DEV diseases at load time and the rest is discarded. The audit never opens `ho-*` items and never uses the `heldout` labels.
- **No new annotation.** Nothing is sent to an LLM. The only new labels are my own audit readings of at most 60 sentences (Q3a); they are an **audit opinion, not ground truth**, are stored apart, and feed no feature or model.
- **Network.** Only Europe PMC, only for Q2 variant G (a direct "drug AND disease" query, first 10 papers), at most 4 requests per second (the project's limiter `engine_score.default_getter`), responses cached under `data/cache_engine/audit/`. No file is downloaded.
- Seed everywhere **20261015**. Every number is computed from the sources; nothing is typed in. Surprises that contradict the claim are reported as such. Constants of the audit (sample size 150, 60 sentences, Europe PMC page size 10, disease-term limits) are UNVALIDATED audit choices, recorded in the output.
- Stack-ml skills applied: `leakage-audit-grouped-split` (DEV only, gold flags used only to split the reporting, never as a feature; annotator knowledge of drugs noted), `weak-label-from-reports` (labels are weak, model-made; agreement and a human re-reading are reported, not assumed), `metricas-y-calibracion` (counts with the denominator always shown, rates not hidden in means, replicate spread given).

## 1. Definitions

- **Gold-positive pair**: a DEV candidate pair whose `g1` flag is true in the candidate file (gold 1, the PrimeKG indication). Other candidates: all remaining DEV pairs (top-30 B0 negatives; also pairs flagged only by gold 2 or gold 3 are in "other"; their number is reported). Expected 465 gold-positive of 4442.
- **Sample** (as the parent): f = 1.0 once; f = 0.3, 0.1, 0.03, 0.01 with 5 replicates, kept papers = `engine_v4_dilution.subsample_indices(n_retrieved, f, rank, fi, rep)`. Rates at f < 1 are means over the 5 replicates, with min and max.
- **Counted sentence** (parent): a verified sentence of the first 3 kept papers (retrieval order) that mention the drug, at most 2 sentences per paper.
- **Annotator relevance** `yes` by *either* annotator / by *both*; **benefit direction** = the parent `consensus()` direction (2 of 3 of Opus, Sonnet, rule extractor) equals `benefit`.

## 2. Questions and metrics

### Q1. Funnel on DEV pairs, gold-positive versus other candidates, at f = 1.0 and at the thinning fractions

Stages (a pair is counted at a stage if it has at least one qualifying item; each stage is a subset of the previous one except where stated):

- S0 candidate pairs;
- S1a >= 1 kept retrieved paper whose title or abstract names the drug **anywhere** (the production drug regex over the full text, no gene gate, no caps);
- S1b >= 1 such paper that also contains a disease gene symbol (the v1 paper gate; a sub-stage of S1a);
- S2 >= 1 sentence extracted (the parent pair index restricted to the first 3 kept papers);
- S3 >= 1 **verified** sentence (the item exists in the DEV shards);
- S4 >= 1 counted sentence with relevance `yes` by **either** annotator;
- S5 >= 1 counted sentence with `yes` by **both**;
- S6 >= 1 counted sentence with `yes` by **both** and benefit direction; S6' the same with `yes` by either.

Reported per stage: count, share of S0, conditional share of the previous stage (retention), for gold-positive and other pairs, per fraction (mean over replicates, min-max). The **largest loss** is the transition with the largest absolute loss of pairs and the transition with the lowest retention, named for the gold-positive pairs at f = 1.0 and at f = 0.03. Also the item-level counts (items, items with `yes`).

### Q2. Retrieval and extraction design, and loss attributable to each rule

(a) From the code and the manifest: how papers are retrieved (query text, field, year, cap, sort), how sentences are selected (drug-name rule, gates, field), the caps, and the unit of annotation, stated with the file and function where each rule lives.

(b) **Relaxation experiment** on a seeded sample of **150 gold-positive DEV pairs** (`random.Random(20261015).sample` of the sorted gold-positive pair list), re-running the extraction **without any labelling**. Variants, each relaxing one rule of the baseline R0 (f = 1.0, first 3 gene-gated papers that mention the drug, first 2 sentences per paper, exact INN with the production regex):

- R0 baseline (must reproduce the stored pair index; the match rate is reported as a check);
- R1 sentence cap lifted (3 papers, all sentences);
- R2 paper cap lifted (all gene-gated retrieved papers, 2 sentences per paper);
- R3 both caps lifted (all gene-gated retrieved papers, all sentences);
- R4 R3 and the gene gate removed (all retrieved papers);
- R5 R3 and drug **synonyms** (section 3) in the drug match;
- R6a direct Europe PMC query `drug AND disease` (title/abstract, year <= 2024, abstract present), first 10 papers, INN only, all sentences;
- R6b the same query and match with the synonyms;
- R7 union of R4, R5 and R6b (everything relaxed).

Metrics per variant, over the 150 pairs: sentences naming the drug **and** the disease in the same sentence (disease terms: the benchmark name only, and the name plus the MONDO / Orphanet / OMIM name set stored in the benchmark identifiers) - mean, median, total, and number (share) of pairs with >= 1; sentences naming the drug whatever the disease; papers whose title+abstract names both; sentences naming drug and disease with a rule-extractor direction of `benefit` or `harm` (a label-free proxy of an evaluative sentence); the number of papers. For R0 to R5 also at f = 0.03 (mean of the 5 replicates). The loss attributable to a rule is the gain in the mean number of drug-and-disease sentences per pair and in the share of pairs with >= 1 when that rule alone is relaxed.

### Q3. Label audit

(a) Gold-positive DEV pairs with at least one annotated sentence where **every** annotated sentence has relevance `no` or `unclear` from **both** annotators (no `yes` by either). Up to **60** sentences, one random sentence per randomly chosen pair (seed 20261015), read by me and classified in a table as: *truly uninformative*; *informative but the label looks too strict*; *label looks right*; one-line reason each. Opinion only. Summary: counts per class.

(b) Cross-tabulations Opus x Sonnet of relevance (3 x 3) and direction (4 x 4) on DEV items of gold-positive pairs versus other pairs, with n, raw agreement and Cohen's kappa for each group, and the share of `yes`.

(c) Sensitivity of coverage to the relevance rule: the share of gold-positive and other pairs with >= 1 counted sentence under `either yes`, `both yes`, `Opus yes`, `Sonnet yes`, each with and without benefit direction, at f = 1.0, 0.1 and 0.03.

(d) v2 pilot: the 800 items of the 8 pilot shards (002, 004, 006, 025, 026, 028, 057, 062). Per annotator, cross-tab v1 versus v2 of relevance (and of direction); the Opus-Sonnet disagreement on relevance under v1 and under v2 (raw disagreement and kappa), the number of items whose disagreement is removed, kept and newly created, and the `yes` counts (either / both) under v1 and v2.

### Q4. Verdict

A plain-language conclusion on whether the low coverage is (i) a property of the literature, (ii) an artefact of retrieval and extraction (naming which rule and how much), or (iii) label strictness, from the numbers of Q1-Q3; plus the cheapest single change expected to raise coverage the most, with an estimated gain computed from the Q2 variants, as a **recommendation only** (nothing implemented, nothing annotated). Decision aids fixed now: a rule is called a major artefact if relaxing it alone at least doubles the mean number of drug-and-disease sentences per pair or raises the share of pairs with >= 1 by at least 20 percentage points; label strictness is called major if "either yes" at least doubles the number of pairs with a counted `yes` sentence relative to "both yes" **and** my reading classifies at least a quarter of the read sentences as informative-but-too-strict. These thresholds are UNVALIDATED audit conventions.

## 3. Drug synonyms (what the project actually holds)

`data/universe.json` stores one name per drug (the ChEMBL preferred name, upper case) and no synonym, brand name or salt field; the drug regex is the exact whole-word match of that name (case-insensitive, names of >= 4 characters). The synonym sources available inside the project are therefore limited, and used as follows, declared before the run: (1) the **other names** of ClinicalTrials.gov interventions in `data/cache_engine/ctgov/` whose intervention name equals a universe drug name; (2) the `molecule_synonyms` of the ChEMBL molecule records that exist in `data/cache/` (a partial set); (3) rule-based **salt and form stripping** of the universe name (a fixed list of salt words, for example "hydrochloride", "sodium", "sulfate") and the stripped name's inverse (the universe name itself); (4) the hyphen/space variant. Brand names are included only where (1) or (2) supply them. The coverage of this synonym set (share of the 150 pairs with at least one synonym beyond the INN) is reported: **R5 and R6b are a lower bound** on what a real synonym table would give. This limitation is stated now and not repaired.

## 4. Leakage and honesty notes

Gold flags split the reporting only; no feature, label or sample is chosen with them except the Q2 sample (gold-positive pairs by definition) and the Q3a pool. The annotators know what drugs treat which disease (parent S4); the audit does not remove that. The Q3a reading was made knowing the pair is gold-positive (it is selected as such), which biases my reading toward finding sentences informative; that bias is stated in the report, and the table includes the reason so it can be re-read. Random thinning is a stress test of volume, not a rare-disease measurement. Q2 relaxations measure the number of **candidate sentences**, not whether they are relevant: they show the capacity of retrieval and extraction, not a gain in signal.

## 5. Additional audit questions (added 2026-10-16 at the request of the user, still before the audit module was run)

Same hard rules: DEV only, no edit of existing files, no HELD-OUT, Europe PMC at most 4 requests per second, responses cached under `data/cache_engine/audit/`, nothing downloaded, nothing labelled. They extend Q2 and Q4.

- **Q2-L (language).** The retrieval query (`engine_score.build_gene_query`) has no `LANG` filter and the drug match is the exact English INN in title+abstract. On a seeded sample of **30 DEV diseases** (`random.Random(20261015)`, diseases that have at least one gold-positive pair) the first 300 records of the project's gene query are re-requested (pages of 100, `cursorMark`, field `language` of the core record) and compared with the cached 300 by id (an agreement check). Reported: records by language and the share non-English among the retrieved; the hit count of the whole gene query with and without `LANG:eng`. For the gold-positive pairs of the 150-pair sample (Q2-Q below): hit count of `drug AND disease` overall and with `NOT LANG:eng`, the languages of the first 20 records, and the share of non-English records among them whose title+abstract matches the INN. Drug names written in a non-English form cannot be searched because the project holds no translation table; this is stated and the number of non-English records is the upper bound of that loss. Language is called **material** if non-English records exceed 5 % of the retrieved records or of the `drug AND disease` hits (UNVALIDATED convention).
- **Q2-Q (recall of the current retrieval against a direct query).** For each of the 150 gold-positive pairs of Q2: one Europe PMC call `drug AND disease (name or synonyms)`, year <= 2024, abstract present, first **20** records (the first 10 serve variant R6) and the hit count. Reported: share of the first-20 record ids that are among the 300 retrieved papers of the disease (recall), its distribution over pairs, and the share of pairs with recall 0. To separate the two possible causes, one further call per pair gives the hit count of the same direct query additionally restricted to the gene clause of the disease (`AND (gene clause)`): hits that the gene query could never retrieve (query-design loss) versus hits that it matches but the 300-record cap may cut (cap loss, bounded by the ratio of the retrieved set size to the gene-query hit count, also reported per disease).
- **Q2-N (drug naming).** Over all gold-positive DEV pairs: the share of gold-positive drugs whose project name is not a plain INN (the name carries a salt or form word, a formulation, a bracket or a digit), the share that has at least one project synonym (section 3), and, in the retrieved 300 papers plus the direct-query records of the sample, the number of sentences that name the drug **only** by a synonym and not by the project name, as a share of all sentences naming the drug by either. Brand names and abbreviations beyond the sources of section 3 cannot be measured: the figure is a lower bound.
- **Q4** additionally states, from these numbers, whether language, query design and drug naming are material losses and which of them the recommendation addresses.

## 6. Run log and deviations (appended after the run)

- Run once on 2026-10-16 (about 47 min, mostly the sentence index and Europe PMC calls, all cached under `data/cache_engine/audit/`); HELD-OUT not touched; 16 offline tests pass.
- Q3a: only 21 of 159 gold-positive pairs with items had no `yes` from either annotator (66 sentences), so "up to 60 sentences" became 29: one per pair, then up to 2 more per pair at random (diversity cap of 3 per pair). The vitamin E pair (nutrient equal to the drug) supplies many of the uninformative sentences.
- Q2-L: records with no `language` field in the retrieved sets are all preprints (source PPR); they are reported apart and not counted as non-English. The `NOT LANG:"eng"` hit counts of the direct query include such records.
- Q2-Q/R6: the disease phrase of the direct query is the benchmark name plus up to five short synonyms; 57 of 150 pairs have zero hits, partly because rare-disease labels differ from the words in papers, so R6 is a lower bound.
- The "major artefact" doubling convention flags R2 to R5 because the baseline is 0.11 sentences per pair; read with the absolute numbers (11 to 16 of 150 pairs).
