# Addendum - benchmark v2 for the literature-grounded evidence engine (documented extension)

Written on 2026-10-15, **before** `src/track2/engine_benchmark_v2.py` was implemented or run, and before any engine score exists for any v2 disease. It extends, and does not replace, `protocol/addendum_engine.md` (benchmark v1). Nothing here modifies `frozen_protocol.md`, `src/track2/engine_benchmark.py`, `out/engine/benchmark.json`, `protocol/addendum_engine.md` or `src/track2/engine_score.py`.

## 0. Status of v1 and v2

- The **v1-only result (17 diseases, gold 1) remains the primary pre-registered result** of the engine validation. The success decision (three thresholds of section 4 of `addendum_engine.md`) is taken on v1 alone.
- v2 is an extension that addresses four limitations of v1. It is reported **in addition** to v1, never instead of it, and never used to rescue or overturn the v1 decision.
- The v2 analysis reports (i) the **union benchmark** (v1 diseases retained plus new diseases) and (ii) the **new diseases alone as a held-out set** (they were not seen when the engine, its components, its thresholds and its endpoints were fixed).
- **All v1 thresholds, endpoints, baselines, bootstrap settings, grouped cross-validation, ablations and reporting rules are unchanged** (sections 2 to 8 of `addendum_engine.md`). v2 adds no threshold and removes none. Gold 1 stays the endpoint gold for the thresholds; gold 3, the union and the intersection of gold 1 and gold 3, and the corrected gold 2 are descriptive sensitivity analyses that do not change any decision.
- **Nothing may be changed after seeing results**: not the matching rules, not the gold definitions, not the mechanism rule table, not the cap, not the seed. Any later deviation goes to `protocol/deviations.md` and does not replace this pre-registration.

## 1. Limitation 1: disease matching by identifiers and synonyms (546 candidates lost in v1)

v1 matched a PrimeKG disease name to HPO or Orphanet only by whole-word name containment and required an exact-name Orphanet match for "rare genetic". v2 adds matching by identifiers and synonyms; **every other inclusion and exclusion rule of v1 is kept** (listed in 1.4).

### 1.1 Identifiers
- PrimeKG disease nodes (`data/sources/primekg/nodes.tab`, `node_type = disease`) carry one MONDO number in `node_id` or several joined by `_` (`MONDO_grouped`); the pool is the same as v1 (curated candidates of `empirical_literature_control.CANDIDATES` plus PrimeKG diseases with at least 2 approved universe drugs linked by an `indication` edge). All PrimeKG nodes with the candidate's name are used.
- MONDO to OMIM / OMIMPS / Orphanet: `xref:` lines of `data/sources/mondo/mondo.obo` (release read from its `data-version` header). The OBO Foundry SSSOM mapping set is a re-serialisation of the same cross-references; it is not used. Only xrefs whose `source` qualifier contains `MONDO:equivalentTo` count (obsolete, related or `equivalentObsolete` xrefs are ignored). Obsolete MONDO terms are ignored.
- HPO (`phenotype.hpoa`, `genes_to_phenotype.txt`) disorders match by identifier when their `OMIM:n` or `ORPHA:n` id equals a mapped xref. Orphanet (Orphadata `en_product6`/`en_product9_ages`) disorders match by ORPHA code. An OMIMPS (phenotypic series) xref is recorded but does not by itself match a disorder.

### 1.2 Synonyms
- Name set of a disease = PrimeKG name, MONDO `name`, MONDO synonyms with scope EXACT (abbreviation-typed synonyms and the `AMBIGUOUS`/`DEPRECATED` types are dropped), and Orphanet synonyms (`en_product1.xml`, `SynonymList`) of the identifier-matched Orphanet disorders. Names shorter than 5 characters are dropped. Names are compared in the v1 normalised form (`norm_name`: lower case, brackets/hyphens/commas removed, anaemia = anemia).
- HPO: the pinned `phenotype.hpoa` carries one name per disease record and no synonym list, so HPO "synonyms" enter only through the identifier mapping above; the name of every identifier-matched HPO or Orphanet disorder is added to the name set of the disease. This is a documented limit.

### 1.3 Matching rule
- Matched disorders = (identifier-matched HPO/Orphanet disorders) union (HPO/Orphanet disorders whose name contains any name of the name set as whole words, the v1 `term_regex`). More than 30 matched disorders = umbrella term, excluded (v1 `MAX_MATCHED`).
- "Identity disorders" = identifier-matched disorders plus disorders whose normalised name EQUALS a normalised name of the name set. **Rare genetic** (replaces the v1 exact-name rule) = at least one identity disorder that is Orphanet-catalogued (ORPHA code in Orphadata, an ORPHA id in HPO, or an Orphanet equivalence xref in MONDO) and has Mendelian, mitochondrial or sex-linked inheritance (HPO aspect I or Orphadata). Curated candidates of the earlier proxy keep the v1 looser alternative (any matched Orphanet-catalogued disorder with Mendelian inheritance).
- Gene list = union of HPO and Orphanet gene symbols of all matched disorders (v1 rule, placeholder symbols removed).

### 1.4 Rules kept from v1 without change
Seed **20261015**; at least 2 genes; at least 2 universe drugs in gold 1; MVA / premature chromatid separation and Fanconi anaemia excluded; a disease is excluded when any gold 1 drug has molecule type Enzyme, Gene or Cell (the rule is applied to gold 1 only, as in v1; biologic types in gold 3 are recorded, not used for exclusion); diseases with an identical gold 1 set are collapsed to the alphabetically first (a new disease whose gold 1 set equals that of a retained v1 disease is dropped in favour of the v1 disease); candidate set of every ranking = all drugs of `data/universe.json`.

### 1.5 Composition and size
- The 17 diseases of v1 are retained unchanged (name, gene list and gold 1 are read from `out/engine/benchmark.json`, read-only; identifiers are extended with the MONDO, OMIM and ORPHA ids found here).
- New diseases = candidates not in v1 that pass all rules. If more than 43 qualify, 43 are drawn with `random.Random(20261015)` from the alphabetically sorted eligible names (total 60); otherwise all qualify. The target is at least 40 diseases in total; if fewer qualify, the shortfall is reported and no rule is relaxed.
- The selection log lists every evaluated candidate with the first failing rule and the counts of exclusions by reason.

## 2. Limitation 2: gold 3, independent of PrimeKG

- **Gold 3** = universe drugs (`data/universe.json`, approved) whose **current FDA label** names the disease. Source: openFDA `drug/label.json`, field `indications_and_usage`; query `openfda.generic_name:"<name>" OR openfda.substance_name:"<name>"`, up to 5 most recent labels per drug (`sort=effective_time:desc`), at most 3 requests per second, responses reduced to name fields, effective time and the indications text and cached under `data/cache_engine/openfda/`.
- **Matching**: a drug-disease pair is in gold 3 when the indications text of at least one retrieved label contains one name of the disease's name set (section 1.2, names of at least 5 characters, plus the PrimeKG name) as a **whole phrase** (v1 `term_regex`: case-insensitive, whole words, space/hyphen tolerant, anaemia = anemia). No negation handling, no stemming, no fuzzy matching: a label that mentions the disease only under limitations of use is counted. Subtypes are covered because a phrase is also matched inside longer labels ("hemophilia" matches "hemophilia A").
- Reported per disease and overall: gold 1, gold 3, their **union**, their **intersection**, and the sizes of each. Gold 1 drugs approved for the indication but absent from labels (e.g. different wording) and gold 3 drugs absent from PrimeKG are not reconciled.

## 3. Limitation 3: gold 2 rebuilt from experimental arms

- Studies: ClinicalTrials.gov API v2, `query.cond` = the PrimeKG disease name, interventional, phase 2 or later (PHASE2, PHASE3, PHASE4; PHASE2/PHASE3 counts), start date strictly after 2018-12-31 (as v1).
- A universe drug counts **only if** it is an intervention that is listed in the `interventionNames` of at least one arm group of type **EXPERIMENTAL** of such a study (match by the v1 drug-name regex on the intervention name and other names; non-drug intervention types device, procedure, behavioural, radiation and diagnostic test are skipped). Drugs that appear only in PLACEBO_COMPARATOR, ACTIVE_COMPARATOR, SHAM_COMPARATOR, NO_INTERVENTION or OTHER arms, or only in studies without arm groups, are excluded and counted in the log.
- Each gold 2 entry records NCT id, start year and arm type (plus arm label and phases). Evidence cutoff year Y = 2018, as in v1.

## 4. Limitation 4: general mechanism classes

- `src/track2/mechanism_general.py` holds a deterministic rule table over UniProt keywords and Gene Ontology term names of the reviewed human UniProt entry (the same annotation source and cache as `disease_landscape`, fetched without editing it). Classes: lysosomal / vesicular trafficking, ion channel / transporter, metabolic enzyme, receptor / signalling, transcription / chromatin regulation, structural / cytoskeletal, mitochondrial, DNA repair / replication, cell-cycle / mitotic checkpoint, RNA processing / translation, extracellular / secreted, other.
- Rule: each rule is a regular expression with a weight (strong 2, weak 1) applied to every distinct term; a class score is the sum of weights of its matched terms, capped at 6 (so that pleiotropic annotation cannot make a generic class win by term count alone); the class with the highest score wins, ties are broken by the fixed order of the table, and a score below 2 gives `other`. The table contains no gene or disease names. Its vocabulary was written from the public UniProt keyword and GO vocabulary and checked once against the cached annotation of genes of earlier landscape runs (term vocabulary only, no gene-level tuning); it is frozen with this file. It is fixed before any v2 output is inspected for scoring and is not tuned afterwards.
- The field `mechanism_class_general` is added to every gene of the v2 benchmark. The v1 field `mechanism_class` (from `disease_landscape.classify_terms`) is kept next to it. `disease_landscape.py` and its outputs are not edited.

## 5. Output and reporting

- `out/engine/benchmark_v2.json` (sorted keys): for every disease identifiers, genes with `mechanism_class_general`, gold 1, gold 3, gold 2 (corrected), the origin (`v1` or `v2_new`), and the selection log with exclusion counts by reason. Every number comes from a source file or API response.
- The v2 analysis, when it is run, will report: the union benchmark; the new diseases alone as a held-out set; gold 1, gold 3, union and intersection; temporal AUC against corrected gold 2; the number of diseases with empty gold 3 or gold 2. Results are reported whatever they are.
- Limits fixed in advance: gold 3 is a text match without negation handling and depends on openFDA coverage of current labels; ClinicalTrials.gov condition search is by disease name only; HPO has no synonym table in the pinned files; the new diseases remain a convenience selection of PrimeKG diseases with at least 2 indicated approved drugs.
