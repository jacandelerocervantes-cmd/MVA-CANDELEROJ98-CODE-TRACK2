# Addendum: application of the frozen engine (B0) to mosaic variegated aneuploidy and Fanconi anaemia

Written before the run. Label of every output: **application of the frozen engine; exploratory; not a candidate list.**

## 1. Nature of this step
This is an APPLICATION of the frozen B0 engine, not a validation. It makes no claim of efficacy and no claim that any listed drug
helps any patient. No threshold is evaluated, no hypothesis is tested, nothing is fed back into the engine. Hypotheses only: no
dosing and no treatment recommendation.

## 2. Diseases and genes
* Mosaic variegated aneuploidy (MVA): the genes of tier `curated` in `out/mva_landscape.json` (BUB1, BUB1B, BUB3, CENATAC, CEP57,
  MAD1L1, TRIP13). Genes of tier `proposed, not curated` are not used for scoring; they are listed in the output.
* Fanconi anaemia: the genes of tier `curated` in `out/landscape_fanconi_anaemia.json`; same rule.
* Mechanism classes are the classes recorded in the two landscape files. Only patient-free information is used (gene symbols and
  disease names). No genomic coordinate, cDNA change or sample identifier is read, printed or stored.

## 3. Scoring
* Model: the frozen head B0 (8 literature and knowledge-graph-prior features, frozen standardisation, coefficients and Platt
  parameters) from `out/engine/engine_v4_frozen.json`, verified by its recorded hash. No refit. The E and MLP heads are not used.
* Features are computed from the retrieved Europe PMC abstracts (gene-symbol query, publication year <= 2024, at most 300 papers,
  cursorMark, <= 4 requests/s, cached under `data/cache_engine/`) exactly as for the benchmark diseases, plus the cached
  knowledge-graph prior. No label (gold 1, 2 or 3) of either disease is read or exists in the benchmark used by the engine (neither
  disease is in the v4 benchmark).
* Score = B0 decision value; the Platt probability is reported beside it and is not a probability of efficacy.

## 4. Candidate set and filter chain
Universe (ChEMBL drugs of `data/universe.json`) -> current FDA label (openFDA label match by generic name, cached; plus a Drugs@FDA
check, marketing status Prescription, for drugs without a label match inside the top 100 of the universe ranking, as in the
Wilms case) -> ChEMBL molecule type `Small molecule`. Counts at each step are computed. A service error is never read as approved.
Drugs of the oncology-cytotoxic group (ATC L01, or a member of the caution-list classes tubulin, mitotic-checkpoint, genotoxic
chemotherapy, proteasome) carry the project's chemotherapy caution, and their label hazard sentences are quoted verbatim.

## 5. Outputs per disease
1. Top 20 of the filtered ranking, with score, Platt value, rank, ties, literature counts, approval state, molecule type, caution flag.
2. Filter-chain counts.
3. Per-drug age-band label support for the bands 0-4, 5-11 and 12-17 years plus unknown age (never assumed), by the project's
   `band_status` rule; the label sections are quoted verbatim (paediatric-use sentence, hazard sentences). `not_labelled`, `verify`
   and `age_not_stated` are valid states; where the classifier is heuristic the entry says `manual review required`.
4. Mechanism class x age-band table: for each landscape mechanism class, the five small-molecule FDA drugs with most retrieved
   papers co-mentioning the drug and a gene of the class (paper counts shown), and their band states.
5. MVA only: the position of the nine drugs of the report (sirolimus, everolimus, metformin, sodium phenylbutyrate, trametinib,
   selumetinib, bortezomib, chloroquine, bosutinib) in the universe, FDA-label and filtered rankings. Comparison only; the report's
   selection used other lines of evidence and this is not a validation of either.
6. Evidence-gap block: case-count tier from Orphanet prevalence (`tier_from_records`, cached file only; `unknown` if the file or the
   code is absent). For MVA the proband is stated to be one known case. Study designs from the project's borrowing logic
   (`DESIGN_BY_TIER`, class-level mean AUC and between-disease variance of the mechanism class pool from
   `validation_v4.json::evidence_gap_hierarchical_borrowing`) that could close each gap: case count, class-level uncertainty, and
   each age band without labelled support. Computed; no dosing; no recommendation. The disease's own AUC does not exist (no
   label), so the weight on it is not computed.
7. Literature-volume covariate: papers retrieved, papers with a gene symbol, Europe PMC hit count and its engine stratum, number of
   universe drugs with at least one co-mention and total co-mentions.
8. Flag `low literature: ranking driven by the prior`, fixed here before the run from the dilution result
   (`validation_v4_dilution.json`, `validation_v4_ablation.json`): when the literature is thinned to about 5-9 papers per disease the
   B0 gain over popularity (0.05-0.06 AUC) is at the level of the prior-only model (0.04 AUC), i.e. the ranking is carried by the
   knowledge-graph prior; the gain grows with volume (0.09 at about 29 papers). Rule: flag when the number of retrieved papers
   carrying a gene symbol is at most the ceiling of the mean volume of the 10 % thinning (29 papers, read from the file) or the Europe PMC
   hit count is below 50 (stratum L of the engine).

## 6. Limits (reproduced in every output)
* The engine was validated on diseases with known drugs (benchmark gold from approved indications and trials); it is exploratory
  where no drug is known.
* On Wilms tumour the engine failed to rank the known drugs in its top 20; a high rank is therefore not evidence of efficacy and a
  low rank is not evidence against.
* Literature counts are co-mentions among at most 300 abstracts; they are not effect evidence. The knowledge-graph prior is
  target-based.
* ChEMBL molecule type and openFDA label presence are used as recorded. Label age support is label text; the heuristic classifier
  is not hand-verified except where `config/pediatric_label_review.json` has an entry.
* Oncology-cytotoxic drugs carry serious label-documented toxicity (marrow suppression, second malignancy, embryo-fetal harm); in
  Fanconi anaemia, DNA-crosslinking agents are of particular concern because of cellular hypersensitivity; this is a flag, not advice.
* The nine drugs comparison is not a validation. Nothing here is a candidate list, a dosing statement or a treatment recommendation.
