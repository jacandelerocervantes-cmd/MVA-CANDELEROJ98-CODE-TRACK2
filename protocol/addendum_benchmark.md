# Addendum — multi-disease benchmark (pre-registered)

Written on 2026-10-06 at about 20:30 UTC and hashed (`addendum_benchmark.sha256`) **before** any proximity was computed for any benchmark disease and before the benchmark disease list was fetched. Before this point only the Open Targets GraphQL schema was inspected (introspection and a few test queries on atopic eczema, MONDO_0004980, to learn how the API behaves). Nothing here modifies `frozen_protocol.md`; all proximity parameters are those of the frozen protocol and nothing is tuned after seeing results. The existing single second-disease check (Fanconi anemia, deviation 3) failed; this benchmark is the broader test requested in response.

## API findings that shape the design (Open Targets Platform v4 GraphQL, https://api.platform.opentargets.org/api/v4/graphql)

- There is no query that lists all diseases. `Query.drugs(chemblIds)` returns, per drug, `indications { count rows { maxClinicalStage disease { id name } } }`; `maxClinicalStage == "APPROVAL"` means approved for that indication. `rows` is complete (count equals the number of rows). Drug ids are ChEMBL ids, the same as `drug_id` in `data/universe.json`.
- `Disease.associatedTargets(datasources, page)` returns targets with an overall association `score` (sorted by score, descending; page size up to 1000). Disease-side `drugAndClinicalCandidates` exists but is not used.
- **Circularity risk.** The default association score includes the datasource `clinical_precedence` (known drugs from ChEMBL and others). Using it would put the targets of the known drugs inside the disease module, making recovery of those drugs trivial. The benchmark therefore sets `{id: "clinical_precedence", weight: 0, propagate: false}` in `datasources`. In the API, datasources not listed keep contributing, with weight 1 (observed: scores differ from the platform default weights, so scores are not the web-portal scores). `propagate: false` means direct associations only (no evidence propagated up from descendant diseases).

## (a) Disease selection rule (deterministic)

1. Pool: every disease `D` that is the `disease.id` of an indication with `maxClinicalStage == "APPROVAL"` of at least one drug in `data/universe.json` that is scorable (see c), excluding ids with prefix `HP_` (phenotype / symptom terms, not diseases).
2. `D` qualifies if all hold, evaluated at the main STRING threshold 700: (i) at least 5 universe drugs scorable at 700 have an APPROVAL indication for `D`; (ii) at least 20 targets have overall score >= 0.5 (association call above); (iii) at least 10 of those targets map to nodes of the STRING 700 giant component (rule b).
3. Take qualifying diseases sorted by ascending id (plain string order of the Open Targets id, e.g. `EFO_...` before `MONDO_...`) and keep the first 30 (all of them if fewer qualify). Because of the string order, the sample is expected to consist mainly of `EFO_` terms; this is a known quirk of the rule and is reported, not corrected.

## (b) Disease module

Targets with overall score >= 0.5 (rule above). Open Targets target `approvedSymbol` is mapped to STRING nodes through the same symbol -> node map as the rest of the project (`Interactome.nodes_for_symbols`). Minimum 10 mapped nodes. For the sensitivity thresholds 400 and 900 the module is re-mapped on that network; a disease with fewer than 10 mapped nodes there is not evaluable at that threshold and is excluded from that threshold's summary.

## (c) Ground truth and scorable universe

Positives: universe drugs with an indication whose `disease.id == D` and `maxClinicalStage == "APPROVAL"` in `Drug.indications` (direct indication only, no ontology propagation). Restricted to drugs scorable on the network at hand (at least one human target mapped in the giant component of that STRING threshold, exactly as `scoring.score_drugs` decides). Negatives: all other scorable drugs (some may be unlabeled true indications; this biases AUC towards 0.5 and is accepted). At thresholds 400/900 a disease is evaluable only if it keeps at least 3 positives.

## Proximity (frozen parameters, unchanged)

STRING v12.0 physical, d_closest, 1000 degree-matched random draws, degree bins of at least 100 nodes, base seed 20261006, drug `seed_offset` = index in the universe list. The null of each disease is built once per threshold (`proximity.build_null`) and drugs are scored with `scoring.score_drugs`. Score of a drug for a disease = -z (more negative z = closer = more likely indicated).

## (d) Metrics

Per disease and threshold: AUC (Mann-Whitney, ties count 1/2) of positives vs. all other scorable drugs; also the mean and median percentile rank (rank by ascending z divided by number of scorable drugs; lower = better) of the positives. Across diseases: mean AUC; 95% percentile bootstrap CI over diseases (10000 resamples of diseases, seed 7); one-sided permutation test (1000 permutations, seed 13): within each disease the positive labels are shuffled among scorable drugs inside degree-matched bins (10 equal-count bins of the drug's mean degree of its targets in the network; each bin keeps its number of positives); the statistic is the mean AUC across diseases and p = (1 + #{permuted >= observed}) / (1 + 1000).

## (e) Decision rule (fixed in advance, applies to the main threshold 700 only)

The method is declared **"supported across diseases"** only if mean AUC >= 0.60 **and** the lower bound of the 95% bootstrap CI > 0.50. Otherwise it is declared **"not supported across diseases"** and the drug ranking stays exploratory. The permutation p-value is reported but is not part of the rule.

## (f) Sensitivity and consensus

The same metrics are reported at STRING thresholds 400 and 900 (descriptive; they do not change the decision). Consensus: for each disease, over drugs scorable at all three thresholds and only for diseases evaluable at all three, each drug's mean of its z-rank at 400, 700 and 900 (rank 1 = most negative z); AUC of the positives vs. the others on that mean rank, then the same cross-disease summary (mean AUC, bootstrap CI, permutation test with degree bins from the 700 network).

## Not done

No parameter, threshold, disease, drug or seed is changed after seeing results. If the method fails the rule, that is the reported result. Known limitations stated in advance: diseases share targets and drugs (not independent); Open Targets associations contain literature-derived evidence (`europepmc`) that can mention known drug targets, a residual leak that is not removed; the disease sample is biased by the id-order rule; proximity does not predict direction of effect.
