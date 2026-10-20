# Audit of the low-coverage claim for LLM-read evidence (DEV only, read-only)

READ-ONLY, DEV-ONLY audit of the low-coverage claim (post hoc); no new annotation; opinion tables are audit readings, not ground truth. Protocol: `protocol\addendum_engine_v4_llm_audit.md`. Seed 20261015. HELD-OUT touched: False.

## Q1. Funnel of DEV candidate pairs

4442 candidate pairs: 465 gold-positive (gold 1), 3977 other (1199 of them gold 2 / gold 3 only). Counts are means over 5 replicates for f < 1 (min-max in the JSON). Each cell: pairs (share of the group).

### f = 1.0

| stage | gold-positive | retention | other | retention |
|---|---|---|---|---|
| S0_candidate_pairs | 465.0 (100.0%) | - | 3977.0 (100.0%) | - |
| S1a_drug_named_in_a_kept_paper | 163.0 (35.1%) | 35.1% | 2279.0 (57.3%) | 57.3% |
| S1b_and_paper_has_disease_gene | 159.0 (34.2%) | 97.5% | 2268.0 (57.0%) | 99.5% |
| S2_sentence_extracted | 159.0 (34.2%) | 100.0% | 2268.0 (57.0%) | 100.0% |
| S3_verified_sentence | 159.0 (34.2%) | 100.0% | 2258.0 (56.8%) | 99.6% |
| S4_relevance_yes_either | 138.0 (29.7%) | 86.8% | 545.0 (13.7%) | 24.1% |
| S5_relevance_yes_both | 77.0 (16.6%) | 55.8% | 156.0 (3.9%) | 28.6% |
| S6_yes_both_and_benefit | 37.0 (8.0%) | 48.1% | 78.0 (2.0%) | 50.0% |
| S6b_yes_either_and_benefit | 74.0 (15.9%) | - | 282.0 (7.1%) | - |
| opus_yes | 78.0 (16.8%) | - | 162.0 (4.1%) | - |
| sonnet_yes | 137.0 (29.5%) | - | 540.0 (13.6%) | - |
| opus_yes_benefit | 37.0 (8.0%) | - | 79.0 (2.0%) | - |
| sonnet_yes_benefit | 74.0 (15.9%) | - | 281.0 (7.1%) | - |
| benefit_direction_only | 82.0 (17.6%) | - | 644.0 (16.2%) | - |

gold_positive: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (302.0 pairs); lowest retention S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (35.1%).

other: largest absolute loss S3_verified_sentence -> S4_relevance_yes_either (1713.0 pairs); lowest retention S3_verified_sentence -> S4_relevance_yes_either (24.1%).

### f = 0.3

| stage | gold-positive | retention | other | retention |
|---|---|---|---|---|
| S0_candidate_pairs | 465.0 (100.0%) | - | 3977.0 (100.0%) | - |
| S1a_drug_named_in_a_kept_paper | 85.4 (18.4%) | 18.4% | 992.4 (25.0%) | 25.0% |
| S1b_and_paper_has_disease_gene | 81.0 (17.4%) | 94.8% | 966.8 (24.3%) | 97.4% |
| S2_sentence_extracted | 81.0 (17.4%) | 100.0% | 966.8 (24.3%) | 100.0% |
| S3_verified_sentence | 81.0 (17.4%) | 100.0% | 963.8 (24.2%) | 99.7% |
| S4_relevance_yes_either | 68.6 (14.8%) | 84.7% | 219.8 (5.5%) | 22.8% |
| S5_relevance_yes_both | 39.0 (8.4%) | 56.9% | 53.2 (1.3%) | 24.2% |
| S6_yes_both_and_benefit | 16.8 (3.6%) | 43.1% | 23.6 (0.6%) | 44.4% |
| S6b_yes_either_and_benefit | 34.0 (7.3%) | - | 110.0 (2.8%) | - |
| opus_yes | 39.2 (8.4%) | - | 54.6 (1.4%) | - |
| sonnet_yes | 68.6 (14.8%) | - | 218.6 (5.5%) | - |
| opus_yes_benefit | 16.8 (3.6%) | - | 23.8 (0.6%) | - |
| sonnet_yes_benefit | 34.0 (7.3%) | - | 110.0 (2.8%) | - |
| benefit_direction_only | 38.8 (8.3%) | - | 250.8 (6.3%) | - |

gold_positive: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (379.6 pairs); lowest retention S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (18.4%).

other: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (2984.6 pairs); lowest retention S3_verified_sentence -> S4_relevance_yes_either (22.8%).

### f = 0.1

| stage | gold-positive | retention | other | retention |
|---|---|---|---|---|
| S0_candidate_pairs | 465.0 (100.0%) | - | 3977.0 (100.0%) | - |
| S1a_drug_named_in_a_kept_paper | 40.2 (8.6%) | 8.6% | 427.8 (10.8%) | 10.8% |
| S1b_and_paper_has_disease_gene | 38.8 (8.3%) | 96.5% | 406.0 (10.2%) | 94.9% |
| S2_sentence_extracted | 38.8 (8.3%) | 100.0% | 406.0 (10.2%) | 100.0% |
| S3_verified_sentence | 38.8 (8.3%) | 100.0% | 404.2 (10.2%) | 99.6% |
| S4_relevance_yes_either | 32.2 (6.9%) | 83.0% | 82.8 (2.1%) | 20.5% |
| S5_relevance_yes_both | 17.0 (3.7%) | 52.8% | 20.4 (0.5%) | 24.6% |
| S6_yes_both_and_benefit | 8.2 (1.8%) | 48.2% | 7.8 (0.2%) | 38.2% |
| S6b_yes_either_and_benefit | 16.0 (3.4%) | - | 41.8 (1.1%) | - |
| opus_yes | 17.2 (3.7%) | - | 21.0 (0.5%) | - |
| sonnet_yes | 32.0 (6.9%) | - | 82.2 (2.1%) | - |
| opus_yes_benefit | 8.2 (1.8%) | - | 8.2 (0.2%) | - |
| sonnet_yes_benefit | 16.0 (3.4%) | - | 41.4 (1.0%) | - |
| benefit_direction_only | 18.0 (3.9%) | - | 90.0 (2.3%) | - |

gold_positive: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (424.8 pairs); lowest retention S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (8.6%).

other: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (3549.2 pairs); lowest retention S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (10.8%).

### f = 0.03

| stage | gold-positive | retention | other | retention |
|---|---|---|---|---|
| S0_candidate_pairs | 465.0 (100.0%) | - | 3977.0 (100.0%) | - |
| S1a_drug_named_in_a_kept_paper | 17.4 (3.7%) | 3.7% | 173.0 (4.4%) | 4.4% |
| S1b_and_paper_has_disease_gene | 16.8 (3.6%) | 96.6% | 164.0 (4.1%) | 94.8% |
| S2_sentence_extracted | 16.8 (3.6%) | 100.0% | 164.0 (4.1%) | 100.0% |
| S3_verified_sentence | 16.8 (3.6%) | 100.0% | 162.2 (4.1%) | 98.9% |
| S4_relevance_yes_either | 13.8 (3.0%) | 82.1% | 36.4 (0.9%) | 22.4% |
| S5_relevance_yes_both | 8.4 (1.8%) | 60.9% | 10.0 (0.3%) | 27.5% |
| S6_yes_both_and_benefit | 3.4 (0.7%) | 40.5% | 4.2 (0.1%) | 42.0% |
| S6b_yes_either_and_benefit | 7.0 (1.5%) | - | 16.8 (0.4%) | - |
| opus_yes | 8.4 (1.8%) | - | 10.2 (0.3%) | - |
| sonnet_yes | 13.8 (3.0%) | - | 36.2 (0.9%) | - |
| opus_yes_benefit | 3.4 (0.7%) | - | 4.2 (0.1%) | - |
| sonnet_yes_benefit | 7.0 (1.5%) | - | 16.8 (0.4%) | - |
| benefit_direction_only | 8.6 (1.8%) | - | 36.6 (0.9%) | - |

gold_positive: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (447.6 pairs); lowest retention S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (3.7%).

other: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (3804.0 pairs); lowest retention S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (4.4%).

### f = 0.01

| stage | gold-positive | retention | other | retention |
|---|---|---|---|---|
| S0_candidate_pairs | 465.0 (100.0%) | - | 3977.0 (100.0%) | - |
| S1a_drug_named_in_a_kept_paper | 11.8 (2.5%) | 2.5% | 119.6 (3.0%) | 3.0% |
| S1b_and_paper_has_disease_gene | 10.6 (2.3%) | 89.8% | 112.6 (2.8%) | 94.1% |
| S2_sentence_extracted | 10.6 (2.3%) | 100.0% | 112.6 (2.8%) | 100.0% |
| S3_verified_sentence | 10.6 (2.3%) | 100.0% | 111.8 (2.8%) | 99.3% |
| S4_relevance_yes_either | 8.8 (1.9%) | 83.0% | 24.2 (0.6%) | 21.6% |
| S5_relevance_yes_both | 3.8 (0.8%) | 43.2% | 5.2 (0.1%) | 21.5% |
| S6_yes_both_and_benefit | 1.8 (0.4%) | 47.4% | 2.0 (0.1%) | 38.5% |
| S6b_yes_either_and_benefit | 3.6 (0.8%) | - | 8.4 (0.2%) | - |
| opus_yes | 3.8 (0.8%) | - | 5.2 (0.1%) | - |
| sonnet_yes | 8.8 (1.9%) | - | 24.2 (0.6%) | - |
| opus_yes_benefit | 1.8 (0.4%) | - | 2.0 (0.1%) | - |
| sonnet_yes_benefit | 3.6 (0.8%) | - | 8.4 (0.2%) | - |
| benefit_direction_only | 3.6 (0.8%) | - | 18.8 (0.5%) | - |

gold_positive: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (453.2 pairs); lowest retention S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (2.5%).

other: largest absolute loss S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (3857.4 pairs); lowest retention S0_candidate_pairs -> S1a_drug_named_in_a_kept_paper (3.0%).

## Q2. Retrieval and extraction design

**paper_retrieval**

- diseases_at_cap: 89
- function: engine_score.fetch_papers / build_gene_query (called by engine_v4_llm._pair_job, cache only)
- mean_papers_kept_by_fraction: `{"0.01": 5.0, "0.03": 8.62135922330097, "0.1": 28.25242718446602, "0.3": 84.59223300970874, "1.0": 281.77669902912623}`
- n_retrieved_per_dev_disease: `{"max": 300.0, "mean": 281.77669902912623, "median": 300.0, "min": 46.0, "n": 103, "p25": 300.0, "p75": 300.0}`
- ordering: Europe PMC default order (no sort parameter), cursorMark paging, first 300 records with an abstract
- query_form: (TITLE_ABS:"ABCA7" OR TITLE_ABS:"APOE" OR TITLE_ABS:"APP") AND PUB_YEAR:[1800 TO 2024] AND HAS_ABSTRACT:y   (example: first 3 gene symbols of one DEV disease)
- query_uses: HGNC gene symbols of the disease only, TITLE_ABS phrase, OR-pooled; year <= 2024; abstract present; NO disease term, NO drug term, NO language filter
- record_cap: 300
- share_of_retrieved_papers_naming_any_candidate_drug: `{"max": 0.7966666666666666, "mean": 0.18963184548445772, "median": 0.14, "min": 0.03, "n": 103, "p25": 0.095, "p75": 0.21778195488721802}`
- share_of_retrieved_papers_with_a_gene_symbol_in_text: `{"max": 1.0, "mean": 0.8682737742969255, "median": 0.9, "min": 0.2633333333333333, "n": 103, "p25": 0.82, "p75": 0.9616666666666667}`
- text_used: title + '. ' + abstract (no full text)

**sentence_selection**

- caps: `{"papers_per_pair": 3, "sentences_per_paper": 2}`
- drug_rule: elc.drug_regex over all 2705 universe names: exact whole-word, case-insensitive match of the single project name (ChEMBL preferred name); names under 4 characters are never matched; longest alternative wins; no synonyms, brand names, salts, plurals or abbreviations; title+abstract sentences only
- function: engine_v4_llm.paper_drug_sentences / pair_papers / union_pair_papers
- name_gate: sentence must contain the first word of some drug name (speed only)
- pair_universe: top-30 B0 drugs union gold of the disease (candidate set P_d)
- paper_gate: paper must contain a disease gene symbol (v1 gate)
- relevance_rule_for_features: both annotators yes; direction 2-of-3 of Opus, Sonnet, rule
- verification: sentence must be a verbatim substring of the Europe PMC title+abstract after one join-period repair; 67 DEV items dropped

Manifest (DEV): `{"n_candidate_pairs": 4442, "n_items": 7637, "n_items_before_verification": 7704, "n_pairs_with_sentences": 2427, "n_papers": 3383}`

### Relaxation experiment on 150 seeded gold-positive pairs (no labelling)

R0 reproduces the stored pair index: paper list identical in 150/150, sentence count identical in 150/150. Pairs with a project synonym beyond the INN: 84 (terms by source: `{"chembl_synonyms": 89, "ctgov_other_names": 639, "hyphen_space": 17, "salt_stripped": 4}`).

| variant | sentences naming drug (mean) | drug+disease sentence, name only (mean) | drug+disease sentence, name+synonyms (mean) | pairs with >=1 | papers drug+disease (mean) | evaluative (rule) drug+disease (mean) | ratio vs R0 | +points pairs>=1 | major |
|---|---|---|---|---|---|---|---|---|---|
| R0_baseline | 0.85 | 0.11 | 0.11 | 11 (7.3%) | 0.35 | 0.03 | - | - | - |
| R1_sentence_cap_lifted | 1.01 | 0.11 | 0.11 | 11 (7.3%) | 0.35 | 0.03 | 1.00 | 0.0 | False |
| R2_paper_cap_lifted | 1.93 | 0.18 | 0.25 | 16 (10.7%) | 0.78 | 0.06 | 2.18 | 3.3 | True |
| R3_both_caps_lifted | 2.54 | 0.18 | 0.25 | 16 (10.7%) | 0.78 | 0.06 | 2.18 | 3.3 | True |
| R4_R3_no_gene_gate | 2.75 | 0.18 | 0.25 | 16 (10.7%) | 0.79 | 0.06 | 2.18 | 3.3 | True |
| R5_R3_drug_synonyms | 3.28 | 0.18 | 0.25 | 16 (10.7%) | 0.89 | 0.06 | 2.18 | 3.3 | True |
| R6a_direct_query_top10_INN | 15.39 | 2.59 | 4.09 | 77 (51.3%) | 5.14 | 0.89 | 36.06 | 44.0 | True |
| R6b_direct_query_top10_synonyms | 16.03 | 2.46 | 3.91 | 80 (53.3%) | 5.33 | 0.83 | 34.47 | 46.0 | True |
| R7_union_R4_R5_R6b | 18.96 | 2.57 | 4.08 | 80 (53.3%) | 5.97 | 0.88 | 36.00 | 46.0 | True |

At f = 0.03 (mean of 5 replicates, same pairs):

| variant | sentences naming drug | drug+disease (name+synonyms) | pairs with >=1 |
|---|---|---|---|
| R0_baseline | 0.05 | 0.01 | 3 (2.0%) |
| R1_sentence_cap_lifted | 0.07 | 0.01 | 3 (2.0%) |
| R2_paper_cap_lifted | 0.05 | 0.01 | 3 (2.0%) |
| R3_both_caps_lifted | 0.07 | 0.01 | 3 (2.0%) |
| R4_R3_no_gene_gate | 0.07 | 0.01 | 3 (2.0%) |
| R5_R3_drug_synonyms | 0.08 | 0.01 | 3 (2.0%) |

### Recall of the current retrieval against a direct `drug AND disease` query

Pairs 150; pairs with zero direct hits 57; direct hits per pair `{"max": 4177.0, "mean": 199.14, "median": 6.5, "min": 0.0, "n": 150, "p25": 0.0, "p75": 138.25}`.

Recall of the first 20 direct records among the 300 retrieved (INN query): `{"max": 1.0, "mean": 0.0376524803469775, "median": 0.0, "min": 0.0, "n": 93, "p25": 0.0, "p75": 0.05}`; pairs with recall 0: 67, recall 1: 1; synonym query: `{"max": 1.0, "mean": 0.04173621338569792, "median": 0.0, "min": 0.0, "n": 97, "p25": 0.0, "p75": 0.05}`.

Share of direct hits that also match the disease gene clause (what any gene-based retrieval can reach): `{"max": 1.0, "mean": 0.15970051220631737, "median": 0.015384615384615385, "min": 0.0, "n": 93, "p25": 0.0, "p75": 0.21052631578947367}`; total hits 29871, with gene clause 3352.

Gene query cap: `{"diseases": 103, "diseases_with_more_hits_than_300": 89, "gene_query_hits": {"max": 241472.0, "mean": 18771.019417475727, "median": 3464.0, "min": 46.0, "n": 103, "p25": 531.5, "p75": 14101.0}, "share_of_gene_query_hits_retrieved": {"max": 1.0, "mean": 0.2914906706787906, "median": 0.08660508083140878, "min": 0.0012423800689106812, "n": 103, "p25": 0.021372502684597623, "p75": 0.5645007116001444}}`

### Language

Retrieved records (30 seeded diseases, language field of the core record): `{"counts": {"chi": 76, "dut": 1, "eng": 7370, "fre": 5, "ger": 6, "heb": 1, "hun": 1, "ita": 4, "jpn": 2, "pol": 1, "por": 1, "rus": 11, "spa": 7, "unknown": 937}, "n": 8423, "non_english": 1053, "non_english_share": 0.12501484031817642}`; by disease `{"max": 0.23, "mean": 0.12115099387007751, "median": 0.11, "min": 0.04975124378109453, "n": 30, "p25": 0.08416666666666667, "p75": 0.13833333333333334}`; share of gene-query hits not English `{"max": 0.11748013620885357, "mean": 0.0766014000541677, "median": 0.07502222047770035, "min": 0.04975124378109453, "n": 30, "p25": 0.06347857860601078, "p75": 0.08763328026245959}`.

Direct query: `{"direct_query_hits_english": 27513, "direct_query_hits_not_english": 2358, "direct_query_hits_total": 29871, "direct_query_not_english_share_of_hits": 0.07893943959023803, "first20_records_language": {"counts": {"chi": 11, "dan": 1, "dut": 1, "eng": 1335, "fre": 11, "ger": 4, "ita": 2, "jpn": 13, "por": 1, "rus": 2, "spa": 7, "unknown": 42}, "n": 1430, "non_english": 95, "non_english_share": 0.06643356643356643}}`

### Drug naming

`{"distinct_gold_positive_drugs": 202, "drugs_with_at_least_one_project_synonym": 104, "drugs_with_at_least_one_synonym_only_sentence": 16, "drugs_with_non_plain_project_name": 2, "gold_positive_pairs": 465, "reading": "lower bound: only synonyms present in the project tables (ClinicalTrials.gov other names, partial ChEMBL records, salt stripping)", "retrieved_300_sentences_naming_drug_by_project_name": 1228, "retrieved_300_sentences_naming_drug_only_by_synonym": 297, "share_drugs_with_non_plain_project_name": 0.009900990099009901, "share_drugs_with_synonym": 0.5148514851485149, "share_of_drug_sentences_only_by_synonym": 0.19475409836065574, "share_of_gold_drugs_with_synonym_only_sentence": 0.07920792079207921, "synonym_terms_by_source_over_gold_drugs": {"chembl_synonyms": 119, "ctgov_other_names": 678, "hyphen_space": 28, "salt_stripped": 2}}`

## Q3. Label audit

### (a) Reading of gold-positive pairs with no `yes` from either annotator

Pool: 21 of 159 gold-positive pairs with items. Reader: the assistant (audit opinion, not ground truth); the pairs are known to be gold-positive, which may bias the reading toward 'informative'. Class counts: `{"informative_label_too_strict": 2, "label_looks_right": 9, "uninformative": 18}`.

| # | drug | disease | sentence | O (rel/dir) | S (rel/dir) | class | reason |
|---|---|---|---|---|---|---|---|
| 1 | methylprednisolone | chondrocalcinosis | Consistently, in rat ONFH models induced by methylprednisolone (MPSL) and imiquimod (IMI), significant bone loss in the femoral head was observed, attributable to increased numbers of TRAP positive osteoclasts. | no/harm | no/harm | uninformative | drug only induces an osteonecrosis rat model |
| 2 | propranolol | tufted angioma (disease) | She had pericardial effusion at the neonatal age that required pericardiocentesis, propranolol, and corticosteroids. | unclear/unclear | unclear/unclear | label_looks_right | treatments listed for a patient but no condition named in the sentence |
| 3 | risedronic acid | juvenile Paget disease | In addition, risedronic acid inhibited osteosarcoma growth at both cellular and molecular levels. | no/benefit | no/benefit | uninformative | drug acts on osteosarcoma growth, not on the disease |
| 4 | doxorubicin | kidney Wilms tumor | This work aims to elucidate the mechanism of thyroid hormone receptor interactor 13 (TRIP13) affecting doxorubicin (DOX) resistance in colorectal cancer (CRC). | no/harm | no/harm | uninformative | doxorubicin resistance in colorectal cancer, other disease |
| 5 | donepezil | Alzheimer disease | APP/PS1 mice of approximately nine months old were randomly divided into the model group, the low, medium, and high(2.5, 5, and 10 g·kg~(-1)) groups of Huanglian Jiedu Decoction, and 0.75 mg·kg~(-1) donepezil hydrochloride group, and the C57BL/6J mice with the same age were taken as the normal group. | unclear/unclear | no/unclear | label_looks_right | donepezil is only the positive-control group of a herbal test |
| 6 | cannabidiol | Dravet syndrome | Here we show that genome-wide miRNA screening of hippocampal tissue from a rat epilepsy model, mice treated with the antiseizure medicine cannabidiol, and plasma from patients with treatment-resistant epilepsy, converge on a single target-miR-335-5p. | unclear/unclear | no/unclear | label_looks_right | drug is a treatment of the mice, but the sentence reports a miRNA screen |
| 7 | apomorphine | parkinsonian-pyramidal syndrome | The mPFF-injected A53T BAC-SNCA Tg mice showed a 38% reduction in tyrosine hydroxylase (TH)-positive neurons in the ipsilateral SNpc, apomorphine-induced rotational behavior, and motor dysfunction at 2 months post-injection. | no/unclear | no/unclear | uninformative | apomorphine used as a tool to induce rotation |
| 8 | atorvastatin | homozygous familial hypercholesterolemia | APOE∗3-Leiden.CETP mice orally treated with NYX-PCSK9i demonstrated a dose-dependent decrease in plasma total cholesterol of up to 57%, while its combination with atorvastatin additively suppressed plasma total cholesterol levels. | unclear/benefit | unclear/benefit | label_looks_right | cholesterol-lowering combination in a mouse, disease not stated |
| 9 | clonazepam | Lennox-Gastaut syndrome | In 2018 she experienced worsening of tremor and developed painful dystonic attacks, resistant to treatments including clonazepam, trihexyphenidyl, baclofen, and levodopa/benserazide. | no/no_effect | no/no_effect | uninformative | patient with another movement disorder resistant to clonazepam |
| 10 | levothyroxine | familial thyroid carcinoma | Levothyroxine (LT4) serves as the standard hypothyroidism treatment, with required dosages influenced by the severity of the individual's disorder, which varies widely among affected individuals. | no/unclear | no/unclear | uninformative | levothyroxine for hypothyroidism, a different disease |
| 11 | vandetanib | thyroid gland carcinoma | Recent clinical trials for RET fusion-positive NSCLC using vandetanib or cabozantinib demonstrated positive clinical response and considerable differential activities for RET inhibitors among fusion partners. | no/benefit | no/benefit | uninformative | RET-fusion lung cancer trials, a different disease |
| 12 | daunorubicin | inherited acute myeloid leukemia | Western blot was used to determine the changes of the expression of CEBPA protein in K562 CEBPA p30, U937 and MOLM-13 cell lines before and after treatments with HHT, daunorubicin (DNR) or cytarabine (Ara-C). | no/unclear | no/unclear | label_looks_right | drug is a perturbation in a Western-blot assay on cell lines |
| 13 | dexamethasone | inherited isolated adrenal insufficiency due to partial CYP11A1 deficiency | Physiological roles of propolis and red ginseng nanoplatforms in alleviating dexamethasone-induced male reproductive challenges in a rat model. | no/harm | no/harm | uninformative | rat reproductive toxicity model |
| 14 | mitoxantrone | inherited acute myeloid leukemia | Interventions Complete remission of APL was achieved after all-trans retinoic acid and arsenic trioxide double induction therapy, followed by 2 cycles of mitoxantrone and cytarabine, and 1 cycle of idarubicin and cytarabine. | unclear/unclear | unclear/benefit | informative_label_too_strict | APL (an AML subtype) treated with mitoxantrone with remission reported; drug evaluated as leukaemia therapy |
| 15 | doxorubicin | inherited acute myeloid leukemia | ACSL4 is essential for HCC cell survival in response to glucose starvation, rather than to hypoxia or to the combination of hypoxia with doxorubicin or cisplatin. | no/unclear | no/unclear | uninformative | liver cancer cells under glucose starvation |
| 16 | vitamin e | familial isolated deficiency of vitamin E | Ataxia with isolated vitamin E deficiency: a Japanese family carrying a novel mutation in the alpha-tocopherol transfer protein gene. | no/unclear | no/unclear | uninformative | title of a family report; vitamin E is the deficient nutrient, no treatment |
| 17 | tocilizumab | temporal arteritis | We recruited 106 patients with active RA who had started abatacept, tocilizumab, or TNF inhibitors as a first-line bDMARDs. | no/unclear | no/unclear | uninformative | rheumatoid arthritis cohort, a different disease |
| 18 | prednisolone | temporal arteritis | Prednisolone treatment was effective for SAT. | unclear/benefit | no/benefit | label_looks_right | treatment effect stated but the condition is the abbreviation SAT, not the item disease |
| 19 | clobazam | Lennox-Gastaut syndrome | Despite treating the seizures with multiple antiepileptic drugs, including phenytoin, sodium valproate, levetiracetam, perampanel, and clobazam, it was very difficult to control the seizures, and genetic testing was suggested. | unclear/no_effect | no/no_effect | informative_label_too_strict | clobazam named among drugs that failed to control the seizures of a patient; a null treatment result in an epilepsy; Sonnet no looks strict |
| 20 | topiramate | Lennox-Gastaut syndrome | Treatment includes antiepileptics such as lamotrigine, valproate or topiramate, calcium blockers such as flunarizine or verapamil and acetazolamide. | unclear/unclear | unclear/unclear | label_looks_right | drug list for an unspecified condition |
| 21 | dexamethasone | X-linked adrenal hypoplasia congenita | Dexamethasone downregulates expression of several genes encoding orphan nuclear receptors that are important to steroidogenesis in stallion testes. | no/unclear | no/unclear | uninformative | stallion testes gene expression |
| 22 | vitamin e | familial isolated deficiency of vitamin E | [Clinical and genetic analysis of a patient with Ataxia and vitamin E deficiency due to homozygous variant of TTPA gene]. | no/unclear | no/unclear | uninformative | title of a case report; vitamin E not given as treatment |
| 23 | vitamin e | familial isolated deficiency of vitamin E | The tocopherol transfer protein mediates vitamin E trafficking between cerebellar astrocytes and neurons. | no/unclear | no/unclear | uninformative | basic biology of vitamin E trafficking |
| 24 | dexamethasone | X-linked adrenal hypoplasia congenita | Stallions treated with the synthetic glucocorticoid dexamethasone had testes collected at 6 or 12 hours postinjection. | no/unclear | no/unclear | uninformative | stallion testes sampling |
| 25 | doxorubicin | kidney Wilms tumor | TRIP13 Activates Glycolysis to Promote Cell Stemness and Strengthen Doxorubicin Resistance of Colorectal Cancer Cells. | no/harm | no/harm | uninformative | colorectal cancer cell stemness |
| 26 | dexamethasone | inherited isolated adrenal insufficiency due to partial CYP11A1 deficiency | Objective This study evaluated the efficiency of red ginseng and propolis, either in powder or as nano-forms against dexamethasone-induced testicular oxidative challenges in adult male albino rats. | no/harm | no/harm | uninformative | dexamethasone as the toxicant in rats |
| 27 | dexamethasone | inherited isolated adrenal insufficiency due to partial CYP11A1 deficiency | We also demonstrated that dexamethasone synchronization in KGN cells enhanced the rhythmic alterations in circadian clock genes. | no/unclear | no/unclear | uninformative | cell-clock experiment, dexamethasone as synchroniser |
| 28 | atorvastatin | homozygous familial hypercholesterolemia | The efficacies of the hypolipidemic effects of allicin (200 µM) and capsaicin (200 µM) were comparable to that of atorvastatin (10 µM) in this study. | no/unclear | no/unclear | label_looks_right | atorvastatin is the comparator in an assay of other compounds |
| 29 | atorvastatin | homozygous familial hypercholesterolemia | Notably, piperines demonstrated synergistic beneficial effects with atorvastatin or ezetimibe, which are two widely used hypocholesterolemic and anti-atherosclerotic drugs. | unclear/benefit | unclear/benefit | label_looks_right | synergy with a statin in a cholesterol model, familial disease not tied |

### (b) Opus x Sonnet cross-tabs (rows Opus, columns Sonnet)

**gold_positive_pairs** (n = 803): relevance agreement 55.3%, kappa 0.350; direction agreement 91.8%, kappa 0.855; share yes `{"opus": 0.3300124533001245, "sonnet": 0.7272727272727273}`, either 73.2%, both 32.5%

- relevance: `{"no": {"no": 147, "unclear": 10, "yes": 12}, "unclear": {"no": 22, "unclear": 36, "yes": 311}, "yes": {"no": 2, "unclear": 2, "yes": 261}}`
- direction: `{"benefit": {"benefit": 254, "harm": 0, "no_effect": 0, "unclear": 5}, "harm": {"benefit": 0, "harm": 47, "no_effect": 1, "unclear": 10}, "no_effect": {"benefit": 0, "harm": 1, "no_effect": 23, "unclear": 2}, "unclear": {"benefit": 44, "harm": 0, "no_effect": 3, "unclear": 413}}`

**other_pairs** (n = 6834): relevance agreement 83.1%, kappa 0.552; direction agreement 93.6%, kappa 0.835; share yes `{"opus": 0.04184957565115598, "sonnet": 0.16637401229148377}`, either 16.8%, both 4.1%

- relevance: `{"no": {"no": 5126, "unclear": 102, "yes": 83}, "unclear": {"no": 186, "unclear": 274, "yes": 777}, "yes": {"no": 6, "unclear": 3, "yes": 277}}`
- direction: `{"benefit": {"benefit": 1028, "harm": 1, "no_effect": 2, "unclear": 60}, "harm": {"benefit": 0, "harm": 304, "no_effect": 2, "unclear": 38}, "no_effect": {"benefit": 0, "harm": 3, "no_effect": 86, "unclear": 17}, "unclear": {"benefit": 237, "harm": 68, "no_effect": 12, "unclear": 4976}}`

### (c) Coverage under different relevance rules (pairs, mean over replicates)

f = 0.03: `{"gold_positive": {"S3_verified_sentence": 16.8, "S4_relevance_yes_either": 13.8, "S5_relevance_yes_both": 8.4, "S6_yes_both_and_benefit": 3.4, "S6b_yes_either_and_benefit": 7.0, "benefit_direction_only": 8.6, "opus_yes": 8.4, "opus_yes_benefit": 3.4, "sonnet_yes": 13.8, "sonnet_yes_benefit": 7.0}, "other": {"S3_verified_sentence": 162.2, "S4_relevance_yes_either": 36.4, "S5_relevance_yes_both": 10.0, "S6_yes_both_and_benefit": 4.2, "S6b_yes_either_and_benefit": 16.8, "benefit_direction_only": 36.6, "opus_yes": 10.2, "opus_yes_benefit": 4.2, "sonnet_yes": 36.2, "sonnet_yes_benefit": 16.8}}`

f = 0.1: `{"gold_positive": {"S3_verified_sentence": 38.8, "S4_relevance_yes_either": 32.2, "S5_relevance_yes_both": 17.0, "S6_yes_both_and_benefit": 8.2, "S6b_yes_either_and_benefit": 16.0, "benefit_direction_only": 18.0, "opus_yes": 17.2, "opus_yes_benefit": 8.2, "sonnet_yes": 32.0, "sonnet_yes_benefit": 16.0}, "other": {"S3_verified_sentence": 404.2, "S4_relevance_yes_either": 82.8, "S5_relevance_yes_both": 20.4, "S6_yes_both_and_benefit": 7.8, "S6b_yes_either_and_benefit": 41.8, "benefit_direction_only": 90.0, "opus_yes": 21.0, "opus_yes_benefit": 8.2, "sonnet_yes": 82.2, "sonnet_yes_benefit": 41.4}}`

f = 1.0: `{"gold_positive": {"S3_verified_sentence": 159.0, "S4_relevance_yes_either": 138.0, "S5_relevance_yes_both": 77.0, "S6_yes_both_and_benefit": 37.0, "S6b_yes_either_and_benefit": 74.0, "benefit_direction_only": 82.0, "opus_yes": 78.0, "opus_yes_benefit": 37.0, "sonnet_yes": 137.0, "sonnet_yes_benefit": 74.0}, "other": {"S3_verified_sentence": 2258.0, "S4_relevance_yes_either": 545.0, "S5_relevance_yes_both": 156.0, "S6_yes_both_and_benefit": 78.0, "S6b_yes_either_and_benefit": 282.0, "benefit_direction_only": 644.0, "opus_yes": 162.0, "opus_yes_benefit": 79.0, "sonnet_yes": 540.0, "sonnet_yes_benefit": 281.0}}`

### (d) v2 pilot on the same 800 items

Disagreement on relevance: `{"kept_in_v2": 23, "newly_created_by_v2": 28, "removed_by_v2": 102, "share_of_v1_disagreement_removed_net": 0.5920000000000001, "v1_items": 125, "v2_items": 51}`; Opus vs Sonnet relevance kappa v1 0.625 -> v2 0.851 (agreement 84.4% -> 93.6%); direction kappa 0.873 -> 0.866; yes counts `{"v1": {"opus_yes": 59, "sonnet_yes": 151, "yes_both": 57, "yes_either": 153}, "v2": {"opus_yes": 45, "sonnet_yes": 66, "yes_both": 41, "yes_either": 70}}`.

- opus relevance v1 (rows) vs v2 (columns): `{"no": {"no": 578, "unclear": 16, "yes": 2}, "unclear": {"no": 4, "unclear": 139, "yes": 2}, "yes": {"no": 2, "unclear": 16, "yes": 41}}`; direction `{"benefit": {"benefit": 130, "harm": 0, "no_effect": 0, "unclear": 3}, "harm": {"benefit": 0, "harm": 49, "no_effect": 0, "unclear": 7}, "no_effect": {"benefit": 0, "harm": 0, "no_effect": 17, "unclear": 0}, "unclear": {"benefit": 5, "harm": 3, "no_effect": 0, "unclear": 586}}`
- sonnet relevance v1 (rows) vs v2 (columns): `{"no": {"no": 567, "unclear": 29, "yes": 4}, "unclear": {"no": 5, "unclear": 40, "yes": 4}, "yes": {"no": 3, "unclear": 90, "yes": 58}}`; direction `{"benefit": {"benefit": 137, "harm": 0, "no_effect": 0, "unclear": 7}, "harm": {"benefit": 0, "harm": 48, "no_effect": 0, "unclear": 10}, "no_effect": {"benefit": 0, "harm": 1, "no_effect": 14, "unclear": 0}, "unclear": {"benefit": 10, "harm": 4, "no_effect": 2, "unclear": 567}}`
- Opus x Sonnet relevance v1 `{"no": {"no": 579, "unclear": 9, "yes": 8}, "unclear": {"no": 20, "unclear": 39, "yes": 86}, "yes": {"no": 1, "unclear": 1, "yes": 57}}` v2 `{"no": {"no": 566, "unclear": 13, "yes": 5}, "unclear": {"no": 9, "unclear": 142, "yes": 20}, "yes": {"no": 0, "unclear": 4, "yes": 41}}`

## Q4. Verdict

Computed summary: `{"gold_positive_at_f1.0": {"largest_absolute_loss": {"from": "S0_candidate_pairs", "lost_pairs": 302.0, "retention": 0.35053763440860214, "to": "S1a_drug_named_in_a_kept_paper"}, "lowest_retention": {"from": "S0_candidate_pairs", "lost_pairs": 302.0, "retention": 0.35053763440860214, "to": "S1a_drug_named_in_a_kept_paper"}}, "gold_positive_at_f0.03": {"largest_absolute_loss": {"from": "S0_candidate_pairs", "lost_pairs": 447.6, "retention": 0.03741935483870967, "to": "S1a_drug_named_in_a_kept_paper"}, "lowest_retention": {"from": "S0_candidate_pairs", "lost_pairs": 447.6, "retention": 0.03741935483870967, "to": "S1a_drug_named_in_a_kept_paper"}}, "gold_positive_share_of_S0": {"1.0": {"S1a_drug_named_in_a_kept_paper": 0.35053763440860214, "S1b_and_paper_has_disease_gene": 0.3419354838709677, "S2_sentence_extracted": 0.3419354838709677, "S3_verified_sentence": 0.3419354838709677, "S4_relevance_yes_either": 0.2967741935483871, "S5_relevance_yes_both": 0.16559139784946236, "S6_yes_both_and_benefit": 0.07956989247311828, "S6b_yes_either_and_benefit": 0.15913978494623657}, "0.1": {"S1a_drug_named_in_a_kept_paper": 0.08645161290322581, "S1b_and_paper_has_disease_gene": 0.08344086021505376, "S2_sentence_extracted": 0.08344086021505376, "S3_verified_sentence": 0.08344086021505376, "S4_relevance_yes_either": 0.069247311827957, "S5_relevance_yes_both": 0.03655913978494624, "S6_yes_both_and_benefit": 0.017634408602150535, "S6b_yes_either_and_benefit": 0.034408602150537634}, "0.03": {"S1a_drug_named_in_a_kept_paper": 0.03741935483870967, "S1b_and_paper_has_disease_gene": 0.03612903225806452, "S2_sentence_extracted": 0.03612903225806452, "S3_verified_sentence": 0.03612903225806452, "S4_relevance_yes_either": 0.02967741935483871, "S5_relevance_yes_both": 0.01806451612903226, "S6_yes_both_and_benefit": 0.007311827956989247, "S6b_yes_either_and_benefit": 0.015053763440860216}}, "variants_meeting_major_artefact_rule": ["R2_paper_cap_lifted", "R3_both_caps_lifted", "R4_R3_no_gene_gate", "R5_R3_drug_synonyms", "R6a_direct_query_top10_INN", "R6b_direct_query_top10_synonyms", "R7_union_R4_R5_R6b"], "relaxation_pairs_with_ge1_drug_disease_sentence": {"R0_baseline": 11, "R1_sentence_cap_lifted": 11, "R2_paper_cap_lifted": 16, "R3_both_caps_lifted": 16, "R4_R3_no_gene_gate": 16, "R5_R3_drug_synonyms": 16, "R6a_direct_query_top10_INN": 77, "R6b_direct_query_top10_synonyms": 80, "R7_union_R4_R5_R6b": 80}, "relaxation_mean_drug_disease_sentences": {"R0_baseline": 0.11333333333333333, "R1_sentence_cap_lifted": 0.11333333333333333, "R2_paper_cap_lifted": 0.24666666666666667, "R3_both_caps_lifted": 0.24666666666666667, "R4_R3_no_gene_gate": 0.24666666666666667, "R5_R3_drug_synonyms": 0.24666666666666667, "R6a_direct_query_top10_INN": 4.086666666666667, "R6b_direct_query_top10_synonyms": 3.9066666666666667, "R7_union_R4_R5_R6b": 4.08}, "recall_mean_first20": 0.0376524803469775, "pairs_zero_direct_hits": 57, "direct_hits_also_matching_gene_clause_share_of_total": 0.1122158615379465, "language": {"explicit_non_english_share_retrieved": 0.01377181526771934, "unknown_language_share_retrieved": 0.11124302505045708, "explicit_non_english_share_first20_direct": 0.03706293706293706, "material_at_5_percent": false, "note": "records with no language field are all preprints (source PPR) in the retrieved set and are not counted as non-English"}, "label_strictness": {"either_over_both_pairs_gold_positive_f1.0": 1.7922077922077921, "either_over_both_pairs_gold_positive_f0.03": 1.6428571428571428, "reading_too_strict_share": 0.06896551724137931, "n_read": 29}}`

**Short answer.** The claim is right about the outcome and wrong about the cause. With the current construction the LLM features are almost empty at 3 % of the literature (gold-positive pairs with a counted relevant benefit sentence: 3.4 of 465, 0.7 %), but they are also nearly empty at 100 % of the literature (37 of 465, 8.0 %; 74 of 465, 15.9 % if either annotator may say yes). So the problem is not 'low volume'; thinning multiplies an already poor base. The base is poor mainly because of how papers are retrieved (ii), partly because some gold pairs really have little drug-and-disease literature (i), and only marginally because of label strictness (iii). Language and drug naming are not material losses.

**Where the pairs are lost (Q1, gold-positive, f = 1.0).** 465 candidate pairs -> 163 (35.1 %) have the drug named in any of the retrieved papers -> 159 have an extracted and verified sentence -> 138 have `yes` from either annotator -> 77 from both -> 37 with benefit direction. The largest single loss is the first step (302 pairs, retention 35 %); the extraction, the paper gate and the verification lose 4 pairs in total. At f = 0.03 the first step keeps 3.7 % of the pairs (17.4 of 465). Non-gold candidates behave the same way at the first step (57 % keep a mention at f = 1.0) but lose most pairs later, at the relevance step (545 of 2258 with a verified sentence have `yes` by either).

**(ii) Retrieval is the artefact.** The paper query pools the gene symbols of the disease only, with no disease term, no drug term and no language filter, and keeps the first 300 records with an abstract; 89 of 103 DEV diseases are at that cap and the 300 are a median 8.7 % of the gene-query hits. A direct `drug AND disease` query shows how far this is from the drug-and-disease literature: of the first 20 direct records only 3.8 % (mean over the 93 pairs that have hits) are among the 300 retrieved, and only 11 % of all direct hits (3352 of 29871) even contain a disease gene symbol, so 89 % cannot be reached by any gene-symbol query whatever the cap. On the 150 seeded gold-positive pairs the baseline gives 0.11 sentences per pair naming drug and disease (11 pairs, 7 %); the top-10 direct-query papers give 3.9 (80 pairs, 53 %), and 0.83 vs 0.03 sentences per pair that the rule extractor reads as evaluative (benefit or harm). Rules relaxed one at a time inside the existing retrieval change little: lifting the sentence cap adds nothing to drug-and-disease sentences; lifting the paper cap gives 0.25 vs 0.11 (16 vs 11 pairs; flagged 'major' by the doubling convention only because the base is tiny); removing the gene gate adds nothing; drug synonyms add nothing to drug-and-disease sentences (0.25).

**(i) Part of it is the literature.** 57 of the 150 gold-positive pairs (38 %) have zero direct hits for the drug AND the disease name or up to five of its synonyms in title or abstract (median 6.5 hits per pair). Some of these are rare diseases whose MONDO labels are not the words used in papers, so this is an upper bound on true absence, but for them no retrieval design will provide sentences; for those pairs the claim 'not testable' is correct. The 93 pairs with hits have a large literature (mean 199 hits over all 150 pairs, maximum 4177).

**(iii) Label strictness is minor for coverage.** On DEV gold-positive items Sonnet says `yes` for 72.7 % and Opus for 33.0 % (kappa 0.35, versus 0.55 for other pairs); `either` vs `both` changes the gold-positive pairs with a counted `yes` sentence from 77 to 138 at f = 1.0 (x1.8) and from 8.4 to 13.8 at f = 0.03 (x1.6), below the doubling convention. In my reading of the gold-positive pairs where neither annotator said `yes` (only 21 pairs; 29 sentences read, audit opinion) 18 sentences were uninformative, 9 labels looked right and 2 looked too strict (7 %), below the quarter threshold. The disagreement is mostly Opus `unclear` against Sonnet `yes` (311 items on gold-positive pairs). The v2 guideline removes 102 of 125 v1 disagreements on the 800 pilot items (kappa 0.62 -> 0.85) but creates 28 and does so by making Sonnet stricter (`yes` by either 153 -> 70, by both 57 -> 41), so it lowers coverage further. The much higher `yes` rate on gold-positive than on other pairs (either 73 % vs 17 %) is compatible with real signal in the sentences and with annotator knowledge of indications (parent S4); this audit cannot separate the two.

**Language and naming.** Explicit non-English records are 1.4 % of the 8423 retrieved records of 30 seeded diseases and 3.7 % of the first 20 direct records; 11 % of retrieved records have no language field and all are preprints, so they are not counted as non-English. Below the 5 % convention, not material. Drug naming: 2 of 202 gold-positive drugs have a non-plain project name; 104 have a project synonym (mostly ClinicalTrials.gov other names, so a lower bound); 16 (7.9 %) have sentences naming the drug only by a synonym, 19.5 % of the drug sentences of pairs with synonyms, but they add nothing to drug-and-disease sentences, so naming is not the bottleneck either.

**Recommendation (not implemented, nothing annotated).** The cheapest change with the largest gain is to retrieve per candidate pair instead of per disease: one Europe PMC query `drug (INN and synonyms) AND disease (name and a few synonyms)` for the first 10 to 20 records, and keep the sentences that name the drug (prefer those that also name the disease). Estimated gain from the 150-pair sample: pairs with at least one drug-and-disease sentence from 7 % to 53 %, sentences per pair from 0.11 to 3.9, i.e. about 7 times the pairs and 35 times the sentences. Cost: about one call per candidate pair (4442 DEV pairs, at most 4 requests per second), and about 4 more sentences per pair to annotate (about 17 thousand DEV items if the other pairs behave like gold-positive ones; fewer in practice). Low volume can then be simulated by keeping a fraction of each pair's own hits, which keeps the test meaningful for pairs with many hits; pairs with none stay empty, which is the true literature limit. Estimates come from gold-positive pairs only and count sentences, not relevant sentences; a second, much smaller lever is accepting `yes` by either annotator (x1.6 to x1.8 pairs), which does not address the 65 % of gold pairs never retrieved.

