# Worked case: kidney Wilms tumor

**worked case on a v2 disease, not part of the pre-registered primary validation.**

> Hypotheses, not advice. No dosing and no treatment recommendation. Every drug below that is an oncology drug (ATC L01) carries serious toxicity (chemotherapy caution): bone-marrow suppression, organ toxicity, second malignancy and embryo-fetal harm are label-documented hazards, not features.

Seed 20261015. Inputs: benchmark_v2 sha256 `42ed7ff06e7b5ca3`, universe sha256 `c4df43ba87d3146b` (ChEMBL_37). Genes (14): BRCA2, CTR9, DICER1, DIS3L2, GPC3, GPC4, H19, IGF2, PAX6, POU6F2, REST, TRIM28, TRIP13, WT1. Known drugs: gold 1 = dactinomycin, doxorubicin, vincristine; gold 3 = dactinomycin, doxorubicin, vincristine; corrected gold 2 has 21 drugs.

## 1. Engine ranking without the disease's own gold

Headline: arm A (prior + popularity correction), refit on the 44 other benchmark-v2 diseases (strict: no label of the target is read), literature <= 2024, applied to the target for all 2705 universe drugs. Check: protocol grouped-CV fold in which the disease is held out (same folds as validation_v2.json) (the fold of 9 diseases that contains the target).
Reproduction of validation_v2.json: fold AUC vs gold 1 = 0.441401 (stored 0.441401; identical: True; folds identical: True). Strict refit AUC vs gold 1 = 0.785097.

Top 20 of all universe drugs (score = logistic decision value; probability = Platt-calibrated):

| # | drug | score | probability | FDA label | type | gold 1 | gold 2 | gold 3 | papers with gene | benefit | human |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | tamoxifen | -1.438516 | 0.145539 | yes | Small molecule |  |  |  | 2 | 1 | 2 |
| 2 | progesterone | -2.24242 | 0.078175 | yes | Small molecule |  |  |  | 2 | 1 | 1 |
| 3 | danazol | -2.931904 | 0.044552 | yes | Small molecule |  |  |  | 1 | 0 | 0 |
| 4 | lenvatinib | -3.056976 | 0.040155 | yes | Small molecule |  |  |  | 1 | 1 | 0 |
| 5 | sevoflurane | -3.079816 | 0.039398 | yes | Small molecule |  |  |  | 1 | 1 | 0 |
| 6 | azacitidine | -3.232289 | 0.034686 | yes | Small molecule |  |  |  | 1 | 0 | 0 |
| 7 | cytarabine | -3.567516 | 0.026162 | yes | Small molecule |  |  |  | 1 | 1 | 1 |
| 8 | testosterone | -3.652527 | 0.024348 | yes | Small molecule |  |  |  | 1 | 0 | 0 |
| 9 | carboplatin | -3.855952 | 0.020489 | yes | Small molecule |  | x |  | 1 | 0 | 0 |
| 10 | rituximab | -3.918948 | 0.019421 | yes | Antibody |  |  |  | 1 | 1 | 1 |
| 11 | zinc chloride | -3.946957 | 0.018964 | yes | Small molecule |  |  |  | 0 | 0 | 0 |
| 12 | cisplatin | -4.19579 | 0.015338 | yes | Small molecule |  |  |  | 1 | 0 | 0 |
| 13 | palbociclib | -4.206648 | 0.015197 | yes | Small molecule |  |  |  | 1 | 0 | 0 |
| 14 | paclitaxel | -4.279512 | 0.014279 | yes | Small molecule |  |  |  | 2 | 1 | 2 |
| 15 | heparin | -4.303735 | 0.013986 | yes | Oligosaccharide |  |  |  | 1 | 0 | 0 |
| 16 | topiramate | -4.394721 | 0.012939 | yes | Small molecule |  |  |  | 0 | 0 | 0 |
| 17 | fostamatinib | -4.396717 | 0.012916 | yes | Small molecule |  |  |  | 0 | 0 | 0 |
| 18 | diethylstilbestrol | -4.402599 | 0.012852 | yes | Small molecule |  |  |  | 0 | 0 | 0 |
| 19 | gemcitabine | -4.448094 | 0.01236 | yes | Small molecule |  |  |  | 2 | 0 | 1 |
| 20 | cannabidiol | -4.558465 | 0.011245 | yes | Small molecule |  |  |  | 0 | 0 | 0 |

## 2. Filter chain and label-based age support

Universe 2705 -> FDA approved 1646 (label rule alone: 1641; added by the Drugs@FDA check of 13 drugs: arsenic trioxide, citric acid, lenvatinib, lurbinectedin, selpercatinib; not approved by either rule: 1059) -> small molecule 1296. Dropped by type after the FDA step: 286 biologics, 11 oligonucleotides, 53 of unknown type.

- Known drugs of gold1 dropped by the filters: dactinomycin (label: True, ChEMBL type Protein).
- Known drugs of gold3 dropped by the filters: dactinomycin (label: True, ChEMBL type Protein).
- Known drugs of gold2 dropped by the filters: atezolizumab (label: True, ChEMBL type Antibody); dactinomycin (label: True, ChEMBL type Protein); ipilimumab (label: True, ChEMBL type Antibody); nivolumab (label: True, ChEMBL type Antibody); oxytetracycline (label: False, ChEMBL type Small molecule); trastuzumab deruxtecan (label: True, ChEMBL type Antibody drug conjugate).

Top 20 after the filters (FDA label, small molecule):

| # | drug | score | Drugs@FDA | 0-4 y | 5-11 y | 12-17 y | basis |
|---|---|---|---|---|---|---|---|
| 1 | tamoxifen | -1.438516 | True | not_labelled | not_labelled | not_labelled | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 2 | progesterone | -2.24242 | True | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 3 | danazol | -2.931904 | True | not_labelled | not_labelled | not_labelled | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 4 | lenvatinib | -3.056976 | True | not_labelled | not_labelled | not_labelled | manual review (config/pediatric_label_review.json) |
| 5 | sevoflurane | -3.079816 | True | partly (from 0.166667 y) | labelled | labelled | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 6 | azacitidine | -3.232289 | True | not_labelled | not_labelled | not_labelled | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 7 | cytarabine | -3.567516 | True | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 8 | testosterone | -3.652527 | True | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 9 | carboplatin | -3.855952 | True | not_labelled | not_labelled | not_labelled | manual review (config/pediatric_label_review.json) |
| 10 | zinc chloride | -3.946957 | True | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 11 | cisplatin | -4.19579 | True | not_labelled | not_labelled | not_labelled | manual review (config/pediatric_label_review.json) |
| 12 | palbociclib | -4.206648 | True | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 13 | paclitaxel | -4.279512 | True | not_labelled | not_labelled | not_labelled | manual review (config/pediatric_label_review.json) |
| 14 | topiramate | -4.394721 | True | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 15 | fostamatinib | -4.396717 | True | not_labelled | not_labelled | not_labelled | manual review (config/pediatric_label_review.json) |
| 16 | diethylstilbestrol | -4.402599 | False | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 17 | gemcitabine | -4.448094 | True | not_labelled | not_labelled | not_labelled | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 18 | cannabidiol | -4.558465 | True | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 19 | cocaine | -4.567293 | True | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand- [...] |
| 20 | spironolactone | -4.572736 | True | not_labelled | not_labelled | not_labelled | heuristic classifier of labels.classify_pediatric (not hand- [...] |

Verbatim label sentences (section pediatric_use unless stated; not dosing information):

- **tamoxifen** (DailyMed, set id 7ee3d3d2-85d1-4018-8e70-5ed8a64ae1f0, section pediatric_use): "Pediatric Use The safety and efficacy of tamoxifen for girls aged 2 to 10 years with McCune-Albright Syndrome and precocious puberty have not been studied beyond one year of treatment."
- **progesterone** (DailyMed, set id d7619751-47b4-47c9-b716-e2a704b350c8, section pediatric_use): "Pediatric Use Progesterone capsules are not indicated in children."
- **danazol** (DailyMed, set id e19acee5-ff33-45a5-bbbf-801756bc59b4, section pediatric_use): "Pediatric Use Safety and effectiveness in pediatric patients have not been established."
- **lenvatinib** (DailyMed, set id f4bedd21-efde-44c6-9d9c-b48b78d7ed1e, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of LENVIMA in pediatric patients have not been established."
- **sevoflurane** (DailyMed, set id c808afe1-ebf1-44cf-9187-09656331bff5, section pediatric_use): "Pediatric Use Induction and maintenance of general anesthesia with sevoflurane have been established in controlled clinical studies in pediatric patients aged 1 to 18 years (see CLINICAL STUDIES, ADVERSE REACTIONS)."
- **azacitidine** (DailyMed, set id a6a5ca37-37e8-ab65-bc0b-9c0682ab52f7, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness of azacitidine in pediatric patients with MDS have not been established."
- **cytarabine** (DailyMed, set id 21e8e8e6-3f85-4fdd-8452-c9d14d2ef69a, section pediatric_use): "Pediatric use: See INDICATIONS AND USAGE."
- **testosterone** (DailyMed, set id a2c73a06-1200-4322-9011-caea5f16db6c, section pediatric_use): "Pediatric use Safety and effectiveness in pediatric patients below the age of 12 years have not been established."
- **carboplatin** (DailyMed, set id 4c89cedc-e48b-445c-ad05-54dbe60e4fb0, section pediatric_use): "Pediatric Use Safety and effectiveness in pediatric patients have not been established (seeWARNINGS : "audiologic toxicity")."
- **zinc chloride** (DailyMed, set id 2bd646d2-2a2a-4376-8da0-c7f08b0511ac, section pediatric_use): "Pediatric Use See DOSAGE AND ADMINISTRATION section."
- **cisplatin** (DailyMed, set id 508496cb-3441-46b3-a4fe-e0d440e6adc6, section pediatric_use): "8.4 Pediatric Use Ototoxic effects may be more severe and detrimental in pediatric patients receiving cisplatin for injection, particularly in patients less than 5 years of age."
- **palbociclib** (DailyMed, set id 89a142d5-7e61-48bf-9e77-e2cd124c4792, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of IBRANCE in pediatric patients have not been established."
- **paclitaxel** (DailyMed, set id 0c12c77c-5fd4-48dc-a318-e96892233992, section pediatric_use): "Pediatric Use The safety and effectiveness of paclitaxel in pediatric patients have not been established."
- **topiramate** (DailyMed, set id 8d5f9f77-403e-432a-8dcb-7995c6947ae0, section pediatric_use): "8.4 Pediatric Use Adjunctive Treatment for Epilepsy Pediatric Patients 2 Years of Age and Older The safety and effectiveness of topiramate extended-release capsules as adjunctive therapy for the treatment of partial-onset seizures, primary generalized tonic-clonic seizures, or seizures associated with Lennox-Gastaut syndrome have been established in pediatric patients 2 years of age and older and is based on controll [...]"
- **fostamatinib** (DailyMed, set id 21149cc3-049b-43e2-b141-c9499160556c, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness in pediatric patients have not been established."
- **diethylstilbestrol**: no readable paediatric-use sentence (heuristic classifier of labels.classify_pediatric (not hand-verified); no readable paediatric-use section); status UNKNOWN
- **gemcitabine** (DailyMed, set id 44bcc16c-341d-49e8-a1c6-95d084667aca, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of Gemcitabine for Injection have not been established in pediatric patients."
- **cannabidiol** (DailyMed, set id 8bf27097-4870-43fb-94f0-f3d0871d1eec, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness of EPIDIOLEX for the treatment of seizures associated with LGS, DS, or TSC have been established in patients 1 year of age and older."
- **cocaine**: no readable paediatric-use sentence (heuristic classifier of labels.classify_pediatric (not hand-verified); no readable paediatric-use section); status UNKNOWN
- **spironolactone** (DailyMed, set id 768563dc-a619-44dc-b97f-1a7e71fd56a0, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness in pediatric patients have not been established."

Known drugs (gold 1 / gold 3) with a current FDA label, label-based age support (any molecule type; not filtered):

| drug | 0-4 y | 5-11 y | 12-17 y | basis | verbatim paediatric-use sentence |
|---|---|---|---|---|---|
| dactinomycin | verify | verify | verify | heuristic classifier of labels.classify_pediatric (not hand-verified) | "8.4 Pediatric Use The safety and effectiveness of dactinomycin have been established in pediatric patients with Wilms tumor, rhabdomyosarcoma, Ewing sarcoma, and metastatic nonseminomatous testicular cancer." (DailyMed, set id 5da55538-2602-4704-a9da-603c9f7ff7dc, section pediatric_use) |
| doxorubicin | age_not_stated | age_not_stated | age_not_stated | manual review (config/pediatric_label_review.json) | "8.4 Pediatric Use Based on postmarketing reports, pediatric patients treated with doxorubicin hydrochloride are at risk for developing late cardiovascular dysfunction." (DailyMed, set id 1fd148fb-0fbc-4b6f-b790-23546fb46a71, section pediatric_use) |
| vincristine | verify | verify | verify | manual review verdict: verify | "Pediatric Use See DOSAGE AND ADMINISTRATION section." (DailyMed, set id b22e4c8e-dbe1-4d4b-a026-e1812974071c, section pediatric_use) |

- Manual-review scope for doxorubicin: "label read: conventional doxorubicin hydrochloride injection (the earlier automatic read used the liposomal product, where paediatric use is not established). The indications include childhood cancers (metastatic Wilms' tumor, metastatic neuroblastoma) and the label gives paediatric safety informati [...]"
- Manual-review scope for vincristine: "the paediatric section only refers to the dosage section; the label gives paediatric dosing (including for patients weighing 10 kg or less) and warns of deaths in children under 13 years after doses ten times the recommended dose; no minimum age and no statement of established effectiveness (implici [...]"

## 3. Concordance with the known drugs

### gold1

among approved (FDA label) candidates: 1646 candidates, 3 of 3 known drugs present, AUC 0.706431; certain hits at 10: 0, at 20: 0 (tie-averaged expected 0.0 / 0.0).

| drug | rank (min-max) | top % |
|---|---|---|
| doxorubicin | 89-89 | 5.407 |
| dactinomycin | 416-416 | 25.273 |
| vincristine | 948-948 | 57.594 |

among approved small molecules: 1296 candidates, 2 of 3 known drugs present, AUC 0.655719; certain hits at 10: 0, at 20: 0 (tie-averaged expected 0.0 / 0.0).

| drug | rank (min-max) | top % |
|---|---|---|
| doxorubicin | 79-79 | 6.096 |
| vincristine | 815-815 | 62.886 |

Not in the candidate set: dactinomycin

### gold3

among approved (FDA label) candidates: 1646 candidates, 3 of 3 known drugs present, AUC 0.706431; certain hits at 10: 0, at 20: 0 (tie-averaged expected 0.0 / 0.0).

| drug | rank (min-max) | top % |
|---|---|---|
| doxorubicin | 89-89 | 5.407 |
| dactinomycin | 416-416 | 25.273 |
| vincristine | 948-948 | 57.594 |

among approved small molecules: 1296 candidates, 2 of 3 known drugs present, AUC 0.655719; certain hits at 10: 0, at 20: 0 (tie-averaged expected 0.0 / 0.0).

| drug | rank (min-max) | top % |
|---|---|---|
| doxorubicin | 79-79 | 6.096 |
| vincristine | 815-815 | 62.886 |

Not in the candidate set: dactinomycin

### gold2

among approved (FDA label) candidates: 1646 candidates, 20 of 21 known drugs present, AUC 0.614468; certain hits at 10: 1, at 20: 1 (tie-averaged expected 1.0 / 1.0).

| drug | rank (min-max) | top % |
|---|---|---|
| carboplatin | 9-9 | 0.547 |
| doxorubicin | 89-89 | 5.407 |
| temozolomide | 146-146 | 8.87 |
| pemetrexed | 389-389 | 23.633 |
| ifosfamide | 406-406 | 24.666 |
| dactinomycin | 416-416 | 25.273 |
| atezolizumab | 501-501 | 30.437 |
| topotecan | 502-502 | 30.498 |
| cabozantinib | 532-532 | 32.321 |
| ipilimumab | 561-561 | 34.083 |
| irinotecan | 649-649 | 39.429 |
| daunorubicin | 704-704 | 42.77 |
| nivolumab | 778-778 | 47.266 |
| trastuzumab deruxtecan | 790-790 | 47.995 |
| etoposide | 913-913 | 55.468 |
| vincristine | 948-948 | 57.594 |
| ivosidenib | 1000-1000 | 60.753 |
| selinexor | 1015-1015 | 61.665 |
| selpercatinib | 1102-1103 | 66.95 |
| lurbinectedin | 1297-1297 | 78.797 |

Not in the candidate set: oxytetracycline

among approved small molecules: 1296 candidates, 15 of 21 known drugs present, AUC 0.567447; certain hits at 10: 1, at 20: 1 (tie-averaged expected 1.0 / 1.0).

| drug | rank (min-max) | top % |
|---|---|---|
| carboplatin | 9-9 | 0.694 |
| doxorubicin | 79-79 | 6.096 |
| temozolomide | 133-133 | 10.262 |
| pemetrexed | 359-359 | 27.701 |
| ifosfamide | 375-375 | 28.935 |
| topotecan | 457-457 | 35.262 |
| cabozantinib | 482-482 | 37.191 |
| irinotecan | 578-578 | 44.599 |
| daunorubicin | 624-624 | 48.148 |
| etoposide | 793-793 | 61.188 |
| vincristine | 815-815 | 62.886 |
| ivosidenib | 855-855 | 65.972 |
| selinexor | 867-867 | 66.898 |
| selpercatinib | 935-936 | 72.145 |
| lurbinectedin | 1070-1070 | 82.562 |

Not in the candidate set: atezolizumab, dactinomycin, ipilimumab, nivolumab, oxytetracycline, trastuzumab deruxtecan

### Temporal check (literature <= 2018, no prior)

Europe PMC records with PUB_YEAR <= 2018 only; knowledge-graph prior excluded (arm B). Labels: corrected gold 2 of the other diseases (v2 temporal protocol; post-2018 trials of OTHER diseases). Drugs first approved after 2018 in the universe: 343.

- approved_fda_label: later-trialled (gold 2) drugs present 20, AUC 0.671141, certain hits at 10/20: 1/2; gold 1 drugs: doxorubicin #4, vincristine #20, dactinomycin #188.
- approved_small_molecule: later-trialled (gold 2) drugs present 15, AUC 0.664039, certain hits at 10/20: 1/2; gold 1 drugs: doxorubicin #4, vincristine #18.

Top 20 with literature up to 2018 (approved small molecules):

| # | drug | score | gold 1 | later trialled (gold 2) | first approval | approved after 2018 |
|---|---|---|---|---|---|---|
| 1 | bortezomib | -2.165857 |  |  | n/a |  |
| 2 | alcohol | -2.246424 |  |  | 1982 |  |
| 3 | fluorouracil | -2.338013 |  |  | 1962 |  |
| 4 | doxorubicin | -2.613504 | x | x | 1974 |  |
| 5 | cabergoline | -2.647593 |  |  | 1996 |  |
| 6 | sorafenib | -2.648065 |  |  | 2005 |  |
| 7 | tamoxifen | -2.707734 |  |  | 1977 |  |
| 8 | vorinostat | -2.785101 |  |  | 2006 |  |
| 9 | folic acid | -2.820425 |  |  | 1982 |  |
| 10 | mitomycin | -2.824817 |  |  | 1981 |  |
| 11 | estradiol | -2.916232 |  |  | 1975 |  |
| 12 | danazol | -2.941429 |  |  | 1976 |  |
| 13 | everolimus | -2.948496 |  |  | 2009 |  |
| 14 | chloroquine | -2.98253 |  |  | 1949 |  |
| 15 | betaine | -3.014835 |  |  | 1996 |  |
| 16 | tyrosine | -3.152998 |  |  | 1982 |  |
| 17 | cyproheptadine | -3.222897 |  |  | 1961 |  |
| 18 | vincristine | -3.981111 | x | x | 1963 |  |
| 19 | helium | -3.981118 |  |  | 2015 |  |
| 20 | propranolol | -3.981148 |  |  | 1967 |  |

Drugs first approved after 2018 with a literature signal or a later trial (approved small molecules): selinexor (#1044, trialled), lurbinectedin (#1170, trialled), selpercatinib (#1251, trialled).

## 4. Mechanism classes of the disease genes

### cell_cycle_mitotic_checkpoint: DIS3L2, TRIP13

Retrieved papers naming a gene of this class: 37.

| drug | papers (class genes) | benefit | human | engine rank (small molecules) | gold 1 | gold 2 |
|---|---|---|---|---|---|---|
| tyrosine | 2 | 1 | 2 | 60 |  |  |
| gefitinib | 1 | 1 | 1 | 552 |  |  |
| cantharidin | 1 | 0 | 0 | 1019 |  |  |
| carboplatin | 1 | 0 | 0 | 9 |  | x |
| cisplatin | 1 | 0 | 0 | 11 |  |  |

### dna_repair_replication: BRCA2, DICER1

Retrieved papers naming a gene of this class: 48.

| drug | papers (class genes) | benefit | human | engine rank (small molecules) | gold 1 | gold 2 |
|---|---|---|---|---|---|---|
| sevoflurane | 1 | 1 | 0 | 5 |  |  |
| tamoxifen | 1 | 1 | 1 | 1 |  |  |
| testosterone | 1 | 0 | 0 | 8 |  |  |
| water | 1 | 0 | 0 | 1288 |  |  |

### extracellular_secreted: GPC4, IGF2

Retrieved papers naming a gene of this class: 60.

| drug | papers (class genes) | benefit | human | engine rank (small molecules) | gold 1 | gold 2 |
|---|---|---|---|---|---|---|
| oxygen | 2 | 0 | 0 | 1002 |  |  |
| danazol | 1 | 0 | 0 | 3 |  |  |
| imatinib | 1 | 0 | 1 | 234 |  |  |

### other: H19

Retrieved papers naming a gene of this class: 82.

| drug | papers (class genes) | benefit | human | engine rank (small molecules) | gold 1 | gold 2 |
|---|---|---|---|---|---|---|
| nitric oxide | 1 | 0 | 0 | 309 |  |  |
| oxygen | 1 | 0 | 0 | 1002 |  |  |
| tyrosine | 1 | 0 | 0 | 60 |  |  |

### receptor_signalling: GPC3

Retrieved papers naming a gene of this class: 22.

| drug | papers (class genes) | benefit | human | engine rank (small molecules) | gold 1 | gold 2 |
|---|---|---|---|---|---|---|
| lenvatinib | 1 | 1 | 0 | 4 |  |  |
| oxygen | 1 | 1 | 1 | 1002 |  |  |

### rna_processing_translation: WT1

Retrieved papers naming a gene of this class: 32.

| drug | papers (class genes) | benefit | human | engine rank (small molecules) | gold 1 | gold 2 |
|---|---|---|---|---|---|---|
| cytarabine | 1 | 1 | 1 | 7 |  |  |
| progesterone | 1 | 1 | 1 | 2 |  |  |
| azacitidine | 1 | 0 | 0 | 6 |  |  |
| gemcitabine | 1 | 0 | 1 | 17 |  |  |
| paclitaxel | 1 | 0 | 1 | 13 |  |  |

### transcription_chromatin_regulation: CTR9, PAX6, POU6F2, REST, TRIM28

Retrieved papers naming a gene of this class: 75.

| drug | papers (class genes) | benefit | human | engine rank (small molecules) | gold 1 | gold 2 |
|---|---|---|---|---|---|---|
| cocaine | 2 | 0 | 2 | 19 |  |  |
| paclitaxel | 1 | 1 | 1 | 13 |  |  |
| duloxetine | 1 | 0 | 0 | 803 |  |  |
| oxygen | 1 | 0 | 0 | 1002 |  |  |
| progesterone | 1 | 0 | 0 | 2 |  |  |

### Mitotic-checkpoint class (explicit; the class relevant to BUB1B in the MVA case)

Genes of kidney Wilms tumor in `cell_cycle_mitotic_checkpoint`: DIS3L2, TRIP13. Universe drugs whose ChEMBL mechanism is a mitotic-checkpoint, Aurora, PLK or TTK mechanism (mechanism list = caution_list class b (ChEMBL mechanism: Aurora, PLK, TTK/MPS1, BUB1, MAD2, mitotic checkpoint, spindle assembly)):

None: the approved-drug universe contains 0 drugs with such a mechanism (no approved mitotic-checkpoint, Aurora, PLK or TTK inhibitor). The approved drugs that act on the mitotic spindle are the tubulin binders:

Tubulin binders (caution_list class a; they act on the mitotic spindle):

| drug | type | small molecule + FDA | Wilms engine rank (all) | papers with a Wilms gene | gold 1 | gold 2 |
|---|---|---|---|---|---|---|
| paclitaxel | Small molecule | True | 14 | 2 |  |  | (label bands 0-4/5-11/12-17: not_labelled / not_labelled / not_labelled)
| vinblastine | Small molecule | True | 164 | 0 |  |  | (label bands 0-4/5-11/12-17: verify / verify / verify)
| ixabepilone | Small molecule | True | 411 | 0 |  |  | (label bands 0-4/5-11/12-17: not_labelled / not_labelled / not_labelled)
| vinorelbine | Small molecule | True | 538 | 0 |  |  | (label bands 0-4/5-11/12-17: not_labelled / not_labelled / not_labelled)
| colchicine | Small molecule | True | 809 | 0 |  |  | (label bands 0-4/5-11/12-17: partly (from 4 y) / labelled / labelled)
| cabazitaxel | Small molecule | True | 959 | 0 |  |  | (label bands 0-4/5-11/12-17: not_labelled / not_labelled / not_labelled)
| eribulin | Small molecule | True | 990 | 0 |  |  | (label bands 0-4/5-11/12-17: not_labelled / not_labelled / not_labelled)
| vindesine | Small molecule | False | 1063 | 0 |  |  |
| vincristine | Small molecule | True | 1170 | 0 | x | x | (label bands 0-4/5-11/12-17: verify / verify / verify)
| vinflunine | Small molecule | False | 1614 | 0 |  |  |
| docetaxel anhydrous | Small molecule | True | 2376 | 0 |  |  | (label bands 0-4/5-11/12-17: not_labelled / not_labelled / not_labelled)

### Mechanism x age-band table (top co-mentioned approved small molecules per class; label support)

| class | 0-4 y | 5-11 y | 12-17 y |
|---|---|---|---|
| cell_cycle_mitotic_checkpoint | not_labelled: carboplatin, cisplatin, gefitinib; verify: cantharidin, tyrosine | not_labelled: carboplatin, cisplatin, gefitinib; verify: cantharidin, tyrosine | not_labelled: carboplatin, cisplatin, gefitinib; verify: cantharidin, tyrosine |
| dna_repair_replication | not_labelled: tamoxifen; partly: sevoflurane; verify: testosterone, water | labelled: sevoflurane; not_labelled: tamoxifen; verify: testosterone, water | labelled: sevoflurane; not_labelled: tamoxifen; verify: testosterone, water |
| extracellular_secreted | not_labelled: danazol; verify: imatinib, oxygen | not_labelled: danazol; verify: imatinib, oxygen | not_labelled: danazol; verify: imatinib, oxygen |
| other | verify: nitric oxide, oxygen, tyrosine | verify: nitric oxide, oxygen, tyrosine | verify: nitric oxide, oxygen, tyrosine |
| receptor_signalling | not_labelled: lenvatinib; verify: oxygen | not_labelled: lenvatinib; verify: oxygen | not_labelled: lenvatinib; verify: oxygen |
| rna_processing_translation | not_labelled: azacitidine, gemcitabine, paclitaxel; verify: cytarabine, progesterone | not_labelled: azacitidine, gemcitabine, paclitaxel; verify: cytarabine, progesterone | not_labelled: azacitidine, gemcitabine, paclitaxel; verify: cytarabine, progesterone |
| transcription_chromatin_regulation | not_labelled: paclitaxel; verify: cocaine, duloxetine, oxygen, progesterone | not_labelled: paclitaxel; verify: cocaine, duloxetine, oxygen, progesterone | not_labelled: paclitaxel; verify: cocaine, duloxetine, oxygen, progesterone |

## 5. Evidence gaps (no dosing, no treatment recommendation)

Evidence rungs: L1 drug co-mentioned with a disease gene in retrieved abstracts; L2 at least one such paper with a sentence the rule lexicon reads as benefit; L3 at least one such paper with a human-patient or human-cell cue; L4 drug is an intervention of an EXPERIMENTAL arm of a registered phase 2+ study of the disease started after 2018; L5 current FDA label names the disease (openFDA indications text).

### #1 tamoxifen - level reached L3

- Missing: L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Other option: a registry of exposure and outcome inside an existing paediatric oncology cooperative-group protocol, to collect the human data that abstracts cannot provide.
- Label sentence (boxed_warning): "WARNING For Women with Ductal Carcinoma in Situ (DCIS) and Women at High Risk for Breast Cancer Serious and life threatening events associated with tamoxifen in the risk reduction setting (women at high risk for cancer and women with DCIS) include uterine malignancies, stroke and pulmonary embolism."
- Label sentence (boxed_warning): "Incidence rates for these events were estimated from the NSABP P-1 trial (see CLINICAL PHARMACOLOGY: Clinical Studies: Reduction in Breast Cancer Incidence in High Risk Women)."
- Label sentence (warnings, toxicidad_medular): "In patients with significant thrombocytopenia, rare hemorrhagic episodes have occurred, but it is uncertain if these episodes are due to tamoxifen therapy."
### #2 progesterone - level reached L3

- Missing: L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Other option: a registry of exposure and outcome inside an existing paediatric oncology cooperative-group protocol, to collect the human data that abstracts cannot provide.
- Label sentence (warnings, segundos_tumores): "Carcinogenesis, Mutagenesis, Impairment of Fertility Progesterone has not been tested for carcinogenicity in animals by the oral route of administration."
- Label sentence (warnings, embrio_fetal): "Pregnancy Category B: Reproductive studies have been performed in mice at doses up to 9 times the human oral dose, in rats at doses up to 44 times the human oral dose, in rabbits at a dose of 10 mcg/day delivered locally within the uterus by an implanted device, in guinea pigs at doses of approximat [...]"
### #3 danazol - level reached L1

- Missing: L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Label sentence (boxed_warning): "Use of danazol in pregnancy is contraindicated."
- Label sentence (boxed_warning): "A sensitive test (e.g., beta subunit test if available) capable of determining early pregnancy is recommended immediately prior to start of therapy."
- Label sentence (warnings, segundos_tumores): "Carcinogenesis, Mutagenesis, Impairment of Fertility Current data are insufficient to assess the carcinogenicity of danazol."
### #4 lenvatinib - level reached L2

- Missing: L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Label sentence (warnings, embrio_fetal): "(5.15) Embryo -F etal T oxicity: Can cause fetal harm."
- Label sentence (warnings, embrio_fetal): "1 6 Embryo -F etal Toxicity Based on its mechanism of action and data from animal reproduction studies, LENVIMA can cause fetal harm when administered to a pregnant woman."
- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.
### #5 sevoflurane - level reached L2

- Missing: L3; L4; L5.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Hazards: no hazard sentence retrieved from the label sections read: not read as safe.
### #6 azacitidine - level reached L1

- Missing: L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Label sentence (warnings, toxicidad_medular): "Anemia, Neutropenia and Thrombocytopenia: Monitor complete blood counts (CBC) frequently (5.2)."
- Label sentence (warnings, toxicidad_medular): "5.2 Anemia, Neutropenia and Thrombocytopenia Azacitidine causes anemia, neutropenia and thrombocytopenia in adult patients with MDS."
- Label sentence (warnings, embrio_fetal): "Embryo-Fetal Toxicity: Azacitidine can cause fetal harm."
- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.
### #7 cytarabine - level reached L3

- Missing: L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a registry of exposure and outcome inside an existing paediatric oncology cooperative-group protocol, to collect the human data that abstracts cannot provide.
- Label sentence (boxed_warning): "WARNING Only physicians experienced in cancer chemotherapy should use Cytarabine Injection."
- Label sentence (boxed_warning): "For induction therapy patients should be treated in a facility with laboratory and supportive resources sufficient to monitor drug tolerance and protect and maintain a patient compromised by drug toxicity."
- Label sentence (boxed_warning, toxicidad_medular): "The main toxic effect of cytarabine is bone marrow suppression with leukopenia, thrombocytopenia and anemia."
- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.
### #8 testosterone - level reached L1

- Missing: L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Hazards: no hazard sentence retrieved from the label sections read: not read as safe.
### #9 carboplatin - level reached L4

- Missing: L2; L3; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: completion and reporting of the registered trial(s) NCT04322318, NCT06401330; a label supplement for the disease needs the sponsor's regulatory submission, which is not something this analysis can supply.
- Label sentence (boxed_warning): "WARNINGS Carboplatin injection should be administered under the supervision of a qualified physician experienced in the use of cancer chemotherapeutic agents."
- Label sentence (boxed_warning): "Appropriate management of therapy and complications is possible only when adequate treatment facilities are readily available."
- Label sentence (boxed_warning, toxicidad_medular): "Bone marrow suppression is dose related and may be severe, resulting in infection and/or bleeding."
- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.
### #10 zinc chloride - level reached L0

- Missing: L1; L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Label sentence (warnings, segundos_tumores): "Carcinogenesis, Mutagenesis, Impairment of Fertility Long-term animal studies to evaluate the carcinogenic potential of Zinc Chloride Injection, USP, 1 mg/mL have not been performed, nor have studies been done to assess mutagenesis or impairment of fertility."
- Label sentence (warnings, embrio_fetal): "It is also not known whether zinc chloride can cause fetal harm when administered to a pregnant woman or can affect reproduction capacity."
### #11 cisplatin - level reached L1

- Missing: L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Label sentence (boxed_warning): "WARNING: NEPHROTOXICITY, PERIPHERAL NEUROPATHY, NAUSEA AND VOMITING and MYELOSUPPRESSION. • Nephrotoxicity: cisplatin for injection can cause severe renal toxicity, including acute renal failure."
- Label sentence (boxed_warning): "Severe renal toxicities are dose-related and cumulative."
- Label sentence (boxed_warning, toxicidad_medular): "WARNING: NEPHROTOXICITY, PERIPHERAL NEUROPATHY, NAUSEA AND VOMITING and MYELOSUPPRESSION. • Nephrotoxicity: cisplatin for injection can cause severe renal toxicity, including acute renal failure."
- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.
### #12 palbociclib - level reached L1

- Missing: L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Label sentence (warnings, toxicidad_medular): "5 WARNINGS AND PRECAUTIONS •Neutropenia: Monitor complete blood count prior to start of IBRANCE therapy and at the beginning of each cycle, as well as on Day 15 of the first 2 cycles, and as clinically indicated."
- Label sentence (warnings, toxicidad_medular): "(5.3, 8.1, 8.3) 5.1 Neutropenia Neutropenia was the most frequently reported adverse reaction in PALOMA-2 with an incidence of 80% and PALOMA-3 with an incidence of 83%."
- Label sentence (warnings, embrio_fetal): "(5.2) •Embryo-Fetal Toxicity: IBRANCE can cause fetal harm."
- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.
### #13 paclitaxel - level reached L3

- Missing: L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a registry of exposure and outcome inside an existing paediatric oncology cooperative-group protocol, to collect the human data that abstracts cannot provide.
- Label sentence (boxed_warning): "WARNING Paclitaxel should be administered under the supervision of a physician experienced in the use of cancer chemotherapeutic agents."
- Label sentence (boxed_warning): "Appropriate management of complications is possible only when adequate diagnostic and treatment facilities are readily available."
- Label sentence (boxed_warning, toxicidad_medular): "In order to monitor the occurrence of bone marrow suppression, primarily neutropenia, which may be severe and result in infection, it is recommended that frequent peripheral blood cell counts be performed on all patients receiving paclitaxel."
- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.
### #14 topiramate - level reached L0

- Missing: L1; L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Label sentence (warnings, embrio_fetal): "5 WARNINGS AND PRECAUTIONS Acute myopia and secondary angle closure glaucoma: can lead to permanent visual loss; discontinue topiramate extended-release capsules as soon as possible (5.1) Visual field defects: consider discontinuation of topiramate extended-release capsules (5.2) Oligohydrosis and h [...]"
- Label sentence (warnings, embrio_fetal): "5.7 Fetal Toxicity Topiramate extended-release capsules can cause fetal harm when administered to a pregnant woman."
### #15 fostamatinib - level reached L0

- Missing: L1; L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Label sentence (warnings, toxicidad_medular): "(5.3) Neutropenia: Monitor ANC monthly, and for infection."
- Label sentence (warnings, toxicidad_medular): "5.4 Neutropenia Neutropenia occurred in 6% of patients treated with TAVALISSE; febrile neutropenia occurred in 1% of patients."
- Label sentence (warnings, segundos_tumores): "In subchronic, chronic, and carcinogenicity studies of TAVALISSE, chondrodystrophy of the femoral head was seen in rodents."
### #16 diethylstilbestrol - level reached L0

- Missing: L1; L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Hazards: no hazard sentence retrieved from the label sections read: not read as safe.
### #17 gemcitabine - level reached L3

- Missing: L2; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a registry of exposure and outcome inside an existing paediatric oncology cooperative-group protocol, to collect the human data that abstracts cannot provide.
- Label sentence (warnings, toxicidad_medular): "(5.1) Myelosuppression: Monitor for myelosuppression prior to each cycle and reduce or withhold dose for severe myelosuppression."
- Label sentence (warnings, toxicidad_medular): "(5.10) 5.1 Schedule-Dependent Toxicity In clinical trials evaluating the maximum tolerated dose of Gemcitabine for Injection, prolongation of the infusion time beyond 60 minutes or more frequent than weekly dosing resulted in an increased incidence of clinically significant hypotension, severe flu-l [...]"
- Label sentence (warnings, embrio_fetal): "(5.6) Embryo-Fetal Toxicity: Can cause fetal harm."
- Chemotherapy caution: oncology drug with serious toxicity; not a recommendation.
### #18 cannabidiol - level reached L0

- Missing: L1; L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Label sentence (warnings, embrio_fetal): "8 USE IN SPECIFIC POPULATIONS Pregnancy: Based on animal data, may cause fetal harm."
- Label sentence (warnings, embrio_fetal): "Administration of cannabidiol to pregnant animals produced evidence of developmental toxicity (increased embryofetal mortality in rats and decreased fetal body weights in rabbits; decreased growth, delayed sexual maturation, long-term neurobehavioral changes, and adverse effects on the reproductive  [...]"
### #19 cocaine - level reached L3

- Missing: L2; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Other option: a registry of exposure and outcome inside an existing paediatric oncology cooperative-group protocol, to collect the human data that abstracts cannot provide.
- Label sentence (boxed_warning): "WARNING: ABUSE AND DEPENDENCE CNS stimulants, including cocaine hydrochloride, have a high potential for abuse and dependence [see Warning and Precautions (5.1).] WARNING: ABUSE AND DEPENDENCE See full prescribing information for complete boxed warning."
- Label sentence (boxed_warning): "CNS stimulants, including cocaine hydrochloride, have a high potential for abuse and dependence."
- Label sentence (warnings, embrio_fetal): "8 USE IN SPECIFIC POPULATIONS • Pregnancy: May cause fetal harm."
### #20 spironolactone - level reached L0

- Missing: L1; L2; L3; L4; L5. Age bands without labelled support: 0-4 y, 5-11 y, 12-17 y.
- Single study or data that would close the gap: a registered, ethics-approved adaptive pilot (phase 1b/2, window or basket design) in relapsed or refractory Wilms tumour, with pre-specified response and toxicity endpoints and an independent safety committee.
- Other option: a single-patient (N-of-1) study design under ethics oversight, only after a mechanism and safety review.
- Label sentence (warnings, embrio_fetal): "8 USE IN SPECIFIC POPULATIONS Pregnancy: Based on animal data, spironolactone may affect sex differentiation of the male during embryogenesis (8.1) 8.1 Pregnancy Risk Summary Based on mechanism of action and findings in animal studies, spironolactone may affect sex differentiation of the male during [...]"
- Label sentence (warnings, embrio_fetal): "Rat embryofetal studies report feminization of male fetuses and endocrine dysfunction in females exposed to spironolactone in utero."

## Limits

- Label support is label text read by a heuristic classifier unless a manual review exists in config; manual_review_required marks the former.
- Co-mention counts are among at most 300 retrieved abstracts for 14 genes together; they are not disease-specific effect evidence.
- ChEMBL molecule type is used as recorded (for example dactinomycin is typed Protein and is therefore dropped by the small-molecule filter).
- openFDA label existence uses current labels by generic-name field and may include combination or homeopathic products that name the substance.
- Wilms tumour is a v2 disease; this is a worked case, not part of the pre-registered primary validation, and no threshold is evaluated here.
