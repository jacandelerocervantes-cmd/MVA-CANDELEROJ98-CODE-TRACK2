# mosaic variegated aneuploidy (MVA): application of the frozen engine (B0)

**application of the frozen engine; exploratory; not a candidate list.** Not a validation, no claim of efficacy, no dosing, no treatment recommendation. Protocol: `protocol/addendum_engine_apply.md`.

**Literature flag: literature volume above the dilution threshold: the literature features carry part of the ranking** (248 retrieved papers carry a gene symbol, threshold <= 29; Europe PMC hit count 2047, stratum H).

## 1. Inputs

* Genes scored (tier curated): BUB1, BUB1B, BUB3, CENATAC, CEP57, MAD1L1, TRIP13
* Genes proposed, not curated, not scored: CEP192, MAD2L1BP, SMC5
* Mechanism classes: centrosome/centriole: CENATAC, CEP57; spindle-assembly checkpoint: BUB1, BUB1B, BUB3, MAD1L1, TRIP13
* Frozen model sha256 `09e3a8e8f5a32ea7`, universe 2705 drugs, literature <= 2024

## 2. Filter chain

| step | drugs |
|---|---|
| universe | 2705 |
| with FDA label (openFDA match 1641, +2 by Drugs@FDA check of 9 asked) | 1643 |
| small molecule (ChEMBL type) | 1293 |

Dropped by type after the label step: {"biologic": 286, "oligonucleotide": 11, "unknown_type": 53}

## 3. Top 20 (FDA-approved small molecules), frozen B0 score

| # | drug | score | Platt | approval | chemo caution | 0-4 y | 5-11 y | 12-17 y | gene papers (benefit/human) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | olaparib | -2.228689 | 0.069514 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 5 (1/4) |
| 2 | tyrosine | -3.065231 | 0.034803 | openFDA label match + Drugs@FDA: no Prescription product | no | verify | verify | verify | 12 (4/4) |
| 3 | carboplatin | -3.566842 | 0.022767 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 2 (0/1) |
| 4 | progesterone | -3.730856 | 0.019797 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 4 (0/0) |
| 5 | zinc chloride | -4.555086 | 0.009757 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |
| 6 | cabazitaxel | -4.666486 | 0.008863 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 1 (0/1) |
| 7 | vorinostat | -4.710194 | 0.008535 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 0 (0/0) |
| 8 | fingolimod | -4.710374 | 0.008533 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |
| 9 | arsenic trioxide | -4.721718 | 0.00845 | Drugs@FDA Prescription | YES | verify | verify | verify | 0 (0/0) |
| 10 | cisplatin | -4.753524 | 0.008221 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 9 (1/3) |
| 11 | doxorubicin | -4.898435 | 0.007254 | openFDA label match + Drugs@FDA Prescription | YES | age_not_stated | age_not_stated | age_not_stated | 2 (0/1) |
| 12 | fostamatinib | -4.903822 | 0.00722 | openFDA label match + Drugs@FDA Prescription | no | not_labelled | not_labelled | not_labelled | 0 (0/0) |
| 13 | belinostat | -4.911879 | 0.00717 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 0 (0/0) |
| 14 | paclitaxel | -4.91242 | 0.007167 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 7 (3/7) |
| 15 | carvedilol | -4.93395 | 0.007034 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |
| 16 | pravastatin | -4.933983 | 0.007034 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |
| 17 | lovastatin | -4.935562 | 0.007025 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |
| 18 | apomorphine | -4.936769 | 0.007017 | openFDA label match + Drugs@FDA Prescription | no | not_labelled | not_labelled | not_labelled | 0 (0/0) |
| 19 | methyldopa | -4.96075 | 0.006873 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |
| 20 | fluvastatin | -4.975177 | 0.006788 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |

Ties share a score; `rank_min` is in the JSON. 'verify' = needs manual verification; 'not_labelled' = label does not establish use in the band; 'age_not_stated' = established use without a minimum age. Unknown age is never assumed.

## 4. Label quotes (verbatim, truncated where marked)

### 1. olaparib (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (openFDA, set id 741ff3e3-dc1a-45a6-84e5-2481b27131aa, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness of Lynparza have not been established in pediatric patients."
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] Monitor complete blood count for cytopenia at baseline and monthly thereafter for clinically significant changes during treatment.
  * [toxicidad_medular, warnings] If the levels have not recovered to Grade 1 or less after 4 weeks, refer the patient to a hematologist for further investigations, including bone marrow analysis and blood sample for cytogenetics.
  * [segundos_tumores, warnings] 5 WARNINGS AND PRECAUTIONS • Myelodysplastic Syndrome/Acute Myeloid Leukemia (MDS/AML): Occurred in approximately 1.2% of patients with various BRCA m, g BRCA m, HRR gene-mutated or HRD-positive cancers exposed to Lynparza and the majority of events had a fatal outcome.
  * [segundos_tumores, warnings] ( 5.5 , 8.1 , 8.3 ) 5.1 Myelodysplastic Syndrome/Acute Myeloid Leukemia Myelodysplastic syndrome (MDS)/Acute Myeloid Leukemia (AML) has occurred in patients treated with Lynparza and some cases were fatal.

### 2. tyrosine (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); no readable paediatric-use section; manual review required)
* Paediatric use: no readable paediatric-use sentence (not labelled or needs manual verification)

### 3. carboplatin (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 4c89cedc-e48b-445c-ad05-54dbe60e4fb0, section pediatric_use): "Pediatric Use Safety and effectiveness in pediatric patients have not been established (seeWARNINGS : "audiologic toxicity")."
* Paediatric warning: Pediatric Use Safety and effectiveness in pediatric patients have not been established (seeWARNINGS : "audiologic toxicity").
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNINGS Carboplatin injection should be administered under the supervision of a qualified physician experienced in the use of cancer chemotherapeutic agents.
  * [boxed_warning, boxed_warning] Appropriate management of therapy and complications is possible only when adequate treatment facilities are readily available.
  * [toxicidad_medular, boxed_warning] Bone marrow suppression is dose related and may be severe, resulting in infection and/or bleeding.
  * [toxicidad_medular, warnings] WARNINGS Bone marrow suppression (leukopenia, neutropenia, and thrombocytopenia) is dose-dependent and is also the dose-limiting toxicity.

### 4. progesterone (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id d7619751-47b4-47c9-b716-e2a704b350c8, section pediatric_use): "Pediatric Use Progesterone capsules are not indicated in children."

### 5. zinc chloride (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 2bd646d2-2a2a-4376-8da0-c7f08b0511ac, section pediatric_use): "Pediatric Use See DOSAGE AND ADMINISTRATION section."

### 6. cabazitaxel (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id ec6892c8-01e1-6a8b-e053-2995a90a956f, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of cabazitaxel in pediatric patients have not been established."
* Paediatric warning: The maximum tolerated dose (MTD) was 30 mg/m 2 intravenously over 1 hour on Day 1 of a 21 day cycle in pediatric patients with solid tumors based on the dose-limiting toxicity (DLT) of febrile neutropenia.
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNING: NEUTROPENIA AND HYPERSENSITIVITY WARNING: NEUTROPENIA AND HYPERSENSITIVITY See full prescribing information for complete boxed warning .
  * [boxed_warning, boxed_warning] Neutropenic deaths have been reported.
  * [toxicidad_medular, boxed_warning] WARNING: NEUTROPENIA AND HYPERSENSITIVITY WARNING: NEUTROPENIA AND HYPERSENSITIVITY See full prescribing information for complete boxed warning .
  * [toxicidad_medular, boxed_warning] Neutropenic deaths have been reported.

### 7. vorinostat (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (openFDA, set id cd86ee78-2781-468b-930c-3c4677bcc092, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of ZOLINZA in pediatric patients have not been established."
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] ( 5.1 ) Myelosuppression: Thrombocytopenia and anemia may require dose modification or discontinuation.
  * [toxicidad_medular, warnings] ( 5.5 ) Severe Thrombocytopenia with Concomitant Use of other HDAC Inhibitors: Severe thrombocytopenia with gastrointestinal bleeding has been reported with concomitant use of ZOLINZA and other HDAC inhibitors (e.g., valproic acid).
  * [embrio_fetal, warnings] ( 5.6 , 7.2 ) Embryo-Fetal Toxicity: Fetal harm can occur when administered to a pregnant woman.
  * [embrio_fetal, warnings] Women should be apprised of the potential harm to the fetus.

### 8. fingolimod (status MIXED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id a545c6fc-7006-52c2-f476-34e9801ac04f, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness of Fingolimod Capsules for the treatment of relapsing forms of multiple sclerosis in pediatric patients 10 to less than 18 years of age were established in one randomized, double-blind clinical study in 215 patients (Fingolimod Capsules n = 107; intramuscular interferon (IFN) beta-1a n = 108) [see Clinical Studies (14.2)]."
* Paediatric warning: Juvenile Animal Toxicity Data In a study in which fingolimod (0.3, 1.5, or 7.5 mg/kg/day) was orally administered to young rats from weaning through sexual maturity, changes in bone mineral density and persistent neurobehavioral impairment (altered auditory startle) were observed at all doses.

### 9. arsenic trioxide (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id c99b01b4-67ff-4793-8375-fe3a21c6c9be, section pediatric_use): "8.4 Pediatric Use The safety and efficacy of arsenic trioxide injection as a single agent for treatment of pediatric patients with relapsed or refractory APL is supported by the pivotal phase 2 study in 40 patients with relapsed or refractory APL."
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] BOXED WARNING Differentiation Syndrome: Patients with acute promyelocytic leukemia (APL) treated with arsenic trioxide injection have experienced differentiation syndrome, which may be life-threatening or fatal.
  * [boxed_warning, boxed_warning] Signs and symptoms may include unexplained fever, dyspnea, hypoxia, acute respiratory distress, pulmonary infiltrates, pleural or pericardial effusions, weight gain, peripheral edema, hypotension, renal insufficiency, hepatopathy, and multi-organ dysfunction, in the presence or absence of leukocytosis.
  * [segundos_tumores, warnings] Monitor patients for the development of second primary malignancies.
  * [segundos_tumores, warnings] Monitor patients for the development of second primary malignancies.

### 10. cisplatin (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 508496cb-3441-46b3-a4fe-e0d440e6adc6, section pediatric_use): "8.4 Pediatric Use Ototoxic effects may be more severe and detrimental in pediatric patients receiving cisplatin for injection, particularly in patients less than 5 years of age."
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNING: NEPHROTOXICITY, PERIPHERAL NEUROPATHY, NAUSEA AND VOMITING and MYELOSUPPRESSION. • Nephrotoxicity: cisplatin for injection can cause severe renal toxicity, including acute renal failure.
  * [boxed_warning, boxed_warning] Severe renal toxicities are dose-related and cumulative.
  * [toxicidad_medular, boxed_warning] WARNING: NEPHROTOXICITY, PERIPHERAL NEUROPATHY, NAUSEA AND VOMITING and MYELOSUPPRESSION. • Nephrotoxicity: cisplatin for injection can cause severe renal toxicity, including acute renal failure.
  * [toxicidad_medular, boxed_warning] Use highly effective antiemetic premedication [see Dosage and Administration (2.1) and Warnings and Precautions (5.3)]. • Myelosuppression: cisplatin for injection can cause severe myelosuppression with fatalities due to infections.

### 11. doxorubicin (status ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 1fd148fb-0fbc-4b6f-b790-23546fb46a71, section pediatric_use): "8.4 Pediatric Use Based on postmarketing reports, pediatric patients treated with doxorubicin hydrochloride are at risk for developing late cardiovascular dysfunction."
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNING: CARDIOMYOPATHY, SECONDARY MALIGNANCIES, EXTRAVASATION AND TISSUE NECROSIS, and SEVERE MYELOSUPPRESSION • Cardiomyopathy: Myocardial damage, including acute left ventricular failure, can occur with doxorubicin hydrochloride.
  * [boxed_warning, boxed_warning] The risk of cardiomyopathy is proportional to the cumulative exposure with incidence rates from 1%–20% for cumulative doses ranging from 300 mg/m2 to 500 mg/m2 when doxorubicin hydrochloride is administered every 3 weeks.
  * [toxicidad_medular, boxed_warning] WARNING: CARDIOMYOPATHY, SECONDARY MALIGNANCIES, EXTRAVASATION AND TISSUE NECROSIS, and SEVERE MYELOSUPPRESSION • Cardiomyopathy: Myocardial damage, including acute left ventricular failure, can occur with doxorubicin hydrochloride.
  * [toxicidad_medular, boxed_warning] Immediately terminate the drug and apply ice to the affected area [see Warnings and Precautions (5.3)]. • Severe myelosuppression resulting in serious infection, septic shock, requirement for transfusions, hospitalization, and death may occur [see Warnings and Precautions (5.4)].

### 12. fostamatinib (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 21149cc3-049b-43e2-b141-c9499160556c, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness in pediatric patients have not been established."
* Paediatric warning: In subchronic, chronic, and carcinogenicity studies of TAVALISSE, chondrodystrophy of the femoral head was seen in rodents.
* Paediatric warning: In a study in juvenile rabbits, growth plate dysplasia was observed in the proximal femur and femoro-tibial joint, and bone marrow cellularity was reduced in the femur and sternum.

### 13. belinostat (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 2e8ef36b-71fa-4492-a16f-577d5f7d111d, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of Beleodaq in pediatric patients have not been established."
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] 5 WARNINGS AND PRECAUTIONS Hematologic Toxicity: Thrombocytopenia, leukopenia (neutropenia and lymphopenia), and anemia: Monitor blood counts and modify dosage for hematologic toxicities.
  * [toxicidad_medular, warnings] (5.6) 5.1 Hematologic Toxicity Beleodaq can cause thrombocytopenia, leukopenia (neutropenia and lymphopenia), and/or anemia; monitor blood counts weekly during treatment, and modify dosage as necessary [see Dosage and Administration (2.2 ) and Adverse Reactions(6.1)].
  * [embrio_fetal, warnings] Embryo-Fetal Toxicity: Can cause fetal harm.
  * [embrio_fetal, warnings] 5.6 Embryo-fetal Toxicity Based on its mechanism of action and findings of genotoxicity, Beleodaq can cause fetal harm when administered to a pregnant woman [see Clinical Pharmacology (12.1) and Nonclinical Toxicology (13.1)].

### 14. paclitaxel (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 0c12c77c-5fd4-48dc-a318-e96892233992, section pediatric_use): "Pediatric Use The safety and effectiveness of paclitaxel in pediatric patients have not been established."
* Paediatric warning: There have been reports of central nervous system (CNS) toxicity (rarely associated with death) in a clinical trial in pediatric patients in which paclitaxel was infused intravenously over 3 hours at doses ranging from 350 mg/m2 to 420 mg/m2.
* Paediatric warning: The toxicity is most likely attributable to the high dose of the ethanol component of the paclitaxel vehicle given over a short infusion time.
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNING Paclitaxel should be administered under the supervision of a physician experienced in the use of cancer chemotherapeutic agents.
  * [boxed_warning, boxed_warning] Appropriate management of complications is possible only when adequate diagnostic and treatment facilities are readily available.
  * [toxicidad_medular, boxed_warning] In order to monitor the occurrence of bone marrow suppression, primarily neutropenia, which may be severe and result in infection, it is recommended that frequent peripheral blood cell counts be performed on all patients receiving paclitaxel.
  * [toxicidad_medular, warnings] Bone marrow suppression (primarily neutropenia) is dose-dependent and is the dose-limiting toxicity.

### 15. carvedilol (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); no readable paediatric-use section; manual review required)
* Paediatric use: no readable paediatric-use sentence (not labelled or needs manual verification)

### 16. pravastatin (status MIXED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 42ccbca7-fa52-4df9-8f1a-ea7da7312883, section pediatric_use): "8.4 Pediatric use The safety and effectiveness of pravastatin as an adjunct to diet to reduce LDL-C have been established in pediatric patients 8 years of age and older with HeFH."

### 17. lovastatin (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); no readable paediatric-use section; manual review required)
* Paediatric use: no readable paediatric-use sentence (not labelled or needs manual verification)

### 18. apomorphine (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id aedf70e4-1949-44bc-bf1b-0f282c4f0705, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness in pediatric patients have not been established."

### 19. methyldopa (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 9e22e1d8-04ec-4b67-92ea-b09c706f2703, section pediatric_use): "Pediatric Use There are no well-controlled clinical trials in pediatric patients."

### 20. fluvastatin (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id ccd6bbe1-e293-4e1c-9ec4-e6c217612e77, section pediatric_use): "8.4 Pediatric Use The safety and efficacy of fluvastatin sodium in children and adolescent patients 9 to 16 years of age with heterozygous familial hypercholesterolemia have been evaluated in open-label, uncontrolled clinical trials for a duration of two years."

## 5. Mechanism class x age band

Per class: the five small-molecule FDA drugs with most papers co-mentioning the drug and a gene of the class.

### centrosome/centriole (genes CENATAC, CEP57; 12 retrieved papers mention a class gene)

| drug | papers | benefit | human | engine rank (small molecule) | 0-4 y | 5-11 y | 12-17 y |
|---|---|---|---|---|---|---|---|
| (no small-molecule FDA drug co-mentioned with a class gene) | | | | | | | |

Band cells by state (drug names):

* 0-4 y: none
* 5-11 y: none
* 12-17 y: none

### spindle-assembly checkpoint (genes BUB1, BUB1B, BUB3, MAD1L1, TRIP13; 241 retrieved papers mention a class gene)

| drug | papers | benefit | human | engine rank (small molecule) | 0-4 y | 5-11 y | 12-17 y |
|---|---|---|---|---|---|---|---|
| tyrosine | 12 | 4 | 4 | 2 | verify | verify | verify |
| cisplatin | 9 | 1 | 3 | 10 | not_labelled | not_labelled | not_labelled |
| paclitaxel | 7 | 3 | 7 | 14 | not_labelled | not_labelled | not_labelled |
| olaparib | 5 | 1 | 4 | 1 | not_labelled | not_labelled | not_labelled |
| progesterone | 4 | 0 | 0 | 4 | verify | verify | verify |

Band cells by state (drug names):

* 0-4 y: not_labelled: cisplatin, olaparib, paclitaxel; verify: progesterone, tyrosine
* 5-11 y: not_labelled: cisplatin, olaparib, paclitaxel; verify: progesterone, tyrosine
* 12-17 y: not_labelled: cisplatin, olaparib, paclitaxel; verify: progesterone, tyrosine

## 6. The nine drugs of the report in this ranking (comparison only, not a validation)

| drug | universe name | score | rank in universe | rank among FDA-label | rank among FDA small molecules | in top 20 | 0-4 y | 5-11 y | 12-17 y |
|---|---|---|---|---|---|---|---|---|---|
| sirolimus | sirolimus | -6.824432 | 1449 of 2705 | 1144 of 1643 | 941 of 1293 | no | not_labelled | not_labelled | partly (from 13 y) |
| everolimus | everolimus | -6.428228 | 648 of 2705 | 557 of 1643 | 484 of 1293 | no | partly (from 1 y) | labelled | labelled |
| metformin | metformin | -7.078504 | 1844 of 2705 | 1363 of 1643 | 1107 of 1293 | no | not_labelled | partly (from 10 y) | labelled |
| sodium phenylbutyrate | phenylbutanoic acid | -7.142514 | 1917 of 2705 | None of 1643 | filtered out of 1293 | no | not read | not read | not read |
| trametinib | trametinib | -5.638962 | 86 of 2705 | 81 of 1643 | 77 of 1293 | no | partly (from 1 y) | labelled | labelled |
| selumetinib | selumetinib | -5.71584 | 108 of 2705 | 100 of 1643 | 93 of 1293 | no | partly (from 1 y) | labelled | labelled |
| bortezomib | bortezomib | -5.321743 | 57 of 2705 | 55 of 1643 | 52 of 1293 | no | not_labelled | not_labelled | not_labelled |
| chloroquine | chloroquine | -6.919698 | 1592 of 2705 | 1226 of 1643 | 1004 of 1293 | no | verify | verify | verify |
| bosutinib | bosutinib | -5.787801 | 126 of 2705 | 116 of 1643 | 107 of 1293 | no | partly (from 1 y) | labelled | labelled |

## 7. Evidence gap (computed; no dosing; no recommendation)

Orphanet case-count tier: >5 (Orphanet prevalence records of ORPHA:1052).

* **number of patients available for any study** (the proband is one known case (stated by the project; not computed from any record)). Design: Single-patient (N-of-1) crossover design with a pre-specified outcome, analysed with a hierarchical model that borrows from the mechanism class.
* **no label to measure the engine on this disease; only the class-level behaviour is known** (class pool cell_cycle_mitotic_checkpoint: mean AUC 0.686759, tau2 0.058352, k=3). Design: Bayesian hierarchical model borrowing from the mechanism class pool; any study result of one disease of the class updates the class level (hypothesis for study design only).
* **age band 0-4 y without labelled support for 20 of the 20 top-ranked drugs** (states verify, not_labelled or age_not_stated). Design: paediatric extrapolation or an age-stratified label supplement is a sponsor/regulatory step; within a study, the N-of-1 or series design of the first gap applies with an independent safety committee.
* **age band 5-11 y without labelled support for 20 of the 20 top-ranked drugs** (states verify, not_labelled or age_not_stated). Design: paediatric extrapolation or an age-stratified label supplement is a sponsor/regulatory step; within a study, the N-of-1 or series design of the first gap applies with an independent safety committee.
* **age band 12-17 y without labelled support for 20 of the 20 top-ranked drugs** (states verify, not_labelled or age_not_stated). Design: paediatric extrapolation or an age-stratified label supplement is a sponsor/regulatory step; within a study, the N-of-1 or series design of the first gap applies with an independent safety committee.

## 8. Literature volume

* Papers retrieved: 300; with a gene symbol: 248; Europe PMC hit count: 2047
* Drugs co-mentioned at least once: 20; sum of paper counts over drugs: 57
* **literature volume above the dilution threshold: the literature features carry part of the ranking**
* Dilution result (benchmark, held-out): B0 gain over popularity by mean papers per disease: 5 papers: 0.051009; 9 papers: 0.064084; 29 papers: 0.093531; 87 papers: 0.119834; 289 papers: 0.127407. Prior-only gain (flat over volume): 5: 0.038025; 9: 0.038025; 29: 0.038025; 87: 0.038025; 289: 0.038025

## 9. Limits

* Application of the frozen B0 engine; exploratory; not a candidate list. It is not a validation and makes no claim of efficacy.
* The engine was validated on diseases with known drugs (benchmark gold from approved indications and trials). Where no drug is known the ranking cannot be checked.
* On Wilms tumour the engine failed to rank the known drugs of the disease in its top 20; a high rank is not evidence that a drug works and a low rank is not evidence against it.
* Literature counts are co-mentions among at most 300 abstracts; they are not effect evidence. The knowledge-graph prior is target based.
* ChEMBL molecule type and openFDA label presence are used as recorded; Drugs@FDA is a confirmation of marketing status only.
* Label age support is label text read by a heuristic classifier unless config/pediatric_label_review.json holds a manual entry (manual review required otherwise); the age of a patient is never assumed.
* Oncology-cytotoxic drugs carry serious label-documented toxicity (marrow suppression, second malignancy, embryo-fetal harm). In Fanconi anaemia DNA-crosslinking agents are of particular concern because of cellular hypersensitivity. This is a flag, not advice.
* No dosing and no treatment recommendation. Hypotheses only.
