# Fanconi anaemia: application of the frozen engine (B0)

**application of the frozen engine; exploratory; not a candidate list.** Not a validation, no claim of efficacy, no dosing, no treatment recommendation. Protocol: `protocol/addendum_engine_apply.md`.

**Literature flag: literature volume above the dilution threshold: the literature features carry part of the ranking** (296 retrieved papers carry a gene symbol, threshold <= 29; Europe PMC hit count 31219, stratum H).

## 1. Inputs

* Genes scored (tier curated): BRCA1, BRCA2, BRIP1, ERCC4, FAAP100, FANCA, FANCB, FANCC, FANCD2, FANCE, FANCF, FANCG, FANCI, FANCL, FANCM, MAD2L2, PALB2, RAD51, RAD51C, RFWD3, SLX4, UBE2T, XRCC2
* Genes proposed, not curated, not scored: COX4I1, ESCO2, NHEJ1
* Mechanism classes: DNA interstrand-crosslink repair: ERCC4, FAAP100, FANCA, FANCB, FANCC, FANCD2, FANCE, FANCF, FANCG, FANCI, FANCL, FANCM, RAD51, RFWD3; DNA repair/replication: BRCA1, BRCA2, BRIP1, PALB2, RAD51C, SLX4, UBE2T, XRCC2; spindle-assembly checkpoint: MAD2L2
* Frozen model sha256 `09e3a8e8f5a32ea7`, universe 2705 drugs, literature <= 2024

## 2. Filter chain

| step | drugs |
|---|---|
| universe | 2705 |
| with FDA label (openFDA match 1641, +3 by Drugs@FDA check of 10 asked) | 1644 |
| small molecule (ChEMBL type) | 1294 |

Dropped by type after the label step: {"biologic": 286, "oligonucleotide": 11, "unknown_type": 53}

## 3. Top 20 (FDA-approved small molecules), frozen B0 score

| # | drug | score | Platt | approval | chemo caution | 0-4 y | 5-11 y | 12-17 y | gene papers (benefit/human) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | niraparib | -0.786029 | 0.207851 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 3 (1/1) |
| 2 | tamoxifen | -0.882718 | 0.194329 | openFDA label match + Drugs@FDA Prescription | no | not_labelled | not_labelled | not_labelled | 3 (2/2) |
| 3 | olaparib | -1.742856 | 0.102374 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 11 (3/5) |
| 4 | carboplatin | -1.879025 | 0.09198 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 3 (1/2) |
| 5 | danazol | -2.123916 | 0.075652 | openFDA label match + Drugs@FDA Prescription | no | not_labelled | not_labelled | not_labelled | 1 (1/1) |
| 6 | talazoparib | -2.603087 | 0.051164 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 1 (1/0) |
| 7 | hydroxyurea | -2.85042 | 0.041664 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 2 (2/0) |
| 8 | anastrozole | -3.207039 | 0.030885 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 1 (1/0) |
| 9 | rucaparib | -3.330839 | 0.027817 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 2 (1/1) |
| 10 | mitomycin | -3.404624 | 0.026131 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 3 (0/1) |
| 11 | cisplatin | -3.631285 | 0.021551 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 4 (1/2) |
| 12 | topotecan | -3.897605 | 0.017167 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 1 (1/0) |
| 13 | paclitaxel | -3.908163 | 0.017013 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 1 (1/1) |
| 14 | fostamatinib | -3.985156 | 0.015927 | openFDA label match + Drugs@FDA Prescription | no | not_labelled | not_labelled | not_labelled | 0 (0/0) |
| 15 | irinotecan | -4.165248 | 0.013647 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 1 (1/0) |
| 16 | zinc chloride | -4.183068 | 0.01344 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |
| 17 | gemcitabine | -4.295755 | 0.012199 | openFDA label match + Drugs@FDA Prescription | YES | not_labelled | not_labelled | not_labelled | 1 (1/1) |
| 18 | estradiol | -4.504182 | 0.010195 | openFDA label match + Drugs@FDA Prescription | no | not_labelled | not_labelled | not_labelled | 1 (0/0) |
| 19 | liothyronine | -4.524504 | 0.010018 | openFDA label match + Drugs@FDA Prescription | no | verify | verify | verify | 0 (0/0) |
| 20 | arsenic trioxide | -4.582091 | 0.009532 | Drugs@FDA Prescription | YES | verify | verify | verify | 0 (0/0) |

Ties share a score; `rank_min` is in the JSON. 'verify' = needs manual verification; 'not_labelled' = label does not establish use in the band; 'age_not_stated' = established use without a minimum age. Unknown age is never assumed.

## 4. Label quotes (verbatim, truncated where marked)

### 1. niraparib (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id b7f675e2-159c-490c-b6f4-3f16d9492b7d, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of ZEJULA have not been established in pediatric patients."
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] (5.1) •Bone Marrow Suppression: Test complete blood counts weekly for the first month, monthly for the next 11 months, and periodically thereafter for clinically significant changes.
  * [toxicidad_medular, warnings] 5.2 Bone Marrow Suppression Hematologic adverse reactions, including thrombocytopenia, anemia, neutropenia, and/or pancytopenia have been reported in patients treated with ZEJULA [see Adverse Reactions (6)].
  * [segundos_tumores, warnings] 5 WARNINGS AND PRECAUTIONS •Myelodysplastic Syndrome/Acute Myeloid Leukemia (MDS/AML): MDS/AML occurred in patients exposed to ZEJULA, and some cases were fatal.
  * [segundos_tumores, warnings] (5.5, 8.1, 8.3) 5.1 Myelodysplastic Syndrome/Acute Myeloid Leukemia Myelodysplastic syndrome/acute myeloid leukemia (MDS/AML), including cases with a fatal outcome, have been reported in patients who received ZEJULA.

### 2. tamoxifen (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 7ee3d3d2-85d1-4018-8e70-5ed8a64ae1f0, section pediatric_use): "Pediatric Use The safety and efficacy of tamoxifen for girls aged 2 to 10 years with McCune-Albright Syndrome and precocious puberty have not been studied beyond one year of treatment."

### 3. olaparib (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (openFDA, set id 741ff3e3-dc1a-45a6-84e5-2481b27131aa, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness of Lynparza have not been established in pediatric patients."
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] Monitor complete blood count for cytopenia at baseline and monthly thereafter for clinically significant changes during treatment.
  * [toxicidad_medular, warnings] If the levels have not recovered to Grade 1 or less after 4 weeks, refer the patient to a hematologist for further investigations, including bone marrow analysis and blood sample for cytogenetics.
  * [segundos_tumores, warnings] 5 WARNINGS AND PRECAUTIONS • Myelodysplastic Syndrome/Acute Myeloid Leukemia (MDS/AML): Occurred in approximately 1.2% of patients with various BRCA m, g BRCA m, HRR gene-mutated or HRD-positive cancers exposed to Lynparza and the majority of events had a fatal outcome.
  * [segundos_tumores, warnings] ( 5.5 , 8.1 , 8.3 ) 5.1 Myelodysplastic Syndrome/Acute Myeloid Leukemia Myelodysplastic syndrome (MDS)/Acute Myeloid Leukemia (AML) has occurred in patients treated with Lynparza and some cases were fatal.

### 4. carboplatin (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 4c89cedc-e48b-445c-ad05-54dbe60e4fb0, section pediatric_use): "Pediatric Use Safety and effectiveness in pediatric patients have not been established (seeWARNINGS : "audiologic toxicity")."
* Paediatric warning: Pediatric Use Safety and effectiveness in pediatric patients have not been established (seeWARNINGS : "audiologic toxicity").
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNINGS Carboplatin injection should be administered under the supervision of a qualified physician experienced in the use of cancer chemotherapeutic agents.
  * [boxed_warning, boxed_warning] Appropriate management of therapy and complications is possible only when adequate treatment facilities are readily available.
  * [toxicidad_medular, boxed_warning] Bone marrow suppression is dose related and may be severe, resulting in infection and/or bleeding.
  * [toxicidad_medular, warnings] WARNINGS Bone marrow suppression (leukopenia, neutropenia, and thrombocytopenia) is dose-dependent and is also the dose-limiting toxicity.

### 5. danazol (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id e19acee5-ff33-45a5-bbbf-801756bc59b4, section pediatric_use): "Pediatric Use Safety and effectiveness in pediatric patients have not been established."

### 6. talazoparib (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 13839d4f-6acf-4ffb-a128-79df59319273, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of TALZENNA have not been established in pediatric patients."
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] (5.1) • Myelosuppression: TALZENNA may affect hematopoiesis and can cause anemia, neutropenia, and/or thrombocytopenia.
  * [toxicidad_medular, warnings] If counts do not recover within 4 weeks, refer the patient to a hematologist for further investigations including bone marrow analysis and blood sample for cytogenetics.
  * [segundos_tumores, warnings] 5 WARNINGS AND PRECAUTIONS • Myelodysplastic Syndrome/Acute Myeloid Leukemia (MDS/AML): MDS/AML occurred in patients exposed to TALZENNA, and some cases were fatal.
  * [segundos_tumores, warnings] (5.3, 8.1, 8.3) 5.1 Myelodysplastic Syndrome/Acute Myeloid Leukemia Myelodysplastic Syndrome/Acute Myeloid Leukemia (MDS/AML), including cases with a fatal outcome, has been reported in patients who received TALZENNA.

### 7. hydroxyurea (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id e76fd60e-7644-48c5-9857-3608a045000b, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness in pediatric patients have not been established."
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] 5 WARNINGS AND PRECAUTIONS •Myelosuppression: Do not give if bone marrow function is markedly depressed.
  * [toxicidad_medular, warnings] (5.8) 5.1 Myelosuppression Hydroxyurea causes severe myelosuppression.
  * [segundos_tumores, warnings] (5.2) •Malignancies: Advise protection from sun exposure and monitor for secondary malignancies.
  * [segundos_tumores, warnings] In patients receiving long-term hydroxyurea for myeloproliferative disorders, secondary leukemia has been reported.

### 8. anastrozole (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 362034e9-2830-19a5-e063-6394a90a8dfb, section pediatric_use): "8.4 Pediatric Use Clinical studies in pediatric patients included a placebo-controlled trial in pubertal boys of adolescent age with gynecomastia and a single-arm trial in girls with McCune-Albright Syndrome and progressive precocious puberty."

### 9. rucaparib (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 0295d202-1cfe-7659-e063-6294a90a476e, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of Rubraca in pediatric patients have not been established."
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] Monitor complete blood counts for cytopenia at baseline and monthly thereafter for clinically significant changes during treatment.
  * [toxicidad_medular, warnings] If the levels have not recovered to Grade 1 or less after 4 weeks or if MDS/AML is suspected, refer the patient to a hematologist for further investigations, including bone marrow analysis and blood sample for cytogenetics.
  * [segundos_tumores, warnings] 5 WARNINGS AND PRECAUTIONS Myelodysplastic Syndrome/Acute Myeloid Leukemia (MDS/AML): MDS/AML occurred in patients exposed to RUBRACA, and some cases were fatal.
  * [segundos_tumores, warnings] ( 5.2, 8.1, 8.3) 5.1 Myelodysplastic Syndrome/Acute Myeloid Leukemia Myelodysplastic Syndrome (MDS)/Acute Myeloid Leukemia (AML) occur in patients treated with RUBRACA, and are potentially fatal adverse reactions.

### 10. mitomycin (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 5917f7a7-cd20-407b-8dda-9fd58c5f47a6, section pediatric_use): "Pediatric Use Safety and effectiveness in pedriatric patients have not been established."
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNINGS Mitomycin should be administered under the supervision of a qualified physician experienced in the use of cancer chemotherapeutic agents.
  * [boxed_warning, boxed_warning] Appropriate management of therapy and complications is possible only when adequate diagnostic and treatment facilities are readily available.
  * [toxicidad_medular, boxed_warning] Bone marrow suppression, notably thrombocytopenia and leukopenia, which may contribute to overwhelming infections in an already compromised patient, is the most common and severe of the toxic effects of mitomycin (see “WARNINGS” and “ADVERSE REACTIONS” Sections).
  * [toxicidad_medular, boxed_warning] Hemolytic Uremic Syndrome (HUS) a serious complication of chemotherapy, consisting primarily of microangiopathic hemolytic anemia, thrombocytopenia, and irreversible renal failure, has been reported in patients receiving systemic mitomycin.

### 11. cisplatin (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 508496cb-3441-46b3-a4fe-e0d440e6adc6, section pediatric_use): "8.4 Pediatric Use Ototoxic effects may be more severe and detrimental in pediatric patients receiving cisplatin for injection, particularly in patients less than 5 years of age."
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNING: NEPHROTOXICITY, PERIPHERAL NEUROPATHY, NAUSEA AND VOMITING and MYELOSUPPRESSION. • Nephrotoxicity: cisplatin for injection can cause severe renal toxicity, including acute renal failure.
  * [boxed_warning, boxed_warning] Severe renal toxicities are dose-related and cumulative.
  * [toxicidad_medular, boxed_warning] WARNING: NEPHROTOXICITY, PERIPHERAL NEUROPATHY, NAUSEA AND VOMITING and MYELOSUPPRESSION. • Nephrotoxicity: cisplatin for injection can cause severe renal toxicity, including acute renal failure.
  * [toxicidad_medular, boxed_warning] Use highly effective antiemetic premedication [see Dosage and Administration (2.1) and Warnings and Precautions (5.3)]. • Myelosuppression: cisplatin for injection can cause severe myelosuppression with fatalities due to infections.

### 12. topotecan (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id b6019088-cd94-43b4-b5ce-ba20b6b8021e, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness in pediatric patients have not been established."
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNING: MYELOSUPPRESSION Topotecan hydrochloride can cause severe myelosuppression.
  * [boxed_warning, boxed_warning] Administer first cycle only to patients with baseline neutrophil counts of greater than or equal to 1,500/mm3 and platelet counts greater than or equal to 100,000/mm3.
  * [toxicidad_medular, boxed_warning] WARNING: MYELOSUPPRESSION Topotecan hydrochloride can cause severe myelosuppression.
  * [toxicidad_medular, boxed_warning] WARNING: MYELOSUPPRESSION See full prescribing information for complete boxed warning.

### 13. paclitaxel (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 0c12c77c-5fd4-48dc-a318-e96892233992, section pediatric_use): "Pediatric Use The safety and effectiveness of paclitaxel in pediatric patients have not been established."
* Paediatric warning: There have been reports of central nervous system (CNS) toxicity (rarely associated with death) in a clinical trial in pediatric patients in which paclitaxel was infused intravenously over 3 hours at doses ranging from 350 mg/m2 to 420 mg/m2.
* Paediatric warning: The toxicity is most likely attributable to the high dose of the ethanol component of the paclitaxel vehicle given over a short infusion time.
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] WARNING Paclitaxel should be administered under the supervision of a physician experienced in the use of cancer chemotherapeutic agents.
  * [boxed_warning, boxed_warning] Appropriate management of complications is possible only when adequate diagnostic and treatment facilities are readily available.
  * [toxicidad_medular, boxed_warning] In order to monitor the occurrence of bone marrow suppression, primarily neutropenia, which may be severe and result in infection, it is recommended that frequent peripheral blood cell counts be performed on all patients receiving paclitaxel.
  * [toxicidad_medular, warnings] Bone marrow suppression (primarily neutropenia) is dose-dependent and is the dose-limiting toxicity.

### 14. fostamatinib (status NOT_ESTABLISHED; manual review (config/pediatric_label_review.json))
* Paediatric use (DailyMed, set id 21149cc3-049b-43e2-b141-c9499160556c, section pediatric_use): "8.4 Pediatric Use Safety and effectiveness in pediatric patients have not been established."
* Paediatric warning: In subchronic, chronic, and carcinogenicity studies of TAVALISSE, chondrodystrophy of the femoral head was seen in rodents.
* Paediatric warning: In a study in juvenile rabbits, growth plate dysplasia was observed in the proximal femur and femoro-tibial joint, and bone marrow cellularity was reduced in the femur and sternum.

### 15. irinotecan (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id b66dbd6d-6ccb-ce9d-e053-2995a90a4c3c, section pediatric_use): "8.4 Pediatric Use The effectiveness of irinotecan in pediatric patients has not been established."
* Paediatric warning: Grade 3 to 4 neutropenia was experienced by 54 (31.8%) patients.
* Paediatric warning: Neutropenia was complicated by fever in 15 (8.8%) patients.
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] BOXED WARNING WARNING: DIARRHEA and MYELOSUPRESSION Early and late forms of diarrhea can occur.
  * [boxed_warning, boxed_warning] Early diarrhea may be accompanied by cholinergic symptoms which may be prevented or ameliorated by atropine.
  * [toxicidad_medular, boxed_warning] Institute antibiotic therapy if patients develop ileus, fever, or severe neutropenia.
  * [toxicidad_medular, boxed_warning] Severe myelosuppression may occur [see Warnings and Precautions (5.2)] WARNING: DIARRHEA and MYELOSUPPRESSION Early and late forms of diarrhea can occur.

### 16. zinc chloride (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 2bd646d2-2a2a-4376-8da0-c7f08b0511ac, section pediatric_use): "Pediatric Use See DOSAGE AND ADMINISTRATION section."

### 17. gemcitabine (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 44bcc16c-341d-49e8-a1c6-95d084667aca, section pediatric_use): "8.4 Pediatric Use The safety and effectiveness of Gemcitabine for Injection have not been established in pediatric patients."
* Paediatric warning: Patients with M1 or M2 bone marrow on Day 28 who did not experience unacceptable toxicity were eligible to receive a maximum of one additional four-week course.
* Paediatric warning: Toxicities observed included myelosuppression, febrile neutropenia, increased serum transaminases, nausea, and rash/desquamation.
* Chemotherapy caution (label-documented hazards): 
  * [toxicidad_medular, warnings] (5.1) Myelosuppression: Monitor for myelosuppression prior to each cycle and reduce or withhold dose for severe myelosuppression.
  * [toxicidad_medular, warnings] (5.10) 5.1 Schedule-Dependent Toxicity In clinical trials evaluating the maximum tolerated dose of Gemcitabine for Injection, prolongation of the infusion time beyond 60 minutes or more frequent than weekly dosing resulted in an increased incidence of clinically significant hypotension, severe flu-like symptoms, myelos [...]
  * [embrio_fetal, warnings] (5.6) Embryo-Fetal Toxicity: Can cause fetal harm.
  * [embrio_fetal, warnings] 5.7 Embryo-Fetal Toxicity Based on animal data and its mechanism of action, Gemcitabine for Injection can cause fetal harm when administered to a pregnant woman.

### 18. estradiol (status NOT_ESTABLISHED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 9333d92f-2585-4b99-aaea-87c01c95125e, section pediatric_use): "Pediatric Use Safety and effectiveness in pediatric patients have not been established."
* Paediatric warning: Large and repeated doses of estrogen over an extended period of time have been shown to accelerate epiphyseal closure, resulting in short adult stature if treatment is initiated before the completion of physiologic puberty in normally developing children.
* Paediatric warning: In patients in whom bone growth is not complete, periodic monitoring of bone maturation and effects on epiphyseal centers is recommended.

### 19. liothyronine (status MIXED; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id 548d684a-86f5-4e57-9844-d8cf0a84e068, section pediatric_use): "8.4 Pediatric Use The initial dose of liothyronine sodium tablets varies with age and body weight."
* Paediatric warning: Overtreatment is associated with craniosynostosis in infants, may adversely affect the tempo of brain maturation, and may accelerate the bone age and result in premature epiphyseal closure and compromised adult stature [see Dosage and Administration (2.6) and Adverse Reactions (6)].
* Paediatric warning: Overtreatment may accelerate the bone age and result in premature epiphyseal closure and compromised adult stature.

### 20. arsenic trioxide (status UNKNOWN; heuristic classifier of labels.classify_pediatric (not hand-verified); manual review required)
* Paediatric use (DailyMed, set id c99b01b4-67ff-4793-8375-fe3a21c6c9be, section pediatric_use): "8.4 Pediatric Use The safety and efficacy of arsenic trioxide injection as a single agent for treatment of pediatric patients with relapsed or refractory APL is supported by the pivotal phase 2 study in 40 patients with relapsed or refractory APL."
* Chemotherapy caution (label-documented hazards): 
  * [boxed_warning, boxed_warning] BOXED WARNING Differentiation Syndrome: Patients with acute promyelocytic leukemia (APL) treated with arsenic trioxide injection have experienced differentiation syndrome, which may be life-threatening or fatal.
  * [boxed_warning, boxed_warning] Signs and symptoms may include unexplained fever, dyspnea, hypoxia, acute respiratory distress, pulmonary infiltrates, pleural or pericardial effusions, weight gain, peripheral edema, hypotension, renal insufficiency, hepatopathy, and multi-organ dysfunction, in the presence or absence of leukocytosis.
  * [segundos_tumores, warnings] Monitor patients for the development of second primary malignancies.
  * [segundos_tumores, warnings] Monitor patients for the development of second primary malignancies.

## 5. Mechanism class x age band

Per class: the five small-molecule FDA drugs with most papers co-mentioning the drug and a gene of the class.

### DNA interstrand-crosslink repair (genes ERCC4, FAAP100, FANCA, FANCB, FANCC, FANCD2, FANCE, FANCF, FANCG, FANCI, FANCL, FANCM, RAD51, RFWD3; 171 retrieved papers mention a class gene)

| drug | papers | benefit | human | engine rank (small molecule) | 0-4 y | 5-11 y | 12-17 y |
|---|---|---|---|---|---|---|---|
| olaparib | 8 | 2 | 3 | 3 | not_labelled | not_labelled | not_labelled |
| carboplatin | 3 | 1 | 2 | 4 | not_labelled | not_labelled | not_labelled |
| cisplatin | 3 | 1 | 2 | 11 | not_labelled | not_labelled | not_labelled |
| niraparib | 3 | 1 | 1 | 1 | not_labelled | not_labelled | not_labelled |
| mitomycin | 3 | 0 | 1 | 10 | not_labelled | not_labelled | not_labelled |

Band cells by state (drug names):

* 0-4 y: not_labelled: carboplatin, cisplatin, mitomycin, niraparib, olaparib
* 5-11 y: not_labelled: carboplatin, cisplatin, mitomycin, niraparib, olaparib
* 12-17 y: not_labelled: carboplatin, cisplatin, mitomycin, niraparib, olaparib

### DNA repair/replication (genes BRCA1, BRCA2, BRIP1, PALB2, RAD51C, SLX4, UBE2T, XRCC2; 262 retrieved papers mention a class gene)

| drug | papers | benefit | human | engine rank (small molecule) | 0-4 y | 5-11 y | 12-17 y |
|---|---|---|---|---|---|---|---|
| olaparib | 11 | 3 | 5 | 3 | not_labelled | not_labelled | not_labelled |
| cisplatin | 4 | 1 | 2 | 11 | not_labelled | not_labelled | not_labelled |
| tamoxifen | 3 | 2 | 2 | 2 | not_labelled | not_labelled | not_labelled |
| carboplatin | 3 | 1 | 2 | 4 | not_labelled | not_labelled | not_labelled |
| niraparib | 3 | 1 | 1 | 1 | not_labelled | not_labelled | not_labelled |

Band cells by state (drug names):

* 0-4 y: not_labelled: carboplatin, cisplatin, niraparib, olaparib, tamoxifen
* 5-11 y: not_labelled: carboplatin, cisplatin, niraparib, olaparib, tamoxifen
* 12-17 y: not_labelled: carboplatin, cisplatin, niraparib, olaparib, tamoxifen

### spindle-assembly checkpoint (genes MAD2L2; 4 retrieved papers mention a class gene)

| drug | papers | benefit | human | engine rank (small molecule) | 0-4 y | 5-11 y | 12-17 y |
|---|---|---|---|---|---|---|---|
| (no small-molecule FDA drug co-mentioned with a class gene) | | | | | | | |

Band cells by state (drug names):

* 0-4 y: none
* 5-11 y: none
* 12-17 y: none

## 7. Evidence gap (computed; no dosing; no recommendation)

Orphanet case-count tier: >5 (Orphanet prevalence records of ORPHA:84).

* **number of patients available for any study** (Orphanet tier >5). Design: Multi-centre single-arm study against a natural-history control, with class-level borrowing as an informative prior.
* **no label to measure the engine on this disease; only the class-level behaviour is known** (class pool dna_repair_replication: mean AUC 0.81901, tau2 0.007521, k=6). Design: Bayesian hierarchical model borrowing from the mechanism class pool; any study result of one disease of the class updates the class level (hypothesis for study design only).
* **age band 0-4 y without labelled support for 20 of the 20 top-ranked drugs** (states verify, not_labelled or age_not_stated). Design: paediatric extrapolation or an age-stratified label supplement is a sponsor/regulatory step; within a study, the N-of-1 or series design of the first gap applies with an independent safety committee.
* **age band 5-11 y without labelled support for 20 of the 20 top-ranked drugs** (states verify, not_labelled or age_not_stated). Design: paediatric extrapolation or an age-stratified label supplement is a sponsor/regulatory step; within a study, the N-of-1 or series design of the first gap applies with an independent safety committee.
* **age band 12-17 y without labelled support for 20 of the 20 top-ranked drugs** (states verify, not_labelled or age_not_stated). Design: paediatric extrapolation or an age-stratified label supplement is a sponsor/regulatory step; within a study, the N-of-1 or series design of the first gap applies with an independent safety committee.

## 8. Literature volume

* Papers retrieved: 300; with a gene symbol: 296; Europe PMC hit count: 31219
* Drugs co-mentioned at least once: 26; sum of paper counts over drugs: 50
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
