# Candidates by line (hypotheses for follow-up), version 2

*hypotheses for follow-up; each line proposes its own candidates; convergence is information, not a requirement; no claim of efficacy, no recommendation, no dose.*

Protocol: `protocol/addendum_multi_evidence.md sections 12 and 13 (version 2: corrected safety screen and age readings)`. Common filters: approved small molecules of the completed pool (224 drugs); hard exclusions (the three hazard drugs, the growth-hormone class, drugs the project's layers exclude, and documented child-specific contraindications quoted from the label; broad label keywords only flag for manual verification); cytotoxic agents are listed as context/hazard under the line that found them. 28 candidates the approval check could not confirm are listed separately: Alfacalcidol, Ataluren, Butabarbital, Carbenoxolone, Chlorprothixene, Copanlisib, Daprodustat, Deserpidine, Dexketoprofen, Didanosine, Dydrogesterone, Glucosamine, Haloprogin, Hexylcaine, Infigratinib, Masoprocol, Mestranol, Methysergide, Niclosamide, Phenacemide, Plicamycin, Racecadotril, Rilmenidine, Tacalcitol, Trichlormethiazide, Tridihexethyl, Trioxsalen, Vernakalant.

## What changed relative to the first output (`candidates_by_line_mva.md`, kept)

Corrections (addendum section 13): the keyword exclusion of the label screen is replaced by the project's hand-read layers as authority (keywords only flag, with the matched sentence), and every candidate has an age-band reading from the project's label functions.

| Line | First output | Version 2 | Added | Removed | Age readings changed |
|---|---|---|---|---|---|
| LIT | Phenylbutanoic Acid, Trametinib, Metformin, Everolimus | Sirolimus, Phenylbutanoic Acid, Trametinib, Metformin, Everolimus | Sirolimus | none | PHENYLBUTANOIC ACID |
| ENG | Olaparib, Tyrosine, Progesterone, Zinc Chloride, Vorinostat | Olaparib, Tyrosine, Progesterone, Zinc Chloride, Vorinostat | none | none | none |
| A | Pomalidomide | Pomalidomide, Disopyramide | Disopyramide | none | none |
| B1 | none | none | none | none | none |
| B2 | Ceritinib, Selinexor, Tirbanibulin | Ceritinib, Selinexor, Tirbanibulin | none | none | TIRBANIBULIN |
| C | Lazertinib, Neratinib, Trametinib | Lazertinib, Neratinib, Prednisolone, Trametinib, Cobimetinib | Prednisolone, Cobimetinib | none | LAZERTINIB, NERATINIB |
| D | Regorafenib, Sunitinib, Nirogacestat, Ritlecitinib, Erdafitinib | Regorafenib, Sunitinib, Nirogacestat, Ritlecitinib, Erdafitinib | none | none | NIROGACESTAT, RITLECITINIB, ERDAFITINIB |
| E | none | none | none | none | none |
| EP | Amphotericin B, Auranofin, Neratinib, Nitazoxanide, Pyrimethamine | Amphotericin B, Auranofin, Fluvastatin, Nelfinavir, Neratinib | Fluvastatin, Nelfinavir | Nitazoxanide, Pyrimethamine | AMPHOTERICIN B, AURANOFIN, NERATINIB |
| ET1 | none | none | none | none | none |

Named drugs under the corrected screen: Sirolimus: screen result candidate; candidate in LIT; also qualifying in no line; Selumetinib: screen result candidate; candidate in no line; also qualifying in LIT; Calcitriol: screen result hard_excluded (documented child-specific contraindication tied to HP:0000121 nephrocalcinosis ('hypercalcemia' in the contraindications section), for human review. Quote: "CONTRAINDICATIONS Calcitriol capsules should not be given to patients wit); candidate in no line; also qualifying in no line.

## Literature (research line)

Native set: 11 drugs; qualifying after the filters: 6 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Sirolimus | rule group constitutional; best help level 1; 7 help papers; net -0.5 | not_labelled / not_labelled / partly (supported in 12-17) | boxed warning in the label; project hand-read label warnings (config/proposal_safety.json); [boxed warning: boxed warning: immunosuppression, lymphoma and other malignancies - "Increased susceptibility to infection and the possible development of lymphoma and other malignancies may result from immunosuppression"]; [warnings 5.18: skin cancer - "Patients on immunosuppressive therapy are at increased risk for skin cancer."] | 2, 3, 4, 5 (Track 1 genes bearing: IGBP1/CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | level 1 literature (7 help, 9 harm papers); contrary: LIT, E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.000 (rank 197) |
| Phenylbutanoic Acid | rule group constitutional; best help level 1; 2 help papers; net 1.0 | labelled / labelled / labelled (supported in 0-4, 5-11, 12-17) | project hand-read label warnings (config/proposal_safety.json); [precautions: phenylacetate neurotoxicity in rat pups - "When given subcutaneously to rat pups, 190–474 mg/kg phenylacetate caused decreased proliferation and increased loss of neurons, and it reduced CNS myelin."]; [warnings: severe heart or kidney failure; sodium retention - "BUPHENYL should be used with great care, if at all, in patients with congestive heart failure or severe renal insufficiency, and in clinical states in which the"] | 2, 3 | level 1 literature (2 help, 0 harm papers); contrary: E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.001 (rank 68) |
| Trametinib | rule group tumour_context; best help level 1; 2 help papers; net 0.75 | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | project hand-read label warnings (config/proposal_safety.json); [warnings 5.1: new malignancies (with dabrafenib) - "In the pooled safety population of MEKINIST administered with dabrafenib, non-cutaneous malignancies occurred in 1% of patients."]; [warnings 5.5: cardiomyopathy - "Cardiomyopathy, including cardiac failure, can occur with MEKINIST."] | 5 (Track 1 genes bearing: CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | level 1 literature (2 help, 0 harm papers); contrary: E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.418 (rank 2) |
| Metformin | rule group constitutional; best help level 1; 1 help papers; net -0.25 | not_labelled / partly / labelled (supported in 5-11, 12-17) | boxed warning in the label; project hand-read label warnings (config/proposal_safety.json); [boxed warning: boxed warning: lactic acidosis - "Postmarketing cases of metformin-associated lactic acidosis have resulted in death, hypothermia, hypotension, and resistant bradyarrhythmias."]; [warnings 5.1: higher risk with renal impairment - "The risk of metformin accumulation and metformin-associated lactic acidosis increases with the severity of renal impairment because metformin is substantially e"] | 2, 3, 4 (Track 1 genes bearing: IGBP1, hypothesis) | level 1 literature (1 help, 2 harm papers); contrary: LIT, E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.000 (rank 87) |
| Everolimus | rule group class_sibling; best help level 0; 2 help papers; net 0.5 | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | boxed warning in the label; project hand-read label warnings (config/proposal_safety.json); [boxed warning: boxed warning (Zortress): immunosuppression, lymphoma and skin cancer - "Increased susceptibility to infection and the possible development of malignancies, such as lymphoma and skin cancer, may result from immunosuppression"]; [warnings, infections: serious infections more frequent below 6 years - "The incidence of serious infections was reported at a higher frequency in patients < 6 years of age"] | 3, 4 (Track 1 genes bearing: IGBP1, hypothesis) | level 0 literature (2 help, 0 harm papers); contrary: E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.004 (rank 50) |

Also qualifying (not shown): Selumetinib.

Found by the line but not proposed (context/hazard or filtered): Chloroquine [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Nicotinamide Mononucleotide [not_eligible: not in the drug universe]; Bortezomib [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Bosutinib [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Dasatinib [not_eligible: not in the drug universe].

## Engine (frozen B0)

Native set: 1293 drugs; qualifying after the filters: 10 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Olaparib | engine rank 1 of 1293, score -2.229 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | 5 (Track 1 genes bearing: CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | literature without a supportive finding (1 harm papers); contrary: LIT, D, E, EP; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.000 (rank 88) |
| Tyrosine | engine rank 2 of 1293, score -3.065 | verify / verify / verify (no band with label support (verify/verify/verify)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.129 (rank 13) |
| Progesterone | engine rank 4 of 1293, score -3.731 | verify / verify / verify (no band with label support (verify/verify/verify)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.013 (rank 34) |
| Zinc Chloride | engine rank 5 of 1293, score -4.555 | verify / verify / verify (no band with label support (verify/verify/verify)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.082 (rank 16) |
| Vorinostat | engine rank 7 of 1293, score -4.710 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.006 (rank 43) |

Also qualifying (not shown): Fingolimod, Arsenic Trioxide, Fostamatinib, Belinostat, Carvedilol.

Found by the line but not proposed (context/hazard or filtered): Carboplatin [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Cabazitaxel [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Cisplatin [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Doxorubicin [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Paclitaxel [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)].

## A knowledge graph

Native set: 3 drugs; qualifying after the filters: 2 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: weak discrimination (excess over a degree baseline 0.022); generic-drug filter applied.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Pomalidomide | graph rank 43; fraction of random queries with the drug in their top 50 0.006 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: weak discrimination (excess over a degree baseline 0.022); generic-drug filter applied | dose-response in BUB1B-deficient versus isogenic control cells for the phenotype the graph links to the drug | 0.003 (rank 53) |
| Disopyramide | graph rank 49; fraction of random queries with the drug in their top 50 0.036 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | boxed warning in the label; needs manual verification: HP:0004322/HP:0001508/HP:0001518 growth 'decreased growth' in the warnings section: "Pregnancy Teratogenic Effects Disopyramide phosphate was associated with decreased numbers of implantation sites and decreased growth and survival of pups when administered to pregnant rats at 250 mg/kg/day (20 or more times the usual daily human dose of 12 mg" | none identified | no verified literature entry; contrary: E, EP; line limit: weak discrimination (excess over a degree baseline 0.022); generic-drug filter applied | dose-response in BUB1B-deficient versus isogenic control cells for the phenotype the graph links to the drug | 0.001 (rank 77) |

Found by the line but not proposed (context/hazard or filtered): Cyclosporine [not_eligible: not a small molecule (Protein)].

## B1 phenotype module

Native set: 0 drugs; qualifying after the filters: 0 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: did not beat a degree baseline; validation not met.

No drug meets the line's native criterion after the filters. This is a result, not a failure of the table.

Below the native threshold, information only (not candidates): Mavacamten (robust in 5 of 9 configurations (needs 7); candidate); Binimetinib (robust in 4 of 9 configurations (needs 7); candidate); Parathyroid Hormone (robust in 4 of 9 configurations (needs 7); not_eligible); Selumetinib (robust in 4 of 9 configurations (needs 7); candidate); Trametinib (robust in 4 of 9 configurations (needs 7); candidate).

## B2 network proximity (mitosis module)

Native set: 4 drugs; qualifying after the filters: 3 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: exploratory; support set of a few drugs, from the pathogenic gene only.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Ceritinib | robust in scheme W3 (pathogenic gene only), two-record links | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: exploratory; support set of a few drugs, from the pathogenic gene only | target engagement and viability or mitotic-fidelity readout in cells with BUB1B loss | 0.238 (rank 8) |
| Selinexor | robust in scheme W3 (pathogenic gene only), two-record links | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | project hand-read label warnings (config/proposal_safety.json); [warnings 5.1: life-threatening thrombocytopenia - "XPOVIO can cause life-threatening thrombocytopenia, potentially leading to hemorrhage."]; [warnings: anorexia and weight loss - "Gastrointestinal Toxicity: Nausea, vomiting, diarrhea, anorexia, and weight loss may occur."] | 5 (Track 1 genes bearing: CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | level 1 literature (1 help, 0 harm papers); contrary: none found; line limit: exploratory; support set of a few drugs, from the pathogenic gene only | target engagement and viability or mitotic-fidelity readout in cells with BUB1B loss | 0.416 (rank 3) |
| Tirbanibulin | robust in scheme W3 (pathogenic gene only), two-record links | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: D; line limit: exploratory; support set of a few drugs, from the pathogenic gene only | target engagement and viability or mitotic-fidelity readout in cells with BUB1B loss | 0.033 (rank 25) |

Found by the line but not proposed (context/hazard or filtered): Natalizumab [not_eligible: not a small molecule (Antibody)].

## C PRISM drug-level vulnerability

Native set: 9 drugs; qualifying after the filters: 5 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Lazertinib | rho 0.209, BH q 0.019, 872 cell lines | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage | isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated | 0.069 (rank 18) |
| Neratinib | rho 0.194, BH q 0.019, 531 cell lines | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage | isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated | 0.211 (rank 9) |
| Prednisolone | rho 0.160, BH q 0.019, 531 cell lines | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | needs manual verification: HP:0003202 muscle atrophy 'myopathy' in the warnings section: "(See DOSAGE AND ADMINISTRATION.) An acute myopathy has been observed with the use of high doses of corticosteroids, most often occurring in patients with disorders of neuromuscular transmission (e.g., myasthenia gravis), or in patients receiving concomitant th" | none identified | no verified literature entry; contrary: D, E, EP; line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage | isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated | 0.001 (rank 79) |
| Trametinib | rho 0.131, BH q 0.047, 550 cell lines | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | project hand-read label warnings (config/proposal_safety.json); [warnings 5.1: new malignancies (with dabrafenib) - "In the pooled safety population of MEKINIST administered with dabrafenib, non-cutaneous malignancies occurred in 1% of patients."]; [warnings 5.5: cardiomyopathy - "Cardiomyopathy, including cardiac failure, can occur with MEKINIST."] | 5 (Track 1 genes bearing: CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | level 1 literature (2 help, 0 harm papers); contrary: E, EP; line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage | isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated | 0.418 (rank 2) |
| Cobimetinib | rho 0.129, BH q 0.053, 546 cell lines | verify / verify / verify (no band with label support (verify/verify/verify)) | needs manual verification: HP:0003202 muscle atrophy 'rhabdomyolysis' in the warnings section: "(5.6, 2.3) Rhabdomyolysis: Monitor creatine phosphokinase periodically and as clinically indicated for signs and symptoms of rhabdomyolysis." | none identified | no verified literature entry; contrary: E, EP; line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage | isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated | 0.048 (rank 19) |

Found by the line but not proposed (context/hazard or filtered): Bosutinib [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)].

## D signed direction

Native set: 335 drugs; qualifying after the filters: 5 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: exploratory; 4 of 6 controls recovered, two by construction.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Regorafenib | strengthens checkpoint, net signed effect 94 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.008 (rank 41) |
| Sunitinib | strengthens checkpoint, net signed effect 57 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.001 (rank 76) |
| Nirogacestat | strengthens checkpoint, net signed effect 49 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.001 (rank 72) |
| Ritlecitinib | strengthens checkpoint, net signed effect 39 | verify / verify / verify (no band with label support (verify/verify/verify)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.000 (rank 191) |
| Erdafitinib | strengthens checkpoint, net signed effect 34 | verify / verify / verify (no band with label support (verify/verify/verify)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.001 (rank 69) |

Found by the line but not proposed (context/hazard or filtered): Bortezomib [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Carfilzomib [not_eligible: not a small molecule (Protein)]; Ixazomib Citrate [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)].

## E signature reversal, BUB1B knockout signature

Native set: 0 drugs; qualifying after the filters: 0 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: knockout signature may not represent partial loss of function; consensus library, no cell-line-specific signatures.

No drug meets the line's native criterion after the filters. This is a result, not a failure of the table.

## E-P signature reversal, patient signature (GSE22206)

Native set: 18 drugs; qualifying after the filters: 5 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Amphotericin B | reversal p 0.0010, FDR q 0.059, R 10 | verify / verify / verify (no band with label support (verify/verify/verify)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.005 (rank 44) |
| Auranofin | reversal p 0.0010, FDR q 0.059, R 8 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.001 (rank 71) |
| Fluvastatin | reversal p 0.0010, FDR q 0.059, R 8 | verify / verify / verify (no band with label support (verify/verify/verify)) | needs manual verification: HP:0003202 muscle atrophy 'myopathy' in the warnings section: "5 WARNINGS AND PRECAUTIONS Skeletal muscle effects (e.g., myopathy and rhabdomyolysis): Risks increase with advanced age (≥ 65), uncontrolled hypothyroidism, renal impairment, and combination use with cyclosporine,or gemfibrozil." | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.018 (rank 30) |
| Nelfinavir | reversal p 0.0010, FDR q 0.059, R 8 | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | needs manual verification: HP:0003202 muscle atrophy 'myopathy' in the contraindications section: "HMG-CoA Reductase Inhibitors Lovastatin, Simvastatin Potential for serious reactions such as myopathy including rhabdomyolysis." | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.090 (rank 15) |
| Neratinib | reversal p 0.0010, FDR q 0.059, R 12 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.211 (rank 9) |

Found by the line but not proposed (context/hazard or filtered): Dexketoprofen [approval_unconfirmed: the approval check could not confirm it (may be approved under another name)]; Glucosamine [approval_unconfirmed: the approval check could not confirm it (may be approved under another name)]; Ixazomib Citrate [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Masoprocol [approval_unconfirmed: the approval check could not confirm it (may be approved under another name)].

## E-T1 signature reversal, Track 1 combined signature

Native set: 0 drugs; qualifying after the filters: 0 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: weighted combination of gene signatures; consensus library.

No drug meets the line's native criterion after the filters. This is a result, not a failure of the table.

## Convergence matrix (candidate shown by the line: x)

Convergence across lines is information and is never required.

| Drug | LIT | ENG | A | B1 | B2 | C | D | E | EP | ET1 | Lines |
|---|---|---|---|---|---|---|---|---|---|---|--:|
| Neratinib |  |  |  |  |  | x |  |  | x |  | 2 |
| Trametinib | x |  |  |  |  | x |  |  |  |  | 2 |
| Amphotericin B |  |  |  |  |  |  |  |  | x |  | 1 |
| Auranofin |  |  |  |  |  |  |  |  | x |  | 1 |
| Ceritinib |  |  |  |  | x |  |  |  |  |  | 1 |
| Cobimetinib |  |  |  |  |  | x |  |  |  |  | 1 |
| Disopyramide |  |  | x |  |  |  |  |  |  |  | 1 |
| Erdafitinib |  |  |  |  |  |  | x |  |  |  | 1 |
| Everolimus | x |  |  |  |  |  |  |  |  |  | 1 |
| Fluvastatin |  |  |  |  |  |  |  |  | x |  | 1 |
| Lazertinib |  |  |  |  |  | x |  |  |  |  | 1 |
| Metformin | x |  |  |  |  |  |  |  |  |  | 1 |
| Nelfinavir |  |  |  |  |  |  |  |  | x |  | 1 |
| Nirogacestat |  |  |  |  |  |  | x |  |  |  | 1 |
| Olaparib |  | x |  |  |  |  |  |  |  |  | 1 |
| Phenylbutanoic Acid | x |  |  |  |  |  |  |  |  |  | 1 |
| Pomalidomide |  |  | x |  |  |  |  |  |  |  | 1 |
| Prednisolone |  |  |  |  |  | x |  |  |  |  | 1 |
| Progesterone |  | x |  |  |  |  |  |  |  |  | 1 |
| Regorafenib |  |  |  |  |  |  | x |  |  |  | 1 |
| Ritlecitinib |  |  |  |  |  |  | x |  |  |  | 1 |
| Selinexor |  |  |  |  | x |  |  |  |  |  | 1 |
| Sirolimus | x |  |  |  |  |  |  |  |  |  | 1 |
| Sunitinib |  |  |  |  |  |  | x |  |  |  | 1 |
| Tirbanibulin |  |  |  |  | x |  |  |  |  |  | 1 |
| Tyrosine |  | x |  |  |  |  |  |  |  |  | 1 |
| Vorinostat |  | x |  |  |  |  |  |  |  |  | 1 |
| Zinc Chloride |  | x |  |  |  |  |  |  |  |  | 1 |

Distinct candidates: 28; unique to one line: 26; in two or more lines: Neratinib, Trametinib.

Overlap between lines against chance (hypergeometric, universe = completed pool of 224 drugs; pairs with candidates in both lines; no multiplicity claim):

| Lines | Sizes | Overlap | Shared | p |
|---|---|--:|---|--:|
| LIT - ENG | 5, 5 | 0 | - | 1 |
| LIT - A | 5, 2 | 0 | - | 1 |
| LIT - B2 | 5, 3 | 0 | - | 1 |
| LIT - C | 5, 5 | 1 | Trametinib | 0.108 |
| LIT - D | 5, 5 | 0 | - | 1 |
| LIT - EP | 5, 5 | 0 | - | 1 |
| ENG - A | 5, 2 | 0 | - | 1 |
| ENG - B2 | 5, 3 | 0 | - | 1 |
| ENG - C | 5, 5 | 0 | - | 1 |
| ENG - D | 5, 5 | 0 | - | 1 |
| ENG - EP | 5, 5 | 0 | - | 1 |
| A - B2 | 2, 3 | 0 | - | 1 |
| A - C | 2, 5 | 0 | - | 1 |
| A - D | 2, 5 | 0 | - | 1 |
| A - EP | 2, 5 | 0 | - | 1 |
| B2 - C | 3, 5 | 0 | - | 1 |
| B2 - D | 3, 5 | 0 | - | 1 |
| B2 - EP | 3, 5 | 0 | - | 1 |
| C - D | 5, 5 | 0 | - | 1 |
| C - EP | 5, 5 | 1 | Neratinib | 0.108 |
| D - EP | 5, 5 | 0 | - | 1 |

## Age-band view (label support, a filter and a column, not evidence of efficacy)

| Drug | Lines | 0-4 | 5-11 | 12-17 | Unknown age (union) |
|---|---|---|---|---|---|
| Amphotericin B | EP | verify | verify | verify | no band with label support (verify/verify/verify) |
| Auranofin | EP | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Ceritinib | B2 | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Cobimetinib | C | verify | verify | verify | no band with label support (verify/verify/verify) |
| Disopyramide | A | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Erdafitinib | D | verify | verify | verify | no band with label support (verify/verify/verify) |
| Everolimus | LIT | partly | labelled | labelled | supported in 0-4, 5-11, 12-17 |
| Fluvastatin | EP | verify | verify | verify | no band with label support (verify/verify/verify) |
| Lazertinib | C | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Metformin | LIT | not_labelled | partly | labelled | supported in 5-11, 12-17 |
| Nelfinavir | EP | partly | labelled | labelled | supported in 0-4, 5-11, 12-17 |
| Neratinib | C, EP | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Nirogacestat | D | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Olaparib | ENG | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Phenylbutanoic Acid | LIT | labelled | labelled | labelled | supported in 0-4, 5-11, 12-17 |
| Pomalidomide | A | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Prednisolone | C | partly | labelled | labelled | supported in 0-4, 5-11, 12-17 |
| Progesterone | ENG | verify | verify | verify | no band with label support (verify/verify/verify) |
| Regorafenib | D | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Ritlecitinib | D | verify | verify | verify | no band with label support (verify/verify/verify) |
| Selinexor | B2 | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Sirolimus | LIT | not_labelled | not_labelled | partly | supported in 12-17 |
| Sunitinib | D | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Tirbanibulin | B2 | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Trametinib | LIT, C | partly | labelled | labelled | supported in 0-4, 5-11, 12-17 |
| Tyrosine | ENG | verify | verify | verify | no band with label support (verify/verify/verify) |
| Vorinostat | ENG | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Zinc Chloride | ENG | verify | verify | verify | no band with label support (verify/verify/verify) |

Age readings come from `engine_case_study.label_support` (age_label_table.row_for, manual review where it exists, otherwise the heuristic classifier). Section used and quote:

- Sirolimus: status ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use Renal Transplant The safety and efficacy of sirolimus in pediatric patients <13 years have not been established."
- Phenylbutanoic Acid: status ESTABLISHED (manual review); section pediatric_use: "Pediatric Use The use of tablets for neonates, infants and children to the weight of 20 kg is not recommended."
- Trametinib: status ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use BRAF V600E Mutation-Positive Unresectable or Metastatic Solid Tumors and LGG The safety and effectiveness of MEKINIST in combination with dabrafenib have been established in pediatric patients 1 year of age and older with unresectable or meta"
- Metformin: status ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of metformin hydrochloride for the treatment of type 2 diabetes mellitus have been established in pediatric patients 10 to 16 years old."
- Everolimus: status ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use TSC-Associated SEGA The safety and effectiveness of AFINITOR/AFINITOR DISPERZ have been established in pediatric patients age 1 year and older with TSC-associated SEGA that requires therapeutic intervention but cannot be curatively resected."
- Olaparib: status NOT_ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use Safety and effectiveness of Lynparza have not been established in pediatric patients."
- Tyrosine: status UNKNOWN (heuristic classifier, not hand-verified); section none found: "no paediatric-use sentence found"
- Progesterone: status UNKNOWN (heuristic classifier, not hand-verified); section pediatric_use: "Pediatric Use Progesterone capsules are not indicated in children."
- Zinc Chloride: status UNKNOWN (heuristic classifier, not hand-verified); section pediatric_use: "Pediatric Use See DOSAGE AND ADMINISTRATION section."
- Vorinostat: status NOT_ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of ZOLINZA in pediatric patients have not been established."
- Pomalidomide: status NOT_ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of pomalidomide have not been established in pediatric patients."
- Disopyramide: status NOT_ESTABLISHED (manual review); section pediatric_use: "Pediatric Use Safety and effectiveness in pediatric patients have not been established (see Dosage and Administration )."
- Ceritinib: status NOT_ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of ZYKADIA in pediatric patients have not been established."
- Selinexor: status NOT_ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of XPOVIO have not been established in pediatric patients."
- Tirbanibulin: status NOT_ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of KLISYRI for actinic keratosis in subjects less than 18 years of age have not been established."
- Lazertinib: status NOT_ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of LAZCLUZE in pediatric patients have not been established."
- Neratinib: status NOT_ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety and efficacy of NERLYNX in pediatric patients has not been established."
- Prednisolone: status ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "PEDIATRIC USE The efficacy and safety of prednisolone in the pediatric population are based on the well-established course of effect of corticosteroids which is similar in pediatric and adult populations."
- Cobimetinib: status MIXED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of COTELLIC have not been established in pediatric patients."
- Regorafenib: status NOT_ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use The safety and efficacy of STIVARGA in pediatric patients less than 18 years of age have not been established."
- Sunitinib: status NOT_ESTABLISHED (manual review); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of sunitinib malate in pediatric patients have not been established."
- Nirogacestat: status NOT_ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of OGSIVEO have not been established in pediatric patients."
- Ritlecitinib: status MIXED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety and effectiveness of LITFULO for the treatment of alopecia areata have been established in pediatric patients ages 12 years and older."
- Erdafitinib: status MIXED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use Safety and effectiveness of BALVERSA in pediatric patients have not been established."
- Amphotericin B: status MIXED (heuristic classifier, not hand-verified); section pediatric_use: "Pediatric Use Pediatric patients, age 1 month to 16 years, with presumed fungal infection (empirical therapy), confirmed systemic fungal infections or with visceral leishmaniasis have been successfully treated with amphotericin B liposome for injection."
- Auranofin: status NOT_ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "Pediatric Use: Auranofin Capsules (auranofin) is not recommended for use in pediatric patients because its safety and effectiveness have not been established."
- Fluvastatin: status UNKNOWN (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety and efficacy of fluvastatin sodium in children and adolescent patients 9 to 16 years of age with heterozygous familial hypercholesterolemia have been evaluated in open-label, uncontrolled clinical trials for a duration of two years"
- Nelfinavir: status ESTABLISHED (heuristic classifier, not hand-verified); section pediatric_use: "8.4 Pediatric Use The safety, tolerability, pharmacokinetic profile and efficacy of VIRACEPT were evaluated in HIV infected pediatric patients from 2 to 13 years of age in multicenter clinical trials, Study 556 and PACTG 337 [see Clinical Pharmacology (12.3) a"

## Cross-check with the multi-evidence SMAA (completed run)

The SMAA rule (P(top 5) >= 0.5) proposed: Calcitriol; trametinib survives: False. The probability of each candidate is in its table row and is not a criterion here. Stability of the SMAA run (completed): 
- label errors 10/20/30 %: baseline candidates retained 10 % 1.00 of 1, 20 % 1.00 of 1, 30 % 1.00 of 1 (trivial when the SMAA candidate has no literature labels).
- equal weights, top 5: Calcitriol, Dexmedetomidine, Trametinib, Selinexor, Acyclovir; literature-only weights, proposed: Eribulin, Trametinib, Selinexor, Vincristine, Temozolomide; gene block off, proposed: Calcitriol.

## Limits

- Each line has its own validity and limits (shown with each table). The lines differ in what they measure, so different candidates are expected; their overlap is reported against chance, not required.
- Candidates are hypotheses for follow-up: no claim of efficacy, no recommendation, no dose. Label keyword exclusions are a screen on label text, not a clinical judgement; a label that could not be read is flagged and not excluded.
- Track 1 genes other than BUB1B are hypotheses to test (karyotype or single-cell sequencing for the megabase events; function tests for IGBP1), not findings.

- Data gap: lines with an empty native set are reported as empty; the B1 near-threshold drugs are information only
