# Candidates by line (hypotheses for follow-up)

*hypotheses for follow-up; each line proposes its own candidates; convergence is information, not a requirement; no claim of efficacy, no recommendation, no dose.*

Protocol: `protocol/addendum_multi_evidence.md section 12`. Common filters: approved small molecules of the completed pool (224 drugs); hard exclusions (the three hazard drugs, the growth-hormone class, and drugs whose label warns about a condition of the child); cytotoxic agents are listed as context/hazard under the line that found them. 28 candidates the approval check could not confirm are listed separately: Alfacalcidol, Ataluren, Butabarbital, Carbenoxolone, Chlorprothixene, Copanlisib, Daprodustat, Deserpidine, Dexketoprofen, Didanosine, Dydrogesterone, Glucosamine, Haloprogin, Hexylcaine, Infigratinib, Masoprocol, Mestranol, Methysergide, Niclosamide, Phenacemide, Plicamycin, Racecadotril, Rilmenidine, Tacalcitol, Trichlormethiazide, Tridihexethyl, Trioxsalen, Vernakalant.

## Literature (research line)

Native set: 11 drugs; qualifying after the filters: 4 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Phenylbutanoic Acid | rule group constitutional; best help level 1; 2 help papers; net 1.0 | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | entry in config/proposal_safety.json; label read, no keyword match for the phenotype terms | 2, 3 | level 1 literature (2 help, 0 harm papers); contrary: E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.001 (rank 68) |
| Trametinib | rule group tumour_context; best help level 1; 2 help papers; net 0.75 | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | entry in config/proposal_safety.json; label read, no keyword match for the phenotype terms | 5 (Track 1 genes bearing: CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | level 1 literature (2 help, 0 harm papers); contrary: E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.418 (rank 2) |
| Metformin | rule group constitutional; best help level 1; 1 help papers; net -0.25 | not_labelled / partly / labelled (supported in 5-11, 12-17) | boxed warning in the label; entry in config/proposal_safety.json; label read, no keyword match for the phenotype terms | 2, 3, 4 (Track 1 genes bearing: IGBP1, hypothesis) | level 1 literature (1 help, 2 harm papers); contrary: LIT, E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.000 (rank 87) |
| Everolimus | rule group class_sibling; best help level 0; 2 help papers; net 0.5 | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | boxed warning in the label; entry in config/proposal_safety.json; label read, no keyword match for the phenotype terms | 3, 4 (Track 1 genes bearing: IGBP1, hypothesis) | level 0 literature (2 help, 0 harm papers); contrary: E, EP; line limit: AI-assigned labels, mostly level 0-1 from mice, yeast and cell lines; rules written after the literature was read | dose-response in BUB1B-hypomorphic or patient-derived cells, reading the consequence the drug targets (growth/mTOR signalling, senescence, aneuploid stress) | 0.004 (rank 50) |

Found by the line but not proposed (context/hazard or filtered): Sirolimus [hard_excluded: label warns about a condition of the child (HP:0003202 muscle atrophy: 'rhabdomyolysis' in the warnings section)]; Chloroquine [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Nicotinamide Mononucleotide [not_eligible: not in the drug universe]; Bortezomib [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Bosutinib [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Dasatinib [not_eligible: not in the drug universe]; Selumetinib [hard_excluded: label warns about a condition of the child (HP:0004322/HP:0001508/HP:0001518 growth: 'reduced growth' in the warnings section)].

## Engine (frozen B0)

Native set: 1293 drugs; qualifying after the filters: 10 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Olaparib | engine rank 1 of 1293, score -2.229 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | 5 (Track 1 genes bearing: CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | literature without a supportive finding (1 harm papers); contrary: LIT, D, E, EP; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.000 (rank 88) |
| Tyrosine | engine rank 2 of 1293, score -3.065 | verify / verify / verify (no band with label support (verify/verify/verify)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.129 (rank 13) |
| Progesterone | engine rank 4 of 1293, score -3.731 | verify / verify / verify (no band with label support (verify/verify/verify)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.013 (rank 34) |
| Zinc Chloride | engine rank 5 of 1293, score -4.555 | verify / verify / verify (no band with label support (verify/verify/verify)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.082 (rank 16) |
| Vorinostat | engine rank 7 of 1293, score -4.710 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: prior from literature co-mention and graph targets; cannot be checked where no drug is known | test the co-mention link directly: drug in BUB1B-depleted versus control cells with checkpoint and proliferation readouts | 0.006 (rank 43) |

Also qualifying (not shown): Arsenic Trioxide, Fostamatinib, Belinostat, Carvedilol, Apomorphine.

Found by the line but not proposed (context/hazard or filtered): Carboplatin [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Cabazitaxel [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Fingolimod [hard_excluded: label warns about a condition of the child (HP:0004322/HP:0001508/HP:0001518 growth: 'growth retardation' in the warnings section)]; Cisplatin [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Doxorubicin [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Paclitaxel [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Pravastatin [hard_excluded: label warns about a condition of the child (HP:0003202 muscle atrophy: 'myopathy' in the warnings section)]; Lovastatin [hard_excluded: label warns about a condition of the child (HP:0003202 muscle atrophy: 'myopathy' in the contraindications section)].

## A knowledge graph

Native set: 3 drugs; qualifying after the filters: 1 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: weak discrimination (excess over a degree baseline 0.022); generic-drug filter applied.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Pomalidomide | graph rank 43; fraction of random queries with the drug in their top 50 0.006 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: weak discrimination (excess over a degree baseline 0.022); generic-drug filter applied | dose-response in BUB1B-deficient versus isogenic control cells for the phenotype the graph links to the drug | 0.003 (rank 53) |

Found by the line but not proposed (context/hazard or filtered): Cyclosporine [not_eligible: not a small molecule (Protein)]; Disopyramide [hard_excluded: label warns about a condition of the child (HP:0004322/HP:0001508/HP:0001518 growth: 'decreased growth' in the warnings section)].

## B1 phenotype module

Native set: 0 drugs; qualifying after the filters: 0 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: did not beat a degree baseline; validation not met.

No drug meets the line's native criterion after the filters. This is a result, not a failure of the table.

Below the native threshold, information only (not candidates): Mavacamten (robust in 5 of 9 configurations (needs 7); candidate); Binimetinib (robust in 4 of 9 configurations (needs 7); hard_excluded); Parathyroid Hormone (robust in 4 of 9 configurations (needs 7); not_eligible); Selumetinib (robust in 4 of 9 configurations (needs 7); hard_excluded); Trametinib (robust in 4 of 9 configurations (needs 7); candidate).

## B2 network proximity (mitosis module)

Native set: 4 drugs; qualifying after the filters: 3 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: exploratory; support set of a few drugs, from the pathogenic gene only.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Ceritinib | robust in scheme W3 (pathogenic gene only), two-record links | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: exploratory; support set of a few drugs, from the pathogenic gene only | target engagement and viability or mitotic-fidelity readout in cells with BUB1B loss | 0.238 (rank 8) |
| Selinexor | robust in scheme W3 (pathogenic gene only), two-record links | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | entry in config/proposal_safety.json; label read, no keyword match for the phenotype terms | 5 (Track 1 genes bearing: CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | level 1 literature (1 help, 0 harm papers); contrary: none found; line limit: exploratory; support set of a few drugs, from the pathogenic gene only | target engagement and viability or mitotic-fidelity readout in cells with BUB1B loss | 0.416 (rank 3) |
| Tirbanibulin | robust in scheme W3 (pathogenic gene only), two-record links | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: D; line limit: exploratory; support set of a few drugs, from the pathogenic gene only | target engagement and viability or mitotic-fidelity readout in cells with BUB1B loss | 0.033 (rank 25) |

Found by the line but not proposed (context/hazard or filtered): Natalizumab [not_eligible: not a small molecule (Antibody)].

## C PRISM drug-level vulnerability

Native set: 9 drugs; qualifying after the filters: 3 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Lazertinib | rho 0.209, BH q 0.019, 872 cell lines | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage | isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated | 0.069 (rank 18) |
| Neratinib | rho 0.194, BH q 0.019, 531 cell lines | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage | isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated | 0.211 (rank 9) |
| Trametinib | rho 0.131, BH q 0.047, 550 cell lines | partly / labelled / labelled (supported in 0-4, 5-11, 12-17) | entry in config/proposal_safety.json; label read, no keyword match for the phenotype terms | 5 (Track 1 genes bearing: CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, hypothesis) | level 1 literature (2 help, 0 harm papers); contrary: E, EP; line limit: correlational, from cancer cell lines, not adjusted for proliferation or lineage | isogenic BUB1B-low versus control lines, controlling for proliferation, with the drug sensitivity of the PRISM screen repeated | 0.418 (rank 2) |

Found by the line but not proposed (context/hazard or filtered): Prednisolone [hard_excluded: label warns about a condition of the child (HP:0003202 muscle atrophy: 'myopathy' in the warnings section)]; Bosutinib [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Cobimetinib [hard_excluded: label warns about a condition of the child (HP:0003202 muscle atrophy: 'rhabdomyolysis' in the warnings section)]; Dasatinib Anhydrous [hard_excluded: label warns about a condition of the child (HP:0004322/HP:0001508/HP:0001518 growth: 'growth retardation' in the warnings section)]; Calcitriol [hard_excluded: label warns about a condition of the child (HP:0000121 nephrocalcinosis: 'hypercalcemia' in the contraindications section)]; Selumetinib [hard_excluded: label warns about a condition of the child (HP:0004322/HP:0001508/HP:0001518 growth: 'reduced growth' in the warnings section)].

## D signed direction

Native set: 335 drugs; qualifying after the filters: 5 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: exploratory; 4 of 6 controls recovered, two by construction.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Regorafenib | strengthens checkpoint, net signed effect 94 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E, EP; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.008 (rank 41) |
| Sunitinib | strengthens checkpoint, net signed effect 57 | not_labelled / not_labelled / not_labelled (no band with label support (not_labelled/not_labelled/not_labelled)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.001 (rank 76) |
| Nirogacestat | strengthens checkpoint, net signed effect 49 | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.001 (rank 72) |
| Ritlecitinib | strengthens checkpoint, net signed effect 39 | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.000 (rank 191) |
| Erdafitinib | strengthens checkpoint, net signed effect 34 | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: none found; line limit: exploratory; 4 of 6 controls recovered, two by construction | spindle-checkpoint readouts (mitotic arrest, BUBR1/MAD2 at kinetochores) after the drug in BUB1B-hypomorph cells | 0.001 (rank 69) |

Found by the line but not proposed (context/hazard or filtered): Bortezomib [hard_excluded: hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)]; Carfilzomib [not_eligible: not a small molecule (Protein)]; Ixazomib Citrate [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)].

## E signature reversal, BUB1B knockout signature

Native set: 0 drugs; qualifying after the filters: 0 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: knockout signature may not represent partial loss of function; consensus library, no cell-line-specific signatures.

No drug meets the line's native criterion after the filters. This is a result, not a failure of the table.

## E-P signature reversal, patient signature (GSE22206)

Native set: 18 drugs; qualifying after the filters: 5 (counted in order of the line's own score and stopped at 5, or at 10 for the engine, to limit label reads). Line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls.

| Drug | Native score | Age 0-4 / 5-11 / 12-17 (unknown age) | Safety flags | Consequences addressed | Evidence tier | Refuting experiment | SMAA P(top 5), cross-check |
|---|---|---|---|---|---|---|--:|
| Amphotericin B | reversal p 0.0010, FDR q 0.059, R 10 | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.005 (rank 44) |
| Auranofin | reversal p 0.0010, FDR q 0.059, R 8 | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | boxed warning in the label; label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.001 (rank 71) |
| Neratinib | reversal p 0.0010, FDR q 0.059, R 12 | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.211 (rank 9) |
| Nitazoxanide | reversal p 0.0010, FDR q 0.059, R 11 | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.003 (rank 51) |
| Pyrimethamine | reversal p 0.0010, FDR q 0.059, R 8 | not_assessed / not_assessed / not_assessed (no band with label support (not_assessed/not_assessed/not_assessed)) | label read, no keyword match for the phenotype terms | none identified | no verified literature entry; contrary: E; line limit: one adult patient lymphoblastoid cell line, 2 arrays against 2 controls | RNA-seq before and after the drug in patient-derived cells, testing whether the patient signature reverses | 0.011 (rank 35) |

Found by the line but not proposed (context/hazard or filtered): Dexketoprofen [approval_unconfirmed: the approval check could not confirm it (may be approved under another name)]; Fluvastatin [hard_excluded: label warns about a condition of the child (HP:0003202 muscle atrophy: 'myopathy' in the warnings section)]; Glucosamine [approval_unconfirmed: the approval check could not confirm it (may be approved under another name)]; Ixazomib Citrate [context_hazard: cytotoxic or caution-list class: not proposed for constitutional use (project caution)]; Masoprocol [approval_unconfirmed: the approval check could not confirm it (may be approved under another name)]; Nelfinavir [hard_excluded: label warns about a condition of the child (HP:0003202 muscle atrophy: 'myopathy' in the contraindications section)]; Niclosamide [approval_unconfirmed: the approval check could not confirm it (may be approved under another name)]; Pimozide [hard_excluded: label warns about a condition of the child (HP:0003202 muscle atrophy: 'rhabdomyolysis' in the warnings section)]; Plicamycin [approval_unconfirmed: the approval check could not confirm it (may be approved under another name)].

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
| Erdafitinib |  |  |  |  |  |  | x |  |  |  | 1 |
| Everolimus | x |  |  |  |  |  |  |  |  |  | 1 |
| Lazertinib |  |  |  |  |  | x |  |  |  |  | 1 |
| Metformin | x |  |  |  |  |  |  |  |  |  | 1 |
| Nirogacestat |  |  |  |  |  |  | x |  |  |  | 1 |
| Nitazoxanide |  |  |  |  |  |  |  |  | x |  | 1 |
| Olaparib |  | x |  |  |  |  |  |  |  |  | 1 |
| Phenylbutanoic Acid | x |  |  |  |  |  |  |  |  |  | 1 |
| Pomalidomide |  |  | x |  |  |  |  |  |  |  | 1 |
| Progesterone |  | x |  |  |  |  |  |  |  |  | 1 |
| Pyrimethamine |  |  |  |  |  |  |  |  | x |  | 1 |
| Regorafenib |  |  |  |  |  |  | x |  |  |  | 1 |
| Ritlecitinib |  |  |  |  |  |  | x |  |  |  | 1 |
| Selinexor |  |  |  |  | x |  |  |  |  |  | 1 |
| Sunitinib |  |  |  |  |  |  | x |  |  |  | 1 |
| Tirbanibulin |  |  |  |  | x |  |  |  |  |  | 1 |
| Tyrosine |  | x |  |  |  |  |  |  |  |  | 1 |
| Vorinostat |  | x |  |  |  |  |  |  |  |  | 1 |
| Zinc Chloride |  | x |  |  |  |  |  |  |  |  | 1 |

Distinct candidates: 24; unique to one line: 22; in two or more lines: Neratinib, Trametinib.

Overlap between lines against chance (hypergeometric, universe = completed pool of 224 drugs; pairs with candidates in both lines; no multiplicity claim):

| Lines | Sizes | Overlap | Shared | p |
|---|---|--:|---|--:|
| LIT - ENG | 4, 5 | 0 | - | 1 |
| LIT - A | 4, 1 | 0 | - | 1 |
| LIT - B2 | 4, 3 | 0 | - | 1 |
| LIT - C | 4, 3 | 1 | Trametinib | 0.0529 |
| LIT - D | 4, 5 | 0 | - | 1 |
| LIT - EP | 4, 5 | 0 | - | 1 |
| ENG - A | 5, 1 | 0 | - | 1 |
| ENG - B2 | 5, 3 | 0 | - | 1 |
| ENG - C | 5, 3 | 0 | - | 1 |
| ENG - D | 5, 5 | 0 | - | 1 |
| ENG - EP | 5, 5 | 0 | - | 1 |
| A - B2 | 1, 3 | 0 | - | 1 |
| A - C | 1, 3 | 0 | - | 1 |
| A - D | 1, 5 | 0 | - | 1 |
| A - EP | 1, 5 | 0 | - | 1 |
| B2 - C | 3, 3 | 0 | - | 1 |
| B2 - D | 3, 5 | 0 | - | 1 |
| B2 - EP | 3, 5 | 0 | - | 1 |
| C - D | 3, 5 | 0 | - | 1 |
| C - EP | 3, 5 | 1 | Neratinib | 0.0658 |
| D - EP | 5, 5 | 0 | - | 1 |

## Age-band view (label support, a filter and a column, not evidence of efficacy)

| Drug | Lines | 0-4 | 5-11 | 12-17 | Unknown age (union) |
|---|---|---|---|---|---|
| Amphotericin B | EP | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Auranofin | EP | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Ceritinib | B2 | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Erdafitinib | D | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Everolimus | LIT | partly | labelled | labelled | supported in 0-4, 5-11, 12-17 |
| Lazertinib | C | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Metformin | LIT | not_labelled | partly | labelled | supported in 5-11, 12-17 |
| Neratinib | C, EP | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Nirogacestat | D | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Nitazoxanide | EP | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Olaparib | ENG | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Phenylbutanoic Acid | LIT | partly | labelled | labelled | supported in 0-4, 5-11, 12-17 |
| Pomalidomide | A | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Progesterone | ENG | verify | verify | verify | no band with label support (verify/verify/verify) |
| Pyrimethamine | EP | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Regorafenib | D | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Ritlecitinib | D | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Selinexor | B2 | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Sunitinib | D | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Tirbanibulin | B2 | not_assessed | not_assessed | not_assessed | no band with label support (not_assessed/not_assessed/not_assessed) |
| Trametinib | LIT, C | partly | labelled | labelled | supported in 0-4, 5-11, 12-17 |
| Tyrosine | ENG | verify | verify | verify | no band with label support (verify/verify/verify) |
| Vorinostat | ENG | not_labelled | not_labelled | not_labelled | no band with label support (not_labelled/not_labelled/not_labelled) |
| Zinc Chloride | ENG | verify | verify | verify | no band with label support (verify/verify/verify) |

## Cross-check with the multi-evidence SMAA (completed run)

The SMAA rule (P(top 5) >= 0.5) proposed: Calcitriol; trametinib survives: False. The probability of each candidate is in its table row and is not a criterion here. Stability of the SMAA run (completed): 
- label errors 10/20/30 %: baseline candidates retained 10 % 1.00 of 1, 20 % 1.00 of 1, 30 % 1.00 of 1 (trivial when the SMAA candidate has no literature labels).
- equal weights, top 5: Calcitriol, Dexmedetomidine, Trametinib, Selinexor, Acyclovir; literature-only weights, proposed: Eribulin, Trametinib, Selinexor, Vincristine, Temozolomide; gene block off, proposed: Calcitriol.

## Limits

- Each line has its own validity and limits (shown with each table). The lines differ in what they measure, so different candidates are expected; their overlap is reported against chance, not required.
- Candidates are hypotheses for follow-up: no claim of efficacy, no recommendation, no dose. Label keyword exclusions are a screen on label text, not a clinical judgement; a label that could not be read is flagged and not excluded.
- Track 1 genes other than BUB1B are hypotheses to test (karyotype or single-cell sequencing for the megabase events; function tests for IGBP1), not findings.

- Data gap: lines with an empty native set are reported as empty; the B1 near-threshold drugs are information only
