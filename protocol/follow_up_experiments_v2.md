# Follow-up plan, version 2: computation first, laboratory only for what computation cannot refute

Status: design. Nothing in sections 3 and 4 has been run. In vitro only. No dose, schedule or treatment recommendation for any person. The case is described only by gene symbol (*BUB1B*) and the categorical profile: a stop-gain predicted to undergo decay, plus a missense variant of uncertain significance, phase unresolved, age unknown. Numbers come from `out/candidates_by_line_mva_v2.md`, `out/line_c_adjusted.md`, `out/line_e_signature_reversal.md`, `out/empirical/power_analysis.json`, `out/caution_list.json` and `config/multi_evidence_track1_gene_block.json`. Items marked DESIGN CHOICE are not derived from project files. This version replaces `protocol/follow_up_experiments.md`, which was written before candidates were organised by line; the old file is kept.

Rule fixed by the user: the multi-parameter analysis serves to propose experiments. A hypothesis that computation can refute is tested by computation; only the rest receives a laboratory design (`protocol/addendum_refutability.md`). "Not supported by computation" is not "refuted".

## 1. From variant to candidates

| Step | Content | Status |
|---|---|---|
| Variant | Biallelic *BUB1B*: stop-gain with predicted decay (loss of function) plus a conserved, buried kinase-domain missense of uncertain significance; phase unresolved | Biallelic state is a hypothesis; the missense effect is unproven |
| Track 1 gene block | BUB1B primary. IGBP1 (hemizygous missense, PP2A regulatory subunit, linked in the literature to AKT dephosphorylation) is a hypothesis-grade modifier of the mTOR/AKT axis. CENPE and MAD1L1: pathway context only. CDC20, MAD2L1BP, TTK, CEP192: megabase events not analysed in Track 1, real or artefact undetermined | Hypotheses to test, not findings |
| Mechanism | Five consequences of reduced BubR1: (1) mosaic aneuploidy, (2) proteotoxic and energy stress, (3) p53/p16 senescence, (4) growth restriction and muscle atrophy (mTORC1 signal in mice), (5) tumour predisposition | Literature, level 0-1 for drugs |
| Lines | LIT (literature), ENG (engine), A (graph), B1, B2 (mitosis module), C (cell-line vulnerability), D (signed direction), E, EP, ET1 (signature reversal). B1, E and ET1 returned no candidate after the filters | Each line has its own candidates |

Result: 28 distinct candidates (30 line entries; Neratinib and Trametinib in two lines), plus four antimitotic cytotoxics found by line C. The 28 coincides with the count of candidates whose approval could not be confirmed; these are different sets, and the unconfirmed ones are not in the triage.

## 2. Triage

Age codes are the label readings for 0-4 / 5-11 / 12-17: L labelled, P partly, N not labelled, V verify. They describe other diseases, not suitability. "Lab" = laboratory design in section 4; "Comp" = computational test in section 3 first. Status: S = supported by a line; U = unsupported; NE = not estimable; H = hazard-class.

| Candidate | Line(s) | Status after checks done | Age | Uncertainty that decides priority | Next step |
|---|---|---|---|---|---|
| Sirolimus | LIT | NE (E: R -6, q 0.925); literature level 1 (7 help, 9 harm papers) | N/N/P | Net literature balance negative (-0.5); boxed warning; IGBP1 link untested | Lab X1; Comp CF8 |
| Everolimus | LIT | NE (E: R -3, q 0.854); level 0 | P/L/L | Class sibling at level 0; must move with sirolimus | Lab X1 |
| Metformin | LIT | NE (E: R 26, q 0.674); level 1 (1 help, 2 harm) | N/P/L | Contradictory papers; AMPK versus mTOR attribution | Lab X1; Comp CF8 |
| Phenylbutanoic Acid | LIT | NE (E: R 5, q 0.777); level 1 (2 help, 0 harm) | L/L/L | Rat-pup neuro-toxicity statement on label | Lab X2 |
| Trametinib | LIT, C | U in C (adjusted q 0.3872); LIT level 1 (2 help) | P/L/L | Tumour-context only; failed adjustment; SMAA rank 2 is not a criterion | Lab X3 (tumour lines) |
| Lazertinib | C | S in C (adjusted rho 0.140, q 0.0316; excl. haem/lymphoid q 0.0380; 15/15 leave-one-lineage-out) | N/N/N | Sign of rho (viability) and mechanism unclear; no literature | Comp CF1, then Lab X4 |
| Neratinib | C, EP | S in C (adjusted rho 0.187, q 0.0316) and EP (q 0.059; one adult line, 2 case vs 4 control arrays) | N/N/N | EP signature from one adult line; sign of C rho | Comp CF1, CF7, then Lab X4 |
| Prednisolone | C | U (adjusted q 0.3227) | P/L/L | Lymphoid suspicion not testable (0 haematopoietic/lymphoid lines in its data) | Closed; no lab |
| Cobimetinib | C | U (adjusted q 0.5892) | V/V/V | Failed adjustment | Closed; no lab |
| Ceritinib | B2 | S in B2 (scheme W3 only, two-record links) | N/N/N | Robust in one scheme; support set of a few drugs | Comp CF5, then Lab X5 |
| Selinexor | B2 | S in B2 (W3 only); level 1 (1 help) | N/N/N | Thrombocytopenia warning; one scheme | Comp CF5, then Lab X5 |
| Tirbanibulin | B2 | S in B2 (W3 only) | N/N/N | Topical label; one scheme | Comp CF5, then Lab X5 |
| Pomalidomide | A | S in A, weak (rank 43; 0.006 of random queries; excess over degree baseline 0.022) | N/N/N | Weak line discrimination; boxed warning | Comp CF4 |
| Disopyramide | A | S in A, weak (rank 49; 0.036 of random queries) | N/N/N | Weak discrimination; growth wording in label | Comp CF4 |
| Olaparib | ENG | NE (engine prior; 1 harm paper) | N/N/N | Co-mention link unverified | Comp CF3 |
| Tyrosine | ENG | NE (no verified entry) | V/V/V | Co-mention only | Comp CF3 |
| Progesterone | ENG | NE (no verified entry) | V/V/V | Co-mention only; boxed warning | Comp CF3 |
| Zinc Chloride | ENG | NE (no verified entry) | V/V/V | Co-mention only | Comp CF3 |
| Vorinostat | ENG | NE (no verified entry) | N/N/N | Co-mention only | Comp CF3 |
| Regorafenib | D | S in D, exploratory (net signed effect 94) | N/N/N | Sign calls; 4 of 6 controls recovered, two by construction | Comp CF6 |
| Sunitinib | D | S in D (57) | N/N/N | Same | Comp CF6 |
| Nirogacestat | D | S in D (49) | N/N/N | Same | Comp CF6 |
| Ritlecitinib | D | S in D (39) | V/V/V | Same; mixed label | Comp CF6 |
| Erdafitinib | D | S in D (34); original C hit did not survive (adjusted q 0.8509) | V/V/V | Same; contradictory C signal | Comp CF6 |
| Amphotericin B | EP | S in EP (q 0.059; R 10) | V/V/V | One adult line, 2 case vs 4 control arrays; boxed warning | Comp CF7 |
| Auranofin | EP | S in EP (R 8) | N/N/N | Same | Comp CF7 |
| Fluvastatin | EP | S in EP (R 8); C original hit did not survive (adjusted q 0.1553) | V/V/V | Same; contradictory C signal | Comp CF7 |
| Nelfinavir | EP | S in EP (R 8) | P/L/L | Same | Comp CF7 |
| Colchicine | C | H: survives (adjusted rho -0.159, q 0.0316); tubulin class in caution list | not read | Score may track mitotic dependency, not selective vulnerability | Lab H1 (toxicity); Comp CF2 |
| Eribulin | C | H: survives (adjusted rho -0.115, q 0.0316) | not read | Same | Lab H1; Comp CF2 |
| Paclitaxel | C | H: survives (adjusted rho -0.141, q 0.0542) | not read | Same; also an engine context hazard | Lab H1; Comp CF2 |
| Vincristine | C | H: survives (adjusted rho -0.158, q 0.0316) | not read | Same | Lab H1; Comp CF2 |

32 rows: 10 with a laboratory design for benefit (LIT four, Trametinib, Lazertinib, Neratinib, Ceritinib, Selinexor, Tirbanibulin), 16 with a computational test first (ENG 5, A 2, D 5, EP 4), 2 closed, 4 hazard-class with toxicity readouts only. Bortezomib, chloroquine and bosutinib are hard-excluded hazards, not rows; they enter the hazard panel. Calcitriol (hard-excluded, documented contraindication) and Selumetinib (in no line) are outside the 28.

## 3. Computational follow-ups that remain (defined, not run)

| ID | Test | Refutation criterion (pre-specified) |
|---|---|---|
| CF1 | Sign audit of line C: state whether positive rho means higher viability (resistance) at a high checkpoint score; recompute with an orientation-explicit sensitivity metric | Lazertinib and Neratinib are dropped as "vulnerabilities" if the sign means resistance |
| CF2 | Class-enrichment test: are the four tubulin-class survivors more than expected among the 190 estimable drugs (hypergeometric)? | If yes, treat the score as a mitotic-dependency control, not a source of candidates; the four stay hazard-class |
| CF3 | Engine specificity: rank of each drug for BUB1B against degree-matched random genes | Drug not above the null: close the row |
| CF4 | Graph-path audit for line A: list the paths drug to disease node and check each edge for a primary source and for generic hub nodes | All paths pass via hubs or lack a source: close |
| CF5 | Scheme sensitivity for B2: recompute support under the other weighting schemes and against random genes | Not in the top set under any other scheme: close |
| CF6 | Line D: sign-permutation null and re-check of the control drugs | Net effect inside the null range: close |
| CF7 | Line E-P: leave-one-array-out (2 case vs 4 control arrays) and reversal in the BUB1B and Track 1 signatures | Reversal q at or above 0.10 in any leave-out, or concordance with the disease signature: close or flag hazard |
| CF8 | IGBP1 axis: network distance IGBP1 to MTOR/AMPK nodes against degree-matched random pairs | Not shorter than null: the IGBP1 rationale is dropped from X1, not the drugs |
| CF9 | Second drug-sensitivity screen for the line C survivors, only if one is already on disk | No same-sign association: close |

## 4. Laboratory designs, by information value

### 4.1 Step 0 diagnostic (E0), clinical-diagnostic setting, with consent

| Test | Purpose |
|---|---|
| Karyotype (metaphase spreads) of blood, fibroblasts if clinically justified; single-cell sequencing of the same material | Confirm and quantify mosaic aneuploidy; the read-depth screen is not a karyotype |
| Same tests for the megabase events of Track 1 (CDC20, MAD2L1BP, TTK, CEP192) | Confirm or discard them |
| Test of the missense variant: BubR1 level, turnover, kinetochore localisation, checkpoint arrest, wild type versus variant in BUB1B-depleted cells | It is of uncertain significance; phase by parental or long-read testing |

Stop rule S0: no excess aneuploidy above the laboratory's control range in every tissue examined, and the missense behaves as wild type or lies in cis with the stop-gain. Under S0 drug work on the child's cells stops; model work may continue as general research, labelled so.

### 4.2 Cell models (M1)

Patient-derived fibroblasts or lymphoblastoid lines only with consent and if a sample exists (low passage; karyotype at start and end, because cultures lose aneuploid clones); isogenic partial-BUB1B-reduction line (titrated, two independent clones, residual BubR1 confirmed); euploid parental line and a donor line; rescue line; tumour-context lines (rhabdomyosarcoma) beside the non-tumour lines. Scoring blinded, plates randomised, vehicle matched.

### 4.3 Experiments

Common readouts: micronuclei; mis-segregation by live imaging; aneuploid fraction by single-cell sequencing or FISH; DNA damage (gamma-H2AX, 53BP1); senescence (SA-beta-gal, p21); growth. In vitro concentrations: 0.1, 0.3, 1, 3 and 10 times the label or published reference exposure (DESIGN CHOICE), with target engagement at the lowest concentration used.

| ID | Candidates | Hypothesis | Specific readouts | Refuted if | Controls |
|---|---|---|---|---|---|
| X1 | Sirolimus, Everolimus, Metformin (IGBP1 axis) | mTOR/AKT modulation lowers growth and senescence burden without raising aneuploidy | pS6, pAKT, p-4E-BP1; AMPK phosphorylation for metformin; IGBP1 status if function test available | Aneuploid fraction or senescence rises, or a relative 25 % rise; rapalogues discordant; no readout improves | Positive: pS6 loss, a second mTOR inhibitor; negative: vehicle, euploid, rescue |
| H1 | Bortezomib, Chloroquine, Bosutinib, Colchicine, Eribulin, Paclitaxel, Vincristine | Harm to constitutional aneuploid cells | Viability, apoptosis, mitotic arrest and slippage, aggregates, LC3/p62 (chloroquine), micronuclei; hypomorph versus euploid; antimitotics also at 0.01 and 0.03 times (DESIGN CHOICE) | Harm if viability is lower in the hypomorph at an equal multiple (a toxicity result, not benefit) | Positive: known aneuploidy-selective cytotoxic; negative: vehicle |
| X4 | Lazertinib, Neratinib | Differential sensitivity tracks BUB1B reduction beyond proliferation | Viability against checkpoint strength, EGFR/HER2 pathway marker, proliferation-matched comparison | No difference after matching proliferation, or hypomorph non-tumour cells equally affected | Resistant line; vehicle |
| X2 | Phenylbutanoic Acid | Lowers aggregation and senescence markers | Aggresome dye, UPR markers, p21, SA-beta-gal | Senescence or aneuploidy rises, none improves | Unrelated chaperone; equal-osmolarity inert salt |
| X5 | Ceritinib, Selinexor, Tirbanibulin | Selective effect on cells with BUB1B loss | Target engagement, mitotic fidelity, viability in hypomorph versus euploid | Selectivity absent | Vehicle; resistant line |
| X3 | Trametinib | Aneuploid tumour cells more sensitive | pERK, viability against aneuploidy score | Sensitivity not tracking aneuploidy, or non-tumour hypomorph equally sensitive | pERK loss; resistant line |

Conditional experiments for rows with a computational step first reuse the X design only if the step does not close the row.

Replicates from `power_analysis.json` (aneuploid fraction, power 0.8, two-sided alpha 0.05, 50 cells per culture, ICC 0.05 to 0.3 assumed; independent cells per arm, then cultures):

| Baseline | Relative reduction | Cells per arm | Cultures per arm |
|---|---|---|---|
| 12 % | 25 % | 1,632 | 113 to 513 |
| 12 % | 50 % | 348 | 25 to 110 |
| 25 % | 25 % | 684 | 48 to 215 |
| 25 % | 50 % | 150 | 11 to 48 |
| 39 % | 50 % | 84 | 6 to 27 |

Use the larger of the arcsine and exact numbers. The file sizes only the aneuploid fraction; viability and imaging endpoints are not sized by the project: at least three independent experiments in triplicate (DESIGN CHOICE), and a null at low power is uninformative, not a refutation.

Order: E0, M1, X1, H1, X4, X2, X5, X3.

## 5. Stop-go and what would not be concluded

| Stage | Go if | Stop if |
|---|---|---|
| E0 | Aneuploidy confirmed or missense impaired | S0 |
| M1 | Hypomorph shows a baseline phenotype that rescue reverses | No phenotype: revise the model, do not test drugs |
| CF1 to CF9 | Row survives its criterion | Criterion met: close the row |
| X1, X2, X4, X5, X3 | No refutation criterion met, reproduced in two clones | Any criterion met; a refuted drug is a recorded dead end |
| H1 | Harm closes the drug | Harm is itself the result |

A drug that survives is not accepted, only not yet refuted. Not concluded from a positive result: clinical benefit, safety or dose for any person including the child; effect in the child's tissues or age group; that label status means suitability; that a boxed warning is outweighed. Line C hits come from cancer lines and do not model partial loss.

## 6. Resources (plain terms)

Cell lines: patient-derived (if consented), isogenic hypomorph with two clones, rescue, euploid controls, tumour lines. Assays: karyotype, single-cell sequencing, immunofluorescence for micronuclei, live-cell imaging, immunoblots, viability and senescence assays. Equipment: cell-culture facility, high-content or confocal microscope with live chamber, sequencing core, cytogenetics laboratory. Staff: cell biologist, cytogeneticist, bioinformatician. Time is not estimated here.

## Inconsistencies found

1. Age bands: old file uses 0-<3, 3-<12 and 12+; v2 outputs use 0-4, 5-11, 12-17.
2. The old file states no drug has two independent lines; Neratinib now has two (C and EP).
3. Line C protocol states viability correlations, so positive rho (Lazertinib, Neratinib) may mean lower sensitivity; unresolved (CF1).
4. The old file's nine drugs include Selumetinib, now in no line.
