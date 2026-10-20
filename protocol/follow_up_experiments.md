# Follow-up experiments for the nine BUB1B hypotheses

Status: design for a cytogenetics or cell-biology laboratory. In vitro only. It gives no dose, no schedule and no treatment recommendation for any person. No experiment has been run. The case is described only by gene symbol (*BUB1B*) and the categorical profile: a stop-gain predicted to undergo nonsense-mediated decay, and a missense variant of uncertain significance at a conserved, buried kinase-domain position; phase unresolved; age unknown. All numbers below come from `submissions/candeleroj98_track2_report.md` (Table 1, S12.4) and `out/empirical/power_analysis.json`. Fold-exposure ranges and replicate counts not in those files are marked DESIGN CHOICE.

Computational context: no drug is supported by two independent lines, and line E (signature reversal) found none. Everything here therefore tests weak, literature-level hypotheses (best level L1; everolimus and selumetinib L0).

## 1. Step 0: diagnostic experiment that conditions everything

Done in a clinical-diagnostic setting, with consent, before any drug work.

| Test | Purpose | Method |
|---|---|---|
| Karyotype (metaphase spreads) of blood lymphocytes, plus skin fibroblasts if a biopsy is clinically justified | Confirm and quantify mosaic aneuploidy; the report states that mosaicism is not shown and that the read-depth screen is not a karyotype | Count at least the cells needed for the fraction expected (see section 3 baselines); report the fraction of aneuploid cells and the chromosome gains and losses |
| Single-cell or shallow single-cell sequencing of the same material | Quantify low-fraction or tissue-restricted mosaic that a karyotype can miss | Per-cell copy number; aneuploid fraction with confidence interval |
| Phase of the two variants | Compound heterozygosity (one variant per copy) is required for the two variants to explain MVA | Parental testing, or long-read or allele-specific sequencing |
| Functional test of the missense variant | It is of uncertain significance | Express wild-type and variant BubR1 in a BUB1B-depleted human line: protein level (the literature describes high turnover of missense mutants), kinetochore localisation, checkpoint-arrest assay under a spindle poison |

**Result that stops the drug experiments (stop rule S0):**
- Karyotype and single-cell sequencing both show no excess aneuploidy above the control range of the testing laboratory, in every tissue examined, and
- the missense variant behaves like wild type in all three functional assays, or the variants are shown to lie on the same copy (in cis).

Under S0 the premise of the report (reduced BubR1 function causing mosaic aneuploidy) is not supported for this child, and sections 2 to 4 are not run for the child's cells. The cell-model experiments on the isogenic model may still be run as general research, labelled as such.
If aneuploidy is confirmed or the missense variant is functionally impaired, proceed. The measured aneuploid fraction selects the planning row in section 3.

## 2. Cell models

| Model | Description | Notes |
|---|---|---|
| Patient-derived fibroblasts, or lymphoblastoid lines | Only with consent and only if a clinically indicated sample exists | Low passage; report passage number, because cultured cells may lose aneuploid clones (caveat of `power_analysis.json`). Karyotype each batch at start and end |
| Isogenic hypomorph | A human fibroblast or myoblast line with partial BUB1B reduction (hypomorph, not full knockout), made by a hypomorphic allele or titrated knockdown | Confirm residual BubR1 protein by immunoblot and kinetochore signal; a full knockout is lethal-prone and does not model the child. Use two independent clones or two independent guides |
| Euploid control | Parental line of the isogenic model, and a euploid donor fibroblast line | Needed for hazard readouts and for patient-versus-control comparisons (the background aneuploid rate in controls was not sourced in S12.4) |
| Rescue control | Hypomorph with wild-type BUB1B re-expressed | Shows that any phenotype is due to BUB1B reduction |
| Tumour-context lines (groups B and C) | Rhabdomyosarcoma lines with graded aneuploidy and spindle-checkpoint strength; the child's tumour cells only if available | Always run beside the non-tumour hypomorph and euploid lines |

All arms are blinded for scoring and randomised on the plate. Vehicle is matched to each compound's solvent.

## 3. Per-drug design

### 3.1 Common readouts

| Readout | Method | Used for |
|---|---|---|
| Micronuclei and lagging chromosomes | Fixed-cell imaging (DAPI, centromere and kinetochore stain, cytokinesis-block) | All |
| Chromosome mis-segregation | Live-cell imaging with a histone reporter, mitotic timing and error frequency | All |
| Aneuploid fraction | Single-cell sequencing or FISH or chromosome spreads | Primary endpoint for group A |
| mTOR pathway activity | Phospho-S6 (and phospho-4E-BP1) by immunoblot or imaging | Sirolimus, everolimus, metformin |
| DNA damage | gamma-H2AX and 53BP1 foci | All |
| Senescence | SA-beta-gal, p21 | All |
| Growth | Proliferation curve, cell count, confluence imaging | All |
| Protein stress | Aggregates (aggresome dye), unfolded-protein-response markers | Phenylbutyrate, hazards |
| Viability and selectivity | Viability against aneuploidy score, hypomorph versus euploid | Groups B and C |

### 3.2 Concentration range (in vitro design only)

For each drug, define the reference exposure from the current label or from a published pharmacokinetic study, and express every concentration as a multiple of the plasma reference. Use at least five points spanning below and above that reference (DESIGN CHOICE: for example 0.1, 0.3, 1, 3 and 10 times), with a vehicle control. Fill the reference concentrations from the label and the literature in the laboratory notebook; this document does not state them. Note the free (unbound) fraction and the protein content of the medium. A mechanism-based readout (phospho-S6 for the mTOR inhibitors) must show target engagement at the lowest concentration used for the phenotype.

### 3.3 Hypotheses, refutation criteria and controls

Common to group A: primary endpoint is the aneuploid fraction, drug against vehicle in the hypomorph. Refutation direction follows the report's rule ("refuted if aneuploidy or senescence rise, or if no readout improves").

| Drug (group, level) | Hypothesis | Refuted if (pre-specified) | Positive control | Negative control |
|---|---|---|---|---|
| Sirolimus (A, L1) | mTORC1 inhibition improves growth or muscle-related readouts and lowers proteotoxic and senescence burden without raising aneuploidy | Aneuploid fraction or senescence markers rise versus vehicle, or none of the readouts improves; also if the aneuploid fraction rises by a relative 25 % or more (the smallest effect sized in S12.4) | Phospho-S6 loss at the lowest effective concentration; a second, structurally unrelated mTOR inhibitor | Vehicle; euploid line; rescue line |
| Everolimus (A, L0, class sibling) | Same as sirolimus | Same criteria. Both rapalogues must move the same way; discordance means a compound-specific effect and no class conclusion | As above | As above |
| Metformin (A, L1) | Reduces aneuploidy-related stress | BubR1 or MAD2 fall further, or aneuploidy rises, or no readout improves (adds AMPK activity and oxidative damage) | AMPK phosphorylation; another AMPK activator | Vehicle; AMPK inhibitor co-treatment to separate the mechanism |
| Sodium phenylbutyrate (A, L1) | Reduces protein-aggregation and senescence markers | Senescence or aneuploidy rise, or none of the readouts improves | A chemical chaperone with a different structure; aggregation induced with a proteasome inhibitor as a reporter of rescue | Vehicle; sodium salt of an inert organic acid at equal osmolarity |
| Trametinib (B, L1) | Aneuploid tumour cells are more sensitive than low-aneuploidy cells | Sensitivity does not track aneuploidy score, or non-tumour hypomorph cells are as sensitive as tumour cells | Phospho-ERK loss; a sensitive reference line | Resistant line; vehicle |
| Selumetinib (B, L0, class sibling) | Same as trametinib | Same criteria; concordance with trametinib required | As above | As above |

Pre-specified numeric thresholds for the non-aneuploidy readouts (growth, markers) are fixed before the first plate is read: a change is counted only if it exceeds the larger of the vehicle replicate range and 25 % relative change (DESIGN CHOICE). Multiplicity: the aneuploid fraction at the pharmacologically active concentration is primary; all other readouts are secondary and are not corrected for, because the power analysis does not correct for multiple readouts or multiple drugs.

### 3.4 Hazard drugs (group C: checked, not proposed)

The aim is to show harm if it exists. A positive result here is a toxicity signal, not a benefit.

| Drug | Hazard readouts | Harm criterion (pre-specified) |
|---|---|---|
| Bortezomib | Viability, apoptosis, aggregates, mitotic arrest and slippage in hypomorph versus euploid | Lower viability or higher death in the hypomorph than in the euploid line at an equal multiple of the reference exposure |
| Chloroquine | Autophagic block (LC3, p62), viability, DNA damage, aneuploidy | Same |
| Bosutinib | Viability, kinase-pathway readouts, mitotic errors | Same |

Positive control for harm: a known aneuploidy-selective cytotoxic compound. Negative control: vehicle. These three drugs also serve as tumour-context probes only in the tumour lines; they are never read as constitutional treatment candidates.

### 3.5 Replicates and power

Power uses `power_analysis.json` (two-sided alpha 0.05; arcsine method; 50 cells per culture; intraclass correlation 0.05 to 0.3, assumed). Choose the row from the baseline measured in Step 0; the published baselines are 12 % (the only BUB1B patient value), 25 % (the lower bound of MVA) and 39 % (patients without BUB1B mutations).

| Baseline | Relative reduction to detect | Cells per arm (power 0.8) | Independent cultures per arm (50 cells each) |
|---|---|---|---|
| 12 % | 25 % | 1,632 | 113 to 513 |
| 12 % | 50 % | 348 | 25 to 110 |
| 25 % | 50 % | 150 | 11 to 48 |
| 25 % | 25 % | 684 | 48 to 215 |
| 39 % | 50 % | 84 | 6 to 27 |

Rules: use the larger of the arcsine and the exact-test numbers (the exact test needs about 5 to 20 % more cells). A two-sided test also detects a rise, which is the failure case for the report. If the lab can reach only a small number of cultures, say in advance that the experiment can detect only large effects and that a null result will then be uninformative, not a refutation. Independent culture means separately seeded and treated, on different days. For the viability and imaging endpoints, at least three independent experiments (DESIGN CHOICE), each with technical triplicates.

## 4. Sequencing, cost, information and stop-go

| Order | Experiment | Cost and information | Go if | Stop if |
|---|---|---|---|---|
| 0 | Step 0 diagnostics | Low cost, highest information; conditions everything | Aneuploidy confirmed or missense variant impaired | Rule S0 |
| 1 | Build and validate isogenic hypomorph and controls (BubR1 level, baseline aneuploid fraction, micronuclei) | Moderate cost; fixes the real baseline and thus the sample size | Hypomorph shows a baseline phenotype beyond the euploid line and the rescue reverses it | No phenotype: the model does not mimic the defect; revise, do not test drugs |
| 2 | Cheap target-engagement and growth screen, all nine drugs, five concentrations | Low cost per drug; discards compounds that do not engage their target | Target engaged at a concentration inside the stated multiples | Target engagement only above the top multiple |
| 3 | Hazard panel (group C and group B in non-tumour cells) | Moderate; can end drugs early | Harm in the hypomorph relative to the euploid line closes that drug | None, harm is itself the result |
| 4 | Group A confirmatory experiments with aneuploid fraction by single-cell sequencing | Highest cost, highest information | Aneuploid fraction not higher and at least one pre-specified readout improved, reproduced in two clones | Any refutation criterion of 3.3 met |
| 5 | Group B in tumour lines against aneuploidy score | Moderate | Selectivity tracks aneuploidy and spares non-tumour cells | Selectivity absent |

A drug that survives a stage is not accepted; it is merely not yet refuted. A drug refuted in two clones is recorded as a documented dead end, which the report treats as a useful result.

## 5. What would not be concluded

- An in-vitro effect is not evidence of clinical benefit, safety or dose for any person, including the child.
- A positive result in an isogenic line does not show the effect in the child's cells, tissues or age group. The age is unknown; label status for the three age bands (0 to under 3, 3 to under 12, 12 and over) describes other diseases, not suitability.
- A null result with low power does not refute a hypothesis.
- Safety statements on the labels (for example, the boxed warning on immunosuppression and malignancy for sirolimus and everolimus) are not outweighed by any in-vitro benefit.
- The nine drugs have no support from two computational lines; passing this protocol would not change that.
- Any clinical use would require a treating physician, a pharmacist reading the full label, and ethics oversight.

## 6. Approximate resources (plain terms, no prices)

| Item | Amount |
|---|---|
| Cell lines | Patient-derived line (if consented), 1 isogenic hypomorph (2 clones), 1 rescue line, 1 to 2 euploid controls, 2 to 3 rhabdomyosarcoma lines |
| Assays | Karyotype, single-cell sequencing, immunofluorescence with micronuclei scoring, live-cell imaging, immunoblots, viability and senescence assays |
| Equipment | Cell-culture facility, confocal or high-content microscope with live-cell chamber, access to a sequencing core and a cytogenetics laboratory |
| Time (indicative) | Step 0: about 4 to 8 weeks; model building and validation: about 8 to 12 weeks; screen and hazard panel: about 6 to 8 weeks; confirmatory group A and group B experiments: about 12 to 16 weeks. The total is about 8 to 12 months (DESIGN CHOICE, not derived from the project files) |
| Staff | A cell biologist, a cytogeneticist, a bioinformatician for single-cell copy number |
