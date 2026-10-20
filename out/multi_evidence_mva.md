# Multi-evidence integration of proposed candidates (case output: `out`)

*hypotheses for follow-up; low evidence; the integration orders hypotheses and creates no evidence; no clinical recommendation, dose or claim of efficacy.*

Protocol: `protocol/addendum_multi_evidence.md (post hoc)`. Seed 20261018, 10000 weight vectors over 11 criteria; rule: P(top 5) >= 0.5, not hard-excluded; primary = top 3.

## 1. Pool

84 approved small molecules entered the pool (81 after hard exclusions). Members by source (a drug can be in several): A 31, B1 13, B2 3, D 14, E 9, ENG 28, EP 7, literature 19. 152 candidates had no approval record in the project outputs or the local openFDA cache (offline: not in the pool, so the pool is incomplete for them) and 3 were checked and are not eligible.

Hard exclusions: BORTEZOMIB (hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)); BOSUTINIB (hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)); CHLOROQUINE (hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)).

## 2. Proposed candidates (hypotheses, low evidence)

| Drug | Status | P(top 5) | P(top 3) | P(rank 1) | P(top 10) | Confidence | Evidence tier |
|---|---|--:|--:|--:|--:|--:|---|
| Trametinib | PRIMARY | 0.502 | 0.397 | 0.210 | 0.667 | 0.82 | level 1 literature (2 help, 0 harm papers); 2 lines supporting, 1 contrary |

### Evidence matrix of the shortlist (grade in [0, 1]; `-` missing; arrows: + supportive, = neutral, x contrary)

| Drug | M | LIT | A | B1 | B2 | C | D | E | EP | ENG | S | Age 0-4 | 5-11 | 12-17 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|---|---|---|
| Trametinib | 1.00+ | 0.68+ | 0.62= | 0.44= | - | - | 1.00+ | 0.07x | 0.00x | 0.94+ | 0.50= | partly | labelled | labelled |

Central weight vectors (weights that put the drug first, or at its best rank):

- Trametinib: M 0.14, LIT 0.12, A 0.09, B1 0.13, B2 0.06, C 0.08, D 0.10, E 0.07, EP 0.05, ENG 0.11, S 0.05

## 3. Rank acceptability, top 15 by P(top 5)

| Drug | P(top 5) | P(top 3) | P(rank 1) | P(top 10) | Mean rank | Missing criteria | Proposed | Evidence tier | Age 0-4 / 5-11 / 12-17 |
|---|--:|--:|--:|--:|--:|--:|---|---|---|
| Trametinib | 0.502 | 0.397 | 0.210 | 0.667 | 9.0 | 2 | yes | level 1 literature (2 help, 0 harm papers); 2 lines supporting, 1 contrary | partly / labelled / labelled |
| Brigatinib | 0.432 | 0.245 | 0.064 | 0.730 | 8.3 | 7 | no | no verified literature entry; 1 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Minocycline | 0.386 | 0.271 | 0.107 | 0.563 | 12.9 | 7 | no | no verified literature entry; 1 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Pindolol | 0.378 | 0.238 | 0.070 | 0.621 | 10.5 | 5 | no | no verified literature entry; 2 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Ceritinib | 0.368 | 0.249 | 0.097 | 0.611 | 10.5 | 4 | no | no verified literature entry; 2 lines supporting, 0 contrary | not_labelled / not_labelled / not_labelled |
| Selinexor | 0.367 | 0.280 | 0.130 | 0.495 | 15.7 | 4 | no | level 1 literature (1 help, 0 harm papers); 2 lines supporting, 0 contrary | not_labelled / not_labelled / not_labelled |
| Tyrosine | 0.322 | 0.246 | 0.115 | 0.451 | 18.5 | 8 | no | no verified literature entry; 1 lines supporting, 0 contrary | verify / verify / verify |
| Fluvastatin | 0.293 | 0.200 | 0.094 | 0.477 | 14.9 | 5 | no | no verified literature entry; 2 lines supporting, 0 contrary | verify / verify / verify |
| Zinc Chloride | 0.233 | 0.148 | 0.000 | 0.352 | 25.8 | 8 | no | no verified literature entry; 1 lines supporting, 0 contrary | verify / verify / verify |
| Lovastatin | 0.211 | 0.088 | 0.008 | 0.484 | 14.0 | 5 | no | no verified literature entry; 1 lines supporting, 0 contrary | verify / verify / verify |
| Entrectinib | 0.162 | 0.055 | 0.000 | 0.433 | 14.5 | 7 | no | no verified literature entry; 1 lines supporting, 0 contrary | partly / labelled / labelled |
| Eribulin | 0.138 | 0.084 | 0.013 | 0.222 | 32.2 | 7 | no | level 1 literature (1 help, 0 harm papers); 0 lines supporting, 0 contrary | not_labelled / not_labelled / not_labelled |
| Belinostat | 0.131 | 0.056 | 0.010 | 0.369 | 14.4 | 4 | no | no verified literature entry; 2 lines supporting, 1 contrary | not_labelled / not_labelled / not_labelled |
| Carvedilol | 0.124 | 0.046 | 0.009 | 0.391 | 14.6 | 4 | no | no verified literature entry; 2 lines supporting, 0 contrary | verify / verify / verify |
| Arsenic Trioxide | 0.123 | 0.070 | 0.025 | 0.244 | 24.5 | 7 | no | no verified literature entry; 2 lines supporting, 0 contrary | verify / verify / verify |

## 4. Age bands (label support is a filter and a column, not evidence of efficacy)

Number of pool drugs with a label reading per band: 0-4 46, 5-11 46, 12-17 46. Unknown age is shown as the union; the age is never assumed.

Headline reading (conservative): a drug is shown for a band only if it is proposed in the full pool and passes the band filter. The restricted-pool reading re-ranks only the drugs that pass the filter; with a small pool the top-5 test is permissive (declared in the addendum, section 10) and it is shown for information only.

- **0-4**: 10 pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: Trametinib. Primary passing: Trametinib. Restricted-pool reading (information only): Trametinib (1.00), Entrectinib (0.96), Everolimus (0.86), Selpercatinib (0.67), Selumetinib (0.65), Larotrectinib (0.55).
- **5-11**: 11 pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: Trametinib. Primary passing: Trametinib. Restricted-pool reading (information only): Trametinib (0.99), Entrectinib (0.93), Everolimus (0.84), Selpercatinib (0.65), Selumetinib (0.64), Larotrectinib (0.54).
- **12-17**: 12 pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: Trametinib. Primary passing: Trametinib. Restricted-pool reading (information only): Trametinib (0.99), Entrectinib (0.93), Everolimus (0.84), Selpercatinib (0.64), Selumetinib (0.62), Larotrectinib (0.51).
- **unknown_age_union**: 12 pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: Trametinib. Primary passing: Trametinib. Restricted-pool reading (information only): Trametinib (0.99), Entrectinib (0.93), Everolimus (0.84), Selpercatinib (0.64), Selumetinib (0.62), Larotrectinib (0.51).

## 5. Sensitivity

Baseline proposed set: Trametinib.

| Variant | Proposed set | Jaccard | Baseline retained |
|---|---|--:|--:|
| leave out M | none | 0.00 | 0 |
| leave out LIT | none | 0.00 | 0 |
| leave out A | Tyrosine, Zinc Chloride, Selinexor | 0.00 | 0 |
| leave out B1 | none | 0.00 | 0 |
| leave out B2 | Trametinib | 1.00 | 1 |
| leave out C | Trametinib | 1.00 | 1 |
| leave out D | Minocycline | 0.00 | 0 |
| leave out E | Minocycline, Trametinib | 0.50 | 1 |
| leave out EP | Trametinib | 1.00 | 1 |
| leave out ENG | Brigatinib, Pindolol, Ceritinib, Minocycline | 0.00 | 0 |
| leave out S | Trametinib, Selinexor | 0.50 | 1 |
| literature-only weights (LIT, M) | Eribulin, Trametinib, Selinexor, Vincristine, Temozolomide | 0.20 | 1 |
| Track 1 non-BUB1B genes off | Trametinib | 1.00 | 1 |

Equal weights (one vector), top 5: Trametinib (0.531), Brigatinib (0.514), Ceritinib (0.511), Pindolol (0.503), Minocycline (0.496). Overlap with the baseline proposed set (Jaccard) 0.20.

Track 1 gene block (cap 0.10): drugs with a bonus: Cabazitaxel (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Carboplatin (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Colchicine (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Eribulin (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Etoposide (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Everolimus (IGBP1, T=1.0), Metformin (IGBP1, T=0.5), Olaparib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Paclitaxel (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Selinexor (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Selpercatinib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Selumetinib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Sirolimus (IGBP1/CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=1.0), Temozolomide (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Temsirolimus (IGBP1/CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=1.0), Trametinib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Treosulfan (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Vincristine (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25). Change in P(top 5) when the genes are off: Temsirolimus +0.014, Brigatinib -0.003, Selpercatinib +0.003, Everolimus +0.002, Carboplatin +0.002, Pindolol -0.002, Belinostat -0.002, Minocycline -0.002. These are hypotheses to test (karyotype or single-cell sequencing for the megabase events; function tests for IGBP1), not findings.

Monte Carlo on label errors (levels or directions of verified literature entries corrupted; pool and exclusions fixed):

| Error rate | Expected baseline candidates retained | Mean set size | Mean Jaccard | Retained fraction by candidate |
|---|--:|--:|--:|---|
| 10 % | 0.62 of 1 | 0.62 | 0.62 | Trametinib 0.62 |
| 20 % | 0.52 of 1 | 0.52 | 0.52 | Trametinib 0.52 |
| 30 % | 0.41 of 1 | 0.41 | 0.41 | Trametinib 0.41 |

## 6. Comparison with the binary rule and with the nine drugs of the report

The earlier rule (at least two usable computational lines) gave candidates: {"b2_counted_as_prespecified": 0, "b2_held_exploratory": 0}. The new rule proposes 1.

| Report drug | Group | Status | Rank by P(top 5) | P(top 5) |
|---|---|---|--:|--:|
| Sirolimus | A constitutional | in pool, not proposed | 73 | 0.0 |
| Everolimus | A constitutional (sibling) | in pool, not proposed | 30 | 0.0119 |
| Metformin | A constitutional | in pool, not proposed | 47 | 0.0002 |
| Phenylbutanoic Acid | A constitutional (sodium phenylbutyrate) | in pool, not proposed | 40 | 0.003 |
| Trametinib | B tumour context | proposed | 1 | 0.5017 |
| Selumetinib | B tumour context (sibling) | in pool, not proposed | 23 | 0.0422 |
| Bortezomib | C hazard | hard-excluded | - | - |
| Chloroquine | C hazard | hard-excluded | - | - |
| Bosutinib | C hazard | hard-excluded | - | - |

Proposed here but not among the nine: none.

## 7. Data gaps and limits

- engine ranking stored only for the top 20 and the nine report drugs: ENG is missing for the other pool drugs and 'top 50 of the engine' is the stored top 20
- line C gives no drug-level result: missing for every drug
- confidences of the five consequences are taken from report section 3.1 (out/differential_diagnosis.json holds the gene-level differential only)
- age-band label support is assessed for the drugs of the label table only; the rest are 'not_assessed'
- The integration orders hypotheses; it does not create evidence. Pool membership comes from the same lines that are integrated, literature labels were assigned by AI search tasks and not reviewed by a specialist, and weights are uniform over a simplex, not elicited.
- Line grades are positions in lists (A, B1, B2, ENG) or exploratory signals (D, E, EP); the patient signature comes from one adult cell line. A rank-based grade is not an effect size.
- Label support by age describes the US label, not benefit. No clinical recommendation, no dosing, no claim of efficacy.
