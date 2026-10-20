# Multi-evidence integration: original and completed runs side by side

*hypotheses for follow-up; low evidence; the integration orders hypotheses and creates no evidence; no clinical recommendation, dose or claim of efficacy.*

Protocol: `protocol/addendum_multi_evidence.md`, section 11 (completion fixed before the rerun). The original outputs are unchanged; the proposal rule, weights, SMAA settings and sensitivity analyses are the same in both runs.

## Answer

- Trametinib does not survive: it is not a proposed candidate in the completed run (P(top 5) 0.418 in the completed run; 0.502 in the original).
- New proposed drugs in the completed run: Calcitriol.
- Calcitriol has no verified literature entry and no mechanism fit (LIT and M missing, 4 of 11 criteria missing, confidence factor 0.64): it is proposed only because of list-position and exploratory lines, so it is the weakest kind of proposal, and its stability under label errors (Monte Carlo) is trivially 1 because it has no literature labels to corrupt. It is a hypothesis to check against the label and the literature, not a finding.
- Calcitriol, flag to verify (not part of the pre-specified safety criterion, which found no hazard flag in the project files): the record of the child documents nephrocalcinosis (report section 3.1) and vitamin D analogues act on calcium handling, so the label must be read for hypercalcaemia and kidney warnings before any follow-up.

## What was completed

- Approval: the project's FDA check was run for 152 candidates without a record; 124 were confirmed (approved small molecules among them entered the pool); 28 were not confirmed.
- Engine: frozen B0 score and rank stored for every drug (consistent with the stored top 20: True); 'top 50 of the engine' now uses the full ranking.
- Line C at drug level: 190 pool drugs in PRISM, 190 tested, 23 with BH q < 0.10. Correlational, from cancer cell lines, not adjusted for proliferation (the score correlates 0.81 with a proliferation score) or lineage.

## Side by side

| | Original | Completed |
|---|---|---|
| Pool size (after exclusions) | 84 (81) | 224 (221) |
| Proposed candidates | Trametinib | Calcitriol |
| Primary | Trametinib | Calcitriol |
| Best by P(top 5) | Trametinib 0.50; Brigatinib 0.43; Minocycline 0.39; Pindolol 0.38; Ceritinib 0.37 | Calcitriol 0.74; Trametinib 0.42; Selinexor 0.42; Dexmedetomidine 0.39; Acyclovir 0.36 |
| Equal weights, top 5 | Trametinib, Brigatinib, Ceritinib, Pindolol, Minocycline | Calcitriol, Dexmedetomidine, Trametinib, Selinexor, Acyclovir |
| Literature-only weights, proposed | Eribulin, Trametinib, Selinexor, Vincristine, Temozolomide | Eribulin, Trametinib, Selinexor, Vincristine, Temozolomide |
| Track 1 non-BUB1B genes off, proposed | Trametinib | Calcitriol |
| Label errors 10 %: baseline candidates retained | 0.62 of 1 | 1.00 of 1 |
| Label errors 20 %: baseline candidates retained | 0.52 of 1 | 1.00 of 1 |
| Label errors 30 %: baseline candidates retained | 0.41 of 1 | 1.00 of 1 |

Leave-one-criterion-out, proposed set:

| Left out | Original | Completed |
|---|---|---|
| M | none | Calcitriol |
| LIT | none | Calcitriol |
| A | Tyrosine, Zinc Chloride, Selinexor | Calcitriol, Selinexor, Tyrosine, Zinc Chloride |
| B1 | none | Calcitriol |
| B2 | Trametinib | Calcitriol |
| C | Trametinib | Dexmedetomidine, Calcitriol |
| D | Minocycline | Acyclovir, Dexmedetomidine, Calcitriol, Minocycline |
| E | Minocycline, Trametinib | Calcitriol |
| EP | Trametinib | Trametinib, Calcitriol |
| ENG | Brigatinib, Pindolol, Ceritinib, Minocycline | Calcitriol |
| S | Trametinib, Selinexor | Calcitriol, Trametinib, Selinexor |

Age bands (proposed in the full pool and passing the label filter):

| Band | Original | Completed |
|---|---|---|
| 0-4 | Trametinib | none |
| 12-17 | Trametinib | none |
| 5-11 | Trametinib | none |
| unknown_age_union | Trametinib | none |

The nine drugs of the report (status, rank by P(top 5), P(top 5)):

| Drug | Original | Completed |
|---|---|---|
| Bortezomib | hard-excluded | hard-excluded |
| Bosutinib | hard-excluded | hard-excluded |
| Chloroquine | hard-excluded | hard-excluded |
| Everolimus | in pool, not proposed, rank 30, P 0.012 | in pool, not proposed, rank 50, P 0.004 |
| Metformin | in pool, not proposed, rank 47, P 0.000 | in pool, not proposed, rank 85, P 0.000 |
| Phenylbutanoic Acid | in pool, not proposed, rank 40, P 0.003 | in pool, not proposed, rank 66, P 0.001 |
| Selumetinib | in pool, not proposed, rank 23, P 0.042 | in pool, not proposed, rank 32, P 0.017 |
| Sirolimus | in pool, not proposed, rank 73, P 0.000 | in pool, not proposed, rank 197, P 0.000 |
| Trametinib | proposed, rank 1, P 0.502 | in pool, not proposed, rank 2, P 0.418 |

---

# Completed run in detail


*hypotheses for follow-up; low evidence; the integration orders hypotheses and creates no evidence; no clinical recommendation, dose or claim of efficacy.*

Protocol: `protocol/addendum_multi_evidence.md (post hoc; section 11 completion)`. Seed 20261018, 10000 weight vectors over 11 criteria; rule: P(top 5) >= 0.5, not hard-excluded; primary = top 3.

## 1. Pool

224 approved small molecules entered the pool (221 after hard exclusions). Members by source (a drug can be in several): A 47, B1 30, B2 3, D 40, E 40, ENG 50, EP 36, literature 28. 30 candidates had no approval record in the project outputs or the local openFDA cache (offline: not in the pool, so the pool is incomplete for them) and 0 were checked and are not eligible.

Hard exclusions: BORTEZOMIB (hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)); BOSUTINIB (hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)); CHLOROQUINE (hazard: verified entries state selective killing of aneuploid cells (hazard for a constitutional mosaic)).

## 2. Proposed candidates (hypotheses, low evidence)

| Drug | Status | P(top 5) | P(top 3) | P(rank 1) | P(top 10) | Confidence | Evidence tier |
|---|---|--:|--:|--:|--:|--:|---|
| Calcitriol | PRIMARY | 0.742 | 0.563 | 0.257 | 0.900 | 0.64 | no verified literature entry; 2 lines supporting, 0 contrary |

### Evidence matrix of the shortlist (grade in [0, 1]; `-` missing; arrows: + supportive, = neutral, x contrary)

| Drug | M | LIT | A | B1 | B2 | C | D | E | EP | ENG | S | Age 0-4 | 5-11 | 12-17 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|---|---|---|
| Calcitriol | - | - | 0.57= | - | - | 0.68= | 1.00+ | 0.33x | 0.80+ | 0.93+ | 1.00+ | not_assessed | not_assessed | not_assessed |

Central weight vectors (weights that put the drug first, or at its best rank):

- Calcitriol: M 0.06, LIT 0.07, A 0.06, B1 0.11, B2 0.06, C 0.11, D 0.11, E 0.09, EP 0.12, ENG 0.10, S 0.12

## 3. Rank acceptability, top 15 by P(top 5)

| Drug | P(top 5) | P(top 3) | P(rank 1) | P(top 10) | Mean rank | Missing criteria | Proposed | Evidence tier | Age 0-4 / 5-11 / 12-17 |
|---|--:|--:|--:|--:|--:|--:|---|---|---|
| Calcitriol | 0.742 | 0.563 | 0.257 | 0.900 | 5.2 | 4 | yes | no verified literature entry; 2 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Trametinib | 0.418 | 0.334 | 0.158 | 0.545 | 20.6 | 1 | no | level 1 literature (2 help, 0 harm papers); 2 lines supporting, 1 contrary | partly / labelled / labelled |
| Selinexor | 0.416 | 0.349 | 0.182 | 0.516 | 27.2 | 2 | no | level 1 literature (1 help, 0 harm papers); 3 lines supporting, 0 contrary | not_labelled / not_labelled / not_labelled |
| Dexmedetomidine | 0.390 | 0.250 | 0.095 | 0.598 | 14.4 | 7 | no | no verified literature entry; 1 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Acyclovir | 0.357 | 0.200 | 0.028 | 0.614 | 15.2 | 5 | no | no verified literature entry; 1 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Minocycline | 0.287 | 0.162 | 0.043 | 0.530 | 17.4 | 5 | no | no verified literature entry; 2 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Pimozide | 0.271 | 0.163 | 0.059 | 0.481 | 18.5 | 3 | no | no verified literature entry; 3 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Ceritinib | 0.238 | 0.170 | 0.061 | 0.390 | 20.2 | 2 | no | no verified literature entry; 2 lines supporting, 0 contrary | not_labelled / not_labelled / not_labelled |
| Neratinib | 0.211 | 0.125 | 0.029 | 0.354 | 28.3 | 3 | no | no verified literature entry; 2 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Crizotinib | 0.172 | 0.067 | 0.002 | 0.421 | 17.9 | 3 | no | no verified literature entry; 1 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Binimetinib | 0.167 | 0.079 | 0.009 | 0.387 | 22.8 | 6 | no | no verified literature entry; 1 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Paliperidone | 0.139 | 0.079 | 0.028 | 0.264 | 36.2 | 5 | no | no verified literature entry; 3 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Tyrosine | 0.129 | 0.072 | 0.013 | 0.233 | 63.2 | 8 | no | no verified literature entry; 1 lines supporting, 0 contrary | verify / verify / verify |
| Pindolol | 0.127 | 0.055 | 0.005 | 0.357 | 20.1 | 3 | no | no verified literature entry; 2 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |
| Nelfinavir | 0.090 | 0.035 | 0.006 | 0.260 | 34.9 | 5 | no | no verified literature entry; 1 lines supporting, 0 contrary | not_assessed / not_assessed / not_assessed |

## 4. Age bands (label support is a filter and a column, not evidence of efficacy)

Number of pool drugs with a label reading per band: 0-4 46, 5-11 46, 12-17 46. Unknown age is shown as the union; the age is never assumed.

Headline reading (conservative): a drug is shown for a band only if it is proposed in the full pool and passes the band filter. The restricted-pool reading re-ranks only the drugs that pass the filter; with a small pool the top-5 test is permissive (declared in the addendum, section 10) and it is shown for information only.

- **0-4**: 10 pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: none. Primary passing: none. Proposed but outside the filter: Calcitriol [not_assessed]. Restricted-pool reading (information only): Trametinib (1.00), Entrectinib (0.98), Everolimus (0.86), Selumetinib (0.73), Selpercatinib (0.67).
- **5-11**: 11 pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: none. Primary passing: none. Proposed but outside the filter: Calcitriol [not_assessed]. Restricted-pool reading (information only): Trametinib (1.00), Entrectinib (0.97), Everolimus (0.83), Selumetinib (0.73), Selpercatinib (0.63).
- **12-17**: 12 pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: none. Primary passing: none. Proposed but outside the filter: Calcitriol [not_assessed]. Restricted-pool reading (information only): Trametinib (1.00), Entrectinib (0.96), Everolimus (0.83), Selumetinib (0.71), Selpercatinib (0.61).
- **unknown_age_union**: 12 pool drugs pass the filter (labelled or partly). Proposed in the full pool and passing: none. Primary passing: none. Proposed but outside the filter: Calcitriol [not_assessed/not_assessed/not_assessed]. Restricted-pool reading (information only): Trametinib (1.00), Entrectinib (0.96), Everolimus (0.83), Selumetinib (0.71), Selpercatinib (0.61).

## 5. Sensitivity

Baseline proposed set: Calcitriol.

| Variant | Proposed set | Jaccard | Baseline retained |
|---|---|--:|--:|
| leave out M | Calcitriol | 1.00 | 1 |
| leave out LIT | Calcitriol | 1.00 | 1 |
| leave out A | Calcitriol, Selinexor, Tyrosine, Zinc Chloride | 0.25 | 1 |
| leave out B1 | Calcitriol | 1.00 | 1 |
| leave out B2 | Calcitriol | 1.00 | 1 |
| leave out C | Dexmedetomidine, Calcitriol | 0.50 | 1 |
| leave out D | Acyclovir, Dexmedetomidine, Calcitriol, Minocycline | 0.25 | 1 |
| leave out E | Calcitriol | 1.00 | 1 |
| leave out EP | Trametinib, Calcitriol | 0.50 | 1 |
| leave out ENG | Calcitriol | 1.00 | 1 |
| leave out S | Calcitriol, Trametinib, Selinexor | 0.33 | 1 |
| literature-only weights (LIT, M) | Eribulin, Trametinib, Selinexor, Vincristine, Temozolomide | 0.00 | 0 |
| Track 1 non-BUB1B genes off | Calcitriol | 1.00 | 1 |

Equal weights (one vector), top 5: Calcitriol (0.620), Dexmedetomidine (0.571), Trametinib (0.570), Selinexor (0.567), Acyclovir (0.566). Overlap with the baseline proposed set (Jaccard) 0.20.

Track 1 gene block (cap 0.10): drugs with a bonus: Cabazitaxel (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Carboplatin (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Colchicine (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Eribulin (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Etoposide (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Everolimus (IGBP1, T=1.0), Metformin (IGBP1, T=0.5), Olaparib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Paclitaxel (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Pralsetinib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Selinexor (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Selpercatinib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.5), Selumetinib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Sirolimus (IGBP1/CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=1.0), Temozolomide (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Temsirolimus (IGBP1/CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=1.0), Trametinib (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Treosulfan (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25), Vincristine (CENPE/MAD1L1/CDC20/MAD2L1BP/TTK/CEP192, T=0.25). Change in P(top 5) when the genes are off: Temsirolimus +0.014, Calcitriol -0.002, Binimetinib -0.002, Minocycline -0.002, Pimozide -0.002, Selumetinib +0.001, Paliperidone -0.001, Pindolol -0.001. These are hypotheses to test (karyotype or single-cell sequencing for the megabase events; function tests for IGBP1), not findings.

Monte Carlo on label errors (levels or directions of verified literature entries corrupted; pool and exclusions fixed):

| Error rate | Expected baseline candidates retained | Mean set size | Mean Jaccard | Retained fraction by candidate |
|---|--:|--:|--:|---|
| 10 % | 1.00 of 1 | 1.00 | 1.00 | Calcitriol 1.00 |
| 20 % | 1.00 of 1 | 1.00 | 1.00 | Calcitriol 1.00 |
| 30 % | 1.00 of 1 | 1.00 | 1.00 | Calcitriol 1.00 |

## 6. Comparison with the binary rule and with the nine drugs of the report

The earlier rule (at least two usable computational lines) gave candidates: {"b2_counted_as_prespecified": 0, "b2_held_exploratory": 0}. The new rule proposes 1.

| Report drug | Group | Status | Rank by P(top 5) | P(top 5) |
|---|---|---|--:|--:|
| Sirolimus | A constitutional | in pool, not proposed | 197 | 0.0 |
| Everolimus | A constitutional (sibling) | in pool, not proposed | 50 | 0.0038 |
| Metformin | A constitutional | in pool, not proposed | 85 | 0.0001 |
| Phenylbutanoic Acid | A constitutional (sodium phenylbutyrate) | in pool, not proposed | 66 | 0.001 |
| Trametinib | B tumour context | in pool, not proposed | 2 | 0.4181 |
| Selumetinib | B tumour context (sibling) | in pool, not proposed | 32 | 0.0166 |
| Bortezomib | C hazard | hard-excluded | - | - |
| Chloroquine | C hazard | hard-excluded | - | - |
| Bosutinib | C hazard | hard-excluded | - | - |

Proposed here but not among the nine: Calcitriol.

## 7. Data gaps and limits

- 28 of the 152 candidates were not confirmed by the approval check (listed in the completion file); the check uses name matching and some may be approved under another name
- line C is correlational, from cancer cell lines, not adjusted for proliferation or lineage, and restricted to pool drugs
- confidences of the five consequences are taken from report section 3.1 (out/differential_diagnosis.json holds the gene-level differential only)
- age-band label support is assessed for the drugs of the label table only; the rest are 'not_assessed'
- The integration orders hypotheses; it does not create evidence. Pool membership comes from the same lines that are integrated, literature labels were assigned by AI search tasks and not reviewed by a specialist, and weights are uniform over a simplex, not elicited.
- Line grades are positions in lists (A, B1, B2, ENG) or exploratory signals (D, E, EP); the patient signature comes from one adult cell line. A rank-based grade is not an effect size.
- Label support by age describes the US label, not benefit. No clinical recommendation, no dosing, no claim of efficacy.
