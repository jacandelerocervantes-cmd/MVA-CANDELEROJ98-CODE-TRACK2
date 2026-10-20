# Line E: signature reversal (fifth line)

Label: exploratory; hypotheses for follow-up; never evidence of efficacy. Protocol: `protocol/addendum_line_e_signature_reversal.md`. The rule of at least two usable lines is unchanged; a candidate here is a hypothesis for follow-up only.

## Data and matching
* Knockout signatures in library: 5212; chemical signatures: 5425.
* Library drugs matched to approved small molecules of the project: 1069 of 5425 library names (2036 project small molecules).
* Curated MVA genes with a knockout signature: BUB1, BUB1B; absent: BUB3, CENATAC, CEP57, MAD1L1, TRIP13.
* Track 1 genes: 7; with a knockout signature: 3 (BUB1B, CENPE, IZUMO4).
* GEO GSE22206: computed, see signature P.

## Top 20 reversal drugs (BUB1B signature, N=150)

| # | Drug | R | q (reversal) | usable | supported by |
|---|---|---|---|---|---|
| 1 | TOFACITINIB | 45 | 0.674 | False | - |
| 2 | MESTRANOL | 41 | 0.534 | False | - |
| 3 | ENZALUTAMIDE | 38 | 0.674 | False | - |
| 4 | LOSARTAN | 38 | 0.674 | False | - |
| 5 | RACECADOTRIL | 36 | 0.674 | False | - |
| 6 | CARBAMAZEPINE | 34 | 0.674 | False | - |
| 7 | SULFACETAMIDE | 34 | 0.534 | False | - |
| 8 | TRIOXSALEN | 34 | 0.534 | False | - |
| 9 | RISPERIDONE | 33 | 0.674 | False | - |
| 10 | QUETIAPINE | 32 | 0.674 | False | - |
| 11 | AMISULPRIDE | 31 | 0.674 | False | - |
| 12 | TEMOZOLOMIDE | 31 | 0.749 | False | - |
| 13 | CEFOPERAZONE | 30 | 0.68 | False | - |
| 14 | LENALIDOMIDE | 30 | 0.674 | False | - |
| 15 | NIACIN | 30 | 0.674 | False | - |
| 16 | SILDENAFIL | 30 | 0.674 | False | - |
| 17 | ESTRIOL | 29 | 0.674 | False | - |
| 18 | LEVOFLOXACIN ANHYDROUS | 28 | 0.674 | False | - |
| 19 | METHYSERGIDE | 28 | 0.534 | False | - |
| 20 | MOLINDONE | 28 | 0.674 | False | - |

## Signature P (patient-derived, GSE22206): 100 Up + 100 Down genes; 17184 genes mapped; usable 18, hazards 141

Limitation: one adult patient's lymphoblastoid cell line (n = 1 patient, 2 arrays) versus 2 controls; cannot represent the child.

| # | Drug | R | q (reversal) | usable |
|---|---|---|---|---|
| 1 | NICLOSAMIDE | 15 | 0.0593 | True |
| 2 | NERATINIB | 12 | 0.0593 | True |
| 3 | NITAZOXANIDE | 11 | 0.0593 | True |
| 4 | MASOPROCOL | 11 | 0.0593 | True |
| 5 | AMPHOTERICIN B | 10 | 0.0593 | True |
| 6 | PIMOZIDE | 9 | 0.0593 | True |
| 7 | TRIFLUOPERAZINE | 8 | 0.0593 | True |
| 8 | NELFINAVIR | 8 | 0.0593 | True |
| 9 | SALMETEROL | 8 | 0.0593 | True |
| 10 | IXAZOMIB CITRATE | 8 | 0.0593 | True |
| 11 | FLUVASTATIN | 8 | 0.0593 | True |
| 12 | PYRIMETHAMINE | 8 | 0.0593 | True |
| 13 | AURANOFIN | 8 | 0.0593 | True |
| 14 | DEXKETOPROFEN | 8 | 0.0593 | True |
| 15 | GLUCOSAMINE | 8 | 0.0593 | True |
| 16 | CRIZOTINIB | 7 | 0.107 | False |
| 17 | AMIODARONE | 7 | 0.153 | False |
| 18 | PLICAMYCIN | 7 | 0.0593 | True |
| 19 | TEMAZEPAM | 7 | 0.0593 | True |
| 20 | BENZONATATE | 6 | 0.128 | False |

## The nine report drugs

| Drug | BUB1B R | q_rev | Track1 combined R | q_rev | P R | P q_rev | line E usable | rank (BUB1B) |
|---|---|---|---|---|---|---|---|---|
| sirolimus | -6 | 0.925 | -8 | 0.948 | 0 | 1 | False | 832 |
| everolimus | -3 | 0.854 | -7 | 0.917 | -4 | 1 | False | 744 |
| metformin | 26 | 0.674 | 24 | 0.607 | -4 | 1 | False | 35 |
| sodium phenylbutyrate | 5 | 0.777 | 8 | 0.748 | -3 | 1 | False | 465 |
| trametinib | -12 | 0.929 | -7 | 0.863 | -5 | 1 | False | 963 |
| selumetinib | -12 | 0.946 | -10 | 0.917 | -7 | 1 | False | 963 |
| bortezomib | -19 | 0.97 | -19 | 0.97 | 3 | 0.869 | False | 1028 |
| chloroquine | 27 | 0.674 | 29 | 0.493 | -8 | 1 | False | 24 |
| bosutinib | -1 | 0.832 | -8 | 0.902 | 5 | 0.289 | False | 680 |

## Controls

Positive controls (inhibitors should mimic their gene knockout; rank by mimicry among matched drugs, top 5% = recovered):

* MTOR: drugs present ['sirolimus', 'everolimus', 'temsirolimus']; mimicry ranks {'Sirolimus': 2, 'Everolimus': 5, 'Temsirolimus': 12}; reversal ranks {'Sirolimus': 1067, 'Everolimus': 1065, 'Temsirolimus': 1057}; recovered (top 5%): True
* MAP2K1: drugs present ['trametinib', 'selumetinib', 'cobimetinib']; mimicry ranks {'Trametinib': 1050, 'Selumetinib': 973, 'Cobimetinib': 854}; reversal ranks {'Trametinib': 17, 'Selumetinib': 86, 'Cobimetinib': 197}; recovered (top 5%): False
* EGFR: drugs present ['erlotinib', 'gefitinib', 'afatinib', 'lapatinib']; mimicry ranks {'Erlotinib': 359, 'Gefitinib': 387, 'Afatinib': 79, 'Lapatinib': 709}; reversal ranks {'Erlotinib': 684, 'Gefitinib': 664, 'Afatinib': 982, 'Lapatinib': 331}; recovered (top 5%): False
* BRAF: drugs present ['vemurafenib', 'dabrafenib']; mimicry ranks {'Vemurafenib': 644, 'Dabrafenib': 208}; reversal ranks {'Vemurafenib': 396, 'Dabrafenib': 835}; recovered (top 5%): False
* HMGCR: drugs present ['atorvastatin', 'simvastatin', 'lovastatin', 'fluvastatin', 'rosuvastatin', 'pitavastatin', 'pravastatin']; mimicry ranks {'Atorvastatin': 7, 'Simvastatin': 18, 'Lovastatin': 30, 'Fluvastatin': 8, 'Rosuvastatin': 79, 'Pitavastatin': 1, 'Pravastatin': 930}; reversal ranks {'Atorvastatin': 1062, 'Simvastatin': 1051, 'Lovastatin': 1039, 'Fluvastatin': 1062, 'Rosuvastatin': 990, 'Pitavastatin': 1069, 'Pravastatin': 129}; recovered (top 5%): True
* ABL1: drugs present ['imatinib', 'dasatinib', 'nilotinib', 'bosutinib']; mimicry ranks {'Imatinib': 546, 'Dasatinib': 291, 'Nilotinib': 561, 'Bosutinib': 600}; reversal ranks {'Imatinib': 500, 'Dasatinib': 745, 'Nilotinib': 500, 'Bosutinib': 450}; recovered (top 5%): False
* Negative control: usable drugs per random pseudo-signature (20): mean fraction 0.0, max 0 drugs.

## Signatures and sensitivity

| Setting | Signature | usable | hazards |
|---|---|---|---|
| N=150 | gene:BUB1 | 0 | 0 |
| N=150 | gene:BUB1B | 0 | 0 |
| N=150 | gene:CENPE | 0 | 0 |
| N=150 | gene:IZUMO4 | 0 | 0 |
| N=150 | mva_set | 0 | 0 |
| N=150 | track1_combined | 0 | 0 |
| N=150 | track1_high_score_only | 0 | 0 |
| N=50 | gene:BUB1 | 0 | 0 |
| N=50 | gene:BUB1B | 0 | 0 |
| N=50 | gene:CENPE | 0 | 0 |
| N=50 | gene:IZUMO4 | 0 | 0 |
| N=50 | mva_set | 0 | 0 |
| N=50 | track1_combined | 0 | 0 |
| N=50 | track1_high_score_only | 0 | 0 |
| N=full | gene:BUB1 | 0 | 0 |
| N=full | gene:BUB1B | 0 | 0 |
| N=full | gene:CENPE | 0 | 0 |
| N=full | gene:IZUMO4 | 0 | 0 |
| N=full | mva_set | 0 | 0 |
| N=full | track1_combined | 0 | 0 |
| N=full | track1_high_score_only | 0 | 0 |
| N=P(100+100) | P | 18 | 141 |

Agreement of usable sets: gene:BUB1B N=50 vs N=150: 0 shared of 0/0; track1_combined N=50 vs N=150: 0 shared of 0/0; mva_set N=50 vs N=150: 0 shared of 0/0; gene:BUB1B N=full vs N=150: 0 shared of 0/0; track1_combined N=full vs N=150: 0 shared of 0/0; mva_set N=full vs N=150: 0 shared of 0/0

## Candidates by the combination rule
Usable sets of the existing lines: {'A': ['CYCLOSPORINE', 'DISOPYRAMIDE', 'POMALIDOMIDE'], 'B2': ['CERITINIB', 'NATALIZUMAB', 'SELINEXOR', 'TIRBANIBULIN']}.
Candidates (B2 as prespecified): 0; with B2 held exploratory: 0.
Drugs supported by line E (BUB1B and/or Track 1 combined): 18; by the Track 1 combined signature: [].

| Drug | existing lines | cytotoxic/caution | FDA label | candidate | age support (per band) |
|---|---|---|---|---|---|
| AMPHOTERICIN B | none | False | False | False | not_assessed |
| MASOPROCOL | none | True | False | False | not_assessed |
| TEMAZEPAM | none | False | False | False | not_assessed |
| NERATINIB | none | True | False | False | not_assessed |
| GLUCOSAMINE | none | False | False | False | not_assessed |
| PIMOZIDE | none | False | False | False | not_assessed |
| DEXKETOPROFEN | none | False | False | False | not_assessed |
| SALMETEROL | none | False | False | False | not_assessed |
| FLUVASTATIN | none | False | True | False | {'0_4': 'verify', '12_17': 'verify', '5_11': 'verify'} |
| NELFINAVIR | none | False | False | False | not_assessed |
| PLICAMYCIN | none | True | False | False | not_assessed |
| AURANOFIN | none | False | False | False | not_assessed |
| IXAZOMIB CITRATE | none | True | True | False | {'0_to_under_3': 'not_labelled', '12_and_older': 'adult_only', '3_to_under_12': 'not_labelled'} |
| TRIFLUOPERAZINE | none | False | False | False | not_assessed |
| PYRIMETHAMINE | none | False | False | False | not_assessed |
| NICLOSAMIDE | none | False | False | False | not_assessed |
| NITAZOXANIDE | none | False | False | False | not_assessed |
| THONZONIUM | none | False | False | False | not_assessed |

## Hazards (concordant, mimic the knockout; FDR < 0.1)

* LOSARTAN: P (BUB1B R = 38)
* RISPERIDONE: P (BUB1B R = 33)
* AMISULPRIDE: P (BUB1B R = 31)
* CEFOPERAZONE: P (BUB1B R = 30)
* LENALIDOMIDE: P (BUB1B R = 30)
* SILDENAFIL: P (BUB1B R = 30)
* ESTRIOL: P (BUB1B R = 29)
* PRAVASTATIN: P (BUB1B R = 28)
* RILMENIDINE: P (BUB1B R = 28)
* CHLOROQUINE: P (BUB1B R = 27)
* ETODOLAC: P (BUB1B R = 27)
* LORAZEPAM: P (BUB1B R = 27)
* PRAZIQUANTEL: P (BUB1B R = 27)
* ESTRADIOL: P (BUB1B R = 26)
* ESTRONE: P (BUB1B R = 26)
* ACAMPROSATE: P (BUB1B R = 25)
* OXAPROZIN: P (BUB1B R = 25)
* VALSARTAN: P (BUB1B R = 25)
* ISRADIPINE: P (BUB1B R = 24)
* VECURONIUM: P (BUB1B R = 24)
* MOXIFLOXACIN: P (BUB1B R = 23)
* NIMODIPINE: P (BUB1B R = 22)
* QUININE: P (BUB1B R = 22)
* ZOLPIDEM: P (BUB1B R = 22)
* BENZTHIAZIDE: P (BUB1B R = 20)
* PEMETREXED: P (BUB1B R = 20)
* PIROXICAM: P (BUB1B R = 20)
* CARTEOLOL: P (BUB1B R = 19)
* DOLUTEGRAVIR: P (BUB1B R = 19)
* MOZAVAPTAN: P (BUB1B R = 19)

Reversal not tied to the mechanism: none

## Limits
A knockout signature is not the partial loss of function of the child; the consensus signatures come from cancer cell lines; the nulls use unrelated gene knockouts; hypotheses only, no dosing, no treatment recommendation.
