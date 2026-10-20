# Transcriptomic reversal layer (SigCom LINCS, L1000)

> Hypothesis-generating evidence. Not medical advice.

## Caveats

- Signatures come from **cancer and immortalised cell lines** (L1000), not from the proband's tissues.
- A CRISPR knockout of BUB1B in a cell line is only a **proxy** for the constitutional *hypomorphic* loss in the proband (partial loss, all tissues, developmental context). The other checkpoint genes are pathway-level proxies only.
- A drug whose signature moves against the knockout signature is **not** shown to correct mosaic variegated aneuploidy; the score is hypothesis-generating, carries no direction of clinical benefit, and is not a dose or treatment statement.
- Consensus chemical signatures merge cell lines, doses and times; cell-line detail comes from a second, descriptive query of the individual-signature library restricted to its top reversers (absence there is not absence of effect).
- Drugs are matched to LINCS perturbagens by cleaned name only; ambiguous matches are flagged `name_match: ambiguous`.

## Method

- Query: top 100 up and 100 down genes of the LINCS CRISPR-KO consensus signature of BUB1B (primary) and of BUB1, TTK, ZWINT, CENPE (proxies).
- Ranking: SigCom LINCS `POST /data-api/api/v1/enrich/ranktwosided` against `l1000_mean_cp` (chemical consensus signatures). Score = best reversal z-sum (-z-sum, larger = stronger reversal) among the signatures of a drug.
- Null: 30 random up/down gene sets of the same size (seed 20260101) through the same pipeline; empirical p pools 34080 drug-level null scores; FDR is Benjamini-Hochberg across scored drugs.
- Direction: `reverses` if FDR <= 0.05 for reversal, `mimics` if FDR <= 0.05 for mimicry, else `none`.

- Structured control: the same pipeline on 20 unrelated CRISPR-KO consensus signatures (BRDN0001147488, BRDN0001148320, BRDN0001482644, BRDN0001494889, BRDN0003228062, BRDN0003578361, DAPK2, GART, HTR2B, IL17RB, LPL, LYRM7, MPZL3, MTA3, PSG1, RPS6KC1, SAMD8, SMAD9, TBC1D7, TXN); `ko_control.p` = (1 + controls scoring >= BUB1B) / (1 + 20). A drug that reverses most knockout signatures is a non-specific reverser.

## Results

- Approved drugs in the universe: 2705; with a LINCS chemical consensus signature (scored): 1136; the rest are `available: false`.
- `reverses` at FDR <= 0.05: 433; `mimics`: 139.
- Null reversal score quantiles (random sets): {'p50': 0.0288, 'p95': 2.3319, 'p99': 3.2427, 'p99.9': 4.2326}; maximum 6.5572.
- Drugs `reverses` AND ko_control p <= 0.05 (specific to BUB1B-like signature): 97.
- Random sets called `reverses` per set (leave-one-out, same rule): mean 0.1 (per set: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]).

## Top 25 reversers

| Drug | Score | p | FDR | Direction | Signatures | Cell lines (top-reverser layer) | In candidates | Proxy sets reversing | p vs unrelated KO controls |
|---|---|---|---|---|---|---|---|---|---|
| SULFACETAMIDE | 17.2714 | 2.9e-05 | 0.000214 | reverses | 1 | A375 | False | 2/4 | 0.0476 |
| CEFOPERAZONE | 14.0797 | 2.9e-05 | 0.000214 | reverses | 1 | HA1E, HCC515, MCF10A | False | 3/4 | 0.0476 |
| ENZALUTAMIDE | 14.0555 | 2.9e-05 | 0.000214 | reverses | 1 | MCLF035SZ, MCLF117SZ, MCLF126CN | False | 3/4 | 0.0476 |
| NIZATIDINE | 13.5293 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 2/4 | 0.0476 |
| TRILOSTANE | 13.2704 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 3/4 | 0.0476 |
| TOFACITINIB | 13.2297 | 2.9e-05 | 0.000214 | reverses | 1 | HELA, HT29, MCF7, THP1, VCAP, YAPC | False | 2/4 | 0.0476 |
| MESTRANOL | 12.6144 | 2.9e-05 | 0.000214 | reverses | 1 | A549, HT29, U2OS | False | 3/4 | 0.0476 |
| TEMOZOLOMIDE | 12.3169 | 2.9e-05 | 0.000214 | reverses | 1 | CAL29, HCC515, XC.P931 | False | 2/4 | 0.0476 |
| AMIKACIN | 12.2429 | 2.9e-05 | 0.000214 | reverses | 1 | HA1E | False | 2/4 | 0.0476 |
| HYDROXYAMPHETAMINE | 12.0297 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 4/4 | 0.0476 |
| TICAGRELOR | 11.7566 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 2/4 | 0.0952 |
| MECLOFENAMIC ACID | 11.3663 | 2.9e-05 | 0.000214 | reverses | 1 | HA1E, YAPC | False | 2/4 | 0.0476 |
| ARGATROBAN | 11.3634 | 2.9e-05 | 0.000214 | reverses | 1 | HA1E | False | 3/4 | 0.0476 |
| PRAZIQUANTEL | 11.2399 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 3/4 | 0.0476 |
| PREDNICARBATE | 11.216 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 1/4 | 0.0476 |
| PARECOXIB | 11.2076 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 2/4 | 0.0476 |
| CEFIXIME | 11.1653 | 2.9e-05 | 0.000214 | reverses | 1 | HCC515, PC3 | False | 3/4 | 0.0476 |
| PRAVASTATIN | 11.1455 | 2.9e-05 | 0.000214 | reverses | 1 | JURKAT, MDAMB231 | False | 2/4 | 0.0952 |
| CARBAMAZEPINE | 11.0585 | 2.9e-05 | 0.000214 | reverses | 1 | A549, HA1E, PC3 | False | 3/4 | 0.0476 |
| ORPHENADRINE | 10.9429 | 2.9e-05 | 0.000214 | reverses | 1 | PC3 | False | 3/4 | 0.0952 |
| NAPROXEN | 10.8383 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 3/4 | 0.0952 |
| ESTRIOL | 10.8159 | 2.9e-05 | 0.000214 | reverses | 1 | A375 | False | 2/4 | 0.0476 |
| IDELALISIB | 10.8153 | 2.9e-05 | 0.000214 | reverses | 1 | K562, OCILY3 | False | 2/4 | 0.0476 |
| IOXAGLIC ACID | 10.8072 | 2.9e-05 | 0.000214 | reverses | 1 | A549 | False | 1/4 | 0.0476 |
| GEMIFLOXACIN | 10.7174 | 2.9e-05 | 0.000214 | reverses | 1 | - | False | 3/4 | 0.0952 |

Full table: `reversal_scores.json`; gene lists and parameters: `reversal_meta.json`.
