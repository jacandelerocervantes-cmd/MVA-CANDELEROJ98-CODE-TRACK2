# Line C at drug level, adjusted for lineage and proliferation

Label: exploratory, correlational, from cancer cell lines; does not model partial BUB1B loss; hypothesis signal only. Protocol: `protocol/addendum_line_c_adjusted.md`.

Statistic: partial rank correlation of PRISM viability with the checkpoint score, lineage fixed effects (OncotreeLineage, lineages < 10 lines pooled) and ranked proliferation covariate. Proliferation covariate: mean z of MKI67, PCNA, TOP2A, MCM2, CCNB1 (no measured doubling time on disk). Null: 1000 permutations of the score within lineage, seed 20261018. Survival rule: adjusted BH q < 0.10 and same sign as the original rho.

Cell lines 895, lineages 29; estimable drugs 190.

## Summary

- Original hits (q < 0.10): 23 (23 estimable); surviving the adjustment: **6**; surviving after also removing haematopoietic and lymphoid lines: 6.
- Drugs with adjusted q < 0.10 among all estimable drugs: 7.
- Original versus adjusted effect size over all 190 estimable drugs: Spearman 0.658, Pearson 0.725; over the original hits: Spearman 0.876; median |adjusted| / |original| for the hits 0.610.
- Prednisolone lymphoid suspicion supported: None (None = not testable: prednisolone has 0 haematopoietic or lymphoid lines in its screen data; 109 such lines exist in the expression set).

## Line C candidates

| Drug | Orig rho | Orig q | n | Adj rho | Adj q | Verdict | Within-lineage pooled rho (p RE; k) | Excl. heme/lymphoid rho, q | LOLO sign kept (p<0.05) | LOLO rho range |
|---|---|---|---|---|---|---|---|---|---|---|
| Lazertinib | 0.209 | 0.0190 | 872 | 0.140 | 0.0316 | survives | 0.160 (0.0001; 15) | 0.148, 0.0380 (survives) | 15/15 (15) | 0.118 to 0.152 |
| Neratinib | 0.194 | 0.0190 | 531 | 0.187 | 0.0316 | survives | 0.161 (0.0144; 7) | 0.187, 0.0380 (survives) | 7/7 (7) | 0.167 to 0.209 |
| Prednisolone | 0.160 | 0.0190 | 531 | 0.097 | 0.3227 | does not survive | 0.048 (0.4105; 7) | 0.097, 0.3227 (does not survive) | 7/7 (5) | 0.081 to 0.119 |
| Trametinib | 0.131 | 0.0475 | 550 | 0.087 | 0.3872 | does not survive | 0.073 (0.2002; 7) | 0.087, 0.3723 (does not survive) | 7/7 (2) | 0.081 to 0.107 |
| Cobimetinib | 0.129 | 0.0527 | 546 | 0.061 | 0.5892 | does not survive | 0.041 (0.6752; 7) | 0.061, 0.5892 (does not survive) | 7/7 (1) | 0.045 to 0.094 |

## All original hits (q < 0.10)

| Drug | Orig rho | Orig q | n | Adj rho | Adj q | Verdict | Within-lineage pooled rho (p RE; k) | Excl. heme/lymphoid rho, q | LOLO sign kept (p<0.05) | LOLO rho range |
|---|---|---|---|---|---|---|---|---|---|---|
| Atorvastatin | -0.098 | 0.0475 | 890 | -0.022 | 0.8509 | does not survive | -0.009 (0.8209; 16) | -0.028, 0.7681 (does not survive) | 16/16 (0) | -0.035 to -0.002 |
| Bosutinib | 0.130 | 0.0475 | 552 | 0.106 | 0.2531 | does not survive | 0.027 (0.6369; 7) | 0.106, 0.2531 (does not survive) | 7/7 (7) | 0.098 to 0.141 |
| Cabazitaxel | -0.110 | 0.0475 | 549 | -0.084 | 0.4307 | does not survive | -0.127 (0.0397; 7) | -0.084, 0.4148 (does not survive) | 7/7 (1) | -0.101 to -0.068 |
| Calcitriol | 0.108 | 0.0949 | 553 | 0.065 | 0.5892 | does not survive | 0.027 (0.7028; 7) | 0.065, 0.5892 (does not survive) | 7/7 (1) | 0.049 to 0.093 |
| Cobimetinib | 0.129 | 0.0527 | 546 | 0.061 | 0.5892 | does not survive | 0.041 (0.6752; 7) | 0.061, 0.5892 (does not survive) | 7/7 (1) | 0.045 to 0.094 |
| Colchicine | -0.150 | 0.0316 | 546 | -0.159 | 0.0316 | survives | -0.154 (0.0284; 7) | -0.159, 0.0380 (survives) | 7/7 (7) | -0.176 to -0.144 |
| Dasatinib Anhydrous | 0.123 | 0.0527 | 544 | 0.012 | 0.9397 | does not survive | 0.037 (0.5852; 7) | 0.012, 0.9339 (does not survive) | 6/7 (0) | -0.004 to 0.041 |
| Erdafitinib | -0.115 | 0.0664 | 533 | -0.030 | 0.8509 | does not survive | -0.064 (0.2695; 7) | -0.030, 0.8302 (does not survive) | 7/7 (0) | -0.043 to -0.007 |
| Eribulin | -0.199 | 0.0190 | 872 | -0.115 | 0.0316 | survives | -0.089 (0.0507; 15) | -0.129, 0.0712 (survives) | 15/15 (15) | -0.126 to -0.100 |
| Etoposide | -0.217 | 0.0190 | 549 | -0.132 | 0.1186 | does not survive | -0.091 (0.1503; 7) | -0.132, 0.1055 (does not survive) | 7/7 (7) | -0.151 to -0.125 |
| Fluvastatin | -0.165 | 0.0190 | 531 | -0.123 | 0.1553 | does not survive | -0.101 (0.0898; 7) | -0.123, 0.1553 (does not survive) | 7/7 (7) | -0.145 to -0.091 |
| Lazertinib | 0.209 | 0.0190 | 872 | 0.140 | 0.0316 | survives | 0.160 (0.0001; 15) | 0.148, 0.0380 (survives) | 15/15 (15) | 0.118 to 0.152 |
| Lovastatin | -0.141 | 0.0190 | 549 | -0.027 | 0.8628 | does not survive | 0.024 (0.6808; 7) | -0.027, 0.8510 (does not survive) | 7/7 (0) | -0.039 to -0.017 |
| Methotrexate | -0.150 | 0.0190 | 553 | -0.075 | 0.5152 | does not survive | -0.074 (0.1922; 7) | -0.075, 0.5009 (does not survive) | 7/7 (0) | -0.082 to -0.068 |
| Neratinib | 0.194 | 0.0190 | 531 | 0.187 | 0.0316 | survives | 0.161 (0.0144; 7) | 0.187, 0.0380 (survives) | 7/7 (7) | 0.167 to 0.209 |
| Olaparib | -0.116 | 0.0599 | 547 | -0.082 | 0.4474 | does not survive | -0.054 (0.3504; 7) | -0.082, 0.4320 (does not survive) | 7/7 (1) | -0.091 to -0.072 |
| Paclitaxel | -0.131 | 0.0316 | 546 | -0.141 | 0.0542 | survives | -0.163 (0.0187; 7) | -0.141, 0.0542 (survives) | 7/7 (7) | -0.158 to -0.096 |
| Prednisolone | 0.160 | 0.0190 | 531 | 0.097 | 0.3227 | does not survive | 0.048 (0.4105; 7) | 0.097, 0.3227 (does not survive) | 7/7 (5) | 0.081 to 0.119 |
| Pyrimethamine | -0.109 | 0.0904 | 551 | -0.036 | 0.7630 | does not survive | -0.034 (0.5779; 7) | -0.036, 0.7681 (does not survive) | 7/7 (0) | -0.047 to -0.022 |
| Selumetinib | 0.102 | 0.0990 | 554 | 0.091 | 0.3227 | does not survive | 0.127 (0.1376; 7) | 0.091, 0.3227 (does not survive) | 7/7 (3) | 0.076 to 0.118 |
| Simvastatin | -0.178 | 0.0190 | 548 | -0.072 | 0.5233 | does not survive | -0.063 (0.2696; 7) | -0.072, 0.5095 (does not survive) | 7/7 (0) | -0.082 to -0.060 |
| Trametinib | 0.131 | 0.0475 | 550 | 0.087 | 0.3872 | does not survive | 0.073 (0.2002; 7) | 0.087, 0.3723 (does not survive) | 7/7 (2) | 0.081 to 0.107 |
| Vincristine | -0.164 | 0.0190 | 546 | -0.158 | 0.0316 | survives | -0.183 (0.0013; 7) | -0.158, 0.0380 (survives) | 7/7 (7) | -0.171 to -0.141 |

## Verdicts

- Lazertinib: survives.
- Neratinib: survives.
- Prednisolone: does not survive.
- Trametinib: does not survive.
- Cobimetinib: does not survive.

## Notes

- LOLO: lineages with at least 30 lines are omitted one at a time; the column gives the omissions that keep the original sign (and, in brackets, also have p < 0.05).
- Within-lineage: DerSimonian-Laird random effects over lineages with at least 30 lines (Fisher z of the proliferation-partialled Spearman correlation); reported, not part of the survival rule.
- Limits: correlational, from cancer cell lines; does not model partial BUB1B loss; the gene-set proliferation proxy leaves possible residual confounding. Hypothesis signal only.
