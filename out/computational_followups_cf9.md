# CF9: independent drug-sensitivity screen (GDSC2) for the line C survivors

Label: exploratory, post hoc, correlational, cancer cell lines; no patient data. Protocol: `protocol/addendum_computational_followups.md (CF9 execution)`.

Screen: GDSC2 fitted dose response, release 8.5 (GDSC2_fitted_dose_response_27Oct23.xlsx); 242036 dose-response rows, 286 drugs. Response: LN_IC50 (lower = more sensitive); AUC reported as secondary. Statistic: line C adjusted partial rank correlation (lineage fixed effects, proliferation covariate, 1,000 within-lineage permutations, seed 20261018). Criterion: original sign and p <= 0.05; not estimable if the drug is absent or fewer than 100 matched lines.

## Cell-line matching coverage

- GDSC lines: 969; matched by Sanger id 966, by COSMIC id 0, by name 1; unmatched 1; ambiguous dropped 1.
- Matched to a DepMap model: 967; retained with expression and lineage: 705.

## Per drug

| Drug | In GDSC2 (name matched) | n lines | Original rho | Adjusted rho (LN_IC50) | p | q (BH over estimable) | AUC rho, p | Verdict |
|---|---|---|---|---|---|---|---|---|
| Lazertinib | absent | 0 | 0.209 | NA | NA | NA | NA | **not estimable** |
| Neratinib | absent | 0 | 0.194 | NA | NA | NA | NA | **not estimable** |
| Paclitaxel | Paclitaxel | 697 | -0.131 | -0.042 | 0.2717 | 0.2717 | 0.025, 0.5235 | **refuted** |
| Vincristine | absent | 0 | -0.164 | NA | NA | NA | NA | **not estimable** |
| Colchicine | absent | 0 | -0.150 | NA | NA | NA | NA | **not estimable** |
| Eribulin | absent | 0 | -0.199 | NA | NA | NA | NA | **not estimable** |

CF9 overall: runnable (1 estimable drug(s)).

## Consequence for each candidate in the plan

- Lazertinib: not estimable; status unchanged: no independent screen covers it.
- Neratinib: not estimable; status unchanged: no independent screen covers it.
- Paclitaxel: refuted; hazard class: stays in the hazard panel and toxicity design H1 (positive/control check only); the line C signal is not replicated in GDSC2 (control not passed); no change to the hazard classification.
- Vincristine: not estimable; hazard class: stays in the hazard panel and toxicity design H1 (positive/control check only); status unchanged: no independent screen covers it.
- Colchicine: not estimable; hazard class: stays in the hazard panel and toxicity design H1 (positive/control check only); status unchanged: no independent screen covers it.
- Eribulin: not estimable; hazard class: stays in the hazard panel and toxicity design H1 (positive/control check only); status unchanged: no independent screen covers it.

## Limits

- Different assay (72 h fitted IC50, Sanger) and a different line panel from PRISM; replication here is partial overlap, not an exact repeat.
- Correlational, post hoc, from cancer cell lines; does not model partial BUB1B loss; a candidate that survives has only not been refuted.
