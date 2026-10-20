# Frozen validation protocol — Track 2

Freeze date: 2026-10-06, **before** running the pipeline on real data (up to this point only tests with synthetic graphs exist). This file is hashed (`protocol/frozen_protocol.sha256`) and is not modified after seeing results. Only code bug fixes that change neither parameters nor decision rules are admitted; any other change is recorded separately as a deviation, with its reason.

## Frozen parameters

- Interactome: STRING v12.0 physical, combined score ≥ 700 (main); sensitivity with 400 and 900.
- Proximity: d_closest, 1000 random iterations, degree matching with bins of at least 100 nodes, base seed 20261006.
- Main proximity threshold: z ≤ −2. In addition, z ≤ −0.15 and z < −4 are reported.
- BUB1B module: seeds from `config/mechanisms/bub1b_sac_lof.json` (Reactome R-HSA-69618 and MVA/checkpoint genes).
- Universe: ChEMBL, `max_phase = 4`, non-withdrawn parent molecules.

## Tests and decision rules (fixed before seeing results)

**P1 — Pathway-level recovery (leave one out).** For each seed gene g of the BUB1B module, the module is rebuilt without g and the z of {g} is computed. Negative controls: 1000 random genes outside the module, matched by degree. Metrics: fraction of seeds with z ≤ −2 and AUC of seeds versus negatives, with a 95 % bootstrap interval (1000 resamples). **Rule:** the method is declared "supported at pathway level" only if AUC ≥ 0.70 and the lower bound of the interval is > 0.50; otherwise it is declared "not supported" and the drug ranking is presented only as exploratory.

**P2 — Drug controls in a second disease (Fanconi anemia; exploratory).** Module: Reactome "Fanconi Anemia Pathway" (R-HSA-6783310). Positive controls: metformin and oxymetholone (percentile of their z among all scorable approved drugs). Conceptual negative controls: cisplatin and mitomycin C (if they have no protein target in ChEMBL they are reported as "not scorable"). **Rule:** the percentiles are reported as they are; nothing is adjusted if the positive controls do not stand out.

**P3 — Robustness.** (a) STRING thresholds 400/700/900; (b) randomly dropping 20 % of the seeds 20 times; (c) three different random seeds. Metric: overlap (Jaccard) of the top 20 drugs and Spearman rank correlation. **Rule:** they are reported; a top 20 with mean Jaccard < 0.5 is declared unstable.

**P4 — Degree-bias control.** Spearman correlation between z and the mean degree of the targets; it is reported (degree normalisation should leave it close to 0).

## What is not done

- No thresholds, seeds or parameters are adjusted after seeing the results.
- No label taken from an external evaluator is used.
- Negative results are reported in the same detail as positive ones.
