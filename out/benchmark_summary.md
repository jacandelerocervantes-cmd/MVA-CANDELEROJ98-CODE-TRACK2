# Multi-disease benchmark (pre-registered, `protocol/addendum_benchmark.md`)

Addendum SHA-256: `348e82a44803a23e1e064219751d4e48579faa48cadcf0f6766902a257e0708b`. 30 diseases selected by the pre-registered rule (Open Targets approved indications of ChEMBL approved drugs; modules from associated targets with score >= 0.5, `clinical_precedence` excluded; first diseases in ascending id order).

## Decision

Rule: mean AUC >= 0.6 and 95% bootstrap CI lower bound > 0.5 at STRING threshold 700.

**Verdict: supported across diseases.** Mean AUC 0.693, 95% CI 0.637 to 0.746.

## Summary across diseases

| Analysis | Diseases | Mean AUC | 95% CI | Median AUC | AUC > 0.5 | Permutation null mean | Permutation p |
|---|---|---|---|---|---|---|---|
| STRING 700 (main) | 30 | 0.693 | 0.637 to 0.746 | 0.724 | 27/30 | 0.532 | 0.0010 |
| STRING 400 | 30 | 0.689 | 0.630 to 0.743 | 0.714 | 25/30 | 0.516 | 0.0010 |
| STRING 900 | 30 | 0.704 | 0.648 to 0.756 | 0.698 | 27/30 | 0.538 | 0.0010 |
| Consensus (mean z-rank of 400/700/900) | 30 | 0.705 | 0.646 to 0.762 | 0.714 | 27/30 | 0.537 | 0.0010 |

AUC: known approved-indication drugs vs all other scorable drugs, score = -z. 0.5 = chance. n/e = not evaluable at that threshold (module < 10 nodes or too few known drugs).

## Per disease

| Disease id | Name | OT targets (score >= 0.5) | Mapped module nodes | Known drugs | AUC 700 | Median pct 700 | AUC 400 | AUC 900 | AUC consensus |
|---|---|---|---|---|---|---|---|---|---|
| EFO_0000544 | infection | 298 | 283 | 7 | 0.407 | 0.492 | 0.340 | 0.504 | 0.461 |
| EFO_0004269 | cardiac arrhythmia | 361 | 271 | 5 | 0.743 | 0.380 | 0.713 | 0.386 | 0.654 |
| EFO_0004272 | anemia (phenotype) | 156 | 131 | 7 | 0.771 | 0.094 | 0.584 | 0.843 | 0.720 |
| EFO_0803536 | multiple sclerosis symptom measurement | 136 | 129 | 23 | 0.679 | 0.306 | 0.645 | 0.688 | 0.704 |
| MONDO_0001056 | gastric cancer | 820 | 684 | 6 | 0.768 | 0.151 | 0.706 | 0.756 | 0.755 |
| MONDO_0001531 | blood coagulation disease | 130 | 111 | 7 | 0.906 | 0.088 | 0.921 | 0.933 | 0.948 |
| MONDO_0002009 | major depressive disorder | 97 | 87 | 22 | 0.654 | 0.251 | 0.644 | 0.661 | 0.655 |
| MONDO_0002050 | depressive disorder | 132 | 118 | 8 | 0.662 | 0.334 | 0.724 | 0.680 | 0.699 |
| MONDO_0002146 | hypogonadism | 91 | 63 | 15 | 0.907 | 0.029 | 0.918 | 0.905 | 0.920 |
| MONDO_0002280 | anemia | 243 | 198 | 22 | 0.840 | 0.101 | 0.797 | 0.853 | 0.831 |
| MONDO_0002406 | dermatitis | 418 | 348 | 7 | 0.504 | 0.503 | 0.470 | 0.593 | 0.520 |
| MONDO_0002525 | inherited lipid metabolism disorder | 383 | 235 | 5 | 0.698 | 0.117 | 0.716 | 0.757 | 0.792 |
| MONDO_0004948 | B-cell chronic lymphocytic leukemia | 123 | 119 | 10 | 0.814 | 0.035 | 0.882 | 0.828 | 0.874 |
| MONDO_0004967 | acute lymphoblastic leukemia | 349 | 327 | 8 | 0.607 | 0.138 | 0.747 | 0.708 | 0.696 |
| MONDO_0004979 | asthma | 171 | 148 | 26 | 0.912 | 0.072 | 0.863 | 0.896 | 0.901 |
| MONDO_0004980 | atopic eczema | 71 | 60 | 22 | 0.761 | 0.212 | 0.759 | 0.780 | 0.766 |
| MONDO_0004981 | atrial fibrillation | 228 | 189 | 5 | 0.487 | 0.640 | 0.498 | 0.341 | 0.415 |
| MONDO_0004989 | breast carcinoma | 1592 | 1325 | 14 | 0.639 | 0.338 | 0.580 | 0.648 | 0.533 |
| MONDO_0004992 | cancer | 8660 | 5683 | 33 | 0.524 | 0.546 | 0.573 | 0.551 | 0.568 |
| MONDO_0004995 | cardiovascular disorder | 2785 | 2101 | 58 | 0.548 | 0.505 | 0.485 | 0.641 | 0.573 |
| MONDO_0005002 | chronic obstructive pulmonary disease | 158 | 135 | 11 | 0.614 | 0.155 | 0.621 | 0.635 | 0.597 |
| MONDO_0005010 | coronary artery disorder | 220 | 184 | 11 | 0.738 | 0.062 | 0.762 | 0.736 | 0.821 |
| MONDO_0005011 | Crohn disease | 21 | 19 | 16 | 0.874 | 0.116 | 0.830 | 0.831 | 0.909 |
| MONDO_0005015 | diabetes mellitus | 758 | 585 | 49 | 0.752 | 0.173 | 0.716 | 0.797 | 0.765 |
| MONDO_0005027 | epilepsy | 531 | 389 | 36 | 0.915 | 0.050 | 0.881 | 0.952 | 0.956 |
| MONDO_0005041 | glaucoma | 167 | 114 | 11 | 0.284 | 0.849 | 0.334 | 0.456 | 0.279 |
| MONDO_0005044 | hypertensive disorder | 324 | 256 | 42 | 0.714 | 0.123 | 0.631 | 0.799 | 0.779 |
| MONDO_0005046 | immune system disorder | 3188 | 2508 | 37 | 0.578 | 0.416 | 0.608 | 0.673 | 0.602 |
| MONDO_0005059 | leukemia | 915 | 821 | 7 | 0.734 | 0.178 | 0.862 | 0.603 | 0.707 |
| MONDO_0005062 | lymphoma | 734 | 672 | 11 | 0.766 | 0.177 | 0.861 | 0.683 | 0.757 |

## Caveats

- Diseases share targets and drugs, so they are not independent; the bootstrap treats them as exchangeable.
- Unlabeled drugs may be true indications; this pulls AUC towards 0.5.
- The id-order rule selects mostly `EFO_` terms; the sample is not representative of all diseases.
- Open Targets literature evidence (`europepmc`) can mention known drug targets; this leak is not removed.
- Proximity does not give the direction of effect. These are hypotheses, not medical advice.
