# Addendum 9: sensitivity without generic complexes (written 2026-10-07, after the first independent review)

Written **after** an independent review found that the pathway-derived modules contain 30 proteasome subunits plus ADRM1 and 24 tubulins, so that drugs
whose targets are those proteins (proteasome inhibitors, taxanes, colchicine, vinca-type agents and the antibody-drug conjugates with tubulin toxins) have distance 0
to the modules **by construction**. This was checked and is true. The strong proximity of these drugs (z of -7 to -10) is therefore pathway membership and not
independent evidence that they act on the disease mechanism. This addendum does not change any earlier result; it adds a sensitivity analysis and a statement.

## Rule (fixed before running)

1. **Core modules:** each submission-derived module with the proteins whose gene symbol matches `^(PSM[ABCDG][0-9]+|ADRM1|POMP)$` (proteasome subunits, adaptor and
   assembly factor) or `^(TUBA\w+|TUBB\w*|TUBG[12])$` (tubulins) removed. Nothing else is removed.
2. The proximity of every approved drug to the core modules is recomputed with the same null model, iterations, seed and expanded links (`e0_e1_e2`) at STRING 400, 700
   and 900, and the nine-configuration robustness count (top 20) is recomputed.
3. **L1 core:** a drug satisfies L1 core if it is among the 20 best in at least 7 of the 9 configurations with the core modules. L1 core is reported next to L1 for every
   candidate. A candidate that satisfies L1 but not L1 core is labelled `proximity_by_generic_complex_only`.
4. The report states, wherever proximity is used as support, that for proteasome and tubulin drugs it is membership; their rationale then rests on the literature, not on proximity.

## What would count as failure

- If most Tier 1 and Tier 2 drugs fail L1 core, the report says that the proximity evidence is mostly pathway membership of generic complexes, and says so in the summary.
