# Addendum 4: new repositioning hypotheses from expanded drug-target links

Written on 2026-10-06 **before** the analysis it describes is run. It does not change the frozen protocol or the earlier addenda.

## Why

The ranking so far uses only the curated mechanism of each drug in ChEMBL (about 1,040 of 2,705 approved drugs have one that maps to a
protein). That cannot find a **new use**: a drug that is not known to act on the disease pathway but binds one of its proteins (the
sildenafil case: designed for one thing, useful for another). Repositioning needs drug-target links beyond the curated table. This addendum adds
them, in three tiers, and recomputes the proximity with the expanded links.

## Sources of links (fixed before running)

- **E0, curated** (current): ChEMBL mechanism table.
- **E1, measured:** any ChEMBL bioactivity of an approved small-molecule drug on a human single-protein target with pChEMBL >= 6.0 (activity at or
  below 1 µM; IC50, Ki, Kd, EC50, Potency) from a **binding assay** (ChEMBL assay type B) with exact relation (=). These are measurements made by someone,
  not predictions. (The ChEMBL activity resource does not carry the assay confidence score, so binding assays on a single-protein target stand in for it.)
- **E2, predicted by ligand similarity:** for each protein of the four modules that has at least 10 reference ligands (pChEMBL >= 6.0, up to the 300
  strongest), the score of a drug is the maximum Tanimoto similarity (Morgan radius 2, 2048 bits) to a reference ligand, excluding the drug itself.
  A link is predicted when the score is >= tau. **tau is the smallest value for which a retrospective test has precision >= 0.5**: known E0 and E1 links of
  approved drugs are the positives, each tested with the drug removed from the reference set, and the negatives are all the other drug-target pairs of the
  same targets. If no tau reaches precision 0.5, E2 is not used and that is reported.
- **E3, docking: not in this addendum.** Docking needs a docking program and PDB structures that are not installed; it is a separate step that would
  need its own validation (redocking of known ligands and enrichment against decoys) before any score is used.

Target universe for E1 and E2: the union of the proteins of the four submission-derived modules (`config/mechanisms_from_submission/`).

## Analysis

1. Recompute the closest-distance proximity z of every approved drug to each module with the target set E0 + E1, and with E0 + E1 + E2 (if E2 is valid),
   using the same null model, iterations, seed and STRING thresholds as before. Drugs not scorable before (no curated mechanism) become scorable.
2. A drug is a **new hypothesis** if its z <= -2 in at least one module with the expanded set and either it was not scorable before or its curated-only
   z was > -2. A drug whose z <= -2 already with E0 alone is **not new**.
3. Robustness: the nine configurations of `addendum_robustness.md` are applied to the expanded ranking (rank within the best 20).
4. Each new hypothesis reports its supporting links (target, source E1 or E2, pChEMBL or similarity), the age-band label status (reviewed where it was
   read), the ATC class and the mechanism-class literature note.

## Classes shown separately (decided now, not after seeing results)

- **Antineoplastic (ATC L01 or L02):** tumour-context hypotheses only.
- **Anaesthetics, sedatives and hypnotics (ATC N01, N05C):** excluded from proposals (acute-use agents; reported only if they pass).
- **Everything else:** hypotheses for the constitutional defect and its consequences.
Whether any drug is appropriate for this child is a clinical decision that this project does not make.

## What would count as failure

- E2 below precision 0.5: reported; only E1 is used.
- No new hypothesis robust in 7 of 9 configurations: reported as a negative result and the list stays exploratory.
- Expanded links produce a ranking dominated by promiscuous drugs (kinase inhibitors with measured activity on hundreds of proteins): reported with the number
  of module proteins each drug touches, and promiscuity (more than 10 module targets) is flagged.
