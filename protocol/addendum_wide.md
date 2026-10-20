# Addendum 6: a wider ring of candidates, for a submission that does not know the child's age

Written on 2026-10-07 **before** the analysis it describes is run. It adds one analysis; no earlier rule changes.

## Why

The plausibility framework (`addendum_plausibility.md`) calls a drug "robustly close" only if it is among the 20 best in at least 7 of 9
configurations. That is a deliberately strict definition and it leaves a short list. Because the age of the child is unknown and an expert panel
reads the output, the list should also contain the drugs that are close but not at the very top, each with its evidence stated, so that nothing
plausible is hidden by a cut-off chosen for convenience. This addendum widens the cut-off in a fixed way and labels the result as a wider ring.

## Rule (fixed before running)

1. **Wide proximity (L1w):** the drug is among the best **50** (instead of 20) in at least 7 of the same 9 configurations (STRING 400, 700, 900 x the
   three weighting schemes), using the ranking with expanded links (variant `e0_e1_e2`).
2. The same gates apply: FDA-eligible (openFDA, currently marketed), not topical-only, not an acute-use anaesthetic or sedative.
3. The same evidence lines L2 (binding to a module protein by curated mechanism or measured activity; predicted counts only as predicted) and L3 (verified
   literature rationale) apply. L4 stays "not assessed".
4. **Tier W:** gates + L1w + L2. **Tier W+:** Tier W + L3. Drugs already in Tier 1 or 2 of the strict framework keep their tier and are not repeated.
5. Antineoplastics (ATC L01 or L02, or mostly neoplasm indications) stay in the tumour-context list.
6. Each drug is reported with what it is approved for today, the gene products it acts on, its closest targets, the three age groups (0-<3, 3-<12, 12 and
   older), and its cautions. The count of drugs per tier is reported whatever it is.

## Complete register

A single table, `out/candidate_register.csv`, lists **every drug that was evaluated and reached at least one evidence line** from any arm of the analysis
(strict tiers, wider ring, extrapolation arm, paediatric search), with which lines it satisfies, which it fails, the FDA status and the age groups. It is the
complete, justified list; the report's main tables are views of it.

## What would count as failure

- The wider ring is dominated by drugs with a broad target profile (more than 20 curated or measured targets) or by network-topology artefacts (for example
  GABA-A modulators): reported with the counts, and those drugs are flagged, not removed.
