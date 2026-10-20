# Addendum 5: FDA scope and the repurposing-plausibility framework

Written on 2026-10-06 **before** the analysis it describes is run. It sits on top of `addendum_dti.md` (expanded drug-target links) and does not
change any earlier rule.

## Scope: FDA-approved drugs only

A drug is eligible for proposal only if openFDA (Drugs@FDA) has a record whose active ingredient matches the drug and at least one of its
products has the marketing status `Prescription` (currently marketed). Discontinued products and drugs approved only by other regulators are not
eligible. The check is made on the active-ingredient name (the salt form is tried first, then the first word of the name). A drug with no match is
reported as `not_fda_eligible`; a drug that cannot be checked because the service failed is reported as an error, never as eligible.

## What "plausible to reuse" means (the framework)

Plausibility is not one score. It is the number of **independent, public lines of evidence** that agree, plus **feasibility gates** that any real use
would have to pass. Nothing is hand-weighted.

Evidence lines (each yes or no):
- **L1, network proximity robust:** the drug is in the best 20 in at least 7 of the 9 configurations (3 STRING thresholds x 3 weighting schemes) of the
  ranking with expanded links (`addendum_dti.md`, `addendum_robustness.md`).
- **L2, binding to a module protein:** a curated mechanism target, or a measured activity (pChEMBL >= 6, ChEMBL binding assay), on at least one protein of the
  modules. A similarity-predicted link counts as `L2_predicted` only if the retrospective validation of `addendum_dti.md` passed, and is labelled as predicted.
- **L3, mechanistic rationale:** a curated, PMID-verified literature class applies to the drug's mechanism with a direction hint other than
  `unknown_possibly_unfavourable_in_weak_checkpoint_context`.
- **L4, independent functional data** (cell-line dependency or drug-sensitivity screens, perturbation signatures): **not assessed in this version.** It is
  reported as not assessed for every drug, never as absent.

Feasibility gates (all must pass to be proposed):
- FDA-eligible (above).
- Not topical-only (ChEMBL route flags are incomplete, so only a drug whose sole recorded route is topical is excluded; this reading of the gate was fixed after finding that bortezomib has no recorded route).
- Class not in {anaesthetics, sedatives and hypnotics (ATC N01, N05C)}.
- Antineoplastic agents (ATC L01, L02) are listed in a separate tumour-context list, never in the constitutional-defect list.
- A boxed warning, promiscuity (more than 10 module proteins linked) and excipient alerts are **shown, not used to exclude**, so that a reader sees them.

## Tiers

- **Tier 1:** gates + L1 + L2 + L3.
- **Tier 2:** gates + L1 + L2.
- **Tier 3:** gates + L1 or L2 (listed, not proposed).
Counts per tier are reported whatever they are, including zero.

## Age

Tiers do not depend on age. Every proposed drug is shown for each age band (neonate, infant, child, adolescent, adult) with its US label status (reviewed
when the label was read, heuristic otherwise). If the age is given, bands that the label does not support are marked and can be filtered. The project does
not give doses.

## Why this is the contribution

The framework needs only a gene list (from any variant-prediction output), public resources and, optionally, an age. Replacing the submission file
reruns the whole pipeline, which is what makes the approach scalable; combining measured binding beyond the curated mechanism, network proximity,
direction-of-effect literature and FDA and age feasibility in one pre-registered rule is what is new.

## What would count as failure

- No Tier 1 or Tier 2 drug: reported as a negative result; the list stays exploratory.
- L3 satisfied for no drug: reported (direction of effect is then the main open gap, not hidden).

## Status notes (2026-10-07, after the first independent review)

- **L3** excludes the class `mva_tumour_case_reports` (case reports are not a mechanistic rationale) as well as `microtubule_targeting`; this is how it was run. Read literally without that exclusion, Tier 1 would have 13 drugs.
- **L2 by measured activity** needs at least two activity records; a single record is shown but does not count.
- **L4** is now assessed (`addendum_prism.md`, with its amendment) and does not change the tiers.
- **L1 core** (`addendum_core.md`) is reported next to L1: for proteasome and tubulin drugs the proximity is membership of generic complexes.
