# Addendum 8: structure-based docking (E3), with its own validation

Written on 2026-10-07 **before** any structure is downloaded or any docking is run. It completes the item "E3, docking" that `addendum_dti.md` left out.

## Why

Measured binding (E1) and ligand-similarity prediction (E2) only see what is already known or resembles something known. Docking can propose a **new** binder for
a protein with a crystal structure. It is only used here if, for that protein, it first reproduces known facts (it re-docks the crystal ligand, and it ranks known
binders above other drugs). A docking score that has not passed this check is not used, and that is reported.

## Targets and structures (selected by a rule, not by hand)

- Candidate proteins: ligandable proteins of the modules and their immediate network neighbours that are kinases or other enzymes with a known drug-like inhibitor:
  CDK1 (P06493), AURKA (O14965), AURKB (Q96GD4), PLK1 (P53350), PLK4 (O00444), TTK (P33981) and NEK7 (Q8TDX7). The list is fixed here.
- For each protein the structure is the **best-resolution human X-ray entry in the RCSB PDB (CC0) that has a drug-like ligand** (heavy-atom count 15 to 60,
  not a metal, solvent or buffer) bound in the active site, with resolution <= 2.8 angstrom. The chain and the ligand are chosen by the same rule (the chain of the
  target protein in contact with the ligand). Waters, other ligands and alternative locations are removed.
- Search box: the centre of the co-crystal ligand, 22 angstrom per side.

## Method

- Program: AutoDock Vina 1.2.7 (official Scripps release; the SHA-256 of the executable used is recorded). Exhaustiveness 8, 9 poses, seed fixed per ligand and
  target. Receptors and ligands are converted to PDBQT with Open Babel (hydrogens at pH 7.4, Gasteiger charges); ligands start from the ChEMBL SMILES (RDKit 3D embedding,
  fixed seed). Covalent binders (for example proteasome inhibitors) cannot be docked this way and are not attempted.
- Ligand set: approved small molecules with at most 50 heavy atoms and at most 15 rotatable bonds.

## Validation required for each target before its scores are used (fixed now)

1. **Redocking:** the co-crystal ligand, prepared the same way, is docked back; the target passes this check if the best pose is within **2.0 angstrom RMSD** of the crystal pose.
2. **Enrichment:** actives = approved drugs with a curated or measured (pChEMBL >= 6) link to that protein; decoys = a fixed random sample of 200 other approved small molecules without such a link.
   The target passes if the **ROC AUC of the Vina affinity (more negative = better) is >= 0.70** and there are at least **4 actives**.
3. A target that fails either check is **not used**; its result is reported.

## Use of the scores

- For each validated target, the hit threshold is the least negative affinity at which the precision on the validation set is >= 0.5. A drug whose best affinity on that target is at or
  below the threshold gets a **predicted docking link (E3)** to the protein, labelled as such (a hypothesis about binding, not a measurement). Links already known (E0, E1) are not repeated.
- E3 links are added to the proximity analysis exactly as E2 was (variant `e0_e1_e2_e3`), and the plausibility tiers are recomputed. E3 never counts as L2 on its own: it is shown as `L2_predicted_only`.
- Separately, **the top candidates** (Tier 1 and 2, and the wider ring) are docked against every validated target and the poses and scores are reported next to them, as context.

## What would count as failure

- No target passes validation: docking is reported as not usable and no E3 link is created.
- A high fraction of hits are very large or very flexible molecules (size bias of the scoring function): reported with the heavy-atom counts.
