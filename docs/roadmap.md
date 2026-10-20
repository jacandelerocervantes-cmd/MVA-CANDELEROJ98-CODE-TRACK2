# Roadmap: what a full development would add

This file lists what is **not** in the current version. Nothing here is claimed as done.

**Status after the second round (2026-10-06).** Built and run: literature and trial retrieval with verified identifiers; pharmacovigilance disproportionality; CPIC allele checks (no diplotypes); Gene Ontology direction attempt (no annotation for any candidate); transcriptomic reversal (non-specific; not used); a pre-registered 30-disease validation benchmark; Dockerfile, CI workflow (not executed here), coverage, data manifest, review page, context-of-use and data-protection documents. **Measured, not predicted:** ChEMBL physicochemical properties, measured affinities on the mechanism targets and measured ADME values are reported per drug (`physchem.py`). **Deliberately not built:** docking and any binding-affinity prediction (no independent validation set), doses and paediatric dosing, and interactions with the patient's medications (no medication list): this is a hypothesis-generation hackathon, not a source of medication advice. The items below remain open.

## 1. Safety and patient fit (highest priority)

- **Real evaluation of the patient requirements.** The rules exist as a structure (`gates.py`) but only the nephrocalcinosis caution can be evaluated today. Needed inputs: current age and weight, tumour status and treatments, renal function (eGFR, calcium/urine calcium), liver function, blood counts, allergy history, current medications, ability to swallow. Each would turn a `NOT_ASSESSABLE` into `PASS`/`FAIL` through explicit rules.
- **Pharmacogenomics.** Star-allele calling for CYP2D6, CYP2C19, CYP2C9, CYP3A5, UGT1A1, DPYD, TPMT, NUDT15 with a dedicated caller and CPIC guideline rules; today only rare coding variants are listed.
- **Drug–drug interactions** against the real medication list (CYP3A4 and transporter substrates, inhibitors, inducers).
- **Pediatric dosing and exposure** (allometric scaling, ontogeny of clearance, renal maturation) before any dose statement; excipient exposure per dose (benzyl alcohol, propylene glycol, cyclodextrin vehicles).
- **Richer toxicity evidence.** Pharmacovigilance with proper disproportionality statistics against a background, organ-specific resources (for example liver and kidney injury databases), pediatric-specific safety data, and cross-reactivity between drug classes for allergy.

## 2. Mechanism, direction and evidence

- **Direction of effect.** Signed interaction networks and mechanism rules so that "inhibits a loss-of-function pathway" is separated from "compensates for it"; today every candidate requires manual direction review.
- **Transcriptomic reversal.** Perturbation signatures (LINCS L1000 / Connectivity Map) and BUB1B-depleted cell signatures to score drugs that reverse the disease signature; cell-line dependency data for vulnerability evidence.
- **Affinity layer.** A validated kinase affinity model with an explicit applicability domain, compared against a simple fingerprint baseline with cluster-wise cross-validation (specified, not implemented).
- **Evidence retrieval.** Literature and clinical-trial retrieval with verified citations (hybrid lexical + dense retrieval, citation-existence checks), so that every "evidence for / against" entry carries a verified source rather than only ChEMBL's curated references.

## 3. Validation

- **Many diseases, not one.** A benchmark of diseases with known approved drug–indication pairs and a pre-registered metric (enrichment with permutation tests and bootstrap intervals), instead of a single Fanconi anemia check.
- **Sensitivity to the interactome.** Alternative interactomes (physical networks of the literature, other STRING versions), stable-ranking criteria and consensus ranking.
- **Calibration and uncertainty** for the combined score (for example conformal intervals) and an explicit decision rule for what is a candidate.
- **Independent review** by a pharmacist and a clinical geneticist of the final shortlist.

## 4. Software and operations

- Containerised, pinned environment; continuous integration running the tests and a small reproducible slice of the pipeline; dependency lock with hashes; coverage reporting.
- Versioned snapshots of the third-party data used, with checksums, so that results do not drift with ChEMBL, DailyMed or openFDA releases.
- A thin web interface for reviewing candidates with their evidence, flags and requirement states.
- Packaging so that another disease needs only a mechanism file (already the design), plus automatic generation of the mechanism file from a gene list.

## 5. Governance and ethics

- A formal statement of context of use and decision risk (credibility framework), data-protection review for any patient-level covariate, and a plan for independent expert review before any clinical use.
- No use of the output as a treatment recommendation: the work is hypothesis generation only.
