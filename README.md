# HSelith — Track 2 (Drug Repurposing), MVA Hackathon 2026

Report and reproducible framework that answer the question of Track 2: **which already approved drugs are plausible candidates for further investigation, given the mechanism disrupted in a proband with Mosaic Variegated Aneuploidy (MVA1, *BUB1B*)?** The report goes from the mechanism (five consequences of reduced BubR1 function, tied to the child's features) to hypotheses, each with verified literature for and against, an evidence level (L0 to L3), the paediatric label read from the label text and the experiment that would refute it. **Nine drugs in three groups (A constitutional: sirolimus, everolimus, metformin, sodium phenylbutyrate; B tumour context: trametinib, selumetinib; C hazards to check, not proposed: bortezomib, chloroquine, bosutinib)** (best evidence level 1; no approved drug reaches level 2 or 3). The selection rules were declared after the computational results (amendment E16, `protocol/addendum_e16.md`). Four computational lines of evidence (a knowledge graph, network proximity, cell-line vulnerability and a direction check), specified with their validations before the final runs, are a filter on the hypotheses: **no drug meets the rule of two usable lines written in the protocol**, and the report shows why with the controls that produced that result. The work also delivers a class-level hypothesis, a caution list (from the literature and from verbatim label sentences), label support by age group (the child's age is unknown) and a pipeline whose first stages (gene modules and network proximity, stages 1 to 16) rerun for any Track 1 submission; the hypothesis layer (lines A and C, the five consequences and the selection of the drugs) is written for this case and would have to be redone for another.

> **Hypotheses, not evidence and not treatment recommendations.** The drugs listed must not be used on the basis of this repository.

This repository contains **Track 2 only**: code, tests, configuration, a frozen validation protocol, the computed outputs and the written report. It contains no patient reads or VCFs, no variant coordinates (the Track 1 submission CSV is not included and the derived files do not carry them) and no third-party datasets (they are downloaded at run time).

**Track 2 submission file:** `submissions/candeleroj98_candidates_track2.csv` lists the 28 distinct candidate medicines (30 entries across the evidence lines), with the lines that propose each one, the label reading by age band (0-4, 5-11, 12-17 years), safety flags and the multi-parameter cross-check rank. It is generated from `out/candidates_by_line_mva_v2.json` and contains no variant coordinates. Every entry is a hypothesis for follow-up.

| | |
|---|---|
| Hugging Face user | `candeleroj98` |
| Report | [`submissions/candeleroj98_track2_report.md`](submissions/candeleroj98_track2_report.md) (and `.pdf`) |
| Methods form | [`submissions/candeleroj98_methods_description_track2.md`](submissions/candeleroj98_methods_description_track2.md) and the official form `submissions/candeleroj98_methods_description_form_track2.xlsx` |
| Final results | [`out/four_lines_review2.json`](out/four_lines_review2.json), [`out/class_analysis.json`](out/class_analysis.json), [`out/caution_list.json`](out/caution_list.json), [`out/age_label_table.json`](out/age_label_table.json), [`out/forking_ledger.json`](out/forking_ledger.json), [`out/proposals_review2.json`](out/proposals_review2.json), [`out/literature_hypotheses/`](out/literature_hypotheses/) |
| Earlier (addendum 10) outputs, kept as sensitivity analyses | [`out/candidates.json`](out/candidates.json), [`out/four_lines.json`](out/four_lines.json) (its single candidate, fostamatinib, is superseded: see the ledger) |
| Protocol | [`protocol/frozen_protocol.md`](protocol/frozen_protocol.md), [`protocol/addendum_four_lines.md`](protocol/addendum_four_lines.md) and the approved [`protocol/addendum_review2.md`](protocol/addendum_review2.md) |

## Headline results (all computed; see the report)

- **Literature-anchored hypotheses, declared post hoc and weakly supported (best level 1):** group A constitutional: sirolimus, everolimus, metformin, sodium phenylbutyrate; group B tumour context: trametinib, selumetinib; group C warnings, **not proposed**: bortezomib, chloroquine, bosutinib. A second drug-agnostic search round (256 queries) found no approved drug that rescues a BubR1-deficient phenotype. Computational candidates by the rule written in the protocol: **none**, under either reading of line B (B2 counted as in the protocol, or held exploratory). The overlap between lines is zero in every pair.
- **Line A (graph)** is valid but modest: AUC 0.802, excess over a degree-only baseline 0.022. Its patient ranking is generic: 43 of the 46 drugs of its top 50 also reach the top 50 of at least 5 % of 500 degree-matched random queries; fostamatinib, the single candidate of the first analysis, is in the top 5 of 98.8 % of them.
- **Line B1 (phenotype module)** does not beat the degree baseline once the causal genes are removed (excess +0.095, interval -0.083 to +0.267). **B2 (BUB1B module)** passes the validation threshold of the protocol (AUC 0.777, 15 diseases) but the check that removes the drugs acting on the causal genes leaves 11 diseases (not testable).
- **Line C (cell lines)**, with TP53, MSI and chromosomal instability in the model and 1,000 null sets: 5 hits against a 95th percentile of 17.0 (22 against 23.05 before). Six MDM2 inhibitors among the first hits were the TP53 signal. Not usable.
- **Class level:** one class passes in two usable lines (ATC L01EX, multi-target kinase inhibitors): evidence level 1, a tumour-context class mostly labelled for adults.
- **Caution list** (tubulin binders, genotoxic chemotherapy, proteasome inhibitors, growth hormone with a tumour history) and **label support by age group** for 70 drugs. Every quoted label sentence (551 across the caution list, the age table and the label readings) is checked to appear in the label it cites (`track2.verify_quotes`: 551 found, 0 not found). Every quoted literature sentence is re-fetched from Europe PMC (`track2.verify_literature`: 355 verified, 6 preprint sentences not checkable and not cited, 2 not found and therefore excluded). The paediatric status of each of the 71 drugs was read from its label and recorded in `config/pediatric_label_review.json` (a reading by the author's assistant, not a pharmacist review); 2 stay unresolved (`verify`: busulfan, cyclosporine).
- All patient requirements that need data we do not have (age, tumour status, organ function, blood counts, allergies, treatments) stay `NOT_ASSESSABLE`.

All 71 drugs of the study are FDA-approved and marketed (not experimental): `python -m track2.fda_approval_check` confirms each one in Drugs@FDA and writes `out/fda_approval_check.json`.

## Lines of evidence and amendments (the final analysis)

Design and amendments are in `protocol/` (addendum 10: four lines; addendum 11: amendments E1 to E15 after two independent reviews by AI agents, plus E16 in `protocol/addendum_e16.md`, with the two conditions that the work is repurposing of approved drugs only and that the child's age is unknown). The numerical thresholds in the code are those of the protocol files. The amendments are post hoc with respect to the results the reviews read, and the protocol rules were written before the runs reported; there is no git history and no external time stamp, so the order of events rests on the author's statement. `out/forking_ledger.json` lists every variant that was run and which is the primary one.

| Module (`src/track2/`) | Role |
|---|---|
| `fetch_sources.py` | downloads HPO, MONDO, PrimeKG, DepMap CRISPR and covariates one source at a time and records SHA-256 in `data/sources/record.json` |
| `kg_validation.py`, `kg_extra.py`, `kg_query.py --e5`, `kg_nullquery.py --e5` | line A: validation, excess over degree, patient query, specificity |
| `phenotype_module.py`, `phenotype_validation.py --leave-causal`, `phenotype_proximity.py --strict-links`, `module_rule_validation.py`, `core_n2.py`, `track1_scenarios.py` | line B (B1 and B2) |
| `vulnerability_scan.py --covariates` | line C with TP53, MSI, CIN and 1,000 null sets |
| `signed_effect.py` | line D (exploratory) |
| `four_lines_review2.py`, `class_analysis.py`, `caution_list.py`, `age_label_table.py`, `forking_ledger.py` | combination and the final outputs |
| `verify_literature.py`, `literature_summary.py`, `report_candidates.py`, `overlap_check.py` | mechanism-guided proposals (E16): Europe PMC re-check of every sentence, selection rules, report sections, expected-overlap control |
| `build_report_review2.py`, `build_methods_review2.py`, `fill_methods_form.py`, `md_to_pdf.py` | report, methods, official form and PDF, all generated from `out/` |

## The earlier pipeline (addendum 10)

Still used for the drug universe, the network, the drug-target links, the labels and the sensitivity analyses; its own ranking is no longer the headline.

1. **Universe** (`chembl_universe.py`): approved parent molecules, curated mechanisms and UniProt targets from ChEMBL.
2. **Network** (`string_net.py`): STRING v12.0 physical interactions, combined score ≥ 700, largest connected component.
3. **Modules** (`submission_input.py`, `auto_modules.py`, `modules.py`): the Track 1 submission CSV gives the genes and weights (best `epcr`); ClinVar is queried for its SNVs; each gene becomes a protein module from its Reactome pathways (STRING neighbours as fallback; minimum size 10). Rule fixed in `protocol/addendum_submission.md`. The hand-defined modules in `config/mechanisms/` are used for validation and as a comparison.
4. **Proximity** (`proximity.py`): d_closest and a z-score against 1000 degree-matched random draws (Guney 2016; Cheng 2018; reimplemented because the original code has no licence).
5. **Scoring** (`scoring.py`): ranks, ChEMBL action, risk flags, `direction_review = REQUIRES_REVIEW` (proximity gives no direction of effect).
6. **Labels, toxicity, allergy and requirements** (`labels.py`, `toxicity.py`, `gates.py`): DailyMed label text → age-band, warning, hypersensitivity and excipient flags (heuristic, verify by hand); an openFDA adverse-event summary (descriptive only, not an incidence); eight requirements with states `PASS`/`FAIL`/`NOT_ASSESSABLE` (a label warning about nephrocalcinosis fails the documented nephrocalcinosis requirement).
7. **Records** (`build_candidates.py`): validated against `schemas/candidate_drug.schema.json`; every PubMed reference is checked to exist.
8. **Validation** (`validation.py`): the frozen protocol (P1–P4).
9. **Genome covariates** (`genome_covariates.py`, optional, needs the patient VCF): sex, aneuploidy signal, variants in pharmacogenes and tumour-predisposition genes (aggregates only are published).
10. **Report** (`build_report.py`, the addendum-10 generator, no longer run by `run_all`): superseded by `build_report_review2.py`.

Additional evidence layers, none of which changes the ranking (`integrate.py` merges them into `out/candidates_enriched.*`):

- `evidence.py`: literature and trials retrieved by keyword (PubMed, ClinicalTrials.gov), every PMID verified to exist.
- `faers_disproportionality.py`: PRR/ROR from FDA adverse-event counts (spontaneous reports; not incidences).
- `pgx.py`: CPIC allele checks for eight genes (no diplotype or phenotype is asserted; needs the private VCF).
- `direction.py`: Gene Ontology orientation of drug targets for the checkpoint; none of the candidates could be oriented.
- `reversal.py`: LINCS L1000 transcriptomic reversal; the result is non-specific and is reported but not used.
- `benchmark.py`: validation over 30 diseases fixed in a protocol file from Open Targets (mean AUC 0.69, biased and optimistic).
- `snapshot.py`, `review_page.py` (`out/review.html`), `delete_data.py`: data manifest, static review page, and the 30-day deletion tool.

Measured ChEMBL binding is used for drug-target links (at least two activity records); docking was run with its own validation and no target passed, so no docking score is used. Doses, paediatric dosing and interactions with the patient's own medications are out of scope for a hypothesis-generation hackathon. See `docs/roadmap.md`.

## Repository layout

```
.
|-- README.md, THIRD_PARTY_LICENSES.md
|-- pytest.ini, requirements-track2.txt / .lock, .python-version, .gitattributes
|-- config/         mechanisms, patient.json, mechanism_screen.json, pediatric_label_review.json, ai_declaration.md
|-- (inputs)        the Track 1 submission CSV is NOT included; see below
|-- protocol/       frozen protocol, addenda (four lines, review2), deviations.md
|-- src/track2/     the pipeline and the final analysis modules
|-- tests/          586 tests collected
|-- out/            computed outputs, forking_ledger.json, checksums.sha256
|-- docs/           data protection, context of use, coverage, roadmap, mechanism screen
|-- notebooks/      Colab notebook
|-- submissions/    report (md, pdf), methods description (md, official xlsx)
`-- data/           third-party downloads and caches (ignored by git; never published)
```

## Tested environment and versions

Executed by the author on **Windows 11 (10.0.26200), Python 3.10.11, pip 23.0.1, setuptools 78.1.0**; 586 tests are collected; 581 passed in the author's environment with the downloaded sources, and in a clean clone without them 561 pass and 25 are skipped. The code uses `pathlib` and writes files with explicit LF.

**Also run on Google Colab (Linux 6.6, glibc 2.39, Python 3.13.16, CPU runtime, 2026-10-07)** with `notebooks/MVA_Hackathon_2026_Track2.ipynb` (all 12 code cells, 59 s, `USE_DRIVE = False`): `requirements-track2.lock` installs without changes, the 169 tests pass, the four sources downloaded again (15 files compared) have the SHA-256 recorded in `data/sources/record.json`, and two results of the final analysis (`out/phenotype_module_e5.json`, `out/kg_query_e5.json`) are regenerated **byte for byte** and match `out/checksums.sha256`. The first attempt showed that the information-content scores of the phenotype module differ in the last digit (about 1e-15, `math.log` differs between operating systems); the written scores are now rounded to 10 decimals, with no change to any gene, variant or ranking.

**Scope of the Colab run.** It was made on 2026-10-07 with the code of that date (169 tests). Later changes (the literature layer, the quote and literature verifiers, the article builder, the safety statements and the documentation) were tested locally only, on Windows with Python 3.10 (171 tests and the recorded checksums); they have not been run on Colab or on Linux. Not tested: the Drive mount (it needs an interactive authorisation), stages 17 to 21 of `run_all` on Colab (several hours), other operating systems, and byte-level reproduction of the remaining deterministic outputs on a different platform: other floating-point outputs may differ in the last digits across systems, so on another platform compare values, not hashes.

| Package | Version | Role |
|---|---|---|
| numpy | 2.2.6 | arrays, random draws |
| scipy | 1.15.3 | sparse graphs, shortest paths, Spearman |
| rdkit | 2026.3.6 | chemical structure handling in the docking and physicochemical layers |
| pillow | 12.3.0 | imaging dependency of rdkit (not imported directly by the code) |
| jsonschema | 4.26.0 | validation against JSON Schema Draft 2020-12 (pulls attrs 26.1.0, jsonschema-specifications 2025.9.1, referencing 0.37.0, rpds-py 0.30.0, typing_extensions 4.16.0) |
| pytest | 9.1.1 | tests (pulls pluggy 1.6.0, iniconfig 2.3.0, packaging 26.3, pygments 2.21.0; colorama 0.4.6 on Windows only; exceptiongroup 1.3.1 and tomli 2.4.1 on Python < 3.11) |

Exact pins are in [`requirements-track2.lock`](requirements-track2.lock). numpy, scipy and rpds-py publish prebuilt wheels for CPython 3.10+ on Windows, macOS (x86_64, arm64) and Linux, so no compiler is needed.

| Data/service | Version | Used for |
|---|---|---|
| ChEMBL | ChEMBL_37 (the manifests record the version; the release date is not recorded) | approved drugs, mechanisms |
| STRING | v12.0 (physical links and aliases) | interaction network |
| Reactome | release 97 | pathway members |
| DailyMed (FDA/NLM) | queried 2026-10-06 | label text |
| NCBI E-utilities | queried 2026-10-06 | PMID existence check |
| PrimeKG, HPO v2026-09-01, MONDO, OmniPath, DepMap 24Q2 | SHA-256 and size in `data/sources/record.json` | the four lines of evidence |
| openFDA adverse events (FAERS) | queried 2026-10-06 | descriptive toxicity summary |

## Run it (Linux, macOS and Windows)

### 1. Set up

Linux / macOS:

```bash
git clone https://github.com/jacandelerocervantes-cmd/MVA-CANDELEROJ98-CODE-TRACK2.git
cd MVA-CANDELEROJ98-CODE-TRACK2
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements-track2.lock
python -m pytest -v              # 561 passed, 25 skipped in a clean clone (586 collected; the 25 need the downloaded sources)
```

Windows (PowerShell):

```powershell
git clone https://github.com/jacandelerocervantes-cmd/MVA-CANDELEROJ98-CODE-TRACK2.git
cd MVA-CANDELEROJ98-CODE-TRACK2
py -3 -m venv .venv
# if activation is blocked: Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements-track2.lock
python -m pytest -v              # 561 passed, 25 skipped in a clean clone (586 collected; the 25 need the downloaded sources)
```

### 2. Reproduce the analysis

**One command.** `python -m track2.run_all --list` shows the 22 stages and their rough duration; `python -m track2.run_all` runs them in order (`--dry-run` prints the commands, `--from-stage N` resumes, `--skip N ...` omits slow stages). Stages 17 to 21 run the four lines and write the report, methods, official form and PDF. The last stage, `python -m track2.verify_outputs`, compares the deterministic outputs with `out/checksums.sha256`. Stages that read live services (DailyMed, openFDA, ChEMBL, PubMed) depend on the recorded cache (`out/data_manifest.json` lists the source hashes); they are listed separately in the checksum file and are never silently accepted. The individual commands behind each stage follow.

The first run downloads ChEMBL (about 13 minutes of API calls), STRING (about 30 MB) and Reactome pathway members into `data/` (ignored by git). Later runs use the cache and finish in about a minute per mechanism. Set `PYTHONPATH` to `src` (or install the package in editable mode).

Linux / macOS:

```bash
export PYTHONPATH=src
python -m track2.from_submission --submission <Track 1 submission CSV> --out out   # modules from the Track 1 CSV + proximity for each gene
for m in bub1b_sac_lof bub1b_down_proteotoxic bub1b_down_p53_senescence bub1b_down_growth_mtor_ampk; do   # hand-defined modules, comparison only
  python -m track2.run_pipeline --mechanism config/mechanisms/$m.json --out out/comparison_manual
done
python -m track2.run_pipeline --mechanism config/mechanisms/fanconi_anemia_demo.json --out out   # control disease of the validation
python -m track2.validation --out out        # frozen protocol P1-P4 (about 30-40 minutes)
python -m track2.build_candidates --out out  # about 3 minutes the first time (DailyMed labels, PubMed checks)
# (the addendum-10 report generator is superseded: run stage 21 of run_all for the report)
# Evidence layers (all optional, run after build_candidates; each caches its API responses)
python -m track2.evidence --candidates out/candidates.json --out out                 # about 10 minutes the first time
python -m track2.faers_disproportionality --candidates out/candidates.json --out out
python -m track2.direction --universe data/universe.json --candidates out/candidates.json --out out
python -m track2.physchem --candidates out/candidates.json --universe data/universe.json --out out   # measured affinity, properties, ADME (ChEMBL)
python -m track2.consensus --out out   # 9 configurations: robust / intermediate / fragile (addendum_robustness.md)
for g in BUB1B CENPE MAD1L1 IGBP1; do python -m track2.validate_submission_modules --mechanism config/mechanisms_from_submission/sub_$g.json --out out; done   # about 30-90 minutes each
python -m track2.literature_context --candidates out/candidates.json --out out   # curated, PMID-verified context
python -m track2.build_candidates --out out/comparison_manual --mechanism-dir config/mechanisms_extrapolation --dest out/arm_extrapolation   # separate extrapolation arm
python -m track2.pediatric_arm --out out   # paediatric-labelled search over the 300 best drugs (addendum_pediatric.md; labels reviewed in config/pediatric_label_review.json)
python -m track2.age_view --out out   # all age bands; add --age-years N when the age is known
python -m track2.depmap_download && python -m track2.prism_l4 --out out   # L4: DepMap/PRISM (about 540 MB of downloads)
python -m track2.get_vina && python -m track2.docking_structures && python -m track2.docking --stage validate   # docking with validation (pip install -r requirements-docking.txt first)
python -m track2.candidate_register --out out   # the complete justified register
python -m track2.reversal --universe data/universe.json --candidates out/candidates.json --out out
python -m track2.benchmark --out out                                                 # about 25 minutes
python -m track2.integrate --out out
python -m track2.review_page --out out && python -m track2.snapshot --out out
# (the addendum-10 report generator is superseded: run stage 21 of run_all for the report)
```

Windows (PowerShell; the layer commands are the same, one per line):

```powershell
$env:PYTHONPATH = "src"
python -m track2.from_submission --submission <Track 1 submission CSV> --out out
foreach ($m in "bub1b_sac_lof","bub1b_down_proteotoxic","bub1b_down_p53_senescence","bub1b_down_growth_mtor_ampk") {
  python -m track2.run_pipeline --mechanism "config/mechanisms/$m.json" --out out/comparison_manual
}
python -m track2.run_pipeline --mechanism config/mechanisms/fanconi_anemia_demo.json --out out
python -m track2.validation --out out
python -m track2.build_candidates --out out
# (the addendum-10 report generator is superseded: run stage 21 of run_all for the report)
```

Another case: point `--submission` to another Track 1-format CSV (genes in the first word of `notes`, weights in `epcr`) and run the same commands, with `--config-dir` set to a separate directory so that the modules of different cases are not mixed in `config/mechanisms_from_submission/` (only the first stages are reusable, see the report); or write a mechanism file in `config/mechanisms/` by hand.

### 3. Optional: genome covariates (requires the patient VCF, not distributed)

```bash
python -m track2.genome_covariates --vcf <sample.vcf.gz> --gtf <MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz> \
  --cds-cache <cds_sequences.json> --gnomad-cache <gnomad_v4_frequencies.json>
```

The files for `--gtf`, `--cds-cache` and `--gnomad-cache` are produced by the Track 1 repository. The per-variant detail goes to `data/private/` (never versioned); only `out/genome_covariates_summary.json` is published.

## Checksums

**Outputs that read live label text.** The label tables, the caution list and the quote checks read DailyMed and openFDA text, and labels change over time, so those outputs can differ between runs. They are listed as network-dependent in the checksums, and the label readings of this study were fixed on 2026-10-08.

The record is `out/checksums.sha256` (184 files at the time of writing: 126 deterministic and 58 that embed live service text, depend on files outside the repository or on the machine, and need the same cache). `python -m track2.verify_outputs` compares the files with it and exits non-zero on any difference in the deterministic group; `python -m track2.verify_outputs --write` records the current state. The hashes of the protocol files at the moment the amendments were approved are in `protocol/hashes_review2.sha256`. Source files downloaded for the four lines are hashed in `data/sources/record.json`.

## Gene landscape, differential diagnosis, variant mechanism and case profile

`src/track2/disease_landscape.py` derives, for any rare disease and by code, the genes associated with it from the HPO annotations and Orphanet, assigns each a mechanism class by a fixed rule and keeps one code-verified Europe PMC sentence per gene (`out/mva_landscape.json`, `out/landscape_fanconi_anaemia.json`; curated genes and genes 'proposed, not curated'). `src/track2/differential.py` classifies the genes of a Track 1 output with rules written before they were applied (`protocol/addendum_differential.md`; `out/differential_diagnosis.json`); `src/track2/variant_level_consequence.py` gives its categorical consequence classes. `src/track2/variant_mechanism.py` measures, from public sources, the mechanism of the two variants at protein level (`out/variant_mechanism.json`), and `src/track2/phase_inference.py` states the phase question as an explicit inference with a sensitivity table (`out/phase_inference.json`). `src/track2/case_profile.py`, `case_reference.py` and `case_contrast.py` define a minimal, privacy-safe case profile with an explicit 'unknown' state for every field and contrast it with published cases (`out/case_profile_child.json`, `out/case_reference_bub1b.json`, `out/case_contrast.json`). `config/track1_gene_level_extract.json` (made by `python -m track2.track1_extract`, which needs the Track 1 repository once) holds the gene-level, categorical subset of the Track 1 output that `differential.py` and `variant_level_consequence.py` read, with no position, allele, genotype call or read metric; with it both modules rerun without the Track 1 repository and `tests/test_track1_extract.py` checks that they reproduce their outputs. The files that fetch literature or annotation (landscape, variant mechanism, case reference, phase) need the network and are listed as network-dependent in the checksums. The text of sections 3.8 to 3.11 and Supplement S13 of the report is generated from them by `src/track2/landscape_section.py`.

## Empirical checks of the method

`out/empirical/` holds measurements of the method itself (not of any drug), each produced by a module `src/track2/empirical_*.py` with its own tests: `label_robustness.json` (stability of the selection when evidence labels are corrupted), `literature_control.json` (a mention-based proxy of rule R2 against known indications of 26 rare diseases; needs Europe PMC and a local cache in `data/cache_empirical/`), `controls.json` (positive control of the caution layer and catch rate of the quote verifier), `variant_calibration.json` (calibrated strength of the computational evidence for the second variant; the cached texts of the two papers are not distributed, rebuild them by running the module), `age_support.json` (fraction of the paediatric age range in which a label supports a drug; label support under an unknown age, not a recommendation), `power_analysis.json` (planning numbers for the refuting experiments, not results) and `scalability.json` (timings and determinism on two synthetic cases; run in separate folders `out_case_*` and `config_case_*`, which are not distributed). The text of section 3.7 and Supplement S12 of the report is generated from these files by `src/track2/empirical_section.py`.

## Field names in the JSON outputs

The label tables (`out/caution_list.json`, `out/age_label_table.json`) keep some field names from the first version of the code, in Spanish; the values and all text are in English. Glossary:

| Field | Meaning |
|---|---|
| `aviso`, `protocolo`, `regla`, `universo`, `nota`, `motivo` | notice, protocol, rule, universe of drugs, note, reason |
| `clases`, `farmacos`, `etiqueta`, `criterio` | classes, drugs, label (or class name), membership criterion |
| `citas`, `texto_literal`, `seccion`, `fuente`, `fecha_lectura`, `truncada` | quotes, verbatim text, label section, source (DailyMed or openFDA), date read, truncated |
| `estado_por_grupo`, `estados`, `labelled_for_group` | status by age group, explanation of the statuses, groups for which the label states use |
| `frase_literal_sustento`, `advertencias_pediatricas` | verbatim sentence that supports the status, paediatric warnings quoted from the label |
| `clasificador_automatico`, `revision_manual`, `revision_manual_en_config` | automatic classifier, recorded manual reading, whether a reading is recorded in `config/pediatric_label_review.json` |
| `sin_cita`, `etiquetas_ilegibles`, `clase_vacia` | drugs without a verifiable quote, unreadable labels, empty class |
| categories `toxicidad_medular`, `embrio_fetal`, `segundos_tumores`, `active_malignancy` | bone-marrow toxicity, embryo-fetal toxicity, second tumours, active malignancy |

## Human Phenotype Ontology

This service/product uses the Human Phenotype Ontology (release 2026-09-01). Find out more at http://www.human-phenotype-ontology.org. The HPO files are used unaltered, under the HPO licence (https://hpo.jax.org/license).

Where each third-party resource was obtained, how access is obtained and how it was used: [`DATA_ACCESS.md`](DATA_ACCESS.md).

## Licence

As the terms of the dataset (Hackathon Rules and Data Transfer Agreement) state, participant submissions and results (predictions, code and reports) are released under CC BY 4.0 and may be reused with attribution; this applies to the submission (reports, tables and the other submitted files). The code is released under the same terms (CC BY 4.0; the full text is in `LICENSE`). Third-party resources keep their own licences (see `THIRD_PARTY_LICENSES.md`); where a resource is CC BY-SA, derived tables that reproduce its content keep that licence.

## Data privacy

Patient-derived detail is not versioned (`data/`, including `data/private/`, is ignored); the Track 1 submission CSV (`candeleroj98_gatk-manta-fusion.csv`, which carries the variant coordinates required by the challenge format) is **not** included in this repository: take it from the Track 1 submission and pass its path to `track2.from_submission --submission`, or set the environment variable `TRACK1_SUBMISSION` before `track2.run_all`. The report names variants at protein level only. **All data must be deleted 30 days after the hackathon ends** (deadline: 23 November 2026, 30 days after the submission deadline of 24 October 2026; confirm by e-mail to cnb@sagebase.org): `python -m track2.delete_data` lists what would be deleted (dry run) and `--execute` deletes it and writes an attestation. See `docs/data_protection.md`. The report uses HPO terms and general descriptions only; birth details from the confidential clinical document are not reproduced.

## Acknowledgments

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

No information that could re-identify the patient or family is included in this repository. We thank the Wilhelm Foundation for their guidance and insight on the project.

## Generative-AI disclosure

See [`config/ai_declaration.md`](config/ai_declaration.md) and section 11 of the report.

## Evidence engine tested on other rare diseases (Section 3.12 and Supplement S14 of the report)

`src/track2/engine_score.py` and its successors (`engine_score_v2.py`, `engine_v3.py`, `engine_v4*.py`, `engine_v5*.py`, `engine_case_study.py`) are deterministic code (no language model in the ranking): Europe PMC retrieval, mention counts of a drug with the disease genes, a PrimeKG-based prior and logistic weights learned by cross-validation grouped by disease. Each stage has a written protocol dated before its scoring (`protocol/addendum_engine*.md`), a frozen model with SHA-256, a lock file for the sealed set and tests. Results are in `out/engine/validation*.json`, the benchmarks in `out/engine/benchmark*.json`, the worked case (Wilms tumour) in `out/engine/case_wilms.md`, the application of the frozen engine to MVA and Fanconi anaemia in `out/engine/apply_mva.md` and `apply_fanconi.md` (`src/track2/engine_apply.py`, exploratory, not a candidate list), the coverage audit in `out/engine/audit_llm_coverage_dev.md`. Language-model labels (two blind annotators, guidelines in `protocol/annotation_guideline.md`, `relevance_guideline.md`, `relevance_guideline_v2.md`) are in `out/engine/llm_labels/` and `out/engine/llm_labels_v5_ul/`. The report text of Section 3.12 and Supplement S14 is generated by `src/track2/engine_section.py`. Three designs of the per-pair retrieval test are documented (full, reduced, ultralight); only the ultralight one was executed. Not in the repository (git-ignored): caches, the two embedding models (about 840 MB, revisions in `data/models/revisions.json`) and the bulky annotation shards. The v1 result is the only pre-registered one; everything after it is post hoc development and says so.

## Candidates by line, cross-check and experiment plan (Sections 3.13 to 3.15 and Supplements S14 and S15 of the report)

Each computational line proposes its own candidates and none is treated as the truth: literature layer, engine, A (knowledge graph), B1 and B2 (network), C (PRISM drug-level vulnerability, `src/track2/line_c_adjusted.py` adjusts it for lineage and proliferation), D (signed direction), E and E-P (signature reversal on LINCS, `src/track2/line_e.py`; the patient signature comes from GSE22206, 2 case and 4 control arrays). `src/track2/candidates_by_line.py` collects them (`out/candidates_by_line_mva_v2.json/.md`) with age-band label readings and safety flags; `src/track2/multi_evidence.py` is the stochastic multicriteria cross-check (`out/multi_evidence_mva*.json`), reusable on another case with `--out`. The protocols are `protocol/addendum_line_e_signature_reversal.md`, `addendum_line_c_adjusted.md`, `addendum_multi_evidence.md` (post hoc, with the three runs and their reasons), `addendum_refutability.md` and the experiment design `protocol/follow_up_experiments_v2.md` (`out/experiment_plan.json`). The text of the new sections is generated by `src/track2/lines_section.py` and `followup_section.py`. Downloaded public data (the two LINCS libraries and one GEO series with its platform file) are recorded with URL, date, size and SHA-256 in `data/lincs/record.json` and `data/geo/record.json` and are git-ignored. Nine computational follow-ups (CF1 to CF9) are defined in the plan and were not run.
