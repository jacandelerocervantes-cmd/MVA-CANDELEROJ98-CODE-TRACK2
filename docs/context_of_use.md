# Context of use and credibility statement

Written in the spirit of risk-informed credibility assessment for computational models (the ASME V&V 40 style of reasoning: question of interest, context of use, model influence, decision consequence, risk, then credibility evidence proportionate to the risk). It is a self-assessment by the authors, not an external review, and it is not a regulatory submission.

## Question of interest

Which already approved drugs are worth a human expert's closer look as hypotheses for a child with mosaic variegated aneuploidy (MVA1, *BUB1B*), given that their targets lie close, in a protein-interaction network, to pathways disrupted by the primary mechanism?

## Context of use

- **Is:** a ranked, explicitly uncertain list of hypotheses ([`out/candidates.json`](../out/candidates.json), and `out/review.html` when generated with `python -m track2.review_page`), with risk flags, label information by age band, toxicity descriptions and requirement states.
- **Is not:** a treatment recommendation, a dosing tool, a safety clearance, or evidence of efficacy. Intended readers are domain experts (clinical geneticist, pharmacist, researcher) deciding what to study further.

## Decision consequence, model influence, risk

| Element | Assessment |
|---|---|
| Decision consequence | Potentially high if misused (a drug given to a child), but the output is not designed to drive a clinical decision; the intended decision is only which hypotheses to investigate further. |
| Model influence | Low: it is one input among literature, clinical judgement and expert review; every candidate carries `direction_review = REQUIRES_REVIEW` and unassessed requirements. |
| Risk level | Moderate if used as intended; unacceptable if read as a recommendation. Hence the banner and the NOT_ASSESSABLE states. |

## Credibility evidence available in this repository

| Evidence | Where | What it shows |
|---|---|---|
| Frozen, pre-registered validation protocol | [`protocol/frozen_protocol.md`](../protocol/frozen_protocol.md), [`protocol/deviations.md`](../protocol/deviations.md) | Rules were fixed before results; deviations are listed |
| Held-out checkpoint genes ranked well | [`out/validation_results.json`](../out/validation_results.json) | AUC 0.92 (95 % CI 0.89-0.95): pathway-level signal |
| Positive control in a second disease | same file | Fanconi anemia control **not** recovered: a negative result, reported |
| Sensitivity to network threshold | same file, [`protocol/deviations.md`](../protocol/deviations.md) | Ranking **unstable** with the interactome threshold (Jaccard 0.48 and 0.38) |
| Degree bias | same | Spearman 0.28 between z and target degree: reduced, not removed |
| Method from the literature | [`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md) | Guney 2016, Cheng 2018, reimplemented |
| Unit tests, locked environment | [`tests/`](../tests), [`requirements-track2.lock`](../requirements-track2.lock), [`coverage.md`](coverage.md) | Core numerics tested; most pipeline logic only end to end (32 % line coverage) |
| Data provenance and checksums | [`README.md`](../README.md), `out/manifest_*.json`, `out/data_manifest.json` (from `python -m track2.snapshot`) | Versions and hashes of inputs and outputs |
| Output schema and citation checks | [`schemas/candidate_drug.schema.json`](../schemas/candidate_drug.schema.json) | Records valid; PubMed IDs confirmed to exist (existence only, not that they support the claim) |
| Generative-AI disclosure | [`config/ai_declaration.md`](../config/ai_declaration.md) | How AI assistance was used |

## Gaps (what the model cannot support)

- **No direction of effect.** Proximity does not say whether a drug helps or worsens chromosomal instability; every candidate needs manual review.
- **Drug-level ranking is exploratory.** Unstable under the network threshold, and the one external positive control failed.
- **No efficacy or safety evidence for this patient.** Age, tumour status, renal, hepatic and haematological function, allergies and current medication are unavailable, so these requirements are NOT_ASSESSABLE; only the documented nephrocalcinosis caution can be evaluated, from label text.
- **Heuristic label parsing.** Age bands, warnings and excipient alerts come from keyword scans of DailyMed text: verify by hand.
- **FAERS summaries are descriptive:** spontaneous reports with no denominator or incidence.
- **No dosing, no pharmacokinetics, no pharmacogenomic star-allele calling, no drug-drug interaction analysis.** The affinity-prediction layer is not implemented.
- **No independent expert review** of the shortlist, no calibration of the combined score, one disease validated with one control.
- Data sources drift with releases; reproducibility is exact only with the cache (see the checksums in the README).

See [`roadmap.md`](roadmap.md) for planned work. Any use beyond hypothesis generation would require, at least, independent expert review and evaluation of the patient requirements.
