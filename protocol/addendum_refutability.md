# Addendum: which hypotheses can be refuted by computation and which need an experiment

Written before the table of Section 3.13 was built. Label of every output: **exploratory; hypotheses for follow-up; never evidence of efficacy; not a treatment or dose recommendation.** No new computation: the table only reads existing outputs.

## 1. Principle
A hypothesis that can be refuted computationally is tested computationally. A hypothesis that cannot gets an experimental design (`protocol/follow_up_experiments.md`).

## 2. Classification rule (fixed before the table)
1. **Computationally refutable.** The hypothesis is a statement about a gene, a pathway or a cell-line effect that public data in the project can measure, that is, that one of lines A to E can score (A: knowledge-graph query; B: phenotype and mitosis modules; C: cell-line vulnerability; D: signed effect on the checkpoint, exploratory; E: signature reversal).
2. **Needs experiment.** The hypothesis requires cells with partial BUB1B loss, a model of the case, or a clinical outcome, none of which the public data contain. This holds for the core claim of every drug in the report (benefit in a hypomorph, selectivity of a tumour-context effect, harm to non-tumour aneuploid cells).
3. **Hazard drugs.** They are tested computationally for harm (concordance with the disease signature in line E, direction in line D, membership of `out/caution_list.json`). The experiment measures toxicity, not benefit.

## 3. Status vocabulary (one per drug)
* **computationally unsupported**: the lines that could score the drug did not support it, and no line shows an effect in the opposite direction.
* **not computationally testable**: line A did not score the drug and line D could not call its direction, so the available lines cannot say anything about it.
* **hazard signal**: the drug is in a caution-list class, or is concordant with a disease signature in line E, or line D says it weakens the checkpoint. A hazard signal comes from a different direction than benefit and is not a refutation of safety either way.

"Not supported by computation" is not "refuted". The word "refuted" is used only if a line shows an effect in the opposite direction; line E concordance (the drug mimics the disease signature) is reported as a hazard signal with the limits of the signature used.

## 4. Sources (existing files only)
`out/four_lines_review2.json` (usable lines and support sets), `out/proposals_review2.json` (group, level, line A rank, line B2 and C flags, line D note, refuting experiment, hazard flag), `out/signed_effect.json` (line D), `out/vulnerability_e6.json` (line C), `out/kg_*` and `out/phenotype_*` (lines A and B, through the files above), `out/line_e_signature_reversal.json` (line E), `out/caution_list.json`, and the literature levels as recorded in `out/proposals_review2.json`.

## 5. What the table cannot say
A row with "unsupported" or "no hazard signal" says only that the public-data lines did not show it. It does not say that the drug is safe or useless. Whether the candidate list holds depends on the experiment, which is not run.
