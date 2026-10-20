# Relevance guideline v2: is the sentence about the drug as a treatment of this disease?

Status: post hoc revision of `protocol/relevance_guideline.md` (v1). It keeps the task and the labels and makes the `yes` / `unclear` / `no` threshold operational. Context and honesty clauses: `protocol/addendum_engine_v4_llm_guideline_v2.md`. All examples below are invented for this guideline; none is taken from any real item. "Drug X", "Drug Z" and "disease Y" are placeholders.

## 1. Task and labels (unchanged from v1)

For each item you receive a drug name, a disease name and one sentence from the title or abstract of a paper. Decide whether the sentence **evaluates the drug as a treatment** of the named disease, or of the mechanism or class of disease that the sentence itself ties to that disease. You judge relevance only, not whether the drug works.

- `yes`: the sentence evaluates, tests, reports on, uses or recommends the drug as treatment, therapy, management, prevention or rescue for the disease or its stated mechanism class. Positive, negative and null results are all `yes`. Safety or toxicity of the drug in patients or models of that disease is `yes`.
- `no`: the drug appears for another purpose (see rule N).
- `unclear`: the sentence names the drug but, read alone, you cannot tell whether it is evaluated as a treatment of this disease.

Output per item: id, relevance, confidence 1-3. Use only the item (sentence, drug name, disease name). No outside knowledge about the drug, the disease or approvals.

## 2. Three conditions

A sentence earns `yes` only if all three hold, each **read from the sentence itself**:

- **C1 Drug.** The sentence names the drug (or a clear synonym or abbreviation of it, or a group that explicitly contains it, see rule L).
- **C2 Treatment aim.** The drug is given, tested, prescribed, used, recommended or compared as a therapy, preventive or rescue measure, and the sentence reports or states that use (an outcome, a regimen, a dose, a recommendation, a patient treated, a status of the therapy).
- **C3 Condition link.** The sentence ties that use to the named disease, to a synonym, subtype or abbreviation of it, or to a mechanism class that the sentence states and that includes it (v1 wording). A gene, a phenotype or a sign that is stated in the sentence as the defining feature of the disease counts as a link only if the sentence itself says it is that disease or its subtype.

## 3. Decision rules

### Rule C: the condition link decides between `yes` and `unclear`

If C1 and C2 hold but C3 does not (no disease name, synonym, subtype or class in the sentence; only "patients", "children", "the cohort", "this condition", "these subjects", a pronoun, or a symptom or organ that is not stated to be the disease), the label is `unclear`, not `yes`. The fact that the paper was retrieved for the disease, or that the item shows a disease name, does not supply C3. The annotator must not fill the gap by guessing which disease the paper is about.

- Synthetic: "Children given drug X for 12 weeks had fewer admissions than before treatment." Disease named in the item: disease Y. The sentence names no condition. Label: `unclear`.
- Synthetic: "Children with disease Y given drug X for 12 weeks had fewer admissions." C1, C2, C3 hold. Label: `yes`.
- Synthetic: "Drug X reduced seizure frequency in this cohort." Disease named in the item: disease Y (not an epilepsy). Seizures are a symptom, not stated as disease Y. Label: `unclear`.

### Rule T: a treatment aim must be stated, not implied by the drug's usual role

If C3 holds but the sentence gives no treatment use (only that a drug exists, or is "available", "in development", or a pure pharmacological property), the label is `unclear` when a treatment use is plausible but not stated, and `no` when the sentence is clearly about something else (rule N).

- Synthetic: "Drug X is a kinase inhibitor with a short half-life." Label: `no` (pharmacology only).
- Synthetic: "Several agents, including drug X, are under investigation for disease Y." Label: `yes` (investigation as treatment of disease Y is stated; the outcome is open).
- Synthetic: "Drug X and disease Y were both discussed at the meeting." Label: `unclear`.

### Rule B: background and review statements about therapy are `yes` when C1-C3 are all in the sentence

A background statement counts as evaluation or recommendation of the drug when it states a treatment role: "is used", "is approved", "is the standard", "is first-line", "improved outcomes" in a population defined as having the disease. A background statement of mechanism or epidemiology without a treatment role is `no` or `unclear` by rule T.

- Synthetic: "Drug X is the first-line treatment of disease Y." Label: `yes`.
- Synthetic: "Disease Y affects about 1 in 50,000 people and is caused by variants in gene G." The drug is not treated at all: `no` (and the sentence would not normally reach annotators without the drug).
- Synthetic: "Drug X was introduced in the 1980s and, together with transplantation, improved the outcome of disease Y." Label: `yes`.

### Rule L: lists

- The drug is in a list and the list is the object of a treatment statement that includes C3: `yes`. Synthetic: "First-line options in disease Y include drug X, drug Z and drug W."
- The drug is in a list that only enumerates compounds with no treatment statement or no C3: `unclear` if the sentence is about therapy of an unnamed or unspecified condition ("approved therapies include drug X and drug Z"), `no` if the list is a catalogue of reagents, targets or side effects.
- The list says "therapy with drug X and drug Z" for disease Y: `yes` for both drugs (the group is evaluated as treatment; do not try to separate their contributions).

### Rule P: patients, trials, regimens

- A patient or case report of disease Y treated with drug X (given as part of a regimen, with or without a stated result): `yes`. Synthetic: "A patient with disease Y was treated with drug X and drug Z, and remitted."
- A treatment named only as part of a regimen for a different, co-existing condition, or as a symptom-control drug in the same patient: `no` when the sentence says the other condition is the target, `unclear` when it does not say.
- Trial design sentences (randomised, placebo-controlled, phase, open-label) that name the drug and the disease as the target: `yes`. If the trial sentence names the drug but the population only as "adults" or "the study group": `unclear` by rule C.
- A population named as having disease Y in which a drug is used for an unrelated indication (for example an anticoagulant for a clot): `no` if the unrelated indication is stated, `unclear` if not stated.

### Rule M: models

- Cell, animal or organoid work: `yes` if the model is presented as a model of disease Y (or its class) and the drug is tested for a therapeutic effect there. Synthetic: "In a mouse model of disease Y, drug X rescued the phenotype."
- A generic experiment (signalling pathway, cell line, inhibitor used to probe biology) with no tie to disease Y: `no`. Synthetic: "Cells were incubated with drug X to block pathway P."
- A mechanism-only sentence (what the drug does to a pathway) that names disease Y but states no treatment effect or aim: `unclear` if it is framed as potential benefit ("may be a candidate", "target for treatment"), `no` if it is framed as basic biology. Synthetic: "Drug X is a candidate for disease Y because it blocks pathway P" is `yes` only when the sentence reports testing or proposing the drug as treatment; as a bare proposal without test it is still `yes` (a proposal is a recommendation), as an unframed pathway statement it is `no`.

### Rule N: `no` conditions

Label `no` when the sentence clearly shows one of these, with confidence 2-3:
- a different disease or indication is the treatment target (synthetic: "Drug X is approved for disease V." with item disease Y);
- the drug is a reagent, tool compound, stain, imaging or diagnostic agent, assay component or biomarker;
- the drug is a comparator or background therapy that is not itself evaluated ("all patients also received drug X" with a different drug evaluated);
- the drug is only discussed for pharmacology, adverse effects unrelated to treating disease Y, or contraindication in an unrelated context;
- a different drug or class is the subject.

### Rule U: `unclear` conditions (the residual class)

Label `unclear` when the case is neither `yes` by C1-C3 nor `no` by rule N. This includes: sentence fragments; a pronoun or "the disease" hiding the condition; a drug in a list with no treatment statement; a bare mention of a therapy outcome with no condition; a title question ("drug X in disease Y?"), where evaluation is the topic, which is `yes` only when it names drug and disease and clearly asks about treatment. Do not use `unclear` as a hedge for sentences that match rule N, and do not use `yes` to mean "the drug probably treats the disease".

### Rule K: confidence

3 = the sentence settles it under these rules; 2 = one condition rests on a judgement; 1 = a guess between two labels. Prefer the lower of two plausible labels in the order `no` < `unclear` < `yes` only when a condition C1-C3 is truly absent from the sentence; if all three are present, say `yes` even if you doubt the paper's quality or the drug's efficacy.

## 4. Decision order (text flowchart)

1. Does the sentence name the drug (C1)? No: `no` (the item should not occur; use `no`).
2. Is the drug clearly a reagent, tool, stain, assay, biomarker, diagnostic, comparator/background, or a drug for a different stated indication (rule N)? Yes: `no`.
3. Does the sentence state a use of the drug as therapy, prevention or rescue, a test of it, a recommendation, or a patient or model treated with it (C2)? No: go to step 6.
4. Does the sentence name disease Y, a synonym or subtype, or a mechanism class that includes it, or a model presented as a model of it (C3)? Yes: `yes`.
5. C3 not satisfied: `unclear`.
6. C2 not satisfied: if the sentence is about an unrelated purpose, `no`; if it is a pure pharmacology/mechanism sentence with no treatment aim, `no`; if a therapy role is plausible but not stated (fragment, bare list, pronoun), `unclear`.
7. Assign a confidence.

## 5. What the sentence alone does not decide

The following are never decided by the sentence and must not be inferred:

- Whether the drug works, is approved, or is standard care for the disease (relevance is not efficacy and not indication knowledge).
- Which disease a paper is "really" about when the sentence does not say; the item's disease name is a retrieval label, not evidence.
- The direction of effect (separate guideline `protocol/annotation_guideline.md`).
- Whether a symptom or phenotype belongs to the disease named in the item, unless the sentence says so.
- The role of a drug in a regimen of several drugs (separate contribution), the dose or timing.
- Anything that needs the abstract or full text: if the context is in a neighbouring sentence, the answer is `unclear`, not `yes`.

## 6. Differences from v1 (summary)

v1 rule 2 said a sentence with no condition at all is `unclear`, but left room to link a sentence to the item's disease by plausibility; v2 makes C3 mandatory, separates `unclear` from `no` by rule N, and states that background treatment statements, regimens, case reports and trial sentences are `yes` when C1-C3 are present. Labels, output format and the no-outside-knowledge rule are unchanged.
