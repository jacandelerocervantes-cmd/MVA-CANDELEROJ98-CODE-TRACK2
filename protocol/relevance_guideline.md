# Relevance guideline: is the sentence about the drug as a treatment of this disease?

Task: for each item you receive a drug name, a disease name and one sentence taken from the title or abstract of a scientific paper. Decide whether the sentence **evaluates the drug as a treatment** of the disease named in the item, or of the biological mechanism or class of disease that the sentence itself ties to that disease. You judge relevance only; you do not judge whether the drug works (that is a different task with a different guideline).

Labels (choose exactly one):
- yes: the sentence evaluates, tests, reports on or recommends the drug as a treatment, therapy, management, prevention or rescue strategy for the named disease, or for the mechanism or class of disease that the sentence states is the one that disease belongs to (for example the sentence says "ciliopathies" or "mitochondrial disorders" and the disease name is a member of that group in the sentence's own wording). The outcome may be positive, negative or null: a trial that found no effect is still `yes`. Safety or toxicity of the drug in patients or models of that disease is `yes`.
- no: the drug is mentioned for a different purpose. Examples: a different disease or condition; a laboratory reagent, inhibitor or tool compound used to study biology; an assay, a staining or imaging use, a biomarker, a diagnostic agent; a comparator or a background therapy listed without being evaluated; a contraindication in an unrelated context; a statement about the drug's pharmacology with no treatment aim; a list of drugs with no evaluation.
- unclear: the sentence names the drug but, read alone, you cannot tell whether it is being evaluated as a treatment of this disease or of its mechanism class (for example a bare list of therapies, a fragment, or a pronoun that hides the condition).

Rules:
1. Use only the item: the sentence, the drug name and the disease name. Do **not** use outside knowledge about what the drug is used for, what the disease is, or what is approved. A drug that you know treats the disease is still `no` or `unclear` if the sentence does not say so.
2. The disease name need not appear in the sentence; the paper was retrieved for the disease, but the sentence must itself support the link to the disease or to its mechanism class. If the sentence only describes a treatment effect with no condition named at all, choose `unclear`.
3. A sentence that evaluates the drug in a model of the disease (cells, animals) counts as `yes` only if the sentence presents the model as a model of this disease or of its mechanism class; a generic cell-biology experiment is `no`.
4. A sentence about a different drug or a drug class (not this drug) is `no`; a sentence that names the drug only inside a list that is evaluated as a group is `yes` if the group is evaluated as treatment of the disease, otherwise `unclear`.
5. If the sentence is negated ("did not improve"), the label is still based on what is being evaluated, not on the result.
6. Give a confidence from 1 (guess) to 3 (certain).

Output per item: id, relevance (yes / no / unclear), confidence (1-3). The direction of the effect (benefit / harm / no_effect / unclear) is labelled separately with `protocol/annotation_guideline.md`.
