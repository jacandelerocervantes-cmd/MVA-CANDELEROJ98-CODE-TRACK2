# Annotation guideline for the direction of a sentence about a drug

Task: for each item you receive a drug name and one sentence from a scientific abstract or full text. Decide what the sentence, read in its own wording, says about the drug in the context of the disease or biological problem the sentence describes.

Labels (choose exactly one):
- benefit: the sentence reports that the drug improved, rescued, reduced or prevented a disease-relevant feature, or killed or sensitised abnormal (for example aneuploid or tumour) cells in a way presented as favourable to the study question.
- harm: the sentence reports toxicity, an adverse reaction, worsening, resistance, genotoxicity, a safety concern or a contraindication linked to the drug.
- no_effect: the sentence reports that the drug did not change the feature, did not rescue, failed, or that no toxicity or no difference was observed.
- unclear: the sentence names the drug but does not, by itself, state a benefit, harm or absence of effect (background, mechanism description, dosing, review statement).

Rules: use only the sentence and the drug name; do not use outside knowledge about the drug; a negated benefit is no_effect, a negated harm is no_effect; if both benefit and harm appear, choose the one the sentence presents as its main finding, otherwise unclear. Also give a confidence from 1 (guess) to 3 (certain).
