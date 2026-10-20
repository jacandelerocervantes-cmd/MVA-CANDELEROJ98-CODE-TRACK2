# Addendum 3: paediatric-labelled arm

Written on 2026-10-06 **before** the analysis it describes is run. It adds one analysis; nothing in the frozen protocol or in
the earlier addenda changes.

## Why

The age of the patient is not documented, but the patient is a child with a paediatric tumour predisposition. None of the 20 robust and
intermediate candidates of the main ranking has a labelled paediatric use (heuristic extraction of the US label). The shortlist was built
from proximity alone, so it never asked "which drugs near the modules *are* used in children?". This arm asks that question over the whole
scorable universe.

## Rule (fixed before running)

1. **Pool.** The 300 drugs with the highest combined score under the main configuration (STRING 700, `epcr` weights) among the scorable
   approved drugs (the universe has 1,040 scorable drugs).
2. **Label.** For each drug in the pool, the US label (DailyMed) is read with the same heuristic as the rest of the project
   (`labels.py`): paediatric status `ESTABLISHED` (the label states established use in some paediatric population) or `MIXED` (verify by
   hand). Drugs whose label says paediatric use is not established, or that have no label, are excluded from this arm.
3. **Reported per drug.** Youngest age stated in the label when parseable, labelled age bands, the label's paediatric excerpt, the rank
   of the drug (of the scorable universe) in each of the nine configurations of `addendum_robustness.md`, and the mechanism-class
   literature note when one applies.
4. **Stable.** A drug is called *stable* if its rank is within the best 100 in at least 7 of the 9 configurations; otherwise it is
   *unstable* (reported, not hidden).
5. **What this is not.** It is not a recommendation or a dose. Proximity of these drugs is weaker than for the main candidates
   (they were chosen by label, not by distance); the report shows their z so a reader can see how weak. The paediatric indication on a label is
   usually unrelated to this condition, and the heuristic extraction must be verified by hand for any drug taken further.

## What would count as failure

- No stable paediatric-labelled drug: reported as such. The conclusion would then be that the proximity signal for this case points only to
  agents without paediatric labels.

## Amendment 1 (2026-10-06, made AFTER seeing the first result; disclosed)

The heuristic status `MIXED` turned out to include labels that state paediatric use is **not** established (docetaxel, cabazitaxel,
natalizumab) next to labels that establish it for one indication (brentuximab vedotin, crizotinib), and it also hid two robust candidates that do
have a labelled paediatric use. Two changes follow, both decided after seeing the result and not before:

1. The labels of the drugs that matter (robust, intermediate and paediatric-arm drugs with a strong proximity) were read and the verdicts recorded in
   `config/pediatric_label_review.json`. The review overrides the heuristic for those drugs. It is a reading by the author's assistant, not a pharmacist
   or regulatory review.
2. A **headline** paediatric candidate is a drug with an *established* paediatric use after review **and** a proximity of z <= -2 (the project's
   primary threshold) in at least one module. Drugs with a paediatric label and weaker proximity are listed but are not candidates: they were selected by
   label, and their distance to the modules is close to chance.

The statement in the first report draft that "none of the robust candidates has a labelled paediatric use" was an artefact of the heuristic and is
corrected here and in the report.
