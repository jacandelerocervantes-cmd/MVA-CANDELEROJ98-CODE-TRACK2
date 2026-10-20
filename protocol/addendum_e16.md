# Amendment E16: literature-anchored proposals (declared post hoc, 2026-10-08)

This amendment is kept in its own file and not mixed into `addendum_review2.md`, which was approved on 2026-10-07.

## Why

The official description of Track 2 asks for approved drugs that are plausible candidates for further investigation **based on the disrupted mechanisms**. The pre-registered rule of addendum 11 (at least two usable computational lines) returned no candidate. Reporting only that result would answer a stricter question than the one the track asks. After the owner approved a realignment on 2026-10-08, the report adds a mechanism-to-hypothesis-to-experiment layer built from the primary literature. The computational lines are kept as a **filter** on it, and their negative result is kept in full.

## What is post hoc

Everything below was specified after the results of the four lines, the class analysis and the caution list were known. It is declared post hoc and is never presented as a finding of the lines.

## Rules (also in the docstring of `src/track2/literature_summary.py`)

- **R1** The drug is approved (in `data/universe.json`, ChEMBL `max_phase = 4`).
- **R2** *Constitutional hypothesis*: at least one verified finding with direction `may_help` and evidence level 1 or more on a feature of the child (growth and muscle, consequence 4) or on a consequence of aneuploidy (proteotoxic or energy stress, consequence 2; senescence, consequence 3), and no verified level-3 harm. Benefit counts from level 1 and harm only from level 3, so level-1 evidence for and against can coexist. *Tumour-context hypothesis*: a verified level-1 finding in an aneuploid or tumour system that is `may_help` or `selectively_kills_aneuploid_cells`; a drug that kills aneuploid cells selectively is also recorded as a hazard for the constitutional mosaic. A verified `may_harm` at level 1 or more marks the evidence as mixed.
- **R3** A class sibling is listed with its class representative, never alone, and keeps its own lower level.
- **R4** A verified harm at level 3 (a case in MVA or BUB1B) or a finding flagged as caution goes to the caution list, whatever the level. For a drug that is also proposed under R2, a verified `may_harm` at level 0 or 1 is not a reason to drop it: it is shown beside the proposal as evidence against and marks the evidence as mixed. (An earlier version of this text sent every `may_harm` to the caution list; the code never did, and the text is corrected here.)

Evidence levels: L0 inferred or indirect; L1 tested in a related system (aneuploid cells, chromosomal instability, yeast, an animal with a checkpoint defect); L2 tested in BUB1B/BUBR1-deficient or MVA-derived human cells, or in a hypomorphic model with rescue; L3 clinical observation in MVA or BUB1B.

## Disclosures about how the proposals were obtained

- **The search was not independent of the final drugs.** Of the 419 queries of the four searches, 143 name a drug that later became a proposal (rapamycin 51, metformin 32, mTOR 30). The searches looked for approved drugs per consequence, but their wording already steered toward these drugs.
- **The rules were written after the findings had been read.** The rules R1 to R4, the class siblings and the link texts were fixed after the literature files existed. File dates cannot show the order of editing and nothing here claims they do.
- **The weight of the evidence is low.** For sirolimus the entries on growth and muscle are, by sentence, 3 `may_help`, 10 `may_harm`, 2 `no_effect` and 1 `unclear` (by paper 3, 5, 2 and 1) in the first search, all level 0 and all verified, and one more `may_harm` (PMID 21380626) added in the third round; its single level-1 help comes from a Drosophila model on the aneuploid-cell consequence. For metformin the only level-1 help is oxidative damage in trisomy-21 fibroblasts, against two level-1 `may_harm`. The proposals are therefore hypotheses with weak and partly contrary support.

## Second search round (2026-10-08) and editorial exclusions

After the first red-team review, four more searches were run, drug-agnostic and on different questions: interventions in BubR1-hypomorph models and in cells from MVA patients; the published clinical management of MVA and BUB1B; approved drugs for the child's own features (rhabdomyosarcoma, nephrocalcinosis, growth, muscle); and vulnerabilities of aneuploid cells beyond the first round. They made 256 queries. Their sentences were checked by the same code (351 verified, 2 not found and therefore excluded, 6 preprints not checkable; 355 verified once the third round below is added). The free-text consequence tags of these searches were mapped to the tags of the rules by keyword (`norm_consequence` in `literature_summary.py`), and names used by the searches were mapped to the names of the drug universe (`ALIAS`; for example SKI606 is bosutinib, which PMID 39924517 states as \"bosutinib (SKI-606)\"). The same rules R1 to R4 were applied to the new entries.

Six drugs pass the mechanical rules and are **not proposed, by editorial decision of the author, each with its reason shown in the last subsection of section 4 of the report** (`EDITORIAL` in `literature_summary.py`): denosumab, somatropin, mecasermin, niacinamide, eribulin and selinexor (the last is one computational line plus literature, which E9 never accepts as a candidate; it stays in the one-line list of section 6.5). This is a judgement after seeing the results and is declared as such. The round found no approved drug that rescues a BubR1-deficient phenotype: its level-2 findings are genetic or use compounds that are not approved.

## Third round (2026-10-09): contrary evidence for sirolimus, added by the AI assistant

After the second red-team review, four sentences were added by the AI assistant, with the author's approval, to `out/literature_hypotheses/round3_audit.json`, all `may_harm` for sirolimus: two from PMID 35141499 (TORC1 inactivation causes mitotic slippage with chromosome instability in yeast, and mitotic slippage in human cells; level 1), one from PMID 21380626 (growth velocity lower in children with a kidney transplant treated with sirolimus than in matched controls; level 0) and one from the full text of PMID 41820377 (the authors suggest that TOR inhibition might have a negative impact on highly aneuploid neuroblasts; level 1). This was not a search (no queries were run): the first two papers were named by the audit and the third was already cited in the text of the report without a verified sentence. The four sentences pass through the same check against Europe PMC (355 verified in all). They were chosen after seeing the results, so they are a correction of an omission and not an unbiased sample of the contrary literature; a search made to find the contrary evidence of every proposal, with the same care as for the supporting evidence, was not done.

## Verification

Four independent literature searches (one per consequence group) wrote the entries. Every quoted sentence is re-fetched from Europe PMC (abstract or full text) by `track2.verify_literature`, whitespace ignored and HTML entities decoded only. Entries not found are excluded. Preprints whose text cannot be re-fetched are reported as not checkable and are not cited. The report cites only verified sentences.

## Limits stated in advance

- No approved drug is expected to reach level 2 or 3; none did.
- The proposals are hypotheses for investigation and for a clinician, not evidence that any drug works. No dose or treatment recommendation is given and the age of the child is unknown.
- The refuting experiment named for each proposal is a design, not a result.
- Approved by the owner in chat on 2026-10-08 ("APRUEBO").
