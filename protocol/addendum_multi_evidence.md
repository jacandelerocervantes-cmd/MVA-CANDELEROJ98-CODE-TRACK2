# Addendum: multi-evidence integration for proposed candidates (post hoc)

Date: 2026-10-18. Code: `src/track2/multi_evidence.py`. Tests: `tests/test_multi_evidence.py`. Outputs: `out/multi_evidence_mva.json`, `out/multi_evidence_mva.md`.

## 0. Status and purpose

This addendum is **post hoc**. The results of every line (knowledge graph A, phenotype module B1, mitosis module B2, cell-line vulnerability C, signed-direction check D, signature reversal E, the ranking engine) and the nine drugs of the report were seen before this text was written. Everything below is fixed **before** the integration is computed. Nothing here changes an earlier rule, result or file.

Reason. The track asks to characterise the mechanism and to propose candidate medications as hypotheses for follow-up. The earlier rule (at least two usable computational lines) gives zero candidates, which answers "no drug is picked out by two independent validated lines" but does not propose anything. No single line is absolute truth. This addendum asks a different question: given all the evidence of every kind, with its uncertainty, which approved drugs are the best-ordered hypotheses, by age band, and how stable is that order?

The integration **orders hypotheses; it does not create evidence**. It makes no clinical recommendation, gives no dose, and claims no efficacy. A proposed candidate is a hypothesis for a refuting experiment or for a clinician to consider, with its evidence tier shown beside it.

## 1. Candidate pool

A drug is in the pool if (a) and (b) hold:

(a) it is an approved small molecule: molecule type "Small molecule" in `data/universe.json` and approved by the project's check (a Prescription record in Drugs@FDA read from the local cache of `track2.fda_approved.lookup`, or confirmed in `out/fda_approval_check.json`, or flagged as an FDA label in `out/line_e_signature_reversal.json` or `out/engine/apply_mva.json`). A drug with no such record is listed as "approval not checked" and is not in the pool.

(b) it appears in the literature layer (a drug with at least one verified entry with an evidence label, `track2.literature_summary.collect`), or in the top 50 of any of lines A, B1, B2, D, E (BUB1B signature), P (patient signature), or in the top 50 of the engine ranking.

Definition of "top 50" per line: A, the 50 best drugs by best graph rank (distinct drugs); B1, the 50 drugs with the highest robustness count (of 9 configurations); B2, all drugs robust in any of the three weighting schemes (a short list); D, the 50 drugs with the largest net strengthening of the checkpoint; E and P, the 50 drugs with the smallest reversal p-value that are not flagged as hazards; engine, every drug present in the stored engine ranking. Line C has no drug-level result (no supported drug), so it contributes no pool member. **Known gap:** `out/engine/apply_mva.json` stores only the top 20 of the engine ranking and the positions of the nine report drugs; the full ranking is not stored, so "top 50 of the engine" is applied to the stored top 20 and the stored positions.

## 2. Evidence vector

Each drug has up to 11 criteria. Each is a grade in [0, 1] or **missing**, with a direction (supportive, neutral, contrary) and a source file. Grade 0.5 is neutral; for signed criteria values above 0.5 are supportive and below 0.5 contrary; for rank-based criteria a low grade means "not supported", not "contrary".

| Id | Criterion | Grade | Source |
|---|---|---|---|
| M | mechanism fit | see 2.1 | `config/repositioning_rationale.json`, literature layer, Track 1 gene block |
| LIT | literature level and direction | see 2.2 | `out/literature_hypotheses/*.json` (verified entries only) |
| A | knowledge graph | 1 - (best rank - 1) / N scored, times (1 - fraction of random queries that put the drug in their top 50) when measured, times 0.5 when in the top 50 and the null is not measured | `out/kg_query_e5.json`, `out/kg_nullquery_e5.json` |
| B1 | phenotype-module proximity | robustness count / 9 when listed; 0 when scored and unlisted; missing when not scored | `out/phenotype_proximity_e5.json`, `out/phenotype_threshold_scores_e5.json` |
| B2 | mitosis-module proximity | number of the three weighting schemes in which the drug is robust / 3 | `out/track1_scenarios.json` (two-record links, no predicted links) |
| C | cell-line vulnerability | 1 if listed as supported, else 0; missing for all if the line gives no drug-level result | `out/vulnerability_e6.json` |
| D | signed direction | strengthens checkpoint 1; weakens checkpoint 0 (contrary); unclear missing | `out/signed_effect.json` |
| E | signature reversal, BUB1B knockout signature | 1 - q(reversal); 0 and contrary if flagged as mimicking; missing if no library entry | `out/line_e_signature_reversal.json` |
| EP | the same for the patient signature P (one adult patient, 2 arrays against 2 controls; kept separate) | as E | same |
| ENG | ranking engine | 1 - (rank among FDA small molecules - 1) / 1293 when a rank is stored; else missing | `out/engine/apply_mva.json` |
| S | label safety | 1 no flag; 0.5 boxed warning or an entry of `config/proposal_safety.json`; 0 member of a class of the caution list | `data/universe.json`, `config/proposal_safety.json`, `out/caution_list.json` |

Hard exclusions, applied before ranking and listed with the reason: the three hazard drugs of the report (bortezomib, chloroquine and bosutinib, derived as the drugs whose verified entries state selective killing of aneuploid cells, which makes them a hazard for a constitutional mosaic); and any drug whose label contraindicates the hypothesised use, here the growth-hormone and IGF-1 class of the caution list (the label contraindicates active malignancy and the record of the child documents a rhabdomyosarcoma, current status unknown).

### 2.1 Mechanism fit (M)

The mechanism is five consequences of reduced BubR1 function (report section 3.1). The confidence of each consequence is the highest verified evidence level of the consequence itself in the primary literature, as in the table of that section, converted by (level + 1) / 4: consequence 1 (missegregation, mosaic aneuploidy) level 2, 0.75; consequence 2 (proteotoxic and energy stress) level 0, 0.25; consequence 3 (p53 and senescence) level 0, 0.25; consequence 4 (growth and muscle) level 0, 0.25; consequence 5 (tumour predisposition) level 3, 1.0. **Data gap:** `out/differential_diagnosis.json` holds the gene-level differential, not the confidences of the five consequences; the values come from the report table and are stored in the code with this provenance.

A drug addresses a consequence if (i) a verified literature entry of the drug is tagged with it (growth_muscle = 4, proteotoxic_energy = 2, senescence_p53 = 3, aneuploidy_tumour or tumour = 5) or (ii) a rule of `config/repositioning_rationale.json` with a non-empty consequence list matches its mechanism text. Support for the consequence is 1.0 with a "may help" entry, 0.5 for a rule-only match or an unclear entry, 0 when its only entries for it are "may harm". M(BUB1B) = max over consequences of confidence x support; missing when the drug addresses none.

### 2.2 Literature (LIT)

For a drug, count each (paper, direction) pair once at its highest level. net = sum over "may help" pairs of (level + 1) / 4 minus the same sum over "may harm" pairs. LIT = 0.5 + 0.5 tanh(net / 2). Entries that are unclear, no effect or "selectively kills aneuploid cells" add nothing (neutral; the last is handled by the hazard exclusion). Missing when the drug has no verified entry. The shown tier is "level k literature" with k the best help level, plus counts of lines supporting and contrary.

### 2.3 Age-band label support (filter and column, not evidence)

Bands: 0-4, 5-11, 12-17 years (intervals [0,5), [5,12), [12,18)). For each band the state comes from the reviewed label intervals and age-specific cautions of `out/empirical/age_support.json`: **labelled** (the label supports use over the whole band with no age-specific caution), **partly** (over part of it), **not_labelled** (paediatric use not established or adult only), **verify** (the label is mixed or unreadable, needs verification), **not_assessed** (drug not in the label table). Unknown age is shown as the union of the three bands and the age is never assumed. The band is not a criterion of the score. It filters the band shortlist (pass = labelled or partly); verify, not_labelled and not_assessed drugs are listed separately in each band. Where the engine's own band reading differs from this one, the difference is reported.

## 3. Integration: stochastic multicriteria acceptability analysis

10,000 weight vectors are drawn uniformly from the simplex over the 11 criteria (Dirichlet(1,...,1)), seed 20261018. For drug i and weights w, with available set V_i and missing weight share m_i = sum of w over missing criteria:

score_i(w) = [ sum over V_i of w_k g_ik / sum over V_i of w_k ] x (1 - 0.5 m_i).

The first factor is the weighted mean over the available criteria (renormalised); the second is the penalty for the share missing (a drug with half of the weight missing loses a quarter of its score). A drug with all criteria missing scores 0. Hard-excluded drugs are removed before ranking.

Per drug: rank acceptability indices (probability of rank 1, of rank within the top 3, 5 and 10); the **central weight vector** (mean weights over the draws in which the drug is ranked first; if never first, over the draws in which it reaches its best observed rank); the **confidence factor** = mean over the draws of (1 - m_i), the expected share of the weight that rests on data.

## 4. Proposal rule (fixed now)

A drug is a **PROPOSED CANDIDATE for follow-up** if P(top 5) >= 0.5 and it is not hard-excluded. **PRIMARY candidates** are the top 3 by P(top 5) (ties by P(top 3), then name) among the proposed. The evidence tier stays visible. If no drug reaches 0.5 the result is reported as zero proposed candidates; the threshold is not changed after seeing results.

## 5. Sensitivity analyses

(1) Leave-one-criterion-out for each of the 11 criteria (the weights are redrawn over the remaining ones). (2) Monte Carlo on label errors of 10, 20 and 30 % using the corruption of `track2.empirical_label_robustness.corrupt` (levels shifted by one, directions flipped, or both, on randomly chosen verified entries), 200 replicates per rate, with 2,000 weight vectors each; LIT and M are recomputed, the pool and the hard exclusions are kept. (3) Equal weights (one vector, 1/11 each; ties reported). (4) Literature-only weights: SMAA over LIT and M only (both come from the literature layer). (5) Track 1 gene block off (cap 0, section 8). Stability = for each baseline proposed drug, the fraction of replicates or variants in which it remains proposed, the Jaccard overlap of the proposed sets, and, for the MC, the expected number of baseline candidates retained.

## 6. Age-band shortlists

Same weight draws, ranking restricted to the drugs that pass the band filter (labelled or partly), for the bands 0-4, 5-11, 12-17; and the union (unknown age) = the drugs that pass at least one band. The same proposal rule is applied inside each pool. The drugs outside the filter, with their state, are listed.

## 7. Comparison

With the earlier binary rule (`out/four_lines_review2.json`: zero candidates under both readings of B2) and with the nine drugs of the report: for each, rank by P(top 5), proposed or not, excluded or not, and why.

## 8. Track 1 gene block (added by instruction of the user on 2026-10-18, before computing)

The mechanism-fit criterion uses **all genes of the Track 1 result**, not only BUB1B. Gene-level facts only (no positions, alleles or the variants column):

| Gene | Role | Evidence tier | Tier weight |
|---|---|---|---|
| BUB1B | recessive candidate, compound heterozygous: stop-gain with predicted decay + missense of uncertain significance; primary | causal candidate | base of M |
| IGBP1 | hemizygous missense, rank 2 genome-wide, phenotype score 0.141; encodes the alpha4 regulatory subunit of PP2A; literature links it to the PP2A switch that dephosphorylates AKT (PMID 34998390, 38419089), translation initiation (PMID 18625885) and autophagy (PMID 40600667); the effect of the variant is untested; not an MVA gene | hypothesis-grade modifier, insufficient evidence | 1.0 |
| CENPE (MCPH13), MAD1L1 (MVA7) | structural calls not overlapping coding sequence | pathway context only | 0.5 |
| CDC20 and MAD2L1BP (duplications), TTK (deletion), CEP192 (deletion) | megabase-scale events not analysed in Track 1 that overlap checkpoint genes; real segmental copy-number change or caller artefact, undetermined | undetermined | 0.25 |
| DNAH10, IZUMO4, ZNF418 | no mechanistic link found in the differential diagnosis | none | 0 |

Which drug hypotheses each gene bears on (declared before computing): IGBP1 (PP2A, AKT, mTOR, translation, autophagy) bears on the mTOR-axis hypotheses (relevance 1.0 when the mechanism text matches mTOR, FKBP, FK506, AKT or PI3K) and, with relevance 0.5, on drugs with a verified literature entry that mentions mTOR, AMPK, AKT, autophagy or translation (e.g. metformin); the checkpoint context genes and the undetermined events bear on the tumour-context and hazard hypotheses (relevance 1.0 when the mechanism text matches tubulin, mitotic, kinesin, checkpoint, aurora, PLK, TTK, CDC20, centrosome or exportin; 0.5 when a literature entry is tagged tumour or aneuploidy_tumour).

Cap. T(drug) = max over the non-BUB1B genes of tier weight x relevance, in [0, 1]. M = min(1, M(BUB1B) + 0.10 x T): the Track 1 non-BUB1B genes can raise the mechanism-fit criterion by at most 10 % of its scale, and never lower it. The effect of the cap is reported as a sensitivity analysis with and without these genes (cap 0.10 versus 0).

These are **hypotheses to test, not findings**: the karyotype or single-cell sequencing for the megabase events; function tests (PP2A-AKT-mTOR readouts) for the IGBP1 variant.

## 9. Limits

The integration orders hypotheses and cannot create evidence. Pool membership depends on the lines (selection on the same signals that are then integrated), the literature labels were assigned by AI search tasks and not reviewed by a specialist, the weights are uniform over a simplex rather than elicited, the lines have very different validity (A weak, B1 not above its degree baseline, C empty, D and E exploratory, P from one adult patient cell line), and a rank-based grade measures position in a list, not an effect. Label support by age is a statement about the US label, not about benefit. No dose, no treatment recommendation, no claim of efficacy.

## 10. Note added on 2026-10-18 after a first exploratory run (declared deviation of presentation, not of the proposal rule)

A first run (20 Monte Carlo replicates, only to test the code) showed that section 6 as written is permissive: when the ranking is restricted to the 10 to 12 drugs that pass an age-band filter, "top 5" means "the better half of a small pool", and drugs with no relation to the mechanism pass. The proposal rule of section 4 and its threshold are **unchanged** and are applied in the full pool. The **headline band shortlist** is therefore the conservative one: the drugs that are proposed in the full pool and pass the band filter. The restricted-pool reading of section 6 is still computed and shown, for information only. The first run also showed that 153 candidates have no approval record in the project outputs or in the local openFDA cache (the machine is offline): they are listed and are not in the pool, so the pool is incomplete for them, and the result must be repeated with the approval check online. The ranking of all pool drugs by P(top 5) is shown whatever the rule says, so a reader can see the order below the threshold; only drugs that meet the rule are called proposed.

## 11. Pool completion (before rerun) - added 2026-10-18, fixed before the completed run is computed

The first run (outputs `out/multi_evidence_mva.json` and `.md`) stays untouched. The user approved completing three data gaps. The completed run writes `out/multi_evidence_mva_completed.json` and `.md`; the proposal rule (section 4), the weight design and SMAA settings (section 3, 10,000 draws, seed 20261018), the criteria, the sensitivity analyses (section 5) and the Track 1 gene block (section 8) are **unchanged**. Both runs are shown side by side in the completed `.md`.

1. **Approval.** The project's FDA approval check (`track2.fda_approval_check.check`: openFDA / Drugs@FDA match on generic name, substance name or active ingredient; confirmed when an NDA, BLA or ANDA record has a product marketed as Prescription or Over-the-counter; the same code and rules as before, unmodified) is run for the 152 candidates that had no approval record. Network use is limited to openFDA (at most about 3 requests per second, enforced by the check's own delay), and the answers are stored in `out/fda_approval_check_pool_completion.json`. Candidates confirmed as approved small molecules are added to the pool; the others are listed.
2. **Engine.** The frozen B0 engine (`engine_v4.predict_frozen_v4`, frozen model hash verified, no refit) is recomputed offline with the functions of `track2.engine_apply` and its inputs reused as they are (Europe PMC answers from the local caches only; a cache miss stops the run and is reported). Score and rank of every drug of the universe are stored in `out/engine/apply_mva_full_ranking.json` (drug name, score, rank among all drugs and rank among FDA small molecules, no other content). The ENG grade is 1 - (rank among FDA small molecules - 1) / N for every pool drug. Because the full ranking now exists, "top 50 of the engine" of section 1 is applied to it, and those drugs enter the pool if they pass the approval rule of the engine's filter chain (openFDA label match or Drugs@FDA). The engine cannot be checked where no drug is known; it is a literature-and-graph prior.
3. **Line C at drug level.** Only data already on disk (`data/depmap`: PRISM Repurposing 24Q2 primary screen, log fold change of viability at 2.5 uM, and expression); nothing is downloaded. Checkpoint-weakness score of a cell line: the mean over cell lines of the z-scored expression (log2 TPM + 1) of the six checkpoint genes BUB1B, BUB1, BUB3, MAD1L1, MAD2L1 and TTK (`track2.prism_l4.sac_scores`, 'mean_all'); a LOW score means a weak checkpoint. For each pool drug found in PRISM (name or synonym match as in `prism_l4.match_rows`; several PRISM rows of one drug are averaged per cell line): the Spearman rank correlation rho between the drug's viability value and the line's score across cell lines with both values. Sign convention: sensitivity = lower viability, so rho > 0 means the drug is more lethal in weak-checkpoint lines (supportive) and rho < 0 means it is more lethal in strong-checkpoint lines (contrary). Null: 1,000 permutations of the score across cell lines (fixed seed 20261018), two-sided permutation p = (1 + number of permuted |rho| >= observed |rho|) / 1001; Benjamini-Hochberg across the tested drugs of the pool. Graded value: if BH q < 0.10, g = 0.5 + 0.5 clip(rho / 0.3, -1, 1); if q >= 0.10 the association is not distinguishable from the null and g = 0.5 (neutral); **missing** when the drug is not in PRISM or has fewer than 100 cell lines with both values. Declared limits: the association is correlational and comes from cancer cell lines, not from patient cells; the score correlates strongly with proliferation (0.84 in the earlier line C), and this statistic is **not** adjusted for proliferation or lineage, so it can reflect general cytotoxicity of fast-growing lines; it is a drug-level hypothesis signal, not a validated line (the earlier lineage- and proliferation-adjusted line C found no supported drug). Line C is restricted to pool drugs (it is not ranked over all PRISM compounds), so it does not add pool members.
4. **Reporting.** The completed run states plainly whether trametinib survives (still proposed) and whether any new drug is proposed, and gives both runs side by side: proposed set, P(top 5), stability numbers (leave-one-out, label-error Monte Carlo, equal and literature-only weights, gene block off), age-band results and the position of the nine report drugs.


## 12. Candidates by line (before computing) - added 2026-10-18; replaces the proposal rule of section 4 as the main output

Change of design requested by the user. The goal is **not** one drug validated on all parameters. **Each line produces its own candidates**, and the sets are expected to differ. Convergence between lines is reported as information and is never required. The SMAA results of the first run and of the completed run (sections 3 to 11) stay untouched and are shown only as a **cross-check column**. Code: `src/track2/candidates_by_line.py`; outputs: `out/candidates_by_line_mva.json` and `.md`. Every candidate is a **hypothesis for follow-up**; nothing here is a claim of efficacy, a recommendation or a dose.

### 12.1 Common filters, applied to every line before its candidates are chosen

1. **Approved small molecule.** The completed pool of 224 drugs of `out/multi_evidence_mva_completed.json` (molecule type small molecule in `data/universe.json` and approved by the project check: label match, Drugs@FDA, or the completion run of `out/fda_approval_check_pool_completion.json`). A drug that a line finds outside this pool is listed under that line as "not eligible" with the reason (not approved, not a small molecule, not in the universe). The 28 candidates that the approval check could not confirm are listed separately (they may be approved under another name).
2. **Hard safety exclusions** (the drug is not proposed by any line; it is listed under the line that found it, with the reason): (a) the three hazard drugs (bortezomib, chloroquine, bosutinib: verified entries state selective killing of aneuploid cells) and the growth-hormone and IGF-1 class of the caution list; (b) any drug whose label (boxed warning, contraindications or warnings section, read from DailyMed or openFDA through the project cache; network only to those services, only for drugs that a line found) warns about a condition in the phenotype of the child, by the fixed keyword map below. The paediatric-use section is not searched for these keywords. Vitamin D analogues such as calcitriol are excluded by this rule (nephrocalcinosis, HP:0000121).
3. **Cytotoxic chemotherapy** (members of the classes of `out/caution_list.json`: tubulin binders, genotoxic chemotherapy, proteasome inhibitors, checkpoint agents): **not proposed for constitutional use** because of the project's cytotoxic-chemotherapy caution; they stay listed under the line that found them as "context/hazard".

Keyword map (case-insensitive, searched in boxed warning, contraindications and warnings only), by documented term of the child. Nephrocalcinosis HP:0000121: nephrocalcinosis, hypercalcemia, hypercalcaemia, hypercalciuria, renal calcification. Short stature HP:0004322, failure to thrive HP:0001508, small for gestational age HP:0001518: growth retardation, growth suppression, growth inhibition, decreased growth, reduced growth, impaired growth, stunted growth, failure to thrive. Skeletal muscle atrophy HP:0003202: muscle atrophy, muscle wasting, myopathy, rhabdomyolysis. Rhabdomyosarcoma HP:0002859: rhabdomyosarcoma anywhere, and, only in the contraindications section, active malignancy or malignant neoplasm. Premature birth HP:0001622 and recurrent spontaneous abortion HP:0200067 (parental history): no keyword (not conditions that a label warns about for the child). A label that cannot be read is flagged "label not read" and is not excluded (it is never read as "no warning").

### 12.2 Native criterion of each line (only the line's own output; fixed now)

- **Literature (research line):** rules R1 to R4 of `track2.empirical_label_robustness.decide` on the verified entries (a constitutional hypothesis needs a verified finding of benefit at level 1 or more on a feature of the child or on consequence 2 or 3; a tumour-context hypothesis a verified level-1 finding in an aneuploid or tumour system; a class sibling is listed with its representative; level-3 harm or a flagged finding goes to the caution list). Ranked by best help level, then number of papers with a help finding, then net (help minus harm) of section 2.2.
- **Engine:** the first 10 drugs of the frozen B0 ranking among FDA small molecules (`out/engine/apply_mva_full_ranking.json`) that pass the common filters, ranked by engine rank (native score = the engine score). A prior from literature co-mention and graph targets, not evidence.
- **A (knowledge graph):** the line's usable set after its generic-drug filter: the top 50 distinct drugs of `out/kg_query_e5.json` whose fraction of random queries placing them in their top 50 is below 5 % (`out/kg_nullquery_e5.json`); drugs whose fraction is not measured are not counted. Ranked by graph rank.
- **B1 (phenotype module):** the line's own support set: drugs robust in at least 7 of 9 configurations (`out/phenotype_proximity_e5.json`). B1 did not pass its validation (no excess over a degree baseline); if the set is empty the line has no candidate and the best drugs below the threshold are shown as information only.
- **B2 (network proximity, mitosis module):** the line's own support set, scheme W3 (pathogenic gene only), two-record links, no predicted links (`out/track1_scenarios.json`). B2 is held exploratory; the set has no score, so it is ordered by name.
- **C (PRISM drug-level vulnerability):** BH q < 0.10 with the correct direction (rho > 0: more lethal in weak-checkpoint lines) from `out/line_c_drug_level.json` (section 11.3), ranked by q then rho. Correlational, from cancer cell lines, not adjusted for proliferation or lineage.
- **D (signed direction):** drugs called "strengthens checkpoint" by `out/signed_effect.json`, ranked by net signed effect. Exploratory (4 of 6 controls recovered, two by construction).
- **E (signature reversal), BUB1B knockout signature; E-P, patient signature (GSE22206); E-T1, Track 1 combined signature:** FDR < 0.1 for reversal and not flagged as mimicking (`usable` in `out/line_e_signature_reversal.json`), ranked by reversal p-value. E-P comes from one adult patient cell line (2 arrays against 2 controls).

When more than 5 candidates qualify in a line after the filters, the top 5 by the line's own score are shown and the rest are listed as also qualifying. A line with no qualifying drug reports an empty set, which is a result.

### 12.3 Per candidate, and the other outputs

Each candidate carries: line, native score and q or rank; evidence tier (best literature level and number of papers for and against, and the line limit in one phrase); age-band label support (0-4, 5-11, 12-17; unknown age = union; not_labelled, verify and not_assessed are valid states); safety flags (boxed warning, `config/proposal_safety.json` entries, cautions, label keyword matches or "label not read"); which of the five consequences it addresses (mechanism fit of section 2.1) and the Track 1 genes that bear on it (section 8, a hypothesis to test, not a finding); the type of experiment that would refute it (one phrase per line, to be expanded in `protocol/follow_up_experiments.md`); and the multi-evidence SMAA probability P(top 5) of the completed run, as a cross-check only.

Outputs: one table per line; the **convergence matrix** (drug by line, candidate or not); the overlap between lines against chance (hypergeometric test over the completed pool of 224, for each pair of lines with at least one candidate; no multiplicity claim is made); the number of candidates unique to one line; the age-band view; the cross-check (SMAA probabilities, and the stability results of the completed run: leave-one-out, label-error Monte Carlo, equal and literature-only weights, gene block off).

### 12.4 Implementation note (added 2026-10-18 after a first test run of the code, before the final run)

The first test run matched the keyword "myopathy" inside "cardiomyopathy" and "encephalomyopathy" (trametinib, metformin). Keywords are therefore matched as whole words, which is the plain reading of the map in 12.1; the keyword lists themselves are unchanged. The label name of phenylbutanoic acid is read as "sodium phenylbutyrate" (the marketed product; the same alias as the approval check). The screen is crude by design: a keyword in a warnings section can refer to an animal finding or to a different population (for example reduced growth rates in pups), and the matched sentence is stored with each exclusion so that a reader can overrule it. No other rule changed.

## 13. Corrections to the candidates-by-line screen (before recomputing) - added 2026-10-18

The first output `out/candidates_by_line_mva.json` and `.md` stays untouched. The corrected output is written with the suffix `_v2`. Per-line criteria (12.2), the top-5 rule, the convergence matrix, the overlap statistics and the SMAA cross-check are unchanged and recomputed. Two corrections, requested by the user.

### 13.1 Safety screen: the project's hand-read layers are the authority; keywords only flag

The whole-word keyword exclusion of 12.1 and 12.4 is replaced. Authority, in this order: `config/proposal_safety.json` (verbatim label statements read by hand), `out/caution_list.json` (classes and quoted label sentences), the verdicts of `src/track2/label_levels.py` (level 6 = not indicated or contraindicated in children) and `config/pediatric_label_review.json`.

**Hard exclusions are only:** (a) the three hazard drugs (bortezomib, chloroquine, bosutinib); (b) drugs that the project's own layers already exclude: the growth-hormone and IGF-1 class of the caution list and drugs with label level 6; (c) a documented child-specific contraindication tied to the child's phenotype terms, namely a **specific** keyword found in the **contraindications or warnings** section of the label: for nephrocalcinosis HP:0000121 the words nephrocalcinosis, hypercalcemia, hypercalcaemia, hypercalciuria, renal calcification; for rhabdomyosarcoma HP:0002859 the word rhabdomyosarcoma anywhere in those sections, and active malignancy or malignant neoplasm in the contraindications section. The quoted sentence is stored for human review. Vitamin D analogues such as calcitriol stay excluded by (c).

**Everything else is a flag, never an exclusion**: the broad keywords of 12.1 (growth suppression and the other growth phrases; muscle atrophy, muscle wasting, myopathy, rhabdomyolysis) in the boxed warning, contraindications or warnings sections are stored with the matched sentence as "needs manual verification". A drug covered by `config/proposal_safety.json` is not screened by keywords at all: its hand-read label warnings are shown instead (this keeps sirolimus and selumetinib as hypotheses with their label warnings shown, as the report treated them, unless a project layer excludes them). Cytotoxic classes of the caution list stay "context/hazard" as in 12.1(3). Pediatric-use sentences are not searched. A label that cannot be read is flagged "label not read" and not excluded.

### 13.2 Age bands for every candidate, with the project's functions

Every candidate of every line gets a reading for 0-4, 5-11 and 12-17 years (unknown age = union) from `track2.engine_case_study.label_support`, which uses `track2.age_label_table.row_for`, the manual review of `config/pediatric_label_review.json` when it exists, and otherwise the heuristic classifier `labels.classify_pediatric`. The label comes from DailyMed (set id of the manual review when there is one) or openFDA through the project cache; network only to those services, only for candidates. The label section used is quoted (the verbatim paediatric-use sentence). The state `not_assessed` remains only when the label cannot be retrieved, and then it is written as "label not retrievable". States: labelled, partly, not_labelled, verify, age_not_stated (shown as labelled, with the note) and the classifier outcome is marked as heuristic where no manual review exists. The age of the patient is never assumed.

### 13.3 Reporting

The `_v2` markdown states what changed relative to the first output, per line (added and removed candidates, changed age readings) and says explicitly what happened to sirolimus and selumetinib.
