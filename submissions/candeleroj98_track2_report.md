# Predicting candidates for repurposing FDA-approved medicines from multiple independent lines of evidence in *BUB1B*-associated mosaic variegated aneuploidy

**MVA Hackathon 2026 · Track 2 – Drug Repurposing** · Juan Antonio Candelero Cervantes (single author) · Code: https://github.com/jacandelerocervantes-cmd/MVA-CANDELEROJ98-CODE-TRACK2 · Licence: CC BY 4.0

---

## Summary

The patient carries two rare changes in *BUB1B*, a gene linked to mosaic variegated aneuploidy, a rare condition in which some of the body's cells carry the wrong number of chromosomes. No approved medicine is known to restore what the gene does, so medicines already approved for other uses are a practical place to look.

Several independent lines of evidence were examined side by side, namely the literature, knowledge networks, drug screens and gene expression, so that a line finding nothing is still informative. Further computational tests discarded about half of the suggestions, and a few medicines were set aside as hazards, because they are selectively toxic to cells with the wrong number of chromosomes or disruptive to cell division.

The result is a set of testable leads, each with its evidence level stated, produced by a method designed to be applied to other diseases; how it behaves when little has been published was tested only by thinning the literature of better-studied ones.

The research that would follow starts with a karyotype or single-cell sequencing, to find out whether cells with the wrong number of chromosomes are present at all. Only if they are would the surviving leads be worth testing in a cell model, and a stop rule ends the sequence if they are not.

---

## 1 Methods

### 1.1 Data

- A gene-level, categorical extract of the Track 1 ranked candidate list: seven genes with their consequence classes, without positions, alleles, genotype calls or read metrics.
- Eight Human Phenotype Ontology (HPO) [1] terms, seven describing the patient's clinical symptomatology and one the parents' history.
- Variant and gene resources: ClinVar [2], UniProt [3], the AlphaFold Protein Structure Database [4], and gene lists for mosaic variegated aneuploidy (MVA) from HPO [1] and Orphanet [5].
- Network and pathway resources: PrimeKG [6]; STRING v12.0 physical interactions [7]; Reactome 97 [8], for the pathway membership used to seed the gene modules; ChEMBL 37 [9], for approved medicines with their curated mechanisms and targets, which define the drug universe; and OmniPath [10].
- Perturbation and screening resources: DepMap 24Q2 CRISPR gene effect [11], PRISM Repurposing 24Q2 [12], LINCS L1000 consensus libraries [13], GDSC2 8.5 [14], and the public expression series GSE22206 with its platform annotation from the Gene Expression Omnibus [15].
- Literature retrieved and verified through Europe PMC [16].
- United States drug labels and approval records from DailyMed [17], openFDA [18] and Drugs@FDA [19].

### 1.2 Mechanism characterisation

The stop-gain variant was assessed against the rule for nonsense-mediated decay summarised by Abou Tayoun *et al.* [20]. The missense variant was characterised from public sources only. Predictor scores [21–25] were read against the ClinGen calibrations [26,28] and interpreted under the ACMG/AMP framework [27]; predictor strengths were not added together.

Allelic configuration was treated as an explicit inference with stated assumptions rather than as an observation.

Reduced BUBR1 function was then described as five downstream consequences, each supported by sentences retrieved from the literature and graded on a four-level scale: level 0 inferred or indirect; level 1 shown in a related system, such as aneuploid cells, yeast or an animal with a checkpoint defect; level 2 shown in human *BUB1B*- or MVA-derived cells, or in a hypomorphic model with rescue; and level 3 shown by clinical observation in MVA or *BUB1B*.

The search tasks assigned the levels from that written definition, and each consequence enters the mechanism-fit criterion with a weight of (level + 1) / 4.

### 1.3 Literature layer

Eight searches assisted by a language model collected sentences about approved medicines in two rounds: four organised by consequence (419 queries) and four worded without reference to any medicine (256 queries). Custom code then re-fetched every quoted sentence from Europe PMC [16] and required a word-for-word match: of 363 sentence slots, 355 were found; 6 preprint sentences could not be checked and were not cited; and 2 were not found and were excluded.

Four selection rules were applied: the medicine must be approved; a constitutional hypothesis requires a benefit finding at level 1 or above on an element of the patient's clinical symptomatology or on the cellular consequences of aneuploidy, while a tumour-context hypothesis requires a verified level-1 finding in an aneuploid or tumour system; a class sibling is listed together with its representative; and a level-3 harm finding sends a medicine to the caution list.

### 1.4 Computational lines

Five lines of evidence were defined, each with its own criterion. A ranking model was applied separately (Section 1.7).

- **Line A, knowledge graph.** Ranking of medicines by personalised PageRank on PrimeKG [6], with side-effect edges removed.
- **Line B, network proximity.** Distance on STRING v12.0 physical interactions [7], following Guney *et al.* [29] and Cheng *et al.* [30], in two forms: B1 against a module built from the patient's clinical symptomatology, and B2 against a module built around BUBR1. The modules were seeded from Reactome pathway membership [8].
- **Line C, cell-line sensitivity.** Partial rank correlation between medicine sensitivity in PRISM Repurposing 24Q2 [12] and a checkpoint-weakness score derived from DepMap 24Q2 [11], adjusted for cell lineage and a proliferation covariate across 895 cell lines and 29 lineages.
- **Line D, signed direction.** Net signed effect over the OmniPath signed network [10].
- **Line E, signature reversal.** Reversal of expression signatures in the LINCS L1000 consensus libraries [13], once against a *BUB1B* knockout signature and once against a signature built from the public series GSE22206 [15].

Lines A and B1 had to show an excess over a degree-matched baseline before being used, and the protocol required a candidate to have support from at least two usable lines.

### 1.5 Candidate assembly and cross-check

The drug universe was taken from ChEMBL [9] and restricted to medicines approved by the Food and Drug Administration: all 71 medicines in the study were confirmed in Drugs@FDA [19], as were 124 of the 152 further candidates checked for the wider pool. The patient's age is unknown and was not assumed, and no dose, schedule or treatment is proposed in this report.

Each line put forward candidates by its own criterion. Every entry carries its evidence level, a reading of the United States label and any safety flag, and agreement between lines is reported as information only.

Label text was read for three fixed age bands, 0–4, 5–11 and 12–17 years, and classified in six levels; code checked that each of the 551 quoted label sentences appears in the label it cites. Because the patient's age is unknown, the union of the three readings is also given. A caution list was built from five classes fixed in advance, and each of its entries requires a verbatim label sentence.

A cross-check, fixed before computation, drew 10,000 weight sets over 11 criteria to order the hypotheses.

### 1.6 Refutation checks

Nine checks were specified before any was run, to test whether each hypothesis could be contradicted with data already available (Annex Table A1). Passing a check establishes only that the available data do not contradict the hypothesis.

### 1.7 Ranking model and its assessment

A literature-grounded ranking model combines counts of co-mentions of a medicine and disease genes, retrieved through Europe PMC [16], with a prior derived from PrimeKG [6]; its logistic weights were fitted by disease-grouped cross-validation.

It was assessed against a popularity baseline on rare diseases whose known medicines were taken from PrimeKG indications [6], drug labels [17,18] and trial arms; the 42 held-out diseases were scored once. Two further analyses thinned the retrieved literature and separated the contribution of the prior from that of the literature.

The frozen model was then applied, without refitting, to MVA and to Fanconi anaemia.

## 2 Results

### 2.1 Mechanism of reduced BUBR1 function

BUBR1 is a component of the mitotic checkpoint complex of the spindle assembly checkpoint; reduced function allows anaphase to begin before all kinetochores are attached. Of the seven genes in the gene-level extract, *BUB1B* was the only one taken forward as the primary mechanism (Table 1).

Four further genes, overlapped by megabase-scale events that Track 1 did not analyse, were added to the mechanism-fit gene block, which brings it to eleven genes; two more genes that Track 1 named in those events are not included (Table 1). Of the ten MVA genes curated or proposed in HPO [1] and Orphanet [5], only *BUB1B* and *MAD1L1* are among the seven of the extract, while two others sit in the gene block.

**Table 1.** Genes entering the analysis, with their consequence class and their role.

| Gene | Consequence class | Role in this analysis |
|---|---|---|
| *BUB1B* | One loss-of-function and one missense variant | Primary mechanism |
| *MAD1L1* | Structural call, not overlapping coding sequence | Pathway context |
| *CENPE* | Structural call, not overlapping coding sequence | Pathway context |
| *IGBP1* | One missense variant | Candidate modifier, no functional test |
| *DNAH10* | Two missense variants | Not modellable in the protein networks used |
| *IZUMO4* | One in-frame and one missense variant | No further role |
| *ZNF418* | One missense variant | No further role |
| *CDC20*, *MAD2L1BP*, *TTK* and *CEP192* in the gene block; *SHLD3* and *CEP63* outside it | Megabase-scale calls | Undetermined until a karyotype or single-cell sequencing |

Both *BUB1B* variants were read at protein level only (Table 2).

**Table 2.** The two *BUB1B* variants at protein level, with the measurements behind each reading.

| Variant | Measurements | Reading |
|---|---|---|
| Stop-gain | Stop in exon 17 of 23, well before the penultimate exon, so decay is predicted [20]; any product would lack 30 % of the protein, including the kinase-like domain [3] | Loss of function, predicted but not tested |
| Missense at a conserved asparagine | Conserved in 116 of 120 aligned orthologues; inside the kinase-like domain [3]; buried in the AlphaFold model [4]; nearest reported missense ten residues away, listed in UniProt [3] as a variant in MVA1, while the only ClinVar [2] missense classified as pathogenic or likely pathogenic in the gene lies far from it; no regional signal in the ±20-residue window (0 pathogenic or likely pathogenic against 2 benign or likely benign, Fisher p 1; 1 against 2, p 0.47, if the UniProt MVA1 variants are counted as pathogenic) | Uncertain significance [27]. Calibrated strengths disagree across predictors. The calibration of Pejaver *et al.* [26] did not include AlphaMissense, and under it alone only PolyPhen-2 0.997 [22] reaches PP3 supporting; the later ClinGen calibration [28] places AlphaMissense 0.923 [21] at PP3 moderate. REVEL 0.472 [23], CADD 24.5 [24] and SIFT 0.01 [25] fall in indeterminate gaps. The tools are correlated, so strengths were not added, and the most favourable reading still leaves the variant of uncertain significance |

Allelic configuration is undetermined. Under stated assumptions, the probability that the two variants lie on different copies of the gene is close to 1, but it falls to 0.5 if *BUB1B* does not explain the presentation. Parental genotyping, or long reads spanning both variants, would decide it.

MVA (OMIM 257300 [31]) was first linked to *BUB1B* by Hanks *et al.* [32]. The reported pattern is a nonsense allele combined with a kinase-domain missense allele, and the patient's combination is consistent with it [33]. Whether the missense variant lowers BUBR1 function is not known; one class of MVA missense alleles located in or near the kinase domain lowers BUBR1 abundance by increasing protein turnover [33].

The five consequences were selected because they bear on the patient's clinical symptomatology (Table 3).

**Table 3.** Consequences of reduced BUBR1 function, with the best evidence level reached and the state of each in the patient.

| Consequence | Best evidence level | State in the patient |
|---|---|---|
| Chromosome missegregation and mosaic aneuploidy | 2 | Not demonstrated: a coarse read-depth screen cannot detect low-fraction mosaicism and is not a karyotype |
| Proteotoxic and energy stress; p53- and p16-dependent senescence; reduced growth and muscle loss | 0 | Not measured |
| Tumour predisposition | 3 | Not assessable |

The three consequences graded level 0 rest on yeast, mouse and cell-line data rather than on human observations; the senescence consequence would carry 0.75 using the level that the same study received in the second search round.

One element of the patient's clinical symptomatology and the parental term have no reported link to *BUB1B*, and the rationale assumes that cells with the wrong number of chromosomes are present, so a karyotype or single-cell sequencing is the first test (Section 3).

Were no excess of aneuploidy found in every tissue examined, and the missense variant to behave as wild type or lie in cis with the stop-gain, stop rule S0 would halt the medicine experiments; the only argument left would be reduced growth and muscle loss, which rests mainly on mouse data.

### 2.2 Candidates proposed by each line

The two-line rule gave no candidate; line C did not pass the rule at gene level, and the drug-level form that Table 4 reports was added later. Under either reading of line B2 the rule gave zero candidates, and the overlap is zero in all six pairs among lines A, B1, B2 and C (1,076 scorable medicines, hypergeometric p 1.0). A class-level test found one class, multi-target kinase inhibitors, present in two usable lines, a result compatible with chance.

**Table 4.** The lines of evidence, their usability and the size of the candidate set each proposes. The criterion of each line is given in Sections 1.3, 1.4 and 1.7.

| Line | Usability | Candidates proposed |
|---|---|---|
| Literature | 355 of 363 sentence slots verified word for word | 5 of 11 |
| A, knowledge graph | Usable: mean area under the curve 0.802 (0.782–0.821) over 200 diseases, excess over a degree-matched baseline 0.022; the patient ranking is generic | 2 of 3 |
| B1, symptom module | Not usable: excess over baseline −0.083 at the lower bound | 0 |
| B2, BUBR1 module | Passes under its module rule: 0.777 (0.687–0.859); too few diseases remained for the specificity check to be conclusive | 3 of 4 |
| C, cell-line sensitivity | 6 of 23 unadjusted associations survive adjustment | 5 of 9 |
| D, signed direction | Exploratory: 4 of 6 control perturbations recovered, two of them by construction | 5 of 335 |
| E, knockout signature | 2 of 6 control targets recovered; no usable medicine | 0 |
| E, patient-cell signature | 18 usable medicines, from lymphoblastoid arrays of one unrelated adult | 5 of 18 |
| Ranking model | Exploratory: it reads mentions, not direction of effect | 5 of 1,293 |

Candidates proposed: the five highest-placed medicines a line puts forward among the medicines it can score. Further medicines met a line's criterion without entering that first five: selumetinib in the literature line, and five more in the ranking model.

For the literature line, the best evidence level reached is 1, and the medicine-agnostic round found no approved medicine that restores a BUBR1-deficient phenotype. In line A, only pomalidomide and disopyramide are specific to this patient's input. In line C, the unadjusted score tracked proliferation closely, and trametinib, prednisolone and cobimetinib did not survive adjustment.

The lines put forward 30 entries, 28 of them distinct medicines, 7 with label support in at least one age band. Two medicines appear in two lines, neratinib (C and E) and trametinib (literature and C), and neither keeps two lines after the checks. Table 5 lists the entries that remain open, and Annex Table A2 those that were closed or kept as hazards.

Each candidate is tied to the part of the mechanism it was selected for. Sirolimus, everolimus, metformin and sodium phenylbutyrate were selected for proteotoxic and energy stress, senescence and reduced growth, the three consequences graded level 0. Trametinib and selumetinib were selected for tumour predisposition, graded level 3. The candidates from the computational lines were selected from the gene and its network rather than from any graded consequence.

The link from mechanism to candidate is therefore weakest where the evidence for the consequence is strongest: missegregation and mosaic aneuploidy, the consequence graded level 2, yields no open candidate of its own, and the medicines that act on cells with the wrong number of chromosomes are the ones held back as hazards.

**Table 5.** Medicines that remain open after the checks, with the line that proposes each one.

| Medicine | Line | Score in that line | Literature level | Label reading 0–4 / 5–11 / 12–17 years | Safety flag | Refutation check |
|---|---|---|---|---|---|---|
| Sirolimus | Literature | Net −0.5 | 1 (7 benefit, 9 harm papers) | N / N / P (union N–P) | Boxed warning; growth, wound-healing and renal warnings | Check 8 withdrew only the *IGBP1* argument |
| Everolimus | Literature, class sibling | Net 0.5 | 0 | P / L / L (union P–L) | Boxed warning | Check 8 withdrew only the *IGBP1* argument |
| Metformin | Literature | Net −0.25 | 1 (1 benefit, 2 harm papers) | N / P / L (union N–L) | Boxed warning | Check 8 withdrew only the *IGBP1* argument |
| Sodium phenylbutyrate | Literature | Net 1.0 | 1 (2 benefit, 0 harm papers) | L / L / L (union L) | Sodium retention; tablet weight limit of 20 kg | None applicable |
| Trametinib | Literature; C unadjusted only | Net 0.75 | 1 (2 benefit, 0 harm papers) | P / L / L (union P–L) | Label warning; concordant hazard flag in the patient-cell signature | Its line C support did not survive adjustment |
| Selumetinib | Literature, class sibling | Qualifies in the literature line | 0 | Labelled from 1 year, not split by band | Preclinical chromosome-level label finding | Not examined by the checks; no laboratory step |
| Lazertinib | C | Adjusted correlation 0.140, q 0.0316 | No entry | N / N / N (union N) | None flagged | Check 1 passed; not estimable in the independent screen |
| Ceritinib | B2 | Robust under one weighting scheme | No entry | N / N / N (union N) | None flagged | Check 5 passed (p 0.0099) |
| Selinexor | B2 | Robust under one weighting scheme | 1 (1 benefit, 0 harm papers) | N / N / N (union N) | Label warning | Check 5 passed (p 0.0099) |
| Tirbanibulin | B2 | Robust under one weighting scheme | No entry | N / N / N (union N) | None flagged | Check 5 passed (p 0.0198) |
| Pomalidomide | A | Graph rank 43; random-query share 0.006 | No entry | N / N / N (union N) | Boxed warning | Check 4 passed, a permissive test |
| Disopyramide | A | Graph rank 49; random-query share 0.036 | No entry | N / N / N (union N) | Boxed warning; growth to verify | Check 4 passed, a permissive test |
| Nirogacestat | D | Net signed effect 49 | No entry | N / N / N (union N) | None flagged | Check 6 passed (p 0.003) |

Graph ranks of pomalidomide and disopyramide are their positions in the line A candidate list; the hazard ranks in Table A2 are positions in the combined graph ranking of 2,128 scorable medicines, where pomalidomide is 36 and disopyramide 41. The two lists differ for reasons that are not recorded, and ranks should not be compared across the two tables.

Label reading: L labelled, P partly labelled, N not labelled, V to verify. Readings were made by a language model and describe the populations studied for other diseases, not suitability for this patient. No entry has been tested in the laboratory.

Bortezomib, chloroquine and bosutinib are not proposed because their sources describe selective killing of cells with the wrong number of chromosomes; together with the four tubulin binders, they serve only to design harm experiments.

### 2.3 Cross-check runs

The cross-check was run three times and gave no single answer; the reasons lie in pool composition and missing-data handling rather than in evidence (Table 6).

**Table 6.** The three cross-check runs.

| Run | Pool | Leading medicine and its probability of ranking in the top five | Reason for the difference |
|---|---|---|---|
| 1 | 84 medicines | Trametinib, 0.50 | Small pool, so a top-five position is permissive |
| 2 | 224 medicines, after approval checks | Calcitriol, 0.74 (excluded; Annex Table A2); trametinib 0.42 | Missing criteria renormalised, which favoured a medicine without verified literature |
| 3 | 224 medicines, reported by line | No single leading medicine, by design | Each line proposes by its own criterion |

### 2.4 Outcome of the refutation checks

Of the 32 entries examined, 16 were closed by these checks, 7 remain open after passing a check, 5 literature medicines had no check that tests the medicine itself, and 4 tubulin binders were kept to design harm experiments. Those 32 are the 28 distinct medicines the lines put forward together with the four tubulin binders, which no line proposed. Annex Table A1 gives the question each check asks and its outcome. Table 5 also includes selumetinib, which qualified in the literature line but was not among the 32 entries.

Three checks changed the reading of the evidence. The enrichment of tubulin binders among the line C survivors showed that the line behaves as a general mitotic-dependency control rather than a *BUB1B*-specific one. The network argument for *IGBP1* was withdrawn, which does not affect the medicines themselves. The one attempt at independent replication failed: only paclitaxel could be estimated in GDSC2 [14], its signal did not replicate, and that check did not pass its own control.

### 2.5 Performance of the ranking model on other rare diseases

The ranking model was assessed on rare diseases other than MVA in five rounds (Table 7).

**Table 7.** Assessment of the ranking model against a popularity baseline, by round.

| Round | Diseases | Gain in area under the curve over the baseline | Reading |
|---|---|---|---|
| 1 | 17, fixed in advance | 0.206 (0.049–0.343) | Criteria met; the reference set of known medicines is partly circular |
| 2 | 28 new, original reference set | 0.00 (−0.09 to +0.08) | The original weights do not transfer cleanly to new diseases |
| 3 | 14 held out | +0.02 (−0.07 to +0.10) against the round-1 weights | No improvement shown; the round-1 weights were kept |
| 4 | 42 held out, of 145 | 0.19 (0.10–0.28) for the simple model | Supported. The embedding variant gained 0.20 but its calibration slope of 1.237 lies outside the accepted range, so its ranking stays exploratory |
| 5 | 42 held out, reduced reading design | Over the prior-only model, 0.005 (0.002–0.009) at 3 % of records; 0.002 at 1 % | Inconclusive: only one of the three designs was affordable, and coverage reached 1.0 % of pairs |

Thinning the retrieved literature to 10, 3 and 1 % of records reduced the advantage from 0.19 to 0.12, 0.09 and 0.07, without taking it below zero. The prior-only model gives 0.06 at every fraction and the literature-only model falls below the baseline, so the advantage at low volume comes from the knowledge-graph prior.

In a worked case of another rare disease with three known medicines, the model placed them at ranks 89, 416 and 948 of 1,646. Applied without refitting to MVA, its ten highest-ranked medicines are olaparib, tyrosine, carboplatin, progesterone, zinc chloride, cabazitaxel, vorinostat, fingolimod, arsenic trioxide and cisplatin, of which six carry a cytotoxic caution and none has paediatric label support; applied to Fanconi anaemia, seven of ten are cytotoxic.

The model reads mentions, not direction of effect, so a mention of a cytotoxic medicine often describes sensitivity or toxicity. It is reported as a ranking procedure and not as a candidate list.

## 3 Experimental design

The study is a single-case, retrospective analysis that starts from a gene-level extract of the Track 1 ranked candidate list and ends in a set of testable leads for laboratory work. Computation comes first, and laboratory work is reserved for what computation cannot contradict.

Five computational lines and the literature were examined side by side, each proposing by its own criterion, and every line is reported whether or not it found anything. Code verified each quoted literature sentence and each quoted label sentence, and the refutation checks were fixed before any of them was run.

**Table 8.** Stages of the analysis, with the resource, the rule applied and the output of each.

| Stage | Resource | Rule and parameters | Output |
|---|---|---|---|
| 1 Input | Gene-level extract of the Track 1 ranked candidate list | Seven genes, plus four in the gene block; modifier weight capped at 0.1 | *BUB1B* as the primary mechanism |
| 2 Variant reading | ClinVar [2], UniProt [3], AlphaFold [4], predictors [21–25] | Decay rule [20]; ClinGen calibrations [26,28]; ACMG/AMP framework [27] | Loss of function predicted; missense of uncertain significance |
| 3 Configuration | — | Explicit inference with prior 0.5, penetrance 0.9, alternative cause 1 × 10⁻⁶ | Posterior close to 1, conditional on the assumptions |
| 4 Consequences | Europe PMC [16] | Five consequences; levels 0 to 3 | Levels 2 and 3 reached |
| 5 Literature layer | Europe PMC [16] | 675 queries; 363 sentence slots; four selection rules | 9 medicines in three groups |
| 6 Label reading | DailyMed [17], openFDA [18] | Three age bands; six levels | 551 of 551 quotes found |
| 7 Approval check | Drugs@FDA [19] | Application and a marketed product required | 71 of 71; 124 of 152 in the wider pool |
| 8 Lines of evidence | PrimeKG [6], STRING [7], PRISM [12], DepMap [11], OmniPath [10], LINCS [13], GEO [15] | Degree-matched baseline for A and B1; two usable lines for a candidate | Zero candidates under the two-line rule |
| 9 Candidates by line | — | Five highest-placed medicines per line | 28 distinct medicines |
| 10 Cross-check | — | Fixed before computation | Three runs, no evidence created |
| 11 Refutation checks | PrimeKG [6], GDSC2 [14], STRING [7], LINCS [13] | Nine checks; 1,000 permutations each where applicable | 16 of 32 entries closed |
| 12 Laboratory plan | — | Stop rule S0 first; power analysis for sizing | Plan in Annex Table A3, none of it carried out |

The parameters of each stage are listed in the repository configuration files.

Sample sizes are planning numbers from a power analysis of independent cells: at 80 % power and a 5 % significance level, 348 cells per arm would detect a halving of a 12 % baseline fraction of aneuploid cells, 1,632 a quarter reduction, and 150 and 84 cells per arm a halving from baselines of 25 % and 39 %.

Those figures are approximate and optimistic. An exact test asks for 380 cells per arm in the first case, and if cells are clustered by culture, 25 to 110 cultures per arm are needed, or 1,250 to 5,500 cells; the larger figure is the one to use.

Concentration ranges are a design choice, and a null result at low power would not contradict a hypothesis.

## 4 Limitations and future work

### 4.1 Limitations of the analysis

- No laboratory experiment of any kind has been carried out, so no entry has been tested.
- This is a single-case analysis, so its performance cannot be generalised.
- The best evidence for any medicine is level 1, a related system; nothing reaches the level of *BUB1B*- or MVA-derived human cells, and no approved medicine is known to restore BUBR1 function.
- The search queries were not independent of the medicines eventually proposed, and the four selection rules were written after the literature had been read.
- The searches were run by language-model search tasks, so the queries, and therefore the literature retrieved, depend on the model used; a different model could have retrieved a different set of papers, and this was not tested.
- The arrangement by which each line proposes its own candidates was adopted after the two-line rule had returned nothing, so the candidate sets order hypotheses rather than test them.
- The knowledge-graph ranking is largely generic.
- The lines are not fully independent: line A and the ranking model share the same knowledge-graph prior, and the literature line and the evidence levels rest on the same language-model reading.
- Line C behaves as a general mitotic-dependency control, and the independent screen does not contain lazertinib, neratinib, vincristine, colchicine or eribulin, so no survivor of that line has independent confirmation.
- The patient-cell signature rests on lymphoblastoid arrays of one unrelated adult, in a cell type that is not the tissue of interest, and no medicine in it survived the leave-one-array-out check.
- The ranking model's advantage at low literature volume comes from the knowledge-graph prior, not from the literature.
- Six medicines that passed the selection rules were set aside editorially, with the reasons recorded in the repository.
- Label and evidence readings were made by language models from written definitions, without review by a pharmacist or a specialist.
- Mosaic aneuploidy, allelic configuration and the effect of the missense variant are untested, the megabase-scale calls remain undetermined, and the *IGBP1* axis has no functional test.
- No laboratory design exists yet for pomalidomide, disopyramide, nirogacestat or selumetinib, and whether to reopen the neratinib entry is a judgement about how much weight the patient-cell signature deserves.
- Twenty-eight of the 152 further candidates checked have unconfirmed approval status and are listed rather than proposed.
- The evidence levels were not applied uniformly. Observational statements in related systems, such as yeast and mice, were graded level 0 although the definition lists such systems at level 1, and one study received different levels in two search rounds. The grades of the three level-0 consequences are conservative and unreviewed.
- The cross-check creates no evidence, and its three runs are not comparable to one another.

### 4.2 Analyses to strengthen the result

Three analyses would strengthen the result most directly. A karyotype or single-cell sequencing would show whether cells with the wrong number of chromosomes are present, and so whether the constitutional rationale applies at all; parental genotyping or long-read sequencing would settle the allelic configuration; and a functional assay would show whether the missense variant lowers BUBR1 abundance, as reported for one class of MVA alleles [33].

Beyond the patient, an independent screen containing the line C survivors would allow the one replication that is missing, and a prospective test of the ranking model on a rare disease with few publications would show whether its literature component contributes anything. Annex Table A4 lists these and further open questions with the test and data each requires.

The biases declared in Section 4.1 could themselves be measured. Repeating the searches with a second language model would show how far the retrieved literature depends on the model. Applying the selection rules to a disease with known medicines, before the result is seen, would test rules that here were written after the literature had been read. Removing the shared prior would separate line A from the ranking model, and a second blind reading of the evidence levels, scored for agreement, would show how far the grades can be trusted.

## 5 Reproducibility and data governance

### 5.1 Reproducibility

The pipeline comprises 150 Python modules with 586 collected tests; SHA-256 checksums are recorded for 184 output files, and the download record lists each source with, where a file was downloaded, its size and SHA-256; the gaps are marked in the repository.

One run on a free cloud runtime installed the locked environment unchanged, passed 169 tests, matched four downloaded sources against their recorded hashes and reproduced two outputs byte for byte; information-content scores differed in the last digit between operating systems, so written scores are rounded. The pipeline also ran unchanged on two synthetic cases, of the Fanconi anaemia and MVA types, with repeat runs producing byte-identical output.

The repository README explains how to install the pipeline and how to re-run it with access to the dataset.

### 5.2 Data governance and licences

The patient's data belong to the patient and their family; they were used under the terms of the Hackathon and are not redistributed. Patient-derived detail is held in a private folder, and only aggregates are published. This report, the supporting tables and the source code are released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0; https://creativecommons.org/licenses/by/4.0/).

Public resources used in the analysis keep their own licences, which are listed in the repository file THIRD_PARTY_LICENSES.md. No raw source file is redistributed.

### 5.3 Use of generative AI

In compliance with the Mandatory Generative AI Disclosure Guidelines, the following commercially available generative-AI tools were used:

- Provider: Anthropic. Tool: Claude Code. Models: Claude Sonnet 5 and Claude Sonnet 5.5 in the main sessions, with Claude Opus and Claude Sonnet as subagents. Plan: Claude Pro (earlier phase), Claude Max (later phase). Setting: model training on my data disabled ("Help improve Claude" off). Use: writing and verifying code; searching and reading the primary literature with a requirement that every reported sentence be found in its source by code; reading drug label text and recording a reading of each label with verbatim quotes; annotating literature sentences for direction and relevance from written definitions; auditing the repository and the reports.
- Provider: Anthropic. Tool: Claude (claude.ai). Model: Claude Opus 5.5. Plan: Claude Max. Setting: model training on my data disabled ("Help improve Claude" off). Use: editing and drafting the final report.
- Provider: Google. Tool: Gemini. Model: Gemini 3.1 Pro. Plan: Google AI Pro. Setting: activity and model training disabled. Use: generating code.

No dataset file (FASTQ, BAM or VCF) was uploaded or attached to any AI assistant, and no AI assistant analysed the data or produced the rankings. The rankings come from the deterministic Python pipeline. The subagents saw only public literature sentences, gene symbols, medicine names and the categorical case profile. The label readings and the evidence levels are the assistants' readings, and are identified as such wherever they appear.

## Acknowledgements

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible. We also thank the Wilhelm Foundation.

---

## Annex A. Supplementary tables

### A1 The refutation checks

**Table A1.** The nine refutation checks and their outcomes.

| Check | Question | Outcome |
|---|---|---|
| 1 | Is the direction of the line C association correct? | Lazertinib and neratinib pass (−0.375 and −0.459 difference in mean viability between the weak- and strong-checkpoint quartiles) |
| 2 | Are the line C survivors enriched for tubulin binders? | 4 of 6 against 0.158 expected, hypergeometric p 1.41 × 10⁻⁶: the line behaves as a general mitotic-dependency control rather than a *BUB1B*-specific one, and the four binders stay in the hazard class |
| 3 | Are the ranking model's top medicines specific to *BUB1B*? | None of the five tested is specific; the co-mention counts fall inside the distribution for degree-matched genes. The check counted co-mentions of the medicine name only, without its synonyms |
| 4 | Do short graph paths exist that avoid hub nodes? | Pomalidomide and disopyramide pass, a permissive test |
| 5 | Does the line B2 result depend on the weighting scheme? | Ceritinib, selinexor and tirbanibulin pass (p 0.0099, 0.0099 and 0.0198) |
| 6 | Is the net signed effect above 1,000 sign-shuffled networks? | Nirogacestat passes (p 0.003); four other medicines are not significant after correction for multiple comparisons (Annex Table A2), and the check's own controls did not rise above chance (p 0.054) |
| 7 | Does the patient-cell signature hold when each array is left out in turn? | None of the five medicines holds, neratinib among them |
| 8 | Is *IGBP1* close to *MTOR*, *PRKAA1* and *PRKAA2* in the network? | Mean distance 4.67 against 4.72 expected, p 0.619: the network argument for *IGBP1* is withdrawn, which does not affect the medicines themselves |
| 9 | Does the line C signal replicate in an independent screen? | Only paclitaxel could be estimated in GDSC2 [14] and its signal did not replicate (adjusted correlation −0.042, p 0.2717); the other five survivors are absent from that screen |

### A2 Entries closed and medicines kept as hazards

**Table A2.** Entries that were closed after the refutation checks, and the medicines retained only for harm experiments.

| Medicine | Line | Score in that line | Label reading 0–4 / 5–11 / 12–17 years | Safety flag | Reason it is closed or retained |
|---|---|---|---|---|---|
| Neratinib | C and E | Adjusted correlation 0.187, q 0.0316; q 0.059 in line E | N / N / N | None flagged | Check 1 passed, but the leave-one-array-out check withdrew its line E support, which closes the entry under the plan rule |
| Prednisolone | C | Correlation 0.160 | P / L / L | Muscle atrophy to verify | Did not survive adjustment (q 0.3227) |
| Cobimetinib | C | Correlation 0.129 | V / V / V | Muscle atrophy to verify | Did not survive adjustment (q 0.5892) |
| Olaparib | Ranking model | Rank 1 of 1,293 | N / N / N | None flagged | Not specific to *BUB1B* in the co-mention check |
| Tyrosine | Ranking model | Rank 2 | V / V / V | None flagged | Not specific to *BUB1B* |
| Progesterone | Ranking model | Rank 4 | V / V / V | Boxed warning | Not specific to *BUB1B* |
| Zinc chloride | Ranking model | Rank 5 | V / V / V | None flagged | Not specific to *BUB1B* |
| Vorinostat | Ranking model | Rank 7 | N / N / N | None flagged | Not specific to *BUB1B* |
| Regorafenib | D | Net 94 | N / N / N | Boxed warning | Not significant after correction for multiple comparisons (q 0.121) |
| Sunitinib | D | Net 57 | N / N / N | Boxed warning | Not significant after correction for multiple comparisons (q 0.121) |
| Ritlecitinib | D | Net 39 | V / V / V | Boxed warning | Not significant after correction for multiple comparisons (q 0.117) |
| Erdafitinib | D | Net 34 | V / V / V | None flagged | Not significant after correction for multiple comparisons (q 0.118) |
| Amphotericin B | E | Reversal score 10 | V / V / V | Boxed warning | Did not hold when arrays were left out (maximum q 0.184) |
| Auranofin | E | Reversal score 8 | N / N / N | Boxed warning | Did not hold when arrays were left out (maximum q 0.391) |
| Fluvastatin | E | Reversal score 8 | V / V / V | Muscle atrophy to verify | Did not hold when arrays were left out (maximum q 0.153) |
| Nelfinavir | E | Reversal score 8 | P / L / L | Muscle atrophy to verify | Did not hold when arrays were left out (maximum q 0.222) |
| Colchicine | C, hazard class | Adjusted correlation −0.159 | Not read | Caution-list class | Tubulin binder; retained for harm experiments only |
| Eribulin | C, hazard class | Adjusted correlation −0.115 | Not read | Caution-list class | Tubulin binder, also excluded editorially; retained for harm experiments only |
| Paclitaxel | C, hazard class | Adjusted correlation −0.141 | Not read | Caution-list class | Tubulin binder; its signal did not replicate in the independent screen |
| Vincristine | C, hazard class | Adjusted correlation −0.158 | Not read | Caution-list class | Tubulin binder; retained for harm experiments only |
| Bortezomib | Literature, hazard group | Graph rank 280 of 2,128 | Paediatric use not established | Selective killing of cells with the wrong number of chromosomes | Hazard, not proposed |
| Chloroquine | Literature, hazard group | Graph rank 102 of 2,128 | Implicit, no age stated | Preclinical genotoxicity statement | Hazard, not proposed |
| Bosutinib | Literature, hazard group, and C | Graph rank 461 of 2,128 | Labelled from 1 year, one indication only | Hazard | Hazard, not proposed; did not survive adjustment in line C (q 0.2531) |
| Calcitriol | Cross-check run 2 only | Probability of a top-five position 0.74 | Not assessed | Contraindication term matching the patient's clinical symptomatology | Excluded; an artefact of that run |

### A3 Planned laboratory work

**Table A3.** Planned experiments in the order in which they would run. None has been carried out.

| Step | Content | Condition that would contradict the hypothesis |
|---|---|---|
| First step | Karyotype and single-cell sequencing of blood, covering the megabase-scale calls; functional test of the missense variant; allelic configuration by parental or long-read testing | Stop rule S0: no excess of aneuploidy in every tissue examined, together with a missense variant that behaves as wild type or lies in cis with the stop-gain. Drug work on the patient's cells then stops, and work in cell models continues only as general research |
| Cell models | Isogenic line with partial *BUB1B* reduction in two clones, euploid controls and rescue | Failure to reproduce a checkpoint defect |
| Medicines, group 1 | Sirolimus, everolimus, metformin | A rise of 25 % or more in the aneuploid fraction or in senescence, or disagreement between the two rapalogues |
| Harm panel | Bortezomib, chloroquine, bosutinib, colchicine, eribulin, paclitaxel, vincristine | Viability no lower in the reduced-function model than in euploid controls at equal multiples |
| Medicines, group 2 | Lazertinib, and neratinib only if that entry is reopened, with proliferation matched | Loss of the association once proliferation is matched |
| Medicines, group 3 | Sodium phenylbutyrate | No effect on proteotoxic stress markers |
| Medicines, group 4 | Ceritinib, selinexor, tirbanibulin | No effect in the reduced-function model |
| Medicines, group 5 | Trametinib, in tumour lines | No effect in the tumour context |

### A4 Open questions

**Table A4.** Open questions, with the proposed test and the data each requires.

| Open question | Proposed test | Data required |
|---|---|---|
| Are cells with the wrong number of chromosomes present? | Karyotype and single-cell sequencing | A blood sample |
| Do the two *BUB1B* variants lie on different copies of the gene? | Parental genotyping, or phasing with reads spanning both variants | Parental DNA; long-read sequencing |
| Does the missense variant lower BUBR1 abundance? | Functional assay with N- and C-terminal antibodies, and transcript recovery under a decay inhibitor | Laboratory cell or protein data |
| Are the megabase-scale calls real, do they affect an MVA gene, and should *SHLD3* and *CEP63* enter the mechanism-fit gene block? | Copy-number analysis of those events | Read depth from the existing alignment |
| Does any line C survivor replicate? | The same association in a screen that contains lazertinib and neratinib | An independent sensitivity screen |
| Does the literature component of the ranking model contribute anything? | Prospective assessment on a rare disease with few publications, with the prior held fixed | A reference set of known medicines for that disease |
| Does the *IGBP1* variant modify the presentation? | Functional test of the variant | Laboratory cell data |
| Would a laboratory design for pomalidomide, disopyramide, nirogacestat and selumetinib be worth adding? | Mechanism review of each against the reduced-function model | Published mechanism data |

## References

1. Gargano MA, Matentzoglu N, Coleman B, et al.: The Human Phenotype Ontology in 2024: phenotypes around the world. Nucleic Acids Res. 2024; 52(D1): D1333–D1346. doi: 10.1093/nar/gkad1005
2. Landrum MJ, Chitipiralla S, Kaur K, et al.: ClinVar: updates to support classifications of both germline and somatic variants. Nucleic Acids Res. 2025; 53(D1): D1313–D1321. doi: 10.1093/nar/gkae1090
3. UniProt Consortium: UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res. 2025; 53(D1): D609–D617. doi: 10.1093/nar/gkae1010
4. Varadi M, Bertoni D, Magana P, et al.: AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Res. 2024; 52(D1): D368–D375. doi: 10.1093/nar/gkad1011
5. Orphanet: Orphadata, free access data from Orphanet. INSERM. https://www.orphadata.com (accessed 18 October 2026).
6. Chandak P, Huang K, Zitnik M: Building a knowledge graph to enable precision medicine. Sci Data. 2023; 10(1): 1–16. doi: 10.1038/s41597-023-01960-3
7. Szklarczyk D, Kirsch R, Koutrouli M, et al.: The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Res. 2023; 51(D1): D638–D646. doi: 10.1093/nar/gkac1000
8. Milacic M, Beavers D, Conley P, et al.: The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res. 2024; 52(D1): D672–D678. doi: 10.1093/nar/gkad1025
9. Zdrazil B, Felix E, Hunter F, et al.: The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Res. 2024; 52(D1): D1180–D1192. doi: 10.1093/nar/gkad1004
10. Türei D, Valdeolivas A, Gul L, et al.: Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. Mol Syst Biol. 2021; 17(3): e9923. doi: 10.15252/msb.20209923
11. Dempster JM, Boyle I, Vazquez F, et al.: Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. Genome Biol. 2021; 22(1): 1–23. doi: 10.1186/s13059-021-02540-7
12. Corsello SM, Nagari RT, Spangler RD, et al.: Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling. Nat Cancer. 2020; 1(2): 235–248. doi: 10.1038/s43018-019-0018-6
13. Subramanian A, Narayan R, Corsello SM, et al.: A Next Generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. Cell. 2017; 171(6): 1437–1452.e17. doi: 10.1016/j.cell.2017.10.049
14. Yang W, Soares J, Greninger P, et al.: Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Res. 2013; 41(D1): D955–D961. doi: 10.1093/nar/gks1111
15. Clough E, Barrett T, Wilhite SE, et al.: NCBI GEO: archive for gene expression and epigenomics data sets: 23-year update. Nucleic Acids Res. 2024; 52(D1): D138–D144. doi: 10.1093/nar/gkad965
16. Ferguson C, Araújo D, Faulk L, et al.: Europe PMC in 2020. Nucleic Acids Res. 2021; 49(D1): D1507–D1514. doi: 10.1093/nar/gkaa994
17. National Library of Medicine: DailyMed. https://dailymed.nlm.nih.gov (accessed 14 October 2026).
18. Food and Drug Administration: openFDA drug label application programming interface. https://open.fda.gov/apis/drug/label/ (accessed 14 October 2026).
19. Food and Drug Administration: Drugs@FDA, FDA-approved drugs. https://www.accessdata.fda.gov/scripts/cder/daf/ (accessed 18 October 2026).
20. Abou Tayoun AN, Pesaran T, DiStefano MT, et al.: Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. Hum Mutat. 2018; 39(11): 1517–1524. doi: 10.1002/humu.23626
21. Cheng J, Novati G, Pan J, et al.: Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023; 381(6664): eadg7492. doi: 10.1126/science.adg7492
22. Adzhubei IA, Schmidt S, Peshkin L, et al.: A method and server for predicting damaging missense mutations. Nat Methods. 2010; 7(4): 248–249. doi: 10.1038/nmeth0410-248
23. Ioannidis NM, Rothstein JH, Pejaver V, et al.: REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. Am J Hum Genet. 2016; 99(4): 877–885. doi: 10.1016/j.ajhg.2016.08.016
24. Schubach M, Maass T, Nazaretyan L, et al.: CADD v1.7: using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions. Nucleic Acids Res. 2024; 52(D1): D1143–D1154. doi: 10.1093/nar/gkad989
25. Ng PC, Henikoff S: SIFT: predicting amino acid changes that affect protein function. Nucleic Acids Res. 2003; 31(13): 3812–3814. doi: 10.1093/nar/gkg509
26. Pejaver V, Byrne AB, Feng BJ, et al.: Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. Am J Hum Genet. 2022; 109(12): 2163–2177. doi: 10.1016/j.ajhg.2022.10.013
27. Richards S, Aziz N, Bale S, et al.: Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 2015; 17(5): 405–424. doi: 10.1038/gim.2015.30
28. Bergquist T, Stenton SL, Nadeau EAW, et al.: Calibration of additional computational tools expands ClinGen recommendation options for variant classification with PP3/BP4 criteria. Genet Med. 2025; 27(6): 101402. doi: 10.1016/j.gim.2025.101402
29. Guney E, Menche J, Vidal M, et al.: Network-based in silico drug efficacy screening. Nat Commun. 2016; 7(1): 10331. doi: 10.1038/ncomms10331
30. Cheng F, Desai RJ, Handy DE, et al.: Network-based approach to prediction and population-based validation of in silico drug repurposing. Nat Commun. 2018; 9(1): 2691. doi: 10.1038/s41467-018-05116-5
31. Amberger JS, Bocchini CA, Scott AF, et al.: OMIM.org: leveraging knowledge across phenotype–gene relationships. Nucleic Acids Res. 2019; 47(D1): D1038–D1043. doi: 10.1093/nar/gky1151
32. Hanks S, Coleman K, Reid S, et al.: Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. Nat Genet. 2004; 36(11): 1159–1161. doi: 10.1038/ng1449
33. Suijkerbuijk SJE, van Osch MHJ, Bos FL, et al.: Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. Cancer Res. 2010; 70(12): 4891–4900. doi: 10.1158/0008-5472.CAN-09-4319
