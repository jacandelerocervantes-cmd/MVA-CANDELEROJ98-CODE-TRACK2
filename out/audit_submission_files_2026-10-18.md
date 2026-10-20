# Audit of the Track 2 submission files (read-only), 2026-10-18

Files audited: `submissions/candeleroj98_track2_report.md` (1,201 lines, about 27.7k words of prose and tables), `submissions/candeleroj98_track2_report.pdf` (47 pages, created 2026-10-18 22:50 UTC, two minutes after the md), `submissions/candeleroj98_methods_description_track2.md` and `submissions/candeleroj98_methods_description_form_track2.xlsx` (both last written 2026-10-15 10:48, i.e. before the 2026-10-18 restructure of the report). Nothing was edited, regenerated or deleted. Line numbers refer to the report md ("R:") and the methods md ("M:"). The xlsx was read cell by cell (both sheets); the PDF through text extraction and page geometry (no raster viewer was available, so visual rendering was not inspected).

Method summary: regex scans of all four files (coordinates, HGVS, ids, paths, emails, dates, dosing words, age words, overclaim words); cell-by-cell comparison of xlsx answers with the methods md; heading and 789-sentence comparison of md against PDF text; reference/section/supplement resolution scripts; 40+ typed numbers compared with JSON outputs in `out/`.

## 1. Result in one paragraph

No privacy or standing-rule violation was found (no genomic coordinates, no cDNA, no rsIDs, no sample ids, no local paths, no emails, no transcript-deletion plan, no assumed age, no dose, only one deletion date, 20 December 2026, in the report). Three High findings concern staleness and disclosure: the CF1 to CF9 follow-ups are described as "defined and not run" although their results were written after the report was built; the methods description and xlsx predate the report restructure; and the AI disclosure omits a model named in the report and mis-states the phase dates. The PDF is the same version as the md (all 56 headings and 789 sampled sentences found; 47 pages).

Findings: 3 High, 8 Medium, 13 Low (24 in total).

## 2. Prioritised findings

| # | Sev | File and location | Quote or value | Why it matters | Suggested fix |
|---|---|---|---|---|---|
| F1 | High | Report R:9 (abstract), R:378-412 (Sec. 3.15), R:447-452 (Sec. 4.7), PDF same; `out/computational_followups.json/.md` | R: "defined and not run: CF1 to CF9"; "Defined and not run (computational follow-ups)". But `out/computational_followups.json` and `.md` (written 2026-10-18 17:47, after the report 16:48 and the PDF 16:50) hold results for CF1 to CF9: CF1 Lazertinib/Neratinib "survives"; CF2 tubulin class "refuted" (hypergeometric p 1.41e-06); CF3 Olaparib, Progesterone, Tyrosine, Vorinostat, Zinc Chloride "refuted"; CF4 Disopyramide/Pomalidomide survive; CF9 not runnable. | The report states a fact (not run) that no longer holds, and omits results that close five ENG candidates and change the "16 to check by computation first, 10 for a laboratory design, 4 hazards, 2 closed" counts (R:449-452). `out/checksums.sha256` (16:50) also predates the follow-ups. Rigor criterion. | Decide whether the follow-ups are part of the submission; if yes, regenerate the report from the JSON (Sec. 3.15, abstract, conclusions, 4.7 counts, checksums); if no, move the files out of `out/` or label them. |
| F2 | High | Methods md M:15-26, M:60-63 and xlsx cells B8, B11, B17 (identical text) | Describes "four computational lines" as the whole filter; no mention of the ranking engine (R Sec. 3.12), of per-line proposals LIT/ENG/EP/E/ET1 (Sec. 3.13), neratinib, the SMAA cross-check (3.14), the case profile (3.11), the differential diagnosis (3.8), CF1-CF9 or the 0-4/5-11/12-17 bands. M:32: "The proposals of section 4 of the report" (the proposals are in Sec. 3.2/3.13 and S4). M:60: "No drug is supported by two lines" without the per-line outcome of the report. | The form is scored with the report (rigor, innovation); a reader of the form gets a different, older design than the 47-page report. The innovation claim (engine, per-line proposals) is invisible in the form. | Regenerate the methods md and xlsx from the final report (same generator, `fill_methods_form.py`) after the report is frozen; fix "section 4". |
| F3 | High | AI disclosure: R:472-479 = M:Section 3 = xlsx B9 (identical) | Models listed: "Claude Sonnet 5 and Claude Sonnet 5.5" and "Gemini 3.1 Pro". Report R:239 says "Two language models (Opus and Sonnet) ... labelled sentences" (files `out/engine/annotations_opus.json`, `annotations_sonnet.json`). Disclosure phase: "Final phase (2026-10-06 to 2026-10-08)"; report and CF work run to 2026-10-18. "executed locally and on private cloud infrastructure" while the text also says Google Colab. "All variant calling, annotation, scoring and the final CSV were produced by ... `track2.run_pipeline` and `track2.build_candidates`" (Track 2 does not call variants; sentence copied from Track 1). | The AI disclosure is the one required item of the form; an unlisted model, a closed date window and a mis-stated infrastructure are the kind of inconsistency a reviewer can quote. | Add every model and plan used (including the Opus annotator and the sub-agents' models), extend the phase dates to cover 2026-10-09 to 2026-10-18, replace "private cloud" by "Google Colab CPU runtime (once)", and remove the variant-calling sentence. Plan and tier statements could not be verified from the files. |
| F4 | Medium | R Table 1 (R:~90-99), R:115 (Sec. 3.4), R:842 (S6.3), R:1127 (S-Appendix A) versus R:9, R:15, R:279-364 (Sec. 3.13), R:362 | Table 1 gives "C no" for trametinib; 3.4 and S6.3: "Line C is not usable and supports no drug"; S-Appendix A: "among the line C hits (none is)". Sec. 3.13 and abstract: five line-C candidates, "trametinib (LIT and C)" and "neratinib (C and E-P)" are in two lines. | Line C is simultaneously "not usable" and a source of proposals; Table 1 and 3.13 disagree on trametinib. The abstract's caveat (each pair includes a weak line) helps but the three places are not reconciled. | State once how line C is used in 3.13 (descriptive, not usable under the protocol) and align Table 1 and Appendix A wording. |
| F5 | Medium | R:470 (Data and code availability); all four files | R: "public repository (link in the submission form; MIT licence)". No "CC BY 4.0" for the submission, no embargo statement (peer-reviewed manuscripts using the dataset until organisers publish a report or preprint) in any of the four files; methods md/xlsx give no licence and no repository link. README.md states CC BY 4.0 for the submission and MIT for the code. | Standing rules require CC BY 4.0 and the embargo statement; the report names only MIT, so the report and README disagree. | Add to Declarations: submission released under CC BY 4.0, code under MIT, and one sentence on the embargo; add the repository URL (the report points to the form, the form has no link). |
| F6 | Medium | M:38 / xlsx B12 versus R:466 | M/xlsx: "all data will be deleted 30 days after the hackathon ends". R: "deleted by 20 December 2026". | Same obligation, two wordings; the form has no date, and `docs/archive_addendum10/README_addendum10.md` still records a different close date. No wrong date appears in the four files. | Use "by 20 December 2026" in the form text. |
| F7 | Medium | R:140, R:281, S10 (R:904-950); M:Q2 item 10 | Two age schemes: 3.6 and S10 use 0 to <3, 3 to <12, 12 and older; 3.13 uses 0-4, 5-11, 12-17 (acknowledged at R:140). Methods md/xlsx describe only the first (plus five finer bands). The abstract says "age bands" without saying which. Age of the child is correctly never assumed. | Two sets of bands in one report invite a reviewer to compare tables that cannot be compared. | Choose one scheme for the main text; mention the second only in the annex or the methods. |
| F8 | Medium | R:858 (S6/S7 text), R:533, R:535-1125 headings | "below level 1 on the scale of Supplement section S2.3": no S2.3 exists. Supplement headings are S4, S5, S6 (6.1-6.6), S7 (only S7.2, no S7.1), S8, S10, S13, S12, S14, S15, S-Appendix A: S1-S3, S9, S11 are absent and S13 precedes S12. R:533 says "Section S4.0 and S6.7 of the full report are in the main text as 3.3 and 4.1" (no "full report" exists). Main text cites S4, S4.1-S4.10, S5, S6.x, S7, S7.2, S8, S10, S12.x, S13, S14, S15: all resolve except S2.3. All "Section n.n" references resolve to existing headings. | Dangling and gappy numbering after the restructure. | Renumber supplements consecutively or add a note; fix S2.3 and R:533. |
| F9 | Medium | Abstract R:9-15; Abbreviations R:491 | Used in the abstract before any definition and absent from the abbreviation list: LIT, ENG, E-P, E, B1, C, CF1 to CF9, SMAA; also W3, R1-R4, X1-X5, E0, REVEL, SIFT, CADD, MONDO, CASMI, NMD, N1. List has 36 entries. | Abstract is not self-contained; undefined codes lower readability for the panel. | Spell out in the abstract (literature layer, ranking engine, ...) and complete the list. |
| F10 | Medium | References R:495-529 | Text says the list is "in order of first citation", but the order of first PMID citation in the main text differs from the list. PMIDs cited in the main text but not in the list: 38987356, 41374403. In the list but not cited by PMID in the main text: 26831545, 30002366 (cited by name), 33264411. 31 entries, 30 with a PMID. Inline citations are PMIDs, not reference numbers. | Reference list does not match the citations; a numbered list with unnumbered citations cannot be cross-checked. | Add the two missing references; use author-year or numbered citations; fix the sentence about order. |
| F11 | Medium | PDF/md structure | 47 pages; the main text with references ends on page 23 (about 22.5 pages), supplements pp. 24-47. | Above the planned 20-page limit (not restructured here; see section 3). | Planning item, see length table. |
| F12 | Low | R:449 (Sec. 4.7) versus R:216-232, R:430 | "the engine was validated on held-out diseases"; Sec. 3.12 says the embeddings model is not supported, ranking stays exploratory, and the low-volume advantage comes from the graph prior. Other uses of "validated" (leave-out validation of lines A and B2) are technical and fine. No "proven", "cure", "safe", "effective" and no "treats" in a claim; "efficacy" appears only in negations; "recommend" only in "no recommendation" or label quotes. | Single wording stronger than the evidence tier. | Replace by "tested on held-out diseases". |
| F13 | Low | R:51 and R:430 | "566 automated tests". `pytest --collect-only` now collects 577 (545 `def test_` lines); README.md says 325/169/171. | Typed number built at build time; stale after new tests. | Regenerate at freeze; update README. |
| F14 | Low | R:15, R:115, R:842 versus R:334, R:382, R:449 | "22 hits" (`out/vulnerability.json` observed_total 22 = 19 PRISM rows + 3 CRISPR genes) versus "6 of 23 original hits survive" (`out/line_c_adjusted.json` n_original_hits 23; n_survive_all_drugs 7). Both are real JSON values for different counting units, but the report never says why the original line-C set is 22 in one place and 23 in another. | Reader sees two totals for one line. | Add one sentence defining 22 (rows/genes) versus 23 (drugs). |
| F15 | Low | R:150, R:39, M:Q5 / xlsx B11 | "363 sentence slots proposed by the AI search tasks, 355 found, 6 not checkable, 2 not found": `verification.json` has 363 rows and 355 verified, but these include the 4 sentences added by the AI assistant (R:39 says 351 + 4 = 355). Methods Q5 says "eight AI-assisted tasks ... a second round of four more searches", which reads as 12 searches; the report says eight in two rounds. Methods gives "143 of its 419 queries" while the report total is 675 (419 is the first round only; 419 + 256 = 675 checked). | Counts are right but wording is ambiguous. | Rephrase "proposed by the AI search tasks" and clarify rounds in the form. |
| F16 | Low | R:434, R:1127, R:307, Table 1 | "only 1,040 of the 2,705 approved drugs can be scored in the network" while line A ranks "of 2128" and S-Appendix A says "946 of the 2,705". 1,040 is the line-B scorable count (`benchmark_results.json`), 2,128 the line-A count. | Easy to read as a contradiction. | Name the line for each denominator. |
| F17 | Low | PDF | 0 link annotations (no clickable PMIDs or cross-references); about 27 % of characters are 7.8 pt (tables) against 10 pt body; PDF title metadata carries literal asterisks ("*BUB1B*") and creator "HeadlessChrome". No text lies outside the page, no replacement characters, no missing glyphs (characters such as the angle quotes, Angstrom, delta are present). Visual check of tables not possible here. | Small tables in a 47-page PDF may be hard to read; cosmetics. | Check pages with wide tables on screen; clean title metadata. |
| F18 | Low | R:9, R:62, R:178, R:152, R:10 occurrences; M/xlsx | Privacy information for judgement: protein-level changes `p.Leu737*` (R 2, M 2, xlsx 2) and `p.Asn1002Lys` (R 7, M 2, xlsx 2), transcript `NM_001211.6` (R 2, M 1, xlsx 1), "exon 17 of 23", "746 nt upstream of the last exon-exon junction", REVEL 0.472, CADD-PHRED 24.5, residue-environment statistics for position 1002; `p.Leu1012Pro` is a published pathogenic variant, not the child's. Absent: chrN:pos, g./c. HGVS, rsIDs, GRCh/hg builds, sample or proband ids (the word "proband's genome" appears once in M/xlsx). | The positional detail is the only item beyond "protein-level without positions"; the user must judge whether exon number and nucleotide distance are acceptable. | Decide; if not, drop "exon 17 of 23" and "746 nt". |
| F19 | Low | R:60 | "premature birth and small for gestational age"; "parental recurrent miscarriage, describes the parents and is not part of any phenotype query. No birth details are published." | Family and perinatal detail; consistent with the report's own phrasing but worth a conscious decision. | Keep or remove the parental miscarriage sentence. |
| F20 | Low | M:Q5 / xlsx B11 | "OMIM numbers and PMIDs from that review are re-validated before publication" (future tense, no outcome); "this reading was done by the author's assistant" repeated. | Open promise inside a final form. | State the outcome or delete. |
| F21 | Low | Report tables | Only "Table 1" and "Table 2" are numbered; about 20 further tables (3.11, 3.12 test table, 3.13 per-line tables, S-sections) have no number or caption. Naming varies: "Phenylbutanoic Acid" (R Sec. 3.13 tables, from ChEMBL) versus "Sodium phenylbutyrate" (Sec. 3.2). | An article needs numbered, captioned tables; the name variation looks like two drugs. | Number tables when restructuring; use one display name. |
| F22 | Low | xlsx | Sheet "Track 1 methods" left unfilled (template); file metadata: created by openpyxl, lastModifiedBy "Verena C" (organiser template), modified 2026-08-27. Cell comments are the template's ("Yellow cells are for your answers"). Answers B7-B17 equal the methods md text (B8, B9, B10, B11, B12, B13, B14, B15, B16 identical; B7 differs by one dash; B17 differs only by the md note "(Abstract length: 499 words.)"). Abstract in B17 has 499 words, within the stated limit of 500 (the only stated word limit). Question order and count (11) match; md headings paraphrase xlsx question text for Q3, Q6, Q8, Q9, Q11. | No defect; reported for completeness. | Optional: clear the organiser name if unwanted. |
| F23 | Low | Cross-file | Consistent across files: author "Juan Antonio Candelero Cervantes", Hugging Face `candeleroj98`, individual submission, ChEMBL_37, 2,705 / 1,040, 200 diseases, AUC 0.802, excess 0.022, 43 of 46, line C 5 hits against 17.0, 551 quotes, 355/2/6 sentences, 22 stages (report only). Not stated in the methods: test count, stage count, licence, embargo. | Confirms the numeric core is coherent; gaps are in F2 and F5. | None. |
| F24 | Low | R:123-129 | Ages "11-year-old", "9-year-old", "38-month-old" appear only for published cases (PMIDs 31053147, 35434947, 31081598), not for the child. The dose regex gave one hit ("474 mg"), whose context was not located (probably inside a quoted sentence); the user should grep it once. Paediatric dosing is mentioned only to say "verify". | No violation. | None. |

Items checked and found correct: all 11 form questions answered in order; no emails; no `C:\`, `/Users/`, scratchpad or AppData strings; no "transcript-deletion plan"; 20 December 2026 is the only calendar deletion date; no "25 December"; dosing only in negations ("No dose, exposure or treatment recommendation is given", "NOT dosing advice"); Track 1 referred to as the input (closed) and nothing in the four files re-opens or re-scores it, but the Track 1 report itself was not compared (not part of this audit).

## 3. Numbers checked against `out/` JSON (all matched unless noted)

| Typed value (report) | Source | Result |
|---|---|---|
| 675 queries (419 first round, 256 second) | `literature_hypotheses/*.json` queries (125+73+75+146 = 419; 38+48+53+117 = 256) | match |
| 363 slots, 355 verified, 6 not checkable, 2 not found | `literature_hypotheses/verification.json` | match (see F15) |
| 551 label quotes, 0 not found | `quote_check.json` | match |
| 2,705 approved drugs; 1,040 scorable; 2,128 line A | manifests, `benchmark_results.json`, `class_analysis.json` | match |
| 200 diseases, mean AUC 0.802, CI 0.782 to 0.821 | `kg_validation.json` | match |
| Excess over degree baseline 0.022 | 0.802 minus 0.780 baseline (text); no stored key | consistent |
| 43 of 46 top-50 drugs generic; fostamatinib 98.8 % | `kg_nullquery_e5.json` | match |
| Line C: 22 hits, null p95 23.05; 5 hits, p95 17.0, p 0.16, 1,183 lines, VIF 4.94, correlation 0.84 | `vulnerability.json`, `vulnerability_e6.json` | match (see F14) |
| 6 of 23 line-C hits survive | `line_c_adjusted.json` | match |
| Lines: 6 LIT, 10 ENG, 2 A, 3 B2, 5 C, 5 D, 5 EP; 28 distinct drugs; 2 in two lines; pool 224; 1,293 native ENG; 38 cytotoxic; 28 unconfirmed approval | `candidates_by_line_mva_v2.json` | match |
| 30 % label corruption: 8.2 of nine drugs keep group | `empirical/label_robustness.json` (8.151) | match |
| Literature proxy AUC 0.75 (0.68 to 0.84), 26 diseases | `empirical/literature_control.json` (0.754; 0.675 to 0.840) | match |
| Kappa 0.84, 0.56, 0.62, 0.85 | `out/engine/*` (0.844, 0.565, 0.625, 0.851) | match |
| Engine gains +0.21 (17 diseases), +0.19 (42 held-out) | table R:~222-226 | not re-derived (only stated in table; consistent with abstract and 4.4) |
| 566 tests | `pytest --collect-only` 577 | stale (F13) |
| 22 stages | `track2.run_all --list` | match |
| 20 December 2026 | `docs/data_protection.md` | match |

Statements about what was run versus defined: CF1 to CF9 are stated as not run but have outputs (F1); the adjustment of line C is stated as done and matches its output; "the early stages were run twice on two synthetic cases" refers to `out_case_*` folders that exist (not re-verified).

## 4. PDF versus markdown

Same version: all 56 headings of the md found in order in the PDF, all 789 sampled md sentences/table cells (longer than 60 characters, alphanumeric match of first and last characters) found in the PDF, the final sentence present. Words: md 27,749 (without table pipes) versus PDF 28,053 extracted tokens (differences are table-layout splits). Pages: 47. No text outside the 612 x 792 pt page; no replacement characters. See F17 for links, font size and metadata.

## 5. Length analysis for the planned article (under 20 pages)

Pages are PDF pages (letter size, 10 pt body, tables 7.8 pt). Start page is where the heading begins; a section ends on the start page of the next one. Words are from the md.

| Top-level section | PDF pages (start to end) | Approx. pages | Words | Candidate |
|---|---|---|---|---|
| Abstract | 1-2 | 1.3 | 909 | main (long: four sub-paragraphs, 909 words) |
| 1. Introduction | 2 | 0.3 | 267 | main |
| 2. Methods | 2-4 | 1.7 | 970 | main |
| 3. Results | 4-18 | 14.0 | 9,133 | main (parts) and annex (parts), see below |
| 4. Discussion | 18-20 | 2.0 | 1,760 | main |
| 5. Conclusions | 20-21 | 0.4 | 172 | main |
| Declarations (data protection, code, AI, interests, acknowledgments, abbreviations) | 21-22 | 1.3 | 816 | main (AI disclosure 372 words) |
| References | 22-23 | 1.3 | 707 | main |
| Supplementary material intro | 23-24 | 0.1 | 55 | annex |
| S4 Proposed candidates (S4.1-S4.10) | 24-33 | 9.0 | 5,957 | annex |
| S5 Design | 33-34 | 0.6 | 380 | annex |
| S6 Results | 34-36 | 1.7 | 1,498 | annex |
| S7 Caution | 36-37 | 0.6 | 356 | annex |
| S8 What would need to be true | 37-38 | 0.8 | 530 | annex |
| S10 Age: label support by group | 38-39 | 0.9 | 563 | annex (3.6 summarises) |
| S13 Gene landscape | 39-41 | 1.3 | 787 | annex |
| S12 Empirical checks | 41-43 | 1.4 | 834 | annex |
| S14 Evidence engine: protocols | 43 | 0.4 | 270 | annex |
| S15 Follow-up experiments | 43-45 | 1.6 | 807 | annex or main (experiments table) |
| S-Appendix A (71 drugs) | 45-47 | 2.2 | 888 | annex |
| Total | 47 | 47 | 27,749 | |

Main body, front matter through references: pages 1-23 (about 22.5 pages, about 14,700 words including references). Supplements: pages 24-47 (about 24 pages, about 12,900 words). To be under 20 pages the main part needs to shrink by roughly 3 pages (about 15 %) with references and declarations included.

Results subsections (3):

| Subsection | Start page | Words | Content type | Candidate |
|---|---|---|---|---|
| 3.1 Mechanism and its consequences | 4 | 1,074 | evidence table | main |
| 3.2 Hypotheses for investigation (Table 1) | 5 | 544 | core result table | main |
| 3.3 Strength of the evidence | 6 | 400 | table | main |
| 3.4 The computational filter | 7 | 303 | core result | main |
| 3.5 Caution | 7 | 764 | caution tables | main summary, detail annex (S7) |
| 3.6 Labels by age group | 8 | 176 | summary | main (detail S10) |
| 3.7 Empirical checks of the method | 8 | 902 | method validation | annex candidate (S12 duplicates it) |
| 3.8 MVA genes and differential diagnosis | 9 | 343 | context | annex candidate (S13 duplicates it) |
| 3.9 Mechanism of the two variants | 10 | 358 | case-specific | main (short) |
| 3.10 Phase of the two variants | 11 | 232 | case-specific | main |
| 3.11 Minimal case profile and comparison | 11 | 723 | table | main or annex |
| 3.12 Ranking engine on other rare diseases | 12 | 1,340 | secondary method (tables) | annex candidate (summary sentence in main) |
| 3.13 Candidates proposed by each line | 14 | 1,132 | per-line tables | main (compressed) or annex |
| 3.14 Multi-parameter cross-check | 16 | 234 | cross-check | annex candidate |
| 3.15 Refutable by computation or experiment | 17 | 606 | table | main (compressed) |

Methods subsections (2.1 to 2.6): 161, 109, 326, 219, 90, 63 words; Discussion subsections: 4.1 287, 4.2 179, 4.3 77, 4.4 551, 4.5 254, 4.6 243, 4.7 167 (4.7 is a restatement of the criteria and could be shortened).
