# Computational follow-ups CF1 to CF9

Label: exploratory, post hoc; cannot prove efficacy; a candidate that survives stays a hypothesis; a refuted one is closed as a computational candidate. Rules fixed before running: `protocol/addendum_computational_followups.md`. Every number is computed by `src/track2/computational_followups.py`; no patient data.

| CF | Candidates concerned | Test | Result | Verdict | Reading |
|---|---|---|---|---|---|
| CF1 | Lazertinib, Neratinib | positive rho: viability against the checkpoint score; orientation-explicit quartile delta | Lazertinib: delta -0.375, adjusted rho 0.14; Neratinib: delta -0.459, adjusted rho 0.187; orientation sanity median LFC Bortezomib -4.17, Colchicine -1.8, Paclitaxel -2.23, Vincristine -1.81 | survives (resolved by documentation) | rho > 0 means lower viability in weak-checkpoint lines, i.e. sensitivity; Lazertinib and Neratinib are not dropped on sign grounds |
| CF2 | Colchicine, Eribulin, Paclitaxel, Vincristine | hypergeometric test among estimable drugs | 4 of 6 survivors are tubulin-class; 5 of 190 estimable drugs are; expected 0.158; hypergeometric p 1.41e-06 | refuted | enriched: the score behaves as a mitotic-dependency control and the four stay hazard-class |
| CF3 | Olaparib, Progesterone, Tyrosine, Vorinostat, Zinc Chloride | BUB1B co-mention count against 99 degree-matched random genes (Europe PMC) | Olaparib 0 vs null median 0 (p 1); Progesterone 4 vs null median 0 (p 0.17); Tyrosine 6 vs null median 1 (p 0.23); Vorinostat 1 vs null median 0 (p 0.14); Zinc Chloride 0 vs null median 0 (p 1) | Olaparib refuted, Progesterone refuted, Tyrosine refuted, Vorinostat refuted, Zinc Chloride refuted | title/abstract co-mention of the drug with BUB1B is not above 99 degree-matched genes; only the co-mention leg is tested and aliases such as BubR1 are not searched |
| CF4 | Pomalidomide, Disopyramide | hub-free paths of length <= 4 from drug to the query nodes in PrimeKG | Disopyramide: hub-free query nodes reached 7; Pomalidomide: hub-free query nodes reached 7 | Disopyramide survives, Pomalidomide survives | hub-free paths of at most 4 edges exist, so line A is not refuted; this is a permissive check (most drugs with targets have such paths) and edge sources are not recorded per edge in PrimeKG |
| CF5 | Ceritinib, Selinexor, Tirbanibulin | top set under W1/W2 from stored scores; rank against 100 random degree-matched modules | Ceritinib: top set W1 True, W2 True, random-module p 0.0099; Selinexor: top set W1 True, W2 True, random-module p 0.0099; Tirbanibulin: top set W1 True, W2 True, random-module p 0.0198 | Ceritinib survives, Selinexor survives, Tirbanibulin survives | support is not specific to scheme W3 if the drug is in the top set under W1 or W2 and above random modules |
| CF6 | Regorafenib, Sunitinib, Nirogacestat, Ritlecitinib, Erdafitinib | net signed effect against 1,000 sign-shuffled networks; control re-check | Erdafitinib: net 34, null range [-10, 39], p 0.0709; Nirogacestat: net 49, null range [-29, 41], p 0.003; Regorafenib: net 94, null range [6.98, 118], p 0.121; Ritlecitinib: net 39, null range [-10, 41], p 0.047; Sunitinib: net 57, null range [-5.02, 69], p 0.105; controls recovered 4/6 (null p 0.0539) | Erdafitinib refuted, Nirogacestat survives, Regorafenib refuted, Ritlecitinib refuted, Sunitinib refuted | signed-path nets are compared with sign-shuffled networks of the same topology; controls not above chance |
| CF7 | Neratinib, Amphotericin B, Auranofin, Fluvastatin, Nelfinavir | six leave-outs of the 2 case vs 4 control arrays; concordance with BUB1B/Track 1 signatures | Amphotericin B: max q 0.184, failed leave-outs 3, concordant False; Auranofin: max q 0.391, failed leave-outs 5, concordant False; Fluvastatin: max q 0.153, failed leave-outs 4, concordant False; Nelfinavir: max q 0.222, failed leave-outs 4, concordant False; Neratinib: max q 0.153, failed leave-outs 3, concordant False | Amphotericin B refuted, Auranofin refuted, Fluvastatin refuted, Nelfinavir refuted, Neratinib refuted | reversal must hold with q < 0.10 in all six leave-outs and the drug must not mimic the BUB1B or Track 1 signature |
| CF8 | Sirolimus, Everolimus, Metformin | STRING distance IGBP1 to MTOR/PRKAA1/PRKAA2 against degree-matched random pairs | mean distance IGBP1 to MTOR/PRKAA1/PRKAA2 4.67 vs null mean 4.72; p 0.619 | refuted | the IGBP1 rationale is dropped from X1; the drugs are not affected |
| CF9 | Lazertinib, Neratinib, Colchicine, Eribulin, Paclitaxel, Vincristine | second independent screen of the same drugs on disk | each line C survivor has exactly one row in the PRISM compound list (one screen), and no GDSC, CTRP or PRISM secondary file is on disk; the DepMap CRISPR matrix is a different assay | not runnable | needs: an independent drug-sensitivity screen covering the same drugs and cell lines (GDSC, CTRP or PRISM secondary/MTS data), which would have to be downloaded |

## Updated status of the triaged candidates

| Candidate | Line(s) | Checks | Status | Next step |
|---|---|---|---|---|
| Sirolimus | LIT | CF8 refuted | unchanged (no per-drug computational check applies to this row); IGBP1 rationale dropped (CF8, affects the rationale only) | laboratory X1 |
| Everolimus | LIT | CF8 refuted | unchanged (no per-drug computational check applies to this row); IGBP1 rationale dropped (CF8, affects the rationale only) | laboratory X1 |
| Metformin | LIT | CF8 refuted | unchanged (no per-drug computational check applies to this row); IGBP1 rationale dropped (CF8, affects the rationale only) | laboratory X1 |
| Phenylbutanoic Acid | LIT | none | unchanged (no per-drug computational check applies to this row) | laboratory X2 |
| Trametinib | LIT, C | none | unchanged (no per-drug computational check applies to this row) | laboratory X3 |
| Lazertinib | C | CF1 survives; CF9 not runnable | still a hypothesis; passed CF1 (CF9 not runnable) | laboratory X4 |
| Neratinib | C, EP | CF1 survives; CF7 refuted; CF9 not runnable | closed as a computational candidate (refuted by CF7); note: CF1 passed, and the refuted check concerns only the other line's support (plan rule: row closed; reopening is a user decision) | closed or flagged hazard (CF7 rule) |
| Prednisolone | C | none | closed (earlier adjusted analysis; not part of the CF set) | none |
| Cobimetinib | C | none | closed (earlier adjusted analysis; not part of the CF set) | none |
| Ceritinib | B2 | CF5 survives | still a hypothesis; passed CF5 | laboratory X5 |
| Selinexor | B2 | CF5 survives | still a hypothesis; passed CF5 | laboratory X5 |
| Tirbanibulin | B2 | CF5 survives | still a hypothesis; passed CF5 | laboratory X5 |
| Pomalidomide | A | CF4 survives | still a hypothesis; passed CF4 | needs a laboratory design (none in the plan yet) |
| Disopyramide | A | CF4 survives | still a hypothesis; passed CF4 | needs a laboratory design (none in the plan yet) |
| Olaparib | ENG | CF3 refuted | closed as a computational candidate (refuted by CF3) | closed |
| Tyrosine | ENG | CF3 refuted | closed as a computational candidate (refuted by CF3) | closed |
| Progesterone | ENG | CF3 refuted | closed as a computational candidate (refuted by CF3) | closed |
| Zinc Chloride | ENG | CF3 refuted | closed as a computational candidate (refuted by CF3) | closed |
| Vorinostat | ENG | CF3 refuted | closed as a computational candidate (refuted by CF3) | closed |
| Regorafenib | D | CF6 refuted | closed as a computational candidate (refuted by CF6) | closed |
| Sunitinib | D | CF6 refuted | closed as a computational candidate (refuted by CF6) | closed |
| Nirogacestat | D | CF6 survives | still a hypothesis; passed CF6 | needs a laboratory design (none in the plan yet) |
| Ritlecitinib | D | CF6 refuted | closed as a computational candidate (refuted by CF6) | closed |
| Erdafitinib | D | CF6 refuted | closed as a computational candidate (refuted by CF6) | closed |
| Amphotericin B | EP | CF7 refuted | closed as a computational candidate (refuted by CF7) | closed or flagged hazard (CF7 rule) |
| Auranofin | EP | CF7 refuted | closed as a computational candidate (refuted by CF7) | closed or flagged hazard (CF7 rule) |
| Fluvastatin | EP | CF7 refuted | closed as a computational candidate (refuted by CF7) | closed or flagged hazard (CF7 rule) |
| Nelfinavir | EP | CF7 refuted | closed as a computational candidate (refuted by CF7) | closed or flagged hazard (CF7 rule) |
| Colchicine | C | CF2 refuted; CF9 not runnable | hazard class; CF2 refuted | laboratory H1 (toxicity readouts only); never a benefit candidate |
| Eribulin | C | CF2 refuted; CF9 not runnable | hazard class; CF2 refuted | laboratory H1 (toxicity readouts only); never a benefit candidate |
| Paclitaxel | C | CF2 refuted; CF9 not runnable | hazard class; CF2 refuted | laboratory H1 (toxicity readouts only); never a benefit candidate |
| Vincristine | C | CF2 refuted; CF9 not runnable | hazard class; CF2 refuted | laboratory H1 (toxicity readouts only); never a benefit candidate |

Not concluded from any row: clinical benefit, safety or dose. A survivor is a hypothesis that was not refuted by a post hoc computational check.

## Not runnable

* CF9: each line C survivor has exactly one row in the PRISM compound list (one screen), and no GDSC, CTRP or PRISM secondary file is on disk; the DepMap CRISPR matrix is a different assay. Needed: an independent drug-sensitivity screen covering the same drugs and cell lines (GDSC, CTRP or PRISM secondary/MTS data), which would have to be downloaded
