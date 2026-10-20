"""Track 2 report laid out as a scientific article: structured abstract, introduction, methods, results, discussion, conclusions, declarations, references, supplementary material.

    python -m track2.build_article --out out --report submissions/candeleroj98_track2_report.md

The detailed material (per-drug cards, the four lines in full, the caution list from labels, the age table, the appendix of drugs that surfaced) is produced by `track2.build_report_review2.build` and placed in the
supplement with its section numbers prefixed by S; the main text is written here, with every number read from `out/`. Nothing about the evidence changes: this is a rewrite of the layout and the prose.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from track2 import brief_section as brief
from track2 import empirical_section as emp
from track2 import engine_section as eng
from track2 import followup_section as fup
from track2 import landscape_section as ls
from track2 import lines_section as lns
from track2.build_report_review2 import _load, build, f3

HEAD = re.compile(r"^## ", re.M)
WORDS = {7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven"}


def split_sections(text: str) -> dict[str, str]:
    """Old report -> {key: block}, key = leading number, 'refs' or 'appendix'."""
    parts = HEAD.split(text)
    out: dict[str, str] = {"front": parts[0]}
    for p in parts[1:]:
        title = p.split("\n", 1)[0]
        m = re.match(r"(\d+)\.", title)
        key = m.group(1) if m else ("refs" if title.startswith("References") else "appendix")
        out[key] = "## " + p
    return out


def body(block: str) -> str:
    return block.split("\n", 1)[1].strip("\n")


def cut_subsection(block: str, number: str) -> tuple[str, str]:
    """Return (subsection body without its heading, block without that subsection)."""
    m = re.search(rf"^### {re.escape(number)} .*$", block, re.M)
    if not m:
        raise SystemExit(f"subsection {number} not found")
    nxt = re.search(r"^#{2,3} ", block[m.end():], re.M)
    end = m.end() + nxt.start() if nxt else len(block)
    return block[m.end():end].strip("\n"), block[:m.start()] + block[end:]


def supplement(block: str) -> str:
    """Prefix section numbers with S, headings and in-text references."""
    block = re.sub(r"^(#{2,3}) (\d+(?:\.\d+)*)\.? ", lambda m: f"{m.group(1)} S{m.group(2)}. " if m.group(1) == "##" else f"{m.group(1)} S{m.group(2)} ", block, flags=re.M)
    block = re.sub(r"\bsections? (\d+(?:\.\d+)*)", lambda m: m.group(0).replace(m.group(1), "S" + m.group(1)), block)
    block = block.replace("sections S", "Supplement sections S").replace("section S", "Supplement section S")
    return block


def main_refs(text: str) -> str:
    """References inside text moved from the old report to the new numbering of the main text (old section -> new section or supplement)."""
    rules = [(r"sections? 4\.0(?![\d])", "section 3.3"), (r"sections? 6\.7(?![\d])", "section 4.1"), (r"sections? 1\.2(?![\d])", "section 4.6"), (r"sections? 1\.1(?![\d])", "section 4.4"),
             (r"sections? 7\.1(?![\d])", "section 3.5"), (r"sections? 7\.2(?![\d])", "Supplement S7.2"), (r"sections? (4\.\d+)", r"Supplement S\1"), (r"sections? (6\.\d+)", r"Supplement S\1"),
             (r"sections? 10(?![\d])", "section 3.6"), (r"sections? 11(?![\d])", "section 2.6"), (r"sections? 9(?![\d.])", "section 4.5"), (r"sections? 8(?![\d.])", "Supplement S8"),
             (r"sections? 7(?![\d.])", "section 3.5"), (r"sections? 5(?![\d.])", "section 2.4"), (r"sections? 4(?![\d.])", "section 3.2"), (r"sections? 3(?![\d.])", "section 3.1")]
    for pat, rep in rules:
        text = re.sub(pat, rep, text)
    return text


ABBREVIATIONS = {
    "AI": "artificial intelligence", "ANDA": "abbreviated new drug application", "ATC": "Anatomical Therapeutic Chemical classification", "ATR": "ATR serine/threonine kinase",
    "AUC": "area under the receiver operating characteristic curve", "BLA": "biologics licence application", "BRAF": "B-Raf proto-oncogene kinase", "ChEMBL": "database of bioactive molecules with drug-like properties",
    "CIN": "chromosomal instability", "CML": "chronic myeloid leukaemia", "ERK": "extracellular signal-regulated kinase", "FASTQ": "text format for sequencing reads",
    "FDA": "United States Food and Drug Administration", "FDR": "false discovery rate", "GH": "growth hormone", "HPO": "Human Phenotype Ontology", "HSCT": "haematopoietic stem cell transplantation",
    "IGF": "insulin-like growth factor", "MANE": "Matched Annotation from NCBI and EMBL-EBI", "MDM2": "MDM2 proto-oncogene, an E3 ubiquitin ligase", "MEK": "mitogen-activated protein kinase kinase",
    "mTOR": "mechanistic target of rapamycin", "mTORC1": "mechanistic target of rapamycin complex 1", "MVA": "mosaic variegated aneuploidy", "NDA": "new drug application", "PARP": "poly(ADP-ribose) polymerase",
    "PLK": "polo-like kinase", "PMC": "PubMed Central", "PMID": "PubMed identifier", "PrimeKG": "precision medicine knowledge graph", "PRISM": "Profiling Relative Inhibition Simultaneously in Mixtures, a drug-sensitivity screen of cancer cell lines",
    "RAF": "RAF kinase", "SAC": "spindle assembly checkpoint", "STRING": "database of protein-protein interactions", "TORC1": "target of rapamycin complex 1", "TTK": "TTK protein kinase (MPS1)",
    "VCF": "variant call format",
}


def abbreviations(text: str) -> str:
    """Abbreviations that occur in the text, in alphabetical order, defined once (article style)."""
    used = sorted((a for a in ABBREVIATIONS if re.search(r"(?<![A-Za-z0-9])" + re.escape(a) + r"(?![A-Za-z0-9])", text)), key=str.lower)
    return "; ".join(f"{a}, {ABBREVIATIONS[a]}" for a in used) + "."


def references(main_text: str, refs_block: str) -> str:
    """Numbered reference list: the method papers, then the papers cited by PMID in the main text in order of first citation, then the data resources."""
    meta = json.loads(Path("config/reference_metadata.json").read_text(encoding="utf-8"))["entries"]
    order = ["26831545", "30002366", "33264411"]
    for m in re.finditer(r"PMIDs? (\d{6,8})", main_text):
        if m.group(1) not in order:
            order.append(m.group(1))
    lines = []
    for n, pmid in enumerate(order, 1):
        e = meta.get(pmid)
        if not e:
            lines.append(f"{n}. PMID {pmid}.")
            continue
        who = e["first_author"] + (" et al." if e["more_authors"] else "")
        where = f"{e['journal']} {e['year']}" + (f";{e['volume']}" if e["volume"] else "") + (f":{e['pages']}" if e["pages"] else "")
        lines.append(f"{n}. {who} {e['title']}. {where}." + (f" doi:{e['doi']}" if e["doi"] else "") + f" (PMID {pmid}).")
    data = re.search(r"^\d+\. (Data resources.*)$", refs_block, re.M)
    if data:
        lines.append(f"{len(order) + 1}. " + data.group(1).replace("Human Phenotype Ontology", "Human Phenotype Ontology (reference 3)", 1))
    note = "Inline citations give the PubMed identifier (PMID); the list gives the full reference of each paper cited in the main text, in order of first citation. Papers cited only in the supplement are identified there by their PMID."
    return "## References\n\n" + note + "\n\n" + "\n".join(lines)


def compact_appendix(block: str) -> str:
    """Appendix of the drugs that surfaced: drop the ChEMBL-mechanism and reason columns (long and repeated) and say so."""
    out_lines = []
    for ln in block.split("\n"):
        if ln.startswith("|"):
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            if len(cells) == 7:
                if set("".join(cells)) <= set("-: "):
                    out_lines.append("|---|---|:-:|---|---|")
                else:
                    cells = [cells[i] for i in (0, 2, 3, 4, 5)]
                    cells[1] = re.sub(r"rank (\d+) \(outside the top 50\)", r"\1", cells[1])
                    cells[2] = cells[2].replace("B2 ", "").replace("C ", "").replace("no / no", "no / no")
                    cells[3] = cells[3].replace(" (exploratory)", "")
                    cells[4] = cells[4].replace("(network-evidence class N1)", "(N1)")
                    out_lines.append("| " + " | ".join(cells) + " |")
                continue
        out_lines.append(ln)
    text = "\n".join(out_lines)
    text = re.sub(r"<!-- widths:[\d,]+ -->", "<!-- widths:18,26,12,18,26 -->", text)
    text = text.replace("its mechanism as recorded in ChEMBL; ", "").replace("| Line A |", "| Line A rank (outside the top 50 unless stated) |").replace("| Line D |", "| Line D (exploratory) |")
    text = re.sub(r"; and the reason it is not a candidate\. The reason is derived.*?checkpoint\.", ". The ChEMBL mechanism and the reason for each row (no individual drug is supported by two usable lines of the protocol rule) were left out of this article for length; they are produced by `track2.build_report_review2.build`. The rule says nothing about whether the drug could ever help this child.", text, flags=re.S)
    return text


def compose(out: Path) -> str:
    full = build(out)
    sec = split_sections(full)
    fl, nq = _load(out / "four_lines_review2.json"), _load(out / "kg_nullquery_e5.json")
    vu, kx = _load(out / "vulnerability_e6.json"), _load(out / "kg_extra.json")["results"]["disease_only"]
    cls, props = _load(out / "class_analysis.json"), _load(out / "proposals_review2.json")
    names = {"PHENYLBUTANOIC ACID": "sodium phenylbutyrate"}
    prop = props["proposed"]
    n_drugs = len(prop)
    def pick(g: str) -> list[str]:
        return [names.get(p["drug"], p["drug"].lower()) for p in prop if p["group"] in g]
    tierA = ["sirolimus", "everolimus", "metformin", "sodium phenylbutyrate"]
    tierB = ["trametinib", "selumetinib"]
    tierC = ["bortezomib", "chloroquine", "bosutinib"]
    n_generic, n_real = nq["n_real_top50_drugs_generic_at_least_5_percent"], nq["real_top50_unique_drugs"]
    pre, exp = fl["readings"]["b2_counted_as_prespecified"], fl["readings"]["b2_held_exploratory"]
    supported = cls["supported_in_two_usable_lines"]
    vu0 = _load(out / "vulnerability.json")
    fda, prov = _load(out / "fda_approval_check.json"), _load(out / "label_provenance.json")
    n_fda, n_lab = fda["n_drugs"], prov["n_labels"]
    no_label = sorted(set(fda["drugs"]) - {k.upper() for k in prov["labels"]})
    pair_txt = lns.pairs_text(out)

    # pieces taken from the old report
    mech_body = body(sec["3"])
    w_body, sec4 = cut_subsection(sec["4"], "4.0")
    m_tab = re.search(r"^\| Drug \| Tier .*?(?=\n\n|\Z)", w_body, re.S | re.M)
    tier_table = m_tab.group(0)
    w_body = (w_body[:m_tab.start()] + w_body[m_tab.end():]).strip()
    w_body = re.sub(r"^The seven drugs.*?$", "", w_body, flags=re.M).strip()
    neg_body, sec6 = cut_subsection(sec["6"], "6.7")
    caut_lit_block = sec["7"]
    caut_lit_body, caut_rest = cut_subsection(caut_lit_block, "7.1")
    open_body, _ = cut_subsection(sec["1"], "1.2")
    limits = [ln[2:] for ln in body(sec["9"]).split("\n") if ln.startswith("- ")]

    t: list[str] = []
    t.append("# Approved drugs to investigate in mosaic variegated aneuploidy with two *BUB1B* variants: a mechanism-first analysis of one case, with a reusable pipeline\n")
    t.append("**Juan Antonio Candelero Cervantes** · Individual submission (no team) · Hugging Face user `candeleroj98` · MVA Hackathon 2026, Track 2 (drug repurposing)\n")
    t.append("*Hypotheses for follow-up and for a clinician. Approved drugs only. Not medical advice: no dose and no recommendation to treat. The age of the child is unknown.*\n")

    t.append("## Abstract\n")
    t.append(brief.answer(Path("out")) + "\n")
    t.append("**Background.** Mosaic variegated aneuploidy (MVA) is a rare disorder of chromosome segregation. We studied a child with two *BUB1B* variants (a nonsense variant and a missense variant of uncertain significance; whether they lie on different copies of the gene is unresolved), short stature, failure to thrive, skeletal muscle atrophy, nephrocalcinosis and rhabdomyosarcoma, and asked which approved drugs are plausible candidates for investigation given the disrupted mechanism.\n")
    t.append("**Methods.** We described the mechanism of reduced BubR1 function as five consequences matched to the child's features. Drug-agnostic literature searches (675 queries in two rounds, plus four sentences of contrary evidence for sirolimus added by the AI assistant after the audit, with the author's approval) were checked by code so that every quoted sentence appears in its source, and written rules selected the hypotheses. The label safety statements of each drug were read. Four computational lines defined in a written protocol (a knowledge graph, network proximity, cell-line vulnerability and a signed-network direction check) served as an independent filter, with random controls. The literature layer, a ranking engine, signature-reversal lines (knockout, combined Track 1 and patient signatures) and these four lines then each proposed candidates separately, and a multi-criteria analysis (SMAA) served as a cross-check.\n")
    t.append(f"**Results.** Nine approved drugs are listed in three tiers: four constitutional hypotheses ({', '.join(tierA)}), two tumour-context hypotheses ({', '.join(tierB)}) and three hazards to check and not propose ({', '.join(tierC)}). The evidence is level 1 at best (related systems); none reaches human *BUB1B*-deficient cells or patients, and for sirolimus the search found more findings of harm than of benefit on growth and muscle. "
             f"No drug met the rule of two usable computational lines written in the protocol, but each line, read separately, proposes its own candidates. " + lns.abstract_sentence(out) + f" The knowledge-graph ranking was generic ({n_generic} of {n_real} top drugs also ranked for random queries) and the cell-line analysis gave {vu0['observed_total']} hits (already below its null) and {vu['observed_total']} after adjusting for confounders. " + eng.abstract_sentence(Path("out/engine")) + "\n")
    t.append("**Conclusions.** The result is a set of weakly supported hypotheses, proposed line by line and not as one ranking (level 1 evidence at best; neratinib is listed by two weak lines and no drug by two strong lines), each with a computational check or an experiment that would refute it, together with a caution list, a label-by-age table and a diagnostic first step (a karyotype or single-cell sequencing). The first stages of the pipeline rerun from another individual's variant file; the hypothesis layer is specific to this case: it rests on AI-assisted literature searches, evidence labels assigned by the AI search tasks from a written definition and not reviewed by a specialist, rules applied by code, and editorial exclusions listed in Supplement S4.10.\n")

    t.append("## 1. Introduction\n")
    t.append("Mosaic variegated aneuploidy is caused by defects in genes that keep chromosome segregation faithful. The best-known cause is biallelic loss of function of *BUB1B*, which encodes BubR1, a component of the mitotic checkpoint complex of the spindle assembly checkpoint. When the checkpoint works poorly, cells enter anaphase with unattached chromosomes and a mosaic of aneuploid cells results. "
             "Affected children show growth restriction, microcephaly in many cases, a progeroid-like picture in mouse models, and a strong predisposition to childhood tumours, mainly Wilms tumour and rhabdomyosarcoma.\n")
    t.append("No approved drug restores BubR1. The question of this report is therefore not which drug corrects the gene, but which drugs that are already approved act on the biological consequences of its loss, and what would be needed to show that they help or harm. "
             "We start from the child: the Track 1 analysis of the genome gave several candidate variants, of which *BUB1B* has the strongest evidence, and the child's clinical features are the ones that the consequences of a weak checkpoint would predict, with one exception, nephrocalcinosis, for which we found no link in the literature.\n")
    t.append("We follow a mechanism-first order: mechanism, then hypotheses, then tests. Computation is used as an independent filter on the hypotheses and not as their source. We report what the filter did not confirm with the same weight as what it did, we state the strength of every piece of evidence, and we propose for each hypothesis the experiment that would refute it. These are hypotheses for follow-up. They are not evidence that any medicine works, and they must not be used to treat.\n")

    t.append("## 2. Methods\n")
    t.append("### 2.1 Case, conditions and data\n")
    t.append("The child's variants are *BUB1B* p.Leu737* (a nonsense allele, ClinVar pathogenic or likely pathogenic) and p.Asn1002Lys (a variant of uncertain significance), in MANE transcript NM_001211.6. Whether the two variants lie on different copies of the gene (compound heterozygosity) is unresolved. The phenotypes used were rhabdomyosarcoma, nephrocalcinosis, short stature, failure to thrive, skeletal muscle atrophy, premature birth and small for gestational age. "
             "Two conditions of the study apply throughout: only approved drugs are considered, and the age of the child is unknown, so label support is read only for fixed age bands and never assumed. Every drug of the study was checked in Drugs@FDA: each has an approved application (NDA, BLA or ANDA) and a product marketed as prescription or over-the-counter, so none is experimental (" + f"{fda['n_confirmed']} of {n_fda}" + " confirmed). No dose, exposure or treatment recommendation is given. The child's reads and VCF are not part of the repository, and the report names variants at protein level only.\n")
    t.append("### 2.2 Mechanism characterisation\n")
    t.append("The genes of the Track 1 submission were screened for function, disease association, mechanism and phenotype match. *BUB1B* was classified as the causal candidate, *MAD1L1* and *CENPE* as pathway context and *IGBP1* as insufficient evidence; the primary analysis uses *BUB1B* alone. Because the phase is unresolved, the combined effect is treated as two direction scenarios, complete loss of BubR1 function and a hypomorph, so that the sign of a drug's effect on the checkpoint is not assumed. "
             "Five consequences of reduced BubR1 function were matched to features of the child, each with verified sentences from the primary literature and the system in which it was measured.\n")
    t.append("### 2.3 Literature search, evidence levels and selection of hypotheses\n")
    t.append("Eight AI-assisted search tasks, in two rounds, collected sentences about approved drugs: four by consequence (growth and muscle; proteotoxic and autophagic stress; senescence and p53; tumour, chromosomal instability and safety) and four drug-agnostic (interventions in BubR1-hypomorph models; clinical management of MVA; the child's own features; vulnerabilities of aneuploid cells), 675 queries in all. "
             "Each entry carries an evidence level: L0, inferred or indirect; L1, tested in a related system (aneuploid cells, chromosomal instability, yeast, an animal with a checkpoint defect); L2, tested in *BUB1B*/BUBR1-deficient or MVA-derived human cells, or in a hypomorphic BubR1 model with rescue; L3, a clinical observation in MVA or *BUB1B*. "
             "A program re-fetched every quoted sentence from Europe PMC, and only sentences found word for word were used: 355 literature sentences in all (351 from the eight searches and 4 added by the AI assistant after the audit, with the author's approval), which is a different count from the " + str(json.loads(Path("out/quote_check.json").read_text(encoding="utf-8"))["n_quotes"]) + " label quotes checked against the label text on DailyMed or openFDA. Written rules then selected the hypotheses: the drug is approved; a constitutional hypothesis needs a verified finding of benefit at level 1 or more on a feature of the child or on a consequence of aneuploidy (consequences 2 and 3), with benefit counted from level 1 and harm only from level 3; a tumour-context hypothesis needs a verified level-1 finding in an aneuploid or tumour system; "
             "a class sibling is listed with its representative at its own level; and a level-3 harm or a flagged finding goes to the caution list. The selection is therefore exploratory: the rules were formulated after the literature had been read, the searches were not independent of the final drugs (143 of the 419 queries of the first round already named a drug that was later selected), and six drugs that satisfy the rules were set aside by judgement, each for a reason given in Supplement S4.10.\n")
    t.append("### 2.4 Computational filter\n")
    t.append("Four lines of evidence were defined in a written protocol, with validations and thresholds, before the final runs (the order rests on the author's statement; there is no external time stamp): line A, a knowledge graph (PrimeKG, personalised PageRank from the disease and the child's phenotypes; validated by leave-out of known indications over 200 diseases and against a degree-only baseline); lines B1 and B2, network proximity of drug targets to a phenotype module and to a *BUB1B* module; "
             "line C, the vulnerability of weak-checkpoint cell lines (DepMap and PRISM, adjusted for lineage, proliferation, TP53, microsatellite instability and chromosomal instability, against 1,000 matched null sets); and line D, a signed-network direction check, exploratory. A drug is a computational candidate if at least two usable lines support it. "
             "Random queries of the same shape as the real one, degree-matched nulls and class-level permutations served as controls. Some analytical choices were refined after two independent reviews by AI agents (in network-pharmacology and validation-audit roles) had seen intermediate results; they are therefore exploratory, and the variants of the analysis that were run are listed in Supplement S5. Later lines (literature layer, ranking engine, signature reversal; `protocol/addendum_multi_evidence.md`) are outside that rule; each proposes candidates with its own score (Section 3.13), and a multicriteria analysis (SMAA) is used only as a cross-check.\n")
    t.append("### 2.5 Labels and safety\n")
    t.append("For every drug in the study, the current US label was read (DailyMed, or openFDA when DailyMed has no readable text), the paediatric status was classified in six levels from label text only, and verbatim sentences were stored with the set identifier and date. For the nine drugs, the label safety statements were quoted. A checker confirms by code that every quoted sentence appears in the label it cites. "
             "The readings were made by an AI assistant; they are not a pharmacist's or a regulatory review.\n")
    t.append("### 2.6 Reproducibility\n")
    t.append("The computational part runs from the Track 1 submission file in 22 stages with fixed seeds, " + str(emp.n_tests()) + " automated tests and recorded checksums, and it was run on Windows; an earlier version was also run on Linux in a cloud notebook. Versions and file hashes of the public sources, and their licences, are recorded with the code (see Data and code availability).\n")

    t.append("## 3. Results\n")
    t.append("### 3.1 Mechanism and its consequences\n")
    t.append(lns.chain_intro(out) + "\n")
    t.append(main_refs(mech_body) + "\n")
    t.append("### 3.2 Hypotheses for investigation\n")
    t.append(f"The selection rules gave {WORDS.get(n_drugs, n_drugs)} approved drugs in three tiers (Table 1; the card of each drug, with its evidence for and against, label safety statements, paediatric label status and refuting experiment, is in Supplement S4). "
             "Tier A holds constitutional hypotheses, supported by a finding of benefit on a feature of the child or on a consequence of aneuploidy: sirolimus, with everolimus as its class sibling (mTOR inhibition, on growth, muscle and aneuploidy-related stress), metformin (oxidative stress in trisomy 21 fibroblasts), and sodium phenylbutyrate (senescence and protein aggregation in trisomic cells). "
             "Tier B holds tumour-context hypotheses: trametinib, with selumetinib as its class sibling, because highly aneuploid non-transformed cells were more sensitive to MEK inhibition. "
             "Tier C lists hazards to check and not hypotheses to treat: bortezomib, chloroquine and bosutinib kill aneuploid cells selectively, which is the reason to study them in a tumour and the reason to fear them in a child whose normal cells are aneuploid.\n")
    t.append("**Table 1.** Hypotheses by tier. Levels (L0 to L3) are defined in section 2.3.\n")
    t.append(tier_table + "\n")
    t.append("### 3.3 Strength of the evidence\n")
    t.append("The evidence for every hypothesis is weak, for five reasons. First, the search was not independent of the final drugs and the selection is exploratory (section 2.3). A second round of four searches that did not start from drugs (256 queries on BubR1-hypomorph models, clinical management of MVA, the child's own features, and vulnerabilities of aneuploid cells) added sodium phenylbutyrate and bosutinib to the list. "
             "It found no approved drug that rescues a BubR1-deficient phenotype: the level-2 findings are genetic (SIRT2 overexpression, clearance of p16-positive cells, loss of p16 or Pak2) or use compounds that are not approved.\n")
    t.append("Second, group A rests on little. For sirolimus the entries on growth and muscle (first search plus one sentence added by the AI assistant after the audit, with the author's approval; all level 0) are, by sentence, 3 suggesting benefit, 11 harm, 2 no effect and 1 unclear (by paper 3, 6, 2 and 1; two papers contribute both a harm and a no-effect sentence), and its one level-1 finding is a Drosophila model of aneuploid neuroblasts, a consequence that is not a feature of the child. "
             "For metformin the only level-1 benefit is less oxidative damage in trisomy 21 fibroblasts, against two level-1 findings that it alters checkpoint genes or is genotoxic in non-cancer cells. Sodium phenylbutyrate rests on two level-1 findings from whole-chromosome trisomy models (fibroblasts and neurons), not from mosaic aneuploidy or *BUB1B*, and its label warns about sodium retention and renal insufficiency in a child with nephrocalcinosis.\n")
    t.append("Third, the same drugs may worsen the child's own problems. Inactivation of TORC1 has been linked to mitotic slippage in yeast and human cells and to chromosomal instability in yeast (PMID 35141499), sirolimus is reported to slow growth in transplanted children (PMID 21380626), and the labels warn about renal function. Group A should be read as hypotheses with a plausible adverse mechanism, to be tested for harm as much as for benefit.\n")
    t.append("Fourth, trametinib shares the hazard of group C: its support is the higher sensitivity of non-transformed aneuploid clones, which also threatens normal aneuploid cells, so group B is a tumour-context hypothesis and not one for the constitutional state.\n")
    t.append("Fifth, \"level 1\" ranks the similarity of the system to MVA and not the strength of a benefit. A reader should look at the \"Best level\" and \"Evidence against or hazard\" columns of Table 1 together.\n")
    t.append("### 3.4 The computational filter\n")
    t.append(f"No drug met the rule written in the protocol (at least two usable lines among A, B and C), under either reading of line B ({pre['n_candidates']} candidates with B2 counted as in the protocol, {exp['n_candidates']} with B2 held exploratory). Read one by one, the lines do propose candidates (Section 3.13), and two drugs are listed by two lines, {pair_txt}, each pair including a weak line. "
             f"The knowledge-graph line is valid but modest: over 200 diseases its mean AUC is {f3(kx['mean_auc'])}, and its excess over a degree-only baseline is {f3(kx['mean_excess_over_degree'])}. Its ranking for this child is generic: {n_generic} of the {n_real} distinct drugs of its top 50 also appear in the top 50 of at least 5 % of 500 degree-matched random queries, and fostamatinib, the single candidate of the first analysis, is in the top 5 of {nq['fostamatinib']['random_top5']:.1%} of them, which exposes it as a connectivity artefact. "
             f"The phenotype module (B1) does not beat the degree baseline once the genes causal for the validation diseases are removed. The *BUB1B* module (B2) passes the validation threshold written in the protocol but cannot be tested for circularity. The cell-line analysis (line C) gave {vu0['observed_total']} hits with lineage and proliferation as covariates, already below its null 95th percentile ({vu0['null']['total_p95']:.2f}); with TP53, microsatellite instability and chromosomal instability added it gave {vu['observed_total']} hits against a null 95th percentile of {vu['null']['total_p95']:.1f}; six of the early hits were MDM2 inhibitors, the signature of TP53 status. "
             "The signed-network line is exploratory and assigns no tier. One class of multi-target kinase inhibitors (ATC L01EX) passes in two lines (network evidence N1, below level 1; tumour context); the result is compatible with chance and none of its members is among the *BUB1B*-module drugs. The detail is in Supplement S5 and S6.\n")
    t.append("### 3.5 Caution\n")
    t.append("**Table 2.** Cautions from the verified literature.\n")
    t.append(main_refs(caut_lit_body) + "\n")
    t.append("The labels of five drug classes (tubulin binders; mitotic-checkpoint, Aurora, PLK and TTK inhibitors; genotoxic chemotherapy; proteasome inhibitors; growth hormone) were also searched for warnings relevant to this child, and each warning is quoted verbatim in Supplement S7.2. The growth hormone label, for example, contraindicates active malignancy, which matters because the child's record documents a rhabdomyosarcoma, although its current status (active or not) is not available. The labels of selumetinib, bortezomib and chloroquine also report chromosome-level findings in preclinical tests (aneugenic activity in the mouse micronucleus assay at high doses, clastogenicity in vitro, and a potential genotoxicity risk in some test systems, respectively). In a syndrome defined by chromosome missegregation these findings are a reason for caution and not evidence of harm in this child; the quotes are in Supplement S4.\n")
    t.append("The vinca and taxane findings are in vitro and from case reports. No clinical data reviewed here show that a child with MVA should avoid standard rhabdomyosarcoma chemotherapy, and every treatment decision, including surveillance for Wilms tumour, rhabdomyosarcoma and leukaemia, belongs to the child's oncology and genetics team; this report makes no statement on intervals or regimens.\n")
    t.append("### 3.6 Labels by age group\n")
    t.append("Because the age of the child is unknown, label support was tabulated for three groups (0 to under 3, 3 to under 12, and 12 years and older) from label text only. A drug is called labelled for a group only when its label says so, and labels that refer to paediatric dosing without a stated age (for example chloroquine) are read as 'verify', not as established use. The table, with the count of drugs per group, is in Supplement S10; the " + str(n_lab) + " label readings, with set identifiers, dates and verbatim quotes, are in the repository (`out/label_provenance.json`" + (f"; the Drugs@FDA check of Section 2.1 covers {n_fda} drugs, and {', '.join(x.lower() for x in no_label)} has no label reading" if no_label else "") + "). The per-line tables of Section 3.13 use the age bands 0-4, 5-11 and 12-17 of the multi-evidence outputs, which are not the groups of this table. Label status describes the populations studied for other diseases. It is not a statement of suitability for this child, and the drugs of S4 are proposed for in-vitro investigation, not for this child.\n")

    t.append(emp.main_section(Path("out/empirical")))
    t.append(ls.main_sections(Path("out")))
    t.append(eng.main_section(Path("out/engine")))
    t.append(lns.main_sections(out))
    t.append(lns.cross_check(out))
    t.append(fup.main_section(out))
    t.append("## 4. Discussion\n")
    t.append("### 4.1 What the negative computational result means\n")
    t.append(main_refs(neg_body).replace(", and every amendment was stricter, not looser", "") + f" Reading the later lines one by one adds a narrower finding: {pair_txt} are the only drugs that two lines list, and each pair includes a weak line (Section 3.13).\n")
    t.append("### 4.2 What this could contribute to understanding and diagnosing MVA\n")
    t.append("If an experiment on *BUB1B*-deficient human cells confirms or refutes one of these hypotheses, the child's clinicians and other families gain either a shortlist with a rationale or a documented dead end, and both are useful. Three contributions do not depend on any drug being right. "
             "First, a diagnostic step: the blood genome of the child shows only a weak aneuploidy signal that cannot be read as mosaicism, so a karyotype or single-cell sequencing is the first test, and the constitutional hypotheses depend on its result. "
             "Second, a caution list for clinicians, built from published cases (a severe reaction to dactinomycin in a child with biallelic *BUB1B* variants and Wilms tumour, graft rejection after a treosulfan-based conditioning), from evidence that a weak checkpoint reduces sensitivity to taxanes, and from the labels. "
             "Third, a table of paediatric label support that can be used whatever the age. The negative results are also information for others: they show where network and cell-line methods lack power for a gene with no targeted drug.\n")
    t.append("### 4.3 What is new in the approach\n")
    t.append("Computation is used to discard as well as to propose. Random queries of the same shape as the real one exposed the first candidate as a connectivity artefact; class-level analysis replaced hub drugs; the unknown age is handled by label bands; every quoted sentence is checked by code against its source, and the exploratory steps and the drugs set aside are written down so that the process can be audited.\n")
    t.append("### 4.4 Applicability to other individuals, diseases and larger datasets\n")
    t.append("The first stages of the pipeline (screening the genes of a Track 1 file for mechanism, building gene modules and network proximity) run from the submission file of any individual and are the reusable part. Three parts are written for this case and would have to be redone for another: the knowledge-graph query and the cell-line analysis (lines A and C), which use a disease node and a checkpoint score chosen for MVA; the five consequences of the mechanism; and the selection of the nine drugs, which applies written rules (R1 to R4, run by code) to evidence labels assigned by the AI search tasks from a written definition and not reviewed by a specialist, on AI-assisted literature searches whose sentences were verified against their sources by code. The rules were written after the literature had been read, and six drugs that satisfy them were set aside by editorial decision (Supplement S4.10). A test with a different submission file (Fanconi anaemia and ATR genes) leaves the report text unchanged. We therefore describe the approach as a method to repeat, not as a program that takes a file and returns hypotheses. It uses public resources only (ChEMBL, STRING, PrimeKG, DepMap, PRISM, the Human Phenotype Ontology, DailyMed), its knowledge-graph line was validated across 200 diseases, and its results are reproducible: fixed seeds, " + str(emp.n_tests()) + " automated tests, recorded checksums and a run of an earlier version in a cloud notebook. The labels must be read for each drug. " + emp.scalability_text(Path("out/empirical")) + " " + ls.scalability_sentence(Path("out")) + " " + lns.engine_scalability(out) + " " + lns.signature_reversal_scalability(out) + " The pipeline has been run on one case, so no claim is made about its performance in general.\n")
    t.append("### 4.5 Limitations\n")
    t.append("This is a single case, so no sensitivity, specificity or clinical significance is claimed. The age of the child is unknown, and the phase and the functional effect of the two variants are unresolved. The current tumour status of the child is not available, so the cautions that depend on it (growth hormone, tumour-context hypotheses) cannot be assessed. The evidence for every hypothesis is indirect (section 3.3), and the hypotheses were selected after the literature had been read. The literature layer was checked with a negative control only through a deterministic proxy based on mentions (Section 3.7); the AI-assisted procedure itself, which also named drugs in its first-round queries, was not tested on drugs without a plausible link to the mechanism. "
             "The computational lines have limited power for a gene without a targeted drug: network and graph methods give no direction of effect, the validations use diseases with known drugs, the *BUB1B*-module line could not be tested for circularity, cell-line evidence is a tumour proxy that does not represent the constitutional picture, and only 1,040 of the 2,705 approved drugs can be scored in the network. "
             "The label readings and the evidence levels were produced with an AI assistant and checked by code only for the existence of each sentence in its source, not for what it is taken to mean; they are not a pharmacist's or a regulatory review and must be re-checked against the current label. Finally, the analysis plan was not registered with an external registry or independently time-stamped.\n")
    t.append("### 4.6 Open questions\n")
    t.append(main_refs(open_body) + "\n")
    t.append(lns.criteria(out))

    t.append("## 5. Conclusions\n")
    t.append("A mechanism-first reading of this case of two *BUB1B* variants gives nine approved drugs to consider, four of them as constitutional hypotheses, two for the tumour context and three only as hazards, with evidence that is weak, partly against and never from human cells with this defect. "
             "The computational filter did not confirm any of them by its own rule and has little power to do so for a gene without a targeted drug; its silence neither supports nor rules out a hypothesis. Read line by line, the other lines propose their own candidates, with neratinib in two weak lines and no drug in two strong ones, and the multi-parameter cross-check is a cross-check and not a result (Sections 3.13 and 3.14). "
             "What the work leaves is a set of refutable experiments, a caution list, a table of label support by age, a diagnostic first step, and a pipeline whose first stages can be rerun for another case. These are hypotheses for follow-up and not evidence that any medicine works.\n")

    t.append("## Declarations\n")
    t.append("### Data protection\n")
    t.append("The child's reads and variant calls are not part of the repository. Per-variant detail is kept in a private folder and is not published, and this article names variants at protein level only. The Track 1 submission file, which contains the variant coordinates required by the challenge format and is the single case-specific input of the pipeline, is not included in the repository. All data derived from the challenge dataset will be deleted by 23 November 2026, as the hackathon rules require.\n")
    t.append("### Data and code availability\n")
    t.append("The code, configuration and derived outputs accompany the submission in a public repository (link in the submission form; MIT licence). The public sources used, with versions, file hashes and licences, are listed in the repository. The challenge dataset is gated and is not redistributed.\n")
    t.append("### Generative-AI disclosure\n")
    t.append(body(sec["13"]) + "\n")
    t.append("### Competing interests and funding\n")
    t.append("The author declares no competing interests. The challenge prizes were sponsored by AWS and Anthropic, and the author used Anthropic's Claude Code as described above. No specific funding; computing and AI-assistant subscriptions were paid by the author.\n")
    t.append("### Acknowledgments\n")
    t.append(body(sec["14"]) + "\n")
    t.append("### Abbreviations\n")
    t.append(abbreviations("\n".join(t)) + "\n")
    t.append(references("\n".join(t), sec["refs"]) + "\n")

    t.append("## Supplementary material\n")
    t.append("The section numbers of the supplement are those of the full report with the prefix S: S4.1 is the card of sirolimus, S6.5 the combination of the lines, S10 the table of labels by age group. Section S4.0 and S6.7 of the full report are in the main text as 3.3 and 4.1.\n")
    for key in ("4", "5", "6", "7", "8", "10"):
        blk = {"4": sec4, "6": sec6, "7": caut_rest}.get(key, sec[key])
        t.append(supplement(blk).strip() + "\n")
    t.append(ls.supplement(Path("out")))
    t.append(emp.supplement(Path("out/empirical")))
    t.append(eng.supplement(Path("out/engine")))
    t.append(fup.supplement(out))
    t.append(compact_appendix(supplement(sec["appendix"])).replace("## Appendix A.", "## S-Appendix A.", 1).strip() + "\n")
    return "\n".join(t)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--report", type=Path, default=Path("submissions/candeleroj98_track2_report.md"))
    args = parser.parse_args()
    text = compose(args.out)
    # one word for the three groups of drugs, in the main text and in the supplement
    for pat, rep in ((r"\btiers\b", "groups"), (r"\bTiers\b", "Groups"), (r"\btier\b", "group"), (r"\bTier\b", "Group")):
        text = re.sub(pat, rep, text)
    args.report.write_text(text, encoding="utf-8", newline="\n")
    print("wrote", args.report, len(text.split()), "words")


if __name__ == "__main__":
    main()
