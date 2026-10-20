"""Sections of the Track 2 report that go from mechanism to hypothesis: the mechanism of the variants, the proposed candidates and the literature-based caution list.

Used by `track2.build_report_review2`. Every quoted sentence is read from the verified literature entries (`out/literature_hypotheses/*.json`, re-checked by `track2.verify_literature`); a source that is not
verified stops the build. The prose around the quotes is written here; the figures (levels, ranks, label status) come from `out/proposals_review2.json`, `out/four_lines_review2.json` and the label review.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

SAFETY_PATH = Path("config/proposal_safety.json")
TIER = {"PHENYLBUTANOIC ACID": "A. constitutional", "SELINEXOR": "B. tumour context", "BOSUTINIB": "C. hazard to check, not proposed", "SIROLIMUS": "A. constitutional", "EVEROLIMUS": "A. constitutional (class sibling)", "METFORMIN": "A. constitutional",
        "TRAMETINIB": "B. tumour context", "SELUMETINIB": "B. tumour context (class sibling)",
        "BORTEZOMIB": "C. hazard to check, not proposed", "CHLOROQUINE": "C. hazard to check, not proposed"}
EXCLUDED_WHY = {"RUXOLITINIB": "its level-1 finding is clinical improvement in Down syndrome with a STAT1 gain-of-function, an interferon signature that is not shown in BUB1B mosaic aneuploidy and is not a feature of the child",
                "CORTICOTROPIN": "the only observation is infantile spasms refractory to ACTH in an MVA boy (level 3); the child has no such feature and there is no evidence that it helps",
                "HYDROCHLOROTHIAZIDE": "halted nephrocalcinosis progression in children with X-linked hypophosphatemia (level 0); no link between nephrocalcinosis and aneuploidy or MVA was found, so it is on-label supportive care at most",
                "CHLOROTHIAZIDE": "did not reverse established nephrocalcinosis in rats (level 0)",
                "SODIUM CITRATE": "did not prevent nephrocalcinosis in preterm neonates (level 0)",
                "TEMSIROLIMUS": "no benefit in rhabdomyosarcoma in a paediatric phase 3 trial (level 0)",
                "HYDROXYCHLOROQUINE": "only level-0 findings",
                "VENETOCLAX": "only level-0 findings"}
REASON = {"MECASERMIN": "IGF-1 axis with the same tumour concern as growth hormone; the evidence is single cases in other syndromes, one with IGF-1 resistance",
          "NIACINAMIDE": "no approved systemic drug label (sold as a vitamin or supplement); one level-1 Drosophila finding, and a mouse study in which nicotinamide supplementation conferred resistance to chemotherapy in tumours",
          "FUROSEMIDE": "the strongest risk factor for nephrocalcinosis in preterm infants, and the child was born prematurely",
          "TEMOZOLOMIDE": "improved survival in relapsed rhabdomyosarcoma but raised grade 3 or higher toxicity; genotoxic in a genome-instability syndrome",
"CARBOPLATIN": "platinum genotoxicity is a safety concern in a genome-instability syndrome; the only finding is an association of chromosomal instability with carboplatin sensitivity in ovarian cancer",
          "DACTINOMYCIN": "a severe adverse reaction in a child with biallelic BUB1B mutations (single case)",
          "DASATINIB ANHYDROUS": "with quercetin as a senolytic it may promote tumorigenesis in vitro; listed with the senolytic caution",
          "ETOPOSIDE": "induces aneuploidy in culture and is genotoxic; no benefit shown in a checkpoint-deficient system",
          "OLAPARIB": "no study shows PARP-inhibitor synthetic lethality with BUBR1 loss, and aneuploid human cells were more resistant",
          "PACLITAXEL": "low BubR1 or a weak checkpoint is associated with taxane resistance, in human cells with germline BUB1B variants",
          "SOMATROPIN": "growth hormone: the label contraindicates active malignancy and the only clinical observation in MVA is an unsatisfactory response (one girl with MVA type 2)",
          "TREOSULFAN": "graft rejection after a treosulfan-based conditioning in a child with MVA1 (single case; the paper does not attribute it to the drug)",
          "VINCRISTINE": "cells with a weak checkpoint slip out of mitosis, which supports reduced sensitivity to vinca alkaloids"}
WHO = {"constitutional": "A laboratory with cell culture and cytogenetics (chromosome spreads or single-cell sequencing). Material: fibroblasts or myoblasts from the child, with consent, or a BUB1B-hypomorphic line, and a euploid control. The study is a standard cell-culture study of weeks; cost and time are not estimated here.",
       "tumour": "A cancer-biology laboratory. Material: rhabdomyosarcoma lines from a public collection scored for aneuploidy and checkpoint strength (no patient material needed), and the child's tumour sample if one exists and the family consents.",
       "hazard": "A cell-biology laboratory with the child's fibroblasts (consented) and a euploid control; only a toxicity comparison is needed."}

TOPICS = ["growth_muscle", "proteotoxic_autophagy", "senescence_p53", "tumour_cin_safety", "context", "round2_models", "round2_clinical", "round2_phenotypes", "round2_vulnerabilities", "round3_audit"]
PRIORITY = ["SIROLIMUS", "EVEROLIMUS", "METFORMIN", "PHENYLBUTANOIC ACID", "TRAMETINIB", "SELUMETINIB", "SELINEXOR", "BORTEZOMIB", "CHLOROQUINE", "BOSUTINIB"]
DISPLAY = {"PHENYLBUTANOIC ACID": "Sodium phenylbutyrate"}


def disp(n: str) -> str:
    return DISPLAY.get(n, n.title())

LEVELS_TEXT = ("**Evidence levels used for the proposals.** L0: inferred or indirect, no experiment with the drug in a relevant system. L1: the drug tested in a related system (aneuploid cells, cells with chromosomal instability, "
               "yeast, an animal with a spindle-checkpoint defect) with a measured effect. L2: tested in human cells with a BUB1B/BUBR1 defect or derived from MVA patients, or in a hypomorphic BubR1 model with rescue of the phenotype. "
               "L3: a clinical observation in patients with MVA or BUB1B. **No approved drug reaches L2 or L3 as a candidate in the literature that was searched.**")

# English description of the system of each source that is quoted below (the research notes mixed languages; these are written once, here).
SYSTEM = {
    "16411201": "cells from patients with premature chromatid separation / MVA (monoallelic BUB1B mutations), 2006",
    "23222101": "budding yeast, 2012",
    "41820377": "Drosophila larval neural stem cells depleted of spindle-checkpoint genes, a model of MVA-type complex aneuploidy, 2026",
    "18516091": "BubR1-insufficient mice, 2008",
    "23602569": "BubR1 progeroid mice lacking p53 or p21, 2013",
    "31738183": "mice carrying MVA-related BubR1 alleles; a correlation, no drug, 2020",
    "40333415": "a 2025 review of congenital aneuploidy syndromes (secondary source, citing the primary literature)",
    "35804254": "a 2022 systematic review of 31 MVA1 cases; ascertainment bias likely",
    "23300461": "BubR1 mutant mice, 2012",
    "16182441": "an embryonal rhabdomyosarcoma from an MVA case with biallelic BUB1B mutations, 2006",
    "31081598": "a 38-month-old child with premature chromatid separation / MVA and biallelic BUB1B mutations, 2019",
    "31312666": "muscle-derived progenitor cells of progeroid Zmpste24-deficient mice, 2019",
    "33484480": "DNA-repair-deficient Ercc1 progeroid mice, 2021",
    "31697050": "mice with genetic inactivation of mTORC1 in muscle, 2020 (no drug)",
    "21380626": "children with a kidney transplant treated with sirolimus against matched controls, 2011; not MVA",
    "23865768": "an introductory statement of a paediatric kidney-transplant study, 2013",
    "31308131": "aged rats with sarcopenia, low-dose rapalog (RAD001), 2019",
    "40313040": "children and young adults with relapsed solid tumours, everolimus with lenvatinib, 2025",
    "40552127": "human trisomy 21 fibroblast lines, 2025",
    "35568294": "human cancer cell lines at 1 to 16 mM, 2022",
    "23900241": "healthy middle-aged male mice fed metformin, 2013",
    "39251587": "non-transformed human RPE1 clones with chemically induced aneuploidy, 2024",
    "39247952": "human cell lines made aneuploid with reversine, and RPE1 clones, 2024",
    "21315436": "trisomic mouse embryonic fibroblasts and tumour cell lines, 2011",
    "26404941": "trisomic mouse fibroblasts (a later paper restating the same finding, so not an independent replication), 2015",
    "34215848": "human constitutional trisomy cell lines, 2021",
    "26144554": "human cells with BubR1 knocked down and a CENP-E inhibitor; MG132 is a research proteasome inhibitor and not an approved drug, 2015",
    "31053147": "an 11-year-old boy with MVA1 (BUB1B) and myelodysplastic syndrome after a transplant, 2019",
    "39264246": "a 2024 review of childhood cancer risk in genomic-instability disorders including MVA",
    "33209710": "a 2020 review of syndromic Wilms tumour (the text on MVA)",
    "17234259": "a paclitaxel-resistant human ovarian carcinoma cell line, 2007",
    "38987356": "the human pancreatic cancer line PANC-1 (a CRISPR screen), 2024",
    "39014450": "human cells and tumour material carrying germline BUB1B variants, 2024",
    "20663956": "mouse models of aneuploidy with a defective spindle checkpoint, 2010",
    "38612842": "cancer-associated fibroblasts in culture, 2024",
    "35434947": "a 9-year-old girl with MVA type 2 (CEP57, not BUB1B), 2022",
}
CURATED = {
    "SIROLIMUS": {"for": [("41820377", "boosting autophagy", 330), ("41820377", "A similar result", 160), ("31738183", "sarcopenia", 240), ("31312666", None, 240), ("31697050", "maintenance of muscle size", 240)],
                  "against": [("41820377", "negative impact", 240), ("41820377", "none of the genetic interventions", 240), ("41820377", "might have a negative impact", 260), ("31697050", "developing and growing muscle", 240), ("33484480", None, 240), ("23865768", "Concerns have been raised", 260),
                              ("21380626", "adverse effect on growth", 440), ("35141499", "chromosome instability", 240), ("35141499", "In human cells", 300)]},
    "EVEROLIMUS": {"for": [("31308131", None, 240), ("40313040", None, 200)], "against": []},
    "METFORMIN": {"for": [("40552127", None, 260)], "against": [("35568294", None, 200), ("23900241", "higher dose", 260)]},
    "TRAMETINIB": {"for": [("39251587", "trametinib", 330)], "against": []},
    "SELUMETINIB": {"for": [("39251587", "selumetinib", 330)], "against": []},
    "BORTEZOMIB": {"for": [("39247952", "proteasome inhibition", 240), ("39247952", "BJ-hTERT", 400)], "against": [("26144554", "MG132", 200)]},
    "PHENYLBUTANOIC ACID": {"for": [("31356639", "premature senescence", 260), ("32820178", "protein aggregates", 300)], "against": []},
    "SELINEXOR": {"for": [("38589567", "KPT-330", 260)], "against": []},
    "BOSUTINIB": {"for": [("31980556", "Furthermore, cells with an ongoing CIN", 330)], "against": []},
    "CHLOROQUINE": {"for": [("21315436", None, 260), ("26404941", None, 240)], "against": [("34215848", "autophagic flux", 260)]},
}
LABEL_SHORT = {"PHENYLBUTANOIC ACID": "urea cycle disorders including neonatal onset (indication); tablets not recommended up to 20 kg", "SELINEXOR": "not established (adult myeloma indications)", "BOSUTINIB": "stated from 1 year (chronic-phase Ph+ CML only)", "SIROLIMUS": "stated from 13 years (renal-transplant prophylaxis); a repackager's label", "EVEROLIMUS": "stated from 1 year (Afinitor, a brain tumour in tuberous sclerosis); not established under 18 years (Zortress)",
               "METFORMIN": "stated from 10 years (type 2 diabetes); a repackager's label", "TRAMETINIB": "stated from 1 year, only with dabrafenib (BRAF V600E tumours); not as a single agent",
               "SELUMETINIB": "stated from 1 year (neurofibromatosis type 1, plexiform neurofibromas)", "BORTEZOMIB": "not established (adult indications); trial data in children",
               "CHLOROQUINE": "implicit (paediatric dosing for malaria); no age stated"}


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


class Literature:
    """Verified literature entries, by PMID."""

    def __init__(self, directory: Path):
        status = {(r["topic"], r["section"], r["index"]): r["status"] for r in _load(directory / "verification.json")["rows"]}
        self.by_pmid: dict[str, list[dict[str, Any]]] = {}
        for topic in TOPICS:
            data = _load(directory / f"{topic}.json")
            for section in ("findings", "background_evidence", "against_or_risk"):
                for i, x in enumerate(data.get(section, [])):
                    if status.get((topic, section, i)) == "verified" and x.get("pmid") and x.get("sentence"):
                        self.by_pmid.setdefault(str(x["pmid"]), []).append({**x, "topic": topic, "section": section})

    def sentence(self, pmid: str, contains: str | None = None, maxlen: int = 230) -> dict[str, Any]:
        rows = self.by_pmid.get(pmid)
        if not rows:
            raise SystemExit(f"PMID {pmid} has no verified entry: the report cannot cite it")
        pick = next((r for r in rows if contains and contains.lower() in r["sentence"].lower()), rows[0])
        text = pick["sentence"].strip()
        if len(text) > maxlen:
            cut = text[:maxlen].rsplit(" ", 1)[0]
            text = cut + " [...]"
        return {"text": text, "level": int(pick.get("evidence_level") or 0), "pmid": pmid, "system": SYSTEM.get(pmid) or str(pick.get("system") or "see out/literature_hypotheses")}


def bullet(ev: dict[str, Any]) -> str:
    return f"- «{ev['text']}» (PMID {ev['pmid']}; {ev['system']}; level {ev['level']})"


def mechanism_section(lit: Literature) -> str:
    ev = lit.sentence
    rows = [
        ("1", "Chromosome missegregation and mosaic aneuploidy", "the defining feature of MVA; the child's genome shows only a weak aneuploidy signal, which is not a karyotype", [ev("16411201")]),
        ("2", "Proteotoxic and energy stress in aneuploid cells", "not a feature of the child; the way aneuploid cells are damaged, and a handle for drugs", [ev("23222101", maxlen=140), ev("41820377", "proteostasis failure", maxlen=220)]),
        ("3", "p53 and p16 senescence, with progeroid features", "growth restriction and muscle loss as ageing-like features", [ev("18516091", maxlen=170), ev("23602569", maxlen=200)]),
        ("4", "Growth restriction and muscle atrophy", "short stature, failure to thrive, small for gestational age, skeletal muscle atrophy", [ev("31738183", "sarcopenia", maxlen=240), ev("31738183", "L1002P", maxlen=260)]),
        ("5", "Tumour predisposition", "rhabdomyosarcoma", [ev("40333415", "high risk of developing", maxlen=200), ev("35804254", "Wilms tumor, rhabdomyosarcoma", maxlen=170), ev("16182441", maxlen=200)]),
    ]
    t = ["## 3. The mechanism of the variants\n"]
    t.append("| Variant (MANE transcript NM_001211.6) | Working interpretation | Status |\n|---|---|---|")
    t.append("| *BUB1B* p.Leu737* | Loss of function (premature stop; ClinVar pathogenic / likely pathogenic) | LoF allele |")
    t.append("| *BUB1B* p.Asn1002Lys | Working hypothesis: partial loss of function; variant of uncertain significance | Effect unproven; phase of the two alleles unresolved |\n")
    t.append("BubR1, encoded by *BUB1B*, is a component of the mitotic checkpoint complex of the spindle assembly checkpoint (SAC), which keeps the anaphase-promoting complex inactive until every kinetochore is attached. "
             "No approved drug restores it. Because the phase of the two variants is unresolved, the combined effect is taken as two **direction scenarios**: complete loss of BubR1 function, and a hypomorph (reduced function). The sign of a drug's effect on the checkpoint is therefore not assumed anywhere in this report. "
             "The mechanism screen of the Track 1 genes classified *BUB1B* as the causal candidate; *MAD1L1* and *CENPE* (non-coding structural calls) as pathway context only, and *IGBP1* (variant of uncertain significance) as insufficient evidence, so the primary analysis uses *BUB1B* alone.\n")
    t.append("Weakened checkpoint function leads to five consequences that matter for drug repurposing. Each is tied below to a feature of the child and to what the primary literature shows. Every sentence is quoted from the source and was re-checked against its text by code; "
             "the system in which it was measured is stated because a result in yeast, mice or tumour lines is not a result in this child.\n")
    t.append("| # | Consequence | Feature of the child it would explain | What is known (verified sentence, source, level) |")
    t.append("|:-:|---|---|---|")
    for n, cons, feat, evs in rows:
        t.append(f"| {n} | {cons} | {feat} | " + "<br>".join(f"«{e['text']}» (PMID {e['pmid']}; {e['system']}; level {e['level']})" for e in evs) + " |")
    t.append("\nNephrocalcinosis and the parental miscarriages are not established consequences of *BUB1B*, and no publication that links nephrocalcinosis with MVA or *BUB1B* was found in the search; they are treated as part of the phenotype, not as explained by the variants. "
             "The evidence of consequences 2 to 5 comes almost entirely from mice, yeast and cell lines; the cases in patients (consequence 5) are single reports.\n")
    return "\n".join(t)


def _label_line(label: dict[str, Any] | None) -> str:
    if not label:
        return "Label not read in this version."
    age = label.get("min_age_years")
    age_txt = "no minimum age stated" if age is None else ("minimum age 1 month" if age < 1 and age > 0 else f"minimum age {age:g} years".replace(" 1 years", " 1 year"))
    verdict = {"established": "paediatric use stated", "not_established": "paediatric use not established", "verify": "unresolved"}.get(label.get("verdict"), label.get("verdict"))
    return f"{verdict} ({age_txt}); {label.get('scope')}."


def candidates_section(lit: Literature, proposals: dict[str, Any], four: dict[str, Any]) -> str:
    by_name = {p["drug"]: p for p in proposals["proposed"]}
    order = [n for n in PRIORITY if n in by_name]
    t = ["## 4. Proposed candidates: approved drugs to investigate, anchored to the mechanism\n"]
    t.append("This is the answer to the question of the track: which approved drugs are plausible candidates for further investigation, given the mechanisms that are disrupted. They are **hypotheses for investigation and for a clinician, not evidence that any drug works and not a treatment recommendation**, "
             "and the age of the child is unknown. Each proposal names the consequence it targets, the evidence for it and **against** it, the paediatric label status, what the computational filters say, and the experiment that would refute it.\n")
    t.append("**How the drugs were chosen.** The rules were fixed after the four computational lines had been run, and are therefore declared as post hoc (amendment E16, `protocol/addendum_e16.md`): "
             "(R1) the drug is approved; (R2) a verified finding of benefit at level 1 or more on a feature of the child or on a consequence of aneuploidy (consequences 2 and 3), with no verified level-3 harm, for a *constitutional* hypothesis, or a verified level-1 finding in an aneuploid or tumour system for a *tumour-context* hypothesis; "
             "(R3) a drug of the same class as a candidate is listed with it and keeps its own lower level; (R4) a verified harm at level 3 or a finding flagged as caution goes to the caution list (section 7), while for a proposed drug a level 0 or 1 harm is shown as evidence against. "
             f"The literature search had eight tasks in two rounds (four by consequence group, then four drug-agnostic) and {proposals['n_drugs_with_verified_entries']} approved drugs with at least one verified entry. " + LEVELS_TEXT + "\n")
    t.append("The nine drugs fall into three tiers. **A, constitutional hypotheses** (sirolimus, with everolimus as its class sibling, metformin and sodium phenylbutyrate): a verified finding of benefit on a feature of the child. **B, tumour-context hypotheses** (trametinib, with selumetinib as its class sibling): relevant to the child's tumour predisposition. "
             "**C, hazards to check, not proposed** (bortezomib, chloroquine and bosutinib): they kill aneuploid cells selectively, which is the property that makes them interesting for a tumour and dangerous for a child whose normal cells are aneuploid. Tier C is listed so that the hazard is explicit and the refuting experiment is a toxicity check.\n")
    t.append("### 4.0 How much weight the evidence carries\n")
    t.append("The proposals are weak. (1) **The search was not independent of the drugs:** 143 of its 419 queries name a drug that later became a proposal, and the selection rules were written after the findings had been read (`protocol/addendum_e16.md`). "
             "(1b) **A second search round** (four drug-agnostic searches on BubR1-hypomorph models, clinical management of MVA, the child's own features and vulnerabilities of aneuploid cells; 256 queries in all) added sodium phenylbutyrate, selinexor and bosutinib; six drugs that pass the mechanical rules were excluded editorially, each with its reason in section §EXCL§ (the six include selinexor, which is one computational line plus literature and so falls under E9). It found no approved drug that rescues a BubR1-deficient phenotype: the level-2 findings are genetic (SIRT2, clearance of p16-positive cells, p16 or Pak2 loss) or use compounds that are not approved. (2) **Tier A rests on very little, and sodium phenylbutyrate on less:** its two level-1 findings come from whole-chromosome trisomy models (fibroblasts and neurons), not from mosaic aneuploidy or BUB1B, and its label warns about sodium retention and renal insufficiency in a child with nephrocalcinosis. (2b) For sirolimus the entries of the growth-and-muscle search (all verified, all level 0) are, by sentence, 3 that suggest benefit, 10 harm, 2 no effect and 1 unclear; by paper, 3, 5, 2 and 1 (two papers contribute both a harm and a no-effect sentence). The count is made by `out/literature_hypotheses/growth_muscle.json`, drug field, direction field; its one level-1 finding is a Drosophila model of aneuploid neuroblasts, a consequence that is not a feature of the child. "
             "For metformin the only level-1 benefit is less oxidative damage in trisomy-21 fibroblasts, against two level-1 findings that it changes checkpoint genes or is genotoxic in non-cancer cells. "
             "(3) **The same drugs may worsen the child's own problems:** mTOR inhibition has been linked to mitotic slippage and chromosomal instability (PMID 35141499), sirolimus is reported to slow growth in transplanted children (PMID 21380626), and the labels warn about renal function, the organ with nephrocalcinosis in this child. Tier A should be read as hypotheses with a plausible adverse mechanism, to be tested for harm as much as for benefit. "
             "(4) **Trametinib has the same hazard as tier C:** its support is the higher sensitivity of non-transformed aneuploid clones, which also threatens normal aneuploid cells; tier B is tumour context only and not for the constitutional state.\n")
    t.append("| Drug | Tier | Consequence targeted | Best level | Evidence against or hazard | Paediatric label | Computational filters |")
    t.append("|---|---|---|:-:|---|---|---|")
    short_cons = {"SIROLIMUS": "4 (growth, muscle) and 2", "EVEROLIMUS": "4 (growth, muscle); class sibling of sirolimus", "METFORMIN": "2 (aneuploidy-related stress)", "TRAMETINIB": "5, tumour context", "SELUMETINIB": "5, tumour context; class sibling of trametinib",
                  "BORTEZOMIB": "2, tumour context", "CHLOROQUINE": "2, tumour context", "PHENYLBUTANOIC ACID": "2 and 3 (proteotoxic stress, senescence)", "SELINEXOR": "5, tumour context", "BOSUTINIB": "5, tumour context"}
    for n in order:
        p = by_name[n]
        s = p["signals"]
        sig = (f"line A rank {s['line_A_rank']} of {s['line_A_scored_drugs']}" if s["line_A_rank"] else "line A: not scored") + f"; B2 {'yes' if s['in_B2_robust_set'] else 'no'}; C {'yes' if s['in_line_C_hits'] else 'no'}"
        if p["group"] == "constitutional":
            hyp = "constitutional" + (" (mixed evidence)" if p["mixed_evidence"] else "")
        elif p["group"] == "class_sibling":
            hyp = "class sibling"
        else:
            hyp = "tumour context only"
        risk = "also a hazard for the constitutional mosaic (kills aneuploid cells selectively)" if p["hazard_for_constitutional_mosaic"] else ("immunosuppression, growth and wound healing" if n in ("SIROLIMUS", "EVEROLIMUS") else ("tumour data point both ways" if n == "METFORMIN" else ("whole-chromosome trisomy models only; rat-pup neurotoxicity of the metabolite" if n == "PHENYLBUTANOIC ACID" else ("thrombocytopenia, weight loss" if n == "SELINEXOR" else "normal aneuploid cells may be more sensitive too"))))
        t.append(f"| {disp(n)} | {TIER[n]} | {short_cons.get(n, '')} | L{p['best_level']} | {risk} | {LABEL_SHORT.get(n, _label_line(p['label'])[:150])} | {sig} |")
    t.append("")
    for n in order:
        p = by_name[n]
        s = p["signals"]
        t.append(f"### 4.{order.index(n) + 1} {disp(n)}" + (f" (class sibling of {p['sibling_of'].title()})" if p["sibling_of"] else "") + "\n")
        if n == "PHENYLBUTANOIC ACID":
            t.append("*(Named phenylbutanoic acid in ChEMBL and in the drug universe used by the pipeline; the approved product is sodium phenylbutyrate.)*\n")
        t.append(f"**Why.** {p['link']}\n")
        cur = CURATED.get(n, {"for": [], "against": []})
        if cur["for"]:
            t.append("**Evidence for (verified sentences).**")
            for pmid, contains, maxlen in cur["for"]:
                t.append(bullet(lit.sentence(pmid, contains, maxlen)))
            t.append("")
        if cur["against"]:
            t.append("**Evidence against or risk (verified sentences).**")
            for pmid, contains, maxlen in cur["against"]:
                t.append(bullet(lit.sentence(pmid, contains, maxlen)))
            t.append("")
        if p["hazard_for_constitutional_mosaic"]:
            t.append("**Hazard.** The property that makes this drug a hypothesis for a tumour, killing aneuploid cells selectively, is also a hazard for the normal aneuploid cells of a child with constitutional mosaic aneuploidy. It is **not** proposed for the constitutional state.\n")
        if n == "EVEROLIMUS":
            t.append("**Rhabdomyosarcoma.** One study reports partial responses in two patients with rhabdomyosarcoma in a trial of everolimus with lenvatinib in children and young adults with relapsed solid tumours (PMID 40313040; level 0 here because it is a tumour observation and not a test of the mechanism). It is mentioned because the child's record documents rhabdomyosarcoma (current tumour status not available); it is not an indication.\n")
        if n == "SIROLIMUS":
            t.append("**Why level 1 and not level 2.** The strongest result (PMID 41820377) is a Drosophila model in which the checkpoint genes *bub3* or *rod* are depleted by RNAi in neural stem cells (the full text does not use a BubR1 allele). The sentence on rapamycin reports a result similar to the genetic Tor and autophagy interventions on the number of aneuploid neuroblasts, and the paper states that the genetic interventions did not rescue the stem-cell marker or the brain size. "
                     "That is a related system with a measured effect (L1). L2 requires BUB1B/BUBR1-deficient or MVA-derived human cells, or a hypomorphic BubR1 model with rescue of the phenotype, and none of these exists for sirolimus. Two limits apply to the rapamycin sentence: the result it reports is shown in a supplementary figure (Supplementary Fig. 6C of PMID 41820377) that was not examined here, so the in-text sentence is the only evidence read; and the inactivation of TORC1 caused chromosome instability only in yeast, while in human cells only mitotic slippage was reported (PMID 35141499).\n")
        if n == "EVEROLIMUS":
            t.append("**Link with the computational class.** Lenvatinib, the other drug of the rhabdomyosarcoma trial above, belongs to ATC L01EX, the one class that passes in two network lines (section 6.6). The trial connects a literature proposal and a computational class by coincidence of drugs, not by mechanism; both are level 1 or lower, and neither is an indication.\n")
        if n == "EVEROLIMUS":
            t.append("**What differs from sirolimus.** The class rationale and the refuting experiment are shared; the evidence of its own is weaker (a rat sarcopenia study and a paediatric trial observation, both level 0), none with a BubR1 defect. Its paediatric label is for a brain tumour in tuberous sclerosis (Afinitor, from 1 year), not for transplantation as with sirolimus, and the transplant label (Zortress) says paediatric use below 18 years is not established. In the experiment, sirolimus is the comparator.\n")
        if n == "METFORMIN":
            t.append("**This is the weakest proposal of tier A.** Its support is one human cell study of oxidative damage in trisomy 21 fibroblasts, and its own evidence against shows that it changes BUBR1 and MAD2 expression in cancer cells, in a direction the abstract does not state. The abstract of that study reports the effect in the trisomy-21 lines at 10 to 50 µM and does not say whether the normal control fibroblasts responded, so it cannot show specificity for aneuploidy, and the cancer-cell study of the evidence against uses 1 to 16 mM. It passes R2 on that single level-1 finding; it is listed also because its label has a paediatric indication from 10 years and its safety profile is well known, which is a reason of convenience and not a rule of E16.\n")
        safety = json.loads(SAFETY_PATH.read_text(encoding="utf-8"))["entries"].get(n, []) if SAFETY_PATH.exists() else []
        if safety:
            t.append("**Label safety statements (verbatim, read 2026-10-08; not a full safety review).**")
            for q in safety:
                t.append(f"- «{q['quote']}» ({q['section']}; {q['product']})")
            if n in ("SIROLIMUS", "METFORMIN"):
                t.append("- These statements come from a repackager's label. The same sentences were found by code in a manufacturer's label of the same drug (" + ("sirolimus oral solution, Mylan" if n == "SIROLIMUS" else "metformin oral solution, VistaPharm") + "); the innovator labels (" + ("Rapamune" if n == "SIROLIMUS" else "Glucophage") + ") were not found in DailyMed on 2026-10-08.")
            t.append("")
        extra = (" A second label was read: " + ", ".join(x["product"] for x in p["label"]["also_read"]) + ".") if p["label"] and p["label"].get("also_read") else ""
        t.append(f"**Paediatric label (read 2026-10-08).** {_label_line(p['label'])}{extra}")
        t.append(f"\n**What the computational filters say.** " + (f"Line A rank {s['line_A_rank']} of {s['line_A_scored_drugs']}" if s["line_A_rank"] else "Line A cannot score it (not in the knowledge graph)") + (f" (in {s['random_query_fraction_in_top50']:.0%} of random queries' top 50)" if s["random_query_fraction_in_top50"] else "") +
                 f"; B2 robust set: {'yes' if s['in_B2_robust_set'] else 'no'}; line C hit: {'yes' if s['in_line_C_hits'] else 'no'}; line D note (exploratory): {s['line_D_note_exploratory'].replace('_', ' ')}; "
                 f"curated protein target in ChEMBL: {'yes' if s['has_curated_target'] else 'no, so the network lines cannot score it'}. None of the proposals is supported by two usable lines, so the case for it rests on the mechanism and the literature above.\n")
        t.append(f"**Experiment that would refute it.** {p['experiment']}\n")
        who = WHO["hazard"] if TIER[n].startswith("C") else (WHO["tumour"] if TIER[n].startswith("B") else WHO["constitutional"])
        t.append(f"**Who could run it, with what.** {who}\n")
    t.append("### §EXCL§ The other drugs with verified entries, and why they are not proposed\n")
    prop = set(by_name)
    cau = [c for c in proposals.get("caution", []) if c["drug"] not in prop]
    cau_names = {c["drug"] for c in cau}
    rows = ["| Drug | Verified levels | Disposition under the rules |", "|---|:-:|---|"]
    for c in sorted(cau, key=lambda c: c["drug"]):
        lv = sorted({int(e["level"]) for e in c["evidence"]})
        rows.append(f"| {c['drug'].title()} | {', '.join('L' + str(x) for x in lv)} | R4, caution list: {REASON.get(c['drug'], 'flagged as a caution by the search')} |")
    for c in sorted(proposals.get("considered_not_proposed", []), key=lambda c: c["drug"]):
        if c["drug"] in cau_names:
            continue
        rows.append(f"| {c['drug'].title()} | {', '.join('L' + str(x) for x in c['levels'])} | not proposed: {EXCLUDED_WHY.get(c['drug'], c['why_not'])} |")
    t.append(f"Of the {proposals['n_drugs_with_verified_entries']} approved drugs with at least one verified entry, {len(by_name)} are in the tiers above and the rest are here, so the post hoc selection can be audited.\n")
    t.extend(rows)
    t.append("")
    return "\n".join(t).replace("§EXCL§", f"4.{len(order) + 1}")


def caution_literature_section(lit: Literature) -> str:
    rows = [
        ("Dactinomycin", "3", "31081598", "dactinomycin", "A case with biallelic BUB1B mutations and bilateral Wilms tumour could not continue chemotherapy."),
        ("Treosulfan", "3", "31053147", "graft rejection", "A child with MVA1 and myelodysplastic syndrome rejected the graft after a treosulfan-based conditioning (one case); the paper does not attribute the rejection to the drug, and treosulfan has a paediatric indication from 1 year in the label review."),
        ("Radiation and chemotherapy in general", "0", "39264246", None, "Not specific to MVA: it describes the group of genomic-instability disorders, and no radiosensitivity experiment in MVA or BUB1B was found."),
        ("Chemotherapy in MVA", "0", "33209710", None, "A review states that the anomalies of MVA may complicate chemotherapeutic regimens."),
        ("Taxanes (paclitaxel, docetaxel, cabazitaxel) and vinca alkaloids", "2", "39014450", None, "Low BubR1 or a weak checkpoint is associated with lower sensitivity (also PMID 17234259 and 38987356); no study was found that shows higher sensitivity with low BubR1."),
        ("Blocking p53, or eliminating senescent cells (for example dasatinib with quercetin)", "0", "20663956", None, "Loss of p53 accelerates tumours in aneuploidy models; BubR1 mutant mice show enhanced growth of carcinogen-induced tumours (PMID 23300461); some senolytics may promote tumorigenesis (PMID 38612842)."),
        ("Growth hormone and mecasermin", "3", "35434947", "unsatisfactory", "The only clinical observation of growth hormone in MVA is one girl with MVA type 2 (CEP57, not BUB1B) with an unsatisfactory response and a high IGF-1; the label of growth hormone contraindicates active malignancy (section 7.2)."),
        ("PARP inhibitors, etoposide", "1", "39251587", "olaparib", "Aneuploid human cells were more resistant to these agents, which argues against a benefit from DNA-damage drugs rather than for one."),
    ]
    t = ["### 7.1 Caution from the literature\n", "| Drug or class | Level | What the verified literature says | Reading |", "|---|:-:|---|---|"]
    for name, level, pmid, contains, reading in rows:
        e = lit.sentence(pmid, contains, maxlen=200)
        t.append(f"| {name} | {level} | «{e['text']}» (PMID {pmid}; {e['system']}) | {reading} |")
    t.append("\nThese are cautions for investigation and for a clinician. Most rest on level-0 or level-1 evidence; the two level-3 entries are single case reports.\n")
    return "\n".join(t)
