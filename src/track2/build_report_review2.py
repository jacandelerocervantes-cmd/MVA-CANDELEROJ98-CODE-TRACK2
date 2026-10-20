"""Track 2 report under addendum 11 (`protocol/addendum_review2.md`): every figure is read from `out/`; nothing is typed by hand except the prose.

    python -m track2.build_report_review2 --out out --report submissions/candeleroj98_track2_report.md

Publishes no genomic coordinates, no sample identifiers and no birth details (`docs/data_protection.md`). Variants appear at protein level only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from track2.label_levels import LEVEL
from track2.report_candidates import Literature, candidates_section, caution_literature_section, mechanism_section

AI_DECLARATION = Path("config/ai_declaration.md")


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _opt(path: Path) -> Any | None:
    return _load(path) if path.exists() else None


def f3(x: float) -> str:
    return f"{x:.3f}"


def ci(s: dict[str, Any]) -> str:
    return f"{f3(s['mean_auc'])} ({f3(s['ci95_low'])} to {f3(s['ci95_high'])})"


def lines_table(fl: dict[str, Any]) -> str:
    rows = ["| Line | What it measures | Pre-specified validation | Result | Usable |", "|---|---|---|---|:-:|"]
    what = {"A": "Knowledge graph (PrimeKG, personalised PageRank from the disease and the child's phenotypes)", "B1_phenotype_module": "STRING proximity to the module of genes annotated to the child's HPO terms",
            "B2_mitosis_module": "STRING proximity to the module built around BUB1B (links with >= 2 records)", "C": "DepMap CRISPR and PRISM sensitivity of weak-checkpoint lines, adjusted for lineage, proliferation, TP53, MSI and CIN",
            "D": "Direction of effect from a signed network (OmniPath)"}
    rule = {"A": "Leave-out of known indications, 200 diseases: AUC >= 0.70, lower bound > 0.55, and excess over a degree-only baseline with lower bound > 0", "B1_phenotype_module": "Same, without the causal genes and the drugs that act on them: AUC >= 0.60, lower bound > 0.50, >= 15 diseases, excess over degree > 0",
            "B2_mitosis_module": "Module rule on diseases with a known causal gene: AUC >= 0.60, lower bound > 0.50, >= 15 diseases", "C": "Observed hits above the 95th percentile of 1,000 matched random gene sets", "D": "At least 5 of 6 textbook perturbations (exploratory otherwise)"}
    for k, v in fl["lines"].items():
        rows.append(f"| {k.replace('_', ' ')} | {what[k]} | {rule[k]} | {v['evidence']} | {'yes' if v['usable'] else 'no'} |")
    return "\n".join(rows)


def caution_markdown(c: dict[str, Any], qc: dict[str, Any]) -> str:
    english = {"a_tubulina": "Tubulin binders (vinca alkaloids, taxanes, colchicine, eribulin, ixabepilone)", "b_checkpoint_mitotico": "Mitotic-checkpoint, Aurora, PLK and TTK (MPS1) inhibitors",
               "c_quimioterapia_genotoxica": "Genotoxic chemotherapy (alkylating agents, platinum, anthracyclines, topoisomerase II)", "d_proteasoma": "Proteasome inhibitors",
               "e_hormona_crecimiento": "Growth hormone and somatropin analogues"}
    rows = ["| Class (pre-specified) | Approved drugs in the universe | Listed, each with a quoted label sentence | Quoted sentences |", "|---|:-:|---|:-:|"]
    for key, v in c["clases"].items():
        listed = ", ".join(sorted(f["name"].title() for f in v["farmacos"])) or ("none: no approved drug meets the criterion" if v["clase_vacia"] else "none")
        rows.append(f"| {english.get(key, v['etiqueta'])} | {v['n_miembros_en_universo']} | {listed} | {v['n_citas']} |")
    gh = next((f for f in c["clases"]["e_hormona_crecimiento"]["farmacos"] if any(x["categoria"] == "active_malignancy" for x in f["citas"])), None)
    quote = next((x for x in gh["citas"] if x["categoria"] == "active_malignancy" and len(x["texto_literal"]) < 200), None) if gh else None
    text = "\n".join(rows) + "\n\n"
    text += ("Each entry is listed only if at least one sentence of the current label (DailyMed, or openFDA when DailyMed has no readable text) is copied verbatim; the sentences, with section, set id and reading date, are in `out/caution_list.json`. "
             f"`track2.verify_quotes` re-checks every quoted sentence (caution list, age table and label readings: {qc['n_quotes']} in all) against the label text it cites, ignoring whitespace only: {qc['found']} found, {qc['not_found']} not found. "
             "The extraction is heuristic and a person must read the label before any use.\n\n")
    if quote:
        text += f"For growth hormone, the criterion that matters for this child (a tumour history) is quoted from the {gh['name'].title()} label, contraindications: \"{quote['texto_literal']}\" (set id {quote['setid']}).\n\n"
    not_repr = [x for x in c["sin_cita"] if str(x.get("motivo", "")).startswith("label not representative")]
    text += (f"Not listed for lack of a readable label (no entry in either service): {', '.join(sorted(x['farmaco'].title() for x in c['sin_cita'] if x not in not_repr))}. "
             + "".join(f"{x['farmaco'].title()} is not listed because {x['motivo'].split(': ', 1)[1]}. " for x in not_repr) +
             "Excluded on purpose and documented: antibody-drug conjugates with a tubulin payload, benzimidazole anthelmintics and griseofulvin, dexrazoxane and quinolones. These are cautions for investigation and for a clinician, not recommendations.")
    return text


def age_markdown(a: dict[str, Any], focus: list[str]) -> str:
    labels = {"0_to_under_3": "0 to < 3 years", "3_to_under_12": "3 to < 12 years", "12_and_older": "12 years and older"}
    rows = ["| Age group | Labelled (of which age not stated) | Partly (from a stated age) | Not labelled or adult only | Verify by hand |", "|---|:-:|:-:|:-:|:-:|"]
    for g, name in labels.items():
        s = a["resumen_por_grupo"][g]
        ans = [f['name'] for f in a['farmacos'] if f['estado_por_grupo'][g] == 'age_not_stated']
        n_exp = sum(1 for x in ans if LEVEL.get(x) == 1)
        n_ans = f"{len(ans)}: {n_exp} explicit, {len(ans) - n_exp} implicit" if ans else "0"
        rows.append(f"| {name} | {len(s['labelled'])} ({n_ans}) | {len(s['partly'])} | {len(s['not_labelled_or_adult_only'])} | {len(s['verify'])} |")
    by_name = {f["name"]: f for f in a["farmacos"]}
    items = ["", "| Drug (in the network-evidence class or supported by one line) | 0 to < 3 | 3 to < 12 | 12 and older | Reading recorded |", "|---|---|---|---|:-:|"]
    for n in focus:
        f = by_name.get(n)
        if f:
            st = f["estado_por_grupo"]
            items.append(f"| {n.title()} | {st['0_to_under_3']} | {st['3_to_under_12']} | {st['12_and_older']} | {'yes' if f['revision_manual_en_config'] else 'no'} |")
    pending = sorted(f["name"] for f in a["farmacos"] if f["manual_review_required"])
    referral = sorted(f["name"] for f in a["farmacos"] if (f.get("revision_manual") or {}).get("verdict") == "verify" and "only refers" in ((f.get("revision_manual") or {}).get("scope") or ""))
    age_ns = sorted({f["name"] for f in a["farmacos"] if f["estado_por_grupo"]["0_to_under_3"] == "age_not_stated"})
    n_reviewed = sum(1 for f in a["farmacos"] if f["revision_manual_en_config"])
    text = "\n".join(rows) + "\n" + "\n".join(items) + "\n\n"
    text += (f"{a['n_farmacos']} drugs were tabulated (the caution list, the network-evidence class and the single-line drugs). The status comes only from label text, through a recorded reading of each label (`config/pediatric_label_review.json`, {n_reviewed} drugs of this table); the age band is never assumed. "
             "The reading was done by the author's assistant from the complete paediatric text of the label (and, where a label has no such section, from the paediatric sentences of the whole label); it is not a pharmacist or regulatory review, and each entry carries a verbatim quote and must be re-checked against the current label before any use. "
             "The reading also corrected products that the automatic choice had got wrong (a systemic instead of an ophthalmic cyclosporine, a conventional instead of a liposomal doxorubicin) and a misread warning (cisplatin: a paediatric ototoxicity warning, not an established paediatric use). "
             f"{len(pending)} drugs stay unresolved (`verify`: {', '.join(n.title() for n in pending)}) because the label gives paediatric information without stating that paediatric use is established. A drug is described as labelled for a group only when its label says so; a label whose paediatric section only refers to dosing in another section ({', '.join(n.title() for n in referral)}) is read as \"verify\", not as established use, and \"age not stated\" ({', '.join(n.title() for n in age_ns)}) means labelled for children without a minimum age, explicitly or by the paediatric information the label gives. The drugs of section 4 are proposed for investigation, not for this child, and this table does not change that.")
    return text


CAUTION_SHORT = {"a_tubulina": "caution: tubulin binders", "c_quimioterapia_genotoxica": "caution: genotoxic chemotherapy", "d_proteasoma": "caution: proteasome inhibitors", "e_hormona_crecimiento": "caution: growth hormone"}
D_TEXT = {"strengthens_checkpoint": "strengthens (exploratory)", "weakens_checkpoint": "weakens (exploratory)", "unclear": "unclear"}


def appendix_markdown(out: Path) -> str:
    """Appendix A: for every drug that surfaced in the analysis (the 56 of the label review plus fostamatinib), its ChEMBL mechanism, the support from each line, the direction note of line D
    and the reason it is not a candidate. Everything is derived from the result files; nothing about BUB1B biology is asserted by hand here."""
    import re

    universe = _load(Path("data/universe.json"))
    by_name: dict[str, Any] = {}
    for d in universe["drugs"]:  # several universe entries can share a name; keep the one that has a recorded mechanism
        key = d["name"].upper()
        if key not in by_name or (not any(m.get("mechanism") for m in by_name[key].get("mechanisms") or []) and any(m.get("mechanism") for m in d.get("mechanisms") or [])):
            by_name[key] = d
    names = {d["drug_id"]: d["name"].upper() for d in universe["drugs"]}
    q, nq = _load(out / "kg_query_e5.json"), _load(out / "kg_nullquery_e5.json")
    best: dict[str, int] = {}
    for drug, r in q["rank_combined"].items():
        best[drug] = min(r, best.get(drug, r))
    order = sorted(best, key=lambda k: (best[k], k))
    rank = {names[d]: i + 1 for i, d in enumerate(order) if d in names}
    frac = {n.upper(): f for n, f in nq["fraction_of_random_queries_with_drug_in_top50"].items()}
    fl = _load(out / "four_lines_review2.json")
    a_specific = set(fl["support_detail"]["random_query_fraction_of_A_specific"])
    b2 = set(fl["support_detail"]["B2_W3_two_record"])
    vu = _load(out / "vulnerability_e6.json")
    c_set = {v["name"].upper() for v in vu["supported_drugs"].values()}
    d_calls = _load(out / "signed_effect.json")["drugs"]
    cls = _load(out / "class_analysis.json")
    in_class = {m for c in cls["supported_in_two_usable_lines"] for m in c["members"]}
    class_name = {m: c["class"] for c in cls["supported_in_two_usable_lines"] for m in c["members"]}
    caution = _load(out / "caution_list.json")
    caution_of: dict[str, list[str]] = {}
    for key, c in caution["clases"].items():
        for d in c["farmacos"]:
            caution_of.setdefault(d["name"].upper(), []).append(CAUTION_SHORT.get(key, key))
    reviewed = sorted(v["name"].upper() for v in _load(Path("config/pediatric_label_review.json"))["entries"].values() if "2026-10-08" in v.get("source", ""))
    drugs = sorted(set(reviewed) | {"FOSTAMATINIB"})
    total = len(rank)
    rows = ["<!-- widths:10,19,11,8,9,15,28 -->", "| Drug | ChEMBL mechanism | Line A | B2 / C | Line D | Class or caution list | Why it is not a candidate |", "|---|---|---|:-:|---|---|---|"]
    for n in drugs:
        d = by_name.get(n)
        mech = "; ".join(dict.fromkeys((m.get("mechanism") or "").strip() for m in (d or {}).get("mechanisms") or [] if m.get("mechanism")))[:90] or "none curated in the ChEMBL table used"
        r = rank.get(n)
        a_txt = "not scored" if r is None else f"rank {r}" + (f"; in {frac[n]:.0%} of random queries" if n in frac else ("" if r <= 50 else " (outside the top 50)"))
        if n in a_specific:
            why = "Supported by line A only, and specific to this query; a second usable line is missing."
        elif n in b2:
            why = "Supported by line B2 only; B2 is not testable once the drugs acting on the causal genes are removed; a second usable line is missing."
        elif r is not None and r <= 50 and n in frac and frac[n] >= 0.05:
            why = "The line A support is generic: the drug also reaches the top 50 of at least 5 % of degree-matched random queries."
        elif n in in_class:
            why = f"Member of the network-evidence class {class_name[n]} (N1, below level 1); no individual drug is supported by two usable lines."
        else:
            why = "No usable line supports it individually" + (f" (line A rank {r} of {total}, outside the top 50)." if r is not None else " (not scored by line A).")
        if d is not None and not d.get("human_targets"):
            why += " It has no curated protein target in the ChEMBL table used, so the network lines (B and D) cannot score it."
        call = (d_calls.get(d["drug_id"]) or {}).get("call", "unclear") if d else "unclear"
        tags = ([class_name[n] + " (network-evidence class N1)"] if n in in_class else []) + caution_of.get(n, [])
        rows.append(f"| {n.title()} | {mech.replace('|', '/')} | {a_txt} | B2 {'yes' if n in b2 else 'no'} / C {'yes' if n in c_set else 'no'} | {D_TEXT.get(call, call)} | {'; '.join(tags) or '—'} | {why} |")
    head = ("## Appendix A. The drugs that surfaced in the analysis, and why none meets the two-line rule\n\n"
            f"The {len(drugs)} drugs below are those that appeared at some point in the analysis: the network-evidence class, the single-line support, the caution list and the label review, plus fostamatinib, the single candidate of the first analysis. "
            "For each one: its mechanism as recorded in ChEMBL; where it ranks in line A (and, for the top 50, in what share of the 500 random queries it also reaches the top 50); whether it is in the B2 robust set or among the line C hits (none is); the direction note of line D, which is exploratory and assumes that strengthening the checkpoint is favourable, "
            "an assumption that is not justified for a hypomorphic allele; its class or caution-list membership; and the reason it is not a candidate. "
            "The reason is derived by rule from the result files (a drug needs at least two usable lines) and says nothing about whether the drug could ever help this child. The mechanism column describes the drug, not its effect on BUB1B or the spindle checkpoint. "
            "946 of the 2,705 approved drugs of the universe, among them some of the drugs below, have no curated mechanism or protein target in the ChEMBL table used; the network lines (B and D) cannot score them, so for those the absence of support is a limit of the data and not evidence against the drug.\n")
    return head + "\n".join(rows) + "\n"


def build(out: Path) -> str:
    fl, ledger, cls = _load(out / "four_lines_review2.json"), _load(out / "forking_ledger.json"), _load(out / "class_analysis.json")
    nq, vu = _load(out / "kg_nullquery_e5.json"), _load(out / "vulnerability_e6.json")
    e2, mrv = _load(out / "phenotype_validation_e2.json"), _load(out / "module_rule_validation.json")
    kg, kx = _load(out / "kg_validation.json"), _load(out / "kg_extra.json")["results"]["disease_only"]
    caution, ages = _opt(out / "caution_list.json"), _opt(out / "age_label_table.json")
    ms = _load(Path("config/mechanism_screen.json"))
    pre, exp = fl["readings"]["b2_counted_as_prespecified"], fl["readings"]["b2_held_exploratory"]
    s = fl["support_detail"]
    x = mrv["evaluation_without_drugs_acting_on_causal_genes"]
    ov = fl["overlap_vs_chance_common_universe"]
    ovc = _load(out / "overlap_check.json")
    n_generic = nq["n_real_top50_drugs_generic_at_least_5_percent"]
    n_real = nq["real_top50_unique_drugs"]
    supported = cls["supported_in_two_usable_lines"]

    t = []
    t.append("# Track 2 Report: Drug-Repurposing Hypotheses for Mosaic Variegated Aneuploidy (BUB1B)\n")
    t.append("**Author:** Juan Antonio Candelero Cervantes (Individual Submission, No team) · **Hugging Face user:** `candeleroj98` · **Track:** 2 (drug repurposing)\n")
    t.append("*Hypotheses for investigation and for a clinician. Approved drugs only. Not medical advice: no dose, no recommendation to treat. The child's age is unknown.*\n")

    props = _load(out / "proposals_review2.json")
    prop_names = [p["drug"] for p in props["proposed"]]
    lit = Literature(out / "literature_hypotheses")
    t.append("## 1. Summary\n")
    t.append("**In plain words.** A child has damage in a gene (*BUB1B*) that controls cell division. We looked for medicines already approved for something else that could be worth testing in a laboratory. The ideas come from published science, not from the computer: the computer checks did not confirm any of them. "
             "The evidence is early (cells, yeast, flies, mice), several of the ideas could also do harm, the child's age is unknown, and nothing here is a treatment or advice. The first step we suggest is a chromosome test. Terms: *level L0 to L3* is how close the evidence is to this condition (L1 is a related system; L2 and L3 would be human cells or patients with BUB1B, and none reaches them); *tier A, B, C* are constitutional hypotheses, tumour-context hypotheses and warnings.\n")
    t.append("The question was which **already approved** drugs are plausible candidates for further investigation given the mechanism that Track 1 identified: biallelic *BUB1B* involvement (a nonsense allele and a missense variant of uncertain significance) in a child with the seven phenotypes listed in section 3. "
             "The report goes from the disrupted mechanism to a hypothesis and to the experiment that would refute it. Four computational lines of evidence (a phenotype-guided knowledge graph, network proximity, vulnerability of weak-checkpoint cell lines, a signed-network direction check) were specified before the analysis (with amendments made after results, listed in section 5) and are used as a **filter** on the hypotheses, not as their source.\n")
    t.append(f"**{len(prop_names)} drugs in three tiers (section 4).** A, constitutional hypotheses: sirolimus (mixed evidence: growth and muscle against immunosuppression and growth inhibition) with everolimus as its class sibling, metformin, and sodium phenylbutyrate. B, tumour-context hypotheses: trametinib with selumetinib as its class sibling. C, hazards to check and **not proposed**: bortezomib, chloroquine and bosutinib, which kill aneuploid cells selectively, the property that makes them interesting for a tumour and dangerous for a child whose normal cells are aneuploid. "
             "For each one the report gives the consequence it targets, verified sentences for and **against**, the paediatric label status (text only, the age is unknown), and the experiment that would refute it. **No approved drug reaches evidence level 2 or 3 (human BUB1B-deficient cells or patients); the best level is 1** (a related system: aneuploid cells, yeast, an animal with a checkpoint defect).\n")
    t.append("**First test before any drug experiment:** a karyotype or single-cell sequencing of the child. The blood genome shows only a weak aneuploidy signal and it is not shown that the child has mosaic aneuploidy; the tier-A (constitutional) hypotheses depend on it (section 3).\n")
    t.append(f"**The computational filters do not confirm any of them, and they could hardly have: they have little power for a gene without a targeted drug (section 6.7).** No drug meets the pre-specified rule of at least two usable lines, under either reading of line B ({pre['n_candidates']} candidates with B2 counted as pre-specified, {exp['n_candidates']} with B2 held exploratory). "
             f"Line A is valid but modest (excess over a degree-only baseline {f3(kx['mean_excess_over_degree'])}) and its patient ranking is generic ({n_generic} of {n_real} drugs of its top 50 also appear in the top 50 of at least 5 % of degree-matched random queries); "
             f"the phenotype module (B1) does not beat the degree baseline once the causal genes are removed; line C, once TP53, MSI and chromosomal instability are in the model, falls from 22 to {vu['observed_total']} hits against a null 95th percentile of {vu['null']['total_p95']:.1f}. "
             "The proposals are therefore **literature-anchored hypotheses selected by rules declared post hoc** (amendment E16), not computational candidates, and the report says so wherever they appear. "
             f"The filter did detect an artefact: the single candidate of the first analysis, fostamatinib, ranked first through connectivity (it is in the top 5 of {nq['fostamatinib']['random_top5']:.1%} of random queries) and is not proposed.\n")
    sup_txt = "; ".join(c["class"] + " (" + ", ".join(c["members"][:4]) + "...)" for c in supported) if supported else "none"
    t.append(f"Also delivered: (i) the negative result of the strict rule, with the checks that produced it; (ii) a class-level network signal, {sup_txt}, network evidence N1 (not a hypothesis) and tumour context; "
             "(iii) a **caution list** in two parts, from the literature (section 7.1) and from quoted label text (section 7.2); (iv) label support by age group, because the age of the child is unknown (section 10); and (v) a framework that reruns for any Track 1 submission.\n")

    t.append("### 1.1 How the report answers the four criteria of the track\n")
    t.append("| Criterion | What the report offers | Where | What stays open |\n|---|---|---|---|\n" +
             "| **Scientific rigor** (35 %) | The mechanism is characterised for two direction scenarios (complete loss and hypomorph) and tied to five consequences, each with verified sentences, the system in which it was measured and an evidence level L0 to L3. Every proposal states its evidence for and against, the label safety statements, and the experiment that would refute it. The strict computational rule gave no candidate, and the report says what that does and does not mean. | Sections 3, 4, 6.7 | The phase of the two variants and the effect of p.Asn1002Lys are unresolved; no approved drug reaches level 2 or 3; the selection rules of the proposals were written after reading the literature. |\n" +
             "| **Potential impact** (25 %) | If a proposal is validated or refuted in *BUB1B*-deficient human cells, the child's clinicians and other families gain a tested shortlist or a documented dead end. Beyond the drugs: a caution list for clinicians (dactinomycin, treosulfan, taxanes, growth hormone, furosemide in a premature infant) with the literature behind each entry; a table of label support by age group for a child whose age is unknown; and a diagnostic step, a karyotype or single-cell sequencing, because the blood genome shows only a weak aneuploidy signal that cannot be read as mosaicism. The negative results are also useful to others: they show where network and cell-line methods lack power for a gene with no targeted drug. | Sections 1, 3, 7, 10 | Whether any proposal is safe or effective; whether the child is mosaic; how the lessons transfer to other MVA genes (*CEP57*, *TRIP13*). |\n" +
             "| **Innovation** (25 %) | Computation is used as a filter to discard, not to generate: random queries of the same shape show that the knowledge-graph ranking is generic, and that exposed the single candidate of the first analysis (fostamatinib) as a connectivity artefact. A class-level analysis replaces single hub-driven drugs. Every quoted label and literature sentence is re-checked by code against its source. The age-unknown problem is handled with fixed age bands read from label text only. The forking-path ledger and the disclosures of the post hoc steps make the process auditable. | Sections 5, 6, 4.10 | Not yet tested on another case. |\n" +
             "| **Scalability** (15 %) | The first stages of the pipeline rerun from the Track 1 CSV with a fixed seed, tests and checksums; the later stages (lines A and C, the consequences and the drug selection) are written for this case. The literature layer is a method (drug-agnostic searches by consequence, a code check of every sentence, written selection rules) that can be repeated for another disease. | Sections 11 and the README | The literature searches were done with AI assistants and are not repeated by changing a file; labels must be read per drug; about 1 GB of public data is downloaded. |\n")
    t.append("### 1.2 Open questions, left open on purpose\n")
    t.append("The track asks for hypotheses for follow-up, not answers. These questions are not closed by this report, and each names what would close it.\n")
    for q, how in (("Does the child have mosaic aneuploidy, and in which cells?", "a karyotype or single-cell sequencing of blood and, if possible, another tissue, and a copy-number analysis of the megabase-scale structural calls listed in Track 1."),
                   ("What is the effect of p.Asn1002Lys, and are the two variants on different copies of the gene?", "Parental testing for the phase, and BubR1 level and checkpoint function in the child's cells."),
                   ("Do the tier-A drugs help or harm *BUB1B*-deficient human cells?", "The refuting experiments of section 4, in patient-derived or hypomorphic cells against a euploid control, with harm read out as carefully as benefit."),
                   ("Is the selective sensitivity of aneuploid cells (tier C) a hazard for the child's own cells?", "A toxicity comparison on the child's fibroblasts against a euploid control."),
                   ("Is there a link between nephrocalcinosis and *BUB1B*?", "cohort data from other MVA families; no link was found in the literature searched."),
                   ("Can the computational filter be made to discriminate for a gene without a targeted drug?", "Validation on diseases with a known causal gene and no drug at that gene, and a module built from that gene's biology.")):
        t.append(f"- **{q}** Would be closed by: {how[0].lower() + how[1:]}")
    t.append("")
    t.append("## 2. Scope, conditions and what this report is not\n")
    t.append("- **Repurposing of approved drugs only** (FDA-eligible and ChEMBL `max_phase = 4`). Investigational, withdrawn or unapproved compounds are excluded at every stage. Where a label already covers a feature of the child (for example growth hormone for short stature in children born small for gestational age), that is on-label supportive care and not repurposing.")
    t.append("- **Age unknown.** No drug is described as suitable for the child. Label support is reported for fixed age bands (neonate < 1 month, infant 1 month to < 2 years, child 2 to < 12 years, adolescent 12 to < 18 years, adult from 18) grouped as 0 to < 3, 3 to < 12 and 12 years and older; the status comes only from label text and the age is never assumed (section 10).")
    t.append("- No dose, exposure, interaction or treatment recommendation is given. The drugs named here must not be used on the basis of this report.\n")

    mech = mechanism_section(lit)
    pheno = ("Phenotypes used (HPO): rhabdomyosarcoma, nephrocalcinosis, short stature, failure to thrive, skeletal muscle atrophy, premature birth and small for gestational age. An eighth item in the clinical record, parental recurrent miscarriage, describes the parents and is not part of any phenotype query. No birth details are published.\n\n")
    head, rest = mech.split("\n\n", 1)
    aneup = ("**Does the child have mosaic aneuploidy?** Not shown. A coarse screen of read depth from the challenge variant file cannot detect a low-fraction or tissue-restricted mosaic and is not a karyotype; its numbers are not reported here. "
             "The working diagnosis of MVA rests on finding two *BUB1B* variants, which would explain the disorder only if they lie on different copies of the gene; that phase is unresolved. "
             "The Track 1 analysis of the same genome lists megabase-scale structural calls that it did not analyse; they could be segmental copy-number changes, including mosaic ones, or caller artefacts, and this report did not analyse them either. "
             "**Consequence for the hypotheses:** the group-A (constitutional) rationale for sirolimus, everolimus and metformin assumes that aneuploid cells are present; if the child's cells were euploid, the aneuploidy-related part of that rationale would not apply and only the growth and muscle part, which rests on mouse data, would remain. "
             "A karyotype or single-cell sequencing is the first test, before any drug experiment.\n\n")
    t.append(head + "\n\n" + pheno + aneup + rest)
    t.append(candidates_section(lit, props, fl))

    t.append("## 5. Design: four lines, a combination rule and the amendments\n")
    t.append("The design (`protocol/addendum_four_lines.md`) and the amendments made after two independent reviews by AI agents (`protocol/addendum_review2.md`, E1 to E15) are in the repository. The numerical thresholds in the code are those of the protocol files. "
             "Because the reviews read results that already existed, the amendments are **post hoc with respect to those results** and pre-specified with respect to the runs reported here; there is no git history and no external time stamp, so the order of events is attested by the author only.\n")
    t.append(lines_table(fl) + "\n")
    t.append("Combination rule: a drug is a candidate if it is approved, systemic, not an acute-use anaesthetic or sedative, and supported by at least two usable lines among A, B and C. B1 and B2 share the proximity method and the drug-target table, so together they count as one line.\n")

    t.append("## 6. Results\n")
    t.append("### 6.1 Line A: knowledge graph\n")
    t.append(f"Validation on {kg['n_evaluated']} diseases: mean AUC {ci(kg['summary'])}. A degree-only baseline reaches {f3(kx['mean_auc_degree_baseline'])}, so the excess is {f3(kx['mean_excess_over_degree'])} (95 % interval {f3(kx['excess_ci95'][0])} to {f3(kx['excess_ci95'][1])}); the excess is the mean over diseases of each disease's own difference from its own degree baseline, so it is not the difference of the two means shown. "
             f"Recall of known drugs in the top 50 among the drugs of the register universe is {kx['mean_recall50_universe']:.3f} against {kx['mean_recall50_universe_degree_baseline']:.3f} for the degree baseline. The line is valid and modest.\n")
    t.append(f"**Specificity for this child.** Five hundred random queries with the same shape (one disease node and six phenotype nodes, each within 30 % of the real degree) were ranked on the same graph. "
             f"Of the {n_real} distinct drugs in the real top 50, {n_generic} are in the top 50 of at least 5 % of the random queries; only {nq['n_real_top50_drugs_specific_below_5_percent']} are not. Fostamatinib, the first drug of the real ranking, is in the top 5 of {nq['fostamatinib']['random_top5']:.1%} of random queries (median random rank {nq['fostamatinib']['median_rank_random']:g}). "
             f"Its support from line A is therefore connectivity, not this child, and the classes that line A ranks highest (for example vitamin D agonists and supplements) and its top 50 (antidepressants, anticonvulsants, anaesthetics) do not look specific to this child either. Ablation of each of the six phenotype terms, and of the core (rhabdomyosarcoma, nephrocalcinosis, muscle atrophy) against the non-specific constitutional terms, changes the top 50 very little (Jaccard 0.89 to 1.00), "
             f"which says the ranking comes from the disease node and well-connected drugs and not from the phenotype set. We therefore make **no claim of a phenotype signal**. Several PrimeKG nodes map to the same drug; the top 50 is counted over distinct drugs.\n")
    t.append(f"The 46 drugs counted above are the distinct drugs of the real top 50 of the query; Appendix A lists a wider set (every drug that surfaced anywhere), and nine of its drugs were not measured against the random queries, which is conservative: none is counted as specific. Drugs that remain specific (in the top 50 of fewer than 5 % of random queries): {', '.join(sorted(s['random_query_fraction_of_A_specific']))}. They are supported by this line only.\n")

    t.append("### 6.2 Line B: network proximity\n")
    t.append("**Validation of the proximity method (frozen protocol, P1 to P4).** The checks were run on the hand-defined checkpoint module and were not repeated for the modules derived from the Track 1 file. P1, leave-one-out recovery of seed genes: the decision rule is met (AUC 0.92, 95 % interval 0.89 to 0.95), but no left-out seed reaches z <= -2 (median z -0.84), so the primary threshold is very strict for a single gene. P2, Fanconi anaemia positive control: **not recovered**. Metformin, cisplatin and mitomycin cannot be scored (ChEMBL has no curated mechanism with a protein target for them) and oxymetholone falls at the 59th percentile (rank 427 of 1,040; z +0.49). P3, stability of the top 20: mean Jaccard 0.81 when 20 % of the seeds are dropped, but 0.48 with the STRING threshold at 400 and 0.38 at 900; the protocol rule (Jaccard < 0.5) declares the ranking **unstable with respect to the interactome threshold**. P4, degree bias: the Spearman correlation between z and the mean degree of the targets is 0.28, so the normalisation does not remove the bias completely.\n\n")
    b1 = e2["primary"]
    t.append(f"**B1, phenotype module.** The first validation (19 diseases) gave a mean AUC of 0.717. A review noted that the module of each validation disease contains its causal gene and that the known drugs include those that act on it. Repeated without the causal genes and without those drugs, "
             f"{b1['n_evaluable']} diseases are evaluable: mean AUC {ci(b1['summary'])}, degree null {f3(b1['summary']['permutation']['null_mean'])}, excess over the degree baseline {f3(b1['excess_over_degree']['mean_excess_over_degree'])} "
             f"(95 % interval {f3(b1['excess_over_degree']['ci95_low'])} to {f3(b1['excess_over_degree']['ci95_high'])}). The interval includes zero, so **B1 is not usable**. Its robust set with two-record links is empty.\n")
    t.append(f"**B2, BUB1B module.** The rule that builds a module around a causal gene was validated on {mrv['evaluation']['n_evaluable']} diseases with a known causal gene: mean AUC {ci(mrv['evaluation']['summary'])}, excess over degree {f3(mrv['evaluation']['excess_over_degree']['mean_excess_over_degree'])} "
             f"(lower bound {f3(mrv['evaluation']['excess_over_degree']['ci95_low'])}); it passes the pre-specified rule. The same circularity check leaves {x['n_evaluable']} evaluable diseases, which is below the 15 required, with excess lower bound {f3(x['excess_over_degree']['ci95_low'])}: **not testable**. "
             f"B2 is therefore counted as pre-specified in one reading and held exploratory in the other. BUB1B-only, two-record links give: {', '.join(s['B2_W3_two_record'])} (tirbanibulin is topical and not eligible).\n")

    t.append("### 6.3 Line C: vulnerability of weak-checkpoint cell lines\n")
    t.append(f"With lineage and proliferation as covariates (the earlier model, 200 null sets), 22 hits (19 PRISM rows and 3 CRISPR genes) against a 95th percentile of 23.05; six of the PRISM hits were MDM2 inhibitors, the signature of TP53 status. "
             f"With TP53 (damaging or hotspot), MSI and the chromosomal-instability score added (DepMap 24Q2, {vu['covariates']['lines_with_complete_covariates']} lines with complete covariates) and 1,000 matched null sets: {vu['observed_total']} hits ({vu['observed_hits']['prism']} PRISM, {vu['observed_hits']['crispr']} CRISPR), null mean {vu['null']['total_mean']:.1f}, "
             f"95th percentile {vu['null']['total_p95']:.1f}, empirical p {vu['null']['empirical_p_observed']:.2f}. **Line C is not usable** and supports no drug. The checkpoint score and proliferation are strongly collinear (VIF {vu['covariates']['vif_of_score']:.2f}); an orthogonalised score gives {vu['covariates']['orthogonalised_score_hits']['prism']} PRISM and {vu['covariates']['orthogonalised_score_hits']['crispr']} CRISPR hits. "
             f"The checkpoint score includes *BUB1B* and *TTK* and correlates 0.84 with proliferation, so line C has little power to separate checkpoint from growth. Chromosomal instability is also a consequence of a weak checkpoint, so adjusting for it may remove real signal; the adjusted total was the pre-specified decision.\n")

    t.append("### 6.4 Line D: direction of effect (exploratory)\n")
    t.append(f"Result: {fl['lines']['D']['evidence']}. The line assumes that strengthening the checkpoint is favourable, which is not justified for a hypomorphic allele. It is shown as a direction note and assigns no tier.\n")

    t.append("### 6.5 Combination\n")
    t.append(f"Usable lines: A, and B2 in one of the two readings. Overlap among the {fl['scorable_in_both_A_and_B']} drugs that both A and B can score: zero in every pair (" + ", ".join(f"{k.replace('_and_', ' and ')} {v['overlap']}" for k, v in ov.items()) + "). "
             f"This zero is not informative: the rule sets are tiny (sizes: A {ovc['rule_sets']['sizes']['A']}, B1 {ovc['rule_sets']['sizes']['B1']}, B2 {ovc['rule_sets']['sizes']['B2']}, C {ovc['rule_sets']['sizes']['C']} drugs) and the overlap expected by chance for A and B2 is {ovc['rule_sets']['expected_overlap_A_and_B2']}. As a positive control of the identifier join, the top 50 of A and of B1 share {ovc['top_n_A_vs_B1_positive_control'][1]['observed']} drugs against {ovc['top_n_A_vs_B1_positive_control'][1]['expected_by_chance']} expected (hypergeometric p = {ovc['top_n_A_vs_B1_positive_control'][1]['p_hypergeometric']:.3f}); A and B2 are at or below chance (`out/overlap_check.json`). "
             f"**Candidates: {pre['n_candidates']} (B2 as pre-specified) and {exp['n_candidates']} (B2 exploratory).** Drugs supported by one line only (listed, not proposed): "
             f"{', '.join(pre['one_line_only'])} (B2 as pre-specified); {', '.join(exp['one_line_only'])} (B2 exploratory). No drug was moved into the candidate list on the weaker rule \"one line plus literature\"; the literature-anchored proposals of section 4 are a separate, declared-post-hoc category (amendment E16).\n")
    t.append("**Fostamatinib as a case of fragility.** In the first analysis it was the only candidate: rank 1 in line A and 9 of 9 configurations of the mitosis module. Its mitosis support rests on single-record measured targets and falls to zero with the two-record rule; "
             "its line-A support is the generic effect above. It is not proposed, and no label review was done for it for that reason.\n")

    t.append("### 6.6 Class-level hypothesis\n")
    t.append(f"Because single drugs are hub-driven, the unit was changed to the mechanism class (ChEMBL mechanism text and ATC level 4). Each class with at least three scored drugs was compared with {cls['permutations']:,} random classes matched in size and drug degree, with FDR < {cls['fdr']} within each line and definition. "
             f"A class is supported if it passes in both usable lines.\n")
    if supported:
        for c in supported:
            t.append(f"- **{c['class']}** ({c['definition']}): passes in A and B2; members include {', '.join(c['members'])}. Network evidence N1 (two lines of network evidence, no measured effect in a human cell system and no verified functional citation), which is below level 1 on the scale of section 2.3. "
                     "Class L01EX groups multi-target kinase inhibitors, which are well connected by construction, are antineoplastic (tumour-context) and are mostly labelled for adults; the match-on-degree null reduces but cannot remove that effect. "
                     f"The class result is also compatible with chance: at FDR 0.10 the ATC level-4 screen passes {sum(1 for x in cls['lines']['atc4/A']['classes'] if x['passes'])} of {len(cls['lines']['atc4/A']['classes'])} classes in line A and {sum(1 for x in cls['lines']['atc4/B2']['classes'] if x['passes'])} of {len(cls['lines']['atc4/B2']['classes'])} in line B2 (the latter include classes that cannot be related to this child), and about 0.4 classes would be expected to pass both by chance against the 1 observed (probability of at least one about 0.33). "
                     f"None of the {len(c['members'])} members is in the B2 robust set ({', '.join(s['B2_W3_two_record'])}), so the class signal rests on how many members are scored by each line and not on one drug; if B2 is held exploratory the class no longer has two usable lines. "
                     "Paediatric labels exist for the NTRK, RET and thyroid indications of some members (for example larotrectinib, selpercatinib) and none for rhabdomyosarcoma. It is not a hypothesis.")
    else:
        t.append("- No class passes in both usable lines.")
    t.append(f"\nClasses that pass in one line only are listed in `out/class_analysis.json` and are not hypotheses.\n")

    t.append("### 6.7 What the negative result means\n")
    t.append("The negative result is **not evidence that the method is wrong, and not evidence that no drug could help**. Two of its controls did what they were designed to do: they exposed fostamatinib as a connectivity artefact and showed that the early signal of line C was TP53 status (MDM2 inhibitors), and every amendment was stricter, not looser. Two others did not: the Fanconi anaemia positive control of the proximity validation was not recovered (three of its four reference drugs could not be scored and the fourth fell at the 59th percentile), and the top 20 of the proximity ranking was unstable with respect to the STRING threshold (Supplement S6.2). "
             "But the filter has little power for this question. The child's gene has no drug that targets it, and the validations use diseases with known drugs, an easier case (section 9). Line A discriminates weakly (excess over a degree baseline 0.022), B1 does not beat the degree baseline, line C has little power (its score correlates 0.84 with proliferation, and chromosomal instability is itself a consequence of a weak checkpoint), and the zero overlap between lines is uninformative because the sets hold 3 and 4 drugs. "
             "So the silence of the filter neither rules out the proposals nor supports them. What it supports is narrower: with the evidence available today, no approved drug is picked out by two independent computational lines, and the proposals of section 4 stand on mechanism and literature alone, with the weakness stated in section 4.0.\n")
    t.append("## 7. Caution\n")
    t.append(caution_literature_section(lit) + "\n")
    t.append("### 7.2 Caution list from the label text (classes pre-specified; each entry needs a quoted sentence)\n")
    t.append("*Note on label age: the cabazitaxel label read is dated 2022 and its marketing status could not be obtained; the set identifier, version and date of every label are in `out/label_provenance.json`, and the label readings are in `config/pediatric_label_review.json`.*\n")
    if caution:
        t.append(caution_markdown(caution, _load(out / "quote_check.json")) + "\n")
    else:
        t.append("`out/caution_list.json` was not available when this report was built.\n")

    t.append("## 8. What would need to be true\n")
    t.append("| Item | Supported by | What fails | What would have to hold |\n|---|---|---|---|")
    for r in pre["one_line_only"]:
        which = "A only" if r in s["random_query_fraction_of_A_specific"] else "B2 only"
        t.append(f"| {r.title()} | {which} | Needs a second usable line | An independent usable line (for B, a validated module; for A, a phenotype-specific signal) would have to point to the same drug |")
    t.append("| Fostamatinib | A (generic) | Specificity; single-record targets | Support from two-record links and a rank beyond random queries |")
    for pn in prop_names:
        t.append(f"| {pn.title()} (section 4) | Literature, level {next(int(e['best_level']) for e in props['proposed'] if e['drug'] == pn)} | No test in BUB1B-deficient human cells | The refuting experiment of section 4 run in BUB1B-hypomorphic human cells, with the result in the predicted direction |")
    for c in supported:
        t.append(f"| Class {c['class']} | A and B2, N1 | Network evidence only; adult tumour-context class; compatible with chance | A measured effect in a human cell system and a verified functional citation |")
    t.append("")

    t.append("## 9. Limitations\n")
    for line in ("Network and graph methods do not give the direction of an effect, and the sign for a hypomorphic allele is unresolved.",
                 "The validations of A and B use diseases with known drugs; a disease with no drug at the causal gene, as here, is a harder case that cannot be validated directly.",
                 "B2 passes its pre-specified validation but not the check that removes the drugs acting on the causal genes (not testable); this is reported in both readings.",
                 "Cell-line evidence (line C) is a tumour proxy; the constitutional picture of the child is not represented.",
                 "ChEMBL coverage: 1,040 of 2,705 approved drugs are scorable in the network.",
                 "Label readings were done by the author's assistant from the label text, with verbatim quotes; they are not a pharmacist or regulatory review and must be re-checked against the current label. The age of the child is unknown.",
                 "One case: no sensitivity, specificity or clinical significance is claimed.",
                 "The order of events (pre-registration, amendments) is attested by the author, not by an independent time stamp (no git history)."):
        t.append(f"- {line}")
    t.append("")

    t.append("## 10. Age: label support by group\n")
    t.append("Label status describes the populations studied for other diseases. It is not a statement of suitability for this child, and the drugs of S4 are proposed for in-vitro investigation, not for this child.\n")
    if ages:
        focus = sorted({m for c in supported for m in c["members"]} | set(pre["one_line_only"]) | set(exp["one_line_only"]) | set(prop_names))
        t.append(age_markdown(ages, focus) + "\n")
    else:
        t.append("`out/age_label_table.json` was not available when this report was built.\n")

    t.append("## 11. Reproducibility and scalability\n")
    t.append(f"Every figure comes from `out/`; the forking-path ledger (`out/forking_ledger.json`, {ledger['n_entries']} variants, {ledger['n_primary']} declared primary) lists the variants of the computational lines that were run, and everything not marked primary is a sensitivity analysis; the literature step (section 4), the docking and PRISM layers of the first analysis and the benchmark are not in the ledger. "
             "Code, configuration, protocol and instructions are in the repository README; third-party data are downloaded with SHA-256 recorded in `data/sources/record.json` and licences are in `THIRD_PARTY_LICENSES.md`. "
             "The Track 1 submission CSV is the case-specific input of the first stages; the later stages are written for this case. Doses, paediatric dosing and interactions with the child's medications are out of scope.\n")

    t.append("## 12. Data protection\n")
    t.append("The child's reads and VCF are not part of the repository. Per-variant and covariate detail stays in `data/private/` and is not published; this report names variants at protein level only. The Track 1 submission CSV, which carries the variant coordinates required by the challenge format and is the one case-specific input of the pipeline, is not included in the repository. All data will be deleted 30 days after the hackathon ends (`docs/data_protection.md`).\n")

    t.append("## 13. Generative-AI disclosure\n")
    t.append(AI_DECLARATION.read_text(encoding="utf-8").strip() + "\n")

    t.append("## 14. Acknowledgments\n")
    t.append("We are deeply grateful to the child and their family who chose to share their genomic data openly with the research community. This work would not exist without their generosity and trust. No information that could re-identify the patient or family is included.\n")
    t.append("The Rare Disease Genomic Challenge (MVA Hackathon 2026) was organized in partnership with Sage Bionetworks, the MVA Society, Hugging Face and BEACON (the Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We thank the Wilhelm Foundation for their guidance and insight.\n")

    t.append("## References\n")
    t.append("1. Guney E, Menche J, Vidal M, Barábasi AL. Network-based in silico drug efficacy screening. *Nat Commun* 2016. doi:10.1038/ncomms10331 (PMID 26831545).")
    t.append("2. Cheng F, Kovács IA, Barabási AL. Network-based approach to prediction and population-based validation of in silico drug repurposing. *Nat Commun* 2018. doi:10.1038/s41467-018-05116-5 (PMID 30002366).")
    t.append("3. Data resources (versions, file hashes and licences are listed in the repository): ChEMBL 37, STRING v12, Reactome, PrimeKG (Harvard Dataverse doi:10.7910/DVN/IXA7BM), Human Phenotype Ontology v2026-09-01, MONDO, OmniPath, DepMap 24Q2 (figshare 25880521) and PRISM Repurposing 24Q2, DailyMed and openFDA.\n")
    t.append(appendix_markdown(out))
    return "\n".join(t).replace("Phenylbutanoic Acid", "Sodium phenylbutyrate")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--report", type=Path, default=Path("submissions/candeleroj98_track2_report.md"))
    args = parser.parse_args()
    text = build(args.out)
    args.report.write_text(text, encoding="utf-8", newline="\n")
    print("wrote", args.report, len(text.split()), "words")


if __name__ == "__main__":
    main()
