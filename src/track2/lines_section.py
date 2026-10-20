"""Sections 3.13 (candidates proposed by each line) and 3.14 (multi-parameter cross-check), the gene-to-mechanism lead of section 3.1,
the answers to the four scoring criteria (4.7) and the sentences of the abstract and the brief, written only from the JSON files of `out/`,
`config/multi_evidence_track1_gene_block.json` and `out/engine/`. Nothing but wording is typed here.

Each line proposes its own candidates with its native score; convergence between lines is information, not a requirement.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from track2 import engine_section as eng

BANDS = ("0-4", "5-11", "12-17")
AGE = {"labelled": "L", "partly": "P", "not_labelled": "N", "verify": "V", "adult_only": "A", "age_not_stated": "S"}
ORDER = ["LIT", "ENG", "A", "B1", "B2", "C", "D", "E", "ET1", "EP"]
SHORT = {"LIT": "LIT, literature layer", "ENG": "ENG, ranking engine", "A": "A, knowledge graph", "B1": "B1, phenotype module", "B2": "B2, mitosis-module proximity", "C": "C, cell-line vulnerability",
         "D": "D, signed direction", "E": "E, BUB1B knockout signature", "ET1": "E with the combined Track 1 signature", "EP": "E-P, patient signature"}


def _load(out: Path, name: str) -> Any | None:
    p = Path(out) / f"{name}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def _have(out: Path) -> bool:
    return all((Path(out) / f"{n}.json").exists() for n in ("candidates_by_line_mva_v2", "multi_evidence_mva", "multi_evidence_mva_completed", "line_c_adjusted", "line_e_signature_reversal"))


def _name(d: str) -> str:
    return d.title()


def _flag(flags: list[str] | None) -> str:
    f = " ".join(flags or [])
    bits = []
    if "boxed warning" in f:
        bits.append("boxed warning")
    if "hand-read label warnings" in f:
        bits.append("label warning")
    m = re.search(r"needs manual verification: HP:\S+ (.+?) '", f)
    if m:
        bits.append("verify: " + m.group(1))
    return "; ".join(bits) or "none flagged"


def _age(s: dict[str, Any]) -> str:
    return "/".join(AGE.get(s["age_bands"][b], "?") for b in BANDS)


def _supported(s: dict[str, Any]) -> bool:
    return str(s.get("age_union", "")).startswith("supported")


def _tier(s: dict[str, Any]) -> str:
    x = s["evidence_tier"].split(";")[0]
    x = re.sub(r"level (\d) literature \((.*?) papers\)", r"L\1 (\2)", x)
    return x.replace("no verified literature entry", "no entry").replace("literature without a supportive finding", "no supportive finding").replace(" papers", "")


def _native(s: dict[str, Any]) -> str:
    x = s["native"].replace("rule group ", "").replace("_", " ")
    for pat, rep in ((r"best help level \d; ", ""), (r"engine rank (\d+) of (\d+), score (-?[\d.]+)", r"rank \1 of \2, score \3"), (r"fraction of random queries with the drug in their top 50", "random-query share"),
                     (r" \(pathogenic gene only\), two-record links", " only"), (r", \d+ cell lines", ""), (r"BH q", "q"), (r"reversal p [\d.]+, FDR ", ""), (r"strengthens checkpoint, ", "")):
        x = re.sub(pat, rep, x)
    return x


def _limit(line: str, l: dict[str, Any], out: Path) -> str:
    """The line's own limit; for C and E-P the text is rebuilt from the adjusted analysis and the array counts so that it cannot be out of date."""
    if line == "C":
        c = _load(out, "line_c_adjusted")
        if c:
            return f"correlational, from cancer cell lines; adjusted for lineage and proliferation below ({c['n_hits_surviving_excluding_heme_lymphoid']} of {c['n_original_hits']} original hits survive when haematopoietic and lymphoid lines are excluded)"
    if line == "EP":
        g = (_load(out, "line_e_signature_reversal") or {}).get("geo")
        if g:
            return f"one adult patient lymphoblastoid cell line, {g['n_cases']} case arrays against {g['n_controls']} control arrays"
    return l["limit"]


def _line_table(line: str, l: dict[str, Any], out: Path) -> list[str]:
    adj = (_load(out, "line_c_adjusted") or {}) if line == "C" else {}
    head = "| Drug | Native score | Evidence | Age " + "/".join(BANDS) + " | Safety flag |" + (" Adjusted q (survives) |" if adj else "")
    t = [head, "|---|---|---|---|---|" + ("---|" if adj else "")]
    for s in l["shown"]:
        extra = ""
        if adj:
            q = (adj["drugs"].get(s["drug"]) or {}).get("adjusted", {}).get("q_bh")
            extra = f" {q:.2f} ({'yes' if s['drug'] in adj['surviving_hits'] else 'no'}) |" if q is not None else " n/a |"
        t.append(f"| {_name(s['drug'])} | {_native(s)} | {_tier(s)} | {_age(s)} | {_flag(s['safety_flags'])} |{extra}")
    return t


def _ln(ls: list[str]) -> str:
    return " and ".join("E-P" if x == "EP" else x for x in ls)


def _pairs(d: dict[str, Any]) -> list[tuple[str, list[str]]]:
    m = d["convergence"]["matrix"]
    return [(k, [x for x in ORDER if m[k].get(x)]) for k in d["convergence"]["in_two_or_more_lines"]]


def _join(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _line_a_ranks(out: Path) -> str:
    pr = _load(out, "proposals_review2")
    n = None
    if pr:
        n = next((x["signals"]["line_A_scored_drugs"] for x in pr["proposed"] if x["signals"].get("line_A_scored_drugs")), None)
    m = _load(out, "multi_evidence_mva_completed")
    pool = m["pool"]["size"] if m else None
    return ("Line A ranks have three sources: the graph rank above (combined rank, `out/kg_query_e5.json`); the position of the nine literature drugs among "
            + (f"the {n} network-scored drugs" if n else "the network-scored drugs") + " in Supplement S4 (`out/proposals_review2.json`); and the cross-check grade, a rank within the "
            + (f"{pool}-drug pool" if pool else "pool") + " damped by the random-query share. They are not interchangeable.")


def main_sections(out: Path) -> str:
    out = Path(out)
    if not _have(out):
        return ""
    d = _load(out, "candidates_by_line_mva_v2")
    lines, filt = d["lines"], d["filters"]
    c = _load(out, "line_c_adjusted")
    e = _load(out, "line_e_signature_reversal")
    t = ["### 3.13 Candidates proposed by each line\n"]
    t.append("The track asks for candidate medications as hypotheses for follow-up. No single line is the truth: each line proposes its own candidates from its own score, with its evidence (literature level, or no entry), the label reading for three age bands "
             "(0-4, 5-11 and 12-17 years: L labelled, P partly, N not labelled, A adults only, S age not stated, V verify; the age of the child is unknown, so a drug counts as label-supported if any band is) and the safety flags of the label screen. "
             f"Pool: {filt['pool_size']} approved small molecules; {len(filt['hard_excluded'])} hard exclusions (hazard drugs, growth-hormone class); {filt['cytotoxic_context_n']} cytotoxic agents are context or hazard and not proposed; "
             f"{len(filt['approval_unconfirmed_28'])} candidates with unconfirmed approval are listed in `out/candidates_by_line_mva_v2.md`. Convergence is information, not a requirement.\n")
    empty = []
    for ln in ORDER:
        l = lines[ln]
        if not l["shown"]:
            empty.append(f"{SHORT[ln]} ({l['n_native_set']} in its native set; {_limit(ln, l, out)})")
            continue
        t.append(f"**{SHORT[ln]}.** {l['n_native_set']} in the native set, {l['n_qualifying']} qualify; limit: {_limit(ln, l, out)}.\n")
        t.extend(_line_table(ln, l, out))
        also = [_name(x) for x in l.get("also_qualifying", [])]
        t.append("")
        if also:
            t.append(f"Also qualifying: {', '.join(also)}.\n")
        if ln == "C" and c:
            extra = [_name(x) for x in c["surviving_hits"] if x not in {s["drug"] for s in l["shown"]}]
            t.append(f"{len(c['surviving_hits'])} of {c['n_original_hits']} original hits survive the adjustment; {_join(extra)} also survive but are cytotoxic or caution-list classes, listed as context or hazard. "
                     f"Whether a positive correlation means vulnerability or resistance is unresolved (check CF1, Section 3.15).\n")
        if ln == "A":
            t.append(_line_a_ranks(out) + "\n")
        if ln == "D" and all(x["native"].startswith("strengthens checkpoint") for x in l["shown"]):
            t.append("All are called as strengthening the checkpoint, an assumption that is not justified for a hypomorphic allele.\n")
    t.append("Lines with no candidate: " + "; ".join(empty) + ". " + _near_b1(d) + "\n")
    if e:
        t.append(_e_note(e))
    t.append(_convergence(d, c, e))
    return "\n".join(t) + "\n"


def _near_b1(d: dict[str, Any]) -> str:
    near = [x for x in d.get("b1_near_threshold_information_only", []) if x["status"] == "candidate"]
    if not near:
        return ""
    top = max(near, key=lambda x: x["robustness_count"])
    return f"The closest B1 drug, information only, was {_name(top['drug'])} ({top['native']})."


def _e_note(e: dict[str, Any]) -> str:
    pc = e["positive_controls"]["results"]
    rec = [g for g, v in pc.items() if v["recovered_top5pct"]]
    nb = e["analysis"]["N=150"]["gene:BUB1B"]["n_usable"]
    p = e["P_signature"]
    g = e["geo"]
    return (f"Line E recovers {len(rec)} of {len(pc)} positive controls; the BUB1B knockout signature gave {nb} usable drugs and the combined Track 1 signature supports none. The patient signature (GSE22206, {g['n_cases']} case and {g['n_controls']} control arrays) "
            f"gave {p['n_usable']} usable drugs and {p['n_hazard']} concordant hazard flags, which are warnings to test, not refutations.\n")


def _convergence(d: dict[str, Any], c: dict[str, Any] | None, e: dict[str, Any] | None) -> str:
    m = d["convergence"]["matrix"]
    pairs = _pairs(d)
    allshown = [s for l in d["lines"].values() for s in l["shown"]]
    sup = {s["drug"] for s in allshown if _supported(s)}
    levels = [s.get("best_help_level") for s in d["lines"]["LIT"]["shown"] if s.get("best_help_level") is not None]
    best = max(levels) if levels else None
    parts = []
    for drug, ls in pairs:
        why = ""
        if c and "C" in ls:
            why = (": its C support survives the lineage and proliferation adjustment" if drug in c["surviving_hits"] else f": its C support does not survive the lineage and proliferation adjustment (adjusted q {c['drugs'][drug]['adjusted']['q_bh']:.2f})")
        if e and drug in {str(h).upper() for h in e["P_signature"]["hazards"]}:
            why += (" and " if why else ": ") + "it is also a concordant hazard for the patient signature"
        parts.append(f"{_name(drug)} ({_ln(ls)}){why}")
    return ("**Convergence matrix.** " + f"{len(m)} distinct drugs appear in the tables above; {len(pairs)} of them appear in two lines: " + ". ".join(parts) + ". "
            "Each pair includes a weak line (cell-line correlation, one patient cell line, AI-assigned literature labels): no drug has independent confirmation by two strong lines. "
            + (f"The best evidence is level {best} (related systems). " if best is not None else "")
            + f"{len(sup)} of the {len({s['drug'] for s in allshown})} drugs in the tables have label support in at least one age band" + (", so most are not labelled for children. " if len(sup) * 2 < len({s["drug"] for s in allshown}) else ". ")
            + "What is proposed is therefore a list per line to test, not a ranking and not a treatment: no claim of efficacy, no dose, no recommendation.\n")


def cross_check(out: Path) -> str:
    out = Path(out)
    if not _have(out):
        return ""
    d = _load(out, "candidates_by_line_mva_v2")
    o, k = _load(out, "multi_evidence_mva"), _load(out, "multi_evidence_mva_completed")
    cc = d["cross_check"]
    po, pk = o["proposed"][0], k["proposed"][0]
    named = d["changes_vs_first_output"]["named_drugs"].get(pk, {})
    kd = k["drugs"][pk]
    t = ["### 3.14 Multi-parameter cross-check\n"]
    t.append(f"A stochastic multicriteria analysis (SMAA, {o['n_draws']:,} weight draws, {len(o['criteria'])} criteria) combines all lines and the label and mechanism criteria into the probability that each drug is in the top five. "
             "It is a cross-check of the per-line tables and not one more line: it needs weights and a rule for missing criteria, and it penalises a missing criterion at "
             f"{o['missing_penalty']}.\n")
    t.append("| Run | Pool | Proposed (P in top 5) | Trametinib (P in top 5) |\n|---|---|---|---|")
    t.append(f"| Original run | {o['pool']['size']} | {_name(po)} ({o['drugs'][po]['p_top5']:.2f}) | {o['drugs']['TRAMETINIB']['p_top5']:.2f} |")
    t.append(f"| Completed with {k['completion']['n_asked']} more drugs ({k['completion']['n_confirmed']} confirmed approved) | {k['pool']['size']} | {_name(pk)} ({kd['p_top5']:.2f}) | {k['drugs']['TRAMETINIB']['p_top5']:.2f} |")
    tram = next((ls for dr, ls in _pairs(d) if dr == "TRAMETINIB"), [])
    t.append(f"| Candidates by line (Section 3.13) | {d['filters']['pool_size']} | no weights; one list per line | listed by {_ln(tram) if tram else 'no two lines'} |")
    t.append("")
    t.append(f"Trametinib is proposed by the original run and falls below the rule in the completed run, where {_name(pk)} is proposed instead. We read this as an artefact: {_name(pk)} has {kd['n_missing']} of {len(k['criteria'])} criteria missing and filled by renormalisation, "
             f"no verified literature entry, a confidence factor of {kd['confidence_factor']:.2f} and no age-band assessment"
             + (", and the label screen of Section 3.13 hard-excludes it for a contraindication tied to the child's nephrocalcinosis, which the SMAA run did not apply" if named.get("screen_result") == "hard_excluded" else "")
             + f". Its stability when each of the {len(cc['stability']['leave_one_out_proposed'])} criteria is left out does not test the missing-data rule, so it is not evidence. "
             f"With literature-only weights the run proposes {_join([_name(x) for x in cc['stability']['literature_only_proposed']])}, which include cytotoxic agents: the support for a drug depends on which lines are counted.\n")
    return "\n".join(t) + "\n"


def chain_intro(out: Path) -> str:
    out = Path(out)
    g = _load(out.parent / "config", "multi_evidence_track1_gene_block")
    if not g:
        return ""
    by = {}
    for x in g["genes"]:
        by.setdefault(x["tier"], []).append("*" + x["gene"] + "*")
    mod = next(x for x in g["genes"] if x["tier"].startswith("hypothesis-grade"))
    undet = by.get("undetermined", [])
    return (f"**The chain.** The primary gene is {_join(by['causal candidate'])} (a stop-gain with predicted decay plus a missense variant of uncertain significance, phase unresolved; Sections 3.9 and 3.10). "
            f"The other Track 1 genes are hypotheses, not findings: {_join(by[mod['tier']])} ({mod['tier']}; PP2A subunit linked in the literature to AKT and mTOR signalling; function tests of the variant; weight capped at {g['cap']}), "
            f"{_join(by['pathway context only'])} (context only) and the megabase-scale events at {_join(undet)} (undetermined: copy-number change or caller artefact; karyotype or single-cell sequencing decides). "
            "After the mechanism below come the candidates of each line (3.2 to 3.4 and 3.13), the cross-check (3.14) and what computation can refute and what needs the laboratory (3.15).\n")


def abstract_sentence(out: Path) -> str:
    out = Path(out)
    if not _have(out):
        return ""
    d = _load(out, "candidates_by_line_mva_v2")
    o, k = _load(out, "multi_evidence_mva"), _load(out, "multi_evidence_mva_completed")
    pairs = _pairs(d)
    n = {ln: d["lines"][ln]["n_qualifying"] for ln in ORDER}
    return (f"The lines proposed {n['LIT']} literature, {n['ENG']} engine, {n['A']} graph, {n['B2']} proximity, {n['C']} cell-line, {n['D']} direction, {n['EP']} patient signature candidates, and none from B1, E or the combined Track 1 signature; "
            f"{pairs_text(out)} are the only drugs in two lines, each pair including a weak line (Section 3.13). "
            f"The multi-parameter cross-check proposed {o['proposed'][0].lower()} in the original run and {k['proposed'][0].lower()} in the completed run, which we attribute to missing-data renormalisation and an unexcluded label hazard (Section 3.14).")


def brief_sentence(out: Path) -> str:
    out = Path(out)
    if not _have(out):
        return ""
    d = _load(out, "candidates_by_line_mva_v2")
    pairs = _pairs(d)
    return (f"Each line proposes its own candidates as hypotheses (Section 3.13); {pairs_text(out)} appear in two lines, always with a weak line, and the best evidence is level 1. "
            "Computation could refute some hypotheses (done: the adjustment of line C; defined and not run: CF1 to CF9) and the rest need a laboratory experiment, beginning with a karyotype or single-cell sequencing (Section 3.15).")


def engine_scalability(out: Path) -> str:
    eo = Path(out) / "engine"
    if not eng._have(eo):
        return ""
    v1, v4 = eng._load(eo, "validation"), eng._load(eo, "validation_v4")
    a = eng._g(v1, "main/primary_gold1/arms/A/auc_minus_popularity")
    b = eng._g(v4, "heldout/evaluation_group_bootstrap/models/B0/auc_minus_popularity")
    return (f"The empirical evidence for scalability and rigor is the ranking engine of Section 3.12: frozen, it beat a popularity baseline in AUC by {a['estimate']:.2f} on 17 pre-registered rare diseases and by {b['estimate']:.2f} on {eng._g(v4, 'heldout/n_diseases')} sealed held-out diseases "
            "and ran without refit on the diseases of this report; its advantage at low literature volume comes from the graph prior. Each line is a separate function with its own native score, so a line can be added or removed without touching the others.")


def criteria(out: Path) -> str:
    out = Path(out)
    if not _have(out):
        return ""
    d = _load(out, "candidates_by_line_mva_v2")
    c = _load(out, "line_c_adjusted")
    q = _load(out, "quote_check")
    ep = _load(out, "experiment_plan")
    surv = f"{c['n_hits_surviving_excluding_heme_lymphoid']} of {c['n_original_hits']}" if c else "some"
    t = ["### 4.7 How the report answers the four criteria of the track\n"]
    t.append("**Scientific rigor.** " + (f"Every one of {q['n_quotes']} label quotes is checked by code against its source; " if q else "") +
             f"the engine was validated on held-out diseases (Section 3.12); a candidate is shown failing (after adjustment for lineage and proliferation {surv} original line-C hits survive); weak and negative results are reported with the same weight.\n")
    t.append("**Potential impact.** Per-line hypotheses to test, "
             + (f"sorted into {ep['counts']['comp_first']} to check by computation first, {ep['counts']['lab_design_benefit']} for a laboratory design, {ep['counts']['hazard_class']} as hazards and {ep['counts']['closed']} closed, " if ep else "")
             + "with a caution list, a label table by age band and a diagnostic first step that bears on every hypothesis.\n")
    t.append("**Innovation.** No line is the truth: each proposes its own candidates, convergence is information, and the multi-criteria analysis is a cross-check whose failure mode is documented (Section 3.14).\n")
    t.append("**Scalability.** The engine validation of Sections 3.12 and 4.4 is the empirical evidence; each line is a separate function with its own score, so one can be added or removed alone.\n")
    return "\n".join(t) + "\n"


def pairs_text(out: Path) -> str:
    """'neratinib (C and EP) and trametinib (LIT and C)': the drugs that two lines propose."""
    d = _load(Path(out), "candidates_by_line_mva_v2")
    if not d:
        return ""
    return _join([p.lower() + " (" + _ln(ls) + ")" for p, ls in _pairs(d)])


def signature_reversal_scalability(out: Path) -> str:
    e = _load(Path(out), "line_e_signature_reversal")
    if not e:
        return ""
    pc = e["positive_controls"]["results"]
    rec = [g for g, v in pc.items() if v["recovered_top5pct"]]
    return (f"Signature reversal was run with public LINCS L1000 consensus knockout signatures (lines E and E-P, Section 3.13); it recovers {len(rec)} of {len(pc)} positive controls, "
            "and a knockout signature may not represent partial loss of function, so for other genes it is a warning line and not a source of candidates.")
