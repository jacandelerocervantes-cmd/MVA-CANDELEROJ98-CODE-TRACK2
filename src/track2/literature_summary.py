"""Mechanism-guided proposals from the verified literature search (addendum 11, E16: selection rules declared after the results of the lines were known).

    python -m track2.literature_summary [--out out]

Input: the four searches `out/literature_hypotheses/*.json` and their independent re-check `verification.json` (`track2.verify_literature`). Only sentences that the re-check confirmed are used. Rules:

R1 The drug is approved (it is in `data/universe.json`, ChEMBL max_phase 4).
R2 *Constitutional hypothesis*: at least one verified finding with direction `may_help` and evidence level >= 1 on a feature of the child (growth and muscle, proteotoxic or energy stress, senescence), and no verified level-3
   harm. *Tumour-context hypothesis*: a verified finding at level >= 1 with direction `may_help` or `selectively_kills_aneuploid_cells` in an aneuploidy or tumour system; for a drug that kills aneuploid cells selectively
   the same property is recorded as a hazard for the constitutional mosaic. A drug with a verified `may_harm` at level >= 1 is marked as having mixed evidence.
R3 A *class sibling* is listed with its class representative, never alone, and keeps its own (lower) level.
R4 *Caution*: a verified harm at level 3 (a case in MVA or BUB1B), a verified `may_harm`, or a finding flagged `use=caution`, whatever the level.
No drug reaches level 2 or 3 as a candidate in the literature found. Levels: L0 inferred or indirect; L1 tested in a related system; L2 tested in BUB1B/BUBR1-deficient or MVA-derived human cells or a hypomorphic model with rescue; L3 clinical observation in MVA or BUB1B.
Writes `out/proposals_review2.json`. The experiments are designs proposed here for each hypothesis, not results.
"""

from __future__ import annotations

import argparse
import re
import json
from pathlib import Path
from typing import Any

TOPICS = ["growth_muscle", "proteotoxic_autophagy", "senescence_p53", "tumour_cin_safety", "round2_models", "round2_clinical", "round2_phenotypes", "round2_vulnerabilities", "round3_audit"]
CONSTITUTIONAL = {"growth_muscle", "proteotoxic_energy", "senescence_p53"}
SIBLINGS = {"EVEROLIMUS": "SIROLIMUS", "SELUMETINIB": "TRAMETINIB"}  # same target class; the sibling has only weaker evidence of its own
LINK = {
    "SIROLIMUS": "Consequence 4 (growth restriction and muscle atrophy) and consequence 2 (aneuploidy-related stress): mTORC1 hyperactivity correlates with sarcopenia in mice with BubR1 alleles, and mTOR inhibition rescued the aneuploid neuroblasts of a Drosophila model of spindle-checkpoint chromosomal instability.",
    "METFORMIN": "Consequence 2 (aneuploidy-related stress): it lowered oxidative damage in cell lines with constitutional trisomy 21; other work reports that it changes BUBR1 and MAD2 expression in tumour cells, so the direction is mixed.",
    "TRAMETINIB": "Consequence 5 (tumour predisposition), tumour context only: highly aneuploid, non-transformed human clones were more sensitive to the MEK inhibitor trametinib; the same property is a hazard for normal aneuploid cells.",
    "PHENYLBUTANOIC ACID": "Consequences 2 and 3 (proteotoxic stress and senescence): the chemical chaperone 4-phenylbutyrate prevented premature senescence in fibroblasts from trisomies 21, 18 and 13 and reduced protein aggregates and apoptosis in trisomic neurons; these are whole-chromosome trisomies, not mosaic aneuploidy.",
    "SELINEXOR": "Consequence 5 (tumour predisposition), tumour context only: the XPO1 inhibitor suppressed TRIP13, which Wilms tumour cell lines need to survive; TRIP13 acts at the spindle checkpoint, and selinexor is also in the BUB1B-module (B2) set of the computational lines.",
    "BOSUTINIB": "Consequence 5 (tumour predisposition), tumour context and a hazard: cells with an ongoing chromosomal instability caused by weakening the spindle checkpoint were more sensitive to the Src inhibitor SKI606 (the code name of bosutinib; the abstract uses the code name, and a separate paper, PMID 39924517, states \"bosutinib (SKI-606)\"), the same selectivity that would threaten normal aneuploid cells.",
    "BORTEZOMIB": "Consequence 2 (proteotoxic stress), tumour context only: aneuploid human cells were more sensitive to the proteasome inhibitor, which is the same property that makes it a hazard for constitutionally aneuploid normal cells.",
    "CHLOROQUINE": "Consequence 2 (proteotoxic and autophagic stress), tumour context only: it was lethal to trisomic mouse fibroblasts and not to euploid ones, which is also a hazard for normal aneuploid cells of the mosaic.",
}
EXPERIMENT = {
    "SIROLIMUS": "In BUB1B-hypomorphic human cells (patient-derived fibroblasts or myoblasts if they exist, otherwise a hypomorphic BUB1B allele made in a myoblast or fibroblast line): phospho-S6 as a measure of mTORC1, aneuploidy by single-cell sequencing or chromosome spreads, senescence (SA-beta-gal, p21), proteotoxic markers and myotube differentiation, with the mTOR inhibitor against vehicle and a second mTOR inhibitor as a control. It is refuted if aneuploidy or senescence rise, or if no readout improves.",
    "METFORMIN": "The same readouts, plus AMPK activity and oxidative damage. It is refuted if BUBR1 or MAD2 fall further or if aneuploidy rises.",
    "TRAMETINIB": "In tumour cells of the child, if available, or in rhabdomyosarcoma lines with graded aneuploidy and spindle-checkpoint strength: sensitivity against aneuploidy score, and the viability of the patient's non-tumour cells in parallel. It is refuted if sensitivity does not track aneuploidy or if non-tumour cells are as sensitive.",
    "PHENYLBUTANOIC ACID": "In patient-derived BUB1B-hypomorphic or MVA fibroblasts against a euploid control: senescence (SA-beta-gal, p21), protein aggregation and aneuploidy, with the compound against vehicle. It is refuted if senescence or aneuploidy rise, or if no readout improves.",
    "SELINEXOR": "In rhabdomyosarcoma and Wilms tumour lines graded for TRIP13 and BUB1B expression: viability against TRIP13 level, and the viability of non-tumour cells in parallel. It is refuted if sensitivity does not track TRIP13 or if non-tumour cells are as sensitive.",
    "BOSUTINIB": "Not a treatment hypothesis for the constitutional state. As a hazard check: its toxicity on the patient's non-tumour fibroblasts against a euploid control; as a tumour hypothesis: sensitivity of tumour lines against checkpoint strength.",
    "BORTEZOMIB": "Not a treatment hypothesis for the constitutional state. As a hazard check: its toxicity on the patient's non-tumour fibroblasts against a euploid control; as a tumour hypothesis: sensitivity of tumour lines against aneuploidy score.",
    "CHLOROQUINE": "Not a treatment hypothesis for the constitutional state. As a hazard check: its toxicity on the patient's non-tumour fibroblasts against a euploid control; as a tumour hypothesis: sensitivity of tumour lines against aneuploidy score.",
}


ALIAS = {"SODIUM PHENYLBUTYRATE": "PHENYLBUTANOIC ACID", "NICOTINAMIDE": "NIACINAMIDE", "KPT-330": "SELINEXOR", "SKI606": "BOSUTINIB", "RAPAMYCIN": "SIROLIMUS"}  # names used by the searches -> name in data/universe.json


EDITORIAL = {
    "SELINEXOR": "one computational line (B2) plus a literature finding, which E9 never accepts as a candidate; the finding is in Wilms tumour and not rhabdomyosarcoma, it ranks 1795 of 2128 in line A, and its label has no paediatric use",  # drugs that pass the mechanical rules but are not proposed; the reason is shown in the report (E16, round 2)
    "DENOSUMAB": "only a murine anti-RANKL antibody in progeroid (Zmpste24) mice, with a bone and muscle readout; no link to the spindle checkpoint or to a feature of the child",
    "SOMATROPIN": "growth hormone: the label contraindicates active malignancy and the child's record documents rhabdomyosarcoma (current tumour status not available); the evidence is case reports in other syndromes (Meier-Gorlin, Seckel) and one unsatisfactory response in MVA type 2; a meta-analysis points to a second-tumour risk",
    "MECASERMIN": "IGF-1 axis with the same tumour concern as growth hormone; the evidence is single cases in other syndromes, one of them with IGF-1 resistance",
    "NIACINAMIDE": "no approved systemic drug label (it is sold as a vitamin or supplement); the one level-1 finding is a Drosophila model, and a mouse study found that nicotinamide supplementation conferred resistance to chemotherapy in tumours",
    "ERIBULIN": "a tubulin binder: BUB1B variants and a weak checkpoint are linked to reduced sensitivity to taxanes (level 2) and cells from MVA infants do not arrest in mitosis with spindle poisons",
}


def norm_consequence(topic: str, consequence: str | None) -> str | None:
    """Round-2 searches used free-text consequence tags; they are mapped here by keyword to the tags of the rules (the mapping is part of E16)."""
    if not topic.startswith("round2"):
        return consequence
    c = (consequence or "").lower()
    if "senescen" in c:
        return "senescence_p53"
    if re.search(r"aggregat|chaperone|proteo|nad\+|atp|energy|oxidative|antioxid|autophag", c):
        return "proteotoxic_energy"
    if re.search(r"growth|muscle|sarcopen|myogen|height|stature", c):
        return "growth_muscle"
    if re.search(r"kill|sensitiv|tumou?r|cancer", c):
        return "tumour"
    return "other_context"


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def collect(directory: Path) -> dict[str, dict[str, list[dict[str, Any]]]]:
    status = {(r["topic"], r["section"], r["index"]): r["status"] for r in _load(directory / "verification.json")["rows"]}
    by_drug: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for topic in TOPICS:
        data = _load(directory / f"{topic}.json")
        for section in ("findings", "against_or_risk"):
            for i, x in enumerate(data.get(section, [])):
                if not x.get("drug") or status.get((topic, section, i)) != "verified":
                    continue
                row = {"topic": topic, "section": section, "level": int(x.get("evidence_level") or 0), "direction": x.get("direction"), "consequence": norm_consequence(topic, x.get("consequence")), "use": x.get("use"),
                       "system": x.get("system"), "finding": x.get("finding"), "sentence": x.get("sentence"), "pmid": x.get("pmid"), "year": x.get("year"), "caveats": x.get("caveats")}
                by_drug.setdefault(ALIAS.get(x["drug"].upper(), x["drug"].upper()), {"entries": []})["entries"].append(row)
    return by_drug


def classify(rows: list[dict[str, Any]]) -> dict[str, Any]:
    helps = [r for r in rows if r["section"] == "findings" and r["direction"] == "may_help"]
    kills = [r for r in rows if r["direction"] == "selectively_kills_aneuploid_cells"]
    harms = [r for r in rows if r["direction"] == "may_harm"]
    l3_harm = [r for r in rows if r["level"] >= 3 and r["direction"] in ("may_harm",)]
    constitutional = [r for r in helps if r["level"] >= 1 and r["consequence"] in CONSTITUTIONAL]
    # selective killing of aneuploid cells is by definition a tumour-context property (and a hazard for the constitutional mosaic), whatever the mechanism that was studied
    tumour = [r for r in helps + kills if r["level"] >= 1 and (r["direction"] == "selectively_kills_aneuploid_cells" or r["consequence"] not in CONSTITUTIONAL | {"other_context"})]
    caution = bool(l3_harm or any(r["use"] == "caution" for r in rows) or any(r["direction"] == "may_harm" and r["level"] >= 3 for r in rows))
    return {"constitutional": constitutional, "tumour": tumour, "kills": kills, "harms": harms, "l3_harm": l3_harm, "caution": caution,
            "best_help_level": max([r["level"] for r in helps] or [-1])}


def net_signals(name: str, out: Path, universe: dict[str, Any]) -> dict[str, Any]:
    q = _load(out / "kg_query_e5.json")
    names = {d["drug_id"]: d["name"].upper() for d in universe["drugs"]}
    best: dict[str, int] = {}
    for drug, r in q["rank_combined"].items():
        best[drug] = min(r, best.get(drug, r))
    order = sorted(best, key=lambda k: (best[k], k))
    rank = {names[d]: i + 1 for i, d in enumerate(order) if d in names}
    nq = _load(out / "kg_nullquery_e5.json")["fraction_of_random_queries_with_drug_in_top50"]
    frac = {n.upper(): f for n, f in nq.items()}
    fl = _load(out / "four_lines_review2.json")
    vu = _load(out / "vulnerability_e6.json")
    d_calls = _load(out / "signed_effect.json")["drugs"]
    by_name = {d["name"].upper(): d for d in universe["drugs"]}
    d = by_name.get(name)
    return {"line_A_rank": rank.get(name), "line_A_scored_drugs": len(rank), "random_query_fraction_in_top50": frac.get(name), "in_B2_robust_set": name in set(fl["support_detail"]["B2_W3_two_record"]),
            "in_line_C_hits": name in {v["name"].upper() for v in vu["supported_drugs"].values()}, "line_D_note_exploratory": (d_calls.get(d["drug_id"]) or {}).get("call", "unclear") if d else "unclear",
            "has_curated_target": bool(d and d.get("human_targets"))}


def label_status(name: str, config: Path) -> dict[str, Any] | None:
    for e in _load(config)["entries"].values():
        if e["name"] == name and "2026-10-08" in e.get("source", ""):
            return {k: e.get(k) for k in ("verdict", "min_age_years", "scope", "route", "source", "quote", "also_read")}
    return None


def best_entries(rows: list[dict[str, Any]], pick: Any, k: int = 4) -> list[dict[str, Any]]:
    chosen = sorted([r for r in rows if pick(r)], key=lambda r: (-r["level"], str(r["pmid"])))
    seen, out = set(), []
    for r in chosen:
        key = (r["pmid"], r["sentence"])
        if key not in seen:
            seen.add(key)
            out.append({x: r[x] for x in ("level", "direction", "consequence", "system", "sentence", "pmid", "year", "caveats")})
    return out[:k]


def build(out: Path, config: Path, universe_path: Path) -> dict[str, Any]:
    universe = _load(universe_path)
    approved = {d["name"].upper() for d in universe["drugs"]}
    data = collect(out / "literature_hypotheses")
    info = {n: classify(v["entries"]) for n, v in data.items() if n in approved}
    proposed, caution, considered = [], [], []
    for name, c in sorted(info.items()):
        rows = data[name]["entries"]
        reps = SIBLINGS.get(name)
        group = "constitutional" if c["constitutional"] else ("tumour_context" if c["tumour"] else None)
        if name in EDITORIAL:
            group = None
            considered.append({"drug": name, "why_not": "editorial exclusion: " + EDITORIAL[name], "levels": sorted({r["level"] for r in rows})})
        if c["caution"] and name not in ("BORTEZOMIB", "CHLOROQUINE"):
            group = None  # a caution flag (level-3 harm or flagged by the search) takes precedence over a proposal
        if name in SIBLINGS:
            group = "class_sibling"
        if group:
            help_rows = best_entries(rows, lambda r: r["direction"] in ("may_help", "selectively_kills_aneuploid_cells") and r["level"] >= (0 if name in SIBLINGS else 1))
            against = best_entries(rows, lambda r: r["direction"] in ("may_harm", "no_effect") or r["section"] == "against_or_risk")
            mixed = bool(c["harms"]) and group == "constitutional"
            proposed.append({"drug": name, "group": group, "sibling_of": reps, "best_level": max([e["level"] for e in help_rows] or [0]), "mixed_evidence": mixed,
                             "hazard_for_constitutional_mosaic": bool(c["kills"]), "link": LINK.get(name) or (LINK.get(reps) if reps else None), "supporting": help_rows, "against_or_risk": against,
                             "signals": net_signals(name, out, universe), "label": label_status(name, config), "experiment": EXPERIMENT.get(name) or EXPERIMENT.get(reps or "", None)})
        if c["caution"] or (c["harms"] and name not in {p["drug"] for p in proposed}):
            caution.append({"drug": name, "level": max([r["level"] for r in rows if r["direction"] == "may_harm" or r["use"] == "caution"] or [0]),
                            "evidence": best_entries(rows, lambda r: r["direction"] == "may_harm" or r["use"] == "caution")})
        if not group and not c["caution"] and name not in EDITORIAL:
            considered.append({"drug": name, "why_not": "no verified finding that meets the rules (only level 0, unclear or no-effect findings)", "levels": sorted({r["level"] for r in rows})})
    order = {"constitutional": 0, "class_sibling": 1, "tumour_context": 2}
    proposed.sort(key=lambda p: (order[p["group"]], -p["best_level"], p["drug"]))
    return {"rules": __doc__, "n_drugs_with_verified_entries": len(info), "proposed": proposed, "caution": caution, "considered_not_proposed": considered}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--config", type=Path, default=Path("config/pediatric_label_review.json"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    args = parser.parse_args()
    result = build(args.out, args.config, args.universe)
    (args.out / "proposals_review2.json").write_text(json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print("proposed:", [(p["drug"], p["group"], f"L{p['best_level']}", "mixed" if p["mixed_evidence"] else "", "hazard" if p["hazard_for_constitutional_mosaic"] else "") for p in result["proposed"]])
    print("caution:", [(c["drug"], f"L{c['level']}") for c in result["caution"]])
    print("considered, not proposed:", [c["drug"] for c in result["considered_not_proposed"]])


if __name__ == "__main__":
    main()
