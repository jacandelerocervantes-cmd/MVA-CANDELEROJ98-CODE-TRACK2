"""Two small empirical analyses (read-only on existing outputs; writes only out/empirical/controls.json).

(1) Positive control of the caution layer: of the drugs with a verified report of harm or reduced sensitivity in MVA/BUB1B
    patients or BUB1B-deficient human cells, how many appear in the final caution output? No rule is tuned to the result.
(2) Catch rate of the quote verifier: counts/rates of verified / not_found / not_checkable per search round/topic, the
    not_found examples, and the label-quote check totals (read from out/quote_check.json, not re-run).

Strict positive set (fixed by hand, BEFORE looking at the caution output; every entry must be a sentence that
out/literature_hypotheses/verification.json marks 'verified', in a human patient with MVA/BUB1B disease or in human cells
with BUB1B deficiency, reporting harm or reduced sensitivity to a drug):
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "out" / "empirical" / "controls.json"
LIT = ROOT / "out" / "literature_hypotheses"

# (drug, pmid, kind, setting, approved_in_ChEMBL_phase4)
POSITIVE_SET: list[dict[str, Any]] = [
    {"drug": "DACTINOMYCIN", "pmid": "31081598", "kind": "harm",
     "setting": "child with biallelic BUB1B PCS/MVA; severe adverse reaction to dactinomycin stopped chemotherapy", "approved": True},
    {"drug": "TREOSULFAN", "pmid": "31053147", "kind": "harm",
     "setting": "MVA1 (BUB1B) boy; treosulfan-based conditioning HSCT followed by graft rejection, pancytopenia, sepsis, death", "approved": True},
    {"drug": "PACLITAXEL", "pmid": "39014450", "kind": "reduced_sensitivity",
     "setting": "human cells/lymphocytes with germline BUB1B variants; reduced paclitaxel sensitivity", "approved": True},
    {"drug": "PACLITAXEL", "pmid": "38987356", "kind": "reduced_sensitivity",
     "setting": "human PANC-1 cells; BUB1B knockdown attenuates paclitaxel-induced cell death", "approved": True},
    {"drug": "COLCEMID", "pmid": "10877982", "kind": "reduced_sensitivity",
     "setting": "fibroblasts of two PCS/MVA infants insensitive to colcemid (no mitotic block)", "approved": False},
]
# Near misses that were examined and left out, with the reason (documented, not tuned).
EXCLUDED_NEAR_MISSES = [
    {"drug": "CYCLOSPORINE", "pmid": "32884756", "why": "sentence is 'not_found' in verification.json (the verifier rejected it); not a verified report"},
    {"drug": "CORTICOTROPIN", "pmid": "23916859", "why": "no benefit (refractory spasms), not harm or reduced sensitivity to a cytotoxic/spindle drug"},
    {"drug": "SOMATROPIN", "pmid": "35434947", "why": "MVA type 2 (CEP57), not BUB1B; unsatisfactory response, not harm"},
    {"drug": "FLUDARABINE / ANTITHYMOCYTE GLOBULIN", "pmid": "31053147", "why": "only listed as conditioning components; no drug-specific harm stated"},
    {"drug": "VINCRISTINE", "pmid": "39698932", "why": "DESI1 knockdown, not BUB1B deficiency"},
    {"drug": "PACLITAXEL", "pmid": "17234259", "why": "ovarian line resistant to paclitaxel, not BUB1B-deficient (kept only as non-BUB1B evidence)"},
]


def _load(p: Path) -> Any:
    return json.loads(p.read_text(encoding="utf-8"))


def _norm(s: str) -> str:
    return " ".join(str(s).upper().split())


def _verified_pmids(rows: list[dict]) -> dict[tuple[str, str], list[str]]:
    d: dict[tuple[str, str], list[str]] = defaultdict(list)
    for r in rows:
        d[(_norm(r.get("drug") or ""), str(r["pmid"]))].append(r["status"])
    return d


def positive_control() -> dict[str, Any]:
    ver = _load(LIT / "verification.json")["rows"]
    vmap = _verified_pmids(ver)
    caution = _load(ROOT / "out" / "caution_list.json")
    cl_names: dict[str, str] = {}
    for cls, v in caution["clases"].items():
        for f in v["farmacos"]:
            cl_names[_norm(f["name"])] = cls
    sin_cita = {_norm(x["farmaco"]): x["motivo"] for x in caution["sin_cita"]}
    p2 = _load(ROOT / "out" / "proposals_review2.json")
    p2_caution = {_norm(c["drug"]): c["level"] for c in p2["caution"]}

    entries = []
    for e in POSITIVE_SET:
        d = _norm(e["drug"])
        statuses = vmap.get((d, e["pmid"]))
        if statuses is None:  # name variants such as 'COLCEMID (DEMECOLCINE)'
            statuses = [s for (dr, pm), ss in vmap.items() if pm == e["pmid"] and dr.startswith(d) for s in ss]
        in_cl = d in cl_names
        in_p2 = d in p2_caution
        if in_cl or in_p2:
            why = None
        elif not e["approved"]:
            why = "not an approved drug (outside data/universe.json, ChEMBL max_phase 4); the caution layer is approved-only by design (E14)"
        elif d in sin_cita:
            why = "class member without a readable label: " + sin_cita[d]
        else:
            why = "drug is not a member of any pre-specified caution class (tubulin, mitotic checkpoint, genotoxic chemotherapy, proteasome, GH)"
        entries.append({**e, "verified_in_verification_json": bool(statuses) and all(s == "verified" for s in statuses),
                        "in_caution_list_json": in_cl, "caution_list_class": cl_names.get(d),
                        "in_proposals_review2_caution": in_p2, "proposals_review2_level": p2_caution.get(d),
                        "in_final_caution_output": in_cl or in_p2, "why_not_in_final_output": why})
    drugs = sorted({_norm(e["drug"]) for e in entries})
    by_drug = {d: {"pmids": sorted({e["pmid"] for e in entries if _norm(e["drug"]) == d}),
                   "in_caution_list_json": any(e["in_caution_list_json"] for e in entries if _norm(e["drug"]) == d),
                   "in_proposals_review2_caution": any(e["in_proposals_review2_caution"] for e in entries if _norm(e["drug"]) == d)}
               for d in drugs}
    rec = sorted(d for d, v in by_drug.items() if v["in_caution_list_json"] or v["in_proposals_review2_caution"])
    rec_cl = sorted(d for d, v in by_drug.items() if v["in_caution_list_json"])
    approved = sorted({_norm(e["drug"]) for e in entries if e["approved"]})
    return {
        "definition": "drug with a sentence marked 'verified' in verification.json reporting harm or reduced sensitivity in MVA/BUB1B patients or BUB1B-deficient human cells (hand-fixed list, see module docstring)",
        "final_caution_output": "union of out/caution_list.json (class-based, label quotes) and the 'caution' block of out/proposals_review2.json (literature-based)",
        "positive_set_drugs": drugs, "M": len(drugs), "N_recovered": len(rec), "recovered": rec,
        "not_recovered": sorted(set(drugs) - set(rec)),
        "recovers": f"recovers {len(rec)} of {len(drugs)}",
        "approved_only": {"M": len(approved), "N_recovered": len([d for d in approved if d in rec]),
                          "recovers": f"recovers {len([d for d in approved if d in rec])} of {len(approved)}"},
        "caution_list_json_only": {"recovered": rec_cl, "N_recovered": len(rec_cl),
                                   "note": "independent of the literature: classes were pre-specified, so this is the fair positive control"},
        "proposals_review2_caution_is_circular": "that block is built from these same verified reports, so recovery there is expected by construction",
        "entries": entries, "excluded_near_misses": EXCLUDED_NEAR_MISSES,
    }


def verifier_catch_rate() -> dict[str, Any]:
    v = _load(LIT / "verification.json")
    rows = v["rows"]
    per: dict[str, Counter] = defaultdict(Counter)
    tot: Counter = Counter()
    for r in rows:
        per[r["topic"]][r["status"]] += 1
        tot[r["status"]] += 1

    def block(c: Counter) -> dict[str, Any]:
        n = sum(c.values())
        out: dict[str, Any] = {"n": n}
        for s in ("verified", "not_found", "not_checkable"):
            out[s] = c.get(s, 0)
            out[s + "_rate"] = round(c.get(s, 0) / n, 4) if n else None
        return out

    nf = [{"drug": r["drug"] or None, "pmid": r["pmid"], "topic": r["topic"], "section": r["section"],
           "claimed_verified_by_ai": r.get("claimed_verified"), "sentence_first_100": r["sentence"][:100]}
          for r in rows if r["status"] == "not_found"]
    nc = [{"drug": r["drug"] or None, "pmid": r["pmid"] or None, "topic": r["topic"], "section": r["section"],
           "sentence_first_100": r["sentence"][:100]} for r in rows if r["status"] == "not_checkable"]
    qc = _load(ROOT / "out" / "quote_check.json")
    return {
        "file_counts_field": v["counts"], "n_sentences": v["n_sentences"],
        "total": block(tot), "per_topic": {t: block(c) for t, c in sorted(per.items())},
        "distinct_not_found_sentences": len({(x["pmid"], x["sentence_first_100"]) for x in nf}),
        "not_found_entries": nf, "not_checkable_entries": nc,
        "caveat": "rows are sentence slots; the same sentence can occupy several sections/topics, so rows overcount distinct sentences",
        "label_quote_check": {"n_quotes": qc["n_quotes"], "found": qc["found"], "not_found": qc["not_found"],
                              "not_checkable": qc["not_checkable"], "by_origin": qc["by_origin"],
                              "source": "out/quote_check.json (not re-run)"},
    }


def run() -> dict[str, Any]:
    res = {"positive_control": positive_control(), "verifier_catch_rate": verifier_catch_rate()}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    return res


if __name__ == "__main__":
    r = run()
    print(r["positive_control"]["recovers"], r["positive_control"]["recovered"], r["positive_control"]["not_recovered"])
    print(r["verifier_catch_rate"]["total"])
