"""CF9: independent GDSC2 drug-sensitivity check of the line C survivors (rule: protocol/addendum_computational_followups.md, section 'CF9 execution').

    PYTHONPATH=src python -I -m track2.cf9 [--root .] [--out out]    (reads a downloaded file: run with -I)

Offline: reads data/gdsc/gdsc2_r8_5, data/depmap and out/line_c_adjusted.json. Correlational, from cancer cell lines.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any

import numpy as np

from track2 import line_c_adjusted as lca
from track2.xlsx_min import read_rows

GDSC_FILE = "GDSC2_fitted_dose_response_27Oct23.xlsx"
SURVIVORS = ("LAZERTINIB", "NERATINIB", "PACLITAXEL", "VINCRISTINE", "COLCHICINE", "ERIBULIN")
BENEFIT = ("LAZERTINIB", "NERATINIB")
MIN_LINES = 100
P_CUT = 0.05


def norm(s: Any) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s or "").upper())


def match_drugs(gdsc_names: list[str], synonyms: dict[str, set[str]]) -> dict[str, list[str]]:
    """Survivor -> GDSC DRUG_NAMEs whose normalised name equals a normalised name or synonym of the survivor (exact, no fuzzy)."""
    by = {}
    for g in gdsc_names:
        by.setdefault(norm(g), []).append(g)
    return {d: sorted({g for k in {norm(d)} | {norm(s) for s in syn} for g in by.get(k, [])}) for d, syn in synonyms.items()}


def match_lines(gdsc: list[dict[str, Any]], models: list[dict[str, str]]) -> tuple[dict[str, str], dict[str, int]]:
    """GDSC line key (COSMIC_ID) -> DepMap ModelID: Sanger id, then COSMIC id, then name. Ambiguous matches are dropped."""
    idx = {"sanger": {}, "cosmic": {}, "name": {}}
    for m in models:
        idx["sanger"].setdefault(m["SangerModelID"], set()).add(m["ModelID"]) if m.get("SangerModelID") else None
        c = (m.get("COSMICID") or "").split(".")[0]
        idx["cosmic"].setdefault(c, set()).add(m["ModelID"]) if c else None
        for k in (norm(m.get("StrippedCellLineName")), norm(m.get("CellLineName"))):
            idx["name"].setdefault(k, set()).add(m["ModelID"]) if k else None
    out: dict[str, str] = {}
    cov = {"gdsc_lines": 0, "by_sanger_id": 0, "by_cosmic_id": 0, "by_name": 0, "unmatched": 0, "ambiguous_dropped": 0}
    for g in gdsc:
        cov["gdsc_lines"] += 1
        keys = [("by_sanger_id", "sanger", g["sanger"]), ("by_cosmic_id", "cosmic", str(int(g["cosmic"])) if g["cosmic"] else ""), ("by_name", "name", norm(g["name"]))]
        for label, table, key in keys:
            s = idx[table].get(key)
            if s:
                if len(s) == 1:
                    out[g["cosmic"]] = next(iter(s))
                    cov[label] += 1
                else:
                    cov["ambiguous_dropped"] += 1
                break
        else:
            cov["unmatched"] += 1
    # reverse ambiguity: two GDSC lines to one model
    seen: dict[str, list[str]] = {}
    for k, v in out.items():
        seen.setdefault(v, []).append(k)
    for v, ks in seen.items():
        if len(ks) > 1:
            for k in ks:
                del out[k]
            cov["ambiguous_dropped"] += len(ks)
    return out, cov


def verdict(adj: dict[str, Any] | None, orig_rho: float) -> str:
    if adj is None or not adj.get("estimable"):
        return "not estimable"
    return "survives" if (np.sign(adj["rho_adj"]) == np.sign(orig_rho) and adj["p_perm"] <= P_CUT) else "refuted"


def run_drug(y: np.ndarray, score: np.ndarray, prolif: np.ndarray, lineage: np.ndarray, min_lines: int = MIN_LINES, n_perm: int = lca.N_PERM) -> dict[str, Any]:
    ok = np.isfinite(y) & np.isfinite(score) & np.isfinite(prolif) & (lineage != "")
    if ok.sum() < min_lines:
        return {"estimable": False, "n_lines": int(ok.sum())}
    r = lca.partial_perm(y[ok], score[ok], prolif[ok], lineage[ok], n_perm, lca.SEED)
    r["estimable"] = True
    r["n_lineages"] = len(set(lca.pool_small(lineage[ok])))
    return r


def load_gdsc(path: Path) -> list[dict[str, Any]]:
    it = read_rows(str(path))
    h = next(it)
    return [dict(zip(h, r)) for r in it]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.out)
    from track2 import prism_l4

    dep = root / "data" / "depmap"
    rows = load_gdsc(root / "data" / "gdsc" / "gdsc2_r8_5" / GDSC_FILE)
    lc = json.loads((out / "line_c_adjusted.json").read_text(encoding="utf-8"))["drugs"]
    keys = prism_l4.compound_keys(dep / "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")
    syn = {d: {k for _, (nm, joined) in keys.items() if d in joined.split("|") for k in joined.split("|")} for d in SURVIVORS}
    names = sorted({r["DRUG_NAME"] for r in rows})
    dmatch = match_drugs(names, syn)
    with open(dep / "Model.csv", encoding="utf-8", newline="") as h:
        models = list(csv.DictReader(h))
    glines = {}
    for r in rows:
        glines[r["COSMIC_ID"]] = {"cosmic": r["COSMIC_ID"], "sanger": r["SANGER_MODEL_ID"] or "", "name": r["CELL_LINE_NAME"]}
    lmap, cov = match_lines(list(glines.values()), models)
    scores = prism_l4.sac_scores(dep / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", prolif_genes=lca.PROLIF_GENES)
    lin = prism_l4.lineages(dep / "Model.csv")
    usable = {c: m for c, m in lmap.items() if m in scores and lin.get(m)}
    cov.update({"matched_to_depmap_model": len(lmap), "retained_with_expression_and_lineage": len(usable)})
    res: dict[str, Any] = {"label": "exploratory, post hoc, correlational, cancer cell lines; no patient data", "protocol": "protocol/addendum_computational_followups.md (CF9 execution)",
                           "screen": "GDSC2 fitted dose response, release 8.5 (" + GDSC_FILE + ")", "n_gdsc_rows": len(rows), "n_gdsc_drugs": len(names), "line_matching_coverage": cov,
                           "min_lines": MIN_LINES, "p_cut": P_CUT, "drugs": {}}
    est: list[str] = []
    for d in SURVIVORS:
        orig = lc[d]["original_rho"]
        gm = dmatch[d]
        row: dict[str, Any] = {"gdsc_names_matched": gm, "original_rho": orig, "orientation": "LN_IC50 low = sensitive; rho > 0 = more lethal in weak-checkpoint lines"}
        if not gm:
            row.update({"adjusted": {"estimable": False, "n_lines": 0}, "verdict": "not estimable", "reason": "drug (name or PRISM synonym) absent from GDSC2"})
        else:
            by: dict[str, dict[str, list[float]]] = {}
            for r in rows:
                if r["DRUG_NAME"] in gm and r["COSMIC_ID"] in usable:
                    by.setdefault(usable[r["COSMIC_ID"]], {"ic50": [], "auc": []})
                    by[usable[r["COSMIC_ID"]]]["ic50"].append(r["LN_IC50"])
                    by[usable[r["COSMIC_ID"]]]["auc"].append(r["AUC"])
            ids = sorted(by)
            sc = np.array([scores[m]["mean_all"] for m in ids])
            pr = np.array([scores[m]["prolif"] for m in ids])
            ln = np.array([lin[m] for m in ids], dtype=object)
            row["adjusted"] = run_drug(np.array([np.mean(by[m]["ic50"]) for m in ids], dtype=float), sc, pr, ln)
            row["adjusted_auc_secondary"] = run_drug(np.array([np.mean(by[m]["auc"]) for m in ids], dtype=float), sc, pr, ln)
            if row["adjusted"]["estimable"]:
                est.append(d)
            row["verdict"] = verdict(row["adjusted"], orig)
        res["drugs"][d] = row
    if est:
        q = lca.bh(np.array([res["drugs"][d]["adjusted"]["p_perm"] for d in est]))
        for d, qv in zip(est, q):
            res["drugs"][d]["adjusted"]["q_bh"] = float(qv)
    res["n_estimable"] = len(est)
    res["cf9_overall"] = "not runnable" if not est else "runnable"
    for d in SURVIVORS:
        v = res["drugs"][d]["verdict"]
        res["drugs"][d]["plan_consequence"] = (
            ("hazard class: stays in the hazard panel and toxicity design H1 (positive/control check only); " if d not in BENEFIT else "")
            + {"survives": "stays a hypothesis; CF9 passed (laboratory design X4 for benefit candidates)", "refuted": ("closed as a computational candidate" if d in BENEFIT else "the line C signal is not replicated in GDSC2 (control not passed); no change to the hazard classification"),
               "not estimable": "status unchanged: no independent screen covers it"}[v])
    (out / "computational_followups_cf9.json").write_text(json.dumps(res, indent=1, sort_keys=True, ensure_ascii=False, default=float) + "\n", encoding="utf-8", newline="\n")
    (out / "computational_followups_cf9.md").write_text(to_markdown(res), encoding="utf-8", newline="\n")


def f(x: Any, nd: int = 3) -> str:
    return "NA" if x is None else f"{x:.{nd}f}"


def to_markdown(r: dict[str, Any]) -> str:
    c = r["line_matching_coverage"]
    L = ["# CF9: independent drug-sensitivity screen (GDSC2) for the line C survivors", "", f"Label: {r['label']}. Protocol: `{r['protocol']}`.", "",
         f"Screen: {r['screen']}; {r['n_gdsc_rows']} dose-response rows, {r['n_gdsc_drugs']} drugs. Response: LN_IC50 (lower = more sensitive); AUC reported as secondary. Statistic: line C adjusted partial rank correlation (lineage fixed effects, proliferation covariate, 1,000 within-lineage permutations, seed {lca.SEED}). Criterion: original sign and p <= {P_CUT}; not estimable if the drug is absent or fewer than {r['min_lines']} matched lines.", "",
         "## Cell-line matching coverage", "",
         f"- GDSC lines: {c['gdsc_lines']}; matched by Sanger id {c['by_sanger_id']}, by COSMIC id {c['by_cosmic_id']}, by name {c['by_name']}; unmatched {c['unmatched']}; ambiguous dropped {c['ambiguous_dropped']}.",
         f"- Matched to a DepMap model: {c['matched_to_depmap_model']}; retained with expression and lineage: {c['retained_with_expression_and_lineage']}.", "",
         "## Per drug", "", "| Drug | In GDSC2 (name matched) | n lines | Original rho | Adjusted rho (LN_IC50) | p | q (BH over estimable) | AUC rho, p | Verdict |", "|---|---|---|---|---|---|---|---|---|"]
    for d, v in r["drugs"].items():
        a, u = v["adjusted"], v.get("adjusted_auc_secondary")
        nm = ", ".join(v["gdsc_names_matched"]) or "absent"
        L.append(f"| {d.title()} | {nm} | {a['n_lines']} | {f(v['original_rho'])} | {f(a.get('rho_adj'))} | {f(a.get('p_perm'), 4)} | {f(a.get('q_bh'), 4)} | "
                 + (f"{f(u.get('rho_adj'))}, {f(u.get('p_perm'), 4)}" if u and u.get("estimable") else "NA") + f" | **{v['verdict']}** |")
    L += ["", f"CF9 overall: {r['cf9_overall']} ({r['n_estimable']} estimable drug(s)).", "", "## Consequence for each candidate in the plan", ""]
    for d, v in r["drugs"].items():
        L.append(f"- {d.title()}: {v['verdict']}; {v['plan_consequence']}.")
    L += ["", "## Limits", "", "- Different assay (72 h fitted IC50, Sanger) and a different line panel from PRISM; replication here is partial overlap, not an exact repeat.",
          "- Correlational, post hoc, from cancer cell lines; does not model partial BUB1B loss; a candidate that survives has only not been refuted.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    main()
