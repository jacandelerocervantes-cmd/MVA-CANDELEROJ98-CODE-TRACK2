"""Line C at drug level adjusted for lineage and proliferation (rule: `protocol/addendum_line_c_adjusted.md`).

    PYTHONPATH=src python -m track2.line_c_adjusted [--root .] [--out out]

Offline: reads data/depmap and out/line_c_drug_level.json only. Correlational, from cancer cell lines; does not model partial BUB1B loss.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

SEED = 20261018
N_PERM = 1000
MIN_LINES = 100
MIN_LINEAGE_DUMMY = 10
MIN_LINEAGE_WITHIN = 30
Q_CUT = 0.10
PROLIF_GENES = ("MKI67", "PCNA", "TOP2A", "MCM2", "CCNB1")
HEME = ("Myeloid", "Lymphoid")
SMALL = "small lineages"


def rank(a: np.ndarray) -> np.ndarray:
    """Average ranks along the last axis."""
    return stats.rankdata(a, axis=-1)


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


def pool_small(lineage: np.ndarray, min_n: int = MIN_LINEAGE_DUMMY) -> np.ndarray:
    """Lineages with fewer than min_n lines (or empty) become one pooled level."""
    vals, counts = np.unique(lineage, return_counts=True)
    keep = {v for v, c in zip(vals, counts) if c >= min_n and v != ""}
    return np.array([v if v in keep else SMALL for v in lineage], dtype=object)


def design(lineage: np.ndarray, prolif_rank: np.ndarray) -> np.ndarray:
    """Fixed-effect design: lineage dummies (full set, so it carries the intercept) and the ranked proliferation covariate."""
    levels = sorted(set(lineage))
    cols = [(lineage == lv).astype(float) for lv in levels]
    return np.column_stack(cols + [prolif_rank - prolif_rank.mean()])


def residualiser(X: np.ndarray) -> np.ndarray:
    """Projection onto the residual space of X (I - hat matrix)."""
    return np.eye(X.shape[0]) - X @ np.linalg.pinv(X)


def within_perms(groups: np.ndarray, n_perm: int, rng: np.random.Generator) -> np.ndarray:
    """(n_perm, n) index matrix; each row permutes positions only inside groups."""
    n = len(groups)
    codes = np.unique(groups, return_inverse=True)[1].astype(float)
    base = np.argsort(codes, kind="stable")
    keys = codes[None, :][:, base] + rng.random((n_perm, n)) * 0.999  # random order inside each block
    order = np.take_along_axis(np.broadcast_to(base, (n_perm, n)), np.argsort(keys, axis=1), axis=1)
    perm = np.empty((n_perm, n), dtype=int)
    perm[:, base] = order
    return perm


def partial_perm(y: np.ndarray, score: np.ndarray, prolif: np.ndarray, lineage: np.ndarray, n_perm: int = N_PERM, seed: int = SEED) -> dict[str, Any]:
    """Partial rank correlation of y with score given lineage fixed effects and prolif, with a within-lineage permutation null."""
    n = len(y)
    lin = pool_small(lineage)
    X = design(lin, rank(prolif))
    M = residualiser(X)
    yr = M @ rank(y)
    yr = yr / np.linalg.norm(yr)
    sr = M @ rank(score)
    rho = float(sr @ yr / np.linalg.norm(sr))
    rng = np.random.default_rng(seed)
    perm = within_perms(lin, n_perm, rng)
    S = rank(score[perm])  # (n_perm, n)
    R = S @ M  # M symmetric
    rp = (R @ yr) / np.linalg.norm(R, axis=1)
    p = (1 + int(np.sum(np.abs(rp) >= abs(rho) - 1e-12))) / (n_perm + 1)
    return {"n_lines": n, "rho_adj": rho, "p_perm": float(p)}


def within_lineage(y: np.ndarray, score: np.ndarray, prolif: np.ndarray, lineage: np.ndarray, min_n: int = MIN_LINEAGE_WITHIN) -> dict[str, Any]:
    """Per-lineage partial rank correlation (proliferation partialled out), combined by DerSimonian-Laird random effects and Stouffer."""
    zs, vs, ns, rhos = [], [], [], []
    for lv in sorted(set(lineage)):
        m = lineage == lv
        if m.sum() < min_n or lv == "":
            continue
        X = np.column_stack([np.ones(m.sum()), rank(prolif[m])])
        M = residualiser(X)
        a, b = M @ rank(y[m]), M @ rank(score[m])
        if np.linalg.norm(a) < 1e-12 or np.linalg.norm(b) < 1e-12:
            continue
        r = float(np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -0.999, 0.999))
        rhos.append(r)
        zs.append(np.arctanh(r))
        vs.append(1.0 / (m.sum() - 4))  # n-3 minus one partialled covariate
        ns.append(int(m.sum()))
    k = len(zs)
    if k < 2:
        return {"n_lineages": k, "pooled_rho": None, "p_random": None, "p_stouffer": None, "i2": None}
    z, v = np.array(zs), np.array(vs)
    w = 1 / v
    fixed = (w * z).sum() / w.sum()
    Q = float((w * (z - fixed) ** 2).sum())
    tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (v + tau2)
    mu = float((wr * z).sum() / wr.sum())
    se = float(np.sqrt(1 / wr.sum()))
    p_re = float(2 * stats.norm.sf(abs(mu / se)))
    wt = np.sqrt(np.array(ns) - 4)
    zst = float((wt * z / np.sqrt(v)).sum() / np.sqrt((wt ** 2).sum()))  # each z / sd is a standard normal
    return {"n_lineages": k, "pooled_rho": float(np.tanh(mu)), "p_random": p_re, "p_stouffer": float(2 * stats.norm.sf(abs(zst))),
            "i2": float(max(0.0, (Q - (k - 1)) / Q)) if Q > 0 else 0.0}


def analyse(Y: dict[str, np.ndarray], score: np.ndarray, prolif: np.ndarray, lineage: np.ndarray, exclude: set[str] | None = None,
            min_lines: int = MIN_LINES, n_perm: int = N_PERM, seed: int = SEED) -> dict[str, dict[str, Any]]:
    """All drugs; lines in `exclude` lineages are dropped. BH across the estimable drugs."""
    keep = np.array([l not in (exclude or set()) for l in lineage])
    rows: dict[str, dict[str, Any]] = {}
    for d in sorted(Y):
        ok = keep & np.isfinite(Y[d]) & np.isfinite(score) & np.isfinite(prolif)
        if ok.sum() < min_lines:
            rows[d] = {"n_lines": int(ok.sum()), "estimable": False}
            continue
        r = partial_perm(Y[d][ok], score[ok], prolif[ok], lineage[ok], n_perm, seed)
        r["estimable"] = True
        r["within"] = within_lineage(Y[d][ok], score[ok], prolif[ok], lineage[ok])
        rows[d] = r
    est = [d for d in rows if rows[d]["estimable"]]
    if est:
        q = bh(np.array([rows[d]["p_perm"] for d in est]))
        for d, qv in zip(est, q):
            rows[d]["q_bh"] = float(qv)
    return rows


def survives(adj: dict[str, Any] | None, orig_rho: float | None) -> str:
    if adj is None or not adj.get("estimable") or orig_rho is None:
        return "not estimable"
    ok = adj["q_bh"] < Q_CUT and np.sign(adj["rho_adj"]) == np.sign(orig_rho)
    return "survives" if ok else "does not survive"


def lolo(y: np.ndarray, score: np.ndarray, prolif: np.ndarray, lineage: np.ndarray, orig_sign: float, min_n: int = MIN_LINEAGE_WITHIN,
         n_perm: int = N_PERM, seed: int = SEED) -> dict[str, Any]:
    out = {}
    for lv in sorted(set(lineage)):
        if (lineage == lv).sum() < min_n or lv == "":
            continue
        m = lineage != lv
        out[lv] = partial_perm(y[m], score[m], prolif[m], lineage[m], n_perm, seed)
    rhos = [v["rho_adj"] for v in out.values()]
    return {"n_omissions": len(out), "n_sign_kept": sum(1 for r in rhos if np.sign(r) == orig_sign),
            "n_sign_kept_and_p_below_0_05": sum(1 for v in out.values() if np.sign(v["rho_adj"]) == orig_sign and v["p_perm"] < 0.05),
            "rho_min": min(rhos) if rhos else None, "rho_max": max(rhos) if rhos else None,
            "worst_omitted_lineage": (min(out, key=lambda k: orig_sign * out[k]["rho_adj"]) if out else None)}


# ------------------------------------------------------------------ data and driver

def load(root: Path, pool_names: list[str]) -> dict[str, Any]:
    from track2 import prism_l4

    data = root / "data" / "depmap"
    need = ("OmicsExpressionProteinCodingGenesTPMLogp1.csv", "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv",
            "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv", "Model.csv")
    for f in need:
        if not (data / f).exists():
            raise FileNotFoundError(f"{f} is not on disk; nothing is downloaded")
    scores = prism_l4.sac_scores(data / need[0], prolif_genes=PROLIF_GENES)
    keys = prism_l4.compound_keys(data / need[2])
    matches = prism_l4.match_rows({n: n for n in pool_names}, keys)
    lin_map = prism_l4.lineages(data / need[3])
    with open(data / need[1], encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        col = {c: i + 1 for i, c in enumerate(header[1:]) if c in scores}
        lines = list(col)
        wanted = {r for rows in matches.values() for r in rows}
        mat = {row[0]: np.array([float(row[col[c]]) if row[col[c]] else np.nan for c in lines]) for row in reader if row[0] in wanted}
    Y = {}
    for d, rows in matches.items():
        with np.errstate(all="ignore"):
            Y[d] = np.nanmean(np.vstack([mat[r] for r in rows]), axis=0)
    return {"Y": Y, "lines": lines, "score": np.array([scores[c]["mean_all"] for c in lines]), "prolif": np.array([scores[c]["prolif"] for c in lines]),
            "lineage": np.array([lin_map.get(c, "") or "" for c in lines], dtype=object)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.out)
    orig = json.loads((out / "line_c_drug_level.json").read_text(encoding="utf-8"))
    pool = sorted(list(orig["drugs"]) + list(orig["not_in_prism"]))
    D = load(root, pool)
    Y, score, prolif, lineage = D["Y"], D["score"], D["prolif"], D["lineage"]
    full = analyse(Y, score, prolif, lineage)
    noheme = analyse(Y, score, prolif, lineage, exclude=set(HEME))
    o = orig["drugs"]
    hits = sorted(d for d, v in o.items() if v.get("q_bh") is not None and v["q_bh"] < Q_CUT)
    cand = ["LAZERTINIB", "NERATINIB", "PREDNISOLONE", "TRAMETINIB", "COBIMETINIB"]
    est = [d for d in full if full[d]["estimable"] and o.get(d, {}).get("rho") is not None]
    rho_o = np.array([o[d]["rho"] for d in est])
    rho_a = np.array([full[d]["rho_adj"] for d in est])
    hit_est = [d for d in hits if d in full and full[d]["estimable"]]
    table: dict[str, Any] = {}
    for d in sorted(set(hits) | set(cand)):
        orig_rho = o.get(d, {}).get("rho")
        row = {"original_rho": orig_rho, "original_q": o.get(d, {}).get("q_bh"), "original_n_lines": o.get(d, {}).get("n_lines"),
               "in_prism": d in Y, "n_heme_lymphoid_lines_in_drug_data": (int(sum(1 for i in np.where(np.isfinite(Y[d]))[0] if lineage[i] in HEME)) if d in Y else None), "adjusted": full.get(d), "excluding_heme_lymphoid": noheme.get(d), "verdict": survives(full.get(d), orig_rho),
               "verdict_excluding_heme_lymphoid": survives(noheme.get(d), orig_rho)}
        if d in Y and full[d]["estimable"] and orig_rho is not None:
            ok = np.isfinite(Y[d])
            row["leave_one_lineage_out"] = lolo(Y[d][ok], score[ok], prolif[ok], lineage[ok], float(np.sign(orig_rho)))
        table[d] = row
    n_surv = sum(1 for d in hits if table[d]["verdict"] == "survives")
    pred = table.get("PREDNISOLONE", {})
    supported = None
    if pred.get("adjusted") and pred.get("excluding_heme_lymphoid") and pred["adjusted"].get("estimable"):
        e = pred["excluding_heme_lymphoid"]
        if not e.get("estimable") or not pred["n_heme_lymphoid_lines_in_drug_data"]:
            supported = None  # not testable: the drug's screen contains no haematopoietic or lymphoid lines
        else:
            supported = bool((pred["original_q"] or 1) < Q_CUT and (e["q_bh"] >= Q_CUT or np.sign(e["rho_adj"]) != np.sign(pred["original_rho"])
                                                                       or abs(e["rho_adj"]) < 0.5 * abs(pred["adjusted"]["rho_adj"])))
    lineage_counts = {k: int((lineage == k).sum()) for k in sorted(set(lineage))}
    res = {"label": "exploratory, correlational, from cancer cell lines; does not model partial BUB1B loss; hypothesis signal only",
           "protocol": "protocol/addendum_line_c_adjusted.md",
           "definition": {"statistic": "partial rank correlation of PRISM viability with the checkpoint score, lineage fixed effects (OncotreeLineage, lineages < 10 lines pooled) and ranked proliferation covariate",
                          "proliferation_covariate": f"mean z of {', '.join(PROLIF_GENES)} (no measured doubling time on disk)", "null": f"{N_PERM} permutations of the score within lineage, seed {SEED}",
                          "survival_rule": "adjusted BH q < 0.10 and same sign as the original rho", "within_lineage_min_lines": MIN_LINEAGE_WITHIN},
           "n_cell_lines": len(D["lines"]), "n_lineages": len(lineage_counts), "lineage_counts": lineage_counts,
           "n_estimable": len(est), "n_original_hits": len(hits), "n_original_hits_estimable": len(hit_est), "n_original_hits_surviving": n_surv,
           "surviving_hits": [d for d in hits if table[d]["verdict"] == "survives"], "non_surviving_hits": [d for d in hits if table[d]["verdict"] != "survives"],
           "n_survive_all_drugs_adjusted_q_below_0_10": sum(1 for d in est if full[d]["q_bh"] < Q_CUT),
           "n_hits_surviving_excluding_heme_lymphoid": sum(1 for d in hits if table[d]["verdict_excluding_heme_lymphoid"] == "survives"),
           "spearman_original_vs_adjusted_all": float(stats.spearmanr(rho_o, rho_a).statistic), "pearson_original_vs_adjusted_all": float(np.corrcoef(rho_o, rho_a)[0, 1]),
           "spearman_original_vs_adjusted_hits": float(stats.spearmanr([o[d]["rho"] for d in hit_est], [full[d]["rho_adj"] for d in hit_est]).statistic) if len(hit_est) > 2 else None,
           "median_abs_ratio_adjusted_to_original_hits": float(np.median([abs(full[d]["rho_adj"]) / abs(o[d]["rho"]) for d in hit_est])) if hit_est else None,
           "line_c_candidates": cand, "drugs": table,
           "prednisolone_lymphoid_suspicion_supported": supported,
           "heme_lymphoid_lines_excluded": int(sum(1 for l in lineage if l in HEME))}
    (out / "line_c_adjusted.json").write_text(json.dumps(res, indent=1, sort_keys=True, ensure_ascii=False, default=float) + "\n", encoding="utf-8", newline="\n")
    (out / "line_c_adjusted.md").write_text(to_markdown(res), encoding="utf-8", newline="\n")


def fmt(x: Any, nd: int = 3) -> str:
    return "NA" if x is None else f"{x:.{nd}f}"


def to_markdown(r: dict[str, Any]) -> str:
    d = r["definition"]
    L = ["# Line C at drug level, adjusted for lineage and proliferation", "",
         f"Label: {r['label']}. Protocol: `{r['protocol']}`.", "",
         f"Statistic: {d['statistic']}. Proliferation covariate: {d['proliferation_covariate']}. Null: {d['null']}. Survival rule: {d['survival_rule']}.", "",
         f"Cell lines {r['n_cell_lines']}, lineages {r['n_lineages']}; estimable drugs {r['n_estimable']}.", "",
         "## Summary", "",
         f"- Original hits (q < 0.10): {r['n_original_hits']} ({r['n_original_hits_estimable']} estimable); surviving the adjustment: **{r['n_original_hits_surviving']}**; surviving after also removing haematopoietic and lymphoid lines: {r['n_hits_surviving_excluding_heme_lymphoid']}.",
         f"- Drugs with adjusted q < 0.10 among all estimable drugs: {r['n_survive_all_drugs_adjusted_q_below_0_10']}.",
         f"- Original versus adjusted effect size over all {r['n_estimable']} estimable drugs: Spearman {fmt(r['spearman_original_vs_adjusted_all'])}, Pearson {fmt(r['pearson_original_vs_adjusted_all'])}; over the original hits: Spearman {fmt(r['spearman_original_vs_adjusted_hits'])}; median |adjusted| / |original| for the hits {fmt(r['median_abs_ratio_adjusted_to_original_hits'])}.",
         f"- Prednisolone lymphoid suspicion supported: {r['prednisolone_lymphoid_suspicion_supported']} (None = not testable: prednisolone has {r['drugs']['PREDNISOLONE']['n_heme_lymphoid_lines_in_drug_data']} haematopoietic or lymphoid lines in its screen data; {r['heme_lymphoid_lines_excluded']} such lines exist in the expression set).", ""]

    def table(drugs: list[str]) -> list[str]:
        T = ["| Drug | Orig rho | Orig q | n | Adj rho | Adj q | Verdict | Within-lineage pooled rho (p RE; k) | Excl. heme/lymphoid rho, q | LOLO sign kept (p<0.05) | LOLO rho range |", "|---|---|---|---|---|---|---|---|---|---|---|"]
        for x in drugs:
            v = r["drugs"][x]
            a, e, lo = v["adjusted"], v["excluding_heme_lymphoid"], v.get("leave_one_lineage_out")
            if a and a.get("estimable"):
                w = a["within"]
                wtxt = f"{fmt(w['pooled_rho'])} ({fmt(w['p_random'], 4)}; {w['n_lineages']})"
                T.append(f"| {x.title()} | {fmt(v['original_rho'])} | {fmt(v['original_q'], 4)} | {v['original_n_lines']} | {fmt(a['rho_adj'])} | {fmt(a['q_bh'], 4)} | {v['verdict']} | {wtxt} | "
                         + (f"{fmt(e['rho_adj'])}, {fmt(e['q_bh'], 4)} ({v['verdict_excluding_heme_lymphoid']})" if e and e.get("estimable") else "NA")
                         + f" | " + (f"{lo['n_sign_kept']}/{lo['n_omissions']} ({lo['n_sign_kept_and_p_below_0_05']})" if lo else "NA") + " | "
                         + (f"{fmt(lo['rho_min'])} to {fmt(lo['rho_max'])}" if lo else "NA") + " |")
            else:
                T.append(f"| {x.title()} | {fmt(v['original_rho'])} | {fmt(v['original_q'], 4)} | {v['original_n_lines']} | NA | NA | {v['verdict']} | NA | NA | NA | NA |")
        return T

    L += ["## Line C candidates", ""] + table(r["line_c_candidates"]) + ["", "## All original hits (q < 0.10)", ""]
    L += table(sorted(r["surviving_hits"] + r["non_surviving_hits"])) + ["", "## Verdicts", ""]
    for x in r["line_c_candidates"]:
        v = r["drugs"][x]
        note = ""
        if v["verdict"] != v["verdict_excluding_heme_lymphoid"]:
            note = f" Without haematopoietic and lymphoid lines: {v['verdict_excluding_heme_lymphoid']}."
        L.append(f"- {x.title()}: {v['verdict']}.{note}")
    L += ["", "## Notes", "",
          "- LOLO: lineages with at least 30 lines are omitted one at a time; the column gives the omissions that keep the original sign (and, in brackets, also have p < 0.05).",
          "- Within-lineage: DerSimonian-Laird random effects over lineages with at least 30 lines (Fisher z of the proliferation-partialled Spearman correlation); reported, not part of the survival rule.",
          "- Limits: correlational, from cancer cell lines; does not model partial BUB1B loss; the gene-set proliferation proxy leaves possible residual confounding. Hypothesis signal only.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    main()
