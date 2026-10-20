"""Independent functional evidence (L4) from DepMap expression and PRISM drug sensitivity (rule: `protocol/addendum_prism.md`).

    python -m track2.prism_l4 --out out

SAC expression score per cell line = mean z-score of log2(TPM+1) of BUB1B, BUB1, BUB3, MAD1L1, MAD2L1 and TTK. Weak-SAC = lowest quartile, strong-SAC =
highest quartile. For each candidate present in PRISM: the difference in mean viability (log fold change) between the groups (negative = the drug kills weak-SAC
cells more), Welch t-test, Benjamini-Hochberg over the tested rows, and a lineage-adjusted regression as the sensitivity analysis. Not clinical evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

from track2.fda_approved import candidate_names

_CFG = json.loads(Path("config/prism_genes.json").read_text(encoding="utf-8")) if Path("config/prism_genes.json").exists() else {}
SAC_GENES = tuple(_CFG.get("score_genes", ("BUB1B", "BUB1", "BUB3", "MAD1L1", "MAD2L1", "TTK")))
PROLIFERATION_GENES = tuple(_CFG.get("proliferation_genes", ("MKI67", "TOP2A", "PCNA", "MCM2", "MCM6", "RRM2")))  # none of them is in the score
MIN_PER_GROUP = 10
FDR = 0.10
DATA = Path("data/depmap")


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


def sac_scores(expression_csv: Path, genes: tuple[str, ...] = SAC_GENES, prolif_genes: tuple[str, ...] = PROLIFERATION_GENES) -> dict[str, dict[str, float]]:
    """Model id -> {'mean_all': mean z of the SAC genes, 'BUB1B': z of BUB1B alone, 'prolif': mean z of the proliferation genes}."""
    allg = tuple(dict.fromkeys(genes + prolif_genes))
    with open(expression_csv, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        index = {g: next(i for i, c in enumerate(header) if c.startswith(g + " (")) for g in allg}
        ids: list[str] = []
        cols: dict[str, list[float]] = {g: [] for g in allg}
        for row in reader:
            ids.append(row[0])
            for g, i in index.items():
                cols[g].append(float(row[i]) if row[i] else np.nan)
    z = {g: (np.array(v) - np.nanmean(v)) / np.nanstd(v) for g, v in cols.items()}
    mean_all = np.nanmean(np.vstack([z[g] for g in genes]), axis=0)
    prolif = np.nanmean(np.vstack([z[g] for g in prolif_genes]), axis=0)
    return {m: {"mean_all": float(mean_all[i]), "BUB1B": float(z["BUB1B"][i]), "prolif": float(prolif[i])} for i, m in enumerate(ids)}


def lineages(model_csv: Path) -> dict[str, str]:
    with open(model_csv, encoding="utf-8", newline="") as handle:
        return {r["ModelID"]: r["OncotreeLineage"] for r in csv.DictReader(handle)}


def compound_keys(compound_csv: Path) -> dict[str, tuple[str, str]]:
    """PRISM row id -> (drug name, all name keys joined with '|'). Keys are upper case."""
    out: dict[str, tuple[str, str]] = {}
    with open(compound_csv, encoding="utf-8", newline="") as handle:
        for r in csv.DictReader(handle):
            names = {r["Drug.Name"].upper()} | {s.strip().upper() for s in (r["Synonyms"] or "").split("|") if s.strip()}
            out[r["IDs"]] = (r["Drug.Name"], "|".join(sorted(names)))
    return out


def match_rows(candidates: dict[str, str], keys: dict[str, tuple[str, str]]) -> dict[str, list[str]]:
    """Candidate drug id -> PRISM row ids whose name or synonym equals the first matching name variant of the candidate."""
    by_key: dict[str, list[str]] = {}
    for row_id, (_, joined) in keys.items():
        for k in joined.split("|"):
            by_key.setdefault(k, []).append(row_id)
    out: dict[str, list[str]] = {}
    for drug_id, name in candidates.items():
        rows: list[str] = []
        for variant in candidate_names(name):
            rows = by_key.get(variant, [])
            if rows:
                break
        if rows:
            out[drug_id] = sorted(set(rows))
    return out


def _ols_coef(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(len(y) - X.shape[1], 1)
    cov = np.linalg.pinv(X.T @ X) * (resid @ resid / dof)
    return beta, np.sqrt(np.diag(cov)), dof


def test_row(viability: np.ndarray, score: np.ndarray, weak: np.ndarray, strong: np.ndarray, lineage: np.ndarray, prolif: np.ndarray | None = None) -> dict[str, Any]:
    a, b = viability[weak], viability[strong]
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < MIN_PER_GROUP or len(b) < MIN_PER_GROUP:
        return {"tested": False}
    t = stats.ttest_ind(a, b, equal_var=False)
    ok = ~np.isnan(viability) & ~np.isnan(score) & (~np.isnan(prolif) if prolif is not None else True)
    levels = sorted(set(lineage[ok]))
    dummies = [(lineage[ok] == lv).astype(float) for lv in levels[1:]]
    y = viability[ok]
    ones = np.ones(int(ok.sum()))
    beta, se, dof = _ols_coef(np.column_stack([score[ok]] + dummies + [ones]), y)  # lineage only
    p_lin = 2 * stats.t.sf(abs(beta[0] / se[0]), dof)
    out = {"tested": True, "n_weak": int(len(a)), "n_strong": int(len(b)), "diff_weak_minus_strong": float(a.mean() - b.mean()), "p": float(t.pvalue),
           "adjusted_coef_per_sd": float(beta[0]), "p_adjusted": float(p_lin), "spearman_all_lines": float(stats.spearmanr(score[ok], y).statistic)}
    if prolif is not None:
        beta2, se2, dof2 = _ols_coef(np.column_stack([score[ok], prolif[ok]] + dummies + [ones]), y)  # lineage and proliferation
        out.update({"coef_lineage_prolif_per_sd": float(beta2[0]), "p_lineage_prolif": float(2 * stats.t.sf(abs(beta2[0] / se2[0]), dof2)), "prolif_coef_per_sd": float(beta2[1])})
    return out


def call_drug(rows: list[dict[str, Any]]) -> str:
    """Drug-level L4 call. A row counts only if it is significant both unadjusted and in the model adjusted for lineage and proliferation, with the same sign
    (protocol/addendum_prism.md, amendment 1); significant only unadjusted is `confounded_by_proliferation_or_lineage`."""
    def adj(r: dict[str, Any]) -> float:
        return r["fdr_lineage_prolif"] if "fdr_lineage_prolif" in r else r["fdr_adjusted"]

    def coef(r: dict[str, Any]) -> float:
        return r["coef_lineage_prolif_per_sd"] if "coef_lineage_prolif_per_sd" in r else r["adjusted_coef_per_sd"]
    confirmed_sens = any(r["fdr"] < FDR and adj(r) < FDR and r["diff_weak_minus_strong"] < 0 and coef(r) > 0 for r in rows)
    confirmed_res = any(r["fdr"] < FDR and adj(r) < FDR and r["diff_weak_minus_strong"] > 0 and coef(r) < 0 for r in rows)
    if confirmed_sens and confirmed_res:
        return "mixed"
    if confirmed_sens:
        return "supports_sensitivity_when_SAC_weak"
    if confirmed_res:
        return "supports_resistance_when_SAC_weak"
    if any(r["fdr"] < FDR for r in rows):
        return "confounded_by_proliferation_or_lineage"
    return "no_association"


def run(candidates: dict[str, str], data: Path = DATA, score_key: str = "mean_all") -> dict[str, Any]:
    scores = sac_scores(data / "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    lin = lineages(data / "Model.csv")
    keys = compound_keys(data / "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv")
    matches = match_rows(candidates, keys)
    with open(data / "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        col = {c: i + 1 for i, c in enumerate(header[1:]) if c in scores}
        lines = list(col)
        wanted = {r for rows in matches.values() for r in rows}
        matrix = {row[0]: np.array([float(row[col[c]]) if row[col[c]] else np.nan for c in lines]) for row in reader if row[0] in wanted}
    score = np.array([scores[c][score_key] for c in lines])
    lo, hi = np.nanquantile(score, 0.25), np.nanquantile(score, 0.75)
    weak, strong = score <= lo, score >= hi
    lineage = np.array([lin.get(c, "Unknown") for c in lines])
    prolif = np.array([scores[c]["prolif"] for c in lines])
    tested: list[dict[str, Any]] = []
    for drug_id, rows in matches.items():
        for r in rows:
            tested.append({"drug_id": drug_id, "prism_row": r, "prism_name": keys[r][0], **test_row(matrix[r], score, weak, strong, lineage, prolif)})
    ok = [t for t in tested if t["tested"]]
    if ok:
        for key, out_key in (("p", "fdr"), ("p_adjusted", "fdr_adjusted"), ("p_lineage_prolif", "fdr_lineage_prolif")):
            for t, v in zip(ok, bh(np.array([t[key] for t in ok]))):
                t[out_key] = float(v)
    per_drug = {}
    for drug_id, name in candidates.items():
        rows = [t for t in ok if t["drug_id"] == drug_id]
        call = call_drug(rows) if rows else ("not_in_PRISM" if drug_id not in matches else "not_testable")
        per_drug[drug_id] = {"name": name, "in_prism": drug_id in matches, "n_rows_tested": len(rows), "call": call, "rows": rows}
    return {"score": score_key, "proliferation_genes": list(PROLIFERATION_GENES), "spearman_score_vs_proliferation": float(stats.spearmanr(score, prolif, nan_policy="omit").statistic), "n_cell_lines": len(lines), "n_weak": int(weak.sum()), "n_strong": int(strong.sum()), "weak_threshold": float(lo),
            "strong_threshold": float(hi), "n_candidates": len(candidates), "n_in_prism": len(matches), "fdr": FDR, "drugs": per_drug}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--data", type=Path, default=DATA)
    args = parser.parse_args()
    candidates: dict[str, str] = {}
    for name in ("proposals.json", "plausibility_e0_e1_e2.json", "plausibility_e0_e1_e2_wide.json"):
        path = args.out / name
        if path.exists():
            data = json.loads(path.read_text(encoding="utf-8"))
            for r in data.get("proposals") or data.get("drugs") or []:
                candidates[r["drug_id"]] = r["name"]
    result = {key: run(candidates, args.data, key) for key in ("mean_all", "BUB1B")}
    (args.out / "prism_l4.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    main_res = result["mean_all"]
    print(main_res["n_cell_lines"], "cell lines;", main_res["n_in_prism"], "of", main_res["n_candidates"], "candidates in PRISM;", dict(Counter(d["call"] for d in main_res["drugs"].values())))


if __name__ == "__main__":
    main()
