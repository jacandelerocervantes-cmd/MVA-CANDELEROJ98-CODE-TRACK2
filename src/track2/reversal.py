"""Transcriptomic-reversal evidence layer (SigCom LINCS, L1000). Hypothesis-generating only.

For the CRISPR-knockout consensus signature of BUB1B (primary) and of other spindle-checkpoint/kinetochore genes
(proxies), take the top 100 up and 100 down genes, ask SigCom LINCS which chemical-perturbation consensus signatures
move against them (reversers) or with them (mimickers), and aggregate per approved drug. A pooled null from random
gene sets of the same size gives empirical p-values. Usage:
    python -m track2.reversal --universe data/universe.json --candidates out/candidates.json --out out
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any

import numpy as np

from track2 import reversal_names, reversal_stats as st, sigcom_client as sc
from track2.reversal_report import write_summary

PRIMARY = "BUB1B"
PROXIES = ["BUB1", "TTK", "ZWINT", "CENPE"]  # the other consensus CRISPR-KO signatures available for checkpoint genes
TOP_N = 100
N_PERMUTATIONS = 30
SEED = 20260101
CELL_TOP = 2000  # reversers requested from the cell-line-level library
ALPHA = 0.05
N_KO_CONTROLS = 20  # structured control: consensus CRISPR-KO signatures of unrelated, randomly chosen genes
FULL = 40000  # larger than the 33,608 consensus signatures: returns the whole library


def query_set(gene: str) -> dict[str, Any]:
    values = sc.consensus_signature("xpr_mean_sigs", gene)
    ordered = sorted(values.items(), key=lambda kv: (kv[1], kv[0]))
    down = [g for g, _ in ordered[:TOP_N]]
    up = [g for g, _ in ordered[-TOP_N:]][::-1]
    ids = sc.gene_ids(up + down)
    return {"gene": gene, "up": up, "down": down, "n_genes_in_signature": len(values),
            "up_ids": [ids[g] for g in up if g in ids], "down_ids": [ids[g] for g in down if g in ids]}


def random_sets(symbols: list[str]) -> list[tuple[list[str], list[str]]]:
    rng = random.Random(SEED)
    out = []
    for _ in range(N_PERMUTATIONS):
        pick = rng.sample(symbols, 2 * TOP_N)
        ids = sc.gene_ids(pick)
        out.append(([ids[g] for g in pick[:TOP_N] if g in ids], [ids[g] for g in pick[TOP_N:] if g in ids]))
    return out


def _calls(rev: np.ndarray, mim: np.ndarray, null_rev: np.ndarray, null_mim: np.ndarray) -> dict[str, np.ndarray]:
    p_rev, p_mim = st.empirical_p(rev, null_rev), st.empirical_p(mim, null_mim)
    return {"rev": rev, "mim": mim, "p_rev": p_rev, "fdr_rev": st.bh_fdr(p_rev), "p_mim": p_mim, "fdr_mim": st.bh_fdr(p_mim)}


def _vec(obs: dict[str, dict[str, Any]], drugs: list[str], key: str) -> np.ndarray:
    return np.array([obs.get(d, {}).get(key, 0.0) for d in drugs])


def cell_level(primary: dict[str, Any], index: dict[str, list[str]]) -> dict[str, list[tuple[str, float]]]:
    """Per drug: (cell line, reversal score) of individual L1000 signatures among the top reversers of BUB1B loss."""
    rows = [r for r in sc.two_sided(primary["up_ids"], primary["down_ids"], sc.LIB_CHEM_CELL, CELL_TOP) if r["z-sum"] < 0]
    meta = sc.signature_meta([r["uuid"] for r in rows], ["meta.pert_name", "meta.cell_line"])
    out: dict[str, list[tuple[str, float]]] = {}
    for r in rows:
        m = meta.get(r["uuid"], {})
        for d in reversal_names.match(m.get("pert_name", ""), index):
            out.setdefault(d, []).append((m.get("cell_line", "?"), -float(r["z-sum"])))
    return out


def run(universe_path: Path, candidates_path: Path, out_dir: Path) -> dict[str, Any]:
    universe = json.loads(universe_path.read_text(encoding="utf-8"))["drugs"]
    candidates = {c["drug_id"] for c in json.loads(candidates_path.read_text(encoding="utf-8"))}
    names = {d["drug_id"]: d["name"] for d in universe}
    index = reversal_names.drug_index(universe)
    sig_meta = sc.library_signature_names(sc.LIB_CHEM_CONSENSUS, ["meta.pert_name"])
    drug_of = {u: reversal_names.match(m["pert_name"], index) for u, m in sig_meta.items()}
    drug_of = {u: d for u, d in drug_of.items() if d}
    scored = sorted({d for v in drug_of.values() for d in v})
    pos = {d: i for i, d in enumerate(scored)}

    sets = {g: query_set(g) for g in [PRIMARY] + PROXIES}
    observed = {g: st.per_drug(sc.two_sided(s["up_ids"], s["down_ids"], sc.LIB_CHEM_CONSENSUS, FULL), drug_of)
                for g, s in sets.items()}
    # Negative control: random gene sets of the same size through the same pipeline (drug level).
    nulls = [st.per_drug(sc.two_sided(u, d, sc.LIB_CHEM_CONSENSUS, FULL), drug_of)
             for u, d in random_sets(sorted(sc.consensus_signature("xpr_mean_sigs", PRIMARY)))]
    null_rev = [_vec(n, scored, "rev") for n in nulls]
    null_mim = [_vec(n, scored, "mim") for n in nulls]
    pooled_rev, pooled_mim = np.concatenate(null_rev), np.concatenate(null_mim)
    stats = {g: _calls(_vec(o, scored, "rev"), _vec(o, scored, "mim"), pooled_rev, pooled_mim) for g, o in observed.items()}

    false_hits = []
    for i in range(len(nulls)):  # leave-one-out: each random set as if it were the observed query
        rest = [j for j in range(len(nulls)) if j != i]
        c = _calls(null_rev[i], null_mim[i], np.concatenate([null_rev[j] for j in rest]), np.concatenate([null_mim[j] for j in rest]))
        false_hits.append(int(sum(st.direction(a, b, ALPHA) == "reverses" for a, b in zip(c["fdr_rev"], c["fdr_mim"]))))

    pool = sorted({m["pert_name"] for m in sc.library_signature_names(sc.LIB_CRISPR_CONSENSUS, ["meta.pert_name"]).values()}
                  - set(sets))
    controls = []
    ctrl_genes = sorted(random.Random(SEED + 1).sample([str(x) for x in pool], N_KO_CONTROLS))
    for g in ctrl_genes:
        q = query_set(g)
        controls.append(_vec(st.per_drug(sc.two_sided(q["up_ids"], q["down_ids"], sc.LIB_CHEM_CONSENSUS, FULL), drug_of), scored, "rev"))
    ctrl = np.array(controls)

    cells = cell_level(sets[PRIMARY], index)
    scores: dict[str, Any] = {}
    for d in sorted(names):
        if d not in pos:
            scores[d] = {"available": False, "name": names[d],
                         "reason": "no LINCS L1000 chemical consensus signature with a matching perturbagen name"}
            continue
        i, s0 = pos[d], stats[PRIMARY]
        sig_ids = sorted(u for u, ds in drug_of.items() if d in ds)
        cl = cells.get(d, [])
        scores[d] = {
            "available": True, "name": names[d], "in_candidates": d in candidates,
            "matched_perturbagens": sorted({sig_meta[u]["pert_name"] for u in sig_ids}),
            "name_match": "ambiguous" if len(index[reversal_names.clean(names[d])]) > 1 else "unique",
            "score": round(float(s0["rev"][i]), 4), "p": round(float(s0["p_rev"][i]), 6), "fdr": round(float(s0["fdr_rev"][i]), 6),
            "direction": st.direction(float(s0["fdr_rev"][i]), float(s0["fdr_mim"][i]), ALPHA),
            "mimic_score": round(float(s0["mim"][i]), 4), "p_mimic": round(float(s0["p_mim"][i]), 6),
            "fdr_mimic": round(float(s0["fdr_mim"][i]), 6),
            "api_fdr_best_reverser": round(float(observed[PRIMARY].get(d, {}).get("api_fdr_rev", 1.0)), 6),
            "ko_control": {"n_controls": N_KO_CONTROLS, "n_controls_with_score_ge_observed": int((ctrl[:, i] >= s0["rev"][i]).sum()),
                           "p": round(float((1 + (ctrl[:, i] >= s0["rev"][i]).sum()) / (1 + N_KO_CONTROLS)), 4),
                           "median_control_score": round(float(np.median(ctrl[:, i])), 4)},
            "n_signatures": len(sig_ids), "source_signature_ids": sig_ids,
            "cell_lines": sorted({c for c, _ in cl}), "n_cell_line_signatures_in_top_reversers": len(cl),
            "best_cell_line_score": round(max((z for _, z in cl), default=0.0), 4),
            "proxy_sac": {g: {"score": round(float(stats[g]["rev"][i]), 4), "p": round(float(stats[g]["p_rev"][i]), 6),
                              "fdr": round(float(stats[g]["fdr_rev"][i]), 6),
                              "direction": st.direction(float(stats[g]["fdr_rev"][i]), float(stats[g]["fdr_mim"][i]), ALPHA)}
                          for g in PROXIES},
        }
    meta = {
        "query_sets": {g: {"up": s["up"], "down": s["down"], "n_genes_in_signature": s["n_genes_in_signature"],
                           "n_up_mapped": len(s["up_ids"]), "n_down_mapped": len(s["down_ids"])} for g, s in sets.items()},
        "libraries": {"chemical_consensus": "l1000_mean_cp " + sc.LIB_CHEM_CONSENSUS, "chemical_cell_level": "l1000_cp " + sc.LIB_CHEM_CELL,
                      "gene_signature_source": "LINCS L1000 CRISPR KO consensus (l1000_mean_xpr) via S3 means/xpr_mean_sigs"},
        "top_n": TOP_N, "n_permutations": N_PERMUTATIONS, "ko_control_genes": ctrl_genes, "seed": SEED, "cell_level_top_reversers": CELL_TOP,
        "n_universe_drugs": len(names), "n_drugs_with_score": len(scored),
        "n_direction_reverses": sum(v.get("direction") == "reverses" for v in scores.values()),
        "n_reverses_and_specific": sum(v.get("direction") == "reverses" and v["ko_control"]["p"] <= ALPHA for v in scores.values()),
        "n_direction_mimics": sum(v.get("direction") == "mimics" for v in scores.values()),
        "null": st.null_summary(pooled_rev, false_hits, len(scored)),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, obj in (("reversal_scores.json", scores), ("reversal_meta.json", meta)):
        (out_dir / name).write_text(json.dumps(obj, sort_keys=True, indent=1) + "\n", encoding="utf-8", newline="\n")
    write_summary(out_dir / "reversal_summary.md", scores, meta)
    return meta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--candidates", type=Path, default=Path("out/candidates.json"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    meta = run(args.universe, args.candidates, args.out)
    print(json.dumps({k: meta[k] for k in ("n_drugs_with_score", "n_direction_reverses", "n_direction_mimics")}))


if __name__ == "__main__":
    main()
