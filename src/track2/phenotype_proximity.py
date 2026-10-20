"""Proximity of approved drugs to the phenotype module (line B, `protocol/addendum_four_lines.md`).

    python -m track2.phenotype_proximity --module out/phenotype_module.json --out out

Interpretation of "nine configurations" for a single module (fixed before the run, see `protocol/deviations.md`): the weighting schemes of the
mitosis modules do not apply to one module, so the nine configurations are 3 STRING thresholds (400, 700, 900) x 3 definitions of the module
(genes annotated to at least 2 terms, to at least 3 terms, and the 100 genes with the highest information-content score). A drug is robust if it
is among the 20 best (largest -z, only z < 0) in at least 7 of the 9. Same null model, iterations and seed as the main analysis; drug-target links
are the curated plus measured ones (E1), as in the register. Writes `out/phenotype_threshold_scores.json` and `out/phenotype_proximity.json`.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from track2.augmented_proximity import augmented_universe
from track2.consensus import ROBUST_MIN, THRESHOLDS, TOP_N

ITERATIONS = 1000


def _run(job: tuple[str, str, list[str], dict[str, Any], int]) -> tuple[str, str, dict[str, float]]:
    from track2.run_pipeline import run_mechanism  # imported in the worker: it loads the network there

    variant, score, genes, universe, iterations = job
    mech = {"id": f"phenotype_{variant}", "seeds": {"genes": genes, "reactome_pathways": []}}
    rows = run_mechanism(mech, universe, int(score), iterations)["rows"]
    return variant, score, {r["drug_id"]: round(r["proximity"]["z"], 6) for r in rows}


def robustness(table: dict[str, dict[str, dict[str, float]]], names: dict[str, str]) -> dict[str, Any]:
    """table[variant][score] = {drug: z}. Counts in how many configurations each drug is in the top 20."""
    counts: dict[str, int] = {}
    n_config = 0
    for by_score in table.values():
        for z in by_score.values():
            n_config += 1
            ranked = sorted((d for d, v in z.items() if v < 0), key=lambda d: (z[d], d))[:TOP_N]
            for d in ranked:
                counts[d] = counts.get(d, 0) + 1
    drugs = [{"drug_id": d, "name": names.get(d), "robustness_count": n, "robust": n >= ROBUST_MIN}
             for d, n in sorted(counts.items(), key=lambda kv: (-kv[1], names.get(kv[0]) or ""))]
    return {"rule": {"configurations": n_config, "top_n": TOP_N, "robust_min": ROBUST_MIN}, "n_robust": sum(1 for x in drugs if x["robust"]), "drugs": drugs}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", type=Path, default=Path("out/phenotype_module.json"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--e1", type=Path, default=Path("data/dti_e1_all.json"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--iterations", type=int, default=ITERATIONS)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--strict-links", action="store_true", help="addendum 11, E8: measured links need >= 2 activity records (the primary rule)")
    parser.add_argument("--tag", default="", help="suffix of the output files, e.g. _e5")
    args = parser.parse_args()
    module = json.loads(args.module.read_text(encoding="utf-8"))
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    e1 = json.loads(args.e1.read_text(encoding="utf-8"))["edges"]
    if args.strict_links:
        from track2.core_n2 import restrict_links

        e1 = restrict_links(e1)
    aug = augmented_universe(universe, e1, None)
    jobs = [(v, s, genes, aug, args.iterations) for v, genes in module["variants"].items() for s in THRESHOLDS]
    table: dict[str, dict[str, dict[str, float]]] = {v: {} for v in module["variants"]}
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for variant, score, z in ex.map(_run, jobs):
            table[variant][score] = z
            print(variant, score, "scored:", len(z), "z<=-2:", sum(1 for x in z.values() if x <= -2), flush=True)
    names = {d["drug_id"]: d["name"] for d in universe["drugs"]}
    (args.out / f"phenotype_threshold_scores{args.tag}.json").write_text(json.dumps(table, sort_keys=True), encoding="utf-8", newline="\n")
    result = robustness(table, names)
    result["module"] = {"variants": {v: len(g) for v, g in module["variants"].items()}, "terms": module["terms"]}
    result["links"] = "measured links with >= 2 records, no predicted links" if args.strict_links else "curated plus measured (E1), single-record links included"
    (args.out / f"phenotype_proximity{args.tag}.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(result["n_robust"], "drugs robust;", [d["name"] for d in result["drugs"] if d["robust"]][:30])


if __name__ == "__main__":
    main()
