"""Run the Track 2 pipeline for one mechanism: approved universe -> module -> proximity -> score.

    python -m track2.run_pipeline --mechanism config/mechanisms/bub1b_sac_lof.json --out out
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


from track2.chembl_universe import build_universe
from track2.http_cache import sha256_file
from track2.modules import build_module, load_mechanism, reactome_uniprot, reactome_version
from track2.proximity import build_null, min_distance_to_set
from track2.scoring import PROXIMITY_Z_PRIMARY, score_drugs
from track2.string_net import load_interactome

UNIVERSE_PATH = Path("data/universe.json")
ITERATIONS = 1000
STRING_SCORES = (700, 400, 900)  # the first is the main one; the others are sensitivity


def load_universe() -> dict[str, Any]:
    if not UNIVERSE_PATH.exists():
        UNIVERSE_PATH.parent.mkdir(parents=True, exist_ok=True)
        UNIVERSE_PATH.write_text(json.dumps(build_universe(), sort_keys=True), encoding="utf-8", newline="\n")
    return json.loads(UNIVERSE_PATH.read_text(encoding="utf-8"))


def run_mechanism(mechanism: dict[str, Any], universe: dict[str, Any], min_score: int, iterations: int = ITERATIONS) -> dict[str, Any]:
    interactome = load_interactome(min_score)
    module = build_module(mechanism, interactome)
    module_uniprot = set()
    for pathway in mechanism["seeds"].get("reactome_pathways", []):
        module_uniprot.update(reactome_uniprot(pathway))
    null = build_null(interactome.adjacency, interactome.degree, module["nodes"], repeats=iterations)
    rows = score_drugs(universe["drugs"], interactome, module["nodes"], module_uniprot, null)
    return {
        "mechanism_id": mechanism["id"], "string_min_score": min_score, "string_version": interactome.string_version,
        "n_network_nodes": len(interactome.proteins), "module": {k: v for k, v in module.items() if k != "nodes"},
        "iterations": iterations, "z_primary_threshold": PROXIMITY_Z_PRIMARY, "n_scored": len(rows),
        "n_pass": sum(1 for r in rows if r["passes_primary_threshold"]), "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--iterations", type=int, default=ITERATIONS)
    args = parser.parse_args()
    mechanism = load_mechanism(args.mechanism)
    universe = load_universe()
    args.out.mkdir(parents=True, exist_ok=True)
    runs = {}
    for min_score in STRING_SCORES:
        runs[min_score] = run_mechanism(mechanism, universe, min_score, args.iterations)
    primary = runs[STRING_SCORES[0]]
    (args.out / f"scores_{mechanism['id']}.json").write_text(json.dumps(primary, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    sensitivity = {
        str(score): {"n_scored": r["n_scored"], "n_pass": r["n_pass"], "top10": [x["name"] for x in r["rows"][:10]]}
        for score, r in runs.items()
    }
    manifest = {
        "mechanism": mechanism["id"], "chembl_version": universe["chembl_version"], "n_approved_drugs": universe["n_drugs"],
        "reactome_version": reactome_version(), "string_version": primary["string_version"],
        "universe_sha256": sha256_file(UNIVERSE_PATH), "iterations": args.iterations, "sensitivity": sensitivity,
    }
    (args.out / f"manifest_{mechanism['id']}.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(manifest["sensitivity"], indent=1))


if __name__ == "__main__":
    main()
