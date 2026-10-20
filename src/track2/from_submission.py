"""Run the Track 2 analysis starting from the Track 1 submission (rule: `protocol/addendum_submission.md`).

    python -m track2.from_submission --submission <Track 1 submission CSV> --out out

Writes `config/mechanisms_from_submission/sub_<GENE>.json`, `out/scores_sub_<GENE>.json`,
`out/manifest_sub_<GENE>.json` and `out/submission_genes.json` (genes, weights, ClinVar verification, module sizes).
"""

from __future__ import annotations

import argparse
import json
import sys
import re
from pathlib import Path
from typing import Any

from track2.auto_modules import build_mechanism, modellable, module_size
from track2.http_cache import sha256_file
from track2.modules import reactome_version
from track2.run_pipeline import ITERATIONS, STRING_SCORES, UNIVERSE_PATH, load_universe, run_mechanism
from track2.string_net import load_interactome
from track2.submission_input import read_submission, redact_cdna, summarise_genes


COORDINATES = re.compile(r"chr[0-9XYMT]+:\d+(?: [ACGT]+>[ACGT]+)?")


def _dump(path: Path, data: Any) -> None:
    """Written files never carry genomic coordinates or cDNA notation of the proband; they stay in the Track 1 submission CSV."""
    text = redact_cdna(COORDINATES.sub("[variant; coordinates not published]", json.dumps(data, indent=1, sort_keys=True)))
    path.write_text(text + "\n", encoding="utf-8", newline="\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submission", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--config-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    parser.add_argument("--iterations", type=int, default=ITERATIONS)
    args = parser.parse_args()
    if not args.submission.is_file():
        sys.exit(f"The Track 1 submission CSV was not found: {args.submission}\nThe repository does not include it (it carries the variant coordinates). Take it from the Track 1 submission and pass its path with --submission, "
                 "or set the environment variable TRACK1_SUBMISSION before running track2.run_all.")
    args.out.mkdir(parents=True, exist_ok=True)
    args.config_dir.mkdir(parents=True, exist_ok=True)
    genes = summarise_genes(read_submission(args.submission))
    universe = load_universe()
    interactome = load_interactome(STRING_SCORES[0])
    summary = []
    for g in genes:
        evidence = {"rows": g["rows"], "finding_types": g["finding_types"], "clinvar": g["clinvar"]}
        mech = build_mechanism(g["gene"], g["weight"], interactome, evidence)
        size = module_size(mech, interactome)
        entry = {"gene": g["gene"], "weight_epcr": g["weight"], "rows": g["rows"], "finding_types": g["finding_types"],
                 "clinvar": g["clinvar"], "module_rule": mech["module_rule"], "n_pathways": len(mech["pathways"]),
                 "module_size_in_network": size, "modellable": modellable(mech, interactome)}
        summary.append(entry)
        print(g["gene"], g["weight"], mech["module_rule"], size, flush=True)
        if not entry["modellable"]:
            continue
        _dump(args.config_dir / f"{mech['id']}.json", mech)
        runs = {s: run_mechanism(mech, universe, s, args.iterations) for s in STRING_SCORES}
        primary = runs[STRING_SCORES[0]]
        _dump(args.out / f"scores_{mech['id']}.json", primary)
        _dump(args.out / f"threshold_scores_{mech['id']}.json",
              {str(s): {r["drug_id"]: round(r["proximity"]["z"], 6) for r in run["rows"]} for s, run in runs.items()})
        _dump(args.out / f"manifest_{mech['id']}.json", {
            "mechanism": mech["id"], "chembl_version": universe["chembl_version"], "n_approved_drugs": universe["n_drugs"],
            "reactome_version": reactome_version(), "string_version": primary["string_version"],
            "universe_sha256": sha256_file(UNIVERSE_PATH), "iterations": args.iterations,
            "sensitivity": {str(s): {"n_scored": r["n_scored"], "n_pass": r["n_pass"], "top10": [x["name"] for x in r["rows"][:10]]}
                            for s, r in runs.items()}})
    _dump(args.out / "submission_genes.json", {"submission_sha256": sha256_file(args.submission), "genes": summary,
                                               "n_modellable": sum(1 for e in summary if e["modellable"])})
    n_ok = sum(1 for e in summary if e["modellable"])
    print("modellable genes:", n_ok, "of", len(summary))
    if n_ok == 0:
        raise SystemExit("No gene of the submission gives a module of at least 10 proteins in the network: nothing to score (see out/submission_genes.json).")


if __name__ == "__main__":
    main()
