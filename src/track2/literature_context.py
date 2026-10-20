"""Attach curated, verified literature context to candidate drugs (`config/literature_context.json`).

    python -m track2.literature_context --candidates out/candidates.json --out out

Matching is by the ChEMBL mechanism-of-action text (lower case substring), so it is a transparent rule and not a model.
Every PMID is re-checked to exist at run time; a PMID that cannot be found is dropped and reported. This is context for a
human reviewer: it does not change any score and it is not evidence that a drug works.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from track2.build_candidates import verify_pmids


def classify(candidate: dict[str, Any], config: dict[str, Any]) -> list[dict[str, Any]]:
    text = " ".join(candidate.get("mechanism_of_action", [])).lower()
    hits = [c for c in config["classes"] if any(k in text for k in c["match_any"])]
    out = [{"class": c["id"], "direction_hint": c["direction_hint"], "claim_level": c["claim_level"], "note": c["note"], "pmids": c["pmids"]} for c in hits]
    reports = config["tumour_context_reports"]
    if reports["applies_to_flag"] in candidate.get("risk_flags", []):
        out.append({"class": "mva_tumour_case_reports", "direction_hint": "reported_not_recommended", "claim_level": "case reports",
                    "note": reports["note"], "pmids": reports["pmids"]})
    return out


def run(candidates: list[dict[str, Any]], config: dict[str, Any], titles: dict[str, str]) -> dict[str, Any]:
    result = {}
    for c in candidates:
        notes = classify(c, config)
        for n in notes:
            n["pmids_verified"] = {p: titles.get(p, "NOT FOUND") for p in n["pmids"]}
        result[c["drug_id"]] = notes
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=Path("out/candidates.json"))
    parser.add_argument("--config", type=Path, default=Path("config/literature_context.json"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--name", default="literature_context.json")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    candidates = json.loads(args.candidates.read_text(encoding="utf-8"))
    pmids = sorted({p for c in config["classes"] for p in c["pmids"]} | set(config["tumour_context_reports"]["pmids"]))
    titles = verify_pmids(pmids)
    missing = [p for p in pmids if p not in titles]
    result = run(candidates, config, titles)
    (args.out / args.name).write_text(json.dumps({"pmids_not_found": missing, "drugs": result}, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(sum(1 for v in result.values() if v), "of", len(result), "drugs with context; PMIDs not found:", missing)


if __name__ == "__main__":
    main()
