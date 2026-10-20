"""Provisional direction of effect on the spindle checkpoint module, from Gene Ontology annotations.

For each (drug, target) pair the ChEMBL action type is combined with the QuickGO annotation of the
target (UniProt) to the checkpoint GO terms listed in `config/direction_go_terms.json`.
[heuristic] Provisional, annotation-based, spindle checkpoint module only. It is a hypothesis
generator for human review, never a statement about dose, treatment or recommendation.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

from track2.http_cache import get_json

QUICKGO = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
BATCH = 20
PAGE = 200
CAVEAT = (
    "provisional; GO annotation-based; applies only to the spindle checkpoint module; "
    "other modules are not assessed"
)
INHIBITORY = {
    "INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR",
    "INVERSE AGONIST", "RNAI INHIBITOR", "ANTISENSE INHIBITOR", "DEGRADER", "ALLOSTERIC ANTAGONIST",
    "GENE EDITING NEGATIVE MODULATOR",
}
ACTIVATING = {
    "AGONIST", "PARTIAL AGONIST", "ACTIVATOR", "OPENER", "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR", "STABILISER",
}
WEAKEN, STRENGTHEN = "may_weaken_checkpoint", "may_strengthen_checkpoint"
NO_ANN, UNCLEAR = "no_checkpoint_annotation", "unclear"
CLASSES = (WEAKEN, STRENGTHEN, NO_ANN, UNCLEAR)


def load_terms(path: Path) -> list[dict[str, Any]]:
    return list(json.loads(path.read_text(encoding="utf-8"))["terms"])


def fetch_annotations(
    accessions: list[str], terms: list[dict[str, Any]], problems: list[str]
) -> dict[str, list[dict[str, Any]]]:
    """UniProt -> annotations (GO descendants included) with root term, role and evidence."""
    out: dict[str, dict[tuple[str, str, str, str], dict[str, Any]]] = {a: {} for a in accessions}
    for term in terms:
        for i in range(0, len(accessions), BATCH):
            batch = accessions[i : i + BATCH]
            page = 1
            while True:
                url = (
                    f"{QUICKGO}?geneProductId={','.join(batch)}&goId={term['id']}"
                    f"&goUsage=descendants&taxonId=9606&limit={PAGE}&page={page}"
                )
                try:
                    data = get_json(url)
                except (LookupError, RuntimeError) as error:
                    problems.append(f"{term['id']} batch {i // BATCH}: {error}")
                    break
                for r in data.get("results", []):
                    if str(r.get("qualifier", "")).startswith("NOT"):
                        continue  # negated annotation: the protein does NOT have this role
                    acc = str(r["geneProductId"]).split(":", 1)[-1]
                    if acc not in out:
                        continue
                    ev = {
                        "go_id": r["goId"], "root_term": term["id"], "role": term["role"],
                        "qualifier": r.get("qualifier"), "evidence": r.get("goEvidence"),
                        "eco": r.get("evidenceCode"), "reference": r.get("reference"),
                    }
                    key = (ev["go_id"], ev["root_term"], str(ev["evidence"]), str(ev["reference"]))
                    out[acc][key] = ev
                if page >= int(data.get("pageInfo", {}).get("total", 1)):
                    break
                page += 1
    return {a: [v[k] for k in sorted(v)] for a, v in out.items()}


def target_role(annotations: list[dict[str, Any]]) -> str:
    if not annotations:
        return "none"
    roles = {a["role"] for a in annotations if a["role"] != "neutral"}
    if roles == {"positive"}:
        return "positive"
    if roles == {"negative"}:
        return "negative"
    return "conflicting" if roles else "neutral"


def classify_pair(action: str | None, role: str) -> str:
    if role == "none":
        return NO_ANN
    if role not in ("positive", "negative"):
        return UNCLEAR  # annotated but direction of the regulator unspecified, or conflicting
    if action in INHIBITORY:
        return WEAKEN if role == "positive" else STRENGTHEN
    if action in ACTIVATING:
        return STRENGTHEN if role == "positive" else WEAKEN
    return UNCLEAR


def drug_pairs(drug: dict[str, Any], targets: dict[str, Any]) -> list[tuple[str, str | None]]:
    seen: set[tuple[str, str | None]] = set()
    for m in drug.get("mechanisms", []):
        for acc in (targets.get(m.get("target_id")) or {}).get("uniprot", []):
            seen.add((acc, m.get("action")))
    return sorted(seen, key=lambda p: (p[0], str(p[1])))


def build_review(universe: dict[str, Any], ids: list[str], terms: list[dict[str, Any]]) -> dict[str, Any]:
    by_id = {d["drug_id"]: d for d in universe["drugs"]}
    pairs = {i: drug_pairs(by_id[i], universe["targets"]) for i in ids}
    accessions = sorted({acc for p in pairs.values() for acc, _ in p})
    problems: list[str] = []
    ann = fetch_annotations(accessions, terms, problems)
    drugs: dict[str, Any] = {}
    for i in ids:
        items = []
        for acc, action in pairs[i]:
            role = target_role(ann[acc])
            items.append({
                "uniprot": acc, "action": action, "target_role": role,
                "classification": classify_pair(action, role), "evidence": ann[acc],
            })
        classes = sorted({x["classification"] for x in items})
        drugs[i] = {
            "name": by_id[i].get("name"), "caveat": CAVEAT, "pairs": items,
            "classifications": classes or [NO_ANN],
        }
    counts = Counter(c for d in drugs.values() for c in d["classifications"])
    return {
        "caveat": CAVEAT, "api_problems": sorted(problems), "drugs": drugs,
        "summary": {
            "n_drugs": len(drugs), "drugs_per_class": {c: counts.get(c, 0) for c in CLASSES},
            "n_targets_queried": len(accessions), "go_terms": [t["id"] for t in terms],
        },
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--universe", default="data/universe.json")
    p.add_argument("--candidates", default="out/candidates.json")
    p.add_argument("--terms", default="config/direction_go_terms.json")
    p.add_argument("--out", default="out")
    a = p.parse_args()
    universe = json.loads(Path(a.universe).read_text(encoding="utf-8"))
    ids = sorted(c["drug_id"] for c in json.loads(Path(a.candidates).read_text(encoding="utf-8")))
    review = build_review(universe, ids, load_terms(Path(a.terms)))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "direction_review.json").write_text(
        json.dumps(review, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    print(json.dumps(review["summary"], indent=2, sort_keys=True))
    print("api problems:", len(review["api_problems"]))


if __name__ == "__main__":
    main()
