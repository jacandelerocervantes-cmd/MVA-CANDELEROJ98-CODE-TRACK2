"""Merge the per-mechanism scores into `CandidateDrug` records validated against the schema.

    python -m track2.build_candidates --out out

Requires `track2.run_pipeline` to have been run for each mechanism with status `mechanism_candidate`.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import jsonschema

from track2.gates import evaluate_gates
from track2.http_cache import get_json
from track2.labels import label_info
from track2.run_pipeline import load_universe
from track2.toxicity import toxicity_profile

SHORTLIST_TOP_N = 15
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
RISK_TEXT = {
    "ANTINEOPLASTIC_AGENT": "Antineoplastic agent (ATC L01): its use requires a tumour context and may increase toxicity in healthy tissue.",
    "INHIBITS_MODULE_COMPONENT": "Inhibits a component of the module: in a weakened-checkpoint context it could worsen chromosomal instability; the direction must be reviewed.",
    "BOXED_WARNING": "Has a boxed warning in ChEMBL.",
    "MICROTUBULE_OR_PROTEIN_STABILISER": "Stabiliser (for example of microtubules): depends on a functional checkpoint for its antimitotic effect.",
}


def verify_pmids(pmids: list[str]) -> dict[str, str]:
    """PMID -> title, only for those that exist in PubMed (batched E-utilities query)."""
    found: dict[str, str] = {}
    unique = sorted(set(pmids))
    for i in range(0, len(unique), 150):
        chunk = unique[i:i + 150]
        data = get_json(ESUMMARY + ",".join(chunk))["result"]
        for pmid in chunk:
            entry = data.get(pmid)
            if entry and "error" not in entry and entry.get("title"):
                found[pmid] = entry["title"]
    return found


def build(out_dir: Path, mechanism_dir: Path = Path("config/mechanisms_from_submission")) -> list[dict[str, Any]]:
    universe = load_universe()
    drugs = {d["drug_id"]: d for d in universe["drugs"]}
    patient = json.loads(Path("config/patient.json").read_text(encoding="utf-8"))
    schema = json.loads(Path("schemas/candidate_drug.schema.json").read_text(encoding="utf-8"))
    mechanisms = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(mechanism_dir.glob("*.json"))]
    mechanisms = [m for m in mechanisms if m.get("status") == "mechanism_candidate"]
    scores: dict[str, dict[str, dict[str, Any]]] = {}
    weights = {}
    for m in mechanisms:
        weights[m["id"]] = float(m["weight"])
        data = json.loads((out_dir / f"scores_{m['id']}.json").read_text(encoding="utf-8"))
        scores[m["id"]] = {r["drug_id"]: r for r in data["rows"]}
    shortlist: set[str] = set()
    for mid, rows in scores.items():
        ordered = sorted(rows.values(), key=lambda r: r["rank"])
        shortlist |= {r["drug_id"] for r in ordered if r["passes_primary_threshold"] or r["rank"] <= SHORTLIST_TOP_N}
    pmids = [ref["id"] for d in shortlist for mech in drugs[d]["mechanisms"] for ref in mech["refs"] if ref["type"] == "PubMed"]
    titles = verify_pmids(pmids)
    records = []
    for drug_id in sorted(shortlist):
        drug = drugs[drug_id]
        label = label_info(drug["name"] or "")
        mech_scores = []
        combined = 0.0
        flags: set[str] = set()
        actions: set[str] = set()
        for mid, rows in scores.items():
            row = rows.get(drug_id)
            if not row:
                continue
            prox = row["proximity"]
            mech_scores.append({"mechanism_id": mid, "z": prox["z"], "d": prox["d"], "p_empirical": prox["p_empirical"], "rank": row["rank"],
                                "passes_primary_threshold": row["passes_primary_threshold"]})
            combined += weights[mid] * max(0.0, -prox["z"])
            flags |= set(row["risk_flags"])
            actions |= set(row["actions"])
        evidence_for = []
        for mech in drug["mechanisms"]:
            for ref in mech["refs"]:
                if ref["type"] == "PubMed":
                    evidence_for.append({"type": "PubMed (ChEMBL mechanism reference)", "id": ref["id"], "title": titles.get(ref["id"], ""),
                                         "citation_verified": ref["id"] in titles})
                elif ref["type"] == "DailyMed":
                    evidence_for.append({"type": "DailyMed (ChEMBL mechanism reference)", "id": ref["id"], "citation_verified": False})
        age = dict(label["age_profile"])
        age.update({"pediatric_label_status": label["pediatric"]["status"], "pediatric_excerpt": label["pediatric"]["excerpt"], "heuristic": True})
        record = {
            "drug_id": drug_id, "name": drug["name"],
            "approval": {"source": "ChEMBL", "chembl_version": universe["chembl_version"], "max_phase": 4, "first_approval": drug["first_approval"],
                         "atc": drug["atc"], "label": {k: label.get(k) for k in ("found", "setid", "title", "has_boxed_warning")}},
            "targets": drug["human_targets"], "mechanism_of_action": sorted({m["mechanism"] for m in drug["mechanisms"] if m["mechanism"]}),
            "mechanism_scores": sorted(mech_scores, key=lambda s: s["mechanism_id"]),
            "combined_score": round(combined, 4), "actions": sorted(actions), "risk_flags": sorted(flags), "direction_review": "REQUIRES_REVIEW",
            "age_profile": age, "gates": evaluate_gates(patient, label, drug), "evidence_for": evidence_for,
            "evidence_against": [{"source": "rule", "text": RISK_TEXT[f]} for f in sorted(flags) if f in RISK_TEXT]
            + [{"source": "DailyMed keyword scan (heuristic)", "text": f"The label mentions {k.lower()}."} for k in label.get("keyword_flags", [])]
            + [{"source": "excipient review (author-curated list)", "text": f"{a['excipient']}: {a['concern']}."} for a in label.get("excipient_alerts", [])]
            + ([{"source": "DailyMed contraindications (heuristic)", "text": "The label lists hypersensitivity among its contraindications."}] if label.get("contraindicated_hypersensitivity") else []),
            "toxicity": toxicity_profile(drug["name"] or ""),
            "excipients": {"inactive_ingredients": label.get("inactive_ingredients", []), "alerts": label.get("excipient_alerts", [])},
        }
        jsonschema.validate(record, schema)
        records.append(record)
    records.sort(key=lambda r: (-r["combined_score"], r["drug_id"]))
    return records


def write_outputs(records: list[dict[str, Any]], out_dir: Path) -> None:
    (out_dir / "candidates.json").write_text(json.dumps(records, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(out_dir / "candidates.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["rank", "drug_id", "name", "combined_score", "best_mechanism", "best_z", "risk_flags", "adult_label", "pediatric_label_status", "nephrocalcinosis_caution", "hypersensitivity_contraindication", "excipient_alerts", "faers_reports"])
        for i, r in enumerate(records, 1):
            best = min(r["mechanism_scores"], key=lambda s: s["z"])
            writer.writerow([i, r["drug_id"], r["name"], r["combined_score"], best["mechanism_id"], round(best["z"], 3),
                             ";".join(r["risk_flags"]), r["age_profile"]["adult"], r["age_profile"]["pediatric_label_status"],
                             next(g["status"] for g in r["gates"] if g["requirement"] == "NEPHROCALCINOSIS_CAUTION"),
                             any("hypersensitivity among its contraindications" in e["text"] for e in r["evidence_against"]),
                             ";".join(a["excipient"] for a in r["excipients"]["alerts"]),
                             r["toxicity"].get("total_reports", "")])
    with open(out_dir / "candidates_by_age_band.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["drug_id", "name", "neonate", "infant", "child", "adolescent", "adult", "heuristic_note"])
        for r in records:
            a = r["age_profile"]
            writer.writerow([r["drug_id"], r["name"], a["neonate"], a["infant"], a["child"], a["adolescent"], a["adult"], "label extracted heuristically; verify by hand"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--dest", type=Path, default=None, help="where to write candidates.* (default: --out)")
    parser.add_argument("--mechanism-dir", type=Path, default=Path("config/mechanisms_from_submission"))
    args = parser.parse_args()
    records = build(args.out, args.mechanism_dir)
    dest = args.dest or args.out
    dest.mkdir(parents=True, exist_ok=True)
    write_outputs(records, dest)
    print(len(records), "candidates written")


if __name__ == "__main__":
    main()
