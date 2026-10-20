"""Merge the optional evidence layers into one enriched candidate table.

    python -m track2.integrate --out out

Reads `out/candidates.json` plus, when present, `evidence_retrieval.json`, `faers_disproportionality.json`,
`pgx_drug_relevance.json`, `direction_review.json`, `reversal_scores.json`, `physchem_adme.json` and `consensus_ranking.json`, and writes
`out/candidates_enriched.json` and `out/candidates_enriched.csv`. No layer changes the ranking: they are
context for a human reviewer. All content is hypothesis support, not medical advice.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

SIGNAL_GROUPS = ("renal", "hepatic", "haematological", "cardiac", "neurological", "hypersensitivity")


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def _evidence(entry: dict[str, Any] | None) -> dict[str, Any]:
    if not entry:
        return {"available": False}
    pubmed = {k: v.get("count") for k, v in entry["pubmed"].items() if v.get("available")}
    trials = {k: v.get("total_count") for k, v in entry["trials"].items() if v.get("available")}
    return {"available": True, "pubmed_counts": pubmed, "trial_counts": trials}


def _faers(entry: dict[str, Any] | None) -> dict[str, Any]:
    if not entry or not entry.get("available"):
        return {"available": False}
    signals = sorted(g for g, v in entry["groups"].items() if v.get("signal"))
    return {"available": True, "signal_groups": signals,
            "prr": {g: v.get("prr") for g, v in entry["groups"].items() if v.get("estimable")}}


def _reversal(entry: dict[str, Any] | None) -> dict[str, Any]:
    if not entry or not entry.get("available"):
        return {"available": False}
    return {"available": True, "direction": entry["direction"], "fdr": entry["fdr"], "ko_control_p": entry["ko_control"]["p"]}


def _physchem(entry: dict[str, Any] | None) -> dict[str, Any]:
    if not entry:
        return {"available": False}
    props = entry["properties"]
    aff = entry["measured_affinity_to_mechanism_targets"]
    return {"available": True, "molecule_type": entry["molecule_type"], "has_small_molecule_properties": entry["has_small_molecule_properties"],
            "mw": props.get("full_mwt"), "alogp": props.get("alogp"), "psa": props.get("psa"), "ro5_violations": props.get("num_ro5_violations"),
            "n_targets_with_measured_affinity": len(aff), "best_median_pactivity": max((v["median_pactivity"] for v in aff.values()), default=None),
            "n_adme_measurement_types": len(entry["measured_adme"]), "routes": [k for k, v in entry["routes"].items() if v]}


def integrate(out: Path) -> list[dict[str, Any]]:
    candidates = json.loads((out / "candidates.json").read_text(encoding="utf-8"))
    evidence, faers = _load(out / "evidence_retrieval.json"), _load(out / "faers_disproportionality.json")
    pgx, direction, reversal = _load(out / "pgx_drug_relevance.json"), _load(out / "direction_review.json"), _load(out / "reversal_scores.json")
    physchem = _load(out / "physchem_adme.json")
    consensus = {d["drug_id"]: d for d in (_load(out / "consensus_ranking.json") or {}).get("drugs", [])}
    rows = []
    for record in candidates:
        drug_id = record["drug_id"]
        row = dict(record)
        row["layers"] = {
            "evidence_retrieved": _evidence((evidence or {}).get(drug_id)),
            "faers_disproportionality": _faers(((faers or {}).get("drugs") or {}).get(drug_id)),
            "pharmacogenomics_cpic": [{"gene": p["gene"], "cpic_level": p.get("cpic_level")} for p in (pgx or {}).get(drug_id, [])],
            "checkpoint_direction": ((direction or {}).get("drugs") or {}).get(drug_id, {}).get("classifications", ["not_run"]),
            "transcriptomic_reversal": _reversal((reversal or {}).get(drug_id)),
            "physchem_affinity_adme": _physchem((physchem or {}).get(drug_id)),
            "robustness": {"robustness_count": consensus[drug_id]["robustness_count"], "of": consensus[drug_id]["of"], "class": consensus[drug_id]["class"]}
            if drug_id in consensus else {"robustness_count": 0, "of": 9, "class": "not_in_top20_of_any_configuration"},
        }
        rows.append(row)
    return rows


def write(rows: list[dict[str, Any]], out: Path) -> None:
    (out / "candidates_enriched.json").write_text(json.dumps(rows, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    with open(out / "candidates_enriched.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["rank", "drug_id", "name", "combined_score", "pubmed_aneuploidy", "pubmed_tumour", "trials_rhabdomyosarcoma",
                         "trials_wilms", "trials_mva", "faers_signal_groups", "cpic_genes", "checkpoint_direction", "reversal_direction",
                         "reversal_ko_control_p", "molecule_type", "mw", "alogp", "psa", "ro5_violations",
                         "best_measured_pactivity", "n_adme_types", "routes", "robustness_count", "robustness_class"])
        for i, r in enumerate(rows, 1):
            layers = r["layers"]
            ev = layers["evidence_retrieved"]
            pm, tr = ev.get("pubmed_counts", {}), ev.get("trial_counts", {})
            rev = layers["transcriptomic_reversal"]
            pc = layers["physchem_affinity_adme"]
            writer.writerow([i, r["drug_id"], r["name"], r["combined_score"], pm.get("aneuploidy_cin", ""), pm.get("tumour_context", ""),
                             tr.get("rhabdomyosarcoma", ""), tr.get("wilms_tumor", ""), tr.get("mosaic_variegated_aneuploidy", ""),
                             ";".join(layers["faers_disproportionality"].get("signal_groups", [])),
                             ";".join(sorted({p["gene"] for p in layers["pharmacogenomics_cpic"]})),
                             ";".join(layers["checkpoint_direction"]), rev.get("direction", ""), rev.get("ko_control_p", ""),
                             pc.get("molecule_type", ""), pc.get("mw", ""), pc.get("alogp", ""), pc.get("psa", ""), pc.get("ro5_violations", ""),
                             pc.get("best_median_pactivity", ""), pc.get("n_adme_measurement_types", ""), ";".join(pc.get("routes", [])),
                             layers["robustness"]["robustness_count"], layers["robustness"]["class"]])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    rows = integrate(args.out)
    write(rows, args.out)
    print(len(rows), "enriched candidates written")


if __name__ == "__main__":
    main()
