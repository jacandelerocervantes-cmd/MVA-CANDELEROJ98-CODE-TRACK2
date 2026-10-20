"""Measured binding, physicochemical properties and measured ADME values of the shortlisted drugs (ChEMBL).

    python -m track2.physchem --candidates out/candidates.json --universe data/universe.json --out out

Everything here is data that someone measured or computed and deposited in ChEMBL; nothing is predicted by this
repository (no docking, no affinity model: see the report). The layer is context for a human reviewer and does
not change the ranking. Biologics (antibodies, antibody-drug conjugates) have no small-molecule properties and
are reported as such, not filled with invented values.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path
from typing import Any, Callable

from track2.http_cache import get_json

BASE = "https://www.ebi.ac.uk/chembl/api/data"
AFFINITY_TYPES = ("IC50", "Ki", "Kd", "EC50", "Potency")
ADME_TYPES = ("Solubility", "LogD", "LogP", "Clearance", "T1/2", "Papp", "Fu", "Bioavailability", "F")
PROPERTY_KEYS = ("full_mwt", "alogp", "cx_logp", "cx_logd", "psa", "hba", "hbd", "rtb", "aromatic_rings", "heavy_atoms",
                 "num_ro5_violations", "qed_weighted", "molecular_species", "cx_most_apka", "cx_most_bpka")


def properties(molecule: dict[str, Any]) -> dict[str, Any]:
    props = molecule.get("molecule_properties") or {}
    return {k: props.get(k) for k in PROPERTY_KEYS if props.get(k) is not None}


def affinity_summary(activities: list[dict[str, Any]]) -> dict[str, Any]:
    """Median pActivity (-log10 of the nM value / 1e9) per target, from exact (=) nanomolar measurements only."""
    per_target: dict[str, list[float]] = {}
    for a in activities:
        try:
            value = float(a["standard_value"])
        except (KeyError, TypeError, ValueError):
            continue
        if a.get("standard_relation") != "=" or a.get("standard_units") != "nM" or value <= 0 or not a.get("target_chembl_id"):
            continue
        per_target.setdefault(a["target_chembl_id"], []).append(9.0 - math.log10(value))
    return {t: {"n": len(v), "median_pactivity": round(statistics.median(v), 2)} for t, v in sorted(per_target.items())}


def adme_summary(activities: list[dict[str, Any]]) -> dict[str, Any]:
    """Per measurement type: number of values, median and units (values with a different unit are not mixed)."""
    groups: dict[tuple[str, str], list[float]] = {}
    for a in activities:
        try:
            value = float(a["standard_value"])
        except (KeyError, TypeError, ValueError):
            continue
        groups.setdefault((a.get("standard_type") or "?", a.get("standard_units") or ""), []).append(value)
    return {f"{t} ({u or 'no unit'})": {"n": len(v), "median": round(statistics.median(v), 4)} for (t, u), v in sorted(groups.items())}


def profile(drug: dict[str, Any], get: Callable[[str], Any] = get_json) -> dict[str, Any]:
    drug_id = drug["drug_id"]
    molecule = get(f"{BASE}/molecule/{drug_id}.json")
    targets = sorted({m["target_id"] for m in drug.get("mechanisms", []) if m.get("target_id")})
    affinity: dict[str, Any] = {}
    if targets:
        url = (f"{BASE}/activity.json?molecule_chembl_id={drug_id}&target_chembl_id__in={','.join(targets)}"
               f"&standard_type__in={','.join(AFFINITY_TYPES)}&limit=1000")
        affinity = affinity_summary(get(url)["activities"])
    adme_url = (f"{BASE}/activity.json?molecule_chembl_id={drug_id}&assay_type=A&standard_type__in="
                f"{','.join(t.replace('/', '%2F') for t in ADME_TYPES)}&limit=1000")
    adme = adme_summary(get(adme_url)["activities"])
    props = properties(molecule)
    return {
        "drug_id": drug_id, "name": drug["name"], "molecule_type": molecule.get("molecule_type"),
        "routes": {k: bool(molecule.get(k)) for k in ("oral", "parenteral", "topical")},
        "properties": props, "has_small_molecule_properties": bool(props),
        "measured_affinity_to_mechanism_targets": affinity, "measured_adme": adme,
        "lipinski_violations": props.get("num_ro5_violations"),
    }


def run(candidates: list[dict[str, Any]], universe: dict[str, Any], get: Callable[[str], Any] = get_json) -> dict[str, Any]:
    drugs = {d["drug_id"]: d for d in universe["drugs"]}
    return {c["drug_id"]: profile(drugs[c["drug_id"]], get) for c in candidates}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=Path("out/candidates.json"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    candidates = json.loads(args.candidates.read_text(encoding="utf-8"))
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    result = run(candidates, universe)
    (args.out / "physchem_adme.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    n_props = sum(1 for r in result.values() if r["has_small_molecule_properties"])
    n_aff = sum(1 for r in result.values() if r["measured_affinity_to_mechanism_targets"])
    n_adme = sum(1 for r in result.values() if r["measured_adme"])
    print(f"{len(result)} drugs: {n_props} with properties, {n_aff} with measured affinity, {n_adme} with measured ADME")


if __name__ == "__main__":
    main()
