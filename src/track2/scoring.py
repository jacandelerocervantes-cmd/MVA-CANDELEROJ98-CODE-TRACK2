"""Score of each approved drug against a disease module: proximity, action and risk flags.

No drug or gene is hard-coded here. Proximity does not give the direction of the effect, so each drug carries
`direction_review = REQUIRES_REVIEW` and flags computed from the ChEMBL action and the ATC class.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from track2.proximity import NullModel, min_distance_to_set, proximity
from track2.string_net import Interactome

INHIBITORY_ACTIONS = {
    "INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR",
    "INVERSE AGONIST", "DISRUPTING AGENT", "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "CROSS-LINKING AGENT",
}
STABILISING_ACTIONS = {"STABILISER"}
PROXIMITY_Z_PRIMARY = -2.0


def _flags(drug: dict[str, Any], module_uniprot: set[str]) -> dict[str, Any]:
    in_module = []
    for mechanism in drug["mechanisms"]:
        in_module.append(mechanism["action"])
    direct_hits = sorted(set(drug["human_targets"]) & module_uniprot)
    inhibits_module = bool(direct_hits) and any(
        (m["action"] or "") in INHIBITORY_ACTIONS for m in drug["mechanisms"]
    )
    flags = []
    if any(code.startswith("L01") for code in drug["atc"]):
        flags.append("ANTINEOPLASTIC_AGENT")
    if inhibits_module:
        flags.append("INHIBITS_MODULE_COMPONENT")
    if (drug.get("black_box_warning") or 0) == 1:
        flags.append("BOXED_WARNING")
    if any((m["action"] or "") in STABILISING_ACTIONS for m in drug["mechanisms"]):
        flags.append("MICROTUBULE_OR_PROTEIN_STABILISER")
    return {
        "actions": sorted({a for a in in_module if a}),
        "direct_module_targets": direct_hits,
        "risk_flags": flags,
        "direction_review": "REQUIRES_REVIEW",
    }


def score_drugs(
    drugs: list[dict[str, Any]], interactome: Interactome, module_nodes: np.ndarray, module_uniprot: set[str], null: NullModel,
) -> list[dict[str, Any]]:
    """One row per drug with targets inside the interactome; sorted by z (more negative = closer)."""
    distance = min_distance_to_set(interactome.adjacency, module_nodes)
    rows = []
    for index, drug in enumerate(drugs):
        nodes = np.array(interactome.nodes_for_uniprot(drug["human_targets"]), dtype=np.int64)
        if nodes.size == 0:
            continue
        stats = proximity(distance, nodes, null, seed_offset=index)
        rows.append({
            "drug_id": drug["drug_id"], "name": drug["name"], "atc": drug["atc"], "first_approval": drug["first_approval"],
            "n_targets_in_network": int(nodes.size), "targets": drug["human_targets"],
            "proximity": stats, **_flags(drug, module_uniprot),
        })
    rows.sort(key=lambda r: (r["proximity"]["z"], r["drug_id"]))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
        row["passes_primary_threshold"] = row["proximity"]["z"] <= PROXIMITY_Z_PRIMARY
    return rows
