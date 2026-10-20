"""Universe of approved drugs from ChEMBL, with mechanism, targets (UniProt) and references.

Approved = `max_phase == 4` and not withdrawn. The parent molecule is used (salts and forms are
deduplicated). Targets come from the ChEMBL mechanism table (curated action), not from raw activities.
Nothing is fixed by hand: if ChEMBL changes version, the version is recorded alongside the data.
"""

from __future__ import annotations

from typing import Any

from track2.http_cache import get_json

BASE = "https://www.ebi.ac.uk/chembl/api/data"
PAGE = 1000
BATCH = 40


def chembl_version() -> str:
    return str(get_json(f"{BASE}/status.json", use_cache=False)["chembl_db_version"])


def _paged(endpoint: str, key: str, query: str) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    offset = 0
    while True:
        data = get_json(f"{BASE}/{endpoint}.json?{query}&limit={PAGE}&offset={offset}")
        out.extend(data[key])
        if not data["page_meta"]["next"]:
            return out
        offset += PAGE


def fetch_approved_molecules() -> list[dict[str, Any]]:
    """Parent molecules with max_phase = 4, not withdrawn."""
    raw = _paged("molecule", "molecules", "max_phase=4")
    parents: dict[str, dict[str, Any]] = {}
    for m in raw:
        hierarchy = m.get("molecule_hierarchy") or {}
        pid = hierarchy.get("parent_chembl_id") or m["molecule_chembl_id"]
        if m["molecule_chembl_id"] != pid or m.get("withdrawn_flag"):
            continue
        structures = m.get("molecule_structures") or {}
        parents[pid] = {
            "drug_id": pid,
            "name": m.get("pref_name"),
            "molecule_type": m.get("molecule_type"),
            "first_approval": m.get("first_approval"),
            "oral": m.get("oral"),
            "parenteral": m.get("parenteral"),
            "topical": m.get("topical"),
            "black_box_warning": m.get("black_box_warning"),
            "atc": sorted({a for a in (m.get("atc_classifications") or [])}),
            "smiles": structures.get("canonical_smiles"),
            "inchikey": structures.get("standard_inchi_key"),
        }
    return sorted(parents.values(), key=lambda d: d["drug_id"])


def fetch_mechanisms(drug_ids: list[str]) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {d: [] for d in drug_ids}
    for i in range(0, len(drug_ids), BATCH):
        chunk = drug_ids[i:i + BATCH]
        for m in _paged("mechanism", "mechanisms", "molecule_chembl_id__in=" + ",".join(chunk)):
            out.setdefault(m["molecule_chembl_id"], []).append({
                "target_id": m.get("target_chembl_id"),
                "action": m.get("action_type"),
                "mechanism": m.get("mechanism_of_action"),
                "direct": m.get("direct_interaction"),
                "refs": [{"type": r["ref_type"], "id": r["ref_id"]} for r in (m.get("mechanism_refs") or [])],
            })
    return out


def fetch_target_uniprot(target_ids: list[str]) -> dict[str, dict[str, Any]]:
    """target_chembl_id -> type and human UniProt accessions of its components."""
    out: dict[str, dict[str, Any]] = {}
    unique = sorted(set(target_ids))
    for i in range(0, len(unique), BATCH):
        chunk = unique[i:i + BATCH]
        for t in _paged("target", "targets", "target_chembl_id__in=" + ",".join(chunk)):
            accessions = sorted({c["accession"] for c in t.get("target_components", []) if c.get("accession")})
            out[t["target_chembl_id"]] = {
                "name": t.get("pref_name"), "type": t.get("target_type"),
                "organism": t.get("organism"), "uniprot": accessions,
            }
    return out


def build_universe() -> dict[str, Any]:
    """Full universe: approved drugs + mechanisms + human targets, with the ChEMBL version."""
    drugs = fetch_approved_molecules()
    mechanisms = fetch_mechanisms([d["drug_id"] for d in drugs])
    targets = fetch_target_uniprot([m["target_id"] for ms in mechanisms.values() for m in ms if m["target_id"]])
    for d in drugs:
        d["mechanisms"] = mechanisms.get(d["drug_id"], [])
        human = []
        for m in d["mechanisms"]:
            t = targets.get(m["target_id"] or "", {})
            if t.get("organism") == "Homo sapiens":
                human.extend(t["uniprot"])
        d["human_targets"] = sorted(set(human))
    return {"chembl_version": chembl_version(), "n_drugs": len(drugs), "drugs": drugs, "targets": targets}
