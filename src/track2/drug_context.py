"""What a drug is approved for today and what it does to which gene (generated from data, never typed by hand).

- Approved indications: ChEMBL `drug_indication` rows with `max_phase_for_ind = 4` (MeSH headings).
- Target roles: for each module protein the drug is linked to, the gene symbol, the ChEMBL action type for curated targets (inhibits, activates, blocks...), or
  "binds (measured, pChEMBL x)" / "predicted binder" where ChEMBL gives no action. The submission gene whose module contains the protein is named.
- A one-sentence hypothesis built from a fixed template, so it carries no information that is not in the data.
"""

from __future__ import annotations

from typing import Any, Callable

from track2.http_cache import get_json

BASE = "https://www.ebi.ac.uk/chembl/api/data"
MAX_INDICATIONS = 8
VERB = {"INHIBITOR": "inhibits", "ANTAGONIST": "blocks", "BLOCKER": "blocks", "AGONIST": "activates", "PARTIAL AGONIST": "partially activates",
        "ACTIVATOR": "activates", "POSITIVE ALLOSTERIC MODULATOR": "positively modulates", "NEGATIVE ALLOSTERIC MODULATOR": "negatively modulates",
        "MODULATOR": "modulates", "STABILISER": "stabilises", "BINDING AGENT": "binds"}


def approved_indications(drug_id: str, get: Callable[[str], Any] = get_json) -> list[str]:
    heads: set[str] = set()
    offset = 0
    while True:
        data = get(f"{BASE}/drug_indication.json?molecule_chembl_id={drug_id}&limit=1000&offset={offset}")
        for r in data["drug_indications"]:
            if str(r.get("max_phase_for_ind")) in ("4", "4.0") and r.get("mesh_heading"):
                heads.add(r["mesh_heading"])
        if not data["page_meta"]["next"]:
            break
        offset += 1000
    return sorted(heads)[:MAX_INDICATIONS]


TUMOUR_WORDS = ("neoplasm", "carcinoma", "lymphoma", "leukemia", "myeloma", "sarcoma", "melanoma", "glioma", "tumor", "blastoma", "hodgkin", "cancer")


def is_tumour_indication(indications: list[str]) -> bool:
    """True when most approved indications are neoplasms (ChEMBL often lacks an ATC code for new antineoplastics)."""
    if not indications:
        return False
    hits = sum(1 for i in indications if any(w in i.lower() for w in TUMOUR_WORDS))
    return hits * 2 > len(indications)


def role_text(link: dict[str, Any]) -> str:
    if link["source"] == "curated":
        verb = VERB.get((link.get("action") or "").upper(), f"acts on ({(link.get('action') or 'action not stated').lower()})")
        return f"{verb} {link['gene']} (curated mechanism)"
    if link["source"] == "measured":
        single = " from a single measurement" if link.get("n_measurements") == 1 else ""
        return f"binds {link['gene']} (measured{single}, median pChEMBL {link['median_pchembl']}; ChEMBL gives no action type)"
    return f"is predicted to bind {link['gene']} (ligand similarity {link['similarity']}; not measured)"


def target_roles(links_by_module: dict[str, list[dict[str, Any]]], max_items: int = 6) -> list[dict[str, Any]]:
    """One entry per linked protein: the gene, the role text and the submission genes (modules) that contain it."""
    seen: dict[str, dict[str, Any]] = {}
    for module, links in links_by_module.items():
        for link in links:
            entry = seen.setdefault(link["uniprot"], {"gene": link["gene"], "role": role_text(link), "source": link["source"], "modules": []})
            entry["modules"].append(module.removeprefix("sub_"))
    order = {"curated": 0, "measured": 1, "predicted": 2}
    return sorted(seen.values(), key=lambda e: (order[e["source"]], e["gene"]))[:max_items]


def proximal_text(proximal: list[dict[str, Any]]) -> str:
    return "; ".join(f"{p['gene']} ({p['source']}{', pChEMBL ' + str(p['median_pchembl']) if p.get('median_pchembl') else ''}; {p['distance']} step(s) from the {p['module']} module)" for p in proximal[:3])


def hypothesis_sentence(name: str, roles: list[dict[str, Any]], indications: list[str], proximal: list[dict[str, Any]] | None = None) -> str:
    approved = f"approved today for {', '.join(i.lower() for i in indications[:3])}" if indications else "no approved indication listed in ChEMBL"
    if not roles and proximal:
        return (f"{name.title()} ({approved}) is not linked to a module protein directly; its closest targets are {proximal_text(proximal)}. "
                f"A hypothesis based on network closeness only; the direction of effect for this patient is not established.")
    if not roles:
        return f"{name.title()} ({approved}); no link to a module protein."
    first = roles[0]
    genes = sorted({m for r in roles for m in r["modules"]})
    extra = f" (and {len(roles) - 1} more module proteins)" if len(roles) > 1 else ""
    return (f"{name.title()} ({approved}) {first['role']}{extra}, a protein of the module of {', '.join(genes)}: a hypothesis that it acts on the pathway "
            f"disrupted by the submission's genes; the direction of effect for this patient is not established.")
