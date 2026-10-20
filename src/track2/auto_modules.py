"""Build a disease module for a gene automatically, with a fixed rule (no hand-picked pathways).

Rule (pre-registered in `protocol/addendum_submission.md`):
1. Resolve the gene symbol to a reviewed human UniProt accession (UniProt search).
2. Take every Reactome pathway that Reactome maps to that accession (human). Keep those with 10 to 300 UniProt
   participants. The module is the union of their participants plus the gene itself.
3. If no pathway qualifies, use the gene plus its STRING physical neighbours (combined score >= 700); if that gives fewer
   than 10 proteins, add the neighbours of those neighbours (second order; amendment 1 of the addendum).
4. A module needs at least 10 proteins in the network; otherwise the gene is reported as "not modellable".
"""

from __future__ import annotations

from typing import Any

from track2.http_cache import get_json
from track2.modules import MIN_MODULE_GENES, reactome_uniprot
from track2.string_net import Interactome

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
REACTOME_MAP = "https://reactome.org/ContentService/data/mapping/UniProt/{acc}/pathways?species=9606"
MIN_PATHWAY, MAX_PATHWAY = 10, 300


def uniprot_accession(symbol: str) -> str | None:
    query = f"gene_exact:{symbol}+AND+organism_id:9606+AND+reviewed:true"
    results = get_json(f"{UNIPROT_SEARCH}?query={query}&fields=accession&format=json&size=1")["results"]
    return results[0]["primaryAccession"] if results else None


def pathways_for(accession: str) -> list[dict[str, Any]]:
    """Reactome pathways mapped to the accession, with their number of UniProt participants."""
    out = []
    try:
        mapped = get_json(REACTOME_MAP.format(acc=accession))
    except LookupError:  # Reactome answers 404 when the protein is in no pathway
        return []
    for p in mapped:
        out.append({"id": p["stId"], "name": p["displayName"], "n_participants": len(reactome_uniprot(p["stId"]))})
    return sorted(out, key=lambda p: p["id"])


def _network_module(gene: str, interactome: Interactome) -> tuple[dict[str, Any], str]:
    """Gene plus its STRING physical neighbours: first-order, then second-order if first-order has fewer than 10 proteins."""
    node = interactome.symbol_to_node.get(gene.upper())
    symbols: dict[int, str] = {}
    for sym, n in interactome.symbol_to_node.items():
        symbols.setdefault(n, sym)
    if node is None:
        return {"reactome_pathways": [], "genes": [gene]}, "not_in_network"
    first = {node} | set(int(n) for n in interactome.adjacency[node].indices)
    members, rule = first, "string_physical_neighbours"
    if len(first) < MIN_MODULE_GENES:
        members, rule = set(first), "string_physical_neighbours_2nd_order"
        for n in first:
            members |= set(int(x) for x in interactome.adjacency[n].indices)
    return {"reactome_pathways": [], "genes": sorted(symbols[n] for n in members if n in symbols)}, rule


def build_mechanism(gene: str, weight: float, interactome: Interactome, evidence: dict[str, Any]) -> dict[str, Any]:
    accession = uniprot_accession(gene)
    chosen = [p for p in (pathways_for(accession) if accession else []) if MIN_PATHWAY <= p["n_participants"] <= MAX_PATHWAY]
    seeds: dict[str, Any] = {"reactome_pathways": [p["id"] for p in chosen], "genes": [gene]}
    rule = "reactome_pathways_of_the_gene"
    if not chosen:
        seeds, rule = _network_module(gene, interactome)
    return {
        "id": f"sub_{gene}", "gene": gene, "uniprot": accession, "status": "mechanism_candidate", "weight": weight,
        "weight_note": "Best epcr of the gene's rows in the Track 1 submission (the submission's own ranking value).",
        "module_rule": rule, "pathways": chosen, "seeds": seeds, "evidence_from_submission": evidence,
    }


def module_size(mechanism: dict[str, Any], interactome: Interactome) -> int:
    from track2.modules import build_module
    try:
        return int(build_module(mechanism, interactome)["n_in_interactome"])
    except ValueError:
        return 0


def modellable(mechanism: dict[str, Any], interactome: Interactome) -> bool:
    return module_size(mechanism, interactome) >= MIN_MODULE_GENES
