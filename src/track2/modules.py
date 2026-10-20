"""Disease module (seed genes) of a mechanism, built from its configuration.

A mechanism is declared in `config/mechanisms*/*.json` (data, not code): Reactome pathways and additional
genes (for example the disease genes in OMIM). Here they are resolved to interactome nodes and the minimum
module size is checked (10 genes, Cheng et al. 2018).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np

from track2.http_cache import get_json
from track2.string_net import Interactome

REACTOME = "https://reactome.org/ContentService"
MIN_MODULE_GENES = 10


def reactome_version() -> str:
    from track2.http_cache import get_bytes
    return get_bytes(f"{REACTOME}/data/database/version", accept="text/plain").decode().strip()


def reactome_uniprot(pathway_id: str) -> list[str]:
    entities = get_json(f"{REACTOME}/data/participants/{pathway_id}/referenceEntities")
    return sorted({e["identifier"] for e in entities if e.get("databaseName") == "UniProt"})


def load_mechanism(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_module(mechanism: dict[str, Any], interactome: Interactome) -> dict[str, Any]:
    """Module nodes plus a record of which seeds were resolved and which were not."""
    seeds = mechanism["seeds"]
    uniprots: list[str] = []
    for pathway in seeds.get("reactome_pathways", []):
        uniprots.extend(reactome_uniprot(pathway))
    nodes = set(interactome.nodes_for_uniprot(uniprots)) | set(interactome.nodes_for_symbols(seeds.get("genes", [])))
    if seeds.get("exclude_symbol_regex"):  # sensitivity analysis: drop generic complexes (see protocol/addendum_core.md)
        rx = re.compile(seeds["exclude_symbol_regex"])
        drop = {n for sym, n in interactome.symbol_to_node.items() if rx.search(sym)}
        nodes -= drop
    missing_symbols = [g for g in seeds.get("genes", []) if g.upper() not in interactome.symbol_to_node]
    arr = np.array(sorted(nodes), dtype=np.int64)
    if len(arr) < MIN_MODULE_GENES:
        raise ValueError(f"Module '{mechanism['id']}' has {len(arr)} nodes: below the minimum of {MIN_MODULE_GENES}.")
    return {
        "id": mechanism["id"], "nodes": arr, "n_seed_uniprot": len(set(uniprots)),
        "n_in_interactome": int(len(arr)), "symbols_not_in_interactome": missing_symbols,
    }
