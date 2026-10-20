"""Physical protein interaction network (STRING v12.0, Homo sapiens) as a sparse graph.

The *physical* subnetwork is used (the interactomes of the proximity literature are physical) with a
configurable combined-score threshold (default 700, high confidence); the sensitivity analysis varies the threshold.
The nodes are proteins of the largest connected component. The UniProt -> node and symbol -> node mappings come
from the alias file of the same version.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from track2.http_cache import get_file

STRING_VERSION = "v12.0"
LINKS_URL = f"https://stringdb-downloads.org/download/protein.physical.links.{STRING_VERSION}/9606.protein.physical.links.{STRING_VERSION}.txt.gz"
ALIASES_URL = f"https://stringdb-downloads.org/download/protein.aliases.{STRING_VERSION}/9606.protein.aliases.{STRING_VERSION}.txt.gz"
UNIPROT_SOURCES = {"UniProt_AC", "Ensembl_UniProt"}
SYMBOL_SOURCES = {"BioMart_HUGO", "Ensembl_HGNC_symbol"}


@dataclass
class Interactome:
    proteins: list[str]                      # ENSP, one element per node
    adjacency: sparse.csr_matrix
    degree: np.ndarray
    uniprot_to_node: dict[str, int]
    symbol_to_node: dict[str, int]
    min_score: int
    string_version: str

    def nodes_for_uniprot(self, accessions: list[str]) -> list[int]:
        return sorted({self.uniprot_to_node[a] for a in accessions if a in self.uniprot_to_node})

    def nodes_for_symbols(self, symbols: list[str]) -> list[int]:
        return sorted({self.symbol_to_node[s] for s in symbols if s in self.symbol_to_node})


def _read_aliases(path: Path, wanted: set[str]) -> tuple[dict[str, str], dict[str, str]]:
    uni: dict[str, str] = {}
    sym: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3 or parts[0] not in wanted:
                continue
            if parts[2] in UNIPROT_SOURCES:
                uni.setdefault(parts[1], parts[0])
            elif parts[2] in SYMBOL_SOURCES:
                sym.setdefault(parts[1].upper(), parts[0])
    return uni, sym


def load_interactome(min_score: int = 700) -> Interactome:
    links = get_file(LINKS_URL, f"9606.protein.physical.links.{STRING_VERSION}.txt.gz")
    aliases = get_file(ALIASES_URL, f"9606.protein.aliases.{STRING_VERSION}.txt.gz")
    index: dict[str, int] = {}
    rows: list[int] = []
    cols: list[int] = []
    with gzip.open(links, "rt", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            a, b, score = line.split()
            if int(score) < min_score or a == b:
                continue
            ia = index.setdefault(a, len(index))
            ib = index.setdefault(b, len(index))
            rows.append(ia)
            cols.append(ib)
    n = len(index)
    adjacency = sparse.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n, n)).tocsr()
    adjacency = ((adjacency + adjacency.T) > 0).astype(np.int8).tocsr()
    _, labels = connected_components(adjacency, directed=False)
    giant = np.flatnonzero(labels == np.bincount(labels).argmax())
    remap = {int(old): new for new, old in enumerate(giant)}
    adjacency = adjacency[giant][:, giant].tocsr()
    proteins = [""] * len(giant)
    for protein, old in index.items():
        if old in remap:
            proteins[remap[old]] = protein
    uni, sym = _read_aliases(aliases, set(proteins))
    protein_to_node = {p: i for i, p in enumerate(proteins)}
    return Interactome(
        proteins=proteins, adjacency=adjacency, degree=np.asarray(adjacency.sum(axis=1)).ravel(),
        uniprot_to_node={u: protein_to_node[p] for u, p in uni.items()},
        symbol_to_node={s: protein_to_node[p] for s, p in sym.items()},
        min_score=min_score, string_version=STRING_VERSION,
    )
