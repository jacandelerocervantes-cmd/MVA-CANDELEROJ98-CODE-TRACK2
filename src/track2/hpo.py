"""Read the Human Phenotype Ontology (hp.obo) and the gene-to-phenotype annotations (genes_to_phenotype.txt).

Both files come from one HPO release (`track2.fetch_sources`, source `hpo_ontology`). Only the structure needed by the
phenotype module is read: term ids, alternative ids, `is_a` parents, and the gene symbol of each annotation.
"""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path


def parse_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Return (term -> direct parents, alternative id -> current id). Obsolete terms are skipped."""
    parents: dict[str, set[str]] = {}
    alt: dict[str, str] = {}
    current: str | None = None
    obsolete = False
    pending: list[str] = []
    alts: list[str] = []

    def close() -> None:
        if current and not obsolete:
            parents[current] = set(pending)
            for a in alts:
                alt[a] = current

    in_term = False
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line == "[Term]":
            if in_term:
                close()
            in_term, current, obsolete, pending, alts = True, None, False, [], []
        elif line.startswith("[") and line.endswith("]"):  # [Typedef] and the like
            if in_term:
                close()
            in_term, current = False, None
        elif in_term and line.startswith("id: "):
            current = line[4:].strip()
        elif in_term and line.startswith("is_a: "):
            pending.append(line[6:].split("!")[0].strip())
        elif in_term and line.startswith("alt_id: "):
            alts.append(line[8:].strip())
        elif in_term and line == "is_obsolete: true":
            obsolete = True
    if in_term:
        close()
    return parents, alt


def ancestors(term: str, parents: dict[str, set[str]]) -> set[str]:
    """The term and all its ancestors."""
    seen, stack = set(), [term]
    while stack:
        t = stack.pop()
        if t not in seen:
            seen.add(t)
            stack.extend(parents.get(t, ()))
    return seen


def parse_gene_annotations(path: Path, alt: dict[str, str]) -> dict[str, set[str]]:
    """gene symbol -> set of directly annotated HPO ids (alternative ids replaced by the current id)."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].lstrip("#").split("\t")
    i_gene, i_hpo = header.index("gene_symbol"), header.index("hpo_id")
    out: dict[str, set[str]] = defaultdict(set)
    for line in lines[1:]:
        cols = line.split("\t")
        if len(cols) > max(i_gene, i_hpo) and cols[i_gene] not in ("", "-"):  # "-" marks an annotation with no gene symbol
            out[cols[i_gene]].add(alt.get(cols[i_hpo], cols[i_hpo]))
    return dict(out)


def propagate(annotations: dict[str, set[str]], parents: dict[str, set[str]]) -> dict[str, set[str]]:
    """Add the ancestors of every annotated term: a gene annotated to a term is also annotated to its parents."""
    cache: dict[str, set[str]] = {}
    out = {}
    for gene, terms in annotations.items():
        full: set[str] = set()
        for t in terms:
            if t not in cache:
                cache[t] = ancestors(t, parents)
            full |= cache[t]
        out[gene] = full
    return out


def information_content(propagated: dict[str, set[str]]) -> dict[str, float]:
    """IC(t) = -ln(fraction of annotated genes carrying t, after propagation). Rare terms weigh more."""
    n = len(propagated)
    count: dict[str, int] = defaultdict(int)
    for terms in propagated.values():
        for t in terms:
            count[t] += 1
    return {t: -math.log(c / n) for t, c in count.items()}
