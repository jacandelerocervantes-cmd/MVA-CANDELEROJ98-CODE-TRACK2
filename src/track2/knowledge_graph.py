"""Line A: phenotype-guided knowledge graph (PrimeKG) with personalised PageRank (rule: `protocol/addendum_four_lines.md`, notes for line A).

PrimeKG (Chandak, Huang and Zitnik, Sci Data 2023; CC0) links diseases, HPO phenotypes, genes, pathways, drugs and more. The pretrained TxGNN model needs DGL with CUDA and is not used;
the fallback fixed in the addendum is personalised PageRank (restart probability 0.15) on the graph **without any drug-disease edge** (indication, contraindication, off-label use), so
that a disease can only reach a drug through genes, phenotypes, pathways and other shared neighbours, and **without the side-effect edges** (`drug_effect`: the drug causes the phenotype).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse

RESTART = 0.15
ITERATIONS, TOLERANCE = 100, 1e-9
DRUG_DISEASE = ("indication", "contraindication", "off-label use")
SIDE_EFFECT = ("drug_effect",)  # a drug linked to a phenotype by this relation causes it (amendment 2 of the addendum)


@dataclass
class Graph:
    node_type: np.ndarray            # (n,) str
    node_id: np.ndarray              # (n,) str, the source identifier (DrugBank id, HPO number, ...)
    node_name: np.ndarray            # (n,) str
    adj: sparse.csr_matrix           # symmetric, binary, without drug-disease edges
    indications: dict[int, set[int]]  # disease node -> drug nodes with an `indication` edge (ground truth, never part of `adj`)


def load_graph(nodes_csv: Path, edges_csv: Path) -> Graph:
    with open(nodes_csv, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        head = next(reader)
        rows = list(reader)
    index = np.array([int(r[head.index("node_index")]) for r in rows])
    if not np.array_equal(index, np.arange(len(rows))):
        raise ValueError("node_index is not 0..n-1 in file order")
    col = lambda name: np.array([r[head.index(name)] for r in rows])  # noqa: E731
    node_type, node_id, node_name = col("node_type"), col("node_id"), col("node_name")
    xs, ys = [], []
    indications: dict[int, set[int]] = {}
    with open(edges_csv, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        h = next(reader)
        ir, ix, iy = h.index("relation"), h.index("x_index"), h.index("y_index")
        for row in reader:
            x, y, rel = int(row[ix]), int(row[iy]), row[ir]
            if rel in DRUG_DISEASE:
                if rel == "indication":
                    drug, disease = (x, y) if node_type[x] == "drug" else (y, x)
                    indications.setdefault(disease, set()).add(drug)
                continue
            if rel in SIDE_EFFECT:
                continue
            xs.append(x)
            ys.append(y)
    n = len(rows)
    a = sparse.coo_matrix((np.ones(len(xs)), (xs, ys)), shape=(n, n)).tocsr()
    a = a.maximum(a.T)
    a.data[:] = 1.0
    return Graph(node_type, node_id, node_name, a.tocsr(), indications)


def personalised_pagerank(adj: sparse.csr_matrix, restart: np.ndarray, alpha: float = RESTART) -> np.ndarray:
    """Columns of `restart` (n x k) are restart distributions (each sums to 1). Returns the stationary scores (n x k)."""
    degree = np.asarray(adj.sum(axis=1)).ravel()
    inv = np.where(degree > 0, 1.0 / np.maximum(degree, 1), 0.0)
    walk = sparse.diags(inv).dot(adj).T.tocsr()          # column-stochastic transition: p_next = walk @ p
    p = restart.copy()
    for _ in range(ITERATIONS):
        nxt = (1 - alpha) * walk.dot(p) + alpha * restart
        done = np.abs(nxt - p).sum(axis=0).max() < TOLERANCE
        p = nxt
        if done:
            break
    return p


def restart_vector(n: int, nodes: dict[int, float]) -> np.ndarray:
    v = np.zeros(n)
    for node, weight in nodes.items():
        v[node] += weight
    return v / v.sum()
