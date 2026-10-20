"""Network proximity between a drug's targets and a disease module (d_closest + z-score).

Reimplementation of the method of Guney et al., Nat Commun 2016 (DOI 10.1038/ncomms10331), with the
extension of Cheng et al., Nat Commun 2018 (DOI 10.1038/s41467-018-05116-5):
- d_closest(S, T): average, over the targets t of T, of the minimum distance from t to any protein of S.
- z = (d - mu) / sigma against R random sets with the same number of nodes and the same degree
  (degree bins with at least 100 nodes); both the module and the targets are randomised.
Proximity does NOT predict the direction of the effect (benefit or harm): whoever uses it must review that separately.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import dijkstra

MIN_BIN_SIZE = 100


def degree_bins(degree: np.ndarray, min_size: int = MIN_BIN_SIZE) -> tuple[np.ndarray, list[np.ndarray]]:
    """Group nodes by degree into bins of at least `min_size` nodes without splitting a single degree."""
    order = np.argsort(degree, kind="stable")
    bins: list[np.ndarray] = []
    start = 0
    n = len(order)
    while start < n:
        end = min(start + min_size, n)
        while end < n and degree[order[end]] == degree[order[end - 1]]:
            end += 1
        if n - end < min_size:  # the last bin does not fall below the minimum
            end = n
        bins.append(order[start:end])
        start = end
    bin_of = np.empty(n, dtype=np.int32)
    for i, members in enumerate(bins):
        bin_of[members] = i
    return bin_of, bins


def min_distance_to_set(adjacency: sparse.csr_matrix, sources: np.ndarray) -> np.ndarray:
    """Distance (in hops) from each node to the nearest node of `sources`."""
    return dijkstra(adjacency, directed=False, unweighted=True, indices=sources, min_only=True).astype(np.float32)


def d_closest(distance_to_module: np.ndarray, targets: np.ndarray) -> float:
    return float(distance_to_module[targets].mean())


def _sample_matched(rng: np.random.Generator, nodes: np.ndarray, bin_of: np.ndarray, bins: list[np.ndarray], repeats: int) -> np.ndarray:
    """(repeats, len(nodes)): for each node, a random one from the same degree bin."""
    out = np.empty((repeats, len(nodes)), dtype=np.int64)
    for j, node in enumerate(nodes):
        members = bins[bin_of[node]]
        out[:, j] = members[rng.integers(0, len(members), size=repeats)]
    return out


@dataclass
class NullModel:
    """Distances to R random modules of equal size and degree; computed once per module."""

    distances: np.ndarray        # (R, n_nodes) float32
    bin_of: np.ndarray
    bins: list[np.ndarray]
    repeats: int
    seed: int


def build_null(adjacency: sparse.csr_matrix, degree: np.ndarray, module: np.ndarray, repeats: int = 1000, seed: int = 20261006) -> NullModel:
    rng = np.random.default_rng(seed)
    bin_of, bins = degree_bins(degree)
    random_modules = _sample_matched(rng, module, bin_of, bins, repeats)
    distances = np.stack([min_distance_to_set(adjacency, np.unique(random_modules[r])) for r in range(repeats)])
    return NullModel(distances=distances, bin_of=bin_of, bins=bins, repeats=repeats, seed=seed)


def proximity(distance_to_module: np.ndarray, targets: np.ndarray, null: NullModel, seed_offset: int = 0) -> dict[str, float]:
    """Observed d, z-score and empirical p of a set of targets against the module."""
    rng = np.random.default_rng(null.seed + 1 + seed_offset)
    d_obs = d_closest(distance_to_module, targets)
    sampled = _sample_matched(rng, targets, null.bin_of, null.bins, null.repeats)
    d_null = null.distances[np.arange(null.repeats)[:, None], sampled].mean(axis=1)
    sigma = float(d_null.std())
    z = (d_obs - float(d_null.mean())) / sigma if sigma > 0 else 0.0
    return {"d": d_obs, "mu": float(d_null.mean()), "sigma": sigma, "z": z, "p_empirical": float((d_null <= d_obs).mean())}
