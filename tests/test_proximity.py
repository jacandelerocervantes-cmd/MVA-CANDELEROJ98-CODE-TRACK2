"""Network proximity tests with a synthetic graph (no network access and no patient data)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from track2.proximity import build_null, d_closest, degree_bins, min_distance_to_set, proximity


def _two_communities(n: int = 600, seed: int = 1) -> sparse.csr_matrix:
    """Two dense communities (0..n/2-1 and n/2..n-1) joined by a few links."""
    rng = np.random.default_rng(seed)
    half = n // 2
    rows, cols = [], []
    for start in (0, half):
        for i in range(start, start + half):
            for j in rng.choice(np.arange(start, start + half), size=6, replace=False):
                if i != j:
                    rows.append(i)
                    cols.append(int(j))
    for k in range(3):  # bridges
        rows.append(k)
        cols.append(half + k)
    a = sparse.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n, n)).tocsr()
    return ((a + a.T) > 0).astype(np.int8).tocsr()


def test_min_distance_is_multisource_and_zero_on_sources() -> None:
    adjacency = sparse.csr_matrix(np.array([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]], dtype=np.int8))
    distance = min_distance_to_set(adjacency, np.array([0, 3]))
    assert distance.tolist() == [0.0, 1.0, 1.0, 0.0]
    assert d_closest(distance, np.array([1, 2])) == 1.0


def test_degree_bins_have_minimum_size_and_keep_equal_degrees_together() -> None:
    degree = np.array([1] * 150 + [2] * 30 + [3] * 130 + [9] * 5)
    bin_of, bins = degree_bins(degree, min_size=100)
    assert all(len(b) >= 100 for b in bins)
    for value in (1, 2, 3, 9):
        assert len(set(bin_of[degree == value])) == 1


def test_targets_inside_the_module_get_strongly_negative_z_and_far_ones_do_not() -> None:
    adjacency = _two_communities()
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    module = np.arange(0, 20)                      # seeds in community A
    null = build_null(adjacency, degree, module, repeats=200, seed=7)
    dist = min_distance_to_set(adjacency, module)
    near = proximity(dist, np.array([2, 5, 11]), null)
    far = proximity(dist, np.array([450, 480, 520]), null)
    assert near["z"] < -2 and far["z"] > near["z"]
    assert near["d"] < far["d"]


def test_proximity_is_deterministic_for_a_fixed_seed() -> None:
    adjacency = _two_communities()
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    module = np.arange(0, 15)
    first = proximity(min_distance_to_set(adjacency, module), np.array([3, 40]), build_null(adjacency, degree, module, 100, 3))
    second = proximity(min_distance_to_set(adjacency, module), np.array([3, 40]), build_null(adjacency, degree, module, 100, 3))
    assert first == pytest.approx(second)
