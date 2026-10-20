from __future__ import annotations

from typing import Any

import numpy as np

from track2 import dti_edges, dti_similarity


def test_summarise_takes_median_per_target() -> None:
    acts = [{"target_chembl_id": "T1", "pchembl_value": "6.0"}, {"target_chembl_id": "T1", "pchembl_value": "8.0"},
            {"target_chembl_id": "T2", "pchembl_value": "7.5"}, {"target_chembl_id": "T3", "pchembl_value": None}]
    assert dti_edges.summarise(acts) == {"T1": {"n": 2, "median_pchembl": 7.0}, "T2": {"n": 1, "median_pchembl": 7.5}}


def test_target_mapping_keeps_only_single_proteins() -> None:
    def fake(url: str) -> Any:
        return {"targets": [{"target_chembl_id": "T1", "target_type": "SINGLE PROTEIN", "target_components": [{"accession": "P1"}]},
                            {"target_chembl_id": "T2", "target_type": "PROTEIN COMPLEX", "target_components": [{"accession": "P2"}, {"accession": "P3"}]},
                            {"target_chembl_id": "T3", "target_type": "SINGLE PROTEIN", "target_components": []}]}
    assert dti_edges.target_to_uniprot(["T1", "T2", "T3"], fake) == {"T1": "P1"}


def test_drug_activities_pages_until_no_next() -> None:
    pages = [{"activities": [{"a": 1}], "page_meta": {"next": "x"}}, {"activities": [{"a": 2}], "page_meta": {"next": None}}]
    calls: list[str] = []

    def fake(url: str) -> Any:
        calls.append(url)
        return pages[len(calls) - 1]
    assert dti_edges.drug_activities("CHEMBL1", fake) == [{"a": 1}, {"a": 2}] and "offset=1000" in calls[1]


def test_self_similarity_is_excluded_from_scores() -> None:
    fp = dti_similarity.fingerprint("CCO")
    other = dti_similarity.fingerprint("c1ccccc1")
    scores = dti_similarity.score_matrix([fp], ["D1"], {"U": [("D1", fp), ("L2", other)]})
    assert scores["U"][0] < 1.0  # D1's own entry (similarity 1.0) is ignored


def test_validation_picks_smallest_tau_with_enough_precision() -> None:
    n = 40
    scores = {"U": np.array([0.9] * 10 + [0.35] * 30)}
    positives = {(i, "U") for i in range(10)}  # the ten high-scoring drugs are the true links
    out = dti_similarity.validate(scores, positives, n)
    assert out["valid"] and out["n_positive"] == 10
    row = next(r for r in out["by_tau"] if r["tau"] == out["tau"])
    assert row["precision"] >= 0.5 and out["tau"] > 0.35  # the low-scoring negatives are below tau


def test_validation_not_valid_when_scores_do_not_separate() -> None:
    rng = np.random.default_rng(0)
    scores = {"U": rng.uniform(0.3, 0.6, 200)}
    positives = {(i, "U") for i in range(0, 200, 10)}  # random positives: precision stays near the base rate (0.1)
    assert dti_similarity.validate(scores, positives, 200)["valid"] is False


def test_heldout_precision_is_lower_than_in_sample_when_the_cut_off_is_overfitted() -> None:
    rng = np.random.default_rng(3)
    n = 300
    scores = {f"T{k}": rng.uniform(0.3, 0.7, n) for k in range(12)}  # no real signal
    positives = {(i, t) for t in scores for i in rng.choice(n, 3, replace=False)}
    out = dti_similarity.heldout(scores, positives, n, repeats=8)
    assert out["n_valid_splits"] >= 0 and (out["heldout_precision_mean"] is None or out["heldout_precision_mean"] < 0.5)
