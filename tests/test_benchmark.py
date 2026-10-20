"""Benchmark tests: offline, with a mocked Open Targets client and a synthetic network."""

from __future__ import annotations

import numpy as np

from track2 import ot_client
from track2.benchmark_stats import DiseaseData, auc_from_scores, bootstrap_mean_ci, drug_bins, percentile_ranks, permutation_test
from track2.benchmark import decide, select_diseases


class FakeInteractome:
    def __init__(self, symbols: list[str], uniprots: list[str]) -> None:
        self.s = symbols
        self.u = uniprots

    def nodes_for_uniprot(self, accessions: list[str]) -> list[int]:
        return sorted({self.u.index(a) for a in accessions if a in self.u})

    def nodes_for_symbols(self, symbols: list[str]) -> list[int]:
        return sorted({self.s.index(a) for a in symbols if a in self.s})


def test_auc_perfect_chance_and_ties() -> None:
    positive = np.array([True, True, False, False])
    assert auc_from_scores(np.array([3.0, 4.0, 1.0, 2.0]), positive) == 1.0
    assert auc_from_scores(np.array([1.0, 2.0, 3.0, 4.0]), positive) == 0.0
    assert auc_from_scores(np.ones(4), positive) == 0.5


def test_percentile_rank_best_is_lowest() -> None:
    pct = percentile_ranks(np.array([5.0, 1.0, 3.0, 2.0]), np.array([True, False, False, False]))
    assert pct[0] == 0.25


def test_drug_bins_equal_count() -> None:
    bins = drug_bins(np.arange(100.0))
    assert set(np.bincount(bins)) == {10}


def test_bootstrap_ci_brackets_mean_and_is_deterministic() -> None:
    values = np.array([0.4, 0.5, 0.6, 0.7, 0.8])
    low, high = bootstrap_mean_ci(values, repeats=500)
    assert low <= values.mean() <= high
    assert (low, high) == bootstrap_mean_ci(values, repeats=500)


def _items(signal: bool) -> list[DiseaseData]:
    rng = np.random.default_rng(0)
    out = []
    for _ in range(6):
        score = rng.normal(size=200)
        positive = np.zeros(200, dtype=bool)
        positive[:8] = True
        if signal:
            score[:8] += 3
        out.append(DiseaseData(score, positive, drug_bins(rng.random(200))))
    return out


def test_permutation_detects_signal_and_not_noise() -> None:
    assert permutation_test(_items(True), repeats=200)["p_one_sided"] < 0.02
    assert permutation_test(_items(False), repeats=200)["p_one_sided"] > 0.05


def test_decision_rule_boundaries() -> None:
    assert decide({"mean_auc": 0.60, "ci95_low": 0.51})["supported_across_diseases"]
    assert not decide({"mean_auc": 0.59, "ci95_low": 0.55})["supported_across_diseases"]
    assert not decide({"mean_auc": 0.70, "ci95_low": 0.50})["supported_across_diseases"]
    assert not decide(None)["supported_across_diseases"]


def test_selection_rule_with_mocked_open_targets(monkeypatch) -> None:
    genes = [f"G{i}" for i in range(40)]
    drugs = [{"drug_id": f"CHEMBL{i}", "human_targets": [f"P{i}"]} for i in range(12)]
    universe = {"drugs": drugs}
    interactome = FakeInteractome(genes, [f"P{i}" for i in range(12)])
    approved = {
        "EFO_2": {f"CHEMBL{i}" for i in range(6)},
        "EFO_1": {f"CHEMBL{i}" for i in range(5)},
        "EFO_3": {"CHEMBL0"},                       # too few drugs
        "EFO_4": {f"CHEMBL{i}" for i in range(6)},  # too few associated targets
        "HP_1": {f"CHEMBL{i}" for i in range(6)},   # phenotype, excluded
    }
    monkeypatch.setattr(ot_client, "approved_indications", lambda ids: (approved, {k: k.lower() for k in approved}))
    monkeypatch.setattr(ot_client, "associated_symbols", lambda d, s: genes[:5] if d == "EFO_4" else genes[:25])
    chosen = select_diseases(universe, interactome, limit=2)  # type: ignore[arg-type]
    assert [c["efo_id"] for c in chosen] == ["EFO_1", "EFO_2"]
    assert len(select_diseases(universe, interactome, limit=30)) == 2  # type: ignore[arg-type]
