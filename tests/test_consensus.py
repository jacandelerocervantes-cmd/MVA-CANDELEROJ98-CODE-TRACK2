from __future__ import annotations

from track2 import consensus


def test_weights_schemes_and_pathogenic_gene() -> None:
    mechs = [{"id": "a", "weight": 1.0, "evidence_from_submission": {"clinvar": [{"classification": "Pathogenic/Likely pathogenic"}]}},
             {"id": "b", "weight": 0.1, "evidence_from_submission": {"clinvar": [{"classification": "Uncertain significance"}]}}]
    w = consensus.weights(mechs)
    assert w["W1_epcr"] == {"a": 1.0, "b": 0.1} and w["W2_equal"] == {"a": 1.0, "b": 1.0} and w["W3_pathogenic_gene_only"] == {"a": 1.0}


def test_combined_ignores_positive_z_and_sums_weighted() -> None:
    z = {"a": {"d1": -3.0, "d2": 2.0}, "b": {"d1": -1.0}}
    out = consensus.combined(z, {"a": 1.0, "b": 0.5})
    assert out == {"d1": 3.5, "d2": 0.0}
    assert consensus.top(out) == ["d1"]  # zero scores are not ranked


def test_w3_falls_back_to_w1_when_no_pathogenic_gene() -> None:
    mechs = [{"id": "a", "weight": 0.4, "evidence_from_submission": {"clinvar": []}}]
    assert consensus.weights(mechs)["W3_pathogenic_gene_only"] == {"a": 0.4}
