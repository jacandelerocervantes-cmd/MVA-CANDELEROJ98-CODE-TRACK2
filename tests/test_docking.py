from __future__ import annotations

from track2 import docking


def test_auc_counts_more_negative_as_better_and_ties_half() -> None:
    assert docking.auc([-9.0, -8.5], [-6.0, -7.0]) == 1.0
    assert docking.auc([-6.0], [-9.0]) == 0.0
    assert docking.auc([-7.0], [-7.0]) == 0.5


def test_threshold_requires_three_actives_and_precision() -> None:
    actives = [-10.0, -9.8, -9.5, -8.0]
    decoys = [-9.0, -8.5, -7.0, -6.5, -6.0]
    t = docking.threshold(actives, decoys)
    assert t == -6.5  # the least negative cut-off that keeps precision >= 0.5 (4 actives and 4 decoys at or below it)
    assert docking.threshold([-9.0, -8.0], decoys) is None  # fewer than three actives can never qualify


def test_heavy_atom_and_rotatable_bond_limits() -> None:
    assert docking.heavy_ok("CC(=O)Oc1ccccc1C(=O)O")
    assert not docking.heavy_ok("C" * 60)
    assert not docking.heavy_ok("")


def test_validation_sets_are_reproducible_and_exclude_actives_from_decoys() -> None:
    universe = {"drugs": [{"drug_id": f"D{i}", "smiles": "CCO", "molecule_type": "Small molecule", "human_targets": ["P1"] if i < 3 else []} for i in range(30)]}
    a1, d1 = docking.validation_sets(universe, {}, "P1")
    a2, d2 = docking.validation_sets(universe, {}, "P1")
    assert set(a1) == {"D0", "D1", "D2"} and d1 == d2 and not set(a1) & set(d1)
