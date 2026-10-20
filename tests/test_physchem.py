from __future__ import annotations

from typing import Any

from track2 import physchem


def test_affinity_uses_only_exact_nanomolar_values() -> None:
    acts = [{"target_chembl_id": "T1", "standard_value": "10", "standard_relation": "=", "standard_units": "nM"},
            {"target_chembl_id": "T1", "standard_value": "1000", "standard_relation": "=", "standard_units": "nM"},
            {"target_chembl_id": "T1", "standard_value": "5", "standard_relation": ">", "standard_units": "nM"},
            {"target_chembl_id": "T2", "standard_value": "1", "standard_relation": "=", "standard_units": "uM"}]
    out = physchem.affinity_summary(acts)
    assert out == {"T1": {"n": 2, "median_pactivity": 7.0}}  # median of pActivity 8 and 6


def test_adme_does_not_mix_units() -> None:
    acts = [{"standard_type": "Solubility", "standard_units": "ug.mL-1", "standard_value": "2"},
            {"standard_type": "Solubility", "standard_units": "ug.mL-1", "standard_value": "4"},
            {"standard_type": "Solubility", "standard_units": "uM", "standard_value": "50"}]
    out = physchem.adme_summary(acts)
    assert out["Solubility (ug.mL-1)"] == {"n": 2, "median": 3.0} and out["Solubility (uM)"]["n"] == 1


def test_biologic_has_no_properties_and_no_invented_values() -> None:
    def fake(url: str) -> Any:
        if "/molecule/" in url:
            return {"molecule_type": "Antibody", "molecule_properties": None, "oral": False, "parenteral": True, "topical": False}
        return {"activities": []}
    out = physchem.profile({"drug_id": "CHEMBL1", "name": "X", "mechanisms": [{"target_id": "T1"}]}, fake)
    assert out["has_small_molecule_properties"] is False and out["properties"] == {} and out["measured_adme"] == {}
    assert out["routes"]["parenteral"] is True
