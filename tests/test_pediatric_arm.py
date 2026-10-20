from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from track2 import pediatric_arm


def test_ranks_are_dense_and_tie_broken_by_id() -> None:
    assert pediatric_arm.ranks({"b": 2.0, "a": 2.0, "c": 5.0}) == {"c": 1, "a": 2, "b": 3}


def _write(tmp: Path) -> None:
    mech = {"id": "sub_G", "status": "mechanism_candidate", "weight": 1.0, "evidence_from_submission": {"clinvar": []}}
    (tmp / "mech").mkdir()
    (tmp / "mech" / "sub_G.json").write_text(json.dumps(mech), encoding="utf-8")
    z = {"D1": -5.0, "D2": -4.0, "D3": -0.5}
    (tmp / "threshold_scores_sub_G.json").write_text(json.dumps({t: z for t in ("700", "400", "900")}), encoding="utf-8")
    (tmp / "scores_sub_G.json").write_text(json.dumps({"rows": [{"drug_id": d, "risk_flags": ["X"]} for d in z]}), encoding="utf-8")


def test_only_paediatric_labelled_drugs_are_kept_and_stability_is_reported(tmp_path: Path, monkeypatch: Any) -> None:
    _write(tmp_path)
    status = {"A": "ESTABLISHED", "B": "NOT_ESTABLISHED", "C": "MIXED"}
    names = {"D1": "A", "D2": "B", "D3": "C"}

    def fake_label(name: str) -> dict[str, Any]:
        return {"found": True, "pediatric": {"status": status[name], "min_age_years": 2.0, "excerpt": "e"},
                "age_profile": {"child": "LABELLED_USE"}, "excipient_alerts": []}
    monkeypatch.setattr(pediatric_arm, "label_info", fake_label)
    universe = {"drugs": [{"drug_id": d, "name": n, "mechanisms": [{"mechanism": "Tubulin inhibitor"}]} for d, n in names.items()]}
    config = {"classes": [], "tumour_context_reports": {"applies_to_flag": "NONE", "note": "", "pmids": []}}
    result = pediatric_arm.run(tmp_path, tmp_path / "mech", universe, config)
    assert [r["name"] for r in result["drugs"]] == ["A", "C"]  # B (not established) is excluded
    assert all(r["stability"] == "stable" and r["n_configs_rank_le_100"] == 9 for r in result["drugs"])
    assert result["drugs"][0]["rank_700_W1"] == 1
