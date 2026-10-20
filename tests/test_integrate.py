"""The enrichment step merges optional layers without changing the ranking and tolerates missing layers."""

from __future__ import annotations

import json
from pathlib import Path

from track2.integrate import integrate, write


def test_integrate_adds_layers_and_keeps_order(tmp_path: Path) -> None:
    base = [{"drug_id": "CHEMBL1", "name": "A", "combined_score": 2.0}, {"drug_id": "CHEMBL2", "name": "B", "combined_score": 1.0}]
    (tmp_path / "candidates.json").write_text(json.dumps(base), encoding="utf-8")
    (tmp_path / "faers_disproportionality.json").write_text(json.dumps({"drugs": {"CHEMBL1": {"available": True, "groups": {"renal": {"signal": True, "estimable": True, "prr": 3.0}}}}}), encoding="utf-8")
    rows = integrate(tmp_path)
    assert [r["drug_id"] for r in rows] == ["CHEMBL1", "CHEMBL2"]
    assert rows[0]["layers"]["faers_disproportionality"]["signal_groups"] == ["renal"]
    assert rows[1]["layers"]["evidence_retrieved"] == {"available": False} and rows[1]["layers"]["checkpoint_direction"] == ["not_run"]
    write(rows, tmp_path)
    assert (tmp_path / "candidates_enriched.csv").read_bytes().count(b"\r") == 0
