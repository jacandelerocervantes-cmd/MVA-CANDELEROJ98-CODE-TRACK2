"""Data manifest tests on a synthetic data directory (offline, no patient data)."""

from __future__ import annotations

import json
from pathlib import Path

from track2.snapshot import build_manifest, group_digest


def _make(tmp_path: Path) -> tuple[Path, Path]:
    data, out = tmp_path / "data", tmp_path / "out"
    (data / "cache").mkdir(parents=True)
    (data / "private").mkdir()
    out.mkdir()
    (data / "private" / "genome_covariates.json").write_text("SECRET-PATIENT", encoding="utf-8")
    (data / "universe.json").write_text("[]", encoding="utf-8")
    (data / "cache" / ("a" * 24 + ".json")).write_text("{}", encoding="utf-8")
    (data / "cache" / ("b" * 24 + ".json")).write_text("[1]", encoding="utf-8")
    (data / "cache" / "spl_x.xml").write_text("<x/>", encoding="utf-8")
    (data / "cache" / "9606.protein.aliases.v12.0.txt.gz").write_bytes(b"alias")
    (out / "manifest_m1.json").write_text(json.dumps({"chembl_version": "ChEMBL_37", "string_version": "v12.0", "reactome_version": "97"}), encoding="utf-8")
    return data, out


def test_manifest_groups_and_versions(tmp_path: Path) -> None:
    data, out = _make(tmp_path)
    m = build_manifest(data, out)
    assert m["api_json_responses"]["count"] == 2
    assert m["dailymed_spl_xml"]["count"] == 1
    assert list(m["large_downloads"]) == ["9606.protein.aliases.v12.0.txt.gz"]
    assert m["string_files"]["9606.protein.aliases.v12.0.txt.gz"]["bytes"] == 5
    assert m["string_files"]["9606.protein.physical.links.v12.0.txt.gz"] is None
    assert m["universe_json"]["bytes"] == 2
    assert m["recorded_versions"]["chembl"] == ["ChEMBL_37"]
    assert "SECRET-PATIENT" not in json.dumps(m) and "genome_covariates" not in json.dumps(m)


def test_digest_is_order_independent_and_deterministic(tmp_path: Path) -> None:
    assert group_digest([("b", "2"), ("a", "1")]) == group_digest([("a", "1"), ("b", "2")])
    data, out = _make(tmp_path)
    assert build_manifest(data, out) == build_manifest(data, out)


def test_missing_data_is_reported_not_crashing(tmp_path: Path) -> None:
    m = build_manifest(tmp_path / "nodata", tmp_path / "noout")
    assert m["api_json_responses"]["count"] == 0 and m["universe_json"] is None
