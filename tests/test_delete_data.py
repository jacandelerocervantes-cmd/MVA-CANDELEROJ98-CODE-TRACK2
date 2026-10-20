"""The deletion tool is dry-run by default and removes only what it lists."""

from __future__ import annotations

from pathlib import Path

from track2.delete_data import delete, inventory


def test_dry_run_deletes_nothing_and_execute_deletes_and_attests(tmp_path: Path) -> None:
    data = tmp_path / "data"
    (data / "private").mkdir(parents=True)
    (data / "private" / "x.json").write_text("secret", encoding="utf-8")
    keep = tmp_path / "out"
    keep.mkdir()
    (keep / "candidates.csv").write_text("kept", encoding="utf-8")
    attestation = tmp_path / "out" / "att.json"

    report = delete([data], execute=False, attestation=attestation)
    assert data.exists() and not attestation.exists() and report["targets"][0]["files"] == 1

    report = delete([data], execute=True, attestation=attestation)
    assert not data.exists() and (keep / "candidates.csv").exists()
    assert attestation.exists() and "secret" not in attestation.read_text(encoding="utf-8")
    assert inventory([data]) == []
