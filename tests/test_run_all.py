from __future__ import annotations

import importlib
from pathlib import Path

from track2 import run_all, verify_outputs


def test_stage_numbers_are_unique_and_every_module_exists() -> None:
    table = run_all.stages()
    assert [s.number for s in table] == list(range(1, len(table) + 1))
    for stage in table:
        for cmd in stage.commands:
            module = cmd[2]  # [python, -m, module, ...]
            assert importlib.util.find_spec(module) is not None, module


def test_checksum_record_roundtrip_detects_a_change(tmp_path: Path, monkeypatch) -> None:
    (tmp_path / "out").mkdir()
    f = tmp_path / "out" / "consensus_ranking.json"
    f.write_text("{}", encoding="utf-8")
    assert verify_outputs.write(tmp_path) == 0 and verify_outputs.verify(tmp_path) == 0
    f.write_text("{ }", encoding="utf-8")
    assert verify_outputs.verify(tmp_path) == 1  # a deterministic file changed
