"""Offline checks of the scalability harness: tiny synthetic CSV parses, a trivial command is timed, protected paths are refused."""

import sys
from pathlib import Path

import pytest

from track2.empirical_scalability import guard, hash_tree, run_timed, stage_commands, write_case_csv
from track2.submission_input import read_submission


def test_synthetic_csv_parses(tmp_path):
    p = tmp_path / "case.csv"
    write_case_csv(p, [("FANCA", "16", 89000001, "A", "G", 0.9, "primary")])
    rows = read_submission(p)
    assert rows[0]["gene"] == "FANCA" and rows[0]["epcr"] == 0.9


def test_trivial_command_is_timed():
    r = run_timed([sys.executable, "-c", "print('ok')"])
    assert r["returncode"] == 0 and r["wall_seconds"] >= 0 and r["output_tail"] == ["ok"]


def test_guard_and_commands(tmp_path):
    with pytest.raises(SystemExit):
        guard(tmp_path / "out" / "x", tmp_path)
    guard(tmp_path / "out_case_x_run1", tmp_path)
    assert "--config-dir" in stage_commands(Path("a.csv"), Path("o"), Path("c"))[0][1]


def test_hash_tree(tmp_path):
    (tmp_path / "a.txt").write_text("x")
    assert list(hash_tree(tmp_path)) == ["a.txt"]
