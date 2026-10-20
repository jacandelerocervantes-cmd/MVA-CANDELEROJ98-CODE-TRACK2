"""DATA_ACCESS.md covers every resource of out/license_check.json and holds no local path or coordinate."""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOC = (ROOT / "DATA_ACCESS.md").read_text(encoding="utf-8")
ROWS = [line for line in DOC.splitlines() if line.startswith("|")]


def test_every_license_check_resource_has_a_row():
    records = json.loads((ROOT / "out" / "license_check.json").read_text(encoding="utf-8"))["records"]
    missing = []
    for r in records:
        key = r["resource"].split(" (")[0].strip()
        if not any(key in row for row in ROWS):
            missing.append(r["resource"])
    assert not missing, missing


def test_no_local_path_or_coordinate():
    patterns = [
        r"[A-Za-z]:\\", r"/Users/", r"/home/", r"[\\]Users[\\]",
        r"\bchr(?:[0-9]{1,2}|X|Y|M)\b", r"\bchr\w*:\d+", r"\bNC_0000\d+", r"\b\d+:\d{4,}-\d{4,}\b",
        r"\bg\.\d+", r"\bc\.\d+[ACGT>_]",
    ]
    hits = [p for p in patterns if re.search(p, DOC)]
    assert not hits, hits


def test_groups_present():
    for g in ("### A.", "### B.", "### C.", "### D.", "## Reproducing without the files", "## Terms that matter for reuse"):
        assert g in DOC


def test_every_not_recorded_cell_gives_a_reason():
    bad = [c.strip()[:60] for row in ROWS for c in row.split("|") if "not recorded" in c and "no record in the repository" not in c]
    assert not bad, bad


def test_recorded_hashes_are_well_formed():
    for h in re.findall(r"SHA-256 `([0-9a-f]+)`", DOC):
        assert len(h) == 64
