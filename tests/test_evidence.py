from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from track2 import evidence

CONFIG_PATH = Path(__file__).resolve().parent.parent / "config" / "evidence_queries.json"
CONFIG = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def fake_get_json(url: str) -> Any:
    if "esearch" in url:
        return {"esearchresult": {"count": "2", "idlist": ["111", "222"]}}
    if "esummary" in url:
        return {"result": {"111": {"title": "Real paper", "pubdate": "2020 Jan"},
                           "222": {"error": "cannot get document summary"}}}
    return {"totalCount": 7, "studies": [{"protocolSection": {
        "identificationModule": {"nctId": "NCT1", "briefTitle": "T"}, "statusModule": {"overallStatus": "COMPLETED"},
        "designModule": {"phases": ["PHASE2"]}, "eligibilityModule": {"minimumAge": "1 Year"}}}]}


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(evidence, "PAUSE_SECONDS", 0)


def test_clean_name() -> None:
    salts = CONFIG["name_cleaning"]["salt_words"]
    assert evidence.clean_name("DOCETAXEL ANHYDROUS", salts) == "docetaxel"
    assert evidence.clean_name("FIBRINOLYSIN, HUMAN", salts) == "fibrinolysin"


def test_drug_evidence_verifies_pmids_and_trials(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(evidence, "get_json", fake_get_json)
    result = evidence.drug_evidence({"name": "IXAZOMIB CITRATE"}, CONFIG, "2026-01-01")
    group = result["pubmed"]["aneuploidy_cin"]
    assert group["count"] == 2 and group["verified_count"] == 1
    assert group["records"][0] == {"pmid": "111", "title": "Real paper", "verified": True, "year": 2020}
    assert group["records"][1]["verified"] is False
    assert group["query"].startswith('"ixazomib" AND (')
    assert result["trials"]["pediatric_population"]["total_count"] == 7
    assert result["trials"]["rhabdomyosarcoma"]["studies"][0]["nct_id"] == "NCT1"
    assert "not curated" in result["caveat"]


def test_api_failure_is_recorded_not_hidden(monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(url: str) -> Any:
        raise RuntimeError("down")

    monkeypatch.setattr(evidence, "get_json", boom)
    result = evidence.drug_evidence({"name": "ASPIRIN"}, CONFIG, "2026-01-01")
    assert result["pubmed"]["pediatric"]["available"] is False
    assert "down" in result["trials"]["wilms_tumor"]["reason"]


def test_run_writes_sorted_lf(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(evidence, "get_json", fake_get_json)
    cand = tmp_path / "c.json"
    cand.write_text(json.dumps([{"drug_id": "B", "name": "X"}, {"drug_id": "A", "name": "Y"}]), encoding="utf-8")
    path = evidence.run(cand, tmp_path, CONFIG_PATH)
    raw = path.read_bytes()
    assert b"\r" not in raw and list(json.loads(raw)) == ["A", "B"]
