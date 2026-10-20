from __future__ import annotations

from typing import Any

import pytest

from track2 import direction as d


def test_classify_pair() -> None:
    assert d.classify_pair("INHIBITOR", "positive") == d.WEAKEN
    assert d.classify_pair("INHIBITOR", "negative") == d.STRENGTHEN
    assert d.classify_pair("AGONIST", "positive") == d.STRENGTHEN
    assert d.classify_pair("ACTIVATOR", "negative") == d.WEAKEN
    assert d.classify_pair(None, "positive") == d.UNCLEAR
    assert d.classify_pair("INHIBITOR", "none") == d.NO_ANN
    assert d.classify_pair("INHIBITOR", "conflicting") == d.UNCLEAR
    assert d.classify_pair("INHIBITOR", "neutral") == d.UNCLEAR


def _row(acc: str, go: str, q: str = "involved_in") -> dict[str, Any]:
    return {"geneProductId": f"UniProtKB:{acc}", "goId": go, "qualifier": q, "goEvidence": "IMP",
            "evidenceCode": "ECO:0000315", "reference": "PMID:1"}


def test_build_review_offline(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake(url: str) -> dict[str, Any]:
        if "goId=GO:1" in url:
            return {"results": [_row("P1", "GO:11"), _row("P2", "GO:11", "NOT|involved_in")],
                    "pageInfo": {"total": 1}}
        if "goId=GO:2" in url:
            return {"results": [_row("P3", "GO:22")], "pageInfo": {"total": 1}}
        raise RuntimeError("down")

    monkeypatch.setattr(d, "get_json", fake)
    terms = [{"id": "GO:1", "role": "positive"}, {"id": "GO:2", "role": "negative"},
             {"id": "GO:3", "role": "positive"}]
    uni = {"targets": {f"T{n}": {"uniprot": [f"P{n}"]} for n in (1, 2, 3)},
           "drugs": [{"drug_id": "D1", "name": "x", "mechanisms": [
               {"action": "INHIBITOR", "target_id": "T1"}, {"action": "INHIBITOR", "target_id": "T3"},
               {"action": "INHIBITOR", "target_id": "T2"}]}]}
    r = d.build_review(uni, ["D1"], terms)
    cls = {p["uniprot"]: p["classification"] for p in r["drugs"]["D1"]["pairs"]}
    assert cls == {"P1": d.WEAKEN, "P2": d.NO_ANN, "P3": d.STRENGTHEN}
    assert r["drugs"]["D1"]["pairs"][0]["evidence"][0]["reference"] == "PMID:1"
    assert r["drugs"]["D1"]["caveat"] == d.CAVEAT
    assert len(r["api_problems"]) == 1 and "GO:3" in r["api_problems"][0]
