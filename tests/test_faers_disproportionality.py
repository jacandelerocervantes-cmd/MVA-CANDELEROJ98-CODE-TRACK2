"""Offline tests: the openFDA calls are replaced by an in-memory getter."""

from __future__ import annotations

import math
import urllib.parse
from typing import Any

from track2 import faers_disproportionality as fd


def test_clean_name_and_group_terms() -> None:
    assert fd.clean_name("DOCETAXEL ANHYDROUS") == "DOCETAXEL"
    assert fd.clean_name("ixazomib citrate") == "IXAZOMIB"
    terms = fd.build_group_terms([["HEPATIC FAILURE", "NAUSEA"], ["LIVER INJURY", "HEPATIC FAILURE"]])
    assert terms["hepatic"] == ["HEPATIC FAILURE", "LIVER INJURY"]
    assert fd.build_group_terms([["WERNICKE^S ENCEPHALOPATHY", "TOXIC ENCEPHALOPATHY"]])["neurological"] == ["TOXIC ENCEPHALOPATHY"]
    assert "NAUSEA" not in sum(terms.values(), [])


def test_prr_ror_known_table() -> None:
    r = fd.disproportionality(a=30, b=70, c=100, d=800)
    assert math.isclose(r["prr"], (30 / 100) / (100 / 900), abs_tol=1e-3)
    assert math.isclose(r["ror"], (30 * 800) / (70 * 100), abs_tol=1e-3)
    assert r["prr_ci95"][0] < r["prr"] < r["prr_ci95"][1]
    assert r["signal"] is True and r["ror_signal"] is True and not r["haldane_correction_applied"]


def test_signal_rule_needs_a_at_least_3() -> None:
    r = fd.disproportionality(a=2, b=2, c=10, d=10000)
    assert r["prr"] >= 2 and r["chi2_yates"] >= 4
    assert r["signal"] is False


def test_no_association_and_zero_cells() -> None:
    assert fd.disproportionality(a=10, b=90, c=100, d=900)["signal"] is False
    zero_a = fd.disproportionality(a=0, b=10, c=5, d=100)
    assert zero_a["estimable"] is False and "prr" not in zero_a
    corrected = fd.disproportionality(a=5, b=0, c=5, d=1000)
    assert corrected["haldane_correction_applied"] is True and corrected["table"]["b"] == 0


def test_group_clause_and_url_quoting() -> None:
    clause = fd.group_clause(["HEPATIC FAILURE", "LIVER INJURY"])
    assert clause == 'patient.reaction.reactionmeddrapt.exact:("HEPATIC FAILURE"+"LIVER INJURY")'
    url = fd._url(f"{fd.drug_search('ASPIRIN')} AND {clause}", "limit=1")  # noqa: SLF001
    assert "%20AND%20" in url and urllib.parse.unquote(url).count('"') == 6


def _fake_getter(drug_totals: dict[str, int]) -> Any:
    def get(url: str) -> Any:
        query = urllib.parse.unquote(url)
        if "count=" in query:
            if "search=" not in query:
                return {"results": [{"term": "HEPATIC FAILURE", "count": 9}, {"term": "NAUSEA", "count": 50}]}
            return {"results": [{"term": "LIVER INJURY", "count": 4}]}
        if "search=" not in query:
            return {"meta": {"results": {"total": 1000}}}
        for name, total in drug_totals.items():
            if f'generic_name.exact:"{name}"' in query:
                if "reactionmeddrapt" not in query:
                    return {"meta": {"results": {"total": total}}}
                return {"meta": {"results": {"total": 20 if "LIVER INJURY" in query else 5}}}
        if "reactionmeddrapt" in query and "generic_name" not in query:
            return {"meta": {"results": {"total": 30}}}
        raise LookupError(url)
    return get


def test_run_with_missing_drug_and_tables_add_up() -> None:
    candidates = [{"drug_id": "D1", "name": "ALPHA CITRATE"}, {"drug_id": "D2", "name": "UNKNOWN"}]
    result = fd.run(candidates, get=_fake_getter({"ALPHA": 100}), log=lambda m: None)
    assert result["drugs"]["D2"] == {"available": False, "name_queried": "UNKNOWN",
                                     "reason": "openFDA has no reports under this generic name"}
    hepatic = result["drugs"]["D1"]["groups"]["hepatic"]["table"]
    assert hepatic["a"] == 20 and hepatic["a"] + hepatic["b"] == 100
    assert sum(hepatic.values()) == result["faers_total_reports"]
    assert result["caveats"] and result["signal_rule"]["prr_min"] == 2.0
