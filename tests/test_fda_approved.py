from __future__ import annotations

from typing import Any

import pytest

from track2 import fda_approved


def _result(statuses: list[str], app: str = "NDA1") -> dict[str, Any]:
    return {"results": [{"application_number": app, "products": [{"marketing_status": s} for s in statuses]}]}


def test_candidate_names_order_and_salt_stripping() -> None:
    assert fda_approved.candidate_names("Ixazomib Citrate") == ["IXAZOMIB CITRATE", "IXAZOMIB"]
    assert fda_approved.candidate_names("DOCETAXEL ANHYDROUS") == ["DOCETAXEL ANHYDROUS", "DOCETAXEL"]


def test_eligible_needs_a_prescription_product() -> None:
    assert fda_approved.lookup("X", lambda url: _result(["Discontinued", "Prescription"]))["eligible"] is True
    out = fda_approved.lookup("X", lambda url: _result(["Discontinued"]))
    assert out["eligible"] is False and out["reason"] == "no product with marketing status Prescription"


def test_no_match_falls_back_to_shorter_names_then_not_eligible() -> None:
    calls: list[str] = []

    def fake(url: str) -> Any:
        calls.append(url)
        if "SODIUM" in url or "FOO%20BAR" in url:
            raise LookupError(url)
        return _result(["Prescription"])
    assert fda_approved.lookup("FOO SODIUM", fake)["matched_name"] == "FOO"  # salt-stripped name matches after the full name fails

    def never(url: str) -> Any:
        raise LookupError(url)
    assert fda_approved.lookup("NOPE", never)["eligible"] is False


def test_service_error_is_not_read_as_not_eligible() -> None:
    def broken(url: str) -> Any:
        raise RuntimeError("service down")
    with pytest.raises(RuntimeError):
        fda_approved.lookup("X", broken)
