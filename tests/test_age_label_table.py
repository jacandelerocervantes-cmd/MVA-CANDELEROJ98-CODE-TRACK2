from __future__ import annotations

import json
from pathlib import Path

from track2 import age_label_table as alt

PED = {
    "MIXEDDRUG": "Safety and effectiveness have been established in pediatric patients 2 years and older for indication X. Safety and effectiveness have not been established for other uses. Growth plate abnormalities were observed in juvenile animals.",
    "NOTEST": "The safety and effectiveness of NOTEST in pediatric patients have not been established.",
    "REVIEWED": "Safety and effectiveness in pediatric patients below 1 year have not been established.",
    "NOSECTION": "",
}


def _fetch(name: str):
    if name == "GONE":
        return None, ["DailyMed: no label", "openFDA: HTTP 500"]
    return {"source": "DailyMed", "setid": f"set-{name}", "title": name, "label_effective": "20240101",
            "sections": {"pediatric_use": PED[name], "boxed_warning": "", "contraindications": "", "warnings": ""}}, []


UNIVERSE = [{"drug_id": f"C{i}", "name": n} for i, n in enumerate(["MIXEDDRUG", "NOTEST", "REVIEWED", "NOSECTION", "GONE"])]
REVIEW = {"C2": {"name": "REVIEWED", "verdict": "established", "min_age_years": 1.0, "source": "x", "scope": "y"}}


def test_mixed_is_verify_and_needs_manual_review() -> None:
    out = alt.build(UNIVERSE, {"MIXEDDRUG": ["t"]}, REVIEW, "2026-10-07", fetch=_fetch)
    r = out["farmacos"][0]
    assert r["clasificador_automatico"]["estado"] == "MIXED"
    assert set(r["estado_por_grupo"].values()) == {"verify"} and r["manual_review_required"] is True
    assert r["labelled_for_group"] == []
    assert "C0" in out["entradas_propuestas_para_config_NO_ESCRITAS"]


def test_literal_sentences_and_warning_quotes() -> None:
    out = alt.build(UNIVERSE, {"MIXEDDRUG": ["t"], "NOTEST": ["t"]}, REVIEW, "2026-10-07", fetch=_fetch)
    for r in out["farmacos"]:
        assert r["frase_literal_sustento"]["texto_literal"] in PED[r["name"]]
        for w in r["advertencias_pediatricas"]:
            assert w["texto_literal"] in PED[r["name"]]
    mixed = out["farmacos"][0]
    assert any("Growth plate" in w["texto_literal"] for w in mixed["advertencias_pediatricas"])
    notest = out["farmacos"][1]
    assert notest["estado_por_grupo"] == {"0_to_under_3": "not_labelled", "3_to_under_12": "not_labelled", "12_and_older": "adult_only"}


def test_reviewed_entry_sets_status_without_assuming_age() -> None:
    out = alt.build(UNIVERSE, {"REVIEWED": ["t"]}, REVIEW, "2026-10-07", fetch=_fetch)
    r = out["farmacos"][0]
    assert r["manual_review_required"] is False and r["revision_manual_en_config"] is True
    assert r["estado_por_grupo"]["0_to_under_3"] == "partly (from 1 y)"
    assert r["estado_por_grupo"]["3_to_under_12"] == "labelled" and r["labelled_for_group"] == ["12_and_older", "3_to_under_12"]


def test_missing_section_and_unreadable_label_are_verify() -> None:
    out = alt.build(UNIVERSE, {"NOSECTION": ["t"], "GONE": ["t"]}, REVIEW, "2026-10-07", fetch=_fetch)
    by = {r["name"]: r for r in out["farmacos"]}
    assert set(by["NOSECTION"]["estado_por_grupo"].values()) == {"verify"}
    assert by["GONE"]["etiqueta"] is None and by["GONE"]["manual_review_required"] is True
    assert out["sin_etiqueta_legible"] == ["GONE"]


def test_gather_warns_when_class_analysis_missing(tmp_path: Path) -> None:
    c = tmp_path / "c.json"
    c.write_text(json.dumps({"clases": {"k": {"farmacos": [{"name": "NOTEST"}]}}}), encoding="utf-8")
    f = tmp_path / "f.json"
    f.write_text(json.dumps({"readings": {"b2": {"one_line_only": ["REVIEWED", "UNKNOWNDRUG"]}}}), encoding="utf-8")
    origins, warnings, unmatched = alt.gather(c, tmp_path / "none.json", f, UNIVERSE)
    assert set(origins) == {"NOTEST", "REVIEWED"} and unmatched == ["UNKNOWNDRUG"]
    assert any("class_analysis" in w or "none.json" in w for w in warnings)
