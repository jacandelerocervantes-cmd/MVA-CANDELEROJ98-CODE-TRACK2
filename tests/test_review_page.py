"""Review page: offline unit tests on a synthetic record."""

from __future__ import annotations

import json
from pathlib import Path

from track2.review_page import index_by_drug, load_optional, render


def _record(name: str = "DRUG<b>", score: float = 1.5) -> dict:
    return {
        "name": name, "drug_id": "CHEMBL1", "combined_score": score, "direction_review": "REQUIRES_REVIEW", "risk_flags": ["INHIBITS_MODULE_COMPONENT"],
        "mechanism_of_action": ["Some inhibitor"], "actions": ["INHIBITOR"], "targets": ["P1"],
        "mechanism_scores": [{"mechanism_id": "m", "rank": 3, "z": 2.0, "d": 1.0, "p_empirical": 0.01, "passes_primary_threshold": True}],
        "gates": [{"requirement": "RENAL_SAFETY", "status": "NOT_ASSESSABLE", "needs": "renal_function"}],
        "age_profile": {"adult": "LABELLED", "child": "NO_LABELLED_USE", "pediatric_label_status": "NOT_ESTABLISHED"},
        "excipients": {"alerts": [{"excipient": "benzyl alcohol"}], "inactive_ingredients": ["TALC"]},
        "toxicity": {"available": True, "caveat": "no incidence", "total_reports": 10, "top_terms": [{"term": "NAUSEA", "count": 4}]},
        "evidence_for": [{"type": "DailyMed", "id": "setid=x", "citation_verified": False}],
        "evidence_against": [{"source": "rule", "text": "review direction"}],
    }


def test_render_core_content_and_escaping() -> None:
    page = render([_record(), _record("OTHER", 3.0)], {}, ["direction_review.json"])
    assert "Hypotheses for follow-up, not medical advice" in page
    assert "DRUG&lt;b&gt;" in page and "DRUG<b>" not in page
    assert "http://" not in page and "https://" not in page  # self-contained: no external resources
    assert page.index("OTHER") < page.index("DRUG&lt;")  # ordered by score
    assert "direction_review.json" in page and "aria-expanded" in page and "benzyl alcohol" in page


def test_optional_files_are_merged_and_absence_is_handled(tmp_path: Path) -> None:
    (tmp_path / "reversal_scores.json").write_text(json.dumps({"rows": [{"drug_id": "CHEMBL1", "score": -0.4}]}), encoding="utf-8")
    optional, missing = load_optional(tmp_path)
    assert "reversal_scores.json" in optional and len(missing) == 4
    page = render([_record()], optional, missing)
    assert "Transcriptomic reversal" in page and "-0.4" in page


def test_index_by_drug_shapes() -> None:
    assert index_by_drug([{"name": "A"}]) == {"A": {"name": "A"}}
    assert index_by_drug({"X": {"v": 1}}) == {"X": {"v": 1}}
    assert index_by_drug(5) == {}
