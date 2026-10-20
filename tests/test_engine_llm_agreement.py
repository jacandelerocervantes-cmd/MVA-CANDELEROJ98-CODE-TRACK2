import json

import pytest

from track2 import engine_llm_agreement as m


def test_reliability_flags() -> None:
    assert m.reliability(["benefit", "benefit", "benefit"]) == "robust"
    assert m.reliability(["benefit", "benefit", "harm"]) == "majority"
    assert m.reliability(["benefit", "harm", "unclear"]) == "unresolved"


def test_load_annotations_checks_ids_and_labels(tmp_path) -> None:
    p = tmp_path / "annotations_opus.json"
    p.write_text(json.dumps({"annotator": "opus", "labels": [{"id": 0, "label": "benefit"}, {"id": 1, "label": "harm"}]}), encoding="utf-8")
    assert m.load_annotations(tmp_path, 2) == {"opus": ["benefit", "harm"]}
    with pytest.raises(ValueError):
        m.load_annotations(tmp_path, 3)
    p.write_text(json.dumps({"annotator": "opus", "labels": [{"id": 0, "label": "maybe"}, {"id": 1, "label": "harm"}]}), encoding="utf-8")
    with pytest.raises(ValueError):
        m.load_annotations(tmp_path, 2)
