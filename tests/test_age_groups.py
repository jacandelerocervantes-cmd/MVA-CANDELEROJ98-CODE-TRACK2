from __future__ import annotations

from track2 import age_groups


def test_min_age_decides_group_coverage() -> None:
    profile = {"neonate": "NO_LABELLED_USE", "infant": "LABELLED_USE", "child": "LABELLED_USE", "adolescent": "LABELLED_USE", "adult": "LABELLED"}
    out = age_groups.group_status(profile, 1.0)
    assert out == {"0_to_under_3": "partly (from 1 y)", "3_to_under_12": "labelled", "12_and_older": "labelled"}
    out2 = age_groups.group_status({"child": "LABELLED_USE", "adolescent": "LABELLED_USE", "adult": "LABELLED"}, 2.0)
    assert out2["0_to_under_3"] == "partly (from 2 y)" and out2["3_to_under_12"] == "labelled"


def test_established_from_12_is_not_labelled_below() -> None:
    out = age_groups.group_status({"adolescent": "LABELLED_USE", "child": "NO_LABELLED_USE", "adult": "LABELLED"}, 12.0)
    assert out["0_to_under_3"] == "not_labelled" and out["3_to_under_12"] == "not_labelled" and out["12_and_older"] == "labelled"


def test_reviewed_not_established_means_adult_only_for_12_plus() -> None:
    out = age_groups.group_status({"adult": "LABELLED"}, None, "not_established")
    assert out == {"0_to_under_3": "not_labelled", "3_to_under_12": "not_labelled", "12_and_older": "adult_only"}


def test_unknown_labels_are_verify_not_assumed() -> None:
    out = age_groups.group_status({"neonate": "UNKNOWN", "infant": "UNKNOWN", "child": "VERIFY_MANUALLY", "adolescent": "UNKNOWN", "adult": "UNKNOWN"}, None)
    assert set(out.values()) == {"verify"}


def test_summarise_buckets() -> None:
    rows = [{"name": "A", "age_groups": {"0_to_under_3": "labelled", "3_to_under_12": "partly (from 4 y)", "12_and_older": "adult_only"}},
            {"name": "B", "age_groups": {"0_to_under_3": "verify", "3_to_under_12": "not_labelled", "12_and_older": "labelled"}}]
    s = age_groups.summarise(rows)
    assert s["0_to_under_3"]["labelled"] == ["A"] and s["3_to_under_12"]["partly"] == ["A"] and s["12_and_older"]["not_labelled_or_adult_only"] == ["A"]
