from pathlib import Path

from track2 import followup_section as fs

OUT = Path(__file__).resolve().parents[1] / "out"


def test_empty_directory_gives_empty_text(tmp_path) -> None:
    assert fs.main_section(tmp_path) == ""
    assert fs.supplement(tmp_path) == ""
    assert fs.abstract_sentence(tmp_path) == ""


def test_section_is_read_from_the_experiment_plan() -> None:
    if not (OUT / "experiment_plan.json").exists():
        return
    plan = fs._load(OUT, "experiment_plan")
    text = fs.main_section(OUT)
    assert text.startswith("### 3.15")
    for cf in plan["computational_followups"]:
        assert f"| {cf['id']} |" in text
    for ex in plan["experiments"]:
        assert f"| {ex['id']} |" in text
    assert "Step 0" in text and "Line C adjusted" in text


def test_supplement_has_one_row_per_triage_candidate_and_corrects_array_counts() -> None:
    if not (OUT / "experiment_plan.json").exists():
        return
    plan = fs._load(OUT, "experiment_plan")
    s = fs.supplement(OUT)
    assert s.startswith("## S15. Follow-up experiments")
    assert "| Candidate | Lines |" in s
    assert sum(1 for ln in s.splitlines() if ln.startswith("| ") and "Candidate" not in ln) == len(plan["triage"])
    assert "2 vs 2 arrays" not in s
    a = fs.abstract_sentence(OUT)
    assert a.endswith(".") and "Step 0" in a
