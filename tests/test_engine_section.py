from pathlib import Path

from track2 import engine_section as es

ENG = Path(__file__).resolve().parents[1] / "out" / "engine"


def test_empty_directory_gives_empty_text(tmp_path) -> None:
    assert es.main_section(tmp_path) == ""
    assert es.supplement(tmp_path) == ""
    assert es.abstract_sentence(tmp_path) == ""


def test_ci_formatting_has_no_negative_zero() -> None:
    assert es._ci({"estimate": -0.001, "ci95": [-0.09, 0.004]}) == "0.00 [-0.09, 0.00]"
    assert es._ci({"estimate": 0.206, "ci95": [0.049, 0.343]}) == "+0.21 [+0.05, +0.34]"


def test_section_is_generated_from_the_engine_files() -> None:
    if not (ENG / "validation_v4.json").exists():
        return
    text = es.main_section(ENG)
    assert text.startswith("### 3.12")
    assert "Pre-registered primary (17 diseases)" in text
    assert text.count("| 100 % |") == 1 and text.count("| 1 % |") == 1
    assert "not testable with those labels" in text and "not informative" in text
    assert "S14" in es.supplement(ENG)
