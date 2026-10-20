from pathlib import Path

from track2 import lines_section as ls

OUT = Path(__file__).resolve().parents[1] / "out"


def _have() -> bool:
    return ls._have(OUT)


def test_empty_directory_gives_empty_text(tmp_path) -> None:
    for fn in (ls.main_sections, ls.cross_check, ls.chain_intro, ls.abstract_sentence, ls.brief_sentence, ls.criteria, ls.engine_scalability, ls.pairs_text, ls.signature_reversal_scalability):
        assert fn(tmp_path) == ""


def test_flag_and_age_helpers() -> None:
    assert ls._flag(["boxed warning in the label", "label read, no keyword match for the phenotype terms"]) == "boxed warning"
    assert ls._flag(["label read, no keyword match for the phenotype terms"]) == "none flagged"
    assert ls._flag(["needs manual verification: HP:0003202 muscle atrophy 'myopathy' in the warnings"]) == "verify: muscle atrophy"
    assert ls._age({"age_bands": {"0-4": "not_labelled", "5-11": "partly", "12-17": "labelled"}}) == "N/P/L"
    assert ls._ln(["C", "EP"]) == "C and E-P"


def test_one_table_per_line_with_candidates() -> None:
    if not _have():
        return
    d = ls._load(OUT, "candidates_by_line_mva_v2")
    text = ls.main_sections(OUT)
    assert text.startswith("### 3.13")
    with_rows = [k for k, v in d["lines"].items() if v["shown"]]
    assert text.count("| Drug | Native score |") == len(with_rows)
    for k in with_rows:
        assert ls.SHORT[k] in text
        for s in d["lines"][k]["shown"]:
            assert ls._name(s["drug"]) in text
    assert "Lines with no candidate" in text and "Convergence matrix" in text
    assert "2 case and 4 control arrays" in text and "2 arrays against 2" not in text
    assert " tier" not in text.lower()  # the build renames tier to group; the section avoids the word


def test_convergence_and_cross_check_follow_the_data() -> None:
    if not _have():
        return
    d = ls._load(OUT, "candidates_by_line_mva_v2")
    text = ls.main_sections(OUT)
    for drug in d["convergence"]["in_two_or_more_lines"]:
        assert ls._name(drug) in text
    cc = ls.cross_check(OUT)
    k = ls._load(OUT, "multi_evidence_mva_completed")
    o = ls._load(OUT, "multi_evidence_mva")
    assert cc.startswith("### 3.14")
    assert ls._name(k["proposed"][0]) in cc and ls._name(o["proposed"][0]) in cc
    assert f"{k['drugs'][k['proposed'][0]]['p_top5']:.2f}" in cc
    assert "not one more line" in cc and "artefact" in cc


def test_sentences_and_chain() -> None:
    if not _have():
        return
    for fn in (ls.abstract_sentence, ls.brief_sentence):
        s = fn(OUT)
        assert s.endswith(".") and "Section 3." in s
    assert "IGBP1" in ls.chain_intro(OUT) and "CDC20" in ls.chain_intro(OUT)
    c = ls.criteria(OUT)
    for crit in ("Scientific rigor", "Potential impact", "Innovation", "Scalability"):
        assert f"**{crit}.**" in c
