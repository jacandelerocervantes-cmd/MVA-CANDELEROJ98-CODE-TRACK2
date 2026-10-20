import numpy as np

from track2.phenotype_validation import MIN_VALID, evaluate, parse_hpoa, parse_mondo_omim, select_diseases

MONDO = """[Term]
id: MONDO:0000001
xref: OMIM:100100 {source="MONDO:equivalentTo"}
xref: Orphanet:5

[Term]
id: MONDO:0000002
xref: OMIM:100200
is_obsolete: true

[Term]
id: MONDO:0000003
name: no omim
"""

HPOA = ("#description: x\n"
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration\n"
        "OMIM:100100\tD\t\tHP:1\tR\tTAS\t\t\t\t\tP\tb\n"
        "OMIM:100100\tD\tNOT\tHP:2\tR\tTAS\t\t\t\t\tP\tb\n"
        "OMIM:100100\tD\t\tHP:3\tR\tTAS\t\t\t\t\tI\tb\n"
        "ORPHA:5\tD\t\tHP:4\tR\tTAS\t\t\t\t\tP\tb\n")


def test_mondo_and_hpoa_parsing(tmp_path):
    (tmp_path / "m.obo").write_text(MONDO, encoding="utf-8")
    (tmp_path / "h.hpoa").write_text(HPOA, encoding="utf-8")
    assert parse_mondo_omim(tmp_path / "m.obo") == {"MONDO_0000001": {"OMIM:100100"}}  # obsolete and OMIM-less terms are skipped
    assert parse_hpoa(tmp_path / "h.hpoa", {}) == {"OMIM:100100": {"HP:1"}}  # negated, non-phenotype and non-OMIM rows are skipped


def test_selection_requires_terms_and_known_drugs():
    mondo = {"MONDO_1": {"OMIM:1"}, "MONDO_2": {"OMIM:2"}, "MONDO_3": {"OMIM:3"}}
    hpoa = {"OMIM:1": {f"HP:{i}" for i in range(6)}, "OMIM:2": {"HP:1", "HP:2"}, "OMIM:3": {f"HP:{i}" for i in range(6)}}
    known = {"MONDO_1": {"a", "b", "c"}, "MONDO_2": {"a", "b", "c"}, "MONDO_3": {"a", "b"}, "MONDO_4": {"a", "b", "c"}}
    got = select_diseases(known, {"a", "b", "c"}, mondo, hpoa)
    assert [d["disease"] for d in got] == ["MONDO_1"]  # MONDO_2 has too few terms, MONDO_3 too few drugs, MONDO_4 no OMIM


def test_evaluate_marks_small_tests_as_not_valid():
    rng = np.random.default_rng(0)
    ids = [f"d{i}" for i in range(60)]
    degrees = {i: float(k) for k, i in enumerate(ids)}
    known = {f"X{n}": {"d0", "d1", "d2"} for n in range(3)}
    z_by = {(d, "v"): {i: (-5.0 if i in known[d] else float(rng.normal())) for i in ids} for d in known}
    res = evaluate(z_by, known, degrees, "v")
    assert res["n_evaluable"] == 3 < MIN_VALID and not res["test_is_valid"] and not res["validated"]
    assert all(v["auc"] > 0.95 for v in res["per_disease"].values())


def test_obsolete_flag_after_the_xrefs_still_excludes_the_term(tmp_path):
    (tmp_path / "m.obo").write_text("[Term]\nid: MONDO:9\nxref: OMIM:9\nis_obsolete: true\n\n[Term]\nid: MONDO:8\nxref: OMIM:8\n", encoding="utf-8")
    assert parse_mondo_omim(tmp_path / "m.obo") == {"MONDO_8": {"OMIM:8"}}
