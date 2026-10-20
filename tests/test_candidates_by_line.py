"""Tests of the candidates-by-line module (synthetic inputs, offline, fast)."""

from __future__ import annotations

import json

from track2 import candidates_by_line as C
from track2 import multi_evidence as ME


def test_label_flags_whole_words_and_sections():
    secs = {"boxed_warning": "", "contraindications": "Do not use in active malignancy.", "pediatric_use": "growth suppression was seen",
            "warnings": "Cardiomyopathy may occur. Hypercalcemia and nephrocalcinosis were reported."}
    flags = C.label_flags(secs)
    terms = {f["term"] for f in flags}
    assert "HP:0000121 nephrocalcinosis" in terms  # warnings section
    assert not any("muscle atrophy" in t for t in terms)  # 'cardiomyopathy' is not 'myopathy'
    assert "HP:0002859 rhabdomyosarcoma" in terms  # 'active malignancy' counts only in contraindications
    assert not any("growth" in t for t in terms)  # the paediatric-use section is not searched
    assert C.label_flags({"warnings": "malignancy risk", "contraindications": "", "boxed_warning": ""}) == []
    broad = C.label_flags({"warnings": "Myopathy and rhabdomyolysis", "contraindications": "", "boxed_warning": ""})[0]
    assert broad["section"] == "warnings" and broad["exclusion"] is False  # broad keywords only flag
    spec = [f for f in flags if f["term"] == "HP:0000121 nephrocalcinosis"][0]
    assert spec["exclusion"] is True and "nephrocalcinosis" in spec["sentence"].lower()  # specific keyword in warnings: documented contraindication
    boxed = C.label_flags({"boxed_warning": "nephrocalcinosis", "contraindications": "", "warnings": ""})[0]
    assert boxed["exclusion"] is False  # the boxed-warning section alone does not exclude (13.1(c))


def make_ctx(flags_by_drug=None):
    flags_by_drug = flags_by_drug or {}
    return {"by_name": {n: {"molecule_type": "Small molecule", "black_box_warning": 0} for n in ("OK", "BAD", "CYTO", "HAZ", "NOPOOL", "UNREAD")} | {"ANTIBODY": {"molecule_type": "Antibody"}},
            "pool": {"OK", "BAD", "CYTO", "HAZ", "UNREAD", "ANTIBODY"}, "unconfirmed": {"UNCONF"}, "excluded": {"HAZ": "hazard"}, "caution_names": {"CYTO"},
            "labels": lambda d: {"read": d != "UNREAD", "flags": flags_by_drug.get(d, []), "source": "x", "failures": []}, "smaa": {"OK": {"p_top5": 0.3, "rank": 4}}, "safety_entries": set()}


def test_classify_drug_order_and_reasons():
    ctx = make_ctx({"BAD": [{"term": "HP:0000121 nephrocalcinosis", "keyword": "hypercalcemia", "section": "warnings", "sentence": "Do not give to patients with hypercalcemia.", "exclusion": True}],
                    "OK": [{"term": "HP:0003202 muscle atrophy", "keyword": "myopathy", "section": "warnings", "sentence": "Myopathy was seen.", "exclusion": False}]})
    assert C.classify_drug("OK", ctx) == ("candidate", "")  # a broad keyword hit is only a flag
    ctx["level6"] = {"LVL6"}
    ctx["by_name"]["LVL6"] = {"molecule_type": "Small molecule"}
    ctx["pool"].add("LVL6")
    assert C.classify_drug("LVL6", ctx)[0] == "hard_excluded" and "level 6" in C.classify_drug("LVL6", ctx)[1]
    ctx["safety_entries"] = {"BAD"}  # covered by the hand-read layer: no keyword screen, its label warnings are shown instead
    assert C.classify_drug("BAD", ctx) == ("candidate", "")
    ctx["safety_entries"] = set()
    assert C.classify_drug("HAZ", ctx)[0] == "hard_excluded"
    assert C.classify_drug("CYTO", ctx)[0] == "context_hazard"
    st, why = C.classify_drug("BAD", ctx)
    assert st == "hard_excluded" and "nephrocalcinosis" in why and "Quote:" in why
    assert C.classify_drug("UNCONF", ctx)[0] == "not_eligible"  # not in the universe
    ctx["by_name"]["UNCONF"] = {"molecule_type": "Small molecule"}
    assert C.classify_drug("UNCONF", ctx)[0] == "approval_unconfirmed"
    assert C.classify_drug("NOPOOL", ctx)[0] == "not_eligible" and "not in the completed pool" in C.classify_drug("NOPOOL", ctx)[1]
    assert "small molecule" in C.classify_drug("ANTIBODY", {**ctx, "pool": set()})[1]
    assert C.classify_drug("UNREAD", ctx)[0] == "candidate"  # an unreadable label is flagged, not read as a warning


def test_evaluate_line_caps_filters_and_keeps_others():
    ctx = make_ctx()
    native = [{"drug": d, "rank_in_line": i + 1, "native": f"n{i}"} for i, d in enumerate(["HAZ", "OK", "CYTO", "UNREAD", "BAD"])]
    rows = {"OK": {"grades": {c: None for c in ME.CRITERIA}, "lit": {"best_help_level": None, "n_help_papers": 0, "n_harm_papers": 0}, "consequences": {"4": 0.5}, "track1_genes": ["IGBP1"]}}
    ctx["age_fn"] = lambda d: {"bands": {"0-4": "labelled", "5-11": "partly", "12-17": "not_labelled"}, "retrievable": True, "quote": "q", "quote_section": "pediatric_use", "label_status": "ESTABLISHED", "basis": "manual review"}
    res = C.evaluate_line("A", native, ctx, rows, {"OK": {"0-4": "labelled", "5-11": "partly", "12-17": "not_labelled"}}, limit_candidates=1, qualify_cap=5)
    assert [c["drug"] for c in res["shown"]] == ["OK"] and res["also_qualifying"] == ["UNREAD", "BAD"] and res["n_qualifying"] == 3
    ok = res["shown"][0]
    assert ok["age_union"] == "supported in 0-4, 5-11" and ok["smaa_p_top5_cross_check"] == 0.3 and ok["consequences_addressed"] == {"4": 0.5}
    assert ok["age_reading"]["quote_section"] == "pediatric_use"
    assert "label read" in " ".join(ok["safety_flags"]) and ok["refuting_experiment"] == C.EXPERIMENT["A"] and "line limit" in ok["evidence_tier"]
    assert {o["drug"] for o in res["found_but_not_proposed"]} >= {"HAZ", "CYTO"}


def test_native_sets_from_synthetic_files(tmp_path):
    out = tmp_path
    (out / "line_c_drug_level.json").write_text(json.dumps({"drugs": {"X": {"rho": 0.2, "q_bh": 0.02, "n_lines": 300}, "Y": {"rho": -0.3, "q_bh": 0.01, "n_lines": 300},
                                                                       "Z": {"rho": 0.5, "q_bh": 0.4, "n_lines": 300}, "W": {"rho": 0.1, "q_bh": 0.001, "n_lines": 300}}}))
    rows = C.native_C(out)
    assert [r["drug"] for r in rows] == ["W", "X"]  # correct direction and q < 0.10 only, ranked by q
    row = {"name": "DRUGA", "bub1b": None, "P": {"usable": True, "hazard": False, "p_rev": 0.002, "q_rev": 0.05, "R": 9}, "track1_combined": {"usable": False, "hazard": False, "p_rev": 0.1, "q_rev": 0.6, "R": 1}}
    row2 = {"name": "DRUGB", "bub1b": None, "P": {"usable": True, "hazard": True, "p_rev": 0.001, "q_rev": 0.05, "R": 9}, "track1_combined": None}
    (out / "line_e_signature_reversal.json").write_text(json.dumps({"table_all_drugs_n150": [row, row2]}))
    assert [r["drug"] for r in C.native_E(out, "P")] == ["DRUGA"] and C.native_E(out, "track1_combined") == [] and C.native_E(out, "bub1b") == []
    (out / "phenotype_proximity_e5.json").write_text(json.dumps({"rule": {"configurations": 9, "robust_min": 7}, "drugs": [{"name": "P", "robust": False, "robustness_count": 5}, {"name": "Q", "robust": False, "robustness_count": 6}]}))
    robust, near = C.native_B1(out)
    assert robust == [] and [n["drug"] for n in near] == ["Q", "P"]
    assert C.native_B2(tmp_path) == [] and C.native_engine(None) == []
    eng = C.native_engine({"n_fda_small_molecule": 3, "drugs": {"A": {"rank_fda_small_molecule": 2, "score": -1.0}, "B": {"rank_fda_small_molecule": 1, "score": -0.5}, "C": {"rank_fda_small_molecule": None, "score": 0}}})
    assert [r["drug"] for r in eng] == ["B", "A"]


def test_convergence_overlap_and_markdown_wording():
    sets = {"A": {"X", "Y"}, "B": {"Y", "Z"}, "C": set(), "D": {"W"}}
    cv = C.convergence(sets, 224)
    assert cv["n_candidates_distinct"] == 4 and cv["n_unique_to_one_line"] == 3 and cv["in_two_or_more_lines"] == ["Y"]
    pair = next(p for p in cv["pairs"] if p["lines"] == ["A", "B"])
    assert pair["overlap"] == 1 and 0 < pair["p_hypergeometric"] < 1
    assert all(p["overlap"] == 0 and p["p_hypergeometric"] == 1.0 for p in cv["pairs"] if "D" in p["lines"])
    res = {"label": "hypotheses for follow-up", "protocol": "p", "filters": {"pool_size": 224, "approval_unconfirmed_28": ["A1"]},
           "lines": {ln: {"name": ln, "n_native_set": 0, "n_qualifying": 0, "shown": [], "also_qualifying": [], "found_but_not_proposed": [], "limit": "l"} for ln in C.LINES},
           "b1_near_threshold_information_only": [], "convergence": cv, "candidates": {}, "cross_check": {"proposed_by_smaa": ["X"], "trametinib_survives": False, "stability": {}},
           "track1_gene_block_statement": "hypotheses to test, not findings", "data_gaps": ["g"]}
    md = C.to_markdown(res)
    assert "No drug meets the line's native criterion" in md and "never required" in md and "Convergence matrix" in md and "not findings" in md
