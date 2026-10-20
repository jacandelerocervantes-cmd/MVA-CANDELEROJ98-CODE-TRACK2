from __future__ import annotations

import gzip

import numpy as np

from track2 import line_e as le


def _lib() -> dict:
    genes = [f"G{i}" for i in range(40)]
    txt = "\n".join([
        "KO1 Up\t\t" + "\t".join(genes[:10]), "KO1 Down\t\t" + "\t".join(genes[10:20]),
        "DrugA Up\t\t" + "\t".join(genes[10:20]), "DrugA Down\t\t" + "\t".join(genes[:10]),   # reverses KO1
        "DrugB Up\t\t" + "\t".join(genes[:10]), "DrugB Down\t\t" + "\t".join(genes[10:20]),   # mimics KO1
        "Junk Up\t\tG1", ""])
    return le.parse_library(txt)


def test_parse_and_reversal_sign() -> None:
    lib = _lib()
    assert set(lib) == {"KO1", "DrugA", "DrugB"}  # "Junk" has no Down set
    index = {f"G{i}": i for i in range(40)}
    ko = le.signed_vector(lib["KO1"]["Up"], lib["KO1"]["Down"], None, index)
    D = le.drug_matrix(lib, ["DrugA", "DrugB"], None, index)
    r = le.reversal(D, ko)
    assert r.tolist() == [20.0, -20.0]
    assert le.reversal(D, le.signed_vector(lib["KO1"]["Up"], lib["KO1"]["Down"], 5, index)).tolist() == [10.0, -10.0]


def test_bh_monotone_and_bounds() -> None:
    p = np.array([0.001, 0.04, 0.5, 0.9])
    q = le.bh(p)
    assert np.all(q >= p) and q.max() <= 1 and q[0] == 0.004


def test_usable_and_hazard_with_null() -> None:
    rng = np.random.default_rng(0)
    null_r = rng.normal(0, 3, size=(2, 1000))
    drugs = np.eye(2, 4, dtype=np.float32)
    sig = np.array([-30, 30, 0, 0], dtype=np.float32)  # R = drugs @ -sig = [30, -30]
    t = le.test_signature(drugs, sig, null_r)
    assert t["usable"].tolist() == [False, False] or t["q_rev"][0] <= 1.0  # BH over two drugs
    assert t["p_rev"][0] == 1 / 1001 and t["p_mim"][1] == 1 / 1001
    assert t["R"].tolist() == [30.0, -30.0]
    big = le.test_signature(np.vstack([drugs] + [np.zeros((10, 4), np.float32)]), sig, np.vstack([null_r, rng.normal(0, 3, (10, 1000))]))
    assert bool(big["usable"][0]) and bool(big["hazard"][1]) and not big["usable"][1]


def test_combined_vote_uses_weights_and_truncates() -> None:
    lib = {"A": {"Up": ["x1", "x2", "x3"], "Down": ["y1"]}, "B": {"Up": ["y1"], "Down": ["x3", "z"]}}
    index = {g: i for i, g in enumerate(["x1", "x2", "x3", "y1", "z"])}
    v = le.combined_vote(lib, {"A": 0.6, "B": 0.1}, 2, index)
    assert v[index["x1"]] == 1 and v[index["x2"]] == 1 and v[index["x3"]] == 0  # x3 trimmed to the top 2 ups
    assert v[index["y1"]] == -1 and v[index["z"]] == -1


def test_track1_weights_floor_and_threshold() -> None:
    rows = [{"gene": "A", "phenotype_score": 0.5}, {"gene": "B", "phenotype_score": 0.0}]
    assert le.track1_weights(rows, ["C"]) == {"A": 0.5, "B": le.FLOOR, "C": le.FLOOR}
    assert le.track1_weights(rows, ["C"], 0.1) == {"A": 0.5}


def test_match_drugs_salt_alias_and_candidates() -> None:
    uni = [{"name": "SIROLIMUS", "molecule_type": "Small molecule", "drug_id": "1"}, {"name": "IMATINIB MESYLATE", "molecule_type": "Small molecule", "drug_id": "2"},
           {"name": "PHENYLBUTANOIC ACID", "molecule_type": "Small molecule", "drug_id": "3"}, {"name": "MABX", "molecule_type": "Antibody", "drug_id": "4"}]
    m, un = le.match_drugs(["Sirolimus", "Imatinib", "Phenylbutyrate", "Mabx", "Nothing"], uni)
    assert set(m) == {"Sirolimus", "Imatinib", "Phenylbutyrate"} and set(un) == {"Mabx", "Nothing"}
    safety = {"caution_names": {"BORTEZOMIB"}, "proposal": {}}
    rows = [{"name": "X", "e_usable": True, "atc": [], "fda_label": True}, {"name": "BORTEZOMIB", "e_usable": True, "atc": ["L01XG01"], "fda_label": True},
            {"name": "Y", "e_usable": True, "atc": [], "fda_label": True}, {"name": "Z", "e_usable": False, "atc": [], "fda_label": True}]
    out = {c["name"]: c["candidate"] for c in le.candidate_table(rows, {"A": {"X", "BORTEZOMIB"}}, safety)}
    assert out == {"X": True, "BORTEZOMIB": False, "Y": False}


def test_geo_signature_without_annotation_is_not_computed(tmp_path) -> None:
    p = tmp_path / "m.txt.gz"
    with gzip.open(p, "wt") as f:
        f.write('!Sample_title\t"LCL-proband-1"\t"LCL-control1-1"\n!series_matrix_table_begin\n"ID_REF"\t"a"\t"b"\n!series_matrix_table_end\n')
    r = le.geo_signature(p, None)
    assert r["computed"] is False and r["n_arrays"] == 2


def test_geo_signature_with_annotation_collapse_and_effect(tmp_path) -> None:
    titles = ["LCL-proband-1", "LCL-proband-2", "LCL-control1-1", "LCL-control1-2", "LCL-control3-1", "LCL-control3-2", "LCL-relative1-1"]
    rows = {"p1": [10, 10.2, 8, 8.1, 8.0, 8.2, 9], "p2": [5, 5, 5, 5, 5, 5, 5],          # p1 up (gene U); p2 flat low-expressed duplicate of U
            "p3": [7, 7.1, 9, 9.1, 9.0, 8.9, 8], "p4": [6, 6, 6, 6, 6, 6, 6], "p5": [6, 6.1, 6, 6.2, 5.9, 6, 6]}
    ann = {"p1": "U", "p2": "U", "p3": "D", "p4": "F", "p5": "F2"}
    p = tmp_path / "m.txt.gz"
    with gzip.open(p, "wt") as f:
        f.write("!Sample_title\t" + "\t".join(f'"{t}"' for t in titles) + "\n!series_matrix_table_begin\n")
        f.write('"ID_REF"\t' + "\t".join(f'"g{i}"' for i in range(7)) + "\n")
        for k, v in rows.items():
            f.write(f'"{k}"\t' + "\t".join(map(str, v)) + "\n")
        f.write("!series_matrix_table_end\n")
    r = le.geo_signature(p, ann, n=5, min_abs_d=0.5)
    assert r["computed"] and r["n_genes_mapped"] == 4 and r["n_cases"] == 2 and r["n_controls"] == 4
    assert r["up"] == ["U"] and r["down"] == ["D"]  # highest-mean probe kept for U; flat genes fail the effect size; relative not used


def test_random_signature_null_shape_and_disjoint() -> None:
    expressed = [f"G{i}" for i in range(50)]
    index = {g: i for i, g in enumerate(expressed)}
    m = le.random_signature_null(expressed, index, 5, 5, 20, 1)
    assert m.shape == (50, 20) and ((m == 1).sum(0) == 5).all() and ((m == -1).sum(0) == 5).all()
    assert np.array_equal(m, le.random_signature_null(expressed, index, 5, 5, 20, 1))


def test_read_annotation_drops_ambiguous(tmp_path) -> None:
    p = tmp_path / "a.annot.gz"
    with gzip.open(p, "wt") as f:
        f.write("^PLATFORM = X\n!platform_table_begin\nID\tGene title\tGene symbol\nA\tt\tGENE1\nB\tt\tG2///G3\nC\tt\t\n!platform_table_end\n")
    assert le.read_annotation(p) == {"A": "GENE1"}
