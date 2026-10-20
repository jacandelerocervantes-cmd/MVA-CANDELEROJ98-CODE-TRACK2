from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

from track2 import auto_modules, submission_input
from track2.string_net import Interactome

HEADER = "finding_type,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,notes\n"


def _csv(tmp_path: Path) -> Path:
    path = tmp_path / "s.csv"
    path.write_text(HEADER
                    + "primary,chr1,100,T,G,chr1,200,A,C,0.9,GENEA compound_het\n"
                    + "secondary,chr1,100,T,G,,,,,0.5,GENEA other\n"
                    + "sv,chr2,300,N,<DEL>,,,,,0.1,GENEB sv\n", encoding="utf-8")
    return path


def test_read_submission_parses_pairs_and_gene(tmp_path: Path) -> None:
    rows = submission_input.read_submission(_csv(tmp_path))
    assert [r["gene"] for r in rows] == ["GENEA", "GENEA", "GENEB"]
    assert rows[0]["variants"] == [("chr1", 100, "T", "G"), ("chr1", 200, "A", "C")]


def test_weight_is_best_epcr_and_clinvar_deduplicated(tmp_path: Path, monkeypatch: Any) -> None:
    calls = []
    monkeypatch.setattr(submission_input, "clinvar_for_variant", lambda *a: calls.append(a) or {"status": "not_in_clinvar", "classification": None})
    genes = submission_input.summarise_genes(submission_input.read_submission(_csv(tmp_path)))
    assert genes[0]["gene"] == "GENEA" and genes[0]["weight"] == 0.9 and genes[0]["is_primary"]
    assert len(genes[0]["clinvar"]) == 2  # chr1:100 appears in two rows but is listed once
    assert genes[1]["gene"] == "GENEB" and genes[1]["clinvar"] == []  # structural variants are not looked up


def test_clinvar_matches_exact_allele(monkeypatch: Any) -> None:
    def fake(url: str) -> dict[str, Any]:
        if "esearch" in url:
            return {"esearchresult": {"idlist": ["1", "2"]}}
        return {"result": {"1": {"variation_set": [{"canonical_spdi": "NC_1:99:T:A"}]},
                           "2": {"variation_set": [{"canonical_spdi": "NC_1:99:T:G"}], "title": "t",
                                 "germline_classification": {"description": "Pathogenic", "review_status": "criteria provided"}}}}
    monkeypatch.setattr(submission_input, "get_json", fake)
    out = submission_input.clinvar_for_variant("chr1", 100, "T", "G")
    assert out["classification"] == "Pathogenic" and out["clinvar_id"] == "2"
    assert submission_input.clinvar_for_variant("chr1", 100, "TT", "G")["status"] == "not_a_snv"


def _interactome(n: int, edges: list[tuple[int, int]]) -> Interactome:
    m = sparse.lil_matrix((n, n))
    for a, b in edges:
        m[a, b] = m[b, a] = 1
    adj = sparse.csr_matrix(m)
    return Interactome(proteins=[f"P{i}" for i in range(n)], adjacency=adj, degree=np.asarray(adj.sum(axis=1)).ravel(),
                       uniprot_to_node={f"U{i}": i for i in range(n)}, symbol_to_node={f"G{i}": i for i in range(n)},
                       min_score=700, string_version="test")


def test_pathway_lookup_404_means_no_pathways(monkeypatch: Any) -> None:
    def missing(url: str) -> Any:
        raise LookupError(url)
    monkeypatch.setattr(auto_modules, "get_json", missing)
    assert auto_modules.pathways_for("P1") == []


def test_small_neighbourhood_is_not_modellable(monkeypatch: Any) -> None:
    monkeypatch.setattr(auto_modules, "uniprot_accession", lambda s: None)
    net = _interactome(30, [(0, 1), (0, 2)])
    mech = auto_modules.build_mechanism("G0", 0.1, net, {})
    assert mech["module_rule"] == "string_physical_neighbours_2nd_order" and mech["seeds"]["genes"] == ["G0", "G1", "G2"]
    assert not auto_modules.modellable(mech, net)


def test_module_rule_keeps_pathways_in_size_range(monkeypatch: Any) -> None:
    monkeypatch.setattr(auto_modules, "uniprot_accession", lambda s: "U0")
    monkeypatch.setattr(auto_modules, "pathways_for", lambda acc: [
        {"id": "R-A", "name": "small", "n_participants": 5}, {"id": "R-B", "name": "ok", "n_participants": 50},
        {"id": "R-C", "name": "huge", "n_participants": 900}])
    mech = auto_modules.build_mechanism("G0", 1.0, _interactome(30, []), {})
    assert mech["seeds"]["reactome_pathways"] == ["R-B"] and mech["module_rule"] == "reactome_pathways_of_the_gene"


def test_malformed_submission_files_fail_with_a_clear_message(tmp_path: Path) -> None:
    import pytest
    bad_cols = tmp_path / "a.csv"
    bad_cols.write_text("foo,bar\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="lacks the columns"):
        submission_input.read_submission(bad_cols)
    empty_notes = tmp_path / "b.csv"
    empty_notes.write_text(HEADER + "primary,chr1,100,T,G,,,,,0.9,\n", encoding="utf-8")
    with pytest.raises(ValueError, match="notes. is empty"):
        submission_input.read_submission(empty_notes)
    no_rows = tmp_path / "c.csv"
    no_rows.write_text(HEADER, encoding="utf-8")
    with pytest.raises(ValueError, match="no rows"):
        submission_input.read_submission(no_rows)
    bad_num = tmp_path / "d.csv"
    bad_num.write_text(HEADER + "primary,chr1,xx,T,G,,,,,0.9,GENEA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not a number"):
        submission_input.read_submission(bad_num)
