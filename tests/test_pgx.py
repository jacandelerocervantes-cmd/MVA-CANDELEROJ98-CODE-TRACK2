"""Offline tests of the star-allele-lite logic with a synthetic CPIC gene and a tiny VCF."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

from track2 import pgx
from track2.pgx_vcf import scan_vcf

# gene GENEX on chr1: pos 100 A>G (rs1), 200 C>T (rs2), 300 delT (not interpretable)
CPIC_GENE = {
    "gene": {"symbol": "GENEX", "chr": "chr1", "chromosequenceid": "NC_000001.11"},
    "alleles": [
        {"id": 1, "name": "*1", "structuralvariation": False},
        {"id": 2, "name": "*2", "structuralvariation": False},
        {"id": 3, "name": "*3", "structuralvariation": False},
        {"id": 4, "name": "*4", "structuralvariation": False},
        {"id": 5, "name": "*5", "structuralvariation": False},
        {"id": 6, "name": "*6", "structuralvariation": True},
    ],
    "locations": {
        10: {"id": 10, "chromosomelocation": "g.100A>G", "dbsnpid": "rs1"},
        11: {"id": 11, "chromosomelocation": "g.200C>T", "dbsnpid": "rs2"},
        12: {"id": 12, "chromosomelocation": "g.300del", "dbsnpid": "rs3"},
    },
    "values": [
        {"alleledefinitionid": 2, "locationid": 10, "variantallele": "G"},
        {"alleledefinitionid": 3, "locationid": 10, "variantallele": "G"},
        {"alleledefinitionid": 3, "locationid": 11, "variantallele": "T"},
        {"alleledefinitionid": 4, "locationid": 11, "variantallele": "Y"},  # C or T: reference allowed, not required
        {"alleledefinitionid": 5, "locationid": 10, "variantallele": "G"},
        {"alleledefinitionid": 5, "locationid": 12, "variantallele": "delT"},
    ],
}


def _vcf(tmp_path: Path, lines: list[str]) -> Path:
    path = tmp_path / "t.vcf.gz"
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(header + "".join(line + "\n" for line in lines))
    return path


def test_parse_location() -> None:
    assert pgx.parse_location("g.100A>G") == (100, "A", {"G"})
    assert pgx.parse_location("g.5G>A; g.5G>T") == (5, "G", {"A", "T"})
    assert pgx.parse_location("g.7_8insT") is None and pgx.parse_location("g.9TA[7]; g.9TA[8]") is None


def test_scan_accepts_both_chromosome_namings_and_ignores_indels(tmp_path: Path) -> None:
    vcf = _vcf(tmp_path, ["1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1", "chr1\t200\t.\tC\tCA\t50\tPASS\t.\tGT\t0/1"])
    calls = scan_vcf(vcf, {"chr1": [100, 200]})
    assert calls["1"][100][0]["alt_bases_called"] == ["G"]
    assert calls["1"].get(200, [{"alt_bases_called": []}])[0]["alt_bases_called"] == []


def test_alleles_and_classification() -> None:
    _, alleles = pgx.build_alleles(CPIC_GENE)
    assert alleles["*2"]["required"] == {100: {"G"}} and alleles["*4"]["required"] == {}
    assert alleles["*5"]["unassessable"] == 1
    detected = {100: {"G"}}
    status = {n: pgx.classify(a, detected) for n, a in alleles.items()}
    assert status == {"*1": "reference_like", "*2": "possible", "*3": "not_supported", "*4": "reference_like",
                      "*5": "unresolved", "*6": "unresolved"}
    both = {100: {"G"}, 200: {"T"}}
    names = [n for n, a in alleles.items() if pgx.classify(a, both) == "possible"]
    assert names == ["*2", "*3"] and pgx.most_specific(names, alleles) == ["*3"]


def test_no_call_is_not_reference(tmp_path: Path) -> None:
    summary, detail = pgx.evaluate_gene("GENEX", CPIC_GENE, scan_vcf(_vcf(tmp_path, ["1\t999\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]), {"1": [100, 200]}))
    assert summary["n_positions_with_nonreference_call"] == 0 and summary["possible_alleles"] == []
    assert "no diplotype" in summary["limitation"] and detail["calls_at_defining_positions"] == []
    assert summary["n_positions_not_interpretable_indel_or_repeat"] == 1


def test_reference_mismatch_is_not_used(tmp_path: Path) -> None:
    calls = scan_vcf(_vcf(tmp_path, ["1\t100\t.\tC\tG\t50\tPASS\t.\tGT\t1/1"]), {"1": [100]})
    summary, _ = pgx.evaluate_gene("GENEX", CPIC_GENE, calls)
    assert summary["n_vcf_records_with_ref_mismatch_vs_cpic_grch38"] == 1 and summary["possible_alleles"] == []


def test_run_writes_no_coordinates_in_public_files(tmp_path: Path, monkeypatch: Any) -> None:
    pairs = [{"drugname": "Alphadrug", "genesymbol": "GENEX", "cpiclevel": "A", "clinpgxlevel": "1A",
              "guidelinename": "G", "guidelineurl": "https://example.org/g", "provisional": False}]

    def get(url: str) -> Any:
        if "pair_view" in url:
            return pairs
        if "select=symbol" in url:
            return [CPIC_GENE["gene"]]
        if "allele_definition" in url:
            return CPIC_GENE["alleles"]
        if "sequence_location" in url:
            return list(CPIC_GENE["locations"].values())
        return CPIC_GENE["values"]

    monkeypatch.setattr(pgx, "GENES", ["GENEX"])
    vcf = _vcf(tmp_path, ["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"])
    cands = [{"drug_id": "D1", "name": "ALPHADRUG HYDROCHLORIDE"}, {"drug_id": "D2", "name": "NOPAIR"}]
    pgx.run(vcf, cands, tmp_path / "out", tmp_path / "private" / "detail.json", get)
    assert not (tmp_path / "out" / "pgx_summary.json").exists()  # the patient-level summary is private: it goes next to the detail file, never to out/
    public = (tmp_path / "private" / "pgx_summary.json").read_text(encoding="utf-8")
    assert "100" not in json.loads(public)["genes"]["GENEX"] and '"pos"' not in public and "rs1" not in public
    assert json.loads((tmp_path / "private" / "detail.json").read_text(encoding="utf-8"))["GENEX"]["calls_at_defining_positions"][0]["pos"] == 100
    relevance = json.loads((tmp_path / "out" / "pgx_drug_relevance.json").read_text(encoding="utf-8"))
    assert relevance["D1"][0]["gene"] == "GENEX" and relevance["D2"] == []
