"""Reversal layer tests: offline, SigCom LINCS calls are mocked (no network, no patient data)."""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from track2 import reversal, reversal_names, reversal_stats as st, sigcom_client as sc


def test_clean_removes_salts_and_punctuation() -> None:
    assert reversal_names.clean("IXAZOMIB CITRATE") == "IXAZOMIB"
    assert reversal_names.clean("Imatinib mesylate") == "IMATINIB"
    assert reversal_names.clean("5-fluorouracil") == "5FLUOROURACIL"
    assert reversal_names.clean("Foo (hydrate) HCl") == "FOO"
    index = reversal_names.drug_index([{"name": "BORTEZOMIB", "drug_id": "A"}, {"name": "X SODIUM", "drug_id": "B"}])
    assert reversal_names.match("bortezomib", index) == ["A"]
    assert reversal_names.match("x", index) == ["B"]
    assert reversal_names.match("unknown", index) == []


def test_per_drug_takes_best_reversal_and_mimic() -> None:
    rows = [
        {"uuid": "s1", "z-sum": -4.0, "fdr-up": 0.01, "fdr-down": 0.02},
        {"uuid": "s2", "z-sum": -6.0, "fdr-up": 0.001, "fdr-down": 0.03},
        {"uuid": "s3", "z-sum": 5.0, "fdr-up": 0.5, "fdr-down": 0.5},
        {"uuid": "s4", "z-sum": -9.0, "fdr-up": 0.0, "fdr-down": 0.0},  # not an approved drug
    ]
    out = st.per_drug(rows, {"s1": ["D1"], "s2": ["D1"], "s3": ["D1"]})
    assert set(out) == {"D1"}
    assert out["D1"]["rev"] == 6.0 and out["D1"]["mim"] == 5.0
    assert out["D1"]["api_fdr_rev"] == 0.03  # conservative: the larger of up/down FDR


def test_empirical_p_and_bh() -> None:
    null = np.arange(100, dtype=float)
    p = st.empirical_p(np.array([0.0, 50.0, 1000.0]), null)
    assert p[0] == pytest.approx(1.0) and p[2] == pytest.approx(1 / 101)
    assert 0.4 < p[1] < 0.6
    fdr = st.bh_fdr(np.array([0.01, 0.04, 0.03, 0.5]))
    assert fdr[0] == pytest.approx(0.04) and fdr[3] == pytest.approx(0.5) and (fdr <= 1).all()
    assert st.direction(0.01, 0.9) == "reverses" and st.direction(0.9, 0.01) == "mimics" and st.direction(0.5, 0.5) == "none"


def _fake_sigcom(monkeypatch: pytest.MonkeyPatch) -> None:
    genes = [f"G{i}" for i in range(400)]
    monkeypatch.setattr(sc, "consensus_signature", lambda kind, name: {g: (i - 200) / 1000.0 for i, g in enumerate(genes)})
    monkeypatch.setattr(sc, "gene_ids", lambda symbols: {s: "id-" + s for s in symbols})
    monkeypatch.setattr(sc, "library_signature_names",
                        lambda lib, fields: ({f"k{i}": {"pert_name": f"GENE{i}"} for i in range(30)} if lib == "xpr"
                                         else {f"sig{i}": {"pert_name": n} for i, n in enumerate(["drugA", "drugB", "drugC", "noise"])}))
    monkeypatch.setattr(sc, "signature_meta", lambda ids, fields: {u: {"pert_name": "drugA", "cell_line": "MCF7"} for u in ids})

    def two_sided(up: list[str], down: list[str], database: str, limit: int) -> list[dict[str, Any]]:
        if database == sc.LIB_CHEM_CELL:
            return [{"uuid": "cell1", "z-sum": -7.0, "fdr-up": 0.01, "fdr-down": 0.01}]
        planted = up and up[0] == "id-G399"  # the real query: drugA reverses, drugB mimics
        rng = random.Random(len(up) + hash(tuple(up)) % 1000)
        z = [-12.0, 12.0, 0.0] if planted else [rng.gauss(0, 1) for _ in range(3)]
        return [{"uuid": f"sig{i}", "z-sum": z[i], "fdr-up": 0.01, "fdr-down": 0.01} for i in range(3)]

    monkeypatch.setattr(sc, "two_sided", two_sided)
    monkeypatch.setattr(reversal, "N_PERMUTATIONS", 30)
    monkeypatch.setattr(reversal, "N_KO_CONTROLS", 3)
    monkeypatch.setattr(sc, "LIB_CRISPR_CONSENSUS", "xpr")


def test_run_end_to_end_with_mocked_api(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _fake_sigcom(monkeypatch)
    universe = {"drugs": [{"name": n, "drug_id": i} for n, i in [("DRUGA", "A"), ("DRUGB HCL", "B"), ("DRUGC", "C"), ("OTHER", "D")]]}
    (tmp_path / "u.json").write_text(json.dumps(universe), encoding="utf-8")
    (tmp_path / "c.json").write_text(json.dumps([{"drug_id": "A"}]), encoding="utf-8")
    meta = reversal.run(tmp_path / "u.json", tmp_path / "c.json", tmp_path / "out")
    scores = json.loads((tmp_path / "out" / "reversal_scores.json").read_text(encoding="utf-8"))
    assert meta["n_drugs_with_score"] == 3
    assert scores["D"] == {"available": False, "name": "OTHER",
                           "reason": "no LINCS L1000 chemical consensus signature with a matching perturbagen name"}
    assert scores["A"]["available"] and scores["A"]["score"] == 12.0 and scores["A"]["direction"] == "reverses"
    assert scores["A"]["cell_lines"] == ["MCF7"] and scores["A"]["in_candidates"] is True
    assert scores["B"]["direction"] == "mimics" and scores["B"]["score"] == 0.0
    assert scores["C"]["direction"] == "none"
    assert len(meta["null"]["false_reverses_calls_per_random_set_leave_one_out"]) == 30
    assert (tmp_path / "out" / "reversal_summary.md").read_text(encoding="utf-8").startswith("# Transcriptomic reversal")
    assert scores["A"]["ko_control"]["n_controls"] == 3
    assert "proxy_sac" in scores["A"] and set(scores["A"]["proxy_sac"]) == set(reversal.PROXIES)
