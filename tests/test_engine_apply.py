"""Offline, synthetic, fast tests of the frozen-engine application (no network, no real data)."""
import numpy as np

from track2 import engine_apply as ap

NAMES = [f"drug{i:02d}" for i in range(12)] + ["sirolimus", "phenylbutanoic acid"]


def test_curated_and_classes():
    ls = {"genes": [{"symbol": "B", "tier": "curated", "mechanism_class": "x"}, {"symbol": "A", "tier": "curated", "mechanism_class": "x"},
                    {"symbol": "C", "tier": "proposed, not curated", "mechanism_class": "y"}]}
    cur, prop = ap.curated_genes(ls)
    assert [g["symbol"] for g in cur] == ["A", "B"] and [g["symbol"] for g in prop] == ["C"]
    assert ap.genes_by_class(cur) == {"x": ["A", "B"]}


def test_low_literature_flag():
    assert ap.low_literature(10, 400, 29)["flag"]
    assert ap.low_literature(300, 30, 29)["flag"]
    f = ap.low_literature(300, 5000, 29)
    assert not f["flag"] and f["hit_count_stratum"] == "H"
    assert ap.low_literature(29, None, 29)["flag"] and not ap.low_literature(30, None, 29)["flag"]


def test_dilution_threshold():
    d = {"heldout": {"table_disease_level_bootstrap_sensitivity": {"fractions": {
        "0.1": {"n_papers_mean": 29.0, "B0": {"auc_gap_vs_popularity": {"estimate": 0.09}}},
        "1.0": {"n_papers_mean": 288.9, "B0": {"auc_gap_vs_popularity": {"estimate": 0.12}}}}}}}
    thr, gains = ap.dilution_threshold(d)
    assert thr == 29 and gains["1.0"]["B0_auc_gap_vs_popularity"] == 0.12


def test_oncology_flag():
    assert ap.oncology_flag(["L01XE01"], False) and ap.oncology_flag([], True) and not ap.oncology_flag(["A10BA02"], False)


def test_nine_resolution_and_positions():
    nine = ap.resolve_nine(NAMES)
    assert nine["sirolimus"] == "sirolimus" and nine["sodium phenylbutyrate"] == "phenylbutanoic acid" and nine["bosutinib"] is None
    score = np.arange(len(NAMES), 0, -1, dtype=float)
    masks = {"universe": np.ones(len(NAMES), dtype=bool), "fda_small_molecule": np.array([i % 2 == 0 for i in range(len(NAMES))])}
    pos = ap.positions(score, NAMES, masks, nine)
    assert pos["sirolimus"]["rank_universe"] == 13 and pos["sirolimus"]["rank_fda_small_molecule"] == 7 and pos["sodium phenylbutyrate"]["rank_fda_small_molecule"] is None
    assert pos["bosutinib"] == {"universe_name": None}
    assert pos["sodium phenylbutyrate"]["rank_universe"] == 14 and pos["sodium phenylbutyrate"]["of_fda_small_molecule"] == 7


def test_pool_and_gap_block():
    borrowing = {"classes": {"cell_cycle_mitotic_checkpoint": {"k": 3, "mu": 0.69, "tau2": 0.06}}}
    pool = ap.pool_for({"spindle-assembly checkpoint": ["A", "B", "C"], "centrosome/centriole": ["D", "E"]}, borrowing)
    assert pool["pool"] == "cell_cycle_mitotic_checkpoint" and pool["pool_k_diseases"] == 3
    assert ap.pool_for({"centrosome/centriole": ["D"]}, borrowing)["pool"] is None
    g = ap.evidence_gap_block("mva", "unknown", ["ORPHA:1"], pool, {"0_4": 5, "5_11": 0, "12_17": 2}, 20, "one known case")
    assert "N-of-1" in g["gaps"][0]["design_from_borrowing_logic"] and len(g["gaps"]) == 2 + 3
    g2 = ap.evidence_gap_block("fanconi", "unknown", [], {"pool": None, "note": "n"}, {}, 20, None)
    assert "registries" in g2["gaps"][0]["design_from_borrowing_logic"]


def test_approval_state_and_band_gaps():
    assert ap.approval_state(True, None) == "openFDA label match"
    assert "failed" in ap.approval_state(False, {"eligible": None})
    assert ap.approval_state(False, {"eligible": True}) == "Drugs@FDA Prescription"
    sup = [{"bands": {"0_4": "labelled", "5_11": "verify", "12_17": "partly (from 14 y)"}},
           {"bands": {"0_4": "not_labelled", "5_11": "age_not_stated", "12_17": "labelled"}}]
    assert ap.band_gap_counts(sup) == {"0_4": 1, "5_11": 2, "12_17": 0}


def test_orphanet_tier_without_cache(tmp_path):
    t, note = ap.orphanet_tier(["ORPHA:1"], tmp_path)
    assert t == "unknown" and "cache" in note
