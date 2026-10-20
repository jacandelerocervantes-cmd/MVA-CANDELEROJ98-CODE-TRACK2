import numpy as np

from track2 import cf9


def test_match_drugs_exact_and_synonym():
    m = cf9.match_drugs(["Paclitaxel", "Vinblastine", "HKI-272"], {"PACLITAXEL": set(), "NERATINIB": {"HKI272"}, "LAZERTINIB": {"YH25448"}})
    assert m["PACLITAXEL"] == ["Paclitaxel"] and m["NERATINIB"] == ["HKI-272"] and m["LAZERTINIB"] == []


def test_match_lines_priority_and_ambiguity():
    models = [{"ModelID": "ACH-1", "SangerModelID": "SIDM1", "COSMICID": "11", "StrippedCellLineName": "AAA", "CellLineName": "AAA"},
              {"ModelID": "ACH-2", "SangerModelID": "", "COSMICID": "22.0", "StrippedCellLineName": "BBB", "CellLineName": "B-BB"},
              {"ModelID": "ACH-3", "SangerModelID": "", "COSMICID": "", "StrippedCellLineName": "CCC", "CellLineName": "CCC"}]
    g = [{"cosmic": 1.0, "sanger": "SIDM1", "name": "x"}, {"cosmic": 22.0, "sanger": "", "name": "x"}, {"cosmic": 3.0, "sanger": "", "name": "c-cc"}, {"cosmic": 4.0, "sanger": "", "name": "zzz"}]
    out, cov = cf9.match_lines(g, models)
    assert out == {1.0: "ACH-1", 22.0: "ACH-2", 3.0: "ACH-3"} and cov["unmatched"] == 1 and cov["by_sanger_id"] == 1 and cov["by_cosmic_id"] == 1 and cov["by_name"] == 1


def test_verdict_and_not_estimable():
    assert cf9.verdict({"estimable": True, "rho_adj": 0.2, "p_perm": 0.01}, 0.1) == "survives"
    assert cf9.verdict({"estimable": True, "rho_adj": -0.2, "p_perm": 0.01}, 0.1) == "refuted"
    assert cf9.verdict({"estimable": True, "rho_adj": 0.2, "p_perm": 0.2}, 0.1) == "refuted"
    assert cf9.verdict({"estimable": False}, 0.1) == "not estimable"


def test_run_drug_synthetic():
    rng = np.random.default_rng(1)
    n = 150
    score = rng.normal(size=n)
    lin = np.array(["A"] * 75 + ["B"] * 75, dtype=object)
    prolif = rng.normal(size=n)
    r = cf9.run_drug(0.8 * score + rng.normal(size=n), score, prolif, lin, n_perm=99)
    assert r["estimable"] and r["rho_adj"] > 0.3 and r["p_perm"] <= 0.05
    assert not cf9.run_drug(score[:50], score[:50], prolif[:50], lin[:50])["estimable"]
