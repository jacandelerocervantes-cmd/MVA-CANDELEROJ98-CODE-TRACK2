"""Offline synthetic tests of the computational follow-ups (CF1 to CF9)."""
import numpy as np
from scipy import sparse

from track2 import computational_followups as cf


def test_bh_and_jsonable():
    q = cf.bh([0.01, 0.04, 0.5])
    assert q[0] <= q[1] <= q[2] and abs(q[0] - 0.03) < 1e-12
    assert cf.jsonable({"a": np.float64("nan"), "b": np.int64(3), "c": {1, 2}}) == {"a": None, "b": 3, "c": [1, 2]}


def test_cf1_sign_convention():
    rng = np.random.default_rng(0)
    score = rng.normal(size=400)
    y_sens = 0.3 * score + rng.normal(scale=0.5, size=400)       # low score -> lower viability = sensitive; rho > 0
    y_res = -0.3 * score + rng.normal(scale=0.5, size=400)       # opposite
    kill = -2 + rng.normal(scale=0.2, size=400)
    doc = {}
    r = cf.cf1_sign_audit({"SENS": y_sens, "RES": y_res, "KILL": kill}, score, ["SENS", "RES", "NOPE"], {"SENS": 0.3, "RES": -0.3}, ["KILL"], doc)
    assert r["orientation_sanity_pass"]
    assert r["per_drug"]["SENS"]["delta_weak_minus_strong"] < 0 and r["per_drug"]["SENS"]["verdict"] == "survives"
    assert r["per_drug"]["RES"]["verdict"] == "refuted"
    assert r["per_drug"]["NOPE"]["verdict"] == "not estimable"


def test_cf2_hypergeometric():
    est = [f"D{i}" for i in range(100)] + ["PACLITAXEL", "VINCRISTINE", "COLCHICINE", "ERIBULIN"]
    members = ["Paclitaxel", "Vincristine", "Colchicine", "Eribulin", "Cabazitaxel"]
    enriched = cf.cf2_enrichment(est, ["PACLITAXEL", "VINCRISTINE", "COLCHICINE", "ERIBULIN", "D1", "D2"], members, [])
    assert enriched["K_tubulin_class_in_estimable"] == 4 and enriched["verdict"] == "refuted"
    plain = cf.cf2_enrichment(est, ["D1", "D2", "D3", "PACLITAXEL"], members, [])
    assert plain["verdict"] == "survives"
    assert cf.cf2_enrichment(est, ["D1"], members, [])["verdict"] == "not estimable"


def test_cf3_specificity_and_query():
    q = cf.epmc_query("olaparib", "BUB1B")
    assert 'TITLE:"olaparib"' in q and "PUB_YEAR:[1900 TO 2024]" in q and "cursorMark=*" in cf.epmc_url(q)
    deg = {f"GENE{i}": 100 + i for i in range(200)}
    deg["BUB1B"] = 150
    n = cf.matched_null_genes(deg, "BUB1B", 20)
    assert len(n["genes"]) == 20 and "BUB1B" not in n["genes"]
    res = cf.cf3_specificity({"HIGH": 50, "LOW": 0, "NA": None}, {"HIGH": [1] * 99, "LOW": [3] * 99, "NA": []})
    assert res["per_drug"]["HIGH"]["verdict"] == "survives" and res["per_drug"]["LOW"]["verdict"] == "refuted" and res["per_drug"]["NA"]["verdict"] == "not estimable"


def test_cf4_hub_free_paths():
    # 0 drug - 1 hub - 4 query (hub path); 0 - 2 - 3 - 5 query (hub-free); drug 6 only via hub
    n = 30
    edges = [(0, 1), (1, 4), (0, 2), (2, 3), (3, 5), (6, 1)] + [(1, 10 + k) for k in range(15)]
    a = sparse.coo_matrix((np.ones(len(edges)), ([e[0] for e in edges], [e[1] for e in edges])), shape=(n, n))
    adj = (a + a.T).tocsr()
    names = [f"n{i}" for i in range(n)]
    types = ["x"] * n
    r = cf.cf4_path_audit(adj, names, types, {"A": [0], "B": [6], "C": []}, {4: "q1", 5: "q2"}, hub_percentile=99)
    assert r["per_drug"]["A"]["verdict"] == "survives" and r["per_drug"]["A"]["shortest_length_hub_free"] == {"q2": 3}
    assert r["per_drug"]["B"]["verdict"] == "refuted" and r["per_drug"]["B"]["shortest_length_hubs_allowed"]["q1"] == 2
    assert r["per_drug"]["C"]["verdict"] == "not estimable"


def test_cf5_scheme_table_and_verdict():
    scores = {"m1": {"700": {"a": -3.0, "b": -2.0, "c": -1.0}}, "m2": {"700": {"a": 1.0, "b": -3.0, "c": -2.5}}}
    wts = {"W1": {"m1": 1.0, "m2": 0.1}, "W2": {"m1": 1.0, "m2": 1.0}, "W3": {"m1": 1.0}}
    t = cf.scheme_table(scores, wts, ["700"], ["W1", "W2", "W3"], top_n=1)
    assert t["top_set"]["W3"] == {"a"} and t["top_set"]["W2"] == {"b"}
    res = cf.cf5_verdicts({"A": "a", "B": "b", "Z": "z"}, t, ["W1", "W2", "W3"], "W3", {"A": 1, "B": 1, "Z": None}, {"A": [5] * 99, "B": [5] * 99, "Z": []})
    assert res["per_drug"]["A"]["criterion_in_top_set_under_other_scheme"] is True  # a is top under W1
    assert res["per_drug"]["B"]["verdict"] == "survives" and res["per_drug"]["Z"]["verdict"] == "not estimable"


def _net():
    import types

    nodes = ["KIN", "BUB1B", "CDC20", "MID"] + [f"X{i}" for i in range(12)]
    index = {n: i for i, n in enumerate(nodes)}
    edges = [("KIN", "BUB1B", 1), ("MID", "CDC20", -1)] + [(f"X{i}", f"X{i + 1}", 1) for i in range(11)] + [("KIN", "MID", 1)]
    rows = [index[s] for s, _, _ in edges]
    cols = [index[t] for _, t, _ in edges]
    adj = sparse.csr_matrix(([float(v) for _, _, v in edges], (rows, cols)), shape=(len(nodes), len(nodes)))
    return types.SimpleNamespace(nodes=nodes, index=index, adj=adj, readout={"BUB1B": 1.0, "CDC20": -1.0}, symbol={n: n for n in nodes}, by_symbol={n: n for n in nodes})


def test_cf6_sign_permutation_runs_and_controls():
    from track2.signed_effect import STRENGTHEN, WEAKEN, call_of

    net = _net()
    drugs = {"D1": {"uniprot_by_target": {"t1": ["KIN"]}, "mechanisms": [{"target_id": "t1", "action": "INHIBITOR"}], "stored_net": None},
             "D0": {"uniprot_by_target": {"t2": ["X0"]}, "mechanisms": [{"target_id": "t2", "action": "INHIBITOR"}], "stored_net": None}}
    res = cf.cf6_sign_permutation(net, drugs, {"INHIBITOR"}, {"ACTIVATOR"}, 50, controls=(("KIN", -1, WEAKEN), ("CDC20", -1, STRENGTHEN)), call_of=call_of, readout_genes=["CDC20"])
    assert res["per_drug"]["D1"]["observed_net"] != 0 and res["per_drug"]["D1"]["p_two_sided"] <= 1.0
    assert res["per_drug"]["D0"]["verdict"] == "not estimable"
    assert res["controls"]["of"] == 2 and 0 <= res["controls"]["p_null_ge_observed"] <= 1


def test_cf7_signature_leave_out_and_verdicts():
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(300)]
    m = rng.normal(8, 0.2, size=(300, 6))
    m[:20, :2] += 2.0       # up in the 2 case arrays (columns 0, 1)
    m[20:40, :] += 4.0
    m[20:40, :2] -= 2.0
    full = cf.signature_from_arrays(m, genes, [0, 1], [2, 3, 4, 5], n=10)
    assert set(full["up"]) <= {f"G{i}" for i in range(20)} and len(full["up"]) == 10
    one_case = cf.signature_from_arrays(m, genes, [0], [2, 3, 4, 5], n=10)  # a single case array must not give NaN
    assert len(one_case["up"]) == 10 and len(one_case["down"]) == 10
    per = {"lo1": {"A": {"R": 5, "q_rev": 0.05}, "B": {"R": 5, "q_rev": 0.05}, "C": {"R": 5, "q_rev": 0.05}},
           "lo2": {"A": {"R": 5, "q_rev": 0.05}, "B": {"R": 5, "q_rev": 0.2}, "C": {"R": 5, "q_rev": 0.05}}}
    conc = {"A": {"concordant_with_disease_signature": False}, "B": {"concordant_with_disease_signature": False}, "C": {"concordant_with_disease_signature": True}}
    res = cf.cf7_verdicts(["A", "B", "C", "D"], per, conc)["per_drug"]
    assert res["A"]["verdict"] == "survives" and res["B"]["verdict"] == "refuted" and res["C"]["action_if_refuted"] == "flag hazard" and res["D"]["verdict"] == "not estimable"


def test_cf8_distance_close_versus_far():
    n = 200
    rng = np.random.default_rng(2)
    edges = [(i, i + 1) for i in range(n - 1)]                 # a path graph: distances are index differences
    edges += [(int(a), int(b)) for a, b in rng.integers(0, n, size=(100, 2)) if a != b]
    a = sparse.coo_matrix((np.ones(len(edges)), ([e[0] for e in edges], [e[1] for e in edges])), shape=(n, n))
    adj = (a + a.T).tocsr()
    adj.data[:] = 1.0
    degree = np.asarray(adj.sum(axis=1)).ravel()
    res = cf.cf8_distance(adj, degree, 10, {"T1": 11, "T2": 12}, 200, chunk=50)
    assert res["observed_mean_distance"] == 1.5 and res["verdict"] == "survives"
    far = cf.cf8_distance(adj, degree, 10, {"T1": 190, "T2": 195}, 200, chunk=50)
    assert far["verdict"] == "refuted"


def test_cf9_not_runnable_and_second_screen():
    rows = [{"Drug.Name": "DRUGA", "Synonyms": "", "screen": "REP.PRIMARY"}, {"Drug.Name": "DRUGB", "Synonyms": "DRUGB2", "screen": "REP.1M"}]
    r = cf.cf9_second_screen(rows, ["DRUGA", "DRUGB"])
    assert r["verdict"] == "not runnable" and "needed" in r
    rows.append({"Drug.Name": "DRUGA", "Synonyms": "", "screen": "REP.1M"})
    r2 = cf.cf9_second_screen(rows, ["DRUGA", "DRUGB"], rerun=lambda d, s: {"same_sign": True, "p_perm": 0.01})
    assert r2["per_drug_verdict"] == {"DRUGA": "survives"}


def test_candidate_status_rules():
    plan = {"triage": [{"candidate": "Alpha", "lines": ["C"], "next_step": "CF1, then Lab X4", "next_step_type": "lab", "experiment_id": "X4"},
                       {"candidate": "Beta", "lines": ["A"], "next_step": "CF4", "next_step_type": "comp", "experiment_id": None},
                       {"candidate": "Gamma", "lines": ["C"], "next_step": "Lab H1; CF2", "next_step_type": "hazard", "experiment_id": "H1"}]}
    results = {"CF1": {"per_drug": {"ALPHA": {"verdict": "survives"}}}, "CF4": {"per_drug": {"BETA": {"verdict": "refuted"}}}, "CF2": {"verdict": "refuted", "concerns": ["GAMMA"]}}
    s = {r["candidate"]: r for r in cf.candidate_status(plan, results)}
    assert s["Alpha"]["next_step"] == "laboratory X4" and "passed CF1" in s["Alpha"]["status"]
    assert s["Beta"]["status"].startswith("closed as a computational candidate")
    assert "hazard class" in s["Gamma"]["status"]
