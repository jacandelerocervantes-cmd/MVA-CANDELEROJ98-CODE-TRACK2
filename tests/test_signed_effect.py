from track2.signed_effect import STRENGTHEN, UNCLEAR, WEAKEN, SignedNetwork, call_of, drug_calls, validate

HEADER = "source\ttarget\tsource_genesymbol\ttarget_genesymbol\tis_directed\tis_stimulation\tis_inhibition\n"


def _net(tmp_path, edges):
    rows = "".join(f"{s}\t{t}\t{s}\t{t}\tTrue\t{a}\t{b}\n" for s, t, a, b in edges)
    (tmp_path / "o.tsv").write_text(HEADER + rows, encoding="utf-8")
    return SignedNetwork(tmp_path / "o.tsv")


def test_direct_readout_nodes_and_signs(tmp_path):
    net = _net(tmp_path, [("BUB1B", "X", "True", "False"), ("CDC20", "X", "True", "False")])
    assert call_of(net.effect("MAD2L1", -1) if "MAD2L1" in net.index else net.effect("BUB1B", -1)) == WEAKEN  # loss of a checkpoint-complex protein
    assert call_of(net.effect("CDC20", -1)) == STRENGTHEN  # inhibiting the APC/C activator prolongs the arrest


def test_signed_paths_and_shortest_layer(tmp_path):
    # KIN stimulates BUB1B directly (+) and inhibits CDC20 via a 2-step path (-,-)=(+) on APC/C (-1): both say inhibiting KIN weakens or not
    net = _net(tmp_path, [("KIN", "BUB1B", "True", "False"), ("KIN", "MID", "False", "True"), ("MID", "CDC20", "False", "True")])
    # inhibition of KIN: BUB1B (+1 readout) x (-1) x (+1 path) = -1 ; CDC20 (-1 readout) x (-1) x (+1 two inhibitions) = +1 -> net 0
    assert net.effect("KIN", -1) == 0 and call_of(0.0) == UNCLEAR
    net2 = _net(tmp_path, [("KIN", "BUB1B", "True", "False")])
    assert call_of(net2.effect("KIN", -1)) == WEAKEN and call_of(net2.effect("KIN", 1)) == STRENGTHEN


def test_edges_with_both_or_no_sign_are_dropped_and_unknown_nodes_give_none(tmp_path):
    net = _net(tmp_path, [("A", "BUB1B", "True", "True"), ("B", "BUB1B", "False", "False"), ("C", "BUB1B", "True", "False")])
    assert net.effect("A", -1) is None and net.effect("B", -1) is None and net.effect("C", -1) == -1.0


def test_validation_counts_recovered_controls(tmp_path):
    net = _net(tmp_path, [("TTK", "BUB1B", "True", "False"), ("BUB1", "BUB1B", "True", "False"), ("AURKB", "BUB3", "True", "False"), ("CDC20", "X", "True", "False"),
                          ("MAD2L1", "X", "True", "False")])
    v = validate(net)
    by = {r["gene"]: r for r in v["controls"]}
    assert by["TTK"]["recovered"] and by["BUB1"]["recovered"] and by["AURKB"]["recovered"] and not by["TRIP13"]["recovered"]
    assert by["TRIP13"]["got"] == "not_in_network" and v["n_recovered"] == 5 and v["used"]


def test_drug_calls_use_action_type(tmp_path):
    net = _net(tmp_path, [("KIN", "BUB1B", "True", "False")])
    universe = {"targets": {"T1": {"uniprot": ["KIN"]}},
                "drugs": [{"drug_id": "D1", "name": "inh", "mechanisms": [{"action": "INHIBITOR", "target_id": "T1"}]},
                          {"drug_id": "D2", "name": "ago", "mechanisms": [{"action": "AGONIST", "target_id": "T1"}]},
                          {"drug_id": "D3", "name": "none", "mechanisms": [{"action": None, "target_id": "T1"}]}]}
    calls = drug_calls(net, universe)
    assert calls["D1"]["call"] == WEAKEN and calls["D2"]["call"] == STRENGTHEN and calls["D3"]["call"] == UNCLEAR
