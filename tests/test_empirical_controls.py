from track2 import empirical_controls as ec


def test_positive_control_consistent():
    pc = ec.positive_control()
    assert pc["M"] == len(pc["positive_set_drugs"])
    assert set(pc["recovered"]) | set(pc["not_recovered"]) == set(pc["positive_set_drugs"])
    assert pc["N_recovered"] == len(pc["recovered"])
    assert all(e["verified_in_verification_json"] for e in pc["entries"])
    for e in pc["entries"]:
        assert (e["why_not_in_final_output"] is None) == e["in_final_caution_output"]


def test_verifier_counts_add_up():
    r = ec.verifier_catch_rate()
    t = r["total"]
    assert t["verified"] + t["not_found"] + t["not_checkable"] == t["n"] == int(r["n_sentences"])
    assert sum(b["n"] for b in r["per_topic"].values()) == t["n"]
    assert len(r["not_found_entries"]) == t["not_found"]
    assert all(len(x["sentence_first_100"]) <= 100 for x in r["not_found_entries"])
    q = r["label_quote_check"]
    assert q["found"] + q["not_found"] + q["not_checkable"] == q["n_quotes"]
