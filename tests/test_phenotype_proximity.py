from track2.phenotype_proximity import robustness


def test_robustness_counts_top_negative_z_across_configurations():
    # 9 configurations: drug A is the best in all of them, B only in 3, C has z > 0 (never counted)
    table = {f"v{i}": {s: {"A": -3.0, "B": (-2.0 if i == 0 else 1.0), "C": 0.5} for s in ("400", "700", "900")} for i in range(3)}
    res = robustness(table, {"A": "drugA", "B": "drugB", "C": "drugC"})
    by = {d["name"]: d for d in res["drugs"]}
    assert res["rule"]["configurations"] == 9 and by["drugA"]["robust"] and by["drugA"]["robustness_count"] == 9
    assert by["drugB"]["robustness_count"] == 3 and not by["drugB"]["robust"] and "drugC" not in by
