import json
import random

from track2 import empirical_label_robustness as E


def row(level, direction, consequence="growth_muscle", use=None, section="findings"):
    return {"level": level, "direction": direction, "consequence": consequence, "use": use, "section": section}


DATA = {
    "SIROLIMUS": {"entries": [row(1, "may_help")]},
    "VINCRISTINE": {"entries": [row(1, "may_help", "tumour")]},
    "BADDRUG": {"entries": [row(3, "may_harm")]},
    "NOTAPPROVED": {"entries": [row(1, "may_help")]},
}
APPROVED = {"SIROLIMUS", "VINCRISTINE", "BADDRUG"}


def test_decide_baseline():
    prop, caut = E.decide(DATA, APPROVED)
    assert prop == {"SIROLIMUS": "constitutional", "VINCRISTINE": "tumour_context"}
    assert caut == {"BADDRUG"}


def test_corrupt_bounds_and_count():
    rng = random.Random(1)
    new = E.corrupt(DATA, 0.5, rng)
    changed = sum(a != b for n in DATA for a, b in zip(DATA[n]["entries"], new[n]["entries"]))
    assert changed == 2
    assert all(0 <= r["level"] <= 3 for v in new.values() for r in v["entries"])
    assert DATA["SIROLIMUS"]["entries"][0]["level"] == 1  # input not mutated


def test_simulate_deterministic_and_zero_rate():
    a = E.simulate(DATA, APPROVED, reps=50, seed=3, rates=(0.0, 0.5))
    b = E.simulate(DATA, APPROVED, reps=50, seed=3, rates=(0.0, 0.5))
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    z = a["rates"]["0"]
    assert z["expected_survivors"] == 2 and all(f == 1.0 for f in z["survival_fraction_by_drug"].values())
    assert a["rates"]["50"]["expected_survivors"] <= 2
