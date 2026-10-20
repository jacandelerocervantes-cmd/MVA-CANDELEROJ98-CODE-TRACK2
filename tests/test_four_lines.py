from track2.four_lines import combine, gates, overlap_statistics

UNIVERSE = {"drugs": [
    {"drug_id": "D1", "name": "ONE", "atc": ["L01XE"], "oral": True},
    {"drug_id": "D2", "name": "TWO", "atc": ["C01"], "oral": True},
    {"drug_id": "D3", "name": "THREE", "atc": ["N01AB"], "parenteral": True},
    {"drug_id": "D4", "name": "FOUR", "atc": ["D07"], "topical": True},
    {"drug_id": "D5", "name": "FIVE", "atc": ["C02"], "oral": True},
]}
FDA_OK = lambda name: {"eligible": name != "FIVE"}  # noqa: E731


def _status(a=True, b=True, c=False, d=False):
    return {k: {"usable": v, "evidence": ""} for k, v in zip("ABCD", (a, b, c, d))}


def test_candidate_needs_two_usable_lines_and_all_gates():
    sets = {"A": {"D1", "D2", "D3", "D4", "D5"}, "B": {"D1", "D3", "D4", "D5"}, "C": {"D2"}, "B_phenotype_module": {"D1"}, "B_mitosis_module_without_generic_complexes": {"D3", "D4", "D5"}}
    rows = {r["name"]: r for r in combine(UNIVERSE, _status(), sets, {}, FDA_OK)}
    assert rows["ONE"]["candidate"] and rows["ONE"]["tier"] == "B" and rows["ONE"]["list"] == "tumour_context" and rows["ONE"]["B_modules"] == ["B_phenotype_module"]
    assert not rows["TWO"]["candidate"] and rows["TWO"]["tier"] == "C"        # C is not usable: only line A counts for it
    assert not rows["THREE"]["candidate"] and rows["THREE"]["tier"] == "not_eligible"   # acute-use anaesthetic (ATC N01)
    assert rows["FOUR"]["tier"] == "not_eligible" and rows["FIVE"]["tier"] == "not_eligible"  # topical only; not FDA-eligible


def test_direction_changes_the_tier_only_when_line_d_is_usable():
    sets = {"A": {"D1"}, "B": {"D1"}, "C": set(), "B_phenotype_module": {"D1"}, "B_mitosis_module_without_generic_complexes": set()}
    weak = {"D1": {"call": "weakens_checkpoint"}}
    assert combine(UNIVERSE, _status(d=False), sets, weak, FDA_OK)[0]["tier"] == "B"      # exploratory D cannot move a candidate
    assert combine(UNIVERSE, _status(d=True), sets, weak, FDA_OK)[0]["tier"] == "avoid"
    assert combine(UNIVERSE, _status(d=True), sets, {"D1": {"call": "strengthens_checkpoint"}}, FDA_OK)[0]["tier"] == "A"


def test_gates_and_overlap_statistics():
    assert gates({"atc": ["N05CD01"], "oral": True}, {"eligible": True})["not_acute_use_class"] is False
    stats = overlap_statistics({"A": {"a", "b", "c"}, "B": {"c", "d"}, "C": set()}, {"a", "b", "c", "d", "e", "f", "g", "h"})
    assert stats["A_and_B"]["overlap"] == 1 and 0 < stats["A_and_B"]["p_hypergeometric"] < 1


def test_one_line_plus_accepted_literature_makes_a_candidate():
    sets = {"A": {"D1", "D2"}, "B": set(), "C": set(), "B_phenotype_module": set(), "B_mitosis_module_without_generic_complexes": set()}
    rows = {r["name"]: r for r in combine(UNIVERSE, _status(a=True, b=False), sets, {}, FDA_OK, literature={"ONE"})}
    assert rows["ONE"]["candidate"] and rows["ONE"]["evidence"] == "one_line_plus_literature" and rows["ONE"]["tier"] == "B"
    assert not rows["TWO"]["candidate"] and rows["TWO"]["evidence"] == "one_line"
    # literature never makes a candidate out of a drug that supports no usable line, nor out of one that fails a gate
    sets2 = {"A": {"D3"}, "B": set(), "C": set(), "B_phenotype_module": set(), "B_mitosis_module_without_generic_complexes": set()}
    assert not combine(UNIVERSE, _status(a=True, b=False), sets2, {}, FDA_OK, literature={"THREE"})[0]["candidate"]
