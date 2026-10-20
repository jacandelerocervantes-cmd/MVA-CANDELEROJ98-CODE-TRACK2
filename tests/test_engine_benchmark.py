import json
import random
from datetime import date

from track2 import engine_benchmark as eb

UIDX = {
    "sirolimus": {"drug_id": "C1", "name": "SIROLIMUS", "molecule_type": "Small molecule", "atc": ["L04AA10", "L04AH02"]},
    "ivacaftor": {"drug_id": "C2", "name": "IVACAFTOR", "molecule_type": "Small molecule", "atc": ["R07AX02"]},
    "imiglucerase": {"drug_id": "C3", "name": "IMIGLUCERASE", "molecule_type": "Enzyme", "atc": []},
}


def match(**kw):
    base = {"hpo": ["ORPHA:1"], "orphanet": ["ORPHA:1"], "n_matched": 1, "genes": ["A1", "B2"], "inheritance": ["Autosomal recessive"],
            "orphanet_catalogued": True, "mendelian": True, "exact_orphanet_mendelian": True}
    base.update(kw)
    return base


def test_norm_and_atc():
    assert eb.norm_name("Beta-thalassemia (intermedia)") == "beta thalassemia"
    assert eb.norm_name("Fanconi anaemia") == "fanconi anemia"
    assert eb.atc_class(["R07AX02", "A01AA01"]) == "A01A"
    assert eb.atc_class([]) is None


def test_evaluate_candidate_rules():
    gold = {"sirolimus", "ivacaftor"}
    assert eb.evaluate_candidate("x disease", gold, match(), UIDX) == (True, "eligible")
    assert not eb.evaluate_candidate("Fanconi anemia", gold, match(), UIDX)[0]
    assert not eb.evaluate_candidate("mosaic variegated aneuploidy", gold, match(), UIDX)[0]
    assert not eb.evaluate_candidate("x", gold, match(genes=["A1"]), UIDX)[0]
    assert not eb.evaluate_candidate("x", {"sirolimus"}, match(), UIDX)[0]
    ok, why = eb.evaluate_candidate("x", gold | {"imiglucerase"}, match(), UIDX)
    assert not ok and "IMIGLUCERASE".lower() in why
    assert not eb.evaluate_candidate("x", gold, match(n_matched=99), UIDX)[0]
    assert not eb.evaluate_candidate("x", gold, match(exact_orphanet_mendelian=False), UIDX)[0]
    assert eb.evaluate_candidate("x", gold, match(exact_orphanet_mendelian=False), UIDX, curated=True)[0]


def test_match_disorders_and_symbols():
    hpo = {"OMIM:1": {"name": "Gaucher disease, type 1", "inheritance": ["Autosomal recessive inheritance"], "genes": {"GBA1"}}}
    orph = {"ORPHA:355": {"name": "Gaucher disease", "inheritance": ["Autosomal recessive"], "genes": {"GBA1", "-", "PSAP"}}}
    m = eb.match_disorders("Gaucher disease", hpo, orph)
    assert m["genes"] == ["GBA1", "PSAP"]  # placeholder '-' dropped
    assert m["exact_orphanet_mendelian"] and m["orphanet_catalogued"] and m["mendelian"]
    assert m["n_matched"] == 2


def test_draw_is_deterministic_and_order_free():
    pool = [f"d{i}" for i in range(50)]
    a = eb.draw(pool)
    b = eb.draw(random.Random(1).sample(pool, 50))
    assert a == b and len(a) == eb.N_TARGET == 30
    assert eb.SEED == 20261015
    assert eb.draw(pool[:5]) == sorted(pool[:5])


def test_duplicate_gold_collapsed():
    gold = {"b": {"x", "y"}, "a": {"x", "y"}, "c": {"x", "z"}}
    keep, dropped = eb.drop_duplicate_gold(["a", "b", "c"], gold)
    assert keep == ["a", "c"] and dropped == {"b": "a"}


def study(nct, start, phases=("PHASE2",), typ="INTERVENTIONAL", ivs=None):
    return {"protocolSection": {
        "identificationModule": {"nctId": nct}, "designModule": {"studyType": typ, "phases": list(phases)},
        "statusModule": {"startDateStruct": {"date": start}},
        "armsInterventionsModule": {"interventions": ivs or [{"name": "Sirolimus 2 mg", "type": "DRUG", "otherNames": ["Rapamune"]}]}}}


def test_gold2_rules():
    dre = eb.drug_name_regex(["sirolimus", "ivacaftor", "ASA"])
    sts = [
        study("N1", "2019-01"),                                              # first day after cutoff: counts
        study("N2", "2018-12-31"),                                           # not strictly after
        study("N3", "2020-05-01", phases=("PHASE1",)),                       # phase 1
        study("N4", "2020-05-01", phases=("PHASE2", "PHASE3"),
              ivs=[{"name": "placebo", "type": "DRUG", "otherNames": ["IVACAFTOR"]}]),
        study("N5", "2021", typ="OBSERVATIONAL"),
        study("N6", "2021", ivs=[{"name": "sirolimus device", "type": "DEVICE"}]),
        study("N7", "2021", ivs=[{"name": "Presirolimus", "type": "DRUG"}]),  # whole word only
        study("N8", None),
    ]
    g = eb.gold2_from_studies(sts, dre)
    assert set(g) == {"sirolimus", "ivacaftor"}
    assert g["sirolimus"] == [{"nct": "N1", "start_year": 2019, "phases": ["PHASE2"]}]
    assert g["ivacaftor"][0]["nct"] == "N4"
    assert eb._parse_date("2019-03") == date(2019, 3, 1) and eb._parse_date("x") is None


def test_fetch_studies_paging_and_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(eb, "CACHE", tmp_path)
    calls = []

    def getter(url, params):
        calls.append(params.get("pageToken"))
        return {"studies": [study(f"N{len(calls)}", "2020")], "nextPageToken": "t" if len(calls) == 1 else None}

    out = eb.fetch_studies("some disease", getter)
    again = eb.fetch_studies("some disease", getter)  # served from the cache
    assert again == out
    assert len(out) == 2 and calls == [None, "t"]


def test_build_offline(tmp_path):
    root = tmp_path
    (root / "data").mkdir()
    drugs = [dict(v, name=v["name"]) for v in UIDX.values()]
    (root / "data" / "universe.json").write_text(json.dumps({"drugs": drugs}), encoding="utf-8")
    cache = root / "data" / "cache_engine"
    old = eb.CACHE
    eb.CACHE = cache
    try:
        cache.mkdir(parents=True)
        (cache / "primekg_indications.json").write_text(json.dumps({"toy disease": ["sirolimus", "ivacaftor"]}), encoding="utf-8")
        hpo = {"ORPHA:9": {"name": "Toy disease", "inheritance": ["Autosomal dominant"], "genes": {"AAA1", "BBB2"}}}
        res = eb.build(root, hpo=hpo, orph={}, studies_fn=lambda n: [study("N1", "2020")],
                       annotate=lambda s: {"accession": "P1", "terms": ["DNA repair"]})
    finally:
        eb.CACHE = old
    assert res["n_benchmark"] == 1
    d = res["diseases"][0]
    assert [g["name"] for g in d["gold1"]] == ["ivacaftor", "sirolimus"]
    assert d["gold2"][0]["name"] == "sirolimus" and d["gold2"][0]["in_gold1"]
    assert {g["mechanism_class"] for g in d["genes"]} == {"DNA repair/replication"}
    assert res["cutoff_year_Y"] == 2018
