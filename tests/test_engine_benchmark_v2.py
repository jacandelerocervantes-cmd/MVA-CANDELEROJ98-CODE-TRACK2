"""Offline tests of benchmark v2 (synthetic data, injected getters, no network)."""
import json
import re

from track2 import engine_benchmark_v2 as b2

UIDX = {
    "alphadrug": {"drug_id": "C1", "name": "ALPHADRUG", "molecule_type": "Small molecule", "atc": ["A10BA02"]},
    "betadrug": {"drug_id": "C2", "name": "BETADRUG", "molecule_type": "Small molecule", "atc": ["N03AX09"]},
    "gammadrug": {"drug_id": "C3", "name": "GAMMADRUG", "molecule_type": "Small molecule", "atc": []},
    "zetaenzyme": {"drug_id": "C4", "name": "ZETAENZYME", "molecule_type": "Enzyme", "atc": []},
}
DRE = b2.v1.drug_name_regex(sorted(UIDX))

OBO = """format-version: 1.2
data-version: releases/2099-01-01

[Term]
id: MONDO:0000010
name: testosis syndrome
synonym: "tests disorder" EXACT []
synonym: "TS" RELATED ABBREVIATION []
synonym: "typed one" EXACT AMBIGUOUS []
synonym: "related name" RELATED []
xref: OMIM:100100 {source="MONDO:equivalentTo"}
xref: OMIM:100200 {source="MONDO:relatedTo"}
xref: OMIMPS:200200 {source="MONDO:equivalentTo"}
xref: Orphanet:777 {source="MONDO:equivalentTo", source="DOID:1"}
xref: Orphanet:888 {source="MONDO:obsoleteEquivalent"}

[Term]
id: MONDO:0000011
name: old thing
is_obsolete: true
xref: OMIM:555 {source="MONDO:equivalentTo"}

[Typedef]
id: part_of
"""


def disorder(name, genes=(), inh=()):
    return {"name": name, "inheritance": list(inh), "genes": set(genes)}


def test_load_mondo_keeps_equivalent_xrefs_exact_untyped_synonyms_and_drops_obsolete(tmp_path):
    f = tmp_path / "m.obo"
    f.write_text(OBO, encoding="utf-8")
    terms, version = b2.load_mondo(f)
    assert version == "releases/2099-01-01"
    assert set(terms) == {"MONDO:0000010"}
    t = terms["MONDO:0000010"]
    assert t["omim"] == {"100100"} and t["omimps"] == {"200200"} and t["orpha"] == {"777"}
    assert t["synonyms"] == ["tests disorder"]


def test_mondo_id_padding_and_primekg_ids(tmp_path):
    assert b2.mondo_id("9067") == "MONDO:0009067"
    f = tmp_path / "nodes.tab"
    f.write_text("node_index,node_id,node_type,node_name,node_source\n0,9796,gene/protein,X,NCBI\n"
                 "1,9067_19745,disease,cystic thing,MONDO_grouped\n2,11669,disease,other thing,MONDO\n", encoding="utf-8")
    ids = b2.load_primekg_disease_ids(f)
    assert ids == {"cystic thing": ["MONDO:0009067", "MONDO:0019745"], "other thing": ["MONDO:0011669"]}


def test_name_set_drops_short_obsolete_parenthetical_and_duplicates():
    ns = b2.build_name_set(["Testosis (disease)", "testosis", "TS", "obsolete testosis", "Test-osis syndrome", "abc"])
    assert ns == ["Test-osis syndrome", "Testosis"] or ns == ["Testosis", "Test-osis syndrome"]
    assert "TS" not in ns and "abc" not in ns


def test_parse_orphanet_names():
    xml = (b"<r><DisorderList><Disorder id='1'><OrphaCode>777</OrphaCode><Name lang='en'>Testosis</Name>"
           b"<SynonymList count='1'><Synonym lang='en'>Syn name</Synonym></SynonymList></Disorder></DisorderList></r>")
    assert b2.parse_orphanet_names(xml) == {"ORPHA:777": {"name": "Testosis", "synonyms": ["Syn name"]}}


def make_index():
    hpo = {"OMIM:100100": disorder("Testosis type 1", ["G1", "G2"], ["Autosomal recessive inheritance"]),
           "OMIM:999": disorder("Unrelated disease", ["Z9"], [])}
    orph = {"ORPHA:777": disorder("Tests disorder", ["G2", "G3"], ["Autosomal recessive"]),
            "ORPHA:1": disorder("Other", ["Q1"], [])}
    return b2.DisorderIndex(hpo, orph)


def make_mondo():
    return {"MONDO:0000010": {"id": "MONDO:0000010", "name": "testosis syndrome", "omim": {"100100"}, "omimps": set(),
                              "orpha": {"777"}, "synonyms": ["tests disorder"], "obsolete": False}}


def test_match_by_identifier_and_synonym_finds_what_the_name_alone_misses():
    index, mondo = make_index(), make_mondo()
    # PrimeKG name "tests disease" matches nothing by name; the MONDO id maps it to OMIM:100100 and ORPHA:777
    name = "tests disease"
    assert not [i for i, n in index.items if b2.v1.term_regex(name).search(n)]
    m = b2.match_disease_v2(name, ["MONDO:0000010"], mondo, index, {})
    assert m["id_matched"] == ["OMIM:100100", "ORPHA:777"]
    assert m["hpo"] == ["OMIM:100100"] and m["orphanet"] == ["ORPHA:777"]
    assert m["genes"] == ["G1", "G2", "G3"] and m["rare_genetic_identity"] and m["mendelian"]
    assert m["n_matched"] == 2 and m["orphanet_catalogued"]


def test_match_by_orphanet_synonym_of_id_matched_disorder():
    index = make_index()
    mondo = {"MONDO:0000010": dict(make_mondo()["MONDO:0000010"], synonyms=[], omim=set())}
    m = b2.match_disease_v2("zzz", ["MONDO:0000010"], mondo, index,
                            {"ORPHA:777": {"name": "Tests disorder", "synonyms": ["Syndrome zzz2"]}})
    assert "Syndrome zzz2" in m["name_set"] and m["id_matched"] == ["ORPHA:777"]


def test_unmapped_disease_without_name_match_has_no_matched_disorders():
    m = b2.match_disease_v2("no such thing", [], make_mondo(), make_index(), {})
    assert m["n_matched"] == 0 and not m["rare_genetic_identity"] and m["genes"] == []


def test_rare_genetic_requires_catalogued_and_mendelian():
    index = b2.DisorderIndex({"OMIM:5": disorder("Foo disease", ["A", "B"], ["Autosomal dominant"])}, {})
    mondo = {"MONDO:1": {"id": "MONDO:1", "name": "foo disease", "omim": {"5"}, "omimps": set(), "orpha": set(),
                         "synonyms": [], "obsolete": False}}
    m = b2.match_disease_v2("foo disease", ["MONDO:1"], mondo, index, {})
    assert not m["rare_genetic_identity"]  # Mendelian but not Orphanet-catalogued (no ORPHA id, no Orphanet xref)
    mondo["MONDO:1"]["orpha"] = {"42"}  # MONDO equivalence to an Orphanet disorder = catalogued
    m = b2.match_disease_v2("foo disease", ["MONDO:1"], mondo, index, {})
    assert m["rare_genetic_identity"]


def test_evaluate_candidate_v2_order_and_reasons():
    index, mondo = make_index(), make_mondo()
    m = b2.match_disease_v2("tests disease", ["MONDO:0000010"], mondo, index, {})
    ok, why = b2.evaluate_candidate_v2("tests disease", {"alphadrug", "betadrug"}, m, UIDX)
    assert ok and why == "eligible"
    assert b2.evaluate_candidate_v2("Fanconi anemia", {"alphadrug", "betadrug"}, m, UIDX)[1].startswith("application case")
    assert "no HPO or Orphanet" in b2.evaluate_candidate_v2("x", {"alphadrug", "betadrug"}, None, UIDX)[1]
    assert "fewer than 2 approved" in b2.evaluate_candidate_v2("x", {"alphadrug"}, m, UIDX)[1]
    assert "enzyme replacement" in b2.evaluate_candidate_v2("x", {"alphadrug", "zetaenzyme"}, m, UIDX)[1]
    assert "umbrella" in b2.evaluate_candidate_v2("x", {"alphadrug", "betadrug"}, dict(m, n_matched=31), UIDX)[1]
    assert "fewer than 2 genes" in b2.evaluate_candidate_v2("x", {"alphadrug", "betadrug"}, dict(m, genes=["A"]), UIDX)[1]
    assert "not a rare genetic" in b2.evaluate_candidate_v2(
        "x", {"alphadrug", "betadrug"}, dict(m, rare_genetic_identity=False), UIDX)[1]
    # curated candidates keep the looser v1 alternative
    assert b2.evaluate_candidate_v2("x", {"alphadrug", "betadrug"}, dict(m, rare_genetic_identity=False), UIDX, True)[0]


def test_draw_new_is_all_when_few_and_seeded_when_many():
    names = [f"d{i:02d}" for i in range(50)]
    assert b2.draw_new(names[:20], 17) == names[:20]
    a = b2.draw_new(names, 17)
    assert len(a) == 43 and a == b2.draw_new(list(reversed(names)), 17) and a == sorted(a)
    assert a != b2.draw_new(names, 17, seed=1)


def test_drop_duplicate_gold_v2_prefers_v1_disease():
    gold = {"new a": {"x", "y"}, "new b": {"x", "y"}, "new c": {"p", "q"}, "new d": {"p", "q"}}
    keep, dropped = b2.drop_duplicate_gold_v2(["new d", "new c", "new b", "new a"], gold, {"old": frozenset({"x", "y"})})
    assert keep == ["new c"] and dropped == {"new a": "old", "new b": "old", "new d": "new c"}


def test_reduce_openfda_and_gold3_whole_phrase_match():
    raw = {"meta": {"results": {"total": 2}},
           "results": [{"set_id": "S1", "effective_time": "20260101", "openfda": {"generic_name": ["ALPHADRUG"]},
                        "indications_and_usage": ["1 INDICATIONS: for Testosis Syndrome in adults."], "dosage": "x"},
                       {"set_id": "S2", "effective_time": "20250101",
                        "indications_and_usage": ["for pretestosis syndrome and sometests"]}]}
    red = b2.reduce_openfda(raw)
    assert red["total"] == 2 and red["labels"][0]["indications"].startswith("1 INDICATIONS") and "dosage" not in red["labels"][0]
    by_drug = {"alphadrug": red["labels"][:1], "betadrug": red["labels"][1:], "gammadrug": []}
    g3 = b2.gold3_for_disease(["testosis syndrome"], by_drug)
    assert list(g3) == ["alphadrug"]  # "pretestosis syndrome" is not a whole-phrase match
    assert g3["alphadrug"][0]["set_id"] == "S1" and g3["alphadrug"][0]["matched_phrase"] == "testosis syndrome"
    assert b2.gold3_for_disease(["tests disorder", "testosis  syndrome"], by_drug).keys() == {"alphadrug"}  # synonym list, spacing
    assert b2.gold3_for_disease([], by_drug) == {}


def test_fetch_labels_uses_cache_and_getter(tmp_path):
    calls = []

    def getter(url, params):
        calls.append(params["search"])
        return {"total": 1, "labels": [{"set_id": "S", "effective_time": "1", "generic_name": [], "indications": "x"}]}

    a = b2.fetch_labels("alphadrug", getter, tmp_path)
    b = b2.fetch_labels("alphadrug", getter, tmp_path)
    assert a == b and len(calls) == 1 and 'openfda.generic_name:"alphadrug"' in calls[0]


def study(nct, arms, interventions, start="2021-05", phases=("PHASE2",), stype="INTERVENTIONAL"):
    return {"protocolSection": {
        "identificationModule": {"nctId": nct}, "statusModule": {"startDateStruct": {"date": start}},
        "designModule": {"studyType": stype, "phases": list(phases)},
        "armsInterventionsModule": {"armGroups": arms, "interventions": interventions}}}


def test_gold2_counts_only_experimental_arms_and_records_arm_and_nct():
    s1 = study("NCT1", [{"label": "Active", "type": "EXPERIMENTAL", "interventionNames": ["Drug: Alphadrug"]},
                        {"label": "Ctrl", "type": "ACTIVE_COMPARATOR", "interventionNames": ["Drug: Betadrug"]},
                        {"label": "Pbo", "type": "PLACEBO_COMPARATOR", "interventionNames": ["Drug: Placebo"]}],
               [{"type": "DRUG", "name": "Alphadrug"}, {"type": "DRUG", "name": "Betadrug"}, {"type": "DRUG", "name": "Placebo"}])
    s2 = study("NCT2", [{"label": "B", "type": "SHAM_COMPARATOR", "interventionNames": ["Drug: Gammadrug"]}],
               [{"type": "DRUG", "name": "Gammadrug"}])
    s3 = study("NCT3", [{"label": "Arm", "type": "EXPERIMENTAL", "interventionNames": ["Drug: BRAND-1"]}],
               [{"type": "DRUG", "name": "BRAND-1", "otherNames": ["Betadrug"]}])
    r = b2.gold2_from_studies_v2([s1, s2, s3], DRE)
    assert set(r["gold"]) == {"alphadrug", "betadrug"}
    assert r["gold"]["alphadrug"] == [{"nct": "NCT1", "start_year": 2021, "phases": ["PHASE2"], "arm_type": "EXPERIMENTAL",
                                       "arm_label": "Active"}]
    assert [t["nct"] for t in r["gold"]["betadrug"]] == ["NCT3"]  # comparator use in NCT1 is not counted
    assert r["excluded"] == {"gammadrug": {"arm_types": ["SHAM_COMPARATOR"], "ncts": ["NCT2"]}}
    assert r["n_qualifying_studies"] == 3


def test_gold2_drug_only_in_comparators_is_excluded_but_kept_if_also_experimental_elsewhere():
    s1 = study("NCT1", [{"label": "C", "type": "ACTIVE_COMPARATOR", "interventionNames": ["Drug: Alphadrug"]}],
               [{"type": "DRUG", "name": "Alphadrug"}])
    s2 = study("NCT2", [{"label": "E", "type": "EXPERIMENTAL", "interventionNames": ["Drug: Alphadrug"]}],
               [{"type": "DRUG", "name": "Alphadrug"}], start="2019-01-01")
    r = b2.gold2_from_studies_v2([s1], DRE)
    assert r["gold"] == {} and r["excluded"]["alphadrug"]["arm_types"] == ["ACTIVE_COMPARATOR"]
    r = b2.gold2_from_studies_v2([s1, s2], DRE)
    assert list(r["gold"]) == ["alphadrug"] and r["excluded"] == {}


def test_gold2_study_filters_and_missing_arm_groups_and_non_drug_types():
    exp = [{"label": "E", "type": "EXPERIMENTAL", "interventionNames": ["Drug: Alphadrug"]}]
    iv = [{"type": "DRUG", "name": "Alphadrug"}]
    studies = [
        study("NCT_OLD", exp, iv, start="2018-12-31"),  # not strictly after the cutoff
        study("NCT_P1", exp, iv, phases=("PHASE1",)),
        study("NCT_OBS", exp, iv, stype="OBSERVATIONAL"),
        study("NCT_NOARMS", [], iv),
        study("NCT_DEV", [{"label": "E", "type": "EXPERIMENTAL", "interventionNames": ["Device: Alphadrug"]}],
              [{"type": "DEVICE", "name": "Alphadrug"}]),
        study("NCT_OK", exp, iv, start="2019", phases=("PHASE2", "PHASE3")),
    ]
    r = b2.gold2_from_studies_v2(studies, DRE)
    assert [t["nct"] for t in r["gold"]["alphadrug"]] == ["NCT_OK"] and r["gold"]["alphadrug"][0]["start_year"] == 2019
    assert r["n_qualifying_without_arms"] == 1 and r["n_qualifying_studies"] == 3
    assert r["excluded"] == {}  # alphadrug is in gold, so it is not reported as excluded


def test_fetch_studies_v2_pages_and_requests_arm_fields(tmp_path, monkeypatch):
    monkeypatch.setattr(b2, "CACHE", tmp_path)
    monkeypatch.setattr(b2.v1, "CACHE", tmp_path)
    pages = {None: {"studies": [{"a": 1}], "nextPageToken": "T"}, "T": {"studies": [{"a": 2}]}}
    seen = []

    def getter(url, params):
        seen.append(params)
        return pages[params.get("pageToken")]

    out = b2.fetch_studies_v2("some disease", getter)
    assert out == [{"a": 1}, {"a": 2}] and "ArmGroupType" in seen[0]["fields"] and seen[0]["query.cond"] == "some disease"


def test_gene_record_has_both_class_fields():
    ann = {"accession": "P0", "terms": ["DNA repair", "double-strand break repair", "Reference proteome"]}
    g = b2.gene_record("SYM1", lambda s: ann)
    assert g["mechanism_class_general"] == "dna_repair_replication" and "mechanism_class" in g and g["symbol"] == "SYM1"
    bad = b2.gene_record("SYM2", lambda s: (_ for _ in ()).throw(RuntimeError("down")))
    assert bad["mechanism_class_general"] == "other" and "annotation_error" in bad


def test_build_end_to_end_synthetic(tmp_path):
    gold1 = {"old disease": {"alphadrug", "betadrug"}, "new disease": {"alphadrug", "gammadrug"},
             "dup disease": {"alphadrug", "betadrug"}, "unmapped": {"alphadrug", "gammadrug"}, "one drug": {"alphadrug"}}
    v1b = {"diseases": [{"name": "old disease", "genes": [{"symbol": "G1"}], "inheritance": ["Autosomal recessive"],
                         "gold1": [{"name": "alphadrug"}, {"name": "betadrug"}]}],
           "selection_log": [{"disease": "new disease", "reason": "no HPO or Orphanet disorder with this name"}]}
    hpo = {"OMIM:100100": disorder("Testosis type 1", ["G1", "G2"], ["Autosomal recessive inheritance"]),
           "OMIM:1": disorder("Old disease", ["G1", "G2"], ["Autosomal recessive inheritance"]),
           "ORPHA:5": disorder("Old disease orpha", ["G3"], ["Autosomal recessive inheritance"])}
    mondo = {"MONDO:0000010": {"id": "MONDO:0000010", "name": "testosis syndrome", "omim": {"100100"}, "omimps": set(),
                               "orpha": {"5"}, "synonyms": [], "obsolete": False},
             "MONDO:0000020": {"id": "MONDO:0000020", "name": "old disease", "omim": {"1"}, "omimps": set(),
                               "orpha": {"5"}, "synonyms": [], "obsolete": False}}
    pk = {"new disease": ["MONDO:0000010"], "old disease": ["MONDO:0000020"], "dup disease": ["MONDO:0000010"]}
    labels = {"gammadrug": [{"set_id": "S", "effective_time": "2026", "generic_name": [],
                             "indications": "indicated for testosis syndrome"}], "alphadrug": [], "betadrug": []}
    exp = [{"label": "E", "type": "EXPERIMENTAL", "interventionNames": ["Drug: Betadrug"]}]
    ctrl = [{"label": "C", "type": "ACTIVE_COMPARATOR", "interventionNames": ["Drug: Alphadrug"]}]
    st = [study("NCT9", exp + ctrl, [{"type": "DRUG", "name": "Betadrug"}, {"type": "DRUG", "name": "Alphadrug"}])]
    out = b2.build(root=tmp_path, hpo=hpo, orph={}, orph_names={}, mondo=mondo, mondo_version="v", primekg_ids=pk,
                   universe=list(UIDX.values()), gold1=gold1, v1_benchmark=v1b, curated=set(), labels_by_drug=labels,
                   studies_fn=lambda n: st,
                   annotate=lambda s: {"accession": "P", "terms": ["Mitosis", "kinetochore"]})
    names = [d["name"] for d in out["diseases"]]
    assert names == ["old disease", "new disease"]  # v1 first; "dup disease" has the gold 1 set of the v1 disease
    assert out["n_v1_retained"] == 1 and out["n_v2_new"] == 1 and out["target_met"] is False
    new = out["diseases"][1]
    assert new["origin"] == "v2_new" and new["genes"][0]["mechanism_class_general"] == "cell_cycle_mitotic_checkpoint"
    assert [g["name"] for g in new["gold1"]] == ["alphadrug", "gammadrug"]
    assert [g["name"] for g in new["gold3"]] == ["gammadrug"] and new["gold3"][0]["in_gold1"]
    assert new["gold_sizes"] == {"gold1": 2, "gold3": 1, "union": 2, "intersection": 1, "gold2": 1, "gold2_excluded": 1}
    assert [g["name"] for g in new["gold2"]] == ["betadrug"] and new["gold2"][0]["trials"][0]["arm_type"] == "EXPERIMENTAL"
    assert list(new["gold2_excluded_non_experimental"]) == ["alphadrug"]
    reasons = {e["disease"]: e["reason"] for e in out["selection_log"]}
    assert reasons["dup disease"].startswith("duplicate gold 1 set of old disease")
    assert reasons["old disease"] == "retained from v1"
    assert out["selection_reason_counts"]["retained from v1"] == 1
    assert out["mechanism_class_general_distribution"]["unique_genes"]["cell_cycle_mitotic_checkpoint"] >= 1
    json.dumps(out, sort_keys=True)  # serialisable with sorted keys
    assert re.fullmatch(r"v2", out["version"])
