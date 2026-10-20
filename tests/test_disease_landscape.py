"""Offline tests for track2.disease_landscape (tiny synthetic annotation files, no network)."""

import json

import pytest

from track2 import disease_landscape as dl

HPOA = (
    "#description: synthetic\n#version: 2099-01-01\n#hpo-version: http://example/hp/2099-01-01\n"
    "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration\n"
    "OMIM:900001\tFoo syndrome 2\t\tHP:0000007\tPMID:1\tPCS\t\t\t\t\tI\tx\n"
    "OMIM:900001\tFoo syndrome 2\t\tHP:0000001\tPMID:1\tPCS\t\t\t\t\tP\tx\n"
    "OMIM:900002\tFoo anaemia, complementation group C\t\tHP:0000006\tPMID:1\tPCS\t\t\t\t\tI\tx\n"
    "OMIM:900003\tUnrelated disease\t\tHP:0000007\tPMID:1\tPCS\t\t\t\t\tI\tx\n"
)
G2P = (
    "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n"
    "1\tAAA1\tHP:0000001\tp\t-\tOMIM:900001\n"
    "2\tBBB2\tHP:0000001\tp\t-\tOMIM:900002\n"
    "3\tCCC3\tHP:0000001\tp\t-\tOMIM:900003\n"
)
OBO = "id: HP:0000007\nname: Autosomal recessive inheritance\n\nid: HP:0000006\nname: Autosomal dominant inheritance\n"
SPEC = dl.DiseaseSpec("foo", "Foo disease", ("foo syndrome", "foo anemia"), abbreviations=("FOO",), out=dl.ROOT / "out/unused.json", track1_crossref=True)


@pytest.fixture
def files(tmp_path):
    (tmp_path / "phenotype.hpoa").write_text(HPOA, encoding="utf-8")
    (tmp_path / "genes_to_phenotype.txt").write_text(G2P, encoding="utf-8")
    (tmp_path / "hp.obo").write_text(OBO, encoding="utf-8")
    return tmp_path


def test_classify_first_rule_wins_and_other():
    assert dl.classify_terms(["centrosome", "DNA repair"])[0] == "DNA repair/replication"
    assert dl.classify_terms(["centrosome", "interstrand cross-link repair"])[0] == "DNA interstrand-crosslink repair"
    assert dl.classify_terms(["kinetochore", "mitotic spindle assembly checkpoint signaling"])[0] == "spindle-assembly checkpoint"
    assert dl.classify_terms(["centriole replication"])[0] == "centrosome/centriole"
    assert dl.classify_terms(["DNA damage"])[0] == "DNA repair/replication"
    assert dl.classify_terms([]) == ("other", [])
    assert dl.classify_terms(["membrane"]) == ("other", [])


def test_rule_table_has_no_disease_or_gene_names():
    text = " ".join(rx.pattern for _, rx in dl.MECHANISM_RULES).lower()
    for word in ("fanconi", "aneuploidy", "mva", "brca", "bub"):
        assert word not in text


def test_name_regex_is_spelling_tolerant():
    rx = dl.SPECS["fanconi_anaemia"].name_re
    assert rx.search("Fanconi anemia, complementation group C") and rx.search("Fanconi-anaemia")
    assert not rx.search("Fanconi-Bickel syndrome")


def test_subtype_parsing():
    assert dl.subtype_of("Fanconi anemia, complementation group D2") == "group D2"
    assert dl.subtype_of("Mosaic variegated aneuploidy syndrome 7 with inflammation") == "type 7"
    assert dl.subtype_of("Fanconi anemia") is None


def test_split_and_verify_sentence():
    text = "First one.  Second\nsentence about FOO syndrome. Third."
    assert "Second sentence about FOO syndrome." in dl.split_sentences(text)
    dl.verify_sentence("Second sentence   about FOO syndrome.", text)
    with pytest.raises(ValueError):
        dl.verify_sentence("Invented sentence.", text)


def test_best_sentence_needs_gene_and_disease_same_sentence():
    text = "AAA1 is a kinase. Foo syndrome is rare. Mutations in AAA1 cause Foo syndrome 2."
    s = dl.best_sentence(text, dl.gene_regex(["AAA1"]), SPEC.assoc_re)
    assert s == "Mutations in AAA1 cause Foo syndrome 2."
    assert dl.best_sentence("AAA1 is a kinase. Foo syndrome is rare.", dl.gene_regex(["AAA1"]), SPEC.assoc_re) is None
    assert dl.best_sentence("AAA1 in the FOO cohort.", dl.gene_regex(["AAA1"]), SPEC.assoc_re)  # case-sensitive abbreviation
    assert dl.best_sentence("AAA1 in the foo cohort.", dl.gene_regex(["AAA1"]), SPEC.assoc_re) is None


def test_hpo_diseases_inheritance_and_genes(files):
    dis = dl.hpo_diseases(SPEC, files / "phenotype.hpoa", files / "hp.obo")
    assert set(dis) == {"OMIM:900001", "OMIM:900002"}
    assert dis["OMIM:900001"]["inheritance"] == ["Autosomal recessive inheritance"]
    assert dl.hpo_gene_map(files / "genes_to_phenotype.txt", set(dis)) == {"OMIM:900001": {"AAA1"}, "OMIM:900002": {"BBB2"}}
    assert dl.hpo_release(files / "phenotype.hpoa")["version"] == "2099-01-01"


def test_build_end_to_end_with_stubs(files):
    orph = [{"symbol": "OLD9", "synonyms": ["AAA1"], "hgnc_id": "HGNC:1", "orpha": "ORPHA:5", "disease_name": "Foo syndrome", "association_type": "t",
             "association_status": "Assessed", "orphanet_pmids": ["1"], "inheritance": "Autosomal recessive"}]
    terms = {"AAA1": ["mitotic spindle assembly checkpoint signaling"], "BBB2": ["membrane"]}
    sentence = "Mutations in BBB2 cause Foo syndrome."

    def evidence(sym, syn, spec):
        if sym == "BBB2":
            dl.verify_sentence(sentence, "Intro. " + sentence)
            return {"evidence": "verified", "pmid": "7", "sentence": sentence}
        return {"evidence": "unverified", "sentence": None}

    res = dl.build(SPEC, files / "phenotype.hpoa", files / "genes_to_phenotype.txt", files / "hp.obo", lambda s: orph,
                   lambda g: {"accession": "P1", "terms": terms[g]}, evidence, lambda: ["AAA1", "ZZZ"], lambda spec: [])
    by = {g["symbol"]: g for g in res["genes"]}
    assert set(by) == {"AAA1", "BBB2"}  # OLD9 merged into AAA1 through the Orphanet synonym; unrelated disease excluded
    assert by["AAA1"]["n_sources"] == 2 and by["BBB2"]["n_sources"] == 1
    assert by["AAA1"]["mechanism_class"] == "spindle-assembly checkpoint" and by["BBB2"]["mechanism_class"] == "other"
    assert by["AAA1"]["in_track1_output"] is True and by["BBB2"]["in_track1_output"] is False
    assert by["AAA1"]["disease_subtype"] == "type 2" and by["BBB2"]["disease_subtype"] == "group C"
    assert res["n_genes_with_verified_sentence"] == 1 and res["header"]["statement"] == dl.HEADER_STATEMENT
    assert res["track1"]["genes_in_track1_output"] == ["AAA1"]
    json.dumps(res)


def test_build_survives_orphanet_failure(files):
    def boom(spec):
        raise RuntimeError("down")

    res = dl.build(SPEC, files / "phenotype.hpoa", files / "genes_to_phenotype.txt", files / "hp.obo", boom,
                   lambda g: {"accession": None, "terms": []}, lambda *a: {"evidence": "unverified"}, lambda: [], lambda spec: [])
    assert res["n_genes"] == 2 and "down" in res["header"]["sources"]["Orphanet"]["error"]


def test_parse_orphanet_synthetic():
    gene = b"""<r><Disorder><OrphaCode>5</OrphaCode><Name lang="en">Foo syndrome</Name><DisorderGeneAssociationList><DisorderGeneAssociation>
    <SourceOfValidation>123[PMID]</SourceOfValidation><Gene><Symbol>AAA1</Symbol><SynonymList><Synonym lang="en">OLD</Synonym></SynonymList>
    <ExternalReferenceList><ExternalReference><Source>HGNC</Source><Reference>77</Reference></ExternalReference></ExternalReferenceList></Gene>
    <DisorderGeneAssociationType><Name lang="en">Disease-causing</Name></DisorderGeneAssociationType></DisorderGeneAssociation></DisorderGeneAssociationList></Disorder>
    <Disorder><OrphaCode>6</OrphaCode><Name lang="en">Other</Name></Disorder></r>"""
    inh = b"""<r><Disorder><OrphaCode>5</OrphaCode><TypeOfInheritanceList><TypeOfInheritance><Name lang="en">Autosomal recessive</Name></TypeOfInheritance></TypeOfInheritanceList></Disorder></r>"""
    (rec,) = dl.parse_orphanet(gene, inh, SPEC)
    assert (rec["symbol"], rec["hgnc_id"], rec["orpha"], rec["inheritance"], rec["orphanet_pmids"]) == ("AAA1", "HGNC:77", "ORPHA:5", "Autosomal recessive", ["123"])


PAPERS = [
    {"pmid": "11", "pubYear": "2023", "title": "Biallelic <i>NEW7</i> variants in Foo syndrome", "abstractText": "Foo syndrome is rare. We report NEW7 variants in two patients with Foo syndrome."},
    {"pmid": "12", "pubYear": "2024", "title": "Homozygous ZZZ9 mutation in a family", "abstractText": "A homozygous ZZZ9 variant was found in one family. The cells showed Foo syndrome features."},
    {"pmid": "13", "pubYear": "2022", "title": "Double homozygosity in AAA1 and DDD4", "abstractText": "Two patients carried DDD4 and AAA1 variants and Foo syndrome."},
    {"pmid": "14", "pubYear": "2022", "title": "MVA and WES in a cohort", "abstractText": "Patients with Foo syndrome were studied by WES."},
    {"pmid": "15", "pubYear": "2021", "title": "Structure of NEW7 domain", "abstractText": "The NEW7 fold is described. Foo syndrome is mentioned."},
    {"pmid": "16", "pubYear": "2021", "title": "More NEW7 patients", "abstractText": "Patients with NEW7 mutations have Foo syndrome."},
]


def test_discover_proposed_filters(monkeypatch):
    ok = lambda t: t in {"NEW7", "ZZZ9", "DDD4"}  # noqa: E731 - stand-in for the UniProt symbol check
    out = dl.discover_proposed(SPEC, {"AAA1"}, {"AAA1"}, PAPERS, ok)
    assert set(out) == {"NEW7", "ZZZ9"}  # DDD4 skipped (curated AAA1 in title); MVA/WES rejected as symbols; "More NEW7" has no gene in title
    assert out["NEW7"]["papers"][0]["sentence_basis"].startswith("gene and disease")
    assert out["ZZZ9"]["papers"][0]["sentence_basis"].startswith("gene and variant")  # fallback sentence, flagged
    for pr in out["NEW7"]["papers"]:
        assert "NEW7" in pr["sentence"]


def test_build_adds_proposed_tier(files):
    ok = lambda t: t == "NEW7"  # noqa: E731
    res = dl.build(SPEC, files / "phenotype.hpoa", files / "genes_to_phenotype.txt", files / "hp.obo", lambda s: [],
                   lambda g: {"accession": None, "terms": []}, lambda *a: {"evidence": "unverified", "sentence": None}, lambda: [],
                   lambda spec: PAPERS[:1], ok)
    by = {g["symbol"]: g for g in res["genes"]}
    assert by["AAA1"]["tier"] == "curated" and by["NEW7"]["tier"] == "proposed, not curated"
    assert by["NEW7"]["n_independent_papers"] == 1 and by["NEW7"]["proposing_pmids"] == ["11"]
    assert by["NEW7"]["association_evidence"]["evidence"] == "verified"
    assert res["n_genes_by_tier"] == {"curated": 2, "proposed, not curated": 1}


def test_title_symbol_candidates_strips_markup():
    assert dl.title_symbol_candidates("Recessive &lt;i&gt;SMC5&lt;/i&gt; Variants in MVA") == ["MVA", "SMC5"]
