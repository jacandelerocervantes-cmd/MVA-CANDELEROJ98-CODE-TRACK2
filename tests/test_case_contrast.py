"""Offline tests for track2.case_contrast and the sentence verification of track2.case_reference."""

from track2 import case_contrast as cc
from track2 import case_profile as cp
from track2 import case_reference as cr


def claim(value, pmid="1"):
    return {"value": value, "pmid": pmid, "population": "synthetic", "sentence": "s"}


def feat(name, hpo, claims):
    return {"feature": name, "hpo_ids": hpo, "claims": claims}


def prof(features=(), history="unknown", genes=()):
    return cp.build_profile("syn", phenotypic_features=list(features), tumour_history=history,
                            genomic_findings=list(genes), when="2026-01-01")


def test_parse_frequency():
    assert cc.parse_frequency("93.5%") == 0.935 and cc.parse_frequency("7/8") == 0.875
    assert cc.parse_frequency(None) is None and cc.parse_frequency("many") is None


def test_rule_present_and_reported():
    p = prof([{"hpo_id": "HP:0000252", "status": "present"}])
    f = feat("a", ["HP:0000252"], [claim("87%", "35804254")])
    assert cc.child_status(p, f) == "present"
    assert cc.read_row("present", f["claims"]) == "present in the child and reported in 87% of published cases, PMID 35804254"


def test_rule_unknown_majority_means_check_clinically():
    r = cc.read_row("unknown", [claim("77.4%")])
    assert "not recorded in the child although reported in a majority of cases" in r and "check clinically" in r


def test_rule_unknown_minority_has_no_check_instruction():
    r = cc.read_row("unknown", [claim("25.8%")])
    assert "minority" in r and "check clinically" not in r


def test_rule_present_not_reported():
    assert cc.read_row("present", []) == "present in the child but not reported in published cases"
    assert cc.read_row("unknown", []) == "not recorded in the child; not reported in the sources"


def test_rule_named_without_frequency():
    assert "named without a frequency" in cc.read_row("present", [claim(None)])


def test_absent_status_requires_all_terms_absent():
    f = feat("a", ["HP:0000252", "HP:0001511"], [])
    only_one = prof([{"hpo_id": "HP:0000252", "status": "absent"}])
    both = prof([{"hpo_id": "HP:0000252", "status": "absent"}, {"hpo_id": "HP:0001511", "status": "absent"}])
    assert cc.child_status(only_one, f) == "unknown" and cc.child_status(both, f) == "absent"


def test_any_malignancy_uses_tumour_history():
    f = feat("any malignancy", ["HP:0002664"], [claim("38.7%")])
    assert cc.child_status(prof(history="yes"), f) == "present"
    assert cc.child_status(prof(history="no"), f) == "absent"
    assert cc.child_status(prof(), f) == "unknown"


def gene(classes, pattern="compound heterozygous", phase="unresolved"):
    return {"gene_symbol": "BUB1B", "mechanism_class": "spindle-assembly checkpoint", "consequence_classes": classes,
            "genotype_pattern": pattern, "phase_status": phase}


REF = {"disease": "d", "retrieved": "2026-01-01", "features": [feat("microcephaly", ["HP:0000252"], [claim("87%")])],
       "genotype_claims": [{"pattern": "p", "pmid": "20516114", "sentence": "s"}],
       "genotype_not_reported": {"two truncating variants (frequency)": cc.NOT_REPORTED}}


def test_genotype_truncating_plus_missense_matches_published_pattern():
    g = cc.genotype_contrast(prof(genes=[gene(["loss of function", "missense"])]), REF)
    assert "pattern class matches" in g["reading"] and "phase not resolved" in g["reading"]
    assert g["not_reported"]["two truncating variants (frequency)"] == cc.NOT_REPORTED


def test_genotype_two_missense_and_missing_gene():
    g = cc.genotype_contrast(prof(genes=[gene(["missense"])]), REF)
    assert "pattern class matches" not in g["reading"]
    assert "no BUB1B finding" in cc.genotype_contrast(prof(), REF)["reading"]


def test_build_contrast_table_shape():
    out = cc.build_contrast(prof([{"hpo_id": "HP:0000252", "status": "absent"}]), REF)
    row = out["rows"][0]
    assert row["child"] == "absent" and row["feature"] == "microcephaly" and row["reading"].startswith("absent in the child")
    assert out["genotype"]["gene"] == "BUB1B"


TEXT = "<p>Growth retardation   occurred in 93.5%, and microcephaly in 87%&#160;of cases.</p>"


def test_sentence_verification_is_verbatim_and_value_bound():
    assert cr.verify_claim(TEXT, "Growth retardation occurred in 93.5%, and microcephaly in 87% of cases.", "93.5%")
    assert not cr.verify_claim(TEXT, "Growth retardation occurred in 93.5%, and microcephaly in 87% of cases.", "50%")
    assert not cr.verify_claim(TEXT, "Growth retardation occurred in 90%.", None)


def test_reference_build_fails_on_unverified_sentence():
    texts = {p: {"text": "nothing here", "kind": "abstract"} for p in cr.SOURCES}
    try:
        cr.build_reference(texts)
    except ValueError as exc:
        assert "not verified" in str(exc)
    else:
        raise AssertionError("expected ValueError")
