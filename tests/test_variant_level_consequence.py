"""Offline tests for variant_level_consequence (synthetic text only; no patient data)."""

import csv
import json

from track2 import variant_level_consequence as v


def test_consequence_class_terms_only():
    assert v.consequence_class("x stop_gained AF=1") == v.LOF
    assert v.consequence_class("x missense_variant") == v.MISSENSE
    assert v.consequence_class("x inframe_insertion") == v.INFRAME
    assert v.consequence_class("x frameshift_variant") == v.LOF
    assert v.consequence_class("x splice_donor_variant") == v.LOF
    assert v.consequence_class("no term here") is None


def test_decay_annotation_states():
    assert v.decay_annotation("stop gained") == "not_annotated"
    assert v.decay_annotation("stop gained, predicted NMD") == "predicted"
    assert v.decay_annotation("transcript escapes NMD") == "escape"


def test_structural_class():
    assert v.structural_class("G DEL 5bp. Overlaps MANE Select coding sequence: no.") == v.SV_NONCODING
    assert v.structural_class("G DEL 5bp. Overlaps MANE Select coding sequence: yes.") == v.SV_CODING
    assert v.structural_class("G DEL 5bp.") == v.SV_UNKNOWN


def test_readers_and_build(tmp_path):
    rk = tmp_path / "r.csv"
    with open(rk, "w", encoding="utf8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "variants"])
        w.writerow(["GA", "p1 stop_gained | p2 missense_variant"])
        w.writerow(["GB", "p3 missense_variant"])
    sv = tmp_path / "s.csv"
    with open(sv, "w", encoding="utf8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["finding_type", "notes"])
        w.writerow(["secondary", "GC single_hit DEL 9bp. Overlaps MANE Select coding sequence: no."])
        w.writerow(["secondary", "GC two heterozygous SV calls in the same gene (DEL + INS)."])
        w.writerow(["secondary", "GA stop_gained component listed separately."])
    snv = v.read_ranking_classes(rk)
    svc, decay = v.read_structural(sv)
    doc = v.build(snv, svc, decay)
    ga, gb, gc = (doc["genes"][g] for g in ("GA", "GB", "GC"))
    assert ga["class_counts"] == {v.LOF: 1, v.MISSENSE: 1}
    assert ga["decay_annotation"] == "not_annotated"
    assert gb["decay_annotation"] == "not_applicable"
    assert gc["class_counts"] == {v.SV_NONCODING: 1}
    # no raw text, only categorical keys and counts
    assert set(ga) == {"class_counts", "n_variants_or_calls", "n_loss_of_function", "decay_annotation"}
    assert "stop_gained" not in json.dumps(doc)


def test_c3_from_classes():
    lof = {"class_counts": {v.LOF: 1, v.MISSENSE: 1}, "n_variants_or_calls": 2, "n_loss_of_function": 1,
           "decay_annotation": "not_annotated"}
    assert v.c3_from_classes(lof) == {"c3": "high_impact", "decay_status": "not_annotated"}
    mis = {"class_counts": {v.MISSENSE: 2}, "n_variants_or_calls": 2, "n_loss_of_function": 0}
    assert v.c3_from_classes(mis)["c3"] == "non_high_impact"
    nc = {"class_counts": {v.SV_NONCODING: 1}, "n_variants_or_calls": 1, "n_loss_of_function": 0}
    assert v.c3_from_classes(nc)["c3"] == "structural_noncoding"
    cod = {"class_counts": {v.SV_CODING: 1}, "n_variants_or_calls": 1, "n_loss_of_function": 0}
    assert v.c3_from_classes(cod)["c3"] == "structural_unpaired"
    assert v.c3_from_classes(None)["c3"] == "none"
