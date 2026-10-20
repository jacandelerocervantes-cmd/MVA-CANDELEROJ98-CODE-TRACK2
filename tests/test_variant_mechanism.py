"""Offline tests for track2.variant_mechanism (small synthetic fixtures, no network)."""
import math

import pytest

from track2 import variant_mechanism as vm


def _exons(lengths, strand=1, start=1000):
    """Synthetic exons separated by 100-nt introns, returned in genomic (ascending) order."""
    out, pos = [], start
    for L in lengths:
        out.append({"start": pos, "end": pos + L - 1})
        pos += L + 100
    return out


def test_transcript_order_and_genomic_to_tx_both_strands():
    ex = _exons([50, 60, 70])
    fwd = vm.transcript_order_exons(ex, 1)
    rev = vm.transcript_order_exons(ex, -1)
    assert vm.exon_lengths(fwd) == [50, 60, 70]
    assert vm.exon_lengths(rev) == [70, 60, 50]
    assert vm.genomic_to_tx(fwd, ex[0]["start"], 1) == 1
    assert vm.genomic_to_tx(fwd, ex[1]["start"], 1) == 51
    assert vm.genomic_to_tx(rev, ex[2]["end"], -1) == 1          # first base of the minus-strand transcript
    assert vm.genomic_to_tx(fwd, ex[0]["end"] + 5, 1) is None     # intronic


def test_nmd_predicted_when_far_from_last_junction():
    # exons 300,300,300,100 ; CDS starts at tx 1 ; codon 30 -> nt 88-90 in exon 1
    r = vm.nmd_assessment([300, 300, 300, 100], 1, 30)
    assert r["stop_exon_number"] == 1 and r["n_exons"] == 4
    assert r["distance_to_last_junction_nt"] == 900 - 90
    assert r["nmd_predicted"] is True and r["distance_class"].startswith("NMD_predicted")


def test_nmd_escape_in_last_exon_and_within_50nt():
    r = vm.nmd_assessment([300, 300, 300], 1, 200)                 # nt 598-600 -> exon 2
    assert r["stop_exon_number"] == 2 and r["distance_to_last_junction_nt"] == 0
    assert r["nmd_predicted"] is False
    r = vm.nmd_assessment([300, 300, 300], 1, 150)                 # nt 448-450, exon 2, junction at 600
    assert r["distance_to_last_junction_nt"] == 150 and r["nmd_predicted"] is True
    r = vm.nmd_assessment([300, 300, 300], 1, 190)                 # nt 568-570 -> dist 30
    assert r["distance_to_last_junction_nt"] == 30 and r["nmd_predicted"] is False
    r = vm.nmd_assessment([300, 300, 300], 1, 250)                 # nt 748-750 in last exon
    assert r["stop_in_last_exon"] and r["nmd_predicted"] is False


def test_nmd_borderline_window_and_cds_offset():
    # 5'UTR of 20 nt: CDS nt1 = tx 21; codon 10 -> CDS 28-30 -> tx 48-50; junction at 100 -> distance 50
    r = vm.nmd_assessment([100, 100], 21, 10)
    assert r["distance_to_last_junction_nt"] == 50
    assert r["nmd_predicted"] is None and r["distance_class"].startswith("borderline")
    r = vm.nmd_assessment([100, 100], 21, 9)                        # distance 53 -> still borderline
    assert r["distance_to_last_junction_nt"] == 53 and r["nmd_predicted"] is None
    r = vm.nmd_assessment([100, 100], 21, 8)                        # distance 56 -> predicted
    assert r["distance_to_last_junction_nt"] == 56 and r["nmd_predicted"] is True


def test_lost_features_split():
    feats = [{"type": "Domain", "description": "N", "start": 10, "end": 50, "source": "UniProt"},
             {"type": "Domain", "description": "Protein kinase", "start": 60, "end": 100, "source": "UniProt"},
             {"type": "Region", "description": "straddles", "start": 40, "end": 70, "source": "UniProt"}]
    r = vm.lost_features(feats, 55, 100)
    assert [f["description"] for f in r["retained"]] == ["N"]
    assert [f["description"] for f in r["lost_entirely"]] == ["Protein kinase"]
    assert [f["description"] for f in r["partially_lost"]] == ["straddles"]
    assert r["truncated_product"]["n_lost"] == 46
    assert math.isclose(r["truncated_product"]["fraction_lost"], 0.46)


def test_scan_motif_and_translate():
    assert vm.scan_motif("AKENKENA", "KEN") == [2, 5]
    assert vm.translate("ATGCTAAAATGA") == "MLK*"


def test_window_identity_with_gaps():
    human = "ABC-DEFGH"
    other = "ABCXD-FGH"
    # human residues: A1 B2 C3 D4 E5 F6 G7 H8 ; position 5 (E) aligned to '-'
    r = vm.window_identity(human, other, 5, 2)    # window residues 3..7 = C D E F G
    assert r["aligned_residue"] == "-" and not r["identical_at_position"]
    assert r["window_residues"] == 5 and r["window_identical"] == 4 and r["window_gaps"] == 1
    r = vm.window_identity(human, other, 4, 1)    # D identical
    assert r["aligned_residue"] == "D" and r["identical_at_position"]
    with pytest.raises(ValueError):
        vm.window_identity(human, other, 99, 1)


def test_parse_protein_change_and_classification():
    t = "NM_001211.6(BUB1B):c.742G>T (p.Ala248Ser)"
    assert vm.parse_protein_change(t) == (248, "A", "S")
    assert vm.parse_protein_change("NM_001211.6(BUB1B):c.3086del (p.Ser1029fs)") is None
    assert vm.parse_protein_change("NM_999999.1(OTHER):c.1A>T (p.Ala1Ser)") is None
    assert vm.classify_clinvar("Likely pathogenic") == "P_LP"
    assert vm.classify_clinvar("Benign/Likely benign") == "B_LB"
    assert vm.classify_clinvar("Uncertain significance") == "other"
    assert vm.classify_clinvar("Conflicting classifications of pathogenicity") == "other"


def test_fisher_window_counts_and_p():
    vs = ([{"pos": 1000 + i, "cls": "P_LP"} for i in range(5)] +
          [{"pos": 1010, "cls": "B_LB"}] +
          [{"pos": 100 + i, "cls": "B_LB"} for i in range(20)] +
          [{"pos": 200, "cls": "P_LP"}, {"pos": 300, "cls": "other"}])
    r = vm.fisher_window(vs, 1002, 20)
    assert r["in_window"] == {"pathogenic": 5, "benign": 1, "n": 6}
    assert r["outside_window"] == {"pathogenic": 1, "benign": 20, "n": 21}
    assert r["fisher_p_two_sided"] < 0.01 and r["fisher_odds_ratio"] > 1
    empty = vm.fisher_window([{"pos": 5, "cls": "B_LB"}], 1002, 20)
    assert empty["fisher_p_two_sided"] is None and "fisher_note" in empty


def test_fraction_and_percentile():
    vals = [0.1, 0.2, 0.3, 0.4]
    assert vm.fraction_at_least_as_damaging(vals, 0.3, True)["n_at_least_as_damaging"] == 2
    assert vm.fraction_at_least_as_damaging(vals, 0.3, False)["n_at_least_as_damaging"] == 3
    assert vm.fraction_at_least_as_damaging([], 0.3, True)["n"] == 0
    assert vm.percentile_rank(vals, 0.3) == pytest.approx(62.5)


def test_shrake_rupley_isolated_atom_and_buried_atom():
    p = vm.METHOD_CONSTANTS["sasa_probe_radius_A"]
    s = vm.shrake_rupley([(0.0, 0.0, 0.0)], [1.7], p, 400)
    assert s[0] == pytest.approx(4 * math.pi * (1.7 + p) ** 2, rel=1e-6)
    # an atom enclosed by a shell of overlapping neighbours is fully buried
    shell = [(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1) if (x, y, z) != (0, 0, 0)]
    coords = [(0.0, 0.0, 0.0)] + [(1.2 * x, 1.2 * y, 1.2 * z) for x, y, z in shell]
    s = vm.shrake_rupley(coords, [1.7] * len(coords), p, 200)
    assert s[0] == 0.0


def _pdb_line(serial, name, res, resi, xyz, b, el):
    return (f"ATOM  {serial:5d} {name:<4s} {res:>3s} A{resi:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00{b:6.2f}          {el:>2s}")


def _synthetic_pdb():
    """Ten-residue straight chain along x (3.8 A CA spacing); residue 5 is Asn with CB; 1 Lys elsewhere."""
    lines, k = [], 1
    for i in range(1, 11):
        res = "ASN" if i == 5 else ("LYS" if i == 2 else "ALA")
        x = 3.8 * i
        lines.append(_pdb_line(k, "CA", res, i, (x, 0.0, 0.0), 90.0 if i != 5 else 55.5, "C")); k += 1
        lines.append(_pdb_line(k, "CB", res, i, (x, 1.5, 0.0), 90.0 if i != 5 else 55.5, "C")); k += 1
    return "\n".join(lines)


def test_parse_pdb_and_structure_features():
    pdb = _synthetic_pdb()
    atoms = vm.parse_pdb(pdb)
    assert len(atoms) == 20 and atoms[8]["res"] == "ASN" and atoms[8]["b"] == pytest.approx(55.5)
    known = [{"pos": 7, "label": "p.Ala7Pro", "source": "synthetic"},
             {"pos": 10, "label": "p.Ala10Pro", "source": "synthetic"}]
    f = vm.structure_features(pdb, 5, "Asn", known)
    assert f["plddt"] == 55.5 and f["plddt_confident"] is False
    assert f["secondary_structure_heuristic"].startswith("not assigned")
    d = {x["variant"]: x for x in f["known_pathogenic_missense_distances"]}
    assert d["p.Ala7Pro"]["sequence_distance_residues"] == 2
    assert d["p.Ala7Pro"]["ca_ca_distance_A"] == pytest.approx(7.6, abs=0.01)
    assert d["p.Ala7Pro"]["within_8A_any_atom"] is True
    assert d["p.Ala10Pro"]["within_8A_any_atom"] is False
    assert f["known_pathogenic_missense_distances"][0]["variant"] == "p.Ala7Pro"   # sorted by 3D distance
    assert "Lys2" not in f["residues_within_8A"] and "Ala3" in f["residues_within_8A"]
    with pytest.raises(ValueError):
        vm.structure_features(pdb, 5, "Lys", known)                                  # wrong residue guard


def test_ss_heuristic_helix_and_strand():
    ca = {}
    for i in range(1, 12):  # ideal alpha helix: r=2.3 A, 1.5 A rise, 100 deg per residue
        a = math.radians(100 * i)
        ca[i] = (2.3 * math.cos(a), 2.3 * math.sin(a), 1.5 * i)
    assert vm.ss_heuristic(ca, 5) == "helix-like"
    strand = {i: (3.3 * i, 0.9 * (-1) ** i, 0.0) for i in range(1, 12)}
    assert vm.ss_heuristic(strand, 5) == "strand-like"


def test_verify_sentences_normalises_markup_and_whitespace():
    abstract = "<h4>X</h4>First sentence.  Low protein   abundance is <i>the</i> result. Last."
    r = vm.verify_sentences(abstract, {"ok": "Low protein abundance is the result.",
                                       "bad": "Low protein abundance is not the result."})
    assert r["ok"]["verified_verbatim_in_europepmc_abstract"] is True
    assert r["bad"]["verified_verbatim_in_europepmc_abstract"] is False


def test_uniprot_known_pathogenic_filters_to_missense_mva_pcs():
    feats = [{"type": "Natural variant", "description": "in MVA1; x", "start": 10, "end": 10,
              "alt": {"originalSequence": "L", "alternativeSequences": ["P"]}},
             {"type": "Natural variant", "description": "in PCS", "start": 20, "end": 20,
              "alt": {"originalSequence": "R", "alternativeSequences": ["Q"]}},
             {"type": "Natural variant", "description": "polymorphism", "start": 30, "end": 30,
              "alt": {"originalSequence": "A", "alternativeSequences": ["T"]}},
             {"type": "Domain", "description": "MVA1", "start": 1, "end": 5, "alt": None}]
    r = vm.uniprot_known_pathogenic(feats)
    assert [x["label"] for x in r] == ["p.Leu10Pro", "p.Arg20Gln"]


def test_mechanism_statements_follow_features():
    nmd = vm.nmd_assessment([300, 300, 300, 100], 1, 30)
    lost = vm.lost_features([{"type": "Domain", "description": "Protein kinase", "start": 60, "end": 100,
                              "source": "UniProt"}], 55, 100)
    s1 = vm.nonsense_statement(nmd, lost, True)
    assert "nonsense-mediated decay" in s1["statement"] and "kinase-like domain" in s1["statement"]
    assert all(c["feature"] for c in s1["clauses"]) and s1["what_would_test_it"]
    s1b = vm.nonsense_statement(None, None, False)
    assert "could not be assessed" in s1b["statement"]

    cons = {"all_one2one_orthologs": {"fraction_identical_at_position_of_aligned": 0.9, "n_identical": 9,
                                      "n_aligned": 10, "mean_window_identity": 0.8}}
    st = {"plddt": 92.0, "plddt_confident": True, "sasa_A2": 3.0, "sasa_percentile_among_same_aa": 5.0,
          "n_residues_within_8A": 30, "known_pathogenic_within_8A": ["p.Leu1012Pro"],
          "known_pathogenic_missense_distances": [{"variant": "p.Leu1012Pro", "sequence_distance_residues": 10,
                                                   "ca_ca_distance_A": 9.1}]}
    reg = {"in_window": {"pathogenic": 2, "benign": 0}, "fisher_p_two_sided": 0.2}
    s2 = vm.missense_statement({"conservation": cons, "structure": st, "clinvar_window": reg,
                                "domain": [766, 1050]})
    assert "conserved position" in s2["statement"] and "buried" in s2["statement"]
    assert s2["predicted_stability_effect"].startswith("expected")
    assert "not demonstrable" in s2["statement"]
    st2 = dict(st, plddt=40.0, plddt_confident=False)
    s3 = vm.missense_statement({"conservation": cons, "structure": st2, "clinvar_window": None, "domain": None})
    assert s3["predicted_stability_effect"].startswith("not determinable")
    s4 = vm.missense_statement({"conservation": None, "structure": None, "clinvar_window": None, "domain": None})
    assert "could not be assessed" in s4["statement"]
