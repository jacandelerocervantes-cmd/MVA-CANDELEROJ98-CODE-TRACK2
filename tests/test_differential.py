"""Offline tests for track2.differential (synthetic gene-level inputs)."""

from track2 import differential as df

PHENO = {"scores": {"AAA": 0.5, "BBB": 0.1}, "ranks": {"AAA": 10, "BBB": 3000, "NOHPO": None}}


def land(tier="curated", inh=("Autosomal recessive inheritance",), sources=("HPO (x)",)):
    return {"symbol": "X", "tier": tier, "inheritance_as_annotated": list(inh), "sources": list(sources)}


def rk(typ="compound_het", n=1, rarity=9.0):
    return {"type": typ, "n_high_impact": n, "log10_joint_rarity": rarity, "phenotype_score": 0.5}


def ev(gene, lg, rr, st=None, sg=None):
    return df.evaluate_gene(gene, lg, rr, st, sg, PHENO)


def test_causal_candidate():
    r = ev("AAA", land(), rk(), sg={"clinvar": [], "finding_types": ["primary"]})
    assert r["verdict"] == "causal candidate" and r["c2"] == "consistent" and r["c3"] == "high_impact"


def test_causal_via_clinvar_plp_without_high_impact():
    r = ev("AAA", land(), rk(n=0), sg={"clinvar": [{"classification": "Pathogenic/Likely pathogenic"}]})
    assert r["c5"] == "P/LP" and r["verdict"] == "causal candidate"


def test_proposed_gene_never_causal():
    assert ev("AAA", land(tier="proposed, not curated"), rk())["verdict"] == "contextual"


def test_non_mva_gene_insufficient_even_if_biallelic():
    assert ev("AAA", None, rk())["verdict"] == "insufficient evidence"


def test_unpaired_structural_call_is_contextual():
    r = ev("AAA", land(), None, st={"secondary"})
    assert r["c2"] == "unpaired" and r["c3"] == "structural_unpaired" and r["verdict"] == "contextual"


def test_vus_only_biallelic_is_insufficient():
    r = ev("BBB", land(), rk(n=0), sg={"clinvar": [{"classification": "Uncertain significance"}]})
    assert r["c5"] == "uncertain" and r["verdict"] == "insufficient evidence"


def test_no_hpo_annotation_is_insufficient():
    r = ev("NOHPO", land(), None)
    assert r["c4"]["has_hpo"] is False and r["verdict"] == "insufficient evidence"


def test_not_recessive_annotation_blocks_causal():
    r = ev("AAA", land(inh=("Autosomal dominant inheritance",)), rk())
    assert r["c2"] == "biallelic_inheritance_not_recessive" and r["verdict"] == "contextual"


def test_hemizygous_needs_x_linked_annotation():
    assert df.c2_pattern(rk("hemizygous"), ["Autosomal recessive inheritance"]) == "inconsistent"
    assert df.c2_pattern(rk("hemizygous"), ["X-linked recessive inheritance"]) == "consistent"


def test_excluded_by_common_pair_and_by_artefact_type():
    assert ev("AAA", land(), rk(rarity=1.5))["verdict"] == "excluded"
    assert ev("AAA", land(), rk(), st={"artefact"})["verdict"] == "excluded"


def test_loaders_read_only_gene_level_fields(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text("rank,gene,type,n_high_impact,log10_joint_rarity,phenotype_score,variants\n1,G1,homozygous,2,5.0,0.3,SECRET\n")
    out = df.load_ranking(p)
    assert out == {"G1": {"type": "homozygous", "n_high_impact": 2, "log10_joint_rarity": 5.0, "phenotype_score": 0.3}}
    assert "SECRET" not in str(out)
    q = tmp_path / "s.csv"
    q.write_text("chrom_1,pos_1,finding_type,notes\nchr1,123,secondary,G1 (overlaps pos 99)\n")
    assert df.load_structural(q) == {"G1": {"secondary"}}
