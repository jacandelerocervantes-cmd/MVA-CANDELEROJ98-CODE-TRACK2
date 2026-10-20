from track2.hpo import ancestors, information_content, parse_gene_annotations, parse_obo, propagate
from track2.phenotype_module import MIN_GENES, TOP_N, build_phenotype_module

OBO = """format-version: 1.2

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000100
name: Growth abnormality
is_a: HP:0000001 ! All

[Term]
id: HP:0000200
name: Short stature
alt_id: HP:0009999
is_a: HP:0000100 ! Growth abnormality

[Term]
id: HP:0000201
name: Severe short stature
is_a: HP:0000200 ! Short stature

[Term]
id: HP:0000300
name: Muscle atrophy
is_a: HP:0000001 ! All

[Term]
id: HP:0000400
name: Old term
is_obsolete: true

[Typedef]
id: part_of
"""

G2P = "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n" + "".join(
    f"{i}\t{g}\t{t}\tx\t-\tOMIM:1\n" for i, (g, t) in enumerate(
        [("GENEA", "HP:0000201"), ("GENEA", "HP:0000300"), ("GENEB", "HP:0000200"), ("GENEC", "HP:0000300"), ("GENED", "HP:0009999")]))


def _load(tmp_path):
    (tmp_path / "hp.obo").write_text(OBO, encoding="utf-8")
    (tmp_path / "g.txt").write_text(G2P, encoding="utf-8")
    parents, alt = parse_obo(tmp_path / "hp.obo")
    return parents, alt, parse_gene_annotations(tmp_path / "g.txt", alt)


def test_obo_parsing_skips_obsolete_and_maps_alt_ids(tmp_path):
    parents, alt, ann = _load(tmp_path)
    assert "HP:0000400" not in parents and alt["HP:0009999"] == "HP:0000200"
    assert ancestors("HP:0000201", parents) == {"HP:0000201", "HP:0000200", "HP:0000100", "HP:0000001"}
    assert ann["GENED"] == {"HP:0000200"}  # alternative id replaced


def test_rows_without_a_gene_symbol_are_skipped(tmp_path):
    (tmp_path / "g.txt").write_text("ncbi_gene_id\tgene_symbol\thpo_id\n1\t-\tHP:0000300\n2\tGENEZ\tHP:0000300\n3\t\tHP:0000300\n", encoding="utf-8")
    assert parse_gene_annotations(tmp_path / "g.txt", {}) == {"GENEZ": {"HP:0000300"}}


def test_descendant_annotation_counts_for_the_parent_term(tmp_path):
    parents, alt, ann = _load(tmp_path)
    module = build_phenotype_module(["HP:0000200", "HP:0000300"], ann, parents, alt)
    # only GENEA is annotated to both terms (HP:0000201 is a child of HP:0000200); fewer than MIN_GENES -> fallback rule
    assert module["rule"].startswith("top_") and "GENEA" in module["genes"]
    assert module["gene_detail"]["GENEA"]["n_terms"] == 2
    assert module["gene_detail"]["GENEA"]["score"] > module["gene_detail"]["GENEC"]["score"]  # two terms outweigh one


def test_strict_rule_when_enough_genes(tmp_path):
    parents, alt, _ = _load(tmp_path)
    ann = {f"G{i}": {"HP:0000201", "HP:0000300"} for i in range(MIN_GENES + 2)}
    ann["ONLYONE"] = {"HP:0000300"}
    module = build_phenotype_module(["HP:0000200", "HP:0000300"], ann, parents, alt)
    assert module["rule"] == "annotated_to_at_least_2_terms"
    assert "ONLYONE" not in module["genes"] and len(module["genes"]) == MIN_GENES + 2


def test_fallback_is_capped_and_ranked_by_information_content():
    parents = {"T1": set(), "T2": set(), "RARE": {"T1"}}
    ann = {f"G{i}": {"T1"} for i in range(300)}
    ann["R"] = {"RARE"}
    module = build_phenotype_module(["T1", "RARE"], ann, parents)
    assert len(module["genes"]) <= TOP_N and "R" in module["genes"]  # the gene with the rarer term outranks the rest


def test_terms_without_annotation_are_reported(tmp_path):
    parents, alt, ann = _load(tmp_path)
    parents["HP:0777777"] = set()
    module = build_phenotype_module(["HP:0777777", "HP:0000300"], ann, parents, alt)
    assert module["terms_without_annotation"] == ["HP:0777777"]


def test_information_content_orders_rare_above_common(tmp_path):
    parents, _, ann = _load(tmp_path)
    ic = information_content(propagate(ann, parents))
    assert ic["HP:0000201"] > ic["HP:0000200"] >= ic["HP:0000001"] == 0


def test_variants_are_nested_and_top_n_is_capped():
    parents = {"T1": set(), "T2": set(), "T3": set()}
    ann = {f"G{i}": {"T1", "T2", "T3"} for i in range(5)} | {f"H{i}": {"T1", "T2"} for i in range(200)}
    module = build_phenotype_module(["T1", "T2", "T3"], ann, parents)
    v = module["variants"]
    assert set(v["min3_terms"]) <= set(v["min2_terms"]) and len(v["min2_terms"]) == 205 and len(v["min3_terms"]) == 5
    assert len(v[f"top{TOP_N}_information_content"]) == TOP_N
