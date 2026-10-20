import numpy as np

from track2.core_n2 import compare, restrict_links
from track2.kg_extra import phenotype_neighbours, recall_at
from track2.knowledge_graph import Graph
from track2.vulnerability_scan import n_classes


def test_restrict_links_keeps_only_links_with_two_records_and_drops_empty_drugs():
    e1 = {"D1": {"P1": {"n": 1}, "P2": {"n": 5}}, "D2": {"P3": {"n": 1}}, "D3": {"P4": {"n": 2}}}
    assert restrict_links(e1) == {"D1": {"P2": {"n": 5}}, "D3": {"P4": {"n": 2}}}


def test_compare_marks_drugs_that_lose_robustness():
    original = {"n_robust_core": 2, "drugs": [{"drug_id": "A", "name": "a", "robustness_count_core": 9, "L1_core": True}, {"drug_id": "B", "name": "b", "robustness_count_core": 8, "L1_core": True},
                                              {"drug_id": "C", "name": "c", "robustness_count_core": 3, "L1_core": False}]}
    strict = {"n_robust_core": 1, "rule": {"robust_min": 7}, "drugs": [{"drug_id": "A", "name": "a", "robustness_count_core": 9, "L1_core": True}, {"drug_id": "B", "name": "b", "robustness_count_core": 4, "L1_core": False}]}
    res = compare(original, strict)
    by = {r["name"]: r for r in res["originally_robust"]}
    assert by["a"]["robust_n2_no_E2"] and not by["b"]["robust_n2_no_E2"] and "c" not in by and res["robust_in_strict_run"] == ["a"]


def test_recall_at_and_class_collapse():
    s = np.array([0.9, 0.8, 0.1, 0.05])
    assert recall_at(s, np.array([True, False, True, False]), top=2) == 0.5
    moa = {"r1": "MDM INHIBITOR", "r2": "MDM INHIBITOR", "r3": "", "r4": "EGFR INHIBITOR"}
    assert n_classes(["r1", "r2", "r3", "r4"], moa) == 3  # two MDM inhibitors count once; the unannotated hit is its own class


def test_phenotype_neighbours_reads_the_raw_relation(tmp_path):
    (tmp_path / "e.csv").write_text("relation,display_relation,x_index,y_index\ndisease_phenotype_positive,x,0,1\ndisease_phenotype_positive,x,2,0\ndisease_phenotype_negative,x,0,3\n", encoding="utf-8")
    g = Graph(np.array(["disease", "effect/phenotype", "effect/phenotype", "effect/phenotype"]), np.array(["a"] * 4), np.array(["a"] * 4), None, {})
    assert sorted(phenotype_neighbours(tmp_path / "e.csv", g, {0})[0]) == [1, 2]


def test_degree_matched_random_queries_have_the_same_shape():
    from scipy import sparse

    from track2.kg_nullquery import draw_queries
    from track2.knowledge_graph import Graph

    # a star graph: disease nodes 0-3 have degree 5, phenotype nodes 4-9 degree 5, a hub (10) connects everything
    n = 11
    edges = [(i, 10) for i in range(10)] + [(i, (i + 1) % 10) for i in range(10)]
    a = sparse.coo_matrix((np.ones(len(edges)), ([e[0] for e in edges], [e[1] for e in edges])), shape=(n, n)).tocsr()
    a = a.maximum(a.T)
    types = np.array(["disease"] * 4 + ["effect/phenotype"] * 6 + ["gene/protein"])
    names = np.array([f"n{i}" for i in range(n)])
    g = Graph(types, names, names, a, {})
    qs = draw_queries(g, disease=0, phenos=[4, 5], n=20)
    assert len(qs) == 20 and all(len(q) == 3 for q in qs)
    assert all(0 not in q and 4 not in q and 5 not in q for q in qs)  # the real nodes are excluded from the random draws


def test_causal_genes_and_equal_weight_combination(tmp_path):
    from track2.module_rule_validation import causal_genes, combine_equal

    (tmp_path / "g.txt").write_text("ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n1\tAAA\tHP:1\tx\t-\tOMIM:1\n2\tBBB\tHP:2\tx\t-\tOMIM:1\n3\t-\tHP:3\tx\t-\tOMIM:1\n4\tCCC\tHP:4\tx\t-\tOMIM:2\n", encoding="utf-8")
    assert causal_genes(tmp_path / "g.txt", {"OMIM:1", "OMIM:3"}) == {"OMIM:1": {"AAA", "BBB"}, "OMIM:3": set()}
    z = combine_equal({"AAA": {"d1": -3.0, "d2": 1.0}, "BBB": {"d1": -1.0, "d2": -2.0}})
    assert z == {"d1": -4.0, "d2": -2.0}  # the drug close to both modules gets the lowest pseudo-z, so the highest score
