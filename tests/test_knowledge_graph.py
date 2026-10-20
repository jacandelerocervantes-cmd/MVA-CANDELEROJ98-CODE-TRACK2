import numpy as np

from track2.kg_validation import sample_diseases, validate
from track2.knowledge_graph import Graph, load_graph, personalised_pagerank, restart_vector

NODES = "node_index,node_id,node_type,node_name,node_source\n" + "\n".join([
    "0,D1,disease,disease one,x", "1,G1,gene/protein,gene one,x", "2,G2,gene/protein,gene two,x",
    "3,DB1,drug,drug A,x", "4,DB2,drug,drug B,x", "5,DB3,drug,drug C,x", "6,DB4,drug,drug D,x", "7,P1,effect/phenotype,pheno,x"]) + "\n"
EDGES = "relation,display_relation,x_index,y_index\n" + "\n".join([
    "disease_protein,associated with,0,1", "disease_protein,associated with,0,2", "drug_protein,target,3,1", "drug_protein,target,4,2",
    "drug_protein,target,5,7", "drug_effect,side effect,6,7",
    "indication,indication,3,0", "indication,indication,4,0", "contraindication,contraindication,5,0", "off-label use,off-label use,6,0"]) + "\n"


def _graph(tmp_path):
    (tmp_path / "n.csv").write_text(NODES, encoding="utf-8")
    (tmp_path / "e.csv").write_text(EDGES, encoding="utf-8")
    return load_graph(tmp_path / "n.csv", tmp_path / "e.csv")


def test_drug_disease_edges_are_kept_out_of_the_graph_and_only_indications_are_ground_truth(tmp_path):
    g = _graph(tmp_path)
    assert g.adj[3, 0] == 0 and g.adj[0, 3] == 0 and g.adj[5, 0] == 0 and g.adj[6, 0] == 0  # none of the three drug-disease relations is walkable
    assert g.indications == {0: {3, 4}}                                                     # contraindication and off-label use are not ground truth
    assert g.adj[6, 7] == 0 and g.adj[7, 6] == 0                                            # a side-effect (drug_effect) edge is not walkable
    assert g.adj[1, 3] == 1 and g.adj[3, 1] == 1                                            # symmetric


def test_pagerank_is_a_distribution_and_reaches_drugs_through_genes(tmp_path):
    g = _graph(tmp_path)
    p = personalised_pagerank(g.adj, restart_vector(8, {0: 1.0})[:, None])
    assert abs(p[:, 0].sum() - 1) < 1e-6 and p[3, 0] > p[5, 0] and p[3, 0] > p[6, 0]  # drugs sharing genes with the disease outrank the others


def test_validation_recovers_hidden_indications_in_a_synthetic_graph():
    from scipy import sparse

    n_dis, per_dis, n_other = 6, 3, 30
    # nodes: diseases, one gene per (disease, drug), the drugs of each disease, and distractor drugs with their own genes
    types, edges, indications = [], [], {}
    for d in range(n_dis):
        dn = len(types); types.append("disease")
        indications[dn] = set()
        for _ in range(per_dis):
            gn = len(types); types.append("gene/protein")
            drug = len(types); types.append("drug")
            edges += [(dn, gn), (drug, gn)]
            indications[dn].add(drug)
    for _ in range(n_other):
        gn = len(types); types.append("gene/protein")
        drug = len(types); types.append("drug")
        edges.append((drug, gn))
    n = len(types)
    a = sparse.coo_matrix((np.ones(len(edges)), ([e[0] for e in edges], [e[1] for e in edges])), shape=(n, n)).tocsr()
    a = a.maximum(a.T)
    names = np.array([f"n{i}" for i in range(n)])
    g = Graph(np.array(types), names, names, a, indications)
    diseases = sample_diseases(g)
    assert len(diseases) == n_dis
    res = validate(g, diseases)
    assert res["n_evaluated"] == n_dis and res["summary"]["mean_auc"] > 0.95 and res["validated"]
