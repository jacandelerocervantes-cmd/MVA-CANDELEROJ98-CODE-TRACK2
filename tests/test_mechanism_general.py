"""Offline tests of the general mechanism rule table (synthetic term lists, no network)."""
import re

from track2 import mechanism_general as mg


def cls(terms):
    return mg.classify_general(terms)[0]


def test_empty_and_unmatched_terms_are_other():
    assert mg.classify_general([]) == ("other", [])
    assert cls(["Reference proteome", "Alternative splicing", "cytoplasm"]) == "other"


def test_single_weak_term_is_below_threshold():
    assert cls(["Hydrolase"]) == "other"  # weak term only: score 1 < MIN_SCORE
    assert cls(["Hydrolase", "Transferase"]) == "metabolic_enzyme"  # two weak terms reach the threshold


def test_each_class_has_a_positive_case():
    cases = {
        "lysosomal_vesicular_trafficking": ["Lysosome", "lysosomal membrane", "autophagy"],
        "ion_channel_transporter": ["Ion channel", "Sodium channel", "ion transport"],
        "metabolic_enzyme": ["Lipid metabolism", "oxidoreductase activity", "fatty acid beta-oxidation"],
        "receptor_signalling": ["Receptor", "signal transduction", "G protein-coupled receptor signaling pathway"],
        "transcription_chromatin_regulation": ["Transcription regulation", "chromatin binding", "histone binding"],
        "structural_cytoskeletal": ["Cytoskeleton", "actin binding", "microtubule"],
        "mitochondrial": ["respiratory chain complex I assembly", "mitochondrial translation", "mitochondrion"],
        "dna_repair_replication": ["DNA repair", "double-strand break repair", "DNA replication"],
        "cell_cycle_mitotic_checkpoint": ["Mitosis", "spindle assembly checkpoint signaling", "kinetochore"],
        "rna_processing_translation": ["Ribosomal protein", "structural constituent of ribosome", "rRNA processing"],
        "extracellular_secreted": ["Secreted", "extracellular matrix", "blood coagulation"],
    }
    assert set(cases) == set(mg.ORDER)
    for expected, terms in cases.items():
        assert cls(terms) == expected, expected


def test_classes_constant_lists_all_rule_classes_and_other():
    assert set(mg.RULES) == set(mg.ORDER)
    assert mg.CLASSES == mg.ORDER + ("other",)


def test_term_order_and_duplicates_do_not_matter():
    t = ["DNA repair", "Mitosis", "DNA replication", "Cell cycle", "Mitosis"]
    assert mg.classify_general(t) == mg.classify_general(list(reversed(t)))
    assert mg.score_terms(t) == mg.score_terms(sorted(set(t)))


def test_ties_are_broken_by_the_fixed_order():
    # one strong term for each of two classes -> equal score 2; the class earlier in ORDER wins
    a, b = "DNA repair", "ion transport"
    s = mg.score_terms([a, b])
    assert s["dna_repair_replication"] == s["ion_channel_transporter"] == 2
    assert cls([a, b]) == "dna_repair_replication"
    assert mg.ORDER.index("dna_repair_replication") < mg.ORDER.index("ion_channel_transporter")


def test_class_score_is_capped():
    many = [f"regulation of transcription step {i}" for i in range(30)]
    assert mg.score_terms(many)["transcription_chromatin_regulation"] == mg.MAX_CLASS_SCORE
    # a specific class with the cap-level score wins the tie against a generic class by order
    assert cls(many + ["DNA repair", "DNA replication", "double-strand break repair"]) == "dna_repair_replication"


def test_a_term_counts_once_per_class_with_its_highest_weight():
    assert mg.score_terms(["DNA repair"])["dna_repair_replication"] == 2  # matches a strong and no second pattern
    assert mg.score_terms(["DNA damage"])["dna_repair_replication"] == 1  # weak only


def test_post_translational_is_not_translation():
    assert mg.score_terms(["post-translational modification", "translational"])["rna_processing_translation"] == 0
    assert mg.score_terms(["translation"])["rna_processing_translation"] == 2


def test_structural_constituent_of_ribosome_is_not_structural():
    s = mg.score_terms(["structural constituent of ribosome"])
    assert s["structural_cytoskeletal"] == 0 and s["rna_processing_translation"] == 2


def test_protein_enzyme_activities_are_not_metabolic_enzymes():
    assert mg.score_terms(["protein serine kinase activity"])["metabolic_enzyme"] == 0
    assert mg.score_terms(["ubiquitin-protein transferase activity"])["metabolic_enzyme"] == 0
    assert mg.score_terms(["alcohol dehydrogenase activity"])["metabolic_enzyme"] == 2


def test_matched_terms_are_sorted_and_belong_to_the_class():
    c, m = mg.classify_general(["Mitosis", "kinetochore", "Reference proteome"])
    assert c == "cell_cycle_mitotic_checkpoint" and m == ["Mitosis", "kinetochore"]


def test_rules_contain_no_gene_like_or_disease_tokens():
    # the table must be gene- and disease-agnostic: no all-caps gene-like symbols with digits, no 'syndrome'/'disease'
    for cls_name, rules in mg.RULES.items():
        for pattern, weight in rules:
            assert weight in (1, 2)
            re.compile(pattern)
            assert not re.search(r"\b[A-Z]{2,}[0-9]+[A-Z0-9]*\b", pattern.replace("G1/S", "").replace("G2/M", "").replace("SH2", "").replace("SH3", "")), (cls_name, pattern)
            assert not re.search(r"syndrome|disease|anemia|anaemia|dystrophy", pattern, re.I), (cls_name, pattern)
