"""Offline tests of the claim extractors with synthetic sentences."""

import numpy as np

from track2 import engine_extractors as E


def test_rule_benefit():
    assert E.rule_direction("Treatment with drugx improved muscle function in the animals.", "drugx") == "benefit"
    assert E.rule_direction("DrugX rescued the aneuploid phenotype.", "drugx") == "benefit"
    assert E.rule_direction("DrugX reduced tumour burden and senescence.", "drugx") == "benefit"


def test_rule_harm():
    assert E.rule_direction("DrugX was toxic to normal cells.", "drugx") == "harm"
    assert E.rule_direction("DrugX increased toxicity in treated cells.", "drugx") == "harm"
    assert E.rule_direction("DrugX reduced growth and muscle mass.", "drugx") == "harm"


def test_rule_negation_and_no_effect():
    assert E.rule_direction("DrugX did not improve survival.", "drugx") == "no_effect"
    assert E.rule_direction("DrugX failed to rescue the defect.", "drugx") == "no_effect"
    assert E.rule_direction("There was no significant difference between groups.", "drugx") == "no_effect"
    assert E.rule_direction("DrugX did not cause toxicity.", "drugx") == "no_effect"


def test_rule_unclear_and_scope():
    assert E.rule_direction("Definitive evidence of benefit remains lacking.", "drugx") == "unclear"
    assert E.rule_direction("The weather was fine.", "drugx") == "unclear"
    # scope: the harm belongs to another compound, the clause about the drug is kept
    assert E.rule_direction("DrugX improved function; however, drugy was toxic.", "drugx") == "benefit"


def test_system_and_design():
    assert E.rule_system("Mice treated with the compound lived longer.") == "mouse or rat"
    assert E.rule_system("In Drosophila larvae the effect was seen.") == "invertebrate or yeast"
    assert E.rule_system("Patients receiving the drug had fewer events.") == "human patients"
    assert E.rule_system("HeLa cells were treated.") == "human cells"
    assert E.rule_system("This review summarises the evidence.") == "review or unspecified"
    assert E.rule_design("A randomised, placebo-controlled study.") == "randomised trial"
    assert E.rule_design("An 11-year-old boy was treated.") == "case report"
    assert E.rule_design("Nothing here.") is None


def test_kappa_and_agreement():
    a = np.array(["x", "y", "x", "y"])
    assert E.kappa(a, a) == 1.0
    assert abs(E.kappa(a, np.array(["y", "x", "y", "x"])) + 1.0) < 1e-9
    ref = ["benefit", "harm", "no_effect", "unclear"] * 6
    res = E.agreement(ref, ref, [f"p{i // 2}" for i in range(24)], n_boot=50)
    assert res["accuracy"]["estimate"] == 1.0 and res["kappa"]["ci95"][0] == 1.0


def test_grouped_folds_never_share_pmid():
    groups = [f"p{i // 3}" for i in range(60)]
    y = (["benefit", "harm", "unclear"] * 20)
    for tr, te in E.grouped_folds(y, groups, k=5):
        assert not ({groups[i] for i in tr} & {groups[i] for i in te})
    texts = [("improved rescue beneficial" if v == "benefit" else "toxic harmful impaired" if v == "harm" else "unclear unknown") + f" item{i % 4}" for i, v in enumerate(y)]
    pred = E.out_of_fold(texts, y, groups, k=5)
    assert np.mean(np.array(pred) == np.array(y)) > 0.8


def test_consensus():
    class Fixed:
        def __init__(self, label):
            self.label = label

        def predict(self, texts):
            return [self.label for _ in texts]

    s = "DrugX improved muscle function."
    assert E.consensus(s, "drugx", "may_help", Fixed("harm"))["label"] == "benefit"
    assert E.consensus(s, "drugx", "may_harm", Fixed("no_effect"))["label"] == "unresolved"
    assert E.consensus(s, "drugx", None, Fixed("benefit"))["label"] == "benefit"
