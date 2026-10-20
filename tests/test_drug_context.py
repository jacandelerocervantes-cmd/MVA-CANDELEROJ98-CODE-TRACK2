from __future__ import annotations

from track2 import drug_context


def test_only_phase_4_indications_are_listed() -> None:
    def fake(url: str):
        return {"drug_indications": [{"max_phase_for_ind": "4.0", "mesh_heading": "Epilepsy"}, {"max_phase_for_ind": "1.0", "mesh_heading": "Cancer"},
                                     {"max_phase_for_ind": "4.0", "mesh_heading": "Epilepsy"}], "page_meta": {"next": None}}
    assert drug_context.approved_indications("CHEMBL1", fake) == ["Epilepsy"]


def test_role_text_by_source_and_action() -> None:
    cur = {"source": "curated", "gene": "MTOR", "action": "INHIBITOR", "median_pchembl": None, "similarity": None}
    meas = {"source": "measured", "gene": "PLK1", "action": None, "median_pchembl": 7.4, "similarity": None}
    pred = {"source": "predicted", "gene": "TTK", "action": None, "median_pchembl": None, "similarity": 0.71}
    assert drug_context.role_text(cur) == "inhibits MTOR (curated mechanism)"
    assert "binds PLK1 (measured, median pChEMBL 7.4" in drug_context.role_text(meas)
    assert "predicted to bind TTK" in drug_context.role_text(pred) and "not measured" in drug_context.role_text(pred)


def test_target_roles_merge_modules_and_sort_curated_first() -> None:
    links = {"sub_BUB1B": [{"uniprot": "P2", "gene": "PLK1", "source": "measured", "action": None, "median_pchembl": 7.0, "similarity": None},
                           {"uniprot": "P1", "gene": "MTOR", "source": "curated", "action": "INHIBITOR", "median_pchembl": None, "similarity": None}],
             "sub_CENPE": [{"uniprot": "P1", "gene": "MTOR", "source": "curated", "action": "INHIBITOR", "median_pchembl": None, "similarity": None}]}
    roles = drug_context.target_roles(links)
    assert [r["gene"] for r in roles] == ["MTOR", "PLK1"] and roles[0]["modules"] == ["BUB1B", "CENPE"]


def test_sentence_without_links_does_not_invent_a_mechanism() -> None:
    assert "no link to a module protein" in drug_context.hypothesis_sentence("DRUG X", [], [])
