from __future__ import annotations

from track2 import caution_list as cl

SECTIONS = {
    "boxed_warning": "WARNING: SEVERE TOXICITY. Fatal if given by other routes. For intravenous use only.",
    "contraindications": "Do not use in patients with active malignancy. Hypersensitivity is also a contraindication.",
    "pediatric_use": "Safety and effectiveness in pediatric patients have not been established. Growth plate findings were seen in animals.",
    "warnings": "Severe myelosuppression can occur. Secondary malignancies have been reported. Fetal harm can occur in pregnancy.",
}


def _fetch(name: str):
    if name == "NOCITE":
        return {"source": "DailyMed", "setid": "s0", "title": "NOCITE", "label_effective": "20200101",
                "sections": {"boxed_warning": "", "contraindications": "None.", "pediatric_use": "", "warnings": "Take with food."}}, []
    if name == "BROKEN":
        return None, ["DailyMed: boom", "openFDA: HTTP 500"]
    return {"source": "DailyMed", "setid": "s1", "title": name, "label_effective": "20200101", "sections": SECTIONS}, []


def _drug(i, name, atc, mech=()):
    return {"drug_id": i, "name": name, "atc": atc, "mechanisms": [{"mechanism": m} for m in mech]}


DRUGS = [
    _drug("C1", "VINCRISTINE", ["L01CA02"]), _drug("C2", "NOCITE", ["L01CD01"]), _drug("C3", "BROKEN", ["L01XA01"]),
    _drug("C4", "SOMATROPIN", ["H01AC01"], ["Growth hormone receptor agonist"]), _drug("C5", "MECASERMIN", ["H01AC03"], ["Insulin-like growth factor I receptor agonist"]),
    _drug("C6", "BORTEZOMIB", ["L01XG01"]), _drug("C7", "ALBENDAZOLE", ["P02CA03"], ["Tubulin inhibitor"]),
    _drug("C8", "ERIBULIN", ["L01XX41"]), _drug("C9", "AMSACRINE", ["L01XX01"], ["DNA topoisomerase II inhibitor"]),
    _drug("C10", "DOCETAXEL ANHYDROUS", [], ["Tubulin stabiliser"]),
]


def test_class_membership_and_empty_class() -> None:
    m = cl.classify_universe(DRUGS)
    assert [x["name"] for x in m["a_tubulina"]] == ["DOCETAXEL ANHYDROUS", "ERIBULIN", "NOCITE", "VINCRISTINE"]
    assert [x["name"] for x in m["e_hormona_crecimiento"]] == ["SOMATROPIN"]  # mecasermin is not GH
    assert m["b_checkpoint_mitotico"] == []
    assert {x["name"] for x in m["c_quimioterapia_genotoxica"]} == {"BROKEN", "AMSACRINE"}
    assert [x["name"] for x in m["d_proteasoma"]] == ["BORTEZOMIB"]


def test_every_quote_is_literal_and_uncited_drugs_are_not_listed() -> None:
    out = cl.build(DRUGS, "2026-10-07", fetch=_fetch)
    assert out["aviso"].startswith("hypotheses, not advice")
    assert out["clases"]["b_checkpoint_mitotico"]["clase_vacia"] is True
    tub = out["clases"]["a_tubulina"]
    assert "NOCITE" not in [d["name"] for d in tub["farmacos"]]
    assert {"clase": "a_tubulina", "farmaco": "NOCITE"}.items() <= {k: v for k, v in out["sin_cita"][0].items()}.items() or any(s["farmaco"] == "NOCITE" for s in out["sin_cita"])
    assert [u["farmaco"] for u in out["etiquetas_ilegibles"]] == ["BROKEN"]
    n = 0
    for cls in out["clases"].values():
        for d in cls["farmacos"]:
            assert d["citas"], "listed without citation"
            for q in d["citas"]:
                assert q["texto_literal"] in SECTIONS[q["seccion"]]
                assert q["setid"] == "s1" and q["fecha_lectura"] == "2026-10-07"
                n += 1
    assert n > 0 and cl.revalidate(out["clases"]["a_tubulina"]["farmacos"][0]["citas"], SECTIONS)


def test_growth_hormone_active_malignancy_and_categories() -> None:
    out = cl.build(DRUGS, "2026-10-07", fetch=_fetch)
    gh = out["clases"]["e_hormona_crecimiento"]["farmacos"][0]
    cats = {q["categoria"] for q in gh["citas"]}
    assert "active_malignancy" in cats
    assert any(q["texto_literal"] == "Do not use in patients with active malignancy." for q in gh["citas"])
    vin = out["clases"]["a_tubulina"]["farmacos"][-1]
    assert {"toxicidad_medular", "embrio_fetal", "boxed_warning", "pediatric_use"} <= {q["categoria"] for q in vin["citas"]}


def test_verify_quote_rejects_altered_text() -> None:
    assert cl.verify_quote("Severe myelosuppression can occur.", SECTIONS["warnings"])
    assert not cl.verify_quote("Severe myelosuppression can occur!", SECTIONS["warnings"])
    assert not cl.verify_quote("", SECTIONS["warnings"])
