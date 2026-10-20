"""E11 caution list (protocol/addendum_review2.md): pre-specified drug classes whose use in this patient needs explicit caution.

A drug enters the list only with at least one VERBATIM label sentence (boxed warning, contraindications, paediatric use,
bone-marrow toxicity, second malignancy, embryo-fetal toxicity; for somatropin also 'active malignancy'). Each quote is
a literal substring of the label text read (whitespace normalised exactly as `labels._sections` does), with setid and read date.
A drug with no quote is recorded in `sin_cita` and is NOT listed. Only approved drugs (data/universe.json, ChEMBL max_phase 4).

HYPOTHESIS, NOT ADVICE: no dosing, no recommendation to treat. The age of the patient is unknown; nothing here assumes an age band.
"""

from __future__ import annotations

import datetime as dt
import json
import re
import sys
import urllib.parse
from pathlib import Path
from typing import Any, Callable

from track2 import labels
from track2.http_cache import get_json

NOTICE = "hypotheses, not advice: caution list for research and for a clinician; not a treatment or dose recommendation. Verbatim label sentences; any dose figure is printed as in the label and is not dosing information."
UNIVERSE = Path("data/universe.json")
OUT = Path("out/caution_list.json")

# LOINC codes of the label sections that are read (labels.SECTIONS has a smaller set).
LOINC = {
    "boxed_warning": ("34066-1",),
    "contraindications": ("34070-3",),
    "pediatric_use": ("34081-0",),
    "warnings": ("43685-7", "34071-1", "42232-9", "43684-0"),
}
SECTION_ORDER = ("boxed_warning", "contraindications", "pediatric_use", "warnings")

# Quote categories -> (sections searched, regex). Pre-specified before reading any label text.
CATEGORIES: dict[str, tuple[tuple[str, ...], str]] = {
    "boxed_warning": (("boxed_warning",), r"."),
    "contraindications": (("contraindications",), r"."),
    "pediatric_use": (("pediatric_use",), r"."),
    "toxicidad_medular": (("boxed_warning", "warnings"), r"myelosuppress|bone marrow|neutropeni|pancytopeni|thrombocytopeni|cytopeni"),
    "segundos_tumores": (("boxed_warning", "warnings"), r"secondary (malignan|acute|leukem|cancer)|second (primary )?(malignan|cancer)|(new|secondary) (malignan|primary)|carcinogenic|leukemogenic|myelodysplastic"),
    "embrio_fetal": (("boxed_warning", "contraindications", "warnings"), r"embryo|fetal (harm|toxicit)|harm (to )?(a|the) fetus|cause fetal"),
    "active_malignancy": (("contraindications", "boxed_warning", "warnings"), r"active malignan|active tumou?r|malignan|neoplas"),
}
MAX_QUOTE = 600
MAX_PER_CATEGORY = 2

TUBULIN = "tubulina"
CLASSES: dict[str, dict[str, Any]] = {
    "a_tubulina": {
        "etiqueta": "Tubulin binders",
        "criterio": "ATC L01CA (vinca), L01CD (taxanes), M04AC01 (colchicine); eribulin and ixabepilone by name; ATC L01XX/L01DC with ChEMBL mechanism 'Tubulin ...'; "
                    "docetaxel (taxane, mechanism 'Tubulin stabiliser', no ATC in the universe) by mechanism and the '-taxel' stem.",
        "categorias": ["boxed_warning", "contraindications", "pediatric_use", "toxicidad_medular", "embrio_fetal"],
    },
    "b_checkpoint_mitotico": {
        "etiqueta": "Mitotic-checkpoint, Aurora, PLK and TTK (MPS1) inhibitors",
        "criterio": "ChEMBL mechanism containing Aurora, Polo-like kinase, TTK/MPS1, BUB1, MAD2 or 'mitotic checkpoint' in data/universe.json (approved drugs).",
        "categorias": ["boxed_warning", "contraindications", "pediatric_use", "toxicidad_medular", "embrio_fetal"],
    },
    "c_quimioterapia_genotoxica": {
        "etiqueta": "Genotoxic chemotherapy",
        "criterio": "ATC L01A* (alkylating agents), L01XA* (platinum), L01DB* (anthracyclines), L01CB* (epipodophyllotoxins, topoisomerase II), L01XX with a topoisomerase mechanism; "
                    "etoposide phosphate (no ATC; topoisomerase II mechanism, antineoplastic) by mechanism and the name stem.",
        "categorias": ["boxed_warning", "contraindications", "pediatric_use", "toxicidad_medular", "segundos_tumores", "embrio_fetal"],
    },
    "d_proteasoma": {
        "etiqueta": "Proteasome inhibitors",
        "criterio": "ATC L01XG*; or ChEMBL mechanism 'proteasome inhibitor' without ATC (ixazomib).",
        "categorias": ["boxed_warning", "contraindications", "pediatric_use", "toxicidad_medular", "embrio_fetal"],
    },
    "e_hormona_crecimiento": {
        "etiqueta": "Growth hormone and somatropin analogues",
        "criterio": "ATC H01AC* with mechanism 'Growth hormone receptor agonist' (somatropin and GH analogues). Sermorelin, mecasermin and tesamorelin are excluded (not GH).",
        "categorias": ["active_malignancy", "contraindications", "pediatric_use", "boxed_warning"],
    },
}
NAME_STEMS_OK = {"DOCETAXEL ANHYDROUS": "a_tubulina", "ETOPOSIDE PHOSPHATE": "c_quimioterapia_genotoxica"}
EXPLICIT_TUBULIN_NAMES = {"ERIBULIN", "IXABEPILONE"}
EXCLUDED_DOCUMENTED = {
    "antibody-drug conjugates with a tubulin payload (brentuximab, polatuzumab, enfortumab, trastuzumab emtansine, tisotumab, mirvetuximab)": "not in the pre-specified definition (ATC L01XX); they are conjugated antibodies, L01F*",
    "albendazole, mebendazole, thiabendazole and griseofulvin (tubulin mechanism, antiparasitic or antifungal use)": "not in the pre-specified classes (ATC P02/D01)",
    "dexrazoxane and quinolones (topoisomerase II or bacterial)": "not pre-specified genotoxic chemotherapy",
}
LABEL_NOT_REPRESENTATIVE = {  # drugs whose only US label found describes a different product than the systemic drug of the class: not listed, with the reason recorded
    "MECHLORETHAMINE": "the only US label found is Valchlor, a topical gel for cutaneous T-cell lymphoma; it does not describe the systemic drug",
}
PREFERRED_TITLES = {  # label choice for drugs whose first DailyMed hit is an unrepresentative product (documented, not data-driven)
    "SOMATROPIN": ["humatrope", "genotropin", "norditropin", "omnitrope", "nutropin"],
}
# Added after the manual label review of 2026-10-08, which found that the first DailyMed hit was the wrong product for these drugs: (search term, text the title must contain, text it must not contain)
PRODUCT_OVERRIDES = {
    "CYCLOSPORINE": ("sandimmune", "SANDIMMUNE", None),            # the first hit was an ophthalmic emulsion; the systemic product is the one relevant here
    "DOXORUBICIN": ("doxorubicin hydrochloride", "DOXORUBICIN", "LIPOSOM"),   # the first hit was the liposomal product
    "ETOPOSIDE PHOSPHATE": ("etopophos", "ETOPOPHOS", None),       # the first hit was etoposide injection, a different product
    "MIDOSTAURIN": ("rydapt", "RYDAPT", None),                     # a search by generic name returns nothing
    "CARMUSTINE": ("carmustine", "CARMUSTINE", "WAFER"),           # the first hit was the Gliadel wafer (an implant), not the intravenous product
}


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def classify_universe(drugs: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Members of each pre-specified class (approved drugs of the universe only). Pure function; classes may be empty."""
    out: dict[str, list[dict[str, Any]]] = {k: [] for k in CLASSES}

    def add(cls: str, d: dict[str, Any], why: str) -> None:
        if all(m["drug_id"] != d["drug_id"] for m in out[cls]):
            out[cls].append({"drug_id": d["drug_id"], "name": d["name"], "atc": d.get("atc", []), "motivo": why})

    for d in drugs:
        atc = d.get("atc") or []
        mech = [m.get("mechanism", "") for m in d.get("mechanisms", [])]
        mtext = " | ".join(mech)
        has = lambda *ps: [a for a in atc if any(a.startswith(p) for p in ps)]  # noqa: E731
        if has("L01CA", "L01CD", "M04AC01"):
            add("a_tubulina", d, f"ATC {has('L01CA', 'L01CD', 'M04AC01')[0]}")
        elif d["name"] in EXPLICIT_TUBULIN_NAMES:
            add("a_tubulina", d, "eribulin or ixabepilone (name given in the protocol)")
        elif has("L01XX", "L01DC") and re.search(r"tubulin", mtext, re.I):
            add("a_tubulina", d, f"ATC {has('L01XX', 'L01DC')[0]} + mecanismo ChEMBL tubulina")
        elif d["name"] in NAME_STEMS_OK and NAME_STEMS_OK[d["name"]] == "a_tubulina" and re.search(r"tubulin", mtext, re.I):
            add("a_tubulina", d, "ChEMBL tubulin mechanism; taxane without ATC in the universe")
        if re.search(r"aurora|polo-like|\bPLK|\bTTK\b|\bMPS1\b|\bBUB1|\bMAD2|mitotic checkpoint|spindle assembly", mtext, re.I) and (
                not atc or has("L01")):
            add("b_checkpoint_mitotico", d, "mecanismo ChEMBL: " + mtext[:120])
        if has("L01A"):
            add("c_quimioterapia_genotoxica", d, f"ATC {has('L01A')[0]} (alquilante)")
        elif has("L01XA"):
            add("c_quimioterapia_genotoxica", d, f"ATC {has('L01XA')[0]} (platino)")
        elif has("L01DB"):
            add("c_quimioterapia_genotoxica", d, f"ATC {has('L01DB')[0]} (antraciclina)")
        elif has("L01CB"):
            add("c_quimioterapia_genotoxica", d, f"ATC {has('L01CB')[0]} (topoisomerasa II)")
        elif has("L01XX") and re.search(r"topoisomerase", mtext, re.I):
            add("c_quimioterapia_genotoxica", d, f"ATC {has('L01XX')[0]} + mecanismo topoisomerasa")
        elif NAME_STEMS_OK.get(d["name"]) == "c_quimioterapia_genotoxica" and re.search(r"topoisomerase", mtext, re.I):
            add("c_quimioterapia_genotoxica", d, "topoisomerase II mechanism; antineoplastic without ATC in the universe")
        if has("L01XG") or (not atc and re.search(r"proteasome inhibitor", mtext, re.I)):
            add("d_proteasoma", d, f"ATC {has('L01XG')[0]}" if has("L01XG") else "proteasome mechanism, no ATC")
        if has("H01AC") and re.search(r"growth hormone receptor agonist", mtext, re.I):
            add("e_hormona_crecimiento", d, f"ATC {has('H01AC')[0]} + agonista del receptor de GH")
    for k in out:
        out[k].sort(key=lambda m: m["name"])
    return out


# ---------------------------------------------------------------- label text -------------------------------------------------
def _read_dailymed(name: str) -> dict[str, Any]:
    """Label text from DailyMed through the cache. Raises on a service error or if no label is found."""
    found = None
    for pref in PREFERRED_TITLES.get(name, []):
        data = get_json(f"{labels.DAILYMED}/spls.json?drug_name={urllib.parse.quote(pref)}&pagesize=25")
        for item in data.get("data", []):
            if pref in item["title"].lower() and "somatropin" in item["title"].lower():
                found = (item["setid"], item["title"])
                break
        if found:
            break
    if not found and name in PRODUCT_OVERRIDES:
        query, must, must_not = PRODUCT_OVERRIDES[name]
        data = get_json(f"{labels.DAILYMED}/spls.json?drug_name={urllib.parse.quote(query)}&pagesize=25")
        for item in data.get("data", []):
            title = item["title"].upper()
            if must in title and not (must_not and must_not in title):
                found = (item["setid"], item["title"])
                break
    if not found:
        base = labels._clean_name(name).lower()
        data = get_json(f"{labels.DAILYMED}/spls.json?drug_name={urllib.parse.quote(base)}&pagesize=25")
        items = [i for i in data.get("data", []) if base.split()[0] in i["title"].lower()]
        items.sort(key=lambda i: (" AND " in i["title"].upper() or " / " in i["title"] or "KIT" in i["title"].upper(), not i["title"].lower().startswith(base.split()[0])))
        if items:
            found = (items[0]["setid"], items[0]["title"])
    if not found:
        raise LookupError(f"DailyMed: no label for {name}")
    return _dailymed_by_setid(*found)


def _dailymed_by_setid(setid: str, title: str | None = None) -> dict[str, Any]:
    """Label text of one known DailyMed set id (the product the paediatric review read), through the cache."""
    root = labels._spl_xml(setid)
    if title is None:
        node = root.find(labels.H + "title")
        org = root.find(".//" + labels.H + "representedOrganization/" + labels.H + "name")
        title = (" ".join("".join(node.itertext()).split()) if node is not None else setid)
        # same shape as the DailyMed list title: "NAME FORM [MANUFACTURER]"
        name = root.find(".//" + labels.H + "manufacturedProduct//" + labels.H + "name")
        title = f"{(name.text or '').strip() if name is not None else title} [{org.text.strip().upper() if org is not None and org.text else '?'}]"
    sections: dict[str, str] = {k: "" for k in LOINC}
    for section in root.iter(labels.H + "section"):
        code = section.find(labels.H + "code")
        if code is None:
            continue
        for key, codes in LOINC.items():
            if code.get("code") in codes:
                text = _norm("".join(section.itertext()))
                if text and text not in sections[key]:
                    sections[key] = (sections[key] + " " + text).strip()
    eff = root.find(labels.H + "effectiveTime")
    return {"source": "DailyMed", "setid": setid, "title": title, "label_effective": eff.get("value") if eff is not None else None, "sections": sections}


def fetch_by_setid(setid: str) -> tuple[dict[str, Any] | None, list[str]]:
    try:
        lab = _dailymed_by_setid(setid)
    except Exception as error:  # noqa: BLE001 - reported
        return None, [f"dailymed set id {setid}: {error!r}"[:240]]
    return (lab, []) if any(lab["sections"].values()) else (None, [f"set id {setid}: label found but no readable text"])


def _read_openfda(name: str) -> dict[str, Any]:
    """Fallback: openFDA drug label API (through the cache). Raises on a service error or if nothing matches."""
    base = labels._clean_name(name).lower().split()[0]
    q = urllib.parse.quote(f'openfda.generic_name:"{base}"')
    data = get_json(f"https://api.fda.gov/drug/label.json?search={q}&limit=1")
    res = data["results"][0]
    pick = lambda *keys: _norm(" ".join(" ".join(res.get(k, [])) for k in keys))  # noqa: E731
    sections = {"boxed_warning": pick("boxed_warning"), "contraindications": pick("contraindications"),
                "pediatric_use": pick("pediatric_use"), "warnings": pick("warnings_and_cautions", "warnings", "precautions")}
    return {"source": "openFDA", "setid": (res.get("openfda", {}).get("spl_set_id") or [res.get("set_id")])[0], "title": base,
            "label_effective": res.get("effective_time"), "sections": sections}


def fetch_label(name: str, readers: tuple[Callable[[str], dict[str, Any]], ...] = (_read_dailymed, _read_openfda)) -> tuple[dict[str, Any] | None, list[str]]:
    """First reader that gives a label with some text wins. Returns (label or None, list of failures); nothing is invented."""
    failures: list[str] = []
    for reader in readers:
        try:
            lab = reader(name)
        except Exception as error:  # noqa: BLE001 - reported, never read as "no warning"
            failures.append(f"{getattr(reader, '__name__', 'reader')}: {error!r}"[:240])
            continue
        if any(lab["sections"].values()):
            return lab, failures
        failures.append(f"{getattr(reader, '__name__', 'reader')}: label found but no readable text")
    return None, failures


# ---------------------------------------------------------------- quotes -----------------------------------------------------
def sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?])\s+(?=[A-Z0-9\[(“\"])", text) if s]


def verify_quote(quote: str, source_text: str) -> bool:
    """A quote is valid only if it appears LITERALLY in the source text."""
    return bool(quote) and quote in source_text


def extract_quotes(sections: dict[str, str], categories: list[str], read_on: str, ref: dict[str, Any]) -> list[dict[str, Any]]:
    quotes: list[dict[str, Any]] = []
    for cat in categories:
        secs, pattern = CATEGORIES[cat]
        found = 0
        for sec in secs:
            text = sections.get(sec, "")
            if not text or found >= MAX_PER_CATEGORY:
                continue
            sents = sentences(text)
            if cat in ("boxed_warning", "contraindications", "pediatric_use"):
                cand = [s for s in sents if len(s) >= 25][:MAX_PER_CATEGORY] if cat != "pediatric_use" else [
                    s for s in sents if re.search(r"pediatric|children|neonat|infant|adolescent", s, re.I)][:MAX_PER_CATEGORY]
            else:
                cand = [s for s in sents if re.search(pattern, s, re.I)][:MAX_PER_CATEGORY - found]
            for s in cand:
                quote = s[:MAX_QUOTE]
                if verify_quote(quote, text):
                    quotes.append({"categoria": cat, "texto_literal": quote, "truncada": len(s) > MAX_QUOTE, "seccion": sec,
                                   "setid": ref["setid"], "fuente": ref["source"], "titulo_etiqueta": ref["title"],
                                   "fecha_lectura": read_on})
                    found += 1
    return quotes


def revalidate(entry_quotes: list[dict[str, Any]], sections: dict[str, str]) -> bool:
    return all(verify_quote(q["texto_literal"], sections[q["seccion"]]) for q in entry_quotes)


# ---------------------------------------------------------------- build ------------------------------------------------------
def build(drugs: list[dict[str, Any]], read_on: str, fetch: Callable[[str], tuple[dict[str, Any] | None, list[str]]] = fetch_label,
          fda_check: Callable[[str], Any] | None = None) -> dict[str, Any]:
    members = classify_universe(drugs)
    result: dict[str, Any] = {}
    sin_cita: list[dict[str, Any]] = []
    unreadable: list[dict[str, Any]] = []
    cache: dict[str, tuple[dict[str, Any] | None, list[str]]] = {}
    for key, spec in CLASSES.items():
        listed = []
        n_quotes = 0
        for m in members[key]:
            if m["name"] in LABEL_NOT_REPRESENTATIVE:
                sin_cita.append({"clase": key, "farmaco": m["name"], "drug_id": m["drug_id"], "motivo": "label not representative: " + LABEL_NOT_REPRESENTATIVE[m["name"]]})
                continue
            if m["name"] not in cache:
                cache[m["name"]] = fetch(m["name"])
            lab, failures = cache[m["name"]]
            if lab is None:
                unreadable.append({"clase": key, "farmaco": m["name"], "drug_id": m["drug_id"], "fallos": failures})
                sin_cita.append({"clase": key, "farmaco": m["name"], "drug_id": m["drug_id"], "motivo": "label unreadable or not found", "fallos": failures})
                continue
            ref = {k: lab[k] for k in ("source", "setid", "title")}
            quotes = extract_quotes(lab["sections"], spec["categorias"], read_on, ref)
            if not quotes:
                sin_cita.append({"clase": key, "farmaco": m["name"], "drug_id": m["drug_id"], "motivo": "no verifiable verbatim sentence in the sections read",
                                 "setid": lab["setid"], "fuente": lab["source"]})
                continue
            row = {**m, "etiqueta": {"fuente": lab["source"], "setid": lab["setid"], "titulo": lab["title"], "fecha_etiqueta": lab["label_effective"], "fecha_lectura": read_on},
                   "citas": quotes, "lectura_alternativa_tras_fallo": failures or None}
            if fda_check:
                row["fda_prescription"] = fda_check(m["name"])
            listed.append(row)
            n_quotes += len(quotes)
        result[key] = {"etiqueta": spec["etiqueta"], "criterio_de_pertenencia": spec["criterio"], "categorias_buscadas": spec["categorias"],
                       "n_miembros_en_universo": len(members[key]), "n_farmacos_listados": len(listed), "n_citas": n_quotes,
                       "clase_vacia": len(members[key]) == 0, "farmacos": listed,
                       "nota": "empty class: no approved drug of the universe meets the criterion; no members are invented" if not members[key] else None}
    return {"aviso": NOTICE, "protocolo": "protocol/addendum_review2.md E11 (+E14 approved drugs only, E15 no age assumed)", "fecha_lectura": read_on,
            "universo": "data/universe.json (ChEMBL max_phase 4)", "clases": result, "sin_cita": sin_cita, "etiquetas_ilegibles": unreadable,
            "excluidos_documentados": EXCLUDED_DOCUMENTED,
            "regla": "a drug is listed only with at least one verifiable verbatim sentence (exact substring of the label text); without a quote it is not listed"}


def _fda(name: str) -> Any:
    from track2 import fda_approved
    try:
        r = fda_approved.lookup(name)
        return {"eligible": r["eligible"], "applications": r["applications"][:3]}
    except Exception as error:  # noqa: BLE001
        return {"eligible": None, "error": repr(error)[:160]}


def main(argv: list[str] | None = None) -> int:
    drugs = json.loads(UNIVERSE.read_text(encoding="utf-8"))["drugs"]
    read_on = dt.date.today().isoformat()
    result = build(drugs, read_on, fda_check=_fda)
    # Re-check every stored quote against a fresh read of the source text.
    for cls in result["clases"].values():
        for row in cls["farmacos"]:
            lab, _ = fetch_label(row["name"])
            assert lab and revalidate(row["citas"], lab["sections"]), f"quote verification failed for {row['name']}"
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1, ensure_ascii=False, sort_keys=True), encoding="utf-8", newline="\n")
    for k, c in result["clases"].items():
        print(k, "miembros", c["n_miembros_en_universo"], "listados", c["n_farmacos_listados"], "citas", c["n_citas"], "VACIA" if c["clase_vacia"] else "")
    print("sin_cita:", [(s["clase"], s["farmaco"]) for s in result["sin_cita"]])
    return 0


if __name__ == "__main__":
    sys.exit(main())
