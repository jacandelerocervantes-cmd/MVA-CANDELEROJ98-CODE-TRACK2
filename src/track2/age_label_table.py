"""E15 age-label table (protocol/addendum_review2.md): label support per age group for every cautioned or candidate drug.

Drugs: (i) all drugs of out/caution_list.json, (ii) class_analysis.json -> supported_in_two_usable_lines members (path is an argument;
if the file is missing a warning is printed and the list is empty), (iii) names in four_lines_review2.json -> readings.*.one_line_only.

The status per group (0-<3, 3-<12, 12+) comes ONLY from the label text (labelled / partly (from X y) / not_labelled / adult_only /
age_not_stated / verify). The age of the patient is unknown and no band is assumed. The automatic classifier output MIXED is `verify`.
A drug without a manual review in config/pediatric_label_review.json has manual_review_required = true; this module never writes that file,
it only proposes the missing entries inside the output.

HYPOTHESIS, NOT ADVICE: label support is not a recommendation and not dosing.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
from pathlib import Path
from typing import Any, Callable

from track2 import age_groups, caution_list, labels
from track2.pediatric_review import load_review, reviewed_profile

NOTICE = "hypotheses, not advice: label support by age group; not a treatment or dose recommendation. The age of the patient is unknown and no age band is assumed. Verbatim label sentences; any dose figure is printed as in the label and is not dosing information."
STATUS_EXPLAINED = {"labelled": "the label supports use across the whole group", "partly": "from a stated minimum age inside the group (reported)",
                    "not_labelled": "paediatric use not established", "adult_only": "adults only", "age_not_stated": "labelled for children without a stated minimum age",
                    "verify": "label mixed or unreadable, or its paediatric section only refers to another section without stating established use: verify by hand"}
SUPPORT_RE = r"established|not recommended|not been|years? (of age )?(and|or) (older|above)|aged? \d|ages? \d|\d+ (years|months|weeks)|neonat|infant|pediatric|children"
WARNING_RE = (r"growth plate|epiphys|growth (retardation|impairment|inhibition|suppression)|cardiomyopath|cardiotox|bone marrow|myelosuppress|neutropeni|"
              r"benzyl alcohol|gasping|ethanol|alcohol content|fatal|death|toxicit|secondary malignan|carcinogen|hepatotox|hypoglyc|scoliosis|slipped capital|"
              r"intracranial hypertension|pancreatitis|neuropath|sudden death|do not use in")


def _group_label(status: str) -> str:
    return "partly" if status.startswith("partly") else status


def pick_support(pediatric_text: str, read_on: str, ref: dict[str, Any]) -> dict[str, Any] | None:
    """First label sentence of the paediatric section that states established / not established / an age (literal)."""
    sents = [s for s in caution_list.sentences(pediatric_text) if len(s) >= 20]
    cand = [s for s in sents if re.search(SUPPORT_RE, s, re.I)] or sents
    if not cand:
        return None
    quote = cand[0][: caution_list.MAX_QUOTE]
    return {"texto_literal": quote, "truncada": len(cand[0]) > caution_list.MAX_QUOTE, "seccion": "pediatric_use", "setid": ref["setid"],
            "fuente": ref["source"], "fecha_lectura": read_on}


def pediatric_warnings(pediatric_text: str, read_on: str, ref: dict[str, Any], limit: int = 4) -> list[dict[str, Any]]:
    out = []
    for s in caution_list.sentences(pediatric_text):
        if re.search(WARNING_RE, s, re.I):
            q = s[: caution_list.MAX_QUOTE]
            out.append({"texto_literal": q, "truncada": len(s) > caution_list.MAX_QUOTE, "seccion": "pediatric_use", "setid": ref["setid"],
                        "fuente": ref["source"], "fecha_lectura": read_on})
        if len(out) >= limit:
            break
    return out


def row_for(drug: dict[str, Any], sources: list[str], lab: dict[str, Any] | None, failures: list[str], review: dict[str, Any], read_on: str) -> dict[str, Any]:
    entry = review.get(drug["drug_id"])
    base = {"drug_id": drug["drug_id"], "name": drug["name"], "origen": sorted(set(sources)),
            "revision_manual_en_config": entry is not None}
    if lab is None:
        return {**base, "etiqueta": None, "fallos_lectura": failures, "estado_por_grupo": {g: "verify" for g in age_groups.GROUPS},
                "frase_literal_sustento": None, "advertencias_pediatricas": [], "clasificador_automatico": None, "manual_review_required": True,
                "nota": "label unreadable or not found in DailyMed/openFDA; nothing is inferred"}
    ref = {k: lab[k] for k in ("source", "setid", "title")}
    ped = lab["sections"].get("pediatric_use", "")
    auto = labels.classify_pediatric(ped)
    profile = labels.band_profile(auto, True)
    profile, prof_source = reviewed_profile(profile, entry)
    verdict = entry["verdict"] if entry and entry["verdict"] in ("established", "not_established") else None
    min_age = entry.get("min_age_years") if verdict == "established" else auto["min_age_years"]
    groups = age_groups.group_status(profile, min_age, verdict)
    if entry and entry["verdict"] in ("verify", "no_us_label"):
        groups = {g: "verify" for g in groups}  # the manual review itself says: not resolvable from the label
    if not ped and verdict is None:
        groups = {g: "verify" for g in groups}  # no paediatric section was readable: the age is never assumed
    needs_manual = verdict is None
    return {**base, "etiqueta": {"fuente": lab["source"], "setid": lab["setid"], "titulo": lab["title"], "fecha_etiqueta": lab["label_effective"], "fecha_lectura": read_on},
            "fallos_lectura": failures or None,
            "clasificador_automatico": {"estado": auto["status"], "edad_minima_anios": auto["min_age_years"], "nota": "MIXED is recorded as verify" if auto["status"] == "MIXED" else None},
            "revision_manual": ({k: entry.get(k) for k in ("verdict", "min_age_years", "scope", "source")} if entry else None),
            "origen_del_perfil": prof_source, "edad_minima_usada_anios": min_age,
            "estado_por_grupo": {g: _group_label(s) if not s.startswith("partly") else s for g, s in groups.items()},
            "labelled_for_group": sorted(g for g, s in groups.items() if s == "labelled"),
            "bandas_etiqueta": {b: profile.get(b, "UNKNOWN") for b in ("neonate", "infant", "child", "adolescent", "adult")},
            "frase_literal_sustento": pick_support(ped, read_on, ref) if ped else None,
            "advertencias_pediatricas": pediatric_warnings(ped, read_on, ref),
            "manual_review_required": needs_manual,
            "nota": None if ped else ("the DailyMed label read has no readable paediatric-use section; the status comes from the manual review in config (another reading, e.g. openFDA)" if verdict else "the label read has no readable paediatric-use section: verify")}


def propose_entries(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Skeleton entries for config/pediatric_label_review.json (NOT written): the author must read the label sentence and decide."""
    out = {}
    for r in rows:
        if r["manual_review_required"] and r.get("etiqueta"):
            out[r["drug_id"]] = {"name": r["name"], "verdict": "PENDIENTE (established | not_established | verify)", "min_age_years": r.get("edad_minima_usada_anios"),
                                 "source": f"{r['etiqueta']['fuente']} setid={r['etiqueta']['setid']}, leído {r['etiqueta']['fecha_lectura']}",
                                 "frase_candidata": (r["frase_literal_sustento"] or {}).get("texto_literal"), "clasificador": r["clasificador_automatico"]["estado"]}
    return out


def gather(caution_path: Path, class_path: Path, four_path: Path, universe: list[dict[str, Any]], proposals_path: Path | None = None) -> tuple[dict[str, list[str]], list[str], list[str]]:
    """Name -> list of origins; warnings; unmatched names."""
    by_name = {d["name"].upper(): d for d in universe}
    origins: dict[str, list[str]] = {}
    warnings: list[str] = []
    if caution_path.exists():
        for ckey, cls in json.loads(caution_path.read_text(encoding="utf-8"))["clases"].items():
            for f in cls["farmacos"]:
                origins.setdefault(f["name"].upper(), []).append(f"caution_list:{ckey}")
    else:
        warnings.append(f"AVISO: no existe {caution_path}; ejecute primero python -m track2.caution_list")
    if class_path.exists():
        for c in json.loads(class_path.read_text(encoding="utf-8")).get("supported_in_two_usable_lines", []):
            for n in c.get("members", []):
                origins.setdefault(n.upper(), []).append("class_analysis:supported_in_two_usable_lines")
    else:
        warnings.append(f"WARNING: {class_path} does not exist (supported_in_two_usable_lines); that source is skipped. Pass the path with --class-analysis")
    if four_path.exists():
        for rname, rd in json.loads(four_path.read_text(encoding="utf-8")).get("readings", {}).items():
            for n in rd.get("one_line_only", []):
                origins.setdefault(n.upper(), []).append(f"four_lines_review2:{rname}:one_line_only")
    else:
        warnings.append(f"WARNING: {four_path} does not exist; that source is skipped")
    if proposals_path is not None and proposals_path.exists():
        for p in json.loads(proposals_path.read_text(encoding="utf-8")).get("proposed", []):
            origins.setdefault(p["drug"].upper(), []).append("proposals_review2:proposed")
    unmatched = sorted(n for n in origins if n not in by_name)
    return {n: o for n, o in origins.items() if n in by_name}, warnings, unmatched


def build(universe: list[dict[str, Any]], origins: dict[str, list[str]], review: dict[str, Any], read_on: str,
          fetch: Callable[[str], tuple[dict[str, Any] | None, list[str]]] = caution_list.fetch_label) -> dict[str, Any]:
    by_name = {d["name"].upper(): d for d in universe}
    rows = []
    for name in sorted(origins):
        lab, failures = fetch(by_name[name]["name"])
        rows.append(row_for(by_name[name], origins[name], lab, failures, review, read_on))
    return {"aviso": NOTICE, "protocolo": "protocol/addendum_review2.md E15 (+E14)", "fecha_lectura": read_on, "grupos": list(age_groups.GROUPS),
            "estados": STATUS_EXPLAINED, "n_farmacos": len(rows), "farmacos": rows,
            "resumen_por_grupo": age_groups.summarise([{"name": r["name"], "age_groups": {g: s for g, s in r["estado_por_grupo"].items()}} for r in rows]),
            "entradas_propuestas_para_config_NO_ESCRITAS": propose_entries(rows),
            "sin_etiqueta_legible": sorted(r["name"] for r in rows if r["etiqueta"] is None)}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--caution-list", type=Path, default=Path("out/caution_list.json"))
    p.add_argument("--class-analysis", type=Path, default=Path("out/class_analysis.json"))
    p.add_argument("--four-lines", type=Path, default=Path("out/four_lines_review2.json"))
    p.add_argument("--proposals", type=Path, default=Path("out/proposals_review2.json"))
    p.add_argument("--out", type=Path, default=Path("out/age_label_table.json"))
    a = p.parse_args(argv)
    universe = json.loads(Path("data/universe.json").read_text(encoding="utf-8"))["drugs"]
    origins, warnings, unmatched = gather(a.caution_list, a.class_analysis, a.four_lines, universe, a.proposals)
    for w in warnings:
        print(w, file=sys.stderr)
    review = load_review()
    setids = {}
    for e in (review.get("entries") or review).values() if isinstance(review.get("entries", review), dict) else []:
        m = re.search(r"setid=([0-9a-f-]{36})", (e.get("source") or "") if isinstance(e, dict) else "")
        if m:
            setids[e["name"].upper()] = m.group(1)

    def fetch(name: str) -> tuple[dict[str, Any] | None, list[str]]:
        """The label the paediatric review read (its set id), so the age table, the quotes and the annex describe one product."""
        if name.upper() in setids:
            lab, failures = caution_list.fetch_by_setid(setids[name.upper()])
            if lab:
                return lab, failures
        return caution_list.fetch_label(name)

    result = build(universe, origins, review, dt.date.today().isoformat(), fetch)
    result["avisos"] = warnings
    result["nombres_sin_correspondencia_en_universo"] = unmatched
    a.out.write_text(json.dumps(result, indent=1, ensure_ascii=False, sort_keys=True), encoding="utf-8", newline="\n")
    print("farmacos:", result["n_farmacos"], "sin etiqueta legible:", result["sin_etiqueta_legible"], "sin correspondencia:", unmatched)
    return 0


if __name__ == "__main__":
    sys.exit(main())
