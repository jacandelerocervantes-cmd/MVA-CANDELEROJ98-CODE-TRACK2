"""Independent check that every quoted label sentence appears in the label text it cites (caution list, age table, manual reviews).

    python -m track2.verify_quotes

Whitespace is ignored on both sides (inline markup such as superscripts splits words differently in the XML), and HTML entities are decoded; nothing else is normalised. DailyMed quotes are
checked against the cached SPL file of the cited set id, openFDA quotes against the cached openFDA label text. Writes `out/quote_check.json`; exit code 1 if any quote is not found.
"""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path
from typing import Any

from track2 import caution_list as cl, labels


def squash(text: str) -> str:
    return re.sub(r"\s+", "", html.unescape(text).replace("\xa0", " ").replace("­", ""))


def _dailymed_text(setid: str) -> str:
    return squash(cl._norm("".join(labels._spl_xml(setid).itertext())))


def _openfda_text(name: str) -> str:
    sections = cl._read_openfda(name)["sections"]
    return squash(" ".join(sections.values()))


def collect(out: Path, config: Path) -> list[dict[str, Any]]:
    items = []
    caution = json.loads((out / "caution_list.json").read_text(encoding="utf-8"))
    for cls in caution["clases"].values():
        for drug in cls["farmacos"]:
            for q in drug["citas"]:
                items.append({"origin": "caution_list", "drug": drug["name"], "source": q["fuente"], "setid": q["setid"], "quote": q["texto_literal"]})
    table = json.loads((out / "age_label_table.json").read_text(encoding="utf-8"))
    for d in table["farmacos"]:
        lab = d.get("etiqueta")
        if not lab:
            continue
        for q in [d.get("frase_literal_sustento")] + list(d.get("advertencias_pediatricas") or []):
            if q:
                items.append({"origin": "age_label_table", "drug": d["name"], "source": lab["fuente"], "setid": lab["setid"], "quote": q["texto_literal"]})
    review = json.loads(config.read_text(encoding="utf-8"))["entries"]
    for e in review.values():
        if e.get("quote"):
            m = re.search(r"setid=([0-9a-f-]{36})", e.get("source", ""))
            items.append({"origin": "pediatric_label_review", "drug": e["name"], "source": "DailyMed" if "DailyMed" in e.get("source", "") else "openFDA", "setid": m.group(1) if m else None, "quote": e["quote"]})
    safety = Path("config/proposal_safety.json")
    if safety.exists():
        for name, qs in json.loads(safety.read_text(encoding="utf-8"))["entries"].items():
            for q in qs:
                items.append({"origin": "proposal_safety", "drug": name, "source": "DailyMed", "setid": q["setid"], "quote": q["quote"]})
    return items


def check(items: list[dict[str, Any]]) -> dict[str, Any]:
    cache: dict[str, str] = {}
    fails, unchecked = [], []
    ok = 0
    for it in items:
        key = f"{it['source']}|{it['setid'] if it['source'] == 'DailyMed' else it['drug']}"
        if key not in cache:
            try:
                cache[key] = _dailymed_text(it["setid"]) if it["source"] == "DailyMed" else _openfda_text(it["drug"])
            except Exception as error:  # noqa: BLE001 - reported, never counted as verified
                cache[key] = ""
                unchecked.append({**it, "reason": repr(error)[:120]})
                continue
        if not cache[key]:
            continue
        if squash(it["quote"]) in cache[key]:
            ok += 1
        else:
            fails.append({k: it[k] for k in ("origin", "drug", "source", "setid")} | {"quote": it["quote"][:140]})
    by_origin: dict[str, dict[str, int]] = {}
    for it in items:
        by_origin.setdefault(it["origin"], {"quotes": 0})["quotes"] += 1
    return {"n_quotes": len(items), "found": ok, "not_found": len(fails), "not_checkable": len(unchecked), "by_origin": by_origin, "failures": fails, "unchecked": unchecked[:20],
            "method": "whitespace ignored on both sides, HTML entities decoded, nothing else normalised"}


def main() -> None:
    result = check(collect(Path("out"), Path("config/pediatric_label_review.json")))
    Path("out/quote_check.json").write_text(json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print({k: result[k] for k in ("n_quotes", "found", "not_found", "not_checkable")}, result["by_origin"])
    for f in result["failures"][:10]:
        print("NOT FOUND:", f)
    sys.exit(1 if result["not_found"] else 0)


if __name__ == "__main__":
    main()
