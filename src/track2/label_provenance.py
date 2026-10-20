"""Provenance of every label read in the manual review of 2026-10-08: set id, label version, effective date, title and marketing status.

    python -m track2.label_provenance

For each entry of `config/pediatric_label_review.json` read on 2026-10-08 the SPL is taken from the cache (DailyMed) to get the version and the effective date, and openFDA `drugsfda` is queried
(through the cache) for the marketing status of the products linked to that set id. A status that cannot be obtained is recorded as `not_found`, never assumed. Writes `out/label_provenance.json`.
"""

from __future__ import annotations

import json
import re
import urllib.parse
from pathlib import Path
from typing import Any

from track2 import labels
from track2.http_cache import get_json

DAY = "2026-10-08"


def spl_facts(setid: str) -> dict[str, Any]:
    root = labels._spl_xml(setid)
    version = root.find(labels.H + "versionNumber")
    effective = root.find(labels.H + "effectiveTime")
    title = root.find(labels.H + "title")
    return {"version": version.get("value") if version is not None else None, "effective": effective.get("value") if effective is not None else None,
            "title": " ".join("".join(title.itertext()).split())[:160] if title is not None else None}


def marketing_status(setid: str) -> dict[str, Any]:
    q = urllib.parse.quote(f'openfda.spl_set_id:"{setid}"')
    try:
        data = get_json(f"https://api.fda.gov/drug/drugsfda.json?search={q}&limit=5")
    except LookupError:
        return {"status": "not_found", "products": []}
    statuses: dict[str, int] = {}
    for result in data.get("results", []):
        for product in result.get("products", []):
            s = product.get("marketing_status") or "unknown"
            statuses[s] = statuses.get(s, 0) + 1
    return {"status": ", ".join(f"{k} ({v})" for k, v in sorted(statuses.items())) or "not_found", "products": sorted(statuses)}


def build(config: Path) -> dict[str, Any]:
    entries = {k: v for k, v in json.loads(config.read_text(encoding="utf-8"))["entries"].items() if DAY in v.get("source", "")}
    out: dict[str, Any] = {}
    for drug_id, e in sorted(entries.items(), key=lambda kv: kv[1]["name"]):
        m = re.search(r"setid=([0-9a-f-]{36})", e["source"])
        setid = m.group(1) if m else None
        row: dict[str, Any] = {"drug_id": drug_id, "setid": setid, "origin": "openFDA" if "openFDA" in e["source"] else "DailyMed"}
        if setid:
            try:
                row.update(spl_facts(setid))
            except Exception as error:  # noqa: BLE001 - recorded, not assumed
                row["spl_error"] = repr(error)[:140]
            row["marketing"] = marketing_status(setid)
        for extra in e.get("also_read") or []:  # further labels read for the same entry (e.g. Neoral and Gengraf for cyclosporine)
            item = {"product": extra["product"], "setid": extra["setid"]}
            try:
                item.update(spl_facts(extra["setid"]))
            except Exception as error:  # noqa: BLE001 - recorded, not assumed
                item["spl_error"] = repr(error)[:140]
            item["marketing"] = marketing_status(extra["setid"])
            row.setdefault("also_read", []).append(item)
        out[e["name"]] = row
    return out


def main() -> None:
    result = build(Path("config/pediatric_label_review.json"))
    n_status = sum(1 for v in result.values() if v.get("marketing", {}).get("status") not in (None, "not_found"))
    payload = {"read_on": DAY, "n_labels": len(result), "n_with_marketing_status": n_status, "labels": result,
               "note": "marketing status comes from openFDA drugsfda for the products linked to the set id; not_found means the lookup returned nothing, not that the product is withdrawn"}
    Path("out/label_provenance.json").write_text(json.dumps(payload, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(len(result), "labels;", n_status, "with a marketing status")
    for name, v in list(result.items())[:6]:
        print(name, v.get("version"), v.get("effective"), v.get("marketing", {}).get("status"))


if __name__ == "__main__":
    main()
