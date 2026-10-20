"""Check that every drug of the study is FDA-approved and marketed (not experimental), against Drugs@FDA (openFDA).

    python -m track2.fda_approval_check

For each drug of `out/age_label_table.json` it looks for the generic name in `openfda.generic_name`, then in `openfda.substance_name`, then in the active ingredients of the products. A drug is confirmed when
a record with an NDA, BLA or ANDA has at least one product with marketing status `Prescription` or `Over-the-counter`. Writes `out/fda_approval_check.json` (application numbers and statuses); exit code 1 if a drug
is not confirmed. A network service is used, so the file is listed among the network-dependent outputs.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

BASE = "https://api.fda.gov/drug/drugsfda.json"
UA = {"User-Agent": "hselith-track2/0.1 (research)"}
OK_STATUS = {"Prescription", "Over-the-counter"}
FDA_NAME = {"PHENYLBUTANOIC ACID": "sodium phenylbutyrate"}  # the name of the drug universe is the parent acid; the approved product is the sodium salt


def _get(search: str) -> list[dict[str, Any]]:
    url = f"{BASE}?search={urllib.parse.quote(search)}&limit=100"
    last: Exception | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=40) as r:  # noqa: S310 - fixed https endpoint
                time.sleep(0.3)  # the service rate-limits
                return json.load(r).get("results", [])
        except urllib.error.HTTPError as error:
            if error.code == 404:  # openFDA answers 404 for "no match"
                return []
            last = error
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last = error
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Drugs@FDA did not answer for {search!r}: {last!r}")


def check(name: str) -> dict[str, Any]:
    word = FDA_NAME.get(name, name).lower()
    for field in ("openfda.generic_name", "openfda.substance_name", "products.active_ingredients.name"):
        results = _get(f'{field}:"{word}"') or _get(f"{field}:{word.split()[0]}*")
        if not results:
            continue
        apps = sorted({r.get("application_number", "") for r in results if r.get("application_number", "").startswith(("NDA", "BLA", "ANDA"))})
        statuses = sorted({p.get("marketing_status", "") for r in results for p in r.get("products", [])})
        brands = sorted({p.get("brand_name", "") for r in results for p in r.get("products", []) if p.get("marketing_status") in OK_STATUS})[:4]
        if apps and OK_STATUS & set(statuses):
            return {"confirmed": True, "found_by": field, "n_records": len(results), "applications": apps[:6], "marketing_statuses": statuses, "example_brands": brands}
    return {"confirmed": False, "found_by": None, "n_records": 0, "applications": [], "marketing_statuses": [], "example_brands": []}


def main() -> None:
    table = json.loads(Path("out/age_label_table.json").read_text(encoding="utf-8"))["farmacos"]
    names = sorted({d["name"] for d in table})
    rows = {n: check(n) for n in names}
    bad = [n for n, r in rows.items() if not r["confirmed"]]
    out = {"source": "Drugs@FDA through openFDA (https://api.fda.gov/drug/drugsfda.json)", "rule": "an NDA, BLA or ANDA record with a product marketed as Prescription or Over-the-counter",
           "n_drugs": len(rows), "n_confirmed": len(rows) - len(bad), "not_confirmed": bad, "drugs": rows}
    Path("out/fda_approval_check.json").write_text(json.dumps(out, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(out["n_confirmed"], "of", out["n_drugs"], "confirmed in Drugs@FDA; not confirmed:", bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
