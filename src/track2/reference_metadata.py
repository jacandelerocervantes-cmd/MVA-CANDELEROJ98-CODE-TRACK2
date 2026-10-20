"""Bibliographic data of the papers cited by PMID in the main text of the report, from Europe PMC.

    python -m track2.reference_metadata

Writes `config/reference_metadata.json` (first author, title, journal, year, DOI, as returned by Europe PMC; nothing is typed by hand).
`track2.build_article` reads that file to number the reference list; a PMID missing from it is listed by its identifier only.
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPORT = Path("submissions/candeleroj98_track2_report.md")
OUT = Path("config/reference_metadata.json")
EXTRA = ["26831545", "30002366", "33264411"]  # method papers (network proximity) and the Human Phenotype Ontology paper, cited in the reference list
READ_ON = "2026-10-09"


def main_text_pmids(report: Path) -> list[str]:
    text = report.read_text(encoding="utf-8")
    main = text[: text.find("## Declarations")]
    order: list[str] = []
    for m in re.finditer(r"PMIDs? (\d{6,8})", main):
        if m.group(1) not in order:
            order.append(m.group(1))
    return order


def fetch(pmid: str) -> dict[str, str]:
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": f"EXT_ID:{pmid} AND SRC:MED", "format": "json", "resultType": "core"})
    for attempt in range(6):
        try:
            data = json.loads(urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "ignore"))
            r = data["resultList"]["result"][0]
            authors = r.get("authorString", "") or ""
            first = authors.split(",")[0].strip() if authors else ""
            return {"first_author": first, "more_authors": "," in authors.rstrip("."), "title": re.sub(r"<[^>]+>", "", r.get("title", "")).rstrip("."),
                    "journal": r.get("journalInfo", {}).get("journal", {}).get("isoabbreviation", "") or r.get("journalInfo", {}).get("journal", {}).get("title", "") or r.get("journalTitle", ""), "year": str(r.get("pubYear", "")),
                    "volume": r.get("journalInfo", {}).get("volume", "") or "", "pages": r.get("pageInfo", "") or "", "doi": r.get("doi", "") or ""}
        except Exception:  # noqa: BLE001 - retried, then reported
            time.sleep(3 * (attempt + 1))
    raise SystemExit(f"Europe PMC did not answer for PMID {pmid}")


def main() -> None:
    pmids = EXTRA + [p for p in main_text_pmids(REPORT) if p not in EXTRA]
    entries = {p: fetch(p) for p in pmids}
    OUT.write_text(json.dumps({"source": "Europe PMC (search by EXT_ID)", "read_on": READ_ON, "entries": entries}, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(len(entries), "references written to", OUT)
    missing = [p for p, e in entries.items() if not e["title"] or not e["first_author"]]
    if missing:
        print("incomplete:", missing)
        sys.exit(1)


if __name__ == "__main__":
    main()
