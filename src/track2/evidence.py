"""Retrieved (not curated) literature and trial context per candidate drug.

PubMed (NCBI E-utilities) and ClinicalTrials.gov API v2 are queried by keyword; every returned PMID is
re-checked through esummary. Presence of a paper or trial is NOT support: this is context for human review.
All HTTP goes through `get_json` (disk cache). A failed call is recorded as `available: false` with the reason.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
from datetime import date
from pathlib import Path
from typing import Any

from track2.http_cache import get_json

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CT_API = "https://clinicaltrials.gov/api/v2/studies"
PAUSE_SECONDS = 0.4  # NCBI allows 3 requests/s without an API key
CAVEAT = (
    "Retrieved by keyword, not curated; the presence of a paper or trial does not mean support. "
    "Hypothesis context for human review, not medical advice."
)


def clean_name(name: str, salt_words: list[str]) -> str:
    """Drug name without salt/form words and without anything after a comma (same convention as toxicity.py)."""
    pattern = r"\b(" + "|".join(re.escape(word) for word in salt_words) + r")\b"
    return re.sub(r"\s+", " ", re.sub(pattern, "", name.upper()).split(",")[0]).strip().lower()


def _get(url: str) -> Any:
    time.sleep(PAUSE_SECONDS)
    return get_json(url)


def _failure(error: Exception) -> dict[str, Any]:
    return {"available": False, "reason": f"{type(error).__name__}: {error}"[:300]}


def pubmed_group(name: str, terms: str, listed: int) -> dict[str, Any]:
    query = f'"{name}" AND ({terms})'
    try:
        search = _get(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax={listed}"
                      f"&term={urllib.parse.quote(query)}")["esearchresult"]
        count, ids = int(search["count"]), list(search["idlist"])
        records: list[dict[str, Any]] = []
        if ids:
            summary = _get(f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(ids)}")["result"]
            for pmid in ids:
                entry = summary.get(pmid)
                if not entry or "error" in entry or not entry.get("title"):
                    records.append({"pmid": pmid, "title": None, "verified": False, "year": None})
                    continue
                year = re.match(r"\d{4}", entry.get("pubdate", ""))
                records.append({"pmid": pmid, "title": entry["title"], "verified": True,
                                "year": int(year.group()) if year else None})
    except (LookupError, RuntimeError, KeyError, ValueError) as error:
        return {"query": query, **_failure(error)}
    return {"available": True, "count": count, "query": query, "records": records,
            "verified_count": sum(1 for record in records if record["verified"])}


def _study(raw: dict[str, Any]) -> dict[str, Any]:
    protocol = raw.get("protocolSection", {})
    ident = protocol.get("identificationModule", {})
    return {
        "minimum_age": protocol.get("eligibilityModule", {}).get("minimumAge"),
        "nct_id": ident.get("nctId"),
        "phase": protocol.get("designModule", {}).get("phases"),
        "status": protocol.get("statusModule", {}).get("overallStatus"),
        "title": ident.get("briefTitle"),
    }


def trial_query(name: str, condition: str | None, advanced: str | None, config: dict[str, Any]) -> dict[str, Any]:
    params = {"query.intr": name, "pageSize": str(config["page_size"]), "countTotal": "true"}
    if condition:
        params["query.cond"] = condition
    if advanced:
        params["filter.advanced"] = advanced
    url = f"{CT_API}?{urllib.parse.urlencode(params, quote_via=urllib.parse.quote)}"
    try:
        payload = _get(url)
        studies = [_study(raw) for raw in payload.get("studies", [])[: config["studies_listed"]]]
        return {"available": True, "studies": studies, "total_count": int(payload["totalCount"])}
    except (LookupError, RuntimeError, KeyError, ValueError) as error:
        return _failure(error)


def drug_evidence(drug: dict[str, Any], config: dict[str, Any], today: str) -> dict[str, Any]:
    name = clean_name(drug["name"], config["name_cleaning"]["salt_words"])
    pubmed_cfg, ct_cfg = config["pubmed"], config["clinicaltrials"]
    pubmed = {group: pubmed_group(name, terms, pubmed_cfg["records_listed"])
              for group, terms in pubmed_cfg["term_groups"].items()}
    trials = {key: trial_query(name, condition, None, ct_cfg) for key, condition in ct_cfg["conditions"].items()}
    trials["pediatric_population"] = trial_query(name, None, ct_cfg["pediatric_filter"], ct_cfg)
    return {"caveat": CAVEAT, "name_queried": name, "pubmed": pubmed, "retrieved_on": today, "trials": trials}


def run(candidates_path: Path, out_dir: Path, config_path: Path) -> Path:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    drugs = json.loads(candidates_path.read_text(encoding="utf-8"))
    today = date.today().isoformat()
    result: dict[str, Any] = {}
    for index, drug in enumerate(sorted(drugs, key=lambda item: item["drug_id"]), 1):
        result[drug["drug_id"]] = drug_evidence(drug, config, today)
        print(f"[{index}/{len(drugs)}] {drug['drug_id']} {drug['name']}", flush=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "evidence_retrieval.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=Path("out/candidates.json"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--config", type=Path, default=Path("config/evidence_queries.json"))
    args = parser.parse_args()
    print(f"Written: {run(args.candidates, args.out, args.config)}")


if __name__ == "__main__":
    main()
