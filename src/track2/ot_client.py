"""Open Targets Platform GraphQL client (read-only, disk cached) used by the multi-disease benchmark.

Every request is cached under `data/cache/` by the hash of (query, variables), like `http_cache`.
A network or API error raises; it is never read as "no data".
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from typing import Any

from track2 import http_cache

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
APPROVED = "APPROVAL"
# Drug-evidence datasource excluded from the association score (see protocol/addendum_benchmark.md).
EXCLUDED_DATASOURCES = [{"id": "clinical_precedence", "weight": 0, "propagate": False}]
BATCH = 40
PAGE = 1000

_DRUGS_Q = "query($i:[String!]!){drugs(chemblIds:$i){id indications{count rows{maxClinicalStage disease{id name}}}}}"
_TARGETS_Q = (
    "query($d:String!,$ds:[DatasourceSettingsInput!],$p:Pagination){disease(efoId:$d){"
    "associatedTargets(datasources:$ds,page:$p){count rows{score target{id approvedSymbol}}}}}"
)


def post_graphql(query: str, variables: dict[str, Any], retries: int = 8) -> dict[str, Any]:
    key = hashlib.sha256(json.dumps([query, variables], sort_keys=True).encode()).hexdigest()[:24]
    path = http_cache.CACHE_DIR / f"ot_{key}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    body = json.dumps({"query": query, "variables": variables}).encode()
    headers = {"Content-Type": "application/json", "User-Agent": http_cache.USER_AGENT}
    last: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(OT_URL, data=body, headers=headers)
            with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310 - fixed URL
                payload = json.loads(response.read())
            break
        except urllib.error.HTTPError as error:
            if error.code == 400:
                raise RuntimeError(f"Open Targets rejected the query: {error.read()[:500]!r}") from error
            last = error
            time.sleep(2 + 3 * attempt)
        except Exception as error:  # noqa: BLE001 - transient network / TLS failure
            last = error
            time.sleep(2 + 3 * attempt)
    else:
        raise RuntimeError(f"Open Targets did not respond after {retries} attempts ({last!r})")
    if payload.get("errors"):
        raise RuntimeError(f"Open Targets returned errors: {payload['errors']!r}")
    http_cache.CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8", newline="\n")
    return payload


def approved_indications(drug_ids: list[str]) -> tuple[dict[str, set[str]], dict[str, str]]:
    """disease id -> set of drug ids with an APPROVAL indication; disease id -> name. Direct indications only."""
    by_disease: dict[str, set[str]] = {}
    names: dict[str, str] = {}
    for start in range(0, len(drug_ids), BATCH):
        payload = post_graphql(_DRUGS_Q, {"i": drug_ids[start:start + BATCH]})
        for drug in payload["data"]["drugs"]:
            block = drug["indications"]
            if block["count"] != len(block["rows"]):
                raise RuntimeError(f"Incomplete indications for {drug['id']}")
            for row in block["rows"]:
                if row["maxClinicalStage"] == APPROVED:
                    by_disease.setdefault(row["disease"]["id"], set()).add(drug["id"])
                    names[row["disease"]["id"]] = row["disease"]["name"]
    return by_disease, names


def associated_symbols(disease_id: str, min_score: float) -> list[str]:
    """approvedSymbol of every target with overall score >= min_score (clinical_precedence excluded)."""
    symbols: list[str] = []
    index = 0
    while True:
        payload = post_graphql(_TARGETS_Q, {"d": disease_id, "ds": EXCLUDED_DATASOURCES, "p": {"index": index, "size": PAGE}})
        block = (payload["data"]["disease"] or {}).get("associatedTargets") or {"rows": [], "count": 0}
        rows = block["rows"]
        symbols.extend(r["target"]["approvedSymbol"] for r in rows if r["score"] >= min_score)
        if not rows or rows[-1]["score"] < min_score or (index + 1) * PAGE >= block["count"]:
            return sorted(set(symbols))
        index += 1
