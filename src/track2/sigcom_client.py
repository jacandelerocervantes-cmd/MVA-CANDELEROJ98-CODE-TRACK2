"""Minimal, cached client for the public SigCom LINCS APIs (metadata-api, data-api, S3 consensus signatures).

Endpoints used (all public, no key):
  GET  {META}/entities?filter=...                 gene symbol -> entity UUID (meta.symbol)
  GET  {META}/libraries/{id}/signatures?filter=.. signature UUID -> pert_name / cell_line
  GET  {META}/signatures?filter=...               same, by signature UUID list or by pert_name
  POST {DATA}/enrich/ranktwosided                 two-sided (up/down) reversal ranking of a library
  GET  s3://lincs-dcic/LINCS-sigs-2021/means/...  consensus gene signature (value per gene)
Responses are cached on disk (trimmed, gzip) so a second run needs no network. A failure raises; it is never
read as "no data".
"""

from __future__ import annotations

import gzip
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from track2 import http_cache

META = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
DATA = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"
S3_MEAN = "https://s3.amazonaws.com/lincs-dcic/LINCS-sigs-2021/means"
LIB_CRISPR_CONSENSUS = "98c6a1ef-e3b7-5a96-9913-480acf78a577"  # dataset l1000_mean_xpr
LIB_CHEM_CONSENSUS = "42cd56da-0ad8-5dad-b27c-fe1d135401b2"  # dataset l1000_mean_cp
LIB_CHEM_CELL = "54198d6e-fe17-5ef8-91ac-02b425761653"  # dataset l1000_cp (one signature per cell line/dose/time)
DATASET = {LIB_CHEM_CONSENSUS: "l1000_mean_cp", LIB_CHEM_CELL: "l1000_cp"}  # the data-api takes the dataset name
KEEP = ("uuid", "type", "z-sum", "z-up", "z-down", "p-up", "p-down", "fdr-up", "fdr-down")


class _Keep307(urllib.request.HTTPRedirectHandler):
    """The data-api answers POST with 307 to its cluster host; re-send the same method and body."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[no-untyped-def]
        if code == 307:
            return urllib.request.Request(newurl, data=req.data, headers={k: v for k, v in req.header_items() if k.lower() not in ("host", "content-length")},
                method=req.get_method())
        return super().redirect_request(req, fp, code, msg, headers, newurl)


_OPENER = urllib.request.build_opener(_Keep307)


def _cache(key: str) -> Path:
    return http_cache.CACHE_DIR / f"sigcom_{hashlib.sha256(key.encode()).hexdigest()[:24]}.json.gz"


def _request(url: str, body: dict[str, Any] | None, retries: int = 6, timeout: int = 300) -> bytes:
    data = None if body is None else json.dumps(body, sort_keys=True).encode()
    headers = {"Accept": "application/json", "User-Agent": http_cache.USER_AGENT}
    if data is not None:
        headers["Content-Type"] = "application/json"
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with _OPENER.open(urllib.request.Request(url, data=data, headers=headers), timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code in (400, 404):
                raise LookupError(f"{error.code}: {url}") from error
            last = error
        except Exception as error:  # noqa: BLE001 - transient network failure
            last = error
        time.sleep(3 + 4 * attempt)
    raise RuntimeError(f"SigCom LINCS did not respond after {retries} attempts: {url} ({last!r})")


def cached_call(url: str, body: dict[str, Any] | None = None, keep: tuple[str, ...] | None = None) -> Any:
    """GET (body None) or POST JSON; cache the (optionally trimmed) result under the hash of url+body."""
    path = _cache(url + "|" + json.dumps(body, sort_keys=True))
    if path.exists():
        return json.loads(gzip.decompress(path.read_bytes()))
    payload = json.loads(_request(url, body))
    if keep is not None and isinstance(payload, dict) and "results" in payload:
        payload = {
            "reversers": payload.get("reversers"), "mimickers": payload.get("mimickers"), "maxRank": payload.get("maxRank"),
            "results": [{k: r[k] for k in keep if k in r} for r in payload["results"]],
        }
    http_cache.CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(json.dumps(payload, sort_keys=True).encode(), mtime=0))
    return payload


def meta_find(endpoint: str, flt: dict[str, Any]) -> list[dict[str, Any]]:
    return cached_call(f"{META}/{endpoint}?filter={urllib.parse.quote(json.dumps(flt, sort_keys=True))}")


def gene_ids(symbols: list[str]) -> dict[str, str]:
    """HGNC symbol -> SigCom entity UUID (symbols missing from SigCom are simply absent)."""
    out: dict[str, str] = {}
    for start in range(0, len(symbols), 150):
        part = symbols[start:start + 150]
        rows = meta_find("entities", {"where": {"meta.symbol": {"inq": part}}, "fields": ["id", "meta.symbol"], "limit": 1000})
        out.update({r["meta"]["symbol"]: r["id"] for r in rows})
    return out


def consensus_signature(kind: str, pert_name: str) -> dict[str, float]:
    """Gene -> consensus value of a LINCS consensus signature (kind 'xpr_mean_sigs' = CRISPR KO)."""
    prefix = {"xpr_mean_sigs": "LINCS_L1000_CRISPRKO"}[kind]
    url = f"{S3_MEAN}/{kind}/{prefix}_{pert_name}_mean_sig.tsv.gz"
    path = http_cache.CACHE_DIR / f"sigcom_sig_{kind}_{pert_name}.tsv.gz"
    if not path.exists():
        http_cache.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        path.write_bytes(http_cache.get_bytes(url, accept="*/*"))
    values: dict[str, float] = {}
    for line in gzip.decompress(path.read_bytes()).decode().splitlines()[1:]:
        gene, value = line.split("\t")
        values[gene] = float(value)
    return values


def two_sided(up: list[str], down: list[str], database: str, limit: int) -> list[dict[str, Any]]:
    """Results of enrich/ranktwosided: `limit` mimickers (z-sum>0) and `limit` reversers (z-sum<0)."""
    body = {"up_entities": up, "down_entities": down, "limit": limit, "database": DATASET[database]}
    return cached_call(f"{DATA}/enrich/ranktwosided", body, keep=KEEP)["results"]


def library_signature_names(library: str, fields: list[str], page: int = 2000) -> dict[str, dict[str, Any]]:
    """signature UUID -> selected meta fields for every signature of a library (paged, cached)."""
    out: dict[str, dict[str, Any]] = {}
    skip = 0
    while True:
        rows = cached_call(f"{META}/libraries/{library}/signatures?filter=" + urllib.parse.quote(
            json.dumps({"fields": ["id"] + fields, "limit": page, "skip": skip}, sort_keys=True)))
        out.update({r["id"]: r["meta"] for r in rows})
        if len(rows) < page:
            return out
        skip += page


def signature_meta(ids: list[str], fields: list[str], batch: int = 60) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for start in range(0, len(ids), batch):
        rows = meta_find("signatures", {"where": {"id": {"inq": sorted(ids[start:start + batch])}}, "fields": ["id"] + fields,
                                        "limit": batch})
        out.update({r["id"]: r["meta"] for r in rows})
    return out
