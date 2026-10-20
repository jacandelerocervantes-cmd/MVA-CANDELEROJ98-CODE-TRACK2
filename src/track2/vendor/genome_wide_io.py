"""Access to public sources for the genome-wide search, with a local disk cache.

- Reference CDS sequences per transcript: Ensembl REST (POST /sequence/id, type=cds).
- Population frequencies per variant: gnomAD v4 API (dataset gnomad_r4, fixed version), queried
  directly per allele; VEP is not used for frequencies because it does not match indels well.
- Ensembl VEP consequence terms: only to VALIDATE the local calculation (consequence_calculator);
  they never decide any label or score.
Results are stored as JSON; with the cache present, the pipeline needs no network.
"""

from __future__ import annotations

import json
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

ENSEMBL_REST = "https://rest.ensembl.org"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
_GNOMAD_BATCH = 25
_CDS_BATCH = 50
_VEP_BATCH = 200


def _post_json(url: str, payload: dict[str, Any], retries: int = 8, timeout: int = 180) -> Any:
    request = urllib.request.Request(
        url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json", "Accept": "application/json"}
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 - fixed Ensembl URL
                return json.load(response)
        except Exception:  # noqa: BLE001 - retry on any transient network failure
            time.sleep(5 + 10 * attempt)
    raise RuntimeError(f"The service did not respond after {retries} attempts: {url}")


def _load_cache(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def fetch_cds_sequences(transcript_ids: list[str], cache_path: Path, workers: int = 6) -> dict[str, str]:
    """CDS sequence (includes the stop codon) per transcript; uses and fills in the cache."""
    cache: dict[str, str] = _load_cache(cache_path)
    pending = [t for t in transcript_ids if t not in cache]
    lock = threading.Lock()

    def work(chunk: list[str]) -> None:
        reply = _post_json(f"{ENSEMBL_REST}/sequence/id?type=cds;multiple_sequences=1", {"ids": [t.split(".")[0] for t in chunk]})
        by_query = {r["query"]: r["seq"] for r in reply}
        with lock:
            for tid in chunk:
                if tid.split(".")[0] in by_query:
                    cache[tid] = by_query[tid.split(".")[0]]

    chunks = [pending[i:i + _CDS_BATCH] for i in range(0, len(pending), _CDS_BATCH)]
    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(work, chunks))
    if pending:  # written only if there were new queries: with a complete cache the disk is not touched
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, sort_keys=True), encoding="utf-8")
    return {t: cache[t] for t in transcript_ids if t in cache}


def validate_cds_lengths(transcripts: dict[str, dict[str, Any]], cds_sequences: dict[str, str]) -> list[str]:
    """Transcripts whose downloaded CDS does not have the same length as CDS+stop from the GTF."""
    bad = []
    for tid, seq in cds_sequences.items():
        t = transcripts[tid]
        expected = sum(e - s + 1 for s, e, _ in t["cds"]) + sum(e - s + 1 for s, e in t["stop"])
        if expected != len(seq):
            bad.append(tid)
    return bad


def fetch_vep_terms(variant_keys: list[str], cache_path: Path, workers: int = 5) -> dict[str, dict[str, Any]]:
    """VEP consequence terms (MANE/canonical transcript) per variant, only to validate the local calculation."""
    cache: dict[str, dict[str, Any]] = _load_cache(cache_path)
    pending = [k for k in variant_keys if k not in cache]
    lock = threading.Lock()

    def work(chunk: list[str]) -> None:
        lines = []
        for key in chunk:
            chrom, pos, ref, alt = key.split(":")
            lines.append(f"{chrom.replace('chr', '')} {pos} . {ref} {alt} . . .")
        reply = _post_json(f"{ENSEMBL_REST}/vep/human/region?canonical=1;mane=1", {"variants": lines})
        found: dict[str, dict[str, Any]] = {}
        for item in reply:
            f = item["input"].split()
            terms = sorted({
                term for tc in item.get("transcript_consequences", [])
                if tc.get("mane_select") or tc.get("canonical") for term in tc["consequence_terms"]
            })
            found[f"chr{f[0]}:{f[1]}:{f[3]}:{f[4]}"] = {"vep_terms": terms}
        with lock:
            for key in chunk:
                cache[key] = found.get(key, {"vep_terms": []})

    chunks = [pending[i:i + _VEP_BATCH] for i in range(0, len(pending), _VEP_BATCH)]
    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(work, chunks))
    if pending:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, sort_keys=True), encoding="utf-8")
    return {k: cache[k] for k in variant_keys if k in cache}


def gnomad_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """gnomAD ID (chrom-pos-ref-alt) in minimal representation: common suffix and then prefix, leaving 1 base."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return f"{chrom.replace('chr', '')}-{pos}-{ref}-{alt}"


def _gnomad_query(ids: list[str]) -> str:
    fields = "exome { ac an af homozygote_count } genome { ac an af homozygote_count }"
    body = " ".join(f'v{i}: variant(variantId: "{vid}", dataset: {GNOMAD_DATASET}) {{ {fields} }}' for i, vid in enumerate(ids))
    return "query {" + body + "}"


def _summarize_gnomad(record: dict[str, Any] | None) -> dict[str, Any]:
    """Per-allele summary: present in gnomAD, maximum AF across exomes and genomes, and observed homozygotes."""
    if not record:
        return {"present": False, "af": None, "homozygotes": 0}
    parts = [record.get("exome"), record.get("genome")]
    afs = [p["af"] for p in parts if p and p.get("af") is not None]
    return {"present": any(parts), "af": max(afs) if afs else 0.0,
            "homozygotes": sum(int(p.get("homozygote_count") or 0) for p in parts if p)}


def fetch_gnomad_frequencies(variant_keys: list[str], cache_path: Path, workers: int = 1) -> dict[str, dict[str, Any]]:
    """gnomAD v4 frequency per allele (maximum across exomes and genomes) and number of homozygotes.

    `present` = False means that gnomAD v4 does not contain the variant (absence confirmed by the API), not that
    the query failed: a network or API failure raises an error and is not interpreted as absence.
    """
    cache: dict[str, dict[str, Any]] = _load_cache(cache_path)
    pending = [k for k in variant_keys if k not in cache]
    lock = threading.Lock()

    def work(chunk: list[str]) -> None:
        ids = []
        for key in chunk:
            chrom, pos, ref, alt = key.split(":")
            ids.append(gnomad_variant_id(chrom, int(pos), ref, alt))
        for attempt in range(10):  # "Service overloaded" is transient: it is retried, never interpreted as absence
            reply = _post_json(GNOMAD_API, {"query": _gnomad_query(ids)})
            other = [e for e in reply.get("errors") or [] if "variant not found" not in e.get("message", "").lower()]
            if not other or not all("overloaded" in e.get("message", "").lower() for e in other):
                break
            time.sleep(10 + 15 * attempt)
        data = reply.get("data")
        if other or data is None:
            raise RuntimeError(f"gnomAD returned errors other than 'Variant not found': {other[:2]}")
        with lock:
            for i, key in enumerate(chunk):
                cache[key] = _summarize_gnomad(data.get(f"v{i}"))  # null alias = variant not found in gnomAD v4
        time.sleep(1.0)

    chunks = [pending[i:i + _GNOMAD_BATCH] for i in range(0, len(pending), _GNOMAD_BATCH)]
    try:
        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(work, chunks))
    finally:  # what has already been queried is kept even if a later query fails
        if pending:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(cache, sort_keys=True), encoding="utf-8")
    return {k: cache[k] for k in variant_keys if k in cache}
