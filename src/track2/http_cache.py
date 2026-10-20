"""HTTP downloads with disk cache, retries and sorted keys (stable JSON).

Every query to an external service goes through here: the response is stored in `data/cache/` under the
hash of the URL, so a second run needs no network and produces the same bytes.
A network or API error raises an exception; it is never interpreted as "missing data".
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

CACHE_DIR = Path("data/cache")
USER_AGENT = "hselith-track2/0.1 (research; contact via repository)"


def _cache_path(url: str, suffix: str) -> Path:
    digest = hashlib.sha256(url.encode()).hexdigest()[:24]
    return CACHE_DIR / f"{digest}{suffix}"


def get_bytes(url: str, retries: int = 6, timeout: int = 120, accept: str = "application/json") -> bytes:
    request = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": USER_AGENT})
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 - fixed URL of each module
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code == 404:  # the resource does not exist: retrying cannot help
                raise LookupError(f"Not found: {url}") from error
            last = error
            time.sleep(3 + 5 * attempt)
        except Exception as error:  # noqa: BLE001 - retry on transient failure
            last = error
            time.sleep(3 + 5 * attempt)
    raise RuntimeError(f"The service did not respond after {retries} attempts: {url} ({last!r})")


def get_json(url: str, use_cache: bool = True, accept: str = "application/json") -> Any:
    path = _cache_path(url, ".json")
    if use_cache and path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    payload = json.loads(get_bytes(url, accept=accept))
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8", newline="\n")
    return payload


def get_file(url: str, name: str, accept: str = "*/*") -> Path:
    """Download a large file only once (under a stable name) and return its path."""
    path = CACHE_DIR / name
    if not path.exists():
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        path.write_bytes(get_bytes(url, timeout=600, accept=accept))
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()
