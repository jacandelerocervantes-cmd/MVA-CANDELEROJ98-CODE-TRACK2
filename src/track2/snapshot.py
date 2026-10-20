"""Data manifest: checksums of the third-party inputs, so a later run can tell whether the data drifted.

    python -m track2.snapshot --out out

Writes `out/data_manifest.json`. Small cached API responses (hundreds of JSON files) and DailyMed label
files are summarised as count + total bytes + one SHA-256 over the sorted (name, sha256) pairs; large named
downloads are listed one by one. Only `data/cache/` and `data/universe.json` are read: `data/private/`
(patient-derived) is never touched and nothing patient-derived is written.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

from track2.http_cache import sha256_file

API_JSON = re.compile(r"^[0-9a-f]{24}\.json$")
STRING_FILES = ("9606.protein.aliases.v12.0.txt.gz", "9606.protein.physical.links.v12.0.txt.gz")
VERSION_KEYS = ("chembl_version", "string_version", "reactome_version")


def group_digest(pairs: list[tuple[str, str]]) -> str:
    """One SHA-256 over the sorted (name, sha256) pairs."""
    digest = hashlib.sha256()
    for name, sha in sorted(pairs):
        digest.update(f"{name}\t{sha}\n".encode())
    return digest.hexdigest()


def summarise_group(files: list[Path]) -> dict[str, Any]:
    pairs = [(f.name, sha256_file(f)) for f in files]
    return {"count": len(files), "total_bytes": sum(f.stat().st_size for f in files), "sha256_of_sorted_name_sha256_pairs": group_digest(pairs)}


def classify(cache_dir: Path) -> dict[str, list[Path]]:
    groups: dict[str, list[Path]] = {"api_json_responses": [], "dailymed_spl_xml": [], "large_downloads": []}
    if not cache_dir.is_dir():
        return groups
    for path in sorted(p for p in cache_dir.iterdir() if p.is_file()):
        if API_JSON.match(path.name):
            groups["api_json_responses"].append(path)
        elif path.name.startswith("spl_") and path.suffix == ".xml":
            groups["dailymed_spl_xml"].append(path)
        else:
            groups["large_downloads"].append(path)
    return groups


def recorded_versions(out: Path) -> dict[str, Any]:
    """Versions recorded in out/manifest_*.json; a disagreement between mechanisms is reported, not hidden."""
    seen: dict[str, set[str]] = {k: set() for k in VERSION_KEYS}
    files = sorted(out.glob("manifest_*.json"))
    for path in files:
        data = json.loads(path.read_text(encoding="utf-8"))
        for key in VERSION_KEYS:
            if key in data:
                seen[key].add(str(data[key]))
    result: dict[str, Any] = {k.replace("_version", ""): sorted(v) for k, v in seen.items()}
    result["manifests_read"] = [p.name for p in files]
    return result


def build_manifest(data_dir: Path, out: Path) -> dict[str, Any]:
    cache = data_dir / "cache"
    groups = classify(cache)
    explicit = {}
    for name in STRING_FILES:
        path = cache / name
        explicit[name] = {"sha256": sha256_file(path), "bytes": path.stat().st_size} if path.exists() else None
    universe = data_dir / "universe.json"
    large = {p.name: {"sha256": sha256_file(p), "bytes": p.stat().st_size} for p in groups["large_downloads"]}
    return {
        "note": "Third-party data only. No patient-derived file is read or listed (data/private is excluded).",
        "universe_json": {"sha256": sha256_file(universe), "bytes": universe.stat().st_size} if universe.exists() else None,
        "string_files": explicit,
        "large_downloads": large,
        "api_json_responses": summarise_group(groups["api_json_responses"]),
        "dailymed_spl_xml": summarise_group(groups["dailymed_spl_xml"]),
        "recorded_versions": recorded_versions(out),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--data", type=Path, default=Path("data"))
    args = parser.parse_args()
    manifest = build_manifest(args.data, args.out)
    target = args.out / "data_manifest.json"
    args.out.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
