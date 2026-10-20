"""Delete the local data of this repository (dry run by default).

    python -m track2.delete_data              # lists what would be deleted
    python -m track2.delete_data --execute    # deletes it and writes out/deletion_attestation.json

Covers `data/` (patient-derived detail in `data/private/`, third-party downloads and caches). It does not touch
source code, protocol, report or the aggregate outputs in `out/`. Copies outside the repository must be deleted
by their holder.
"""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

TARGETS = [Path("data")]
ATTESTATION = Path("out/deletion_attestation.json")


def inventory(targets: list[Path]) -> list[dict[str, Any]]:
    """Existing paths with file counts and bytes (names of top-level entries only, no content)."""
    rows = []
    for target in targets:
        if not target.exists():
            continue
        files = [p for p in target.rglob("*") if p.is_file()] if target.is_dir() else [target]
        rows.append({"path": str(target), "files": len(files), "bytes": sum(p.stat().st_size for p in files)})
    return rows


def delete(targets: list[Path], execute: bool, attestation: Path = ATTESTATION) -> dict[str, Any]:
    rows = inventory(targets)
    report = {"executed": execute, "when_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "targets": rows}
    if execute:
        for row in rows:
            path = Path(row["path"])
            shutil.rmtree(path) if path.is_dir() else path.unlink()
        attestation.parent.mkdir(parents=True, exist_ok=True)
        attestation.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="actually delete (the default is a dry run)")
    args = parser.parse_args()
    report = delete(TARGETS, args.execute)
    verb = "Deleted" if args.execute else "Would delete (dry run)"
    for row in report["targets"]:
        print(f"{verb}: {row['path']} ({row['files']} files, {row['bytes']:,} bytes)")
    if not report["targets"]:
        print("Nothing to delete.")


if __name__ == "__main__":
    main()
