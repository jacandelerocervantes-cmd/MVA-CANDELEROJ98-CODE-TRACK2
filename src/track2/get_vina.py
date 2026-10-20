"""Download the official AutoDock Vina 1.2.7 binary for this platform and check it against the SHA-256 recorded for the Windows build.

    python -m track2.get_vina

Source: https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7 (Apache-2.0). The binary is not stored in the repository (`tools/` is ignored by git).
The recorded SHA-256 (`tools/vina_download_record.json`) pins the Windows build used for the results; for Linux and macOS the file name differs and the hash is
printed so that it can be recorded, because the published results were produced with the Windows build.
"""

from __future__ import annotations

import hashlib
import json
import platform
import stat
import sys
import urllib.request
from pathlib import Path

BASE = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/"
RECORD = Path("tools/vina_download_record.json")


def asset() -> tuple[str, str]:
    system, machine = platform.system(), platform.machine().lower()
    if system == "Windows":
        return "vina_1.2.7_win.exe", "vina_1.2.7_win.exe"
    if system == "Darwin":
        return ("vina_1.2.7_mac_aarch64" if "arm" in machine or "aarch64" in machine else "vina_1.2.7_mac_x86_64"), "vina"
    return ("vina_1.2.7_linux_aarch64" if "aarch64" in machine or "arm" in machine else "vina_1.2.7_linux_x86_64"), "vina"


def main() -> None:
    name, local = asset()
    dest = Path("tools") / local
    dest.parent.mkdir(exist_ok=True)
    data = urllib.request.urlopen(urllib.request.Request(BASE + name, headers={"User-Agent": "hselith-track2/0.1"}), timeout=120).read()  # noqa: S310 - fixed https URL
    digest = hashlib.sha256(data).hexdigest()
    dest.write_bytes(data)
    dest.chmod(dest.stat().st_mode | stat.S_IEXEC)
    record = json.loads(RECORD.read_text(encoding="utf-8")) if RECORD.exists() else {}
    if name == record.get("file"):
        if digest != record["sha256"]:
            dest.unlink()
            sys.exit(f"SHA-256 mismatch for {name}: expected {record['sha256']}, got {digest}. File deleted.")
        print("verified against the recorded SHA-256:", digest)
    else:
        print(f"downloaded {name}; SHA-256 {digest} (the recorded results used the Windows build)")


if __name__ == "__main__":
    main()
