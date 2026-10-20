"""Download the DepMap/PRISM files used by `protocol/addendum_prism.md` from figshare and check their MD5.

    python -m track2.depmap_download            # downloads the files that are missing (about 540 MB in total)
    python -m track2.depmap_download --list     # shows what would be downloaded

Source: figshare articles 25917643 (PRISM Repurposing Public 24Q2) and 25880521 (DepMap Public 24Q2), both CC BY 4.0.
Files go to `data/depmap/` (ignored by git). A failed MD5 check deletes the file and stops.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

DEST = Path("data/depmap")
WANTED = {
    25917643: ["README.txt", "Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv", "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv",
               "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv", "Repurposing_Public_24Q2_Treatment_Meta_Data.csv"],
    25880521: ["Model.csv", "OmicsExpressionProteinCodingGenesTPMLogp1.csv"],
}
UA = {"User-Agent": "hselith-track2/0.1 (research)"}


def article(aid: int) -> dict:
    return json.loads(urllib.request.urlopen(urllib.request.Request(f"https://api.figshare.com/v2/articles/{aid}", headers=UA), timeout=60).read())


def md5(path: Path) -> str:
    h = hashlib.md5()  # noqa: S324 - integrity check against the publisher's checksum, not a security use
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    DEST.mkdir(parents=True, exist_ok=True)
    record = {}
    for aid, names in WANTED.items():
        meta = article(aid)
        by_name = {f["name"]: f for f in meta["files"]}
        for name in names:
            f = by_name[name]
            target = DEST / name
            record[name] = {"figshare_article": aid, "title": meta["title"], "license": meta["license"]["name"], "md5": f["computed_md5"], "size_mb": round(f["size"] / 1e6, 1)}
            if args.list:
                print(f"{name}: {record[name]['size_mb']} MB ({meta['title']}, {meta['license']['name']})")
                continue
            if target.exists() and md5(target) == f["computed_md5"]:
                print("already present and verified:", name)
                continue
            print("downloading", name, record[name]["size_mb"], "MB", flush=True)
            urllib.request.urlretrieve(urllib.request.Request(f["download_url"], headers=UA).full_url, target)  # noqa: S310 - fixed https figshare URL
            if md5(target) != f["computed_md5"]:
                target.unlink()
                sys.exit(f"MD5 mismatch for {name}; file deleted.")
    if not args.list:
        (DEST / "download_record.json").write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        print("all files verified")


if __name__ == "__main__":
    main()
