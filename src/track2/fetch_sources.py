"""Download the sources of the four-line analysis (`protocol/addendum_four_lines.md`), one at a time, and record each file's source, size and SHA-256.

    python -m track2.fetch_sources --list                    # what exists, with sizes where the publisher states them (reads only, downloads nothing)
    python -m track2.fetch_sources --only hpo_annotations    # download one source
    python -m track2.fetch_sources --only primekg mgi        # or several

Nothing is downloaded unless a source is named with --only. Files go to `data/sources/<source>/` (ignored by git) and the record to `data/sources/record.json`.
Intended to run on the Linux virtual machine. Data are used for non-commercial research; derived tables of any source whose licence restricts redistribution are not published.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path
from typing import Any

DEST = Path("data/sources")
UA = {"User-Agent": "hselith-track2/0.1 (research)"}
HPO_RELEASE = "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/"
MGI = "https://www.informatics.jax.org/downloads/reports/"
DATAVERSE = "https://dataverse.harvard.edu/api"
PRIMEKG_DOI = "doi:10.7910/DVN/IXA7BM"
FIGSHARE_DEPMAP = 25880521

SOURCES = {
    "hpo_annotations": {"licence": "HPO licence, https://hpo.jax.org/app/license", "files": {"phenotype.hpoa": HPO_RELEASE + "phenotype.hpoa"}},
    "hpo_ontology": {"licence": "HPO licence, https://hpo.jax.org/app/license",
                     "files": {"hp.obo": HPO_RELEASE + "hp.obo", "genes_to_phenotype.txt": HPO_RELEASE + "genes_to_phenotype.txt"}},
    "mondo": {"licence": "CC BY 4.0 (monarch-initiative/mondo; 'latest' release, so the SHA-256 recorded identifies the file used)",
              "files": {"mondo.obo": "https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo"}},
    "mgi": {"licence": "MGI data use terms, https://www.informatics.jax.org/mgihome/other/copyright.shtml",
            "files": {"MGI_GenePheno.rpt": MGI + "MGI_GenePheno.rpt", "HMD_HumanPhenotype.rpt": MGI + "HMD_HumanPhenotype.rpt"}},
    "omnipath": {"licence": "per resource, https://omnipathdb.org/info",
                 "files": {"omnipath_interactions.tsv": "https://omnipathdb.org/interactions?datasets=omnipath,pathwayextra,kinaseextra&genesymbols=yes&fields=sources,references,curation_effort&format=tsv"}},
    "primekg": {"licence": "CC0 1.0 (Harvard Dataverse record " + PRIMEKG_DOI + ")", "dataverse": ["nodes.tab", "edges.csv"]},
    "depmap_crispr": {"licence": "CC BY 4.0 (figshare article 25880521)", "figshare": [FIGSHARE_DEPMAP, ["CRISPRGeneEffect.csv"]]},
    # Covariates of line C (addendum 11, E6): TP53 status (damaging / hotspot matrices), MSI and chromosomal-instability signatures, and the model-profile maps.
    "depmap_covariates": {"licence": "CC BY 4.0 (figshare article 25880521)",
                          "figshare": [FIGSHARE_DEPMAP, ["OmicsSomaticMutationsMatrixDamaging.csv", "OmicsSomaticMutationsMatrixHotspot.csv", "OmicsSignatures.csv",
                                                         "OmicsDefaultModelProfiles.csv", "OmicsProfiles.csv"]]},
}


def _get(url: str, timeout: int = 120) -> bytes:
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout).read()  # noqa: S310 - fixed https URLs listed above


def _resolve(source: dict[str, Any]) -> dict[str, tuple[str, int | None, str | None]]:
    """name -> (url, size in bytes if the publisher states it, md5 if stated)."""
    out: dict[str, tuple[str, int | None, str | None]] = {}
    for name, url in source.get("files", {}).items():
        out[name] = (url, None, None)
    if "dataverse" in source:
        meta = json.loads(_get(f"{DATAVERSE}/datasets/:persistentId/?persistentId={PRIMEKG_DOI}"))["data"]["latestVersion"]["files"]
        by_name = {f["dataFile"]["filename"]: f["dataFile"] for f in meta}
        for name in source["dataverse"]:
            f = by_name[name]
            out[name] = (f"{DATAVERSE}/access/datafile/{f['id']}", f["filesize"], f.get("md5"))
    if "figshare" in source:
        aid, names = source["figshare"]
        files = {f["name"]: f for f in json.loads(_get(f"https://api.figshare.com/v2/articles/{aid}"))["files"]}
        for name in names:
            out[name] = (files[name]["download_url"], files[name]["size"], files[name]["computed_md5"])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--only", nargs="*", choices=sorted(SOURCES), default=[])
    args = parser.parse_args()
    if not args.list and not args.only:
        parser.error("nothing to do: use --list, or name the sources to download with --only")
    record_path = DEST / "record.json"
    record = json.loads(record_path.read_text(encoding="utf-8")) if record_path.exists() else {}
    for key in (sorted(SOURCES) if args.list else args.only):
        src = SOURCES[key]
        files = _resolve(src)
        for name, (url, size, md5) in files.items():
            if args.list:
                print(f"{key}/{name}: {round(size / 1e6, 1) if size else 'size not stated'} MB; licence: {src['licence']}")
                continue
            target = DEST / key / name
            target.parent.mkdir(parents=True, exist_ok=True)
            print("downloading", key, name, f"({round(size / 1e6, 1)} MB)" if size else "", flush=True)
            data = _get(url, timeout=3600)
            if md5 and hashlib.md5(data).hexdigest() != md5 and "dataverse" in src:  # noqa: S324 - Dataverse serves the ingested form of tables; its MD5 is the original's
                url += "?format=original"
                data = _get(url, timeout=3600)
            if md5 and hashlib.md5(data).hexdigest() != md5:  # noqa: S324 - integrity check against the publisher's checksum
                sys.exit(f"MD5 mismatch for {name}; nothing was kept.")
            target.write_bytes(data)
            record[f"{key}/{name}"] = {"url": url if "figshare" not in src else f"figshare article {FIGSHARE_DEPMAP}", "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(),
                                       "licence": src["licence"]}
            record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    if not args.list:
        print("done; recorded in", record_path)


if __name__ == "__main__":
    main()
