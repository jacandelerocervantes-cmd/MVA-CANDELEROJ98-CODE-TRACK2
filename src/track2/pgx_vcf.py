"""Scan of a germline SNV/indel VCF for the calls at a fixed set of GRCh38 positions (pharmacogene defining variants).

The VCF holds variant calls only: a position without a record is NOT evidence of the reference allele (it can also be
a no-coverage or filtered position). Only records with the same length of REF and ALT are read per base (SNV and
MNP); insertions and deletions are not interpreted here. Chromosome names are accepted with or without `chr`.
"""

from __future__ import annotations

import bisect
import gzip
from collections import defaultdict
from pathlib import Path
from typing import Any


def norm_chrom(name: str) -> str:
    return name[3:] if name.startswith("chr") else name


def _called_alleles(genotype: str, n_alts: int) -> set[int]:
    """Indices (1-based into ALT) of the alternate alleles present in a GT string; unphased/phased alike."""
    found = set()
    for token in genotype.replace("|", "/").split("/"):
        if token.isdigit() and 1 <= int(token) <= n_alts:
            found.add(int(token))
    return found


def scan_vcf(vcf: Path, wanted: dict[str, list[int]]) -> dict[str, dict[int, list[dict[str, Any]]]]:
    """chrom -> position -> observations `{ref, alt_bases_called, gt, filter}` for the wanted positions (any chrom naming)."""
    positions = {norm_chrom(c): sorted(set(p)) for c, p in wanted.items()}
    out: dict[str, dict[int, list[dict[str, Any]]]] = {c: defaultdict(list) for c in positions}
    with gzip.open(vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line[0] == "#":
                continue
            cut = line.index("\t")
            chrom = norm_chrom(line[:cut])
            sites = positions.get(chrom)
            if sites is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            pos, ref, alts = int(parts[1]), parts[3], parts[4].split(",")
            lo = bisect.bisect_left(sites, pos)
            hi = bisect.bisect_right(sites, pos + len(ref) - 1)
            if lo == hi:
                continue
            keys = parts[8].split(":")
            values = parts[9].split(":")
            genotype = values[keys.index("GT")] if "GT" in keys and len(values) > keys.index("GT") else "./."
            called = _called_alleles(genotype, len(alts))
            for site in sites[lo:hi]:
                offset = site - pos
                bases = {alts[i - 1][offset] for i in called if len(alts[i - 1]) == len(ref)}
                out[chrom][site].append({"ref": ref[offset], "alt_bases_called": sorted(bases - {ref[offset]}),
                                         "gt": genotype, "filter": parts[6]})
    return {c: dict(v) for c, v in out.items()}
