"""Read the Track 1 submission and turn it into the inputs of the repositioning analysis.

The Track 1 file (one row per proposed causal variant or compound-heterozygous pair, each with an `epcr` value
that ranks the rows) is the only source of "which genes matter here and how much". For every gene in it we keep the
rows, the best `epcr` (the gene weight) and, for single-nucleotide variants, the ClinVar classification looked up
by GRCh38 position. Nothing about the genes is typed by hand.
"""

from __future__ import annotations

import csv
import re
import urllib.parse
from pathlib import Path
from typing import Any

from track2.http_cache import get_json

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Nucleotide-level (cDNA) notation of the proband's variants is not published; protein notation and the transcript are kept.
CDNA = re.compile(r"\bc\.[-*]?\d+(?:[+-]\d+)?(?:_[-*]?\d+(?:[+-]\d+)?)?(?:[ACGT]+>[ACGT]+|del(?:ins)?[ACGT]*|dup[ACGT]*|ins[ACGT]*|[ACGT]*=)")
CDNA_MARKER = "[cDNA notation not published]"


def redact_cdna(text: str) -> str:
    return CDNA.sub(CDNA_MARKER, text)


def read_submission(path: Path) -> list[dict[str, Any]]:
    """Rows of the Track 1 CSV with parsed positions and the gene taken from the first word of the notes."""
    rows = []
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = {"chrom_1", "pos_1", "ref_1", "alt_1", "epcr", "finding_type", "notes"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: the Track 1 submission file lacks the columns {sorted(missing)}")
        for i, r in enumerate(reader, 1):
            if not (r["notes"] or "").split():
                raise ValueError(f"{path}, row {i}: `notes` is empty; its first word must be the gene symbol")
            try:
                variants = [(r["chrom_1"], int(r["pos_1"]), r["ref_1"], r["alt_1"])]
                if r.get("chrom_2"):
                    variants.append((r["chrom_2"], int(r["pos_2"]), r["ref_2"], r["alt_2"]))
                epcr = float(r["epcr"])
            except (ValueError, TypeError) as error:
                raise ValueError(f"{path}, row {i}: a position or the epcr value is not a number ({error})") from error
            rows.append({"row": i, "gene": r["notes"].split()[0], "finding_type": r["finding_type"], "epcr": epcr, "variants": variants, "notes": r["notes"]})
    if not rows:
        raise ValueError(f"{path}: the submission file has no rows")
    return rows


def clinvar_for_variant(chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any]:
    """ClinVar germline classification of the exact GRCh38 variant, or `not_in_clinvar` / `not_a_snv`."""
    if len(ref) != 1 or len(alt) != 1:
        return {"status": "not_a_snv", "classification": None}
    term = f"{chrom.replace('chr', '')}[Chromosome] AND {pos}[Base Position for Assembly GRCh38]"
    ids = get_json(f"{EUTILS}/esearch.fcgi?db=clinvar&retmode=json&term={urllib.parse.quote(term)}")["esearchresult"]["idlist"]
    if not ids:
        return {"status": "not_in_clinvar", "classification": None}
    summary = get_json(f"{EUTILS}/esummary.fcgi?db=clinvar&retmode=json&id={','.join(ids)}")["result"]
    for uid in ids:
        record = summary[uid]
        for variation in record.get("variation_set", []):
            spdi = variation.get("canonical_spdi") or ""
            if spdi.endswith(f":{pos - 1}:{ref}:{alt}"):
                cls = record["germline_classification"]
                return {"status": "in_clinvar", "classification": cls["description"], "review_status": cls["review_status"],
                        "clinvar_id": uid, "title": redact_cdna(record["title"])}
    return {"status": "other_allele_in_clinvar", "classification": None}


def summarise_genes(rows: list[dict[str, Any]], with_clinvar: bool = True) -> list[dict[str, Any]]:
    """One record per gene: rows, weight (best epcr), finding types and ClinVar results of its SNVs."""
    genes: dict[str, dict[str, Any]] = {}
    for r in rows:
        g = genes.setdefault(r["gene"], {"gene": r["gene"], "rows": [], "weight": 0.0, "finding_types": set(), "clinvar": []})
        g["rows"].append(r["row"])
        g["weight"] = max(g["weight"], r["epcr"])
        g["finding_types"].add(r["finding_type"])
        if with_clinvar:
            for chrom, pos, ref, alt in r["variants"]:
                if len(ref) == 1 and len(alt) == 1 and not r["notes"].split()[1:2] == ["single_hit"]:
                    g["clinvar"].append({"variant": f"{chrom}:{pos} {ref}>{alt}", **clinvar_for_variant(chrom, pos, ref, alt)})
    out = []
    for g in genes.values():
        g["clinvar"] = list({c["variant"]: c for c in g["clinvar"]}.values())
        g["finding_types"] = sorted(g["finding_types"])
        g["is_primary"] = "primary" in g["finding_types"]
        out.append(g)
    return sorted(out, key=lambda x: (-x["weight"], x["gene"]))
