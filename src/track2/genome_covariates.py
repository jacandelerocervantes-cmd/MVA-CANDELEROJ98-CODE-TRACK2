"""Covariates derived from the proband's genome (aggregates only; the VCF is not published).

1. Sex: heterozygosity pattern on chrX outside the pseudoautosomal regions.
2. Mosaic aneuploidy signal: relative depth per chromosome and deviation of the allele fraction
   of heterozygotes. It is a SIGNAL, not a karyotype: a low-fraction mosaic or one in another tissue is not seen.
3. Rare coding-effect variants in pharmacogenes and tumour predisposition genes. This is not star-allele
   calling (CYP2D6 and others require dedicated tools); it is a list of variants with their frequency.

    python -m track2.genome_covariates --vcf <vcf.gz> --gtf <MANE.gtf.gz> --cds-cache <cds_sequences.json> \
        --gnomad-cache <gnomad_v4_frequencies.json> --out out
"""

from __future__ import annotations

import argparse
import gzip
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

from track2.vendor.genome_wide_search import (
    compute_effects, infer_hemizygous_x, load_mane_transcripts, scan_vcf_coding, variant_key,
)

PGX_GENES = ["CYP3A4", "CYP3A5", "UGT1A1", "DPYD", "TPMT", "NUDT15", "ABCB1", "CYP2D6", "CYP2C9", "CYP2C19", "SLCO1B1"]
PREDISPOSITION_GENES = ["TP53", "RB1", "WT1", "DICER1", "NF1", "PTCH1", "CDKN2A", "ATM", "BRCA2", "SMARCB1"]
RARE_AF = 1e-3
DEPTH_FLAG_LOW, DEPTH_FLAG_HIGH = 0.95, 1.05
MIN_HET_SITES = 200


def chromosome_signal(vcf: Path) -> dict[str, Any]:
    """Median depth and allele-fraction deviation in PASS heterozygotes, per chromosome."""
    depth: dict[str, list[int]] = defaultdict(list)
    het_dev: dict[str, list[float]] = defaultdict(list)
    with gzip.open(vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line[0] == "#":
                continue
            parts = line.split("\t", 10)
            if parts[6] != "PASS" or len(parts) < 10 or len(parts[3]) != 1 or len(parts[4]) != 1:
                continue
            keys = parts[8].split(":")
            values = parts[9].strip().split(":")
            if "DP" not in keys or "AD" not in keys or len(values) < len(keys):
                continue
            try:
                dp = int(values[keys.index("DP")])
                ref_alt = [int(x) for x in values[keys.index("AD")].split(",")]
            except ValueError:
                continue
            chrom = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
            depth[chrom].append(dp)
            if values[0].replace("|", "/") in ("0/1", "1/0") and sum(ref_alt) >= 10:
                het_dev[chrom].append(abs(ref_alt[1] / sum(ref_alt) - 0.5))
    autosomes = [f"chr{i}" for i in range(1, 23)]
    genome = statistics.median([statistics.median(depth[c]) for c in autosomes if depth[c]])
    rows = {}
    for chrom in autosomes + ["chrX", "chrY"]:
        if not depth[chrom]:
            continue
        ratio = statistics.median(depth[chrom]) / genome
        dev = statistics.median(het_dev[chrom]) if len(het_dev[chrom]) >= MIN_HET_SITES else None
        rows[chrom] = {
            "median_depth": statistics.median(depth[chrom]), "depth_ratio": round(ratio, 3), "n_het_sites": len(het_dev[chrom]),
            "median_het_af_deviation": None if dev is None else round(dev, 4),
            "flag": chrom in autosomes and not (DEPTH_FLAG_LOW <= ratio <= DEPTH_FLAG_HIGH),
        }
    return {"genome_median_depth": genome, "per_chromosome": rows,
            "thresholds": {"depth_ratio_normal_range": [DEPTH_FLAG_LOW, DEPTH_FLAG_HIGH]},
            "limits": "Signal in blood; a low-fraction mosaic or one in another tissue is not detected."}


def variant_panel(effects: list[dict[str, Any]], frequencies: dict[str, Any], genes: list[str]) -> list[dict[str, Any]]:
    out = []
    for e in effects:
        if e["gene"] not in genes or e["impact"] not in ("HIGH", "MODERATE"):
            continue
        freq = frequencies.get(variant_key(e), {})
        af = freq.get("af") if freq.get("present") else None
        if freq.get("present") and af is not None and af >= RARE_AF:
            continue
        out.append({"gene": e["gene"], "chrom": e["chrom"], "pos": e["pos"], "consequence": e["consequence"], "hgvsp": e.get("hgvsp"),
                    "genotype": e["gt"], "gnomad_af": af, "gnomad_checked": bool(freq)})
    return sorted(out, key=lambda r: (r["gene"], r["pos"]))


def _count_by_gene(variants: list[dict[str, Any]], genes: list[str]) -> dict[str, int]:
    counts = {g: 0 for g in genes}
    for v in variants:
        counts[v["gene"]] += 1
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--cds-cache", type=Path, required=True)
    parser.add_argument("--gnomad-cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("data/private"), help="per-variant detail (not versioned)")
    parser.add_argument("--summary-out", type=Path, default=Path("data/private"), help="publishable summary")
    args = parser.parse_args()
    transcripts = load_mane_transcripts(args.gtf)
    hits = scan_vcf_coding(args.vcf, transcripts)
    effects = compute_effects(hits, transcripts, json.loads(args.cds_cache.read_text(encoding="utf-8")))
    frequencies = json.loads(args.gnomad_cache.read_text(encoding="utf-8"))
    result = {
        "sex_inference": {"single_x": infer_hemizygous_x(hits), "method": "fraction of heterozygous genotypes on chrX outside the PAR (<= 25 %)"},
        "aneuploidy_signal": chromosome_signal(args.vcf),
        "pharmacogene_variants": {"note": "rare HIGH/MODERATE variants; not star-allele calling", "variants": variant_panel(effects, frequencies, PGX_GENES)},
        "predisposition_gene_variants": {"note": "rare HIGH/MODERATE variants; no pathogenicity assignment", "variants": variant_panel(effects, frequencies, PREDISPOSITION_GENES)},
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "genome_covariates.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    summary = {
        "sex_inference": result["sex_inference"], "aneuploidy_signal": result["aneuploidy_signal"],
        "pharmacogene_variants_by_gene": _count_by_gene(result["pharmacogene_variants"]["variants"], PGX_GENES),
        "predisposition_variants_by_gene": _count_by_gene(result["predisposition_gene_variants"]["variants"], PREDISPOSITION_GENES),
        "note": "Publishable summary: no coordinates or individual genotypes. The per-variant detail stays in data/private (not versioned).",
    }
    args.summary_out.mkdir(parents=True, exist_ok=True)
    (args.summary_out / "genome_covariates_summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: (v if k == "sex_inference" else "ok") for k, v in result.items()}))


if __name__ == "__main__":
    main()
