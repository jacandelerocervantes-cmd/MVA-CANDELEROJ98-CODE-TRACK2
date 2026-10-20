"""Pharmacogenomics, star-allele-LITE, from the CPIC allele definitions (public API https://api.cpicpgx.org/v1/).

What is done, per gene (CYP2C19, CYP2C9, CYP3A5, TPMT, NUDT15, DPYD, UGT1A1, SLCO1B1):
  1. CPIC `allele_definition` + `allele_location_value` + `sequence_location` give, per named allele, the variant
     nucleotide(s) at GRCh38 positions (`chromosomelocation` is `g.<pos><ref>><alt>` on the chromosome strand;
     the `variantallele` values are on that same strand, IUPAC codes mean "either base" and are not required).
  2. An allele REQUIRES the positions where its defining nucleotide set does not contain the GRCh38 reference base.
  3. The proband's VCF is scanned at those positions (`track2.pgx_vcf`). An allele is `possible` when every required
     position has a call carrying one of its nucleotides; `unresolved` when all its SNV positions are called but it also
     needs an indel/repeat/structural feature that is not interpreted here; `not_supported` when a required SNV has no
     matching call; `reference_like` when it has no required position (it cannot be evidenced by variant calls).

What is NOT concluded (stated in every output): the VCF has variant calls only (no reference calls) and is unphased.
A missing call is not evidence of the reference allele, `possible` alleles may lie on the same or on different
haplotypes, so NO diplotype and NO phenotype are asserted. CYP2D6 is excluded on purpose: it needs copy-number and
hybrid-gene (CYP2D7) calling that a short-variant VCF cannot provide; HLA genes need HLA typing. Output is a flag for
expert review, not a dosing or treatment statement.

Per-variant detail goes to `data/private/pgx_detail.json` (git-ignored); the publishable summary has no coordinates.

    python -m track2.pgx --vcf <vcf.gz> --candidates out/candidates.json --out out
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from track2.faers_disproportionality import clean_name
from track2.http_cache import get_json
from track2.pgx_vcf import norm_chrom, scan_vcf

CPIC = "https://api.cpicpgx.org/v1/"
GENES = ["CYP2C19", "CYP2C9", "CYP3A5", "TPMT", "NUDT15", "DPYD", "UGT1A1", "SLCO1B1"]
EXCLUDED = {
    "CYP2D6": "Structural gene (deletion, duplication, CYP2D7 hybrids): needs a dedicated caller; short-variant VCF is insufficient.",
    "HLA-A/HLA-B": "Highly polymorphic; require HLA typing, not SNV calls.",
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
SNV = re.compile(r"g\.(\d+)([ACGT])>([ACGT])")
LIMITATION = ("The VCF contains variant calls only and is unphased: a missing call is not evidence of the reference allele, "
              "so no diplotype or phenotype is asserted. 'Possible' alleles are consistent with the calls, not confirmed. "
              "Indel/repeat-defined features are not interpreted. Flag for expert review, not a clinical result.")
SPLIT = 150  # ids per CPIC `in.(...)` request


def parse_location(text: str | None) -> tuple[int, str, set[str]] | None:
    """`g.123A>G` or `g.123G>A; g.123G>T` -> (pos, ref, alts); None for indels, duplications and repeats."""
    parts = [SNV.fullmatch(p.strip()) for p in (text or "").split(";")]
    if not parts or any(m is None for m in parts) or len({m.group(1) for m in parts if m}) != 1:
        return None
    first = parts[0]
    assert first is not None
    return int(first.group(1)), first.group(2), {m.group(3) for m in parts if m}


def load_cpic(gene: str, get: Any = get_json) -> dict[str, Any]:
    info = get(f"{CPIC}gene?symbol=eq.{gene}&select=symbol,chr,chromosequenceid")[0]
    defs = get(f"{CPIC}allele_definition?genesymbol=eq.{gene}&limit=1000")
    locs = get(f"{CPIC}sequence_location?genesymbol=eq.{gene}&limit=5000")
    values: list[dict[str, Any]] = []
    for i in range(0, len(defs), SPLIT):
        ids = ",".join(str(d["id"]) for d in defs[i:i + SPLIT])
        values += get(f"{CPIC}allele_location_value?alleledefinitionid=in.({ids})&limit=100000")
    return {"gene": info, "alleles": defs, "locations": {loc["id"]: loc for loc in locs}, "values": values}


def build_alleles(cpic: dict[str, Any]) -> tuple[dict[int, tuple[int, str, set[str]]], dict[str, dict[str, Any]]]:
    """Assessable SNV locations by id, and per allele its required SNV positions and unassessable features."""
    snv = {i: p for i, loc in cpic["locations"].items() if (p := parse_location(loc.get("chromosomelocation")))}
    by_allele: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for v in cpic["values"]:
        by_allele[v["alleledefinitionid"]].append(v)
    alleles: dict[str, dict[str, Any]] = {}
    for d in cpic["alleles"]:
        required: dict[int, set[str]] = {}
        unassessable = 0
        for v in by_allele.get(d["id"], []):
            loc, base = snv.get(v["locationid"]), v["variantallele"]
            if loc is None or not set(base) <= set("ACGTRYSWKMBDHVN") or len(base) != 1:
                unassessable += 1
            elif loc[1] not in IUPAC[base]:  # reference base not allowed -> the allele needs a variant here
                required[loc[0]] = set(IUPAC[base])
        alleles[d["name"]] = {"required": required, "unassessable": unassessable, "structural": bool(d["structuralvariation"])}
    return {i: p for i, p in snv.items()}, alleles


def classify(allele: dict[str, Any], detected: dict[int, set[str]]) -> str:
    if allele["structural"]:
        return "unresolved"
    if any(not (detected.get(pos, set()) & bases) for pos, bases in allele["required"].items()):
        return "not_supported"
    if allele["unassessable"]:
        return "unresolved"
    return "possible" if allele["required"] else "reference_like"


def most_specific(possible: list[str], alleles: dict[str, dict[str, Any]]) -> list[str]:
    """Possible alleles whose required positions are not a strict subset of another possible allele's."""
    req = {n: set(alleles[n]["required"]) for n in possible}
    return sorted(n for n in possible if not any(req[n] < req[m] for m in possible))


def evaluate_gene(gene: str, cpic: dict[str, Any], calls: dict[str, dict[int, list[dict[str, Any]]]]) -> tuple[dict[str, Any], dict[str, Any]]:
    """(public summary, private detail) of one gene."""
    snv, alleles = build_alleles(cpic)
    chrom = norm_chrom(cpic["gene"]["chr"])
    by_pos = {p[0]: (loc_id, p) for loc_id, p in snv.items()}
    observed = calls.get(chrom, {})
    detected: dict[int, set[str]] = {}
    detail_rows, mismatches, nonref = [], 0, 0
    for pos, (loc_id, (_, ref, alts)) in sorted(by_pos.items()):
        obs = observed.get(pos, [])
        good = [o for o in obs if o["ref"] == ref]
        mismatches += len(obs) - len(good)
        called = {b for o in good for b in o["alt_bases_called"]}
        if called:
            nonref += 1
            detected[pos] = called
        if obs:
            detail_rows.append({"chrom": chrom, "pos": pos, "dbsnp": cpic["locations"][loc_id].get("dbsnpid"), "cpic_ref": ref,
                                "cpic_alts": sorted(alts), "observations": obs})
    status = {name: classify(a, detected) for name, a in alleles.items()}
    names = {s: sorted(n for n, v in status.items() if v == s) for s in ("possible", "unresolved", "reference_like", "not_supported")}
    summary = {
        "n_defining_positions_checked": len(by_pos),
        "n_positions_not_interpretable_indel_or_repeat": len(cpic["locations"]) - len(snv),
        "n_positions_with_nonreference_call": nonref,
        "n_vcf_records_with_ref_mismatch_vs_cpic_grch38": mismatches,
        "n_alleles_defined_in_cpic": len(alleles),
        "possible_alleles": names["possible"], "most_specific_possible_alleles": most_specific(names["possible"], alleles),
        "unresolved_alleles_need_indel_repeat_or_structural": names["unresolved"],
        "reference_like_alleles_not_evidenced_by_variant_calls": names["reference_like"],
        "n_alleles_not_supported_by_calls": len(names["not_supported"]),
        "limitation": LIMITATION,
    }
    return summary, {"calls_at_defining_positions": detail_rows,
                     "required_positions_per_possible_allele": {n: sorted(alleles[n]["required"]) for n in names["possible"]}}


def drug_relevance(candidates: list[dict[str, Any]], get: Any = get_json) -> dict[str, list[dict[str, Any]]]:
    """drug_id -> CPIC drug-gene pairs (match by cleaned lower-case name); [] when CPIC has no pair for the drug."""
    pairs = get(f"{CPIC}pair_view?limit=5000")
    by_name: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for p in pairs:
        by_name[p["drugname"].lower()].append({
            "gene": p["genesymbol"], "cpic_level": p.get("cpiclevel"), "clinpgx_level": p.get("clinpgxlevel"),
            "guideline_name": p.get("guidelinename"), "guideline_url": p.get("guidelineurl"), "provisional": p.get("provisional")})
    return {c["drug_id"]: sorted(by_name.get(clean_name(c["name"] or "").lower(), []), key=lambda r: (r["gene"], r["guideline_name"] or ""))
            for c in candidates}


def run(vcf: Path, candidates: list[dict[str, Any]], out: Path, private: Path, get: Any = get_json) -> dict[str, Any]:
    cpic = {g: load_cpic(g, get) for g in GENES}
    wanted: dict[str, list[int]] = defaultdict(list)
    for g, data in cpic.items():
        snv, _ = build_alleles(data)
        wanted[data["gene"]["chr"]] += [p[0] for p in snv.values()]
    calls = scan_vcf(vcf, wanted)
    summary, detail = {}, {}
    for g in GENES:
        summary[g], detail[g] = evaluate_gene(g, cpic[g], calls)
    relevance = drug_relevance(candidates, get)
    result = {"source": "CPIC API v1 (allele_definition, allele_location_value, sequence_location, pair_view); GRCh38",
              "genes": summary, "excluded_genes": EXCLUDED, "limitation": LIMITATION,
              "drugs_with_cpic_pair": sum(1 for v in relevance.values() if v), "drugs_total": len(relevance)}
    private.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    private.write_text(json.dumps(detail, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    private.parent.mkdir(parents=True, exist_ok=True)  # the summary lists possible star alleles of the patient: private, never published
    (private.parent / "pgx_summary.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    (out / "pgx_drug_relevance.json").write_text(json.dumps(relevance, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, default=Path("out/candidates.json"))
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--private", type=Path, default=Path("data/private/pgx_detail.json"))
    args = parser.parse_args()
    result = run(args.vcf, json.loads(args.candidates.read_text(encoding="utf-8")), args.out, args.private)
    print(json.dumps({g: {k: v for k, v in s.items() if k != "limitation"} for g, s in result["genes"].items()}, indent=1))


if __name__ == "__main__":
    main()
