"""Deterministic differential diagnosis of the Track 1 genes (rules: protocol/addendum_differential.md).

Run:  PYTHONPATH=src python -m track2.differential

Gene-level fields only. Genomic positions, alleles, genotypes and read metrics are never read or stored: the Track 1
ranking is read through an explicit column whitelist and the structural-call file through `finding_type` and the first
word of `notes` (a gene symbol).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TRACK1 = ROOT.parent / "HSelith_MVA-Track1" / "submissions" / "track1"
OUT = ROOT / "out" / "differential_diagnosis.json"

RANKING_COLUMNS = ("gene", "type", "n_high_impact", "log10_joint_rarity", "phenotype_score")
EXCLUDED_FINDING_TYPES = {"artefact", "common_variant", "filtered"}
MIN_LOG10_JOINT_RARITY = 2.0  # joint frequency above 1e-2 is a common-variant pair
BIALLELIC_TYPES = {"compound_het", "homozygous"}


# ---------------------------------------------------------------- loaders (gene-level only)
def load_ranking(path: Path) -> dict:
    out = {}
    with open(path, encoding="utf8", newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["gene"]] = {
                "type": row["type"],
                "n_high_impact": int(row["n_high_impact"]),
                "log10_joint_rarity": float(row["log10_joint_rarity"]),
                "phenotype_score": float(row["phenotype_score"]),
            }
    return out


def load_structural(path: Path) -> dict:
    """gene -> set of finding types; only `finding_type` and the first word of `notes` are read."""
    out: dict = {}
    with open(path, encoding="utf8", newline="") as fh:
        for row in csv.DictReader(fh):
            words = (row.get("notes") or "").split()
            if words:
                out.setdefault(words[0], set()).add(row["finding_type"])
    return out


def load_phenotype(path: Path) -> dict:
    d = json.loads(Path(path).read_text(encoding="utf8"))
    return {
        "scores": d.get("candidate_scores", {}),
        "ranks": d.get("candidate_phenotype_rank_among_all_genes", {}),
    }


# ---------------------------------------------------------------- criteria
def c1_status(landscape_gene: dict | None) -> str:
    if not landscape_gene:
        return "none"
    return "curated" if landscape_gene.get("tier") == "curated" else "proposed"


def c2_pattern(ranking_row: dict | None, inheritance: list) -> str:
    if ranking_row is None:
        return "unpaired"
    text = " ".join(inheritance or []).lower()
    typ = ranking_row["type"]
    if typ in BIALLELIC_TYPES:
        return "consistent" if "recessive" in text else "biallelic_inheritance_not_recessive"
    if typ == "hemizygous":
        return "consistent" if "x-linked" in text else "inconsistent"
    return "inconsistent"


def c3_consequence(ranking_row: dict | None) -> str:
    if ranking_row is None:
        return "structural_unpaired"
    return "high_impact" if ranking_row["n_high_impact"] >= 1 else "non_high_impact"


def c5_clinvar(records: list) -> str:
    classes = [(r.get("classification") or "").lower() for r in records or []]
    if any("pathogenic" in c and "uncertain" not in c and "conflicting" not in c for c in classes):
        return "P/LP"
    if any("uncertain" in c or "conflicting" in c for c in classes):
        return "uncertain"
    return "none"


def c4_phenotype(gene: str, score_table: dict, rank_table: dict, landscape_gene: dict | None) -> dict:
    rank = rank_table.get(gene)
    if gene in rank_table:
        has_hpo = rank is not None
    else:
        has_hpo = bool(landscape_gene) and any("HPO" in s for s in landscape_gene.get("sources", []))
    return {"phenotype_score": score_table.get(gene), "phenotype_rank": rank, "has_hpo": has_hpo}


# ---------------------------------------------------------------- verdict
def verdict(c: dict) -> tuple:
    """Rules of the addendum, first match wins. `c` carries c1..c5 values plus finding_types and joint rarity."""
    if (set(c.get("finding_types", [])) & EXCLUDED_FINDING_TYPES) or (
        c.get("log10_joint_rarity") is not None and c["log10_joint_rarity"] < MIN_LOG10_JOINT_RARITY
    ):
        return "excluded", "common variant or artefact family (rule 1)"
    if c["c1"] == "curated" and c["c2"] == "consistent" and (c["c3"] == "high_impact" or c["c5"] == "P/LP"):
        return "causal candidate", (
            "curated MVA gene, two rare alleles with recessive annotation, and a high-impact allele or ClinVar P/LP "
            "(phase not resolved at gene level)"
        )
    if c["c1"] == "none":
        return "insufficient evidence", "not an MVA gene in the landscape"
    if not c["c4"]["has_hpo"]:
        return "insufficient evidence", "no HPO annotation"
    if c["c2"] == "consistent" and c["c3"] != "high_impact" and c["c5"] != "P/LP":
        return "insufficient evidence", "biallelic genotype but no high-impact allele and no ClinVar P/LP (uncertain only)"
    return "contextual", "MVA gene without a qualifying coding genotype (single or unpaired call, or not recessive-consistent)"


def evaluate_gene(gene, landscape_gene, ranking_row, structural_types, submission_gene, pheno) -> dict:
    inheritance = (landscape_gene or {}).get("inheritance_as_annotated", [])
    finding_types = sorted(set(structural_types or []) | set((submission_gene or {}).get("finding_types", [])))
    c = {
        "c1": c1_status(landscape_gene),
        "c2": c2_pattern(ranking_row, inheritance),
        "c3": c3_consequence(ranking_row),
        "c4": c4_phenotype(gene, pheno["scores"], pheno["ranks"], landscape_gene),
        "c5": c5_clinvar((submission_gene or {}).get("clinvar", [])),
        "finding_types": finding_types,
        "log10_joint_rarity": ranking_row["log10_joint_rarity"] if ranking_row else None,
        "genotype_type": ranking_row["type"] if ranking_row else None,
    }
    v, reason = verdict(c)
    c["verdict"], c["reason"] = v, reason
    return c


def run(landscape: dict, submission: dict, ranking: dict, structural: dict, pheno: dict) -> dict:
    land = {g["symbol"]: g for g in landscape["genes"]}
    sub = {g["gene"]: g for g in submission["genes"]}
    genes = sorted(set(ranking) | set(structural))
    return {g: evaluate_gene(g, land.get(g), ranking.get(g), structural.get(g), sub.get(g), pheno) for g in genes}


def main() -> None:
    landscape = json.loads((ROOT / "out" / "mva_landscape.json").read_text(encoding="utf8"))
    submission = json.loads((ROOT / "out" / "submission_genes.json").read_text(encoding="utf8"))
    if (TRACK1 / "candeleroj98_genomewide_ranking.csv").exists():
        inputs = (
            load_ranking(TRACK1 / "candeleroj98_genomewide_ranking.csv"),
            load_structural(TRACK1 / "candeleroj98_gatk-manta-fusion.csv"),
            load_phenotype(TRACK1 / "candeleroj98_phenotype_summary.json"),
        )
    else:  # the Track 1 repository is not available: use the gene-level extract kept in config/
        from track2 import track1_extract

        inputs = track1_extract.differential_inputs()
    result = run(landscape, submission, *inputs)
    doc = {"rules": "protocol/addendum_differential.md", "genes": result}
    OUT.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf8")
    for g, r in result.items():
        print(g, r["verdict"])


if __name__ == "__main__":
    main()
