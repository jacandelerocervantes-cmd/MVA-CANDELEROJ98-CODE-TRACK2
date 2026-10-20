"""Variant-level consequence classes per Track 1 gene (answers limit A of the rule-based differential diagnosis).

Run:  PYTHONPATH=src python -m track2.variant_level_consequence

Limit A: `differential.c3_consequence` uses the gene-level count `n_high_impact`, a proxy, because predicted
nonsense-mediated decay (NMD) is not a gene-level field. Here every variant is reduced, privately and in memory, to a
categorical consequence class, and only classes and counts are written to out/variant_consequence_classes.json.

Privacy: the only substrings ever extracted from the Track 1 files are (i) a whitelisted Sequence Ontology consequence
term (e.g. stop_gained), (ii) a structural-call type word (DEL, INS, DUP, INV, BND), (iii) the yes/no answer of the
text "Overlaps MANE Select coding sequence", (iv) the presence (not the content) of a decay annotation, and (v) the
first word of `notes` (a gene symbol). Position, allele, genotype and read-metric columns/text are never parsed,
stored, or printed. Nothing is written except gene -> class -> count.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TRACK1 = ROOT.parent / "HSelith_MVA-Track1" / "submissions" / "track1"
OUT = ROOT / "out" / "variant_consequence_classes.json"

LOF = "loss of function"
MISSENSE = "missense"
INFRAME = "in-frame"
OTHER = "other coding or non-coding consequence"
SV_NONCODING = "structural call not overlapping coding sequence"
SV_CODING = "structural call overlapping coding sequence"
SV_UNKNOWN = "structural call, coding overlap not stated"

# Sequence Ontology term -> class. Only these exact terms are ever extracted from the text.
SO_TO_CLASS = {
    "stop_gained": LOF,
    "frameshift_variant": LOF,
    "splice_donor_variant": LOF,
    "splice_acceptor_variant": LOF,
    "start_lost": LOF,
    "missense_variant": MISSENSE,
    "inframe_insertion": INFRAME,
    "inframe_deletion": INFRAME,
    "stop_lost": OTHER,
    "synonymous_variant": OTHER,
    "splice_region_variant": OTHER,
    "protein_altering_variant": OTHER,
}
_SO_RE = re.compile(r"\b(" + "|".join(sorted(SO_TO_CLASS, key=len, reverse=True)) + r")\b")
_SV_TYPE_RE = re.compile(r"\b(DEL|INS|DUP|INV|BND)\b")
_MANE_RE = re.compile(r"Overlaps MANE Select coding sequence:\s*(yes|no)\b", re.IGNORECASE)
_DECAY_RE = re.compile(r"(nonsense[- ]mediated|\bNMD\b|predicted[- ]decay|escapes? (?:NMD|decay))", re.IGNORECASE)
_DECAY_NEG_RE = re.compile(r"(NMD[- ]escape|escapes? (?:NMD|decay)|not (?:predicted|subject) to (?:NMD|decay))", re.IGNORECASE)
_PAIR_ROW_RE = re.compile(r"\b(two|pair|reciprocal|alternate representation|single-variant component)\b", re.IGNORECASE)


def consequence_class(segment: str) -> str | None:
    """Class of one variant description; only the whitelisted Sequence Ontology term is read."""
    m = _SO_RE.search(segment or "")
    return SO_TO_CLASS[m.group(1)] if m else None


def decay_annotation(text: str) -> str:
    """'predicted' / 'escape' / 'not_annotated' from free text (presence of the concept only)."""
    if _DECAY_NEG_RE.search(text or ""):
        return "escape"
    return "predicted" if _DECAY_RE.search(text or "") else "not_annotated"


def structural_class(notes: str) -> str:
    m = _MANE_RE.search(notes or "")
    if not m:
        return SV_UNKNOWN
    return SV_CODING if m.group(1).lower() == "yes" else SV_NONCODING


# ---------------------------------------------------------------- readers (in memory only)
def read_ranking_classes(path: Path) -> dict:
    """gene -> list of class strings, one per variant, from the `variants` column (terms only)."""
    out: dict = {}
    with open(path, encoding="utf8", newline="") as fh:
        for row in csv.DictReader(fh):
            for seg in (row.get("variants") or "").split("|"):
                cls = consequence_class(seg)
                if cls is not None:
                    out.setdefault(row["gene"], []).append(cls)
    return out


def read_structural(path: Path) -> tuple[dict, dict]:
    """(gene -> list of structural classes for single SV calls, gene -> decay annotation state from notes)."""
    sv: dict = {}
    decay: dict = {}
    with open(path, encoding="utf8", newline="") as fh:
        for row in csv.DictReader(fh):
            notes = row.get("notes") or ""
            words = notes.split()
            if not words:
                continue
            gene = words[0]
            state = decay_annotation(notes)
            if decay.get(gene, "not_annotated") == "not_annotated":
                decay[gene] = state
            if _SV_TYPE_RE.search(notes) and not _PAIR_ROW_RE.search(notes):
                sv.setdefault(gene, []).append(structural_class(notes))
    return sv, decay


def read_phenotype_decay(path: Path) -> bool:
    """True if any string leaf of the phenotype summary mentions a decay annotation (concept only)."""
    d = json.loads(Path(path).read_text(encoding="utf8"))

    def walk(x) -> bool:
        if isinstance(x, str):
            return bool(_DECAY_RE.search(x))
        if isinstance(x, dict):
            return any(walk(v) for v in list(x.keys()) + list(x.values()))
        if isinstance(x, list):
            return any(walk(v) for v in x)
        return False

    return walk(d)


# ---------------------------------------------------------------- assembly
def build(snv_classes: dict, sv_classes: dict, decay_states: dict, phenotype_has_decay: bool = False) -> dict:
    genes = {}
    for gene in sorted(set(snv_classes) | set(sv_classes)):
        counts: dict = {}
        for c in snv_classes.get(gene, []) + sv_classes.get(gene, []):
            counts[c] = counts.get(c, 0) + 1
        state = decay_states.get(gene, "not_annotated")
        if state == "not_annotated" and phenotype_has_decay:
            state = "predicted"
        n_lof = counts.get(LOF, 0)
        genes[gene] = {
            "class_counts": dict(sorted(counts.items())),
            "n_variants_or_calls": sum(counts.values()),
            "n_loss_of_function": n_lof,
            "decay_annotation": state if n_lof else "not_applicable",
        }
    return {
        "note": (
            "Categorical classes and counts only. Predicted nonsense-mediated decay is not annotated in the Track 1 "
            "files, so loss-of-function variants are reported with decay_annotation = not_annotated; decay is not inferred."
        ),
        "genes": genes,
    }


def c3_from_classes(entry: dict | None) -> dict:
    """Variant-level replacement for the gene-level `n_high_impact` proxy of differential.c3_consequence.

    Returns {"c3": ..., "decay_status": ...}. `c3` keeps the vocabulary used by differential.verdict:
    'high_impact' (at least one loss-of-function variant), 'non_high_impact' (coding SNV/indel classes only),
    'structural_noncoding' (only structural calls that do not overlap the coding sequence; never high impact),
    'structural_unpaired' (structural call overlapping coding sequence or coding overlap not stated), 'none'.
    `decay_status` is the Track 1 decay annotation of the loss-of-function variants (never inferred).
    """
    if not entry or not entry.get("n_variants_or_calls"):
        return {"c3": "none", "decay_status": "not_applicable"}
    counts = entry["class_counts"]
    if entry.get("n_loss_of_function", 0) >= 1:
        return {"c3": "high_impact", "decay_status": entry.get("decay_annotation", "not_annotated")}
    snv_like = sum(counts.get(k, 0) for k in (MISSENSE, INFRAME, OTHER))
    if snv_like:
        return {"c3": "non_high_impact", "decay_status": "not_applicable"}
    if counts.get(SV_CODING, 0) or counts.get(SV_UNKNOWN, 0):
        return {"c3": "structural_unpaired", "decay_status": "not_applicable"}
    return {"c3": "structural_noncoding", "decay_status": "not_applicable"}


def main() -> None:
    if (TRACK1 / "candeleroj98_genomewide_ranking.csv").exists():
        snv = read_ranking_classes(TRACK1 / "candeleroj98_genomewide_ranking.csv")
        sv, decay = read_structural(TRACK1 / "candeleroj98_gatk-manta-fusion.csv")
        pheno_decay = read_phenotype_decay(TRACK1 / "candeleroj98_phenotype_summary.json")
    else:  # the Track 1 repository is not available: use the gene-level extract kept in config/
        from track2 import track1_extract

        snv, sv, decay, pheno_decay = track1_extract.variant_level_inputs()
    doc = build(snv, sv, decay, pheno_decay)
    OUT.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf8")
    for g, e in doc["genes"].items():
        print(g, e["class_counts"], "decay:", e["decay_annotation"])


if __name__ == "__main__":
    main()
