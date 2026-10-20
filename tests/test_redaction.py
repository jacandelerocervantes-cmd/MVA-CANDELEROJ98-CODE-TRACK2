"""Nucleotide-level notation of the proband's variants must not reach the written files."""

from track2.submission_input import CDNA_MARKER, redact_cdna


def test_cdna_notation_is_replaced_and_protein_notation_kept() -> None:
    title = "NM_000000.1(GENEX):c.123A>G (p.Lys41Arg)"
    assert redact_cdna(title) == f"NM_000000.1(GENEX):{CDNA_MARKER} (p.Lys41Arg)"
    for notation in ("c.100del", "c.100_102delinsTT", "c.-12A>G", "c.*5dup", "c.88+2T>C", "c.9ins", "c.5="):
        assert CDNA_MARKER in redact_cdna(f"X:{notation} (p.Ter)") and notation not in redact_cdna(f"X:{notation} (p.Ter)")


def test_text_without_cdna_is_unchanged() -> None:
    assert redact_cdna("a pathogenic variant in the c. terminal region") == "a pathogenic variant in the c. terminal region"
