"""The gene-level Track 1 extract must reproduce the outputs of the two modules that read Track 1, and must carry no patient-level string."""

import json
import re
from pathlib import Path

from track2 import differential as diff
from track2 import track1_extract as ex
from track2 import variant_level_consequence as vlc

ROOT = Path(__file__).resolve().parents[1]


def test_rerun_from_the_extract_reproduces_the_differential_diagnosis() -> None:
    landscape = json.loads((ROOT / "out" / "mva_landscape.json").read_text(encoding="utf8"))
    submission = json.loads((ROOT / "out" / "submission_genes.json").read_text(encoding="utf8"))
    result = diff.run(landscape, submission, *ex.differential_inputs())
    expected = json.loads((ROOT / "out" / "differential_diagnosis.json").read_text(encoding="utf8"))["genes"]
    assert json.loads(json.dumps(result, sort_keys=True, default=list)) == expected


def test_rerun_from_the_extract_reproduces_the_variant_level_classes() -> None:
    doc = vlc.build(*ex.variant_level_inputs())
    expected = json.loads((ROOT / "out" / "variant_consequence_classes.json").read_text(encoding="utf8"))
    assert json.loads(json.dumps(doc, sort_keys=True)) == expected


def test_the_extract_has_no_position_allele_genotype_or_read_metric() -> None:
    text = ex.EXTRACT.read_text(encoding="utf8")
    for pattern in (r"chr?[0-9XYM]{1,2}[:_]\d{4,}", r"(?<![\d.])\d{6,}(?![\d.])", r"\b[01][/|][01]\b(?!\d)", r"\bc\.[-*]?\d+", r"[ACGT]>[ACGT]", r"\b(GQ|DP|QUAL|VAF)\b"):
        assert not re.search(pattern, text), pattern
