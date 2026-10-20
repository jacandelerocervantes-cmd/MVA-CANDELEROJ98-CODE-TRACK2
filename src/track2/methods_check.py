"""Check the Track 2 methods description and write out/methods_check.json (sorted keys).

    python -m track2.methods_check

Counts the abstract words, lists the answered questions, scans the answers for privacy leaks (genomic coordinates, local paths, sample identifiers, cDNA notation)
and compares each answer of the Markdown with the matching cell of the Excel form.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from track2.fill_methods_form import ANSWER_CELLS, answers, read_cells

PRIVACY_PATTERNS = {
    "genomic_coordinate": r"\bchr[0-9XYM]{1,2}\s*[:_-]\s*\d{3,}|\bg\.\d{3,}|\bNC_0\d+",
    "cdna_notation": r"\bc\.\d+[ACGT>_+\-\d]*|\bNM_\d+",
    "local_path": r"[A-Za-z]:\\|/Users/|/home/|\\\\",
    "sample_id": r"\bSRR\d+|\bERR\d+|\bsample[_ ]?id\b|\bHG\d{5}\b",
}
REQUIRED_PHRASES = ["GitHub link: provided with the submission", "23 November 2026", "CC BY 4.0", "embargo", "hypotheses", "age is unknown"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods", type=Path, default=Path("submissions/candeleroj98_methods_description_track2.md"))
    parser.add_argument("--form", type=Path, default=Path("submissions/candeleroj98_methods_description_form_track2.xlsx"))
    parser.add_argument("--out", type=Path, default=Path("out/methods_check.json"))
    args = parser.parse_args()
    md = args.methods.read_text(encoding="utf-8")
    ans = answers(md)
    cells = read_cells(args.form)
    scan = {name: sorted({m.group(0) for a in ans.values() for m in re.finditer(pat, a)}) for name, pat in PRIVACY_PATTERNS.items()}
    equal = {str(q): cells[c] == ans[q] for q, c in ANSWER_CELLS.items()}
    text = " ".join(ans.values())
    report = {
        "abstract_word_count": len(re.findall(r"\S+", ans[11])),
        "answered_questions": sorted(q for q, a in ans.items() if a.strip()),
        "answer_length_characters": {str(q): len(a) for q, a in sorted(ans.items())},
        "md_equals_xlsx": equal,
        "md_equals_xlsx_all": all(equal.values()) and len(equal) == 11,
        "privacy_scan": {"matches": scan, "passed": not any(scan.values())},
        "required_phrases_present": {p: p in text for p in REQUIRED_PHRASES},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
