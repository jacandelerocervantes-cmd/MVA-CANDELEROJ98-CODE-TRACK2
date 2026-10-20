"""Write the answers of the methods description (Markdown, official template order) into the official Excel form, sheet 'Track 2 methods', cells B7 to B17.

    python -m track2.fill_methods_form --methods submissions/candeleroj98_methods_description_track2.md --form submissions/candeleroj98_methods_description_form_track2.xlsx

The answer cells hold inline strings; only those eleven cells (empty in the template, or inline strings) are written, by direct text substitution in the sheet XML (the XML is never re-serialised, so
namespace declarations, styles, comments and the Track 1 sheet stay byte for byte as in the form). Standard library only.
"""

from __future__ import annotations

import argparse
import re
import shutil
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

ANSWER_CELLS = {1: "B7", 2: "B8", 3: "B9", 4: "B10", 5: "B11", 6: "B12", 7: "B13", 8: "B14", 9: "B15", 10: "B16", 11: "B17"}
MAX_CELL = 32000
SHEET = "xl/worksheets/sheet2.xml"


def plain(text: str) -> str:
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"(?<!\*)\*(?!\*)([^*\n]+?)\*(?!\*)", r"\1", text)
    return text.replace("`", "").strip()


def answers(md: str) -> dict[int, str]:
    parts = re.split(r"^### (\d+)\. .*$", md, flags=re.M)
    out = {}
    for i in range(1, len(parts), 2):
        body = re.sub(r"\n\(Abstract length: \d+ words\.\)\s*$", "", parts[i + 1].strip())
        out[int(parts[i])] = plain(body)
    return out


def fill(form: Path, md: str) -> dict[str, int]:
    ans = answers(md)
    with zipfile.ZipFile(form) as z:
        names = z.namelist()
        data = {n: z.read(n) for n in names}
    xml = data[SHEET].decode("utf-8")
    sizes = {}
    for q, cell in ANSWER_CELLS.items():
        text = ans[q]
        if len(text) > MAX_CELL:
            raise SystemExit(f"answer {q} has {len(text)} characters (Excel limit 32767)")
        pattern = re.compile(r'<c r="%s"([^>]*?)(?:/>|>(?:<is>.*?</is>)?</c>)' % cell, re.S)
        if not pattern.search(xml):
            raise SystemExit(f"cell {cell} not found")

        def repl(m, text=text):
            attrs = re.sub(r'\s*t="[^"]*"', "", m.group(1))
            return '<c r="%s"%s t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>' % (cell, attrs, escape(text))

        xml = pattern.sub(repl, xml, count=1)
        sizes[cell] = len(text)
    data[SHEET] = xml.encode("utf-8")
    tmp = form.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
        for n in names:
            zo.writestr(n, data[n])
    shutil.move(str(tmp), str(form))
    return sizes


def read_cells(form: Path) -> dict[str, str]:
    """Read back the text of cells B7 to B17 of the Track 2 sheet."""
    from html import unescape

    with zipfile.ZipFile(form) as z:
        xml = z.read(SHEET).decode("utf-8")
    out = {}
    for cell in ANSWER_CELLS.values():
        m = re.search(r'<c r="%s"[^>]*>(?:<is><t[^>]*>(.*?)</t></is>)?</c>' % cell, xml, re.S)
        out[cell] = unescape(m.group(1)) if m and m.group(1) is not None else ""
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods", type=Path, default=Path("submissions/candeleroj98_methods_description_track2.md"))
    parser.add_argument("--form", type=Path, default=Path("submissions/candeleroj98_methods_description_form_track2.xlsx"))
    parser.add_argument("--template", type=Path, default=None, help="blank official template: copied to --form before filling")
    args = parser.parse_args()
    if args.template:
        shutil.copyfile(args.template, args.form)
    print("filled", fill(args.form, args.methods.read_text(encoding="utf-8")))


if __name__ == "__main__":
    main()
