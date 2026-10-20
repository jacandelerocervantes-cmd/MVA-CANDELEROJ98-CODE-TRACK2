"""Minimal stdlib reader for the first sheet of an .xlsx file (strings and numbers only; no formulas, no dates)."""

from __future__ import annotations

import re
import zipfile
from xml.etree import ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _col(ref: str) -> int:
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_rows(path: str):
    """Yield each row of the first worksheet as a list (None for empty cells)."""
    with zipfile.ZipFile(path) as z:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.parse(z.open("xl/sharedStrings.xml")).getroot().iter(NS + "si"):
                shared.append("".join(t.text or "" for t in si.iter(NS + "t")))
        sheet = sorted(n for n in z.namelist() if n.startswith("xl/worksheets/sheet"))[0]
        for _, el in ET.iterparse(z.open(sheet), events=("end",)):
            if el.tag != NS + "row":
                continue
            row: list = []
            for c in el.iter(NS + "c"):
                i = _col(c.get("r"))
                row.extend([None] * (i + 1 - len(row)))
                v = c.find(NS + "v")
                t = c.get("t")
                if t == "inlineStr":
                    val = "".join(x.text or "" for x in c.iter(NS + "t"))
                elif v is None or v.text is None:
                    val = None
                elif t == "s":
                    val = shared[int(v.text)]
                elif t in ("str", "e", "b"):
                    val = v.text
                else:
                    val = float(v.text)
                row[i] = val
            yield row
            el.clear()
