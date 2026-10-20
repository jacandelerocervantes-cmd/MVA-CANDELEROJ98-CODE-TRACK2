"""Markdown report -> PDF with no third-party Python package: a small converter to HTML and headless Edge (Windows) or Chromium to print it.

    python -m track2.md_to_pdf submissions/candeleroj98_track2_report.md submissions/candeleroj98_track2_report.pdf

Supports what the report uses: headings, paragraphs, bullet and numbered lists, tables, block quotes, bold, italic and inline code. No network is used.
"""

from __future__ import annotations

import argparse
import html
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

CSS = """
@page { size: letter; margin: 18mm 16mm; }
body { font-family: 'Times New Roman', Times, serif; font-size: 10pt; line-height: 1.4; color: #111; }
h1 { font-size: 16pt; text-align: center; border-bottom: 2px solid #222; padding-bottom: 6px; }
h2 { font-size: 12.5pt; border-bottom: 1px solid #555; padding-bottom: 3px; margin-top: 16px; page-break-after: avoid; }
h3 { font-size: 11pt; margin-top: 12px; page-break-after: avoid; }
p { margin: 0 0 7px 0; text-align: justify; }
ul, ol { margin: 0 0 7px 0; padding-left: 22px; }
li { margin-bottom: 3px; }
blockquote { border-left: 3px solid #333; background: #f7f7f7; margin: 8px 0; padding: 6px 12px; font-size: 9.5pt; }
table { width: 100%; border-collapse: collapse; margin: 10px 0; font-size: 7.8pt; table-layout: auto; page-break-inside: auto; }
th, td { border: 1px solid #444; padding: 3px 5px; vertical-align: top; overflow-wrap: anywhere; }
th { background: #eee; }
table.wide { table-layout: fixed; }
table.wide th, table.wide td { overflow-wrap: normal; word-break: normal; hyphens: manual; }
table.wide td:nth-child(3) { overflow-wrap: anywhere; }
table.wide tr { page-break-inside: avoid; }
code { font-family: Consolas, monospace; font-size: 8.5pt; background: #f2f2f2; padding: 0 2px; }
"""


def inline(text: str) -> str:
    codes: list[str] = []

    def keep(m: re.Match[str]) -> str:
        codes.append(m.group(1))
        return f"\u0000{len(codes) - 1}\u0000"

    text = re.sub(r"`([^`]+)`", keep, text)
    text = html.escape(text, quote=False).replace("&lt;br&gt;", "<br>")
    text = re.sub(r"\[([^\]]+)\]\((https?://[^)\s]+)\)", r'<a href="\2">\1</a>', text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<![\w*])\*(?!\*)([^*\n]+?)\*(?![\w*])", r"<i>\1</i>", text)
    return re.sub(r"\u0000(\d+)\u0000", lambda m: "<code>" + html.escape(codes[int(m.group(1))], quote=False) + "</code>", text)


def convert(md: str) -> str:
    out: list[str] = []
    lines = md.splitlines()
    i = 0
    para: list[str] = []
    widths: list[int] | None = None

    def flush() -> None:
        if para:
            out.append("<p>" + inline(" ".join(para)) + "</p>")
            para.clear()

    while i < len(lines):
        line = lines[i]
        if not line.strip():
            flush()
            i += 1
        elif re.match(r"^(?:\[\[widths:[\d,\s]+\]\]|<!--\s*widths:[\d,\s]+-->)\s*$", line):
            flush()
            widths = [int(x) for x in re.findall(r"\d+", line)]
            i += 1
        elif line.startswith("#"):
            flush()
            level = len(line) - len(line.lstrip("#"))
            out.append(f"<h{level}>{inline(line[level:].strip())}</h{level}>")
            i += 1
        elif line.startswith("|") and i + 1 < len(lines) and re.match(r"^\|[\s:|-]+\|\s*$", lines[i + 1]):
            flush()
            head = [c.strip() for c in line.strip().strip("|").split("|")]
            rows = []
            i += 2
            while i < len(lines) and lines[i].startswith("|"):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            wide = widths is not None and len(widths) == len(head)  # a `[[widths:3,9,...]]` line before a table gives fixed column widths (percent) so short columns do not split words
            cols = "<colgroup>" + "".join(f'<col style="width:{w}%">' for w in widths) + "</colgroup>" if wide else ""
            widths = None
            out.append(f'<table class="wide">{cols}<thead><tr>' if wide else "<table><thead><tr>")
            out[-1] += ("".join(f"<th>{inline(c)}</th>" for c in head) + "</tr></thead><tbody>"
                       + "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in rows) + "</tbody></table>")
        elif re.match(r"^\s*[-*] ", line):
            flush()
            items = []
            while i < len(lines) and re.match(r"^\s*[-*] ", lines[i]):
                items.append(re.sub(r"^\s*[-*] ", "", lines[i]))
                i += 1
            out.append("<ul>" + "".join(f"<li>{inline(x)}</li>" for x in items) + "</ul>")
        elif re.match(r"^\d+\. ", line):
            flush()
            items = []
            while i < len(lines) and re.match(r"^\d+\. ", lines[i]):
                items.append(re.sub(r"^\d+\. ", "", lines[i]))
                i += 1
            out.append("<ol>" + "".join(f"<li>{inline(x)}</li>" for x in items) + "</ol>")
        elif line.startswith(">"):
            flush()
            quote = []
            while i < len(lines) and lines[i].startswith(">"):
                quote.append(lines[i].lstrip("> ").strip())
                i += 1
            out.append("<blockquote><p>" + inline(" ".join(quote)) + "</p></blockquote>")
        else:
            para.append(line.strip())
            i += 1
    flush()
    return "\n".join(out)


def find_browser() -> str:
    for p in (r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe", r"C:\Program Files\Microsoft\Edge\Application\msedge.exe", r"C:\Program Files\Google\Chrome\Application\chrome.exe"):
        if Path(p).exists():
            return p
    for name in ("chromium", "chromium-browser", "google-chrome", "microsoft-edge"):
        if shutil.which(name):
            return str(shutil.which(name))
    raise SystemExit("no Edge, Chrome or Chromium found to print the PDF")


def to_pdf(md_path: Path, pdf_path: Path, landscape: bool = False, lang: str = "en") -> None:
    title = re.search(r"^# (.+)$", md_path.read_text(encoding="utf-8"), re.M)
    body = convert(md_path.read_text(encoding="utf-8"))
    css = CSS.replace("size: letter;", "size: letter landscape;") if landscape else CSS
    doc = f"<!DOCTYPE html><html lang=\"{lang}\"><head><meta charset=\"utf-8\"><title>{html.escape(title.group(1) if title else md_path.stem)}</title><style>{css}</style></head><body>{body}</body></html>"
    with tempfile.TemporaryDirectory() as tmp:
        page = Path(tmp) / "report.html"
        page.write_text(doc, encoding="utf-8")
        pdf_path.unlink(missing_ok=True)
        subprocess.run([find_browser(), "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--user-data-dir={tmp}/profile", f"--print-to-pdf={pdf_path.resolve()}", page.as_uri()],
                       check=True, capture_output=True, timeout=180)
        # the browser can return before the file is complete: wait until it exists and its size stops changing
        last, stable = -1, 0
        for _ in range(180):
            size = pdf_path.stat().st_size if pdf_path.exists() else -1
            stable = stable + 1 if size == last and size > 0 else 0
            if stable >= 2:
                break
            last = size
            time.sleep(0.5)
        else:
            raise SystemExit("the PDF was not written")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("markdown", type=Path)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--landscape", action="store_true")
    parser.add_argument("--lang", default="en")
    args = parser.parse_args()
    to_pdf(args.markdown, args.pdf, args.landscape, args.lang)
    print("wrote", args.pdf, args.pdf.stat().st_size, "bytes")


if __name__ == "__main__":
    main()
