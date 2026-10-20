"""Self-contained static review page (no external resources, no server).

    python -m track2.review_page --out out        # writes out/review.html

Reads out/candidates.json and, if present, direction_review.json, evidence_retrieval.json,
faers_disproportionality.json, pgx_drug_relevance.json and reversal_scores.json. A missing optional file is
skipped and said so on the page. Every value is HTML-escaped. Aggregate, non-patient data only.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

OPTIONAL = {
    "direction_review.json": "Direction review",
    "evidence_retrieval.json": "Evidence retrieval",
    "faers_disproportionality.json": "FAERS disproportionality",
    "pgx_drug_relevance.json": "Pharmacogenomic relevance",
    "reversal_scores.json": "Transcriptomic reversal",
}
CSS = """:root{--bg:#fff;--fg:#1b1f23;--mut:#586069;--line:#d0d7de;--ban:#fff4ce;--banfg:#4a3b00;--acc:#0b5cad;--card:#f6f8fa}
@media(prefers-color-scheme:dark){:root{--bg:#0d1117;--fg:#e6edf3;--mut:#9da7b3;--line:#30363d;--ban:#3b2f00;--banfg:#ffe9a0;--acc:#6cb6ff;--card:#161b22}}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif}
main{max-width:1200px;margin:0 auto;padding:0 16px 48px}
.banner{background:var(--ban);color:var(--banfg);padding:12px 16px;font-weight:600;text-align:center;border-bottom:1px solid var(--line)}
.tools{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}
input,select,button{font:inherit;color:var(--fg);background:var(--card);border:1px solid var(--line);border-radius:6px;padding:6px 10px}
table{border-collapse:collapse;width:100%}th,td{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}
th button{border:0;background:none;font-weight:700;padding:0;cursor:pointer}
th[aria-sort=ascending] button::after{content:" \\25B2"}th[aria-sort=descending] button::after{content:" \\25BC"}
.detail td{background:var(--card)}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px}
h3{margin:8px 0 4px;font-size:1em}ul{margin:0;padding-left:18px}.mut{color:var(--mut)}
:focus-visible{outline:2px solid var(--acc);outline-offset:2px}[hidden]{display:none!important}"""
JS = """(function(){var t=document.getElementById('t'),q=document.getElementById('q'),f=document.getElementById('f'),n=document.getElementById('n');
function bodies(){return Array.prototype.slice.call(t.tBodies)}
function filt(){var s=q.value.toLowerCase(),g=f.value,c=0;bodies().forEach(function(b){var ok=(!s||b.dataset.text.indexOf(s)>-1)&&(!g||(' '+b.dataset.flags+' ').indexOf(' '+g+' ')>-1);b.hidden=!ok;if(ok)c++});n.textContent=c+' shown'}
q.addEventListener('input',filt);f.addEventListener('change',filt);
t.querySelectorAll('th button').forEach(function(b){b.addEventListener('click',function(){var th=b.parentNode,i=th.cellIndex,dir=th.getAttribute('aria-sort')==='ascending'?-1:1,num=th.dataset.type==='num';
t.querySelectorAll('th').forEach(function(x){x.removeAttribute('aria-sort')});th.setAttribute('aria-sort',dir===1?'ascending':'descending');
var rows=bodies().sort(function(a,b2){var x=a.rows[0].cells[i].dataset.v,y=b2.rows[0].cells[i].dataset.v;if(num){x=parseFloat(x);y=parseFloat(y);return(isNaN(x)?1e9:x)-(isNaN(y)?1e9:y)||0}return x.localeCompare(y)});
if(dir===-1)rows.reverse();rows.forEach(function(r){t.appendChild(r)})})});
t.addEventListener('click',function(e){var b=e.target.closest('button.tog');if(!b)return;var d=document.getElementById(b.getAttribute('aria-controls')),o=b.getAttribute('aria-expanded')==='true';b.setAttribute('aria-expanded',String(!o));d.hidden=o});
filt()})();"""


def esc(value: Any) -> str:
    return html.escape("" if value is None else str(value))


def load_optional(out: Path) -> tuple[dict[str, Any], list[str]]:
    found, missing = {}, []
    for name in OPTIONAL:
        path = out / name
        if path.exists():
            found[name] = json.loads(path.read_text(encoding="utf-8"))
        else:
            missing.append(name)
    return found, missing


def index_by_drug(data: Any) -> dict[str, Any]:
    """Accept a list of records or a dict (optionally wrapped in rows/drugs/results); key = drug_id or name."""
    if isinstance(data, dict):
        for wrap in ("rows", "drugs", "results", "candidates"):
            if isinstance(data.get(wrap), (list, dict)):
                data = data[wrap]
                break
    if isinstance(data, list):
        return {str(r.get("drug_id") or r.get("name")): r for r in data if isinstance(r, dict) and (r.get("drug_id") or r.get("name"))}
    if isinstance(data, dict):
        return {str(k): v for k, v in data.items()}
    return {}


def _kv(value: Any) -> str:
    if isinstance(value, dict):
        return "<ul>" + "".join(f"<li><b>{esc(k)}</b>: {_kv(v)}</li>" for k, v in value.items()) + "</ul>"
    if isinstance(value, list):
        return "<ul>" + "".join(f"<li>{_kv(v)}</li>" for v in value[:25]) + ("<li class=mut>...</li>" if len(value) > 25 else "") + "</ul>"
    return esc(value)


def _ul(items: list[str]) -> str:
    return "<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>" if items else "<span class=mut>none</span>"


def detail_html(r: dict[str, Any], extras: dict[str, dict[str, Any]]) -> str:
    key_ids = (r["drug_id"], r["name"])
    mech = [f"{esc(m['mechanism_id'])}: rank {m['rank']}, z={m['z']:.2f}, d={m['d']:.2f}, p={m['p_empirical']}, {'passes' if m['passes_primary_threshold'] else 'below'} primary threshold" for m in r["mechanism_scores"]]
    gates = [f"<b>{esc(g['requirement'])}</b>: {esc(g['status'])} <span class=mut>(needs {esc(g['needs'])})</span>" for g in r["gates"]]
    age = r["age_profile"]
    ages = [f"{esc(k)}: {esc(age[k])}" for k in ("neonate", "infant", "child", "adolescent", "adult") if k in age]
    ages.append(f"pediatric label status: {esc(age.get('pediatric_label_status'))} (heuristic, verify by hand)")
    exc = r["excipients"]
    alerts = [esc(a.get("excipient", a) if isinstance(a, dict) else a) for a in exc.get("alerts", [])]
    tox = r["toxicity"]
    if tox.get("available"):
        terms = [f"{esc(t['term'])} ({t['count']:,})" for t in tox["top_terms"][:8]]
        toxh = f"<p>{tox['total_reports']:,} FAERS reports. {esc(tox['caveat'])}</p>" + _ul(terms)
    else:
        toxh = "<span class=mut>not available</span>"
    ev_for = [f"{esc(e.get('type'))}: {esc(e.get('id'))} (citation verified: {esc(e.get('citation_verified'))})" for e in r["evidence_for"]]
    ev_against = [f"{esc(e.get('text'))} <span class=mut>[{esc(e.get('source'))}]</span>" for e in r["evidence_against"]]
    parts = [
        ("Mechanism scores", _ul(mech)),
        ("Flags and direction", _ul([esc(f) for f in r["risk_flags"]] + [f"direction_review: {esc(r['direction_review'])}"])),
        ("Requirement states", _ul(gates)),
        ("Label by age band", _ul(ages)),
        ("Excipient alerts", _ul(alerts) + f"<p class=mut>{len(exc.get('inactive_ingredients', []))} inactive ingredients listed</p>"),
        ("Toxicity summary (descriptive)", toxh),
        ("Evidence for", _ul(ev_for)),
        ("Evidence against", _ul(ev_against)),
    ]
    for fname, title in OPTIONAL.items():
        rec = next((extras[fname][k] for k in key_ids if fname in extras and k in extras[fname]), None)
        if rec is not None:
            parts.append((title, _kv(rec)))
    return '<div class="grid">' + "".join(f"<section><h3>{t}</h3>{b}</section>" for t, b in parts) + "</div>"


def row_html(i: int, r: dict[str, Any], extras: dict[str, dict[str, Any]]) -> str:
    ranks = [m["rank"] for m in r["mechanism_scores"]]
    best = min(ranks) if ranks else ""
    states = [g["status"] for g in r["gates"]]
    adult = r["age_profile"].get("adult", "UNKNOWN")
    flags = " ".join(r["risk_flags"])
    text = f"{r['name']} {r['drug_id']} {' '.join(r['mechanism_of_action'])} {flags}".lower()
    pid = f"d{i}"
    cells = [
        (i, i),
        (r["name"], r["name"]),
        (r["drug_id"], r["drug_id"]),
        (r["combined_score"], f"{r['combined_score']:.2f}"),
        (best, best),
        (flags, ", ".join(r["risk_flags"]) or "-"),
        (adult, adult),
        (states.count("FAIL"), f"{states.count('FAIL')} FAIL / {states.count('NOT_ASSESSABLE')} NOT_ASSESSABLE"),
    ]
    tds = "".join(f'<td data-v="{esc(v)}">{esc(shown)}</td>' for v, shown in cells)
    return (f'<tbody data-text="{esc(text)}" data-flags="{esc(flags)}"><tr>{tds}<td data-v="">'
            f'<button class="tog" type="button" aria-expanded="false" aria-controls="{pid}">Details</button></td></tr>'
            f'<tr class="detail" id="{pid}" hidden><td colspan="9">{detail_html(r, extras)}</td></tr></tbody>')


def render(candidates: list[dict[str, Any]], optional: dict[str, Any], missing: list[str]) -> str:
    extras = {name: index_by_drug(data) for name, data in optional.items()}
    ordered = sorted(candidates, key=lambda r: -r["combined_score"])
    flags = sorted({f for r in ordered for f in r["risk_flags"]})
    heads = [("Rank", "num"), ("Drug", "text"), ("ChEMBL", "text"), ("Score", "num"), ("Best module rank", "num"), ("Flags", "text"), ("Adult label", "text"), ("Requirements", "num")]
    ths = "".join(f'<th scope="col" data-type="{t}"><button type="button">{h}</button></th>' for h, t in heads) + '<th scope="col">Details</th>'
    options = '<option value="">All flags</option>' + "".join(f'<option value="{esc(f)}">{esc(f)}</option>' for f in flags)
    note = f"<p class=mut>Optional inputs not found (section omitted): {esc(', '.join(missing))}.</p>" if missing else ""
    body = "".join(row_html(i, r, extras) for i, r in enumerate(ordered, 1))
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>Track 2 candidate review</title><style>{CSS}</style></head><body>'
            f'<div class="banner" role="note">Hypotheses for follow-up, not medical advice</div><main><h1>Track 2 candidate review</h1>'
            f'<p>{len(ordered)} approved drugs ranked as hypotheses. A requirement state of NOT_ASSESSABLE means the patient data needed were not available. '
            f'Label flags are heuristic keyword scans: verify by hand.</p>{note}'
            f'<div class="tools"><label>Search <input id="q" type="search"></label><label>Flag <select id="f">{options}</select></label>'
            f'<span id="n" role="status" aria-live="polite"></span></div>'
            f'<table id="t"><caption class="mut">Click a column header to sort; use Details to expand a drug.</caption><thead><tr>{ths}</tr></thead>{body}</table>'
            f'</main><script>{JS}</script></body></html>')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    args = parser.parse_args()
    candidates = json.loads((args.out / "candidates.json").read_text(encoding="utf-8"))
    optional, missing = load_optional(args.out)
    target = args.out / "review.html"
    target.write_text(render(candidates, optional, missing), encoding="utf-8", newline="\n")
    print(f"wrote {target} ({len(candidates)} drugs; missing optional: {missing or 'none'})")


if __name__ == "__main__":
    main()
