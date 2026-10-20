"""Independent re-check of every sentence of the literature search (`out/literature_hypotheses/*.json`), written by the research tasks that produced them.

    python -m track2.verify_literature

Each quoted sentence is looked for in the text of the cited source, fetched again from Europe PMC (the abstract, or the full text when the entry says it came from the full text; NCBI efetch is the
fallback for full text). Whitespace is ignored on both sides and HTML entities are decoded; nothing else is normalised. A source whose text cannot be fetched is reported as `not_checkable`, never as
verified. Writes `out/literature_hypotheses/verification.json`; exit code 1 if any sentence is not found.
"""

from __future__ import annotations

import html
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from track2.http_cache import get_json

TOPICS = ["growth_muscle", "proteotoxic_autophagy", "senescence_p53", "tumour_cin_safety", "context", "round2_models", "round2_clinical", "round2_phenotypes", "round2_vulnerabilities", "round3_audit"]
SECTIONS = ("findings", "background_evidence", "against_or_risk", "not_approved")
UA = {"User-Agent": "hselith-track2/0.1 (research)"}


def squash(text: str) -> str:
    text = re.sub(r"</?[A-Za-z][^<>]*>", " ", text)  # only real tags: a bare "<" in the text (for example "p<0.001") must not swallow what follows
    return re.sub(r"\s+", "", html.unescape(text).replace("\xa0", " "))


def _get(url: str, timeout: int = 40) -> bytes | None:
    try:
        return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout).read()  # noqa: S310 - fixed https endpoints
    except (urllib.error.URLError, TimeoutError, OSError):
        return None


def abstract_text(source_id: str) -> str | None:
    """Abstract (and title) from Europe PMC for a PubMed id or a preprint id (PPR...)."""
    src = "PPR" if source_id.startswith("PPR") else "MED"
    q = urllib.parse.quote(f"EXT_ID:{source_id} AND SRC:{src}")
    try:
        data = get_json(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&resultType=core&pageSize=1")
    except Exception:  # noqa: BLE001 - not checkable, reported as such
        return None
    res = (data.get("resultList") or {}).get("result") or []
    if not res:
        return None
    return " ".join([res[0].get("title") or "", res[0].get("abstractText") or ""])


def fulltext_text(pmcid: str) -> str | None:
    pmcid = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
    raw = _get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML")
    if raw is None:
        raw = _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid[3:]}&retmode=xml")
    return raw.decode("utf-8", errors="replace") if raw else None


def entries(directory: Path) -> list[dict[str, Any]]:
    out = []
    for topic in TOPICS:
        data = json.loads((directory / f"{topic}.json").read_text(encoding="utf-8"))
        for section in SECTIONS:
            for i, x in enumerate(data.get(section, [])):
                if isinstance(x, dict) and x.get("sentence"):
                    out.append({"topic": topic, "section": section, "index": i, "drug": x.get("drug"), "pmid": str(x.get("pmid") or ""), "pmcid": str(x.get("pmcid") or ""),
                                "source_section": x.get("source_section") or "abstract", "sentence": x["sentence"], "claimed_verified": bool(x.get("verified"))})
    return out


def check(items: list[dict[str, Any]]) -> dict[str, Any]:
    abstracts: dict[str, str | None] = {}
    fulls: dict[str, str | None] = {}
    ids = sorted({i["pmid"] for i in items if i["pmid"]})
    with ThreadPoolExecutor(max_workers=8) as pool:
        for pid, text in zip(ids, pool.map(abstract_text, ids)):
            abstracts[pid] = text
        pmcs = sorted({i["pmcid"] for i in items if i["pmcid"] and i["source_section"] == "full_text"})
        for pc, text in zip(pmcs, pool.map(fulltext_text, pmcs)):
            fulls[pc] = text
    rows = []
    for it in items:
        texts = []
        if it["pmid"] and abstracts.get(it["pmid"]):
            texts.append(squash(abstracts[it["pmid"]] or ""))
        if it["pmcid"] and fulls.get(it["pmcid"]):
            texts.append(squash(fulls[it["pmcid"]] or ""))
        if not texts:
            status = "not_checkable"
        else:
            status = "verified" if any(squash(it["sentence"]) in t for t in texts) else "not_found"
        rows.append({**{k: it[k] for k in ("topic", "section", "index", "drug", "pmid", "pmcid", "source_section", "claimed_verified")}, "status": status, "sentence": it["sentence"][:160]})
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    return {"n_sentences": len(rows), "counts": counts, "rows": rows,
            "method": "abstract or full text re-fetched from Europe PMC (NCBI efetch as fallback); whitespace ignored, HTML entities decoded, nothing else normalised"}


def main() -> None:
    directory = Path("out/literature_hypotheses")
    result = check(entries(directory))
    (directory / "verification.json").write_text(json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(result["n_sentences"], "sentences;", result["counts"])
    for r in result["rows"]:
        if r["status"] == "not_found":
            print("NOT FOUND:", r["topic"], r["section"], r["drug"], r["pmid"], "|", r["sentence"][:90])
    sys.exit(1 if result["counts"].get("not_found") else 0)


if __name__ == "__main__":
    main()
