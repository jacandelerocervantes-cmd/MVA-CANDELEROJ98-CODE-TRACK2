"""Reference set of published BUB1B (MVA1) cases from Europe PMC, with code-verified verbatim sentences.

Every frequency stored in out/case_reference_bub1b.json comes from a sentence that was found verbatim
(after whitespace and entity normalisation) in the text fetched from Europe PMC; the number is also
required to occur inside that sentence. A feature that none of the sources quantifies is stored as
"not reported in the sources". Network errors raise; they are never read as missing data.
"""

from __future__ import annotations

import html
import json
import re
import time
import urllib.request
from datetime import date
from pathlib import Path
from typing import Callable

REST = "https://www.ebi.ac.uk/europepmc/webservices/rest"
NOT_REPORTED = "not reported in the sources"

SOURCES = {
    "35804254": {"pmcid": "PMC9616775", "kind": "full text", "label": "Systematic review of MVA1 cases (Pavone et al. 2022)"},
    "18548531": {"pmcid": None, "kind": "abstract", "label": "Clinical and genetic heterogeneity in MVA (Garcia-Castillo et al. 2008)"},
    "20516114": {"pmcid": "PMC2887387", "kind": "abstract", "label": "Molecular causes for BUBR1 dysfunction in MVA (Suijkerbuijk et al. 2010)"},
    "39264246": {"pmcid": "PMC11705613", "kind": "abstract", "label": "AACR update on genomic instability disorders (2024)"},
}

# feature -> list of claims. A claim is a verbatim sentence, its PMID, the frequency token that must
# occur inside the sentence ("value") and the population it refers to. value=None means the sentence
# names the feature without any frequency.
_REVIEW = "35804254"
_GC = "18548531"
S_GROWTH = "In the group of 30 MVA patients, pre-postnatal growth retardation, microcephaly, ID/DD, and epileptic seizures occurred in 93.5%, 87%, 77.4%, and 45.1% of the cases, respectively."
S_DWM = "Furthermore, DWM was present in 25.8% of the cases and other cerebral anomalies including the corpus callosum hypo-aplasia (CAA) were reported in 32.2% of the cases."
S_RENAL = "Less frequent but serious presentations include genital anomalies (12.9%), renal dysfunction (12.9%), and endocrine dysfunction (12.9%)."
S_MALIG = "Malignancies had noticeably a high incidence (38.7%) with Wilms tumor, rhabdomyosarcoma and leukemia being the most frequently reported and recently polycystic ovary syndrome"
S_LESS = "As suggested by this study, less frequent features included minor facial features, bilateral cataracts, microphthalmia, cleft palate, ambiguous genitalia in males, skin abnormalities, and malignancies."
S_IUGR = "MVA1 is clinically characterized by intrauterine growth restriction, post-natal growth retardation, and severe neurologic impairment including microcephaly, developmental delay/intellectual disability, epileptic seizures, and generalized hypotonia."
S_GC = "Patients with monoallelic BUB1B mutations were severely affected with Dandy-Walker complex (7/8), cataracts (6/6), and Wilms' tumor (7/8); premature chromatid separation (PCS) was observed in 8/8 propositi and 7/7 carrier parents."
S_GENO = "In patients with biallelic mutations, a missense mutation pairs with a truncating mutation."
S_KIN = "In this group of missense mutants, the amino acid change consistently occurs in or near the BUBR1 kinase domain."
S_HOMCH = "MVA1 type 1 (OMIM 257,300) is caused by homozygous or compound heterozygous mutations in the BUB1B gene (OMIM 602,860) located in chromosome 15q1 and is clinically characterized by microcephaly, developmental delay/intellectual disability (DD/ID), epileptic seizures, generalized hypotonia, intrauterine growth restriction (IUGR), and postnatal growth retardation."
S_DOMAIN = "In the protein structure can be identified an N-terminal region required for BUB1B binding, an intermediate region containing the BUB3 and CDC20 binding domains, and a catalytic serine/threonine kinase domain in the C-terminal region."

POP_REVIEW = "published MVA cases in the systematic review (text: 'the group of 30 MVA patients'; percentages as published)"
POP_GC = "patients with BUB1B mutations screened in the 2008 review (denominator as published)"

FEATURES: list[dict] = [
    {"feature": "microcephaly", "hpo": ["HP:0000252"], "claims": [{"pmid": _REVIEW, "sentence": S_GROWTH, "value": "87%", "population": POP_REVIEW}]},
    {"feature": "pre- or postnatal growth restriction", "hpo": ["HP:0001510", "HP:0004322", "HP:0001508", "HP:0001518", "HP:0001511"], "claims": [{"pmid": _REVIEW, "sentence": S_GROWTH, "value": "93.5%", "population": POP_REVIEW}]},
    {"feature": "intrauterine growth restriction (separately)", "hpo": ["HP:0001511"], "claims": [{"pmid": _REVIEW, "sentence": S_IUGR, "value": None, "population": "named as a characteristic feature without a separate frequency"}]},
    {"feature": "developmental delay / intellectual disability", "hpo": ["HP:0001263"], "claims": [{"pmid": _REVIEW, "sentence": S_GROWTH, "value": "77.4%", "population": POP_REVIEW}]},
    {"feature": "epileptic seizures", "hpo": ["HP:0001250"], "claims": [{"pmid": _REVIEW, "sentence": S_GROWTH, "value": "45.1%", "population": POP_REVIEW}]},
    {"feature": "Dandy-Walker malformation", "hpo": ["HP:0001305"], "claims": [
        {"pmid": _REVIEW, "sentence": S_DWM, "value": "25.8%", "population": POP_REVIEW},
        {"pmid": _GC, "sentence": S_GC, "value": "7/8", "population": POP_GC}]},
    {"feature": "other brain malformations (incl. corpus callosum hypo-aplasia)", "hpo": ["HP:0012443"], "claims": [{"pmid": _REVIEW, "sentence": S_DWM, "value": "32.2%", "population": POP_REVIEW}]},
    {"feature": "any malignancy", "hpo": ["HP:0002664"], "claims": [{"pmid": _REVIEW, "sentence": S_MALIG, "value": "38.7%", "population": POP_REVIEW}]},
    {"feature": "rhabdomyosarcoma", "hpo": ["HP:0002859"], "claims": [{"pmid": _REVIEW, "sentence": S_MALIG, "value": None, "population": "named among the most frequently reported malignancies; no individual frequency"}]},
    {"feature": "Wilms tumour", "hpo": ["HP:0002667"], "claims": [
        {"pmid": _REVIEW, "sentence": S_MALIG, "value": None, "population": "named among the most frequently reported malignancies; no individual frequency"},
        {"pmid": _GC, "sentence": S_GC, "value": "7/8", "population": POP_GC}]},
    {"feature": "leukaemia", "hpo": ["HP:0001909"], "claims": [{"pmid": _REVIEW, "sentence": S_MALIG, "value": None, "population": "named among the most frequently reported malignancies; no individual frequency"}]},
    {"feature": "renal dysfunction", "hpo": ["HP:0000077"], "claims": [{"pmid": _REVIEW, "sentence": S_RENAL, "value": "12.9%", "population": POP_REVIEW}]},
    {"feature": "nephrocalcinosis", "hpo": ["HP:0000121"], "claims": []},
    {"feature": "cataract", "hpo": ["HP:0000518"], "claims": [
        {"pmid": _GC, "sentence": S_GC, "value": "6/6", "population": POP_GC},
        {"pmid": _REVIEW, "sentence": S_LESS, "value": None, "population": "named among less frequent features; no frequency"}]},
    {"feature": "microphthalmia", "hpo": ["HP:0000568"], "claims": [{"pmid": _REVIEW, "sentence": S_LESS, "value": None, "population": "named among less frequent features; no frequency"}]},
]

GENOTYPE_CLAIMS = [
    {"pattern": "biallelic with one truncating and one missense variant", "pmid": "20516114", "sentence": S_GENO, "value": None},
    {"pattern": "missense variants lie in or near the kinase domain (protein level)", "pmid": "20516114", "sentence": S_KIN, "value": None},
    {"pattern": "homozygous or compound heterozygous", "pmid": _REVIEW, "sentence": S_HOMCH, "value": None},
    {"pattern": "protein domain organisation (N-terminal kinetochore-binding region, BUB3/CDC20-binding region, C-terminal kinase domain)", "pmid": _REVIEW, "sentence": S_DOMAIN, "value": None},
]
GENOTYPE_NOT_REPORTED = {
    "two truncating variants (frequency)": NOT_REPORTED,
    "counts of truncating plus missense versus two truncating": NOT_REPORTED,
}

# Absence checks (regex over each fetched text) for child features that no source mentions.
ABSENCE_PROBES = {
    "nephrocalcinosis": r"nephrocalcinosis",
    "premature birth": r"preterm|premature birth|prematurity",
    "skeletal muscle atrophy": r"muscle atrophy|muscular atrophy|sarcopenia",
}


def normalise(text: str) -> str:
    text = html.unescape(re.sub(r"<[^>]+>", " ", text))
    text = text.replace("−", "-").replace("–", "-").replace("’", "'").replace(" ", " ")
    return re.sub(r"\s+", " ", text).strip()


def verify_claim(text: str, sentence: str, value: str | None) -> bool:
    """True when the sentence occurs verbatim in the normalised text and contains the value token."""
    norm = normalise(text)
    return sentence in norm and (value is None or value in sentence)


def _get(url: str, retries: int = 5) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "hselith-track2/0.1"}), timeout=60) as r:
                data = r.read()
            time.sleep(0.3)  # at most ~3 requests/s
            return data
        except Exception as exc:  # noqa: BLE001 - retried then re-raised
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Europe PMC request failed: {url}: {last}")


def fetch_texts() -> dict[str, dict]:
    """PMID -> {'text', 'kind', 'retrieved'|'error'} (full text where open access, else abstract)."""
    out: dict[str, dict] = {}
    for pmid, meta in SOURCES.items():
        entry: dict = {"kind": None, "text": ""}
        if meta["kind"] == "full text":
            entry = {"kind": "full text", "text": normalise(_get(f"{REST}/{meta['pmcid']}/fullTextXML").decode("utf-8"))}
        else:
            res = json.loads(_get(f"{REST}/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json&resultType=core"))
            rec = res["resultList"]["result"][0]
            entry = {"kind": "abstract", "text": normalise(rec.get("abstractText") or ""), "title": normalise(rec.get("title") or "")}
            try:  # full text of the AACR chapter is not open access: record the outcome only
                if pmid == "39264246":
                    _get(f"{REST}/{meta['pmcid']}/fullTextXML", retries=2)
                    entry["full_text_note"] = "full text endpoint answered; abstract only was used"
            except RuntimeError:
                entry["full_text_note"] = "full text not retrievable from Europe PMC (not open access); abstract only"
        out[pmid] = entry
    return out


def build_reference(texts: dict[str, dict], when: str | None = None) -> dict:
    """Assemble the reference set; raises if any sentence or number is not verified in the fetched text."""
    def claim_out(c: dict) -> dict:
        t = texts[c["pmid"]]["text"]
        if not verify_claim(t, c["sentence"], c.get("value")):
            raise ValueError(f"sentence not verified for PMID {c['pmid']}: {c['sentence'][:70]}")
        d = {k: c[k] for k in c if k != "sentence"}
        d.update({"sentence": c["sentence"], "verified_verbatim": True, "source_kind": texts[c["pmid"]]["kind"]})
        return d
    feats = []
    for f in FEATURES:
        claims = [claim_out(c) for c in f["claims"]]
        feats.append({"feature": f["feature"], "hpo_ids": f["hpo"], "claims": claims,
                      "frequency_status": "quantified" if any(c["value"] for c in claims)
                      else "named without frequency" if claims else NOT_REPORTED})
    absence = {k: {"pattern": rx, "found_in_pmids": [p for p, e in texts.items() if re.search(rx, e["text"], re.I)]}
               for k, rx in ABSENCE_PROBES.items()}
    return {
        "disease": "mosaic variegated aneuploidy type 1 (BUB1B)",
        "retrieved": when or date.today().isoformat(),
        "sources": {p: {**SOURCES[p], "text_used": texts[p]["kind"], "note": texts[p].get("full_text_note")} for p in SOURCES},
        "features": feats,
        "genotype_claims": [claim_out(c) for c in GENOTYPE_CLAIMS],
        "genotype_not_reported": GENOTYPE_NOT_REPORTED,
        "absence_probes": absence,
        "notes": [
            "PMID 18407612 given in the task resolves to an unrelated article; the Garcia-Castillo MVA review is PMID 18548531 (abstract only; no PMC full text).",
            "The AACR 2024 genomic instability update (PMID 39264246) full text is not retrievable from Europe PMC; its abstract carries no MVA frequencies.",
            "The group label 'monoallelic BUB1B mutations' is quoted as published in PMID 18548531.",
            "The systematic review (PMID 35804254) says 'group of 30 MVA patients' for its percentages; frequencies are reproduced as published and not recomputed.",
        ],
    }


def main(root: Path = Path("."), fetch: Callable[[], dict] = fetch_texts) -> Path:
    ref = build_reference(fetch())
    out = root / "out/case_reference_bub1b.json"
    out.write_text(json.dumps(ref, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    return out


if __name__ == "__main__":
    print(main())
