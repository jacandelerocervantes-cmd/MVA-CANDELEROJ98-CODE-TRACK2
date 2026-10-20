"""Code-derived, disease-agnostic landscape of the genes associated with a rare disease.

Run:  PYTHONPATH=src python -m track2.disease_landscape mva|fanconi_anaemia
      PYTHONPATH=src python -m track2.disease_landscape custom --terms "term one" "term two" --ids OMIM:123456 ORPHA:99 --out out/x.json

A disease is a `DiseaseSpec`: search terms and/or OMIM/ORPHA identifiers (nothing else is disease specific; the
mechanism rule table and the evidence logic are shared by every disease).

Gene list (structured sources, no hand-made list)
  * HPO (release read from the phenotype.hpoa header): diseases of phenotype.hpoa whose name matches a search term, or
    whose identifier is listed; gene symbols from genes_to_phenotype.txt (disease id -> gene); mode of inheritance =
    aspect-I HPO term(s) annotated to the disease, labelled from hp.obo. Files are looked up under data/sources, then data/.
  * Orphanet (Orphadata en_product6 gene-disease, en_product9_ages inheritance; downloaded, cached) for the disorders
    whose name matches a term or whose ORPHA code is listed. If Orphadata is unreachable the run continues with HPO only
    and records the failure in the header.
  Genes of the two sources are merged by HGNC symbol (an Orphanet symbol that is an Orphanet synonym of an HPO symbol
  is merged into the HPO symbol). Per gene: symbol, disease names/identifiers, subtype label (parsed generically from
  'complementation group X' / 'syndrome N' / 'type N'), inheritance as annotated, sources listing it.

Mechanism class (deterministic; public functional annotation of the reviewed human UniProt entry: keywords + Gene
Ontology term names, fetched from UniProt REST because no GO/Reactome file is under data/sources). The FIRST rule whose
regular expression matches any term wins; the order is fixed:
  1 DNA interstrand-crosslink repair      'interstrand cross-link repair'
  2 spindle-assembly checkpoint           'spindle assembly checkpoint', 'spindle checkpoint', 'mitotic checkpoint'
  3 DNA repair/replication                'DNA repair', 'DNA replication', 'double-strand break repair',
                                          'homologous recombination'
  4 kinetochore/microtubule attachment    'kinetochore', 'attachment of spindle microtubules'
  5 centrosome/centriole                  'centrosom', 'centriole'
  6 sister-chromatid cohesion             'sister chromatid cohesion', 'cohesin'
  7 DNA repair/replication (weak term)    'DNA damage'
  8 other                                 no rule matched, or no UniProt entry
Specific rules precede the generic ones (e.g. DNA repair before centrosome/cohesion, because BRCA1/BRCA2/RAD51C also carry
centrosome or cohesion terms). The table contains no gene or disease name, so it needs no edit for another disease. Classes are coarse
functional annotation, not a statement about disease mechanism.

Association evidence: up to one verbatim sentence per gene from Europe PMC (abstracts of up to 25 hits for
'gene AND disease term', then up to 3 open-access PMC full texts). The sentence must name the gene (symbol or
Orphanet synonym) and a disease term/abbreviation in the same sentence. It is verified by code as a
whitespace-normalised substring of a freshly fetched abstract or full text; ValueError is raised if a kept sentence is
not found. If no sentence verifies the gene is marked 'unverified'.

Track 1 cross-reference (only when spec.track1_crossref): gene symbols from out/submission_genes.json, read-only;
no coordinates or variants are read into the output.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data/cache/disease_landscape"
SUBMISSION_GENES = ROOT / "out/submission_genes.json"
ORPHADATA = {
    "gene_disease": "https://www.orphadata.com/data/xml/en_product6.xml",
    "inheritance": "https://www.orphadata.com/data/xml/en_product9_ages.xml",
}
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
HEADER_STATEMENT = "derived by code from public sources; association evidence is that of the cited sources"


# ---------------------------------------------------------------- disease specification
def term_regex(term: str) -> str:
    """Regex for a search term: case-insensitive, space/hyphen tolerant, anemia/anaemia tolerant."""
    parts = [re.escape(w) for w in re.split(r"[\s-]+", term.strip()) if w]
    return r"[\s-]+".join(re.sub(r"an(a)?emia", "an(a)?emia", p, flags=re.I) if re.search("anaemia|anemia", p, re.I) else p for p in parts)


@dataclass(frozen=True)
class DiseaseSpec:
    key: str
    label: str
    terms: tuple[str, ...]  # disease search terms (name match in HPO / Orphanet / Europe PMC)
    ids: tuple[str, ...] = ()  # optional 'OMIM:257300' / 'ORPHA:84' identifiers
    abbreviations: tuple[str, ...] = ()  # case-sensitive short forms accepted in an evidence sentence
    out: Path = ROOT / "out/landscape.json"
    track1_crossref: bool = False

    @property
    def name_re(self) -> re.Pattern[str]:
        return re.compile("|".join(term_regex(t) for t in self.terms), re.I)

    @property
    def assoc_re(self) -> re.Pattern[str]:
        pats = [term_regex(t) for t in self.terms] + [rf"(?-i:\b{re.escape(a)}\b)" for a in self.abbreviations]
        return re.compile("|".join(pats), re.I)


SPECS = {
    "mva": DiseaseSpec(
        "mva",
        "mosaic variegated aneuploidy / premature chromatid separation syndrome",
        ("mosaic variegated aneuploidy", "premature chromatid separation"),
        abbreviations=("MVA", "PCS"),
        out=ROOT / "out/mva_landscape.json",
        track1_crossref=True,
    ),
    "fanconi_anaemia": DiseaseSpec(
        "fanconi_anaemia",
        "Fanconi anaemia",
        ("Fanconi anaemia", "Fanconi anemia"),
        abbreviations=("FA",),
        out=ROOT / "out/landscape_fanconi_anaemia.json",
    ),
}


# ---------------------------------------------------------------- network (cached, throttled, retried)
def _fetch(url: str, retries: int = 6, timeout: int = 180) -> bytes:
    path = CACHE / (hashlib.sha256(url.encode()).hexdigest()[:24] + ".bin")
    if path.exists():
        return path.read_bytes()
    last: Exception | None = None
    for attempt in range(retries):
        time.sleep(0.25)  # <= 4 requests per second
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "hselith-track2/0.2 (research)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310 - fixed public endpoints
                data = r.read()
            CACHE.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            return data
        except Exception as error:  # noqa: BLE001 - retry (SSL, transient)
            last = error
            time.sleep(2 + 4 * attempt)
    raise RuntimeError(f"No response after {retries} attempts: {url} ({last!r})")


# ---------------------------------------------------------------- HPO
def find_source(name: str, roots: tuple[Path, ...] = (ROOT / "data/sources", ROOT / "data")) -> Path:
    for base in roots:
        hits = sorted(base.rglob(name)) if base.exists() else []
        if hits:
            return hits[0]
    raise FileNotFoundError(f"{name} not found under {[str(r) for r in roots]}")


def hpo_labels(obo: Path, wanted: set[str]) -> dict[str, str]:
    names, cur = {}, None
    for line in obo.read_text(encoding="utf-8").splitlines():
        if line.startswith("id: "):
            cur = line[4:].strip()
        elif line.startswith("name: ") and cur in wanted:
            names[cur] = line[6:].strip()
    return names


def hpo_diseases(spec: DiseaseSpec, hpoa: Path, obo: Path) -> dict[str, dict[str, Any]]:
    """Diseases of phenotype.hpoa matching the spec by name or identifier, with aspect-I inheritance labels."""
    rx, ids = spec.name_re, set(spec.ids)
    diseases: dict[str, dict[str, Any]] = {}
    for line in hpoa.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or line.startswith("database_id") or not line.strip():
            continue
        f = line.split("\t")
        if rx.search(f[1]) or f[0] in ids:
            d = diseases.setdefault(f[0], {"name": f[1], "inheritance_ids": set()})
            if len(f) > 10 and f[10] == "I":
                d["inheritance_ids"].add(f[3])
    labels = hpo_labels(obo, {i for d in diseases.values() for i in d["inheritance_ids"]})
    for d in diseases.values():
        d["inheritance"] = sorted(labels.get(i, i) for i in d.pop("inheritance_ids"))
    return diseases


def hpo_gene_map(g2p: Path, disease_ids: set[str]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for line in g2p.read_text(encoding="utf-8").splitlines()[1:]:
        f = line.split("\t")
        if len(f) > 5 and f[5] in disease_ids:
            out[f[5]].add(f[1])
    return out


def hpo_release(hpoa: Path) -> dict[str, str]:
    meta = {}
    for line in hpoa.read_text(encoding="utf-8").splitlines()[:8]:
        for key in ("version", "hpo-version"):
            if line.startswith(f"#{key}:"):
                meta[key] = line.split(":", 1)[1].strip()
    return meta


SUBTYPE_RES = (
    re.compile(r"complementation group ([A-Z0-9]+)", re.I),
    re.compile(r"\bsyndrome (\d+[A-Z]?)\b", re.I),
    re.compile(r"\btype (\d+[A-Z]?)\b", re.I),
)


def subtype_of(name: str) -> str | None:
    """Generic subtype label parsed from a disease name, e.g. 'complementation group C' -> 'group C', 'syndrome 2' -> 'type 2'."""
    for i, rx in enumerate(SUBTYPE_RES):
        m = rx.search(name)
        if m:
            return f"group {m.group(1).upper()}" if i == 0 else f"type {m.group(1)}"
    return None


# ---------------------------------------------------------------- Orphanet
def _t(el: ET.Element | None, path: str) -> str | None:
    x = el.find(path) if el is not None else None
    return x.text.strip() if x is not None and x.text else None


def parse_orphanet(gene_xml: bytes, inh_xml: bytes | None, spec: DiseaseSpec) -> list[dict[str, Any]]:
    rx, ids = spec.name_re, {i.split(":", 1)[1] for i in spec.ids if i.startswith("ORPHA:")}
    inheritance: dict[str, str] = {}
    if inh_xml:
        for d in ET.fromstring(inh_xml).iter("Disorder"):
            code = _t(d, "OrphaCode")
            inheritance[code] = "; ".join(x.text.strip() for x in d.findall("./TypeOfInheritanceList/TypeOfInheritance/Name") if x.text) or "not annotated"
    recs: list[dict[str, Any]] = []
    for d in ET.fromstring(gene_xml).iter("Disorder"):
        code, name = _t(d, "OrphaCode"), _t(d, "Name") or ""
        if not (rx.search(name) or code in ids):
            continue
        for a in d.findall("./DisorderGeneAssociationList/DisorderGeneAssociation"):
            g = a.find("Gene")
            hgnc = [r.findtext("Reference") for r in g.findall("./ExternalReferenceList/ExternalReference") if r.findtext("Source") == "HGNC"]
            recs.append(
                {
                    "symbol": _t(g, "Symbol"),
                    "synonyms": sorted(s.text.strip() for s in g.findall("./SynonymList/Synonym") if s.text),
                    "hgnc_id": f"HGNC:{hgnc[0]}" if hgnc else None,
                    "orpha": f"ORPHA:{code}",
                    "disease_name": name,
                    "association_type": _t(a, "DisorderGeneAssociationType/Name"),
                    "association_status": _t(a, "DisorderGeneAssociationStatus/Name"),
                    "orphanet_pmids": re.findall(r"(\d+)\[PMID\]", _t(a, "SourceOfValidation") or ""),
                    "inheritance": inheritance.get(code, "not annotated"),
                }
            )
    return recs


def fetch_orphanet(spec: DiseaseSpec) -> list[dict[str, Any]]:
    return parse_orphanet(_fetch(ORPHADATA["gene_disease"]), _fetch(ORPHADATA["inheritance"]), spec)


# ---------------------------------------------------------------- mechanism class by rule
DNA_REPAIR = "DNA repair/replication"
MECHANISM_RULES: list[tuple[str, re.Pattern[str]]] = [
    ("DNA interstrand-crosslink repair", re.compile(r"interstrand cross-?link repair", re.I)),
    ("spindle-assembly checkpoint", re.compile(r"spindle assembly checkpoint|spindle checkpoint|mitotic checkpoint", re.I)),
    (DNA_REPAIR, re.compile(r"DNA repair|DNA replication|double-strand break repair|homologous recombination", re.I)),
    ("kinetochore/microtubule attachment", re.compile(r"kinetochore|attachment of spindle microtubules", re.I)),
    ("centrosome/centriole", re.compile(r"centrosom|centriole", re.I)),
    ("sister-chromatid cohesion", re.compile(r"sister chromatid cohesion|cohesin", re.I)),
    (DNA_REPAIR, re.compile(r"DNA damage", re.I)),  # weak generic term, only if nothing more specific matched
]


def classify_terms(terms: list[str]) -> tuple[str, list[str]]:
    """(class, matched terms of the winning rule); first matching rule wins; ('other', []) otherwise."""
    for label, rx in MECHANISM_RULES:
        hit = sorted({t for t in terms if rx.search(t)})
        if hit:
            return label, hit
    return "other", []


def uniprot_annotation(symbol: str) -> dict[str, Any]:
    for field in ("gene_exact", "gene"):
        q = f"({field}:{symbol}) AND (organism_id:9606) AND (reviewed:true)"
        url = f"{UNIPROT}?{urllib.parse.urlencode({'query': q, 'fields': 'accession,gene_primary,keyword,go', 'format': 'tsv'})}"
        rows = _fetch(url).decode("utf-8").splitlines()[1:]
        if rows:
            f = (rows[0].split("\t") + ["", "", "", ""])[:4]
            kws = [x.strip() for x in f[2].split(";") if x.strip()]
            gos = [re.sub(r"\s*\[GO:\d+\]", "", x).strip() for x in f[3].split(";") if x.strip()]
            return {"accession": f[0], "terms": sorted(set(kws + gos))}
    return {"accession": None, "terms": []}


# ---------------------------------------------------------------- Europe PMC evidence
def norm(text: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html.unescape(text))).strip()


def split_sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(\[])", norm(text)) if s.strip()]


def gene_regex(names: list[str]) -> re.Pattern[str]:
    return re.compile(r"(?<![A-Za-z0-9])(" + "|".join(re.escape(n) for n in names) + r")(?![A-Za-z0-9])")


MAX_SENTENCE_CHARS = 350  # longer 'sentences' are usually table or list fragments
CAUSAL_RE = re.compile(r"mutation|variant|biallelic|homozygous|compound heterozygous|cause|identif|associated|deficien|loss|patients?", re.I)


def best_sentence(text: str, gene_rx: re.Pattern[str], assoc_rx: re.Pattern[str]) -> str | None:
    cands = [s for s in split_sentences(text) if len(s) <= MAX_SENTENCE_CHARS and gene_rx.search(s) and assoc_rx.search(s)]
    if not cands:
        return None
    return sorted(cands, key=lambda s: (-len(CAUSAL_RE.findall(s)), len(s), s))[0]


def verify_sentence(sentence: str, source_text: str, where: str = "") -> None:
    """Raise ValueError unless the sentence is a whitespace-normalised substring of the source text."""
    if norm(sentence) not in norm(source_text):
        raise ValueError(f"Sentence not found in source text {where}: {sentence[:80]!r}")


def _epmc_json(path: str, **params: Any) -> dict[str, Any]:
    return json.loads(_fetch(f"{EPMC}/{path}?{urllib.parse.urlencode({**params, 'format': 'json'})}"))


def europepmc_evidence(symbol: str, synonyms: list[str], spec: DiseaseSpec) -> dict[str, Any]:
    names = sorted({symbol, *[s for s in synonyms if len(s) >= 4]})
    gene_rx, assoc_rx = gene_regex(names), spec.assoc_re
    gene_q = " OR ".join(f'"{n}"' for n in names)
    dis_q = " OR ".join(f'"{t}"' for t in spec.terms)
    query = f"({gene_q}) AND ({dis_q})"
    results = _epmc_json("search", query=query, resultType="core", pageSize=25)["resultList"]["result"]
    found: dict[str, Any] | None = None
    checked: list[str] = []
    for r in results:
        pmid, abstract = r.get("pmid"), r.get("abstractText") or ""
        if not pmid or not abstract:
            continue
        s = best_sentence(abstract, gene_rx, assoc_rx)
        if not s:
            continue
        # independent verification against a fresh fetch of this PMID's record
        chk = _epmc_json("search", query=f"EXT_ID:{pmid} AND SRC:MED", resultType="core")["resultList"]["result"]
        verify_sentence(s, (chk[0].get("abstractText") or "") if chk else "", f"PMID {pmid}")
        found = {"pmid": pmid, "sentence": s, "title": norm(r.get("title", "")), "year": r.get("pubYear"), "source_text": "abstract"}
        checked.append(pmid)
        break
    if found is None:  # fallback: open-access PMC full text, at most 3 articles
        for r in [x for x in results if x.get("pmcid") and x.get("isOpenAccess") == "Y"][:3]:
            try:
                xml_text = _fetch(f"{EPMC}/{r['pmcid']}/fullTextXML").decode("utf-8", "replace")
            except RuntimeError:
                continue
            body = norm(xml_text)
            s = best_sentence(body, gene_rx, assoc_rx)
            if s and r.get("pmid"):
                verify_sentence(s, body, r["pmcid"])
                found = {"pmid": r["pmid"], "sentence": s, "title": norm(r.get("title", "")), "year": r.get("pubYear"), "source_text": f"PMC full text {r['pmcid']}"}
                checked.append(r["pmid"])
                break
    if found is None:
        return {"evidence": "unverified", "query": query, "sentence": None, "n_hits": len(results)}
    return {"evidence": "verified", "query": query, "n_hits": len(results), **found}



# ---------------------------------------------------------------- second tier: proposed, not curated
CASE_RE = re.compile(r"patients?|famil(y|ies)|proband|individuals?|siblings?|fetus(es)?|child(ren)?|cases?|cohort|infant|boy|girl", re.I)
NONHUMAN_RE = re.compile(r"arabidopsis|plants?|yeast|drosophila|zebrafish|elegans", re.I)
VARIANT_RE = re.compile(r"biallelic|homozygous|compound heterozygous", re.I)
TITLE_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9])([A-Z][A-Z0-9]+(?:-[A-Z0-9]+)?)(?![A-Za-z0-9])")


def title_symbol_candidates(title: str) -> list[str]:
    """Upper-case tokens of a title that look like human gene symbols (validated afterwards against a symbol source)."""
    return sorted({t for t in TITLE_TOKEN_RE.findall(norm(title)) if len(t) >= 3 and not t.isdigit()})


def uniprot_symbol_ok(token: str) -> bool:
    """True if the token is the approved primary gene symbol of a reviewed human UniProt entry (no free-text guessing)."""
    q = f"(gene_exact:{token}) AND (organism_id:9606) AND (reviewed:true)"
    url = f"{UNIPROT}?{urllib.parse.urlencode({'query': q, 'fields': 'accession,gene_primary', 'format': 'tsv'})}"
    return any(r.split(chr(9))[1:2] == [token] for r in _fetch(url).decode("utf-8").splitlines()[1:])


def europepmc_papers(spec: DiseaseSpec, years: int = 10) -> list[dict[str, Any]]:
    dis_q = " OR ".join(f'TITLE_ABS:"{t}"' for t in spec.terms)
    first = date.today().year - years
    query = f'({dis_q}) AND (biallelic OR homozygous OR "compound heterozygous") AND mosaic AND PUB_YEAR:[{first} TO {date.today().year}]'
    return _epmc_json("search", query=query, resultType="core", pageSize=100)["resultList"]["result"]


def discover_proposed(
    spec: DiseaseSpec,
    known_symbols: set[str],
    curated: set[str],
    papers: list[dict[str, Any]],
    symbol_ok: Callable[[str], bool] = uniprot_symbol_ok,
) -> dict[str, dict[str, Any]]:
    """Genes named in the title of recent case papers but absent from the curated list.

    A title token counts as a gene only if it is a symbol known from the HPO gene file or the approved primary symbol of a
    reviewed human UniProt entry. The paper is kept only if its abstract describes cases (CASE_RE) and contains a sentence
    naming the gene and a disease term (verified as a substring of the abstract). Genes with no such sentence are dropped.
    """
    out: dict[str, dict[str, Any]] = {}
    seen_ok: dict[str, bool] = {}
    for r in sorted(papers, key=lambda r: (str(r.get("pmid")), r.get("title", ""))):
        pmid, abstract = r.get("pmid"), r.get("abstractText") or ""
        if not pmid or not abstract or not CASE_RE.search(abstract) or NONHUMAN_RE.search(abstract):
            continue
        cands = title_symbol_candidates(r.get("title", ""))
        if curated.intersection(cands) or not spec.assoc_re.search(abstract):
            continue  # a curated gene already explains the paper (e.g. a second, unrelated gene in a double-homozygosity title)
        for tok in cands:
            if tok not in known_symbols:
                if tok not in seen_ok:
                    seen_ok[tok] = symbol_ok(tok)
                if not seen_ok[tok]:
                    continue
            rx = gene_regex([tok])
            s = best_sentence(abstract, rx, spec.assoc_re)
            basis = "gene and disease term in the same sentence"
            if not s:  # disease named elsewhere in the abstract (the search required it): accept gene + variant-type sentence
                s = best_sentence(abstract, rx, VARIANT_RE)
                basis = "gene and variant type in the sentence; disease term elsewhere in the same abstract"
            if not s:
                continue
            verify_sentence(s, abstract, f"PMID {pmid}")
            e = out.setdefault(tok, {"papers": []})
            e["papers"].append({"pmid": pmid, "year": r.get("pubYear"), "title": norm(r.get("title", "")), "sentence": s, "sentence_basis": basis})
    return out


# ---------------------------------------------------------------- Track 1 cross-reference (symbols only)
def track1_symbols(path: Path = SUBMISSION_GENES) -> list[str]:
    return sorted({g["gene"] for g in json.loads(path.read_text(encoding="utf-8"))["genes"]})


# ---------------------------------------------------------------- assembly
def build(
    spec: DiseaseSpec,
    hpoa: Path | None = None,
    g2p: Path | None = None,
    obo: Path | None = None,
    orphanet_fn: Callable[[DiseaseSpec], list[dict[str, Any]]] = fetch_orphanet,
    annotate_fn: Callable[[str], dict[str, Any]] = uniprot_annotation,
    evidence_fn: Callable[[str, list[str], DiseaseSpec], dict[str, Any]] = europepmc_evidence,
    track1_fn: Callable[[], list[str]] = track1_symbols,
    papers_fn: Callable[[DiseaseSpec], list[dict[str, Any]]] = europepmc_papers,
    symbol_ok: Callable[[str], bool] = uniprot_symbol_ok,
) -> dict[str, Any]:
    hpoa = hpoa or find_source("phenotype.hpoa")
    g2p = g2p or find_source("genes_to_phenotype.txt")
    obo = obo or find_source("hp.obo")
    hpo_dis = hpo_diseases(spec, hpoa, obo)
    gmap = hpo_gene_map(g2p, set(hpo_dis))
    orphanet_error = None
    try:
        orph = orphanet_fn(spec)
    except Exception as error:  # noqa: BLE001 - Orphadata unreachable: continue with HPO only
        orph, orphanet_error = [], repr(error)

    genes: dict[str, dict[str, Any]] = {}

    def entry(sym: str) -> dict[str, Any]:
        return genes.setdefault(sym, {"symbol": sym, "diseases": [], "sources": set(), "synonyms": set(), "hgnc_id": None})

    for did, d in sorted(hpo_dis.items()):
        for sym in sorted(gmap.get(did, ())):
            e = entry(sym)
            e["sources"].add("HPO (genes_to_phenotype.txt + phenotype.hpoa)")
            e["diseases"].append({"source": "HPO", "id": did, "name": d["name"], "subtype": subtype_of(d["name"]), "inheritance": "; ".join(d["inheritance"]) or "not annotated"})
    for r in sorted(orph, key=lambda r: (r["symbol"], r["orpha"])):
        sym = r["symbol"]
        if sym not in genes:  # merge an outdated Orphanet symbol into the HPO symbol when it is a listed synonym
            alias = [s for s in r["synonyms"] if s in genes]
            sym = alias[0] if alias else sym
        e = entry(sym)
        e["sources"].add("Orphanet (Orphadata en_product6 + en_product9_ages)")
        e["synonyms"].update(r["synonyms"] + ([r["symbol"]] if r["symbol"] != sym else []))
        e["hgnc_id"] = r["hgnc_id"]
        e["diseases"].append(
            {
                "source": "Orphanet", "id": r["orpha"], "name": r["disease_name"], "subtype": subtype_of(r["disease_name"]), "inheritance": r["inheritance"],
                "association_type": r["association_type"], "association_status": r["association_status"], "orphanet_pmids": r["orphanet_pmids"],
            }
        )
    known = {line.split(chr(9))[1] for line in g2p.read_text(encoding="utf-8").splitlines()[1:] if chr(9) in line}
    proposed = discover_proposed(spec, known, set(genes), papers_fn(spec), symbol_ok)
    for sym in sorted(proposed):
        e = entry(sym)
        e["tier"] = "proposed, not curated"
        e["sources"].add("Europe PMC literature (title/abstract search)")
        e["diseases"].append({"source": "Europe PMC", "id": None, "name": f"{spec.label} (proposed association)", "subtype": None, "inheritance": "not annotated"})
    t1 = set(track1_fn()) if spec.track1_crossref else set()
    out_genes = []
    for sym in sorted(genes):
        e = genes[sym]
        ann = annotate_fn(sym)
        cls, matched = classify_terms(ann["terms"])
        subtypes = sorted({d["subtype"] for d in e["diseases"] if d["subtype"]})
        g = {
            "symbol": sym,
            "hgnc_id": e["hgnc_id"],
            "diseases": sorted(e["diseases"], key=lambda d: (d["source"], d["id"] or "")),
            "disease_subtype": ", ".join(subtypes) if subtypes else "not assigned in sources",
            "inheritance_as_annotated": sorted({d["inheritance"] for d in e["diseases"]}),
            "sources": sorted(e["sources"]),
            "n_sources": len(e["sources"]),
            "mechanism_class": cls,
            "mechanism_basis": {"uniprot_accession": ann["accession"], "matched_terms": matched},
            "tier": e.get("tier", "curated"),
            "association_evidence": evidence_fn(sym, sorted(e["synonyms"]), spec),
        }
        if sym in proposed:  # the proposing paper's own verified sentence is the evidence; counted independent papers
            p0 = sorted(proposed[sym]["papers"], key=lambda x: (-len(CAUSAL_RE.findall(x["sentence"])), len(x["sentence"]), x["pmid"]))[0]
            g["association_evidence"] = {"evidence": "verified", "pmid": p0["pmid"], "sentence": p0["sentence"], "year": p0["year"], "title": p0["title"], "source_text": "abstract", "sentence_basis": p0["sentence_basis"]}
            g["n_independent_papers"] = len({x["pmid"] for x in proposed[sym]["papers"]})
            g["proposing_pmids"] = sorted({x["pmid"] for x in proposed[sym]["papers"]})
        if spec.track1_crossref:
            g["in_track1_output"] = sym in t1
        out_genes.append(g)
    by_class: dict[str, int] = defaultdict(int)
    for g in out_genes:
        by_class[g["mechanism_class"]] += 1
    result: dict[str, Any] = {
        "header": {
            "statement": HEADER_STATEMENT,
            "access_date": date.today().isoformat(),
            "disease": spec.label,
            "search_terms": list(spec.terms),
            "identifiers": list(spec.ids),
            "sources": {
                "HPO": {**hpo_release(hpoa), "files": [hpoa.name, g2p.name, obo.name]},
                "Orphanet": {"urls": ORPHADATA, "error": orphanet_error},
                "UniProt": UNIPROT,
                "Europe PMC": EPMC,
            },
            "hpo_diseases_matched": {k: v["name"] for k, v in sorted(hpo_dis.items())},
            "orphanet_disorders_matched": sorted({f"{r['orpha']} {r['disease_name']}" for r in orph}),
            "mechanism_rule": [[label, rx.pattern] for label, rx in MECHANISM_RULES] + [["other", "no rule matched"]],
        },
        "n_genes": len(out_genes),
        "n_genes_by_mechanism_class": dict(sorted(by_class.items())),
        "n_genes_by_tier": {t: sum(g["tier"] == t for g in out_genes) for t in ("curated", "proposed, not curated")},
        "n_genes_with_verified_sentence": sum(g["association_evidence"]["evidence"] == "verified" for g in out_genes),
        "genes": out_genes,
    }
    if spec.track1_crossref:
        result["track1"] = {
            "note": "gene symbols only, from out/submission_genes.json",
            "genes_in_track1_output": [g["symbol"] for g in out_genes if g["in_track1_output"]],
            "genes_not_in_track1_output": [g["symbol"] for g in out_genes if not g["in_track1_output"]],
        }
    return result


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("disease", help="mva | fanconi_anaemia | custom")
    ap.add_argument("--terms", nargs="*", default=[])
    ap.add_argument("--ids", nargs="*", default=[])
    ap.add_argument("--abbreviations", nargs="*", default=[])
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    spec = SPECS.get(a.disease) or DiseaseSpec(a.disease, " / ".join(a.terms or a.ids), tuple(a.terms), tuple(a.ids), tuple(a.abbreviations), ROOT / (a.out or f"out/landscape_{a.disease}.json"))
    result = build(spec)
    spec.out.write_text(json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    for g in result["genes"]:
        print(g["symbol"], g["tier"], g.get("n_independent_papers", ""), g["disease_subtype"], g["mechanism_class"], g["n_sources"], g["association_evidence"]["evidence"], g.get("in_track1_output", ""))


if __name__ == "__main__":
    main()
