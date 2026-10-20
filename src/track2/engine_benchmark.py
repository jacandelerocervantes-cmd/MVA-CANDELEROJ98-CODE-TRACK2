"""Validation benchmark for the literature-grounded evidence engine (built BEFORE any scorer exists).

Run: PYTHONPATH=src python -m track2.engine_benchmark      (writes out/engine/benchmark.json)

Pre-registration: protocol/addendum_engine.md. Nothing here scores a drug.

Disease pool and selection (deterministic, seed 20261015)
  * pool = the curated rare-disease candidates of track2.empirical_literature_control.CANDIDATES (reused) plus every
    PrimeKG disease with >= 2 approved universe drugs linked by an `indication` edge;
  * a candidate is matched to HPO (phenotype.hpoa + genes_to_phenotype.txt) and Orphanet (Orphadata en_product6/9)
    disorders whose name contains the candidate name as whole words (same term_regex as disease_landscape); more than
    MAX_MATCHED matched disorders = umbrella term, excluded;
  * rare + genetic = a disorder whose normalised name EQUALS the candidate name, is catalogued by Orphanet (Orphadata, or
    an ORPHA id in HPO) and has a Mendelian/mitochondrial/sex-linked mode of inheritance (HPO aspect I or Orphadata); the
    earlier-proxy curated candidates (rare genetic by construction) may instead pass with any matched Orphanet-catalogued
    disorder of Mendelian inheritance (e.g. umbrella names such as "mucopolysaccharidosis");
  * (a) >= 2 gene symbols (union of HPO and Orphanet) ; (b) gold 1 >= 2 universe drugs ; (c) MVA and Fanconi anaemia
    excluded (application cases) ; (c) a disease is excluded when ANY gold 1 drug is an enzyme (replacement), gene or cell
    therapy (universe molecule_type Enzyme/Gene/Cell): conservative reading of "mechanism is the drug target itself";
  * diseases with an identical gold 1 set are not independent: only the alphabetically first is kept;
  * N_TARGET diseases are drawn with random.Random(SEED) from the eligible ones (all of them if fewer qualify).

Gold 1 (primary): universe drugs with a PrimeKG `indication` edge to the disease (as in the earlier proxy).
Gold 2 (temporal): universe drugs named (whole word, case-insensitive; names and other names of the intervention) as an
  intervention of an interventional ClinicalTrials.gov study whose condition search matches the disease name, with phase
  2 or later and START DATE after 2018-12-31; evidence cutoff year Y = 2018 for the temporal analysis.
Mechanism class per gene: disease_landscape.classify_terms on the reviewed UniProt entry (same fixed rules).

All network responses are cached under data/cache_engine/ (git-ignored by data/*).
"""
from __future__ import annotations

import hashlib
import json
import random
import re
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261015
N_TARGET = 30
MIN_GENES = 2
MIN_GOLD = 2
MAX_MATCHED = 30
CUTOFF_YEAR = 2018
GOLD2_AFTER = date(2018, 12, 31)  # start date must be strictly after this day
CACHE = ROOT / "data" / "cache_engine"
CTGOV = "https://clinicaltrials.gov/api/v2/studies"
ORPHADATA = {
    "gene_disease": "https://www.orphadata.com/data/xml/en_product6.xml",
    "inheritance": "https://www.orphadata.com/data/xml/en_product9_ages.xml",
}
APPLICATION_CASES = re.compile(r"mosaic[\s-]+variegated|premature[\s-]+chromatid|fanconi", re.I)
BIOLOGIC_THERAPY_TYPES = frozenset({"Enzyme", "Gene", "Cell"})  # ERT, gene therapy, cell therapy
MENDELIAN = re.compile(r"autosomal|x-linked|y-linked|mitochondrial|semidominant|digenic|oligogenic", re.I)
SYMBOL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.@-]+$")  # drops placeholders such as "-"
NON_DRUG_TYPES = frozenset({"DEVICE", "PROCEDURE", "BEHAVIORAL", "RADIATION", "DIAGNOSTIC_TEST"})
PHASE2_PLUS = frozenset({"PHASE2", "PHASE3", "PHASE4"})  # PHASE2/PHASE3 studies list both


# ----------------------------------------------------------------------------- generic helpers
def term_regex(name: str) -> "re.Pattern[str]":
    """Whole-word, case-insensitive, space/hyphen tolerant regex of a disease name (as disease_landscape.term_regex)."""
    parts = [re.escape(w) for w in re.split(r"[\s-]+", name.strip()) if w]
    body = r"[\s-]+".join(re.sub(r"an(a)?emia", "an(a)?emia", p, flags=re.I) if re.search("anaemia|anemia", p, re.I) else p
                          for p in parts)
    return re.compile(r"(?<![\w])" + body + r"(?![\w])", re.I)


def norm_name(name: str) -> str:
    """Comparison form of a disease name: lower case, hyphens/commas/brackets removed, anaemia = anemia."""
    n = re.sub(r"\([^)]*\)", " ", name.lower()).replace("anaemia", "anemia")
    return re.sub(r"\s+", " ", re.sub(r"[-,]", " ", n)).strip()


def atc_class(atc_codes: Sequence[str], level: int = 4) -> Optional[str]:
    """Pharmacological class = first `level` characters (ATC level 3 by default) of the lexicographically first code."""
    codes = sorted(c for c in atc_codes if c and len(c) >= level)
    return codes[0][:level] if codes else None


def http_get_json(url: str, params: dict, cache_dir: Path = CACHE / "http",
                  getter: Optional[Callable[[str, dict], Any]] = None) -> Any:
    """Cached GET returning JSON; retry with exponential backoff on SSL/transient errors; <= 4 requests/s."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256((url + "?" + json.dumps(params, sort_keys=True)).encode()).hexdigest()[:24]
    f = cache_dir / f"{key}.json"
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    if getter is None:
        import requests

        def getter(u, p):  # type: ignore[no-redef]
            for attempt in range(6):
                try:
                    time.sleep(0.25)
                    r = requests.get(u, params=p, timeout=60, headers={"User-Agent": "hselith-track2/0.3 (research)"})
                    r.raise_for_status()
                    return r.json()
                except Exception:  # noqa: BLE001 - SSL / transient errors
                    if attempt == 5:
                        raise
                    time.sleep(2 ** attempt)
    data = getter(url, params)
    f.write_text(json.dumps(data), encoding="utf-8")
    return data


# ----------------------------------------------------------------------------- universe and PrimeKG gold 1
def load_universe(path: Path) -> List[dict]:
    return json.loads(path.read_text(encoding="utf-8"))["drugs"]


def universe_index(drugs: Sequence[dict]) -> Dict[str, dict]:
    """lower-case name -> record (first record when a name repeats; deterministic by drug_id)."""
    out: Dict[str, dict] = {}
    for d in sorted(drugs, key=lambda d: d["drug_id"]):
        out.setdefault(d["name"].lower(), d)
    return out


def primekg_gold1(root: Path, drug_names: Set[str]) -> Dict[str, Set[str]]:
    """PrimeKG disease name -> lower-case universe drug names (cached derivation of the earlier proxy's function)."""
    f = CACHE / "primekg_indications.json"
    if f.exists():
        raw = json.loads(f.read_text(encoding="utf-8"))
        return {k: set(v) & drug_names for k, v in raw.items() if set(v) & drug_names}
    from track2 import empirical_literature_control as elc
    av = elc.load_primekg_indications(root / "data" / "sources" / "primekg", drug_names)
    CACHE.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps({k: sorted(v) for k, v in av.items()}, sort_keys=True), encoding="utf-8")
    return av


# ----------------------------------------------------------------------------- HPO / Orphanet gene sources
def load_hpo(hpoa: Path, g2p: Path, obo: Path) -> Dict[str, dict]:
    """HPO disease id -> {name, inheritance(list of labels), genes(set)}."""
    dis: Dict[str, dict] = {}
    inh_ids: Dict[str, Set[str]] = defaultdict(set)
    for line in hpoa.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or line.startswith("database_id") or not line.strip():
            continue
        f = line.split("\t")
        dis.setdefault(f[0], {"name": f[1], "inheritance": [], "genes": set()})
        if len(f) > 10 and f[10] == "I":
            inh_ids[f[0]].add(f[3])
    wanted = {i for s in inh_ids.values() for i in s}
    labels, cur = {}, None
    for line in obo.read_text(encoding="utf-8").splitlines():
        if line.startswith("id: "):
            cur = line[4:].strip()
        elif line.startswith("name: ") and cur in wanted:
            labels[cur] = line[6:].strip()
    for did, ids in inh_ids.items():
        dis[did]["inheritance"] = sorted(labels.get(i, i) for i in ids)
    for line in g2p.read_text(encoding="utf-8").splitlines()[1:]:
        f = line.split("\t")
        if len(f) > 5 and f[5] in dis:
            dis[f[5]]["genes"].add(f[1])
    return dis


def parse_orphanet_all(gene_xml: bytes, inh_xml: Optional[bytes]) -> Dict[str, dict]:
    """ORPHA code -> {name, inheritance(list), genes(set of symbols), hgnc{symbol: id}} for every Orphadata disorder."""
    inheritance: Dict[str, List[str]] = {}
    if inh_xml:
        for d in ET.fromstring(inh_xml).iter("Disorder"):
            code = (d.findtext("OrphaCode") or "").strip()
            inheritance[code] = sorted(x.text.strip() for x in d.findall("./TypeOfInheritanceList/TypeOfInheritance/Name") if x.text)
    out: Dict[str, dict] = {}
    for d in ET.fromstring(gene_xml).iter("Disorder"):
        code = (d.findtext("OrphaCode") or "").strip()
        rec = out.setdefault(f"ORPHA:{code}", {"name": (d.findtext("Name") or "").strip(),
                                               "inheritance": inheritance.get(code, []), "genes": set(), "hgnc": {}})
        for a in d.findall("./DisorderGeneAssociationList/DisorderGeneAssociation"):
            g = a.find("Gene")
            sym = (g.findtext("Symbol") or "").strip() if g is not None else ""
            if sym:
                rec["genes"].add(sym)
                for r in g.findall("./ExternalReferenceList/ExternalReference"):
                    if r.findtext("Source") == "HGNC":
                        rec["hgnc"][sym] = f"HGNC:{r.findtext('Reference')}"
    return out


def fetch_orphanet(fetch: Optional[Callable[[str], bytes]] = None) -> Dict[str, dict]:
    fetch = fetch or _cached_bytes
    return parse_orphanet_all(fetch(ORPHADATA["gene_disease"]), fetch(ORPHADATA["inheritance"]))


def _cached_bytes(url: str) -> bytes:
    from track2 import disease_landscape as dl
    old = dl.CACHE
    dl.CACHE = CACHE / "landscape"  # disease_landscape._fetch caches (and retries) under the module-level CACHE
    try:
        return dl._fetch(url)
    finally:
        dl.CACHE = old


def match_disorders(name: str, hpo: Dict[str, dict], orph: Dict[str, dict]) -> Dict[str, Any]:
    """Disorders of HPO and Orphanet whose name contains `name` as whole words; union of genes; rarity/genetic flags."""
    rx = term_regex(name)
    h = {i: d for i, d in hpo.items() if rx.search(d["name"])}
    o = {i: d for i, d in orph.items() if rx.search(d["name"])}
    genes: Set[str] = set()
    for d in list(h.values()) + list(o.values()):
        genes |= {g for g in d["genes"] if SYMBOL_RE.match(g)}
    inh = sorted({x for d in list(h.values()) + list(o.values()) for x in d["inheritance"]})
    target = norm_name(name)
    exact = [d for i, d in list(o.items()) + list(h.items())
             if norm_name(d["name"]) == target and (i in o or i.startswith("ORPHA:"))
             and any(MENDELIAN.search(x) for x in d["inheritance"])]
    return {"exact_orphanet_mendelian": bool(exact), "hpo": sorted(h), "orphanet": sorted(o), "n_matched": len(set(h) | set(o)), "genes": sorted(genes),
            "inheritance": inh,
            "orphanet_catalogued": bool(o) or any(i.startswith("ORPHA:") for i in h),
            "mendelian": any(MENDELIAN.search(x) for x in inh)}


# ----------------------------------------------------------------------------- selection
def evaluate_candidate(name: str, gold_names: Set[str], match: Optional[Dict[str, Any]],
                       uidx: Dict[str, dict], curated: bool = False) -> Tuple[bool, str]:
    """(eligible, reason). Checks in a fixed order; the first failing check is the logged reason."""
    if APPLICATION_CASES.search(name):
        return False, "application case (MVA / Fanconi anaemia) excluded by design"
    if match is None or (not match["hpo"] and not match["orphanet"]):
        return False, "no HPO or Orphanet disorder with this name"
    if match["n_matched"] > MAX_MATCHED:
        return False, f"umbrella term: {match['n_matched']} matched disorders (> {MAX_MATCHED})"
    if not match["exact_orphanet_mendelian"] and not (curated and match["orphanet_catalogued"] and match["mendelian"]):
        return False, ("not a rare genetic disease: no Orphanet-catalogued disorder of exactly this name with Mendelian/"
                       "mitochondrial/sex-linked inheritance (and not a curated rare-disease candidate)")
    if len(match["genes"]) < MIN_GENES:
        return False, f"fewer than {MIN_GENES} genes ({len(match['genes'])})"
    gold = {g for g in gold_names if g in uidx}
    if len(gold) < MIN_GOLD:
        return False, f"fewer than {MIN_GOLD} approved universe drugs indicated ({len(gold)})"
    bad = sorted(g for g in gold if uidx[g].get("molecule_type") in BIOLOGIC_THERAPY_TYPES)
    if bad:
        return False, "gold drug is an enzyme replacement / gene / cell therapy: " + ", ".join(bad)
    return True, "eligible"


def drop_duplicate_gold(eligible: Sequence[str], gold: Dict[str, Set[str]]) -> Tuple[List[str], Dict[str, str]]:
    """Diseases with an identical gold 1 set are one PrimeKG indication cluster (e.g. Parkinson variants): keep the
    alphabetically first, return the dropped ones with the name they duplicate."""
    keep, seen, dropped = [], {}, {}
    for n in sorted(eligible):
        key = frozenset(gold[n])
        if key in seen:
            dropped[n] = seen[key]
        else:
            seen[key] = n
            keep.append(n)
    return keep, dropped


def draw(eligible: Sequence[str], seed: int = SEED, n: int = N_TARGET) -> List[str]:
    """Deterministic draw: names are sorted first, so the result does not depend on discovery order."""
    pool = sorted(eligible)
    return sorted(random.Random(seed).sample(pool, min(n, len(pool))))


# ----------------------------------------------------------------------------- ClinicalTrials.gov gold 2
def _parse_date(s: Optional[str]) -> Optional[date]:
    """'2021', '2021-03' or '2021-03-15' -> date (missing parts = first day)."""
    if not s:
        return None
    p = (s.split("-") + ["01", "01"])[:3]
    try:
        return date(int(p[0]), int(p[1]), int(p[2]))
    except (ValueError, IndexError):
        return None


def parse_study(study: dict) -> dict:
    ps = study.get("protocolSection", {})
    ivs = (ps.get("armsInterventionsModule") or {}).get("interventions") or []
    return {
        "nct": (ps.get("identificationModule") or {}).get("nctId"),
        "study_type": (ps.get("designModule") or {}).get("studyType"),
        "phases": list((ps.get("designModule") or {}).get("phases") or []),
        "start": _parse_date(((ps.get("statusModule") or {}).get("startDateStruct") or {}).get("date")),
        "interventions": [{"name": i.get("name", ""), "type": i.get("type", ""), "other": list(i.get("otherNames") or [])}
                          for i in ivs],
    }


def study_qualifies(p: dict) -> bool:
    return (p["study_type"] == "INTERVENTIONAL" and bool(PHASE2_PLUS & set(p["phases"]))
            and p["start"] is not None and p["start"] > GOLD2_AFTER)


def drug_name_regex(names: Sequence[str]) -> "re.Pattern[str]":
    clean = sorted({n.lower() for n in names if n and len(n) >= 4}, key=len, reverse=True)
    return re.compile("|".join(r"(?<![\w])" + re.escape(t) + r"(?![\w])" for t in clean), re.I)


def gold2_from_studies(studies: Sequence[dict], dre: "re.Pattern[str]") -> Dict[str, List[dict]]:
    """lower-case universe drug name -> sorted list of {nct, start_year, phases}; only qualifying studies count."""
    out: Dict[str, Dict[str, dict]] = defaultdict(dict)
    for s in studies:
        p = parse_study(s)
        if not p["nct"] or not study_qualifies(p):
            continue
        for iv in p["interventions"]:
            if iv["type"] in NON_DRUG_TYPES:
                continue
            for text in [iv["name"]] + iv["other"]:
                for m in dre.finditer(text or ""):
                    out[m.group(0).lower()][p["nct"]] = {"nct": p["nct"], "start_year": p["start"].year,
                                                         "phases": sorted(p["phases"])}
    return {d: sorted(v.values(), key=lambda t: t["nct"]) for d, v in sorted(out.items())}


def fetch_studies(condition: str, getter: Optional[Callable[[str, dict], Any]] = None) -> List[dict]:
    """All interventional phase 2+ studies starting after the cutoff for a condition search (paged, cached)."""
    studies: List[dict] = []
    token = None
    while True:
        params = {
            "query.cond": condition, "pageSize": 1000, "format": "json",
            "filter.advanced": "AREA[StudyType]INTERVENTIONAL AND AREA[Phase](PHASE2 OR PHASE3 OR PHASE4) AND "
                               f"AREA[StartDate]RANGE[{(GOLD2_AFTER.replace(day=31)).isoformat()},MAX]",
            "fields": "NCTId|StudyType|Phase|StartDate|InterventionName|InterventionOtherName|InterventionType",
        }
        if token:
            params["pageToken"] = token
        data = http_get_json(CTGOV, params, CACHE / "ctgov", getter)
        studies.extend(data.get("studies") or [])
        token = data.get("nextPageToken")
        if not token:
            return studies


# ----------------------------------------------------------------------------- mechanism class of genes
def gene_mechanism(symbol: str, annotate: Optional[Callable[[str], dict]] = None) -> dict:
    from track2 import disease_landscape as dl
    if annotate is None:
        def annotate(sym):  # type: ignore[no-redef]
            old = dl.CACHE
            dl.CACHE = CACHE / "landscape"
            try:
                return dl.uniprot_annotation(sym)
            finally:
                dl.CACHE = old
    try:
        ann = annotate(symbol)
    except Exception as e:  # noqa: BLE001 - keep the gene, class 'other', failure recorded
        ann = {"terms": [], "accession": None, "error": repr(e)[:200]}
    cls, matched = dl.classify_terms(ann["terms"])
    return {"symbol": symbol, "mechanism_class": cls, "matched_terms": matched, "uniprot_accession": ann["accession"]}


# ----------------------------------------------------------------------------- driver
def build(root: Path = ROOT, hpo: Optional[Dict[str, dict]] = None, orph: Optional[Dict[str, dict]] = None,
          studies_fn: Callable[[str], List[dict]] = fetch_studies,
          annotate: Optional[Callable[[str], dict]] = None) -> Dict[str, Any]:
    from track2 import empirical_literature_control as elc
    drugs = load_universe(root / "data" / "universe.json")
    uidx = universe_index(drugs)
    gold1 = primekg_gold1(root, set(uidx))
    if hpo is None:
        src = root / "data" / "sources"
        hpo = load_hpo(next(src.rglob("phenotype.hpoa")), next(src.rglob("genes_to_phenotype.txt")), next(src.rglob("hp.obo")))
    orph_error = None
    if orph is None:
        try:
            orph = fetch_orphanet()
        except Exception as e:  # noqa: BLE001 - continue with HPO only, recorded
            orph, orph_error = {}, repr(e)
    curated = {c[0] for c in elc.CANDIDATES}
    pool = sorted({d for d, g in gold1.items() if len(g) >= MIN_GOLD} | (curated & set(gold1)))
    log, eligible, matches = [], [], {}
    for name in pool:
        m = match_disorders(name, hpo, orph)
        ok, why = evaluate_candidate(name, gold1[name], m, uidx, name in curated)
        matches[name] = m
        if ok:
            eligible.append(name)
        log.append({"disease": name, "eligible": ok, "reason": why, "curated_candidate_of_earlier_proxy": name in curated})
    eligible, dup = drop_duplicate_gold(eligible, gold1)
    chosen = draw(eligible)
    for e in log:
        if e["disease"] in dup:
            e["eligible"], e["reason"] = False, f"duplicate gold 1 set of {dup[e['disease']]} (not independent)"
        elif e["eligible"]:
            e["reason"] = "selected (seeded draw)" if e["disease"] in chosen else "eligible, not drawn"
    dre = drug_name_regex(sorted(uidx))
    diseases = []
    for name in chosen:
        m = matches[name]
        g1 = sorted(gold1[name])
        studies = studies_fn(name)
        g2 = gold2_from_studies(studies, dre)
        g2 = {d: t for d, t in g2.items() if d in uidx}
        diseases.append({
            "name": name,
            "identifiers": {"hpo": m["hpo"], "orphanet": m["orphanet"], "primekg_name": name},
            "inheritance": m["inheritance"],
            "genes": [gene_mechanism(s, annotate) for s in m["genes"]],
            "gold1": [{"name": d, "drug_id": uidx[d]["drug_id"], "molecule_type": uidx[d].get("molecule_type"),
                       "atc": uidx[d].get("atc", []), "atc_class": atc_class(uidx[d].get("atc", []))} for d in g1],
            "gold2": [{"name": d, "drug_id": uidx[d]["drug_id"], "in_gold1": d in gold1[name], "trials": t}
                      for d, t in g2.items()],
            "gold2_query": {"endpoint": CTGOV, "query.cond": name, "n_studies_returned": len(studies)},
        })
    reasons: Dict[str, int] = defaultdict(int)
    for e in log:
        reasons[re.sub(r"[:(].*$", "", e["reason"]).strip()] += 1
    return {
        "seed": SEED, "n_target": N_TARGET, "n_benchmark": len(diseases), "cutoff_year_Y": CUTOFF_YEAR,
        "gold2_rule": "universe drug is an intervention of an interventional ClinicalTrials.gov study with phase 2 or later "
                      "and start date after 2018-12-31",
        "gold1_rule": "universe drug with a PrimeKG `indication` edge to the disease",
        "n_universe_drugs": len(uidx), "n_candidates_evaluated": len(pool), "n_eligible": len(eligible),
        "orphanet_error": orph_error,
        "gold1_sizes": {d["name"]: len(d["gold1"]) for d in diseases},
        "gold2_sizes": {d["name"]: len(d["gold2"]) for d in diseases},
        "diseases": diseases,
        "selection_log": sorted(log, key=lambda e: e["disease"]),
        "selection_reason_counts": dict(sorted(reasons.items())),
        "limitations": [
            "Gold 1 (PrimeKG/DrugBank/DrugCentral indications) is incomplete and partly circular for drugs developed against the disease gene.",
            "Gold 2 matches universe drug names in intervention names (no synonym table); a drug in a trial may be a comparator or background therapy.",
            "ClinicalTrials.gov condition search is by disease name only; trials registered under synonyms are missed.",
            "Rare genetic = exact-name Orphanet-catalogued disorder with Mendelian/mitochondrial/sex-linked inheritance (curated candidates of the earlier proxy: looser match); gene lists are the union over all whole-word-matched subtypes.",
            "Enzyme, gene and cell therapies in gold 1 exclude the whole disease (conservative reading).",
        ],
    }


def main() -> None:
    out = build()
    p = ROOT / "out" / "engine" / "benchmark.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: out[k] for k in ("n_benchmark", "n_eligible", "gold1_sizes", "gold2_sizes", "selection_reason_counts")},
                     indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
