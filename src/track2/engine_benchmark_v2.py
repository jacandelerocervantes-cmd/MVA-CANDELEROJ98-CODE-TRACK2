"""Benchmark v2 of the literature-grounded evidence engine: a documented extension of benchmark v1.

Run: PYTHONPATH=src python -m track2.engine_benchmark_v2      (writes out/engine/benchmark_v2.json)

Pre-registration: protocol/addendum_engine_v2.md (written before this module). v1 (engine_benchmark.py, benchmark.json,
addendum_engine.md) is imported/read only and never modified; the v1-only result stays the primary result. Nothing here
scores a drug. All thresholds and endpoints of v1 are unchanged.

What v2 changes (the four limitations of v1)
  1. Disease matching by identifiers and synonyms: PrimeKG disease node (MONDO ids) -> OMIM / OMIMPS / Orphanet via the
     `xref:` lines of mondo.obo (only `MONDO:equivalentTo` xrefs; the SSSOM sets are a re-serialisation of them) ->
     HPO and Orphadata disorders by id, plus whole-word name containment of the PrimeKG name, MONDO name, MONDO EXACT
     synonyms (no abbreviation type), Orphanet synonyms (en_product1) and the names of the id-matched disorders.
     "Rare genetic" = the identifier-matched disorders of the disease (one entity by MONDO equivalence) are
     Orphanet-catalogued and Mendelian, OR (as v1) a disorder of exactly a name of the disease is catalogued and Mendelian.
  2. Gold 3 (independent of PrimeKG): universe drugs whose current FDA label (openFDA drug/label.json,
     indications_and_usage) names the disease or one of its synonyms as a whole phrase.
  3. Gold 2 rebuilt: a drug counts only if it is an intervention of an EXPERIMENTAL arm group of a qualifying study.
  4. `mechanism_class_general` per gene (track2.mechanism_general), next to the v1 `mechanism_class`.
Every other inclusion/exclusion rule and the seed of v1 are kept (see evaluate_candidate_v2 and build).
"""
from __future__ import annotations

import json
import random
import re
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from track2 import engine_benchmark as v1
from track2.mechanism_general import CLASSES, classify_general

ROOT = v1.ROOT
SEED = v1.SEED  # 20261015
CAP_TOTAL = 60  # v1 diseases + new diseases; new diseases are drawn only when more than CAP_TOTAL - n_v1 qualify
TARGET_TOTAL = 40
MIN_GENES, MIN_GOLD, MAX_MATCHED = v1.MIN_GENES, v1.MIN_GOLD, v1.MAX_MATCHED
CACHE = v1.CACHE
MIN_NAME_CHARS = 5
LABELS_PER_DRUG = 5
OPENFDA = "https://api.fda.gov/drug/label.json"
ORPHADATA_NAMES = "https://www.orphadata.com/data/xml/en_product1.xml"
MONDO_ROOT = ROOT / "data" / "sources" / "mondo" / "mondo.obo"
EXPERIMENTAL = "EXPERIMENTAL"
GOLD2_SKIP_TYPES = v1.NON_DRUG_TYPES
XREF_RE = re.compile(r"^xref: (OMIM|OMIMPS|Orphanet):(\d+)(?: \{(.*)\})?\s*$")
SYN_RE = re.compile(r'^synonym: "((?:[^"\\]|\\.)*)" (EXACT|RELATED|BROAD|NARROW)(?: ([A-Za-z_:0-9]+))? \[')


# ----------------------------------------------------------------------------- names
def clean_name(name: str) -> str:
    """Display form of a disease name: parenthetical qualifiers such as '(disease)' removed, spaces collapsed."""
    return re.sub(r"\s+", " ", re.sub(r"\s*\([^)]*\)", "", name)).strip()


def build_name_set(names: Iterable[str]) -> List[str]:
    """Unique (by v1 normalised form) names of at least MIN_NAME_CHARS characters; obsolete names dropped; sorted."""
    seen: Dict[str, str] = {}
    for n in sorted({clean_name(x) for x in names if x}):
        k = v1.norm_name(n)
        if len(k) >= MIN_NAME_CHARS and not k.startswith("obsolete") and k not in seen:
            seen[k] = n
    return sorted(seen.values())


def names_regex(names: Sequence[str]) -> Optional["re.Pattern[str]"]:
    """Whole-phrase, case-insensitive, space/hyphen tolerant alternation of the names (v1 term_regex); None if empty."""
    pats = [v1.term_regex(n).pattern for n in names if n]
    return re.compile("|".join(pats), re.I) if pats else None


# ----------------------------------------------------------------------------- MONDO / PrimeKG / Orphanet names
def load_mondo(path: Path) -> Tuple[Dict[str, dict], Optional[str]]:
    """mondo.obo -> ({MONDO:id: {name, omim, omimps, orpha (equivalent xrefs only), synonyms (EXACT, untyped)}}, data-version).
    Obsolete terms are dropped."""
    terms: Dict[str, dict] = {}
    version = None
    cur: Optional[dict] = None

    def close() -> None:
        if cur and cur["id"].startswith("MONDO:") and not cur["obsolete"]:
            terms[cur["id"]] = cur

    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("data-version:"):
            version = line.split(":", 1)[1].strip()
        elif line.startswith("["):
            close()
            cur = {"id": "", "name": "", "obsolete": False, "omim": set(), "omimps": set(), "orpha": set(), "synonyms": []} \
                if line.strip() == "[Term]" else None
        elif cur is None:
            continue
        elif line.startswith("id: "):
            cur["id"] = line[4:].strip()
        elif line.startswith("name: "):
            cur["name"] = line[6:].strip()
        elif line.startswith("is_obsolete: true"):
            cur["obsolete"] = True
        elif line.startswith("synonym: "):
            m = SYN_RE.match(line)
            if m and m.group(2) == "EXACT" and not m.group(3):  # typed synonyms (abbreviation, ambiguous...) are dropped
                cur["synonyms"].append(m.group(1).replace('\\"', '"'))
        elif line.startswith("xref: "):
            m = XREF_RE.match(line)
            if m and "MONDO:equivalentTo" in (m.group(3) or ""):
                cur[{"OMIM": "omim", "OMIMPS": "omimps", "Orphanet": "orpha"}[m.group(1)]].add(m.group(2))
    close()
    return terms, version


def mondo_id(n: str) -> str:
    return "MONDO:" + n.zfill(7)


def load_primekg_disease_ids(nodes_tab: Path) -> Dict[str, List[str]]:
    """PrimeKG disease node name -> sorted MONDO ids (node_id is one number or several joined by '_')."""
    import csv
    out: Dict[str, Set[str]] = defaultdict(set)
    with nodes_tab.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["node_type"] == "disease":
                out[r["node_name"]].update(mondo_id(x) for x in r["node_id"].split("_") if x.isdigit())
    return {k: sorted(v) for k, v in out.items()}


def parse_orphanet_names(xml: bytes) -> Dict[str, dict]:
    """Orphadata en_product1: 'ORPHA:code' -> {name, synonyms}."""
    out: Dict[str, dict] = {}
    for d in ET.fromstring(xml).iter("Disorder"):
        code = (d.findtext("OrphaCode") or "").strip()
        if code:
            out[f"ORPHA:{code}"] = {"name": (d.findtext("Name") or "").strip(),
                                    "synonyms": [s.text.strip() for s in d.findall("./SynonymList/Synonym") if s.text]}
    return out


# ----------------------------------------------------------------------------- matching (limitation 1)
class DisorderIndex:
    """Precomputed (id, name) lists of HPO and Orphanet disorders plus a normalised-name index."""

    def __init__(self, hpo: Dict[str, dict], orph: Dict[str, dict]):
        self.hpo, self.orph = hpo, orph
        self.items: List[Tuple[str, str]] = [(i, d["name"]) for i, d in sorted(hpo.items())] + \
                                             [(i, d["name"]) for i, d in sorted(orph.items())]
        self.by_norm: Dict[str, Set[str]] = defaultdict(set)
        for i, n in self.items:
            self.by_norm[v1.norm_name(n)].add(i)

    def inheritance(self, i: str) -> List[str]:
        return list(self.hpo.get(i, {}).get("inheritance", [])) + list(self.orph.get(i, {}).get("inheritance", []))

    def genes(self, i: str) -> Set[str]:
        return set(self.hpo.get(i, {}).get("genes", set())) | set(self.orph.get(i, {}).get("genes", set()))

    def name(self, i: str) -> str:
        return (self.hpo.get(i) or self.orph.get(i) or {}).get("name", "")


def disease_xrefs(mondo_ids: Sequence[str], mondo: Dict[str, dict]) -> Dict[str, Any]:
    """Identifiers and names of a disease from its MONDO ids (ids absent from `mondo`, e.g. obsolete, are reported)."""
    omim, omimps, orpha, names, missing = set(), set(), set(), [], []
    for m in mondo_ids:
        t = mondo.get(m)
        if t is None:
            missing.append(m)
            continue
        omim |= t["omim"]
        omimps |= t["omimps"]
        orpha |= t["orpha"]
        names += [t["name"]] + list(t["synonyms"])
    return {"mondo": sorted(mondo_ids), "mondo_missing": sorted(missing), "omim": sorted(omim), "omimps": sorted(omimps),
            "orphanet": sorted(orpha), "names": names}


def match_disease_v2(name: str, mondo_ids: Sequence[str], mondo: Dict[str, dict], index: DisorderIndex,
                     orph_names: Dict[str, dict]) -> Dict[str, Any]:
    """Matched HPO/Orphanet disorders of a disease by identifiers and synonyms (protocol section 1)."""
    x = disease_xrefs(mondo_ids, mondo)
    ids = {f"OMIM:{n}" for n in x["omim"]} | {f"ORPHA:{n}" for n in x["orphanet"]}
    id_matched = sorted(i for i in ids if i in index.hpo or i in index.orph)
    names = [name] + x["names"]
    for o in x["orphanet"]:
        rec = orph_names.get(f"ORPHA:{o}")
        if rec:
            names += [rec["name"]] + rec["synonyms"]
    names += [index.name(i) for i in id_matched]
    name_set = build_name_set(names)
    rx = names_regex(name_set)
    contained = {i for i, n in index.items if rx and rx.search(n)}
    matched = sorted(contained | set(id_matched))
    # identity disorders: id-matched (one entity by MONDO equivalence) + disorders of exactly a name of the disease
    exact_ids = sorted({i for n in name_set for i in index.by_norm.get(v1.norm_name(n), ())})
    mondo_orpha = bool(x["orphanet"])
    mend = lambda inh: any(v1.MENDELIAN.search(s) for s in inh)  # noqa: E731
    id_inh = [s for i in id_matched for s in index.inheritance(i)]
    rare_by_id = bool(id_matched) and (mondo_orpha or any(i.startswith("ORPHA:") for i in id_matched)) and mend(id_inh)
    rare_by_name = any(i.startswith("ORPHA:") and mend(index.inheritance(i)) for i in exact_ids)
    inh = sorted({s for i in matched for s in index.inheritance(i)})
    genes = sorted({g for i in matched for g in index.genes(i) if v1.SYMBOL_RE.match(g)})
    return {
        "mondo": x["mondo"], "mondo_missing": x["mondo_missing"], "omim": x["omim"], "omimps": x["omimps"],
        "orphanet_xref": x["orphanet"], "name_set": name_set,
        "hpo": sorted(i for i in matched if i in index.hpo), "orphanet": sorted(i for i in matched if i in index.orph),
        "id_matched": id_matched, "n_matched": len(matched), "genes": genes, "inheritance": inh,
        "rare_genetic_identity": rare_by_id or rare_by_name,
        "orphanet_catalogued": mondo_orpha or any(i.startswith("ORPHA:") for i in matched),
        "mendelian": mend(inh),
    }


# ----------------------------------------------------------------------------- selection
def evaluate_candidate_v2(name: str, gold1_names: Set[str], match: Optional[Dict[str, Any]], uidx: Dict[str, dict],
                          curated: bool = False) -> Tuple[bool, str]:
    """(eligible, reason): the v1 checks in the v1 order; only the matching and the rare-genetic test differ."""
    if v1.APPLICATION_CASES.search(name):
        return False, "application case (MVA / Fanconi anaemia) excluded by design"
    if match is None or match["n_matched"] == 0:
        return False, "no HPO or Orphanet disorder matched by identifier, name or synonym"
    if match["n_matched"] > MAX_MATCHED:
        return False, f"umbrella term: {match['n_matched']} matched disorders (> {MAX_MATCHED})"
    if not match["rare_genetic_identity"] and not (curated and match["orphanet_catalogued"] and match["mendelian"]):
        return False, ("not a rare genetic disease: no Orphanet-catalogued identity disorder (identifier or exact "
                       "name/synonym) with Mendelian/mitochondrial/sex-linked inheritance (and not a curated candidate)")
    if len(match["genes"]) < MIN_GENES:
        return False, f"fewer than {MIN_GENES} genes ({len(match['genes'])})"
    gold = {g for g in gold1_names if g in uidx}
    if len(gold) < MIN_GOLD:
        return False, f"fewer than {MIN_GOLD} approved universe drugs indicated ({len(gold)})"
    bad = sorted(g for g in gold if uidx[g].get("molecule_type") in v1.BIOLOGIC_THERAPY_TYPES)
    if bad:
        return False, "gold drug is an enzyme replacement / gene / cell therapy: " + ", ".join(bad)
    return True, "eligible"


def draw_new(eligible: Sequence[str], n_v1: int, seed: int = SEED, cap_total: int = CAP_TOTAL) -> List[str]:
    """All eligible new diseases, or a seeded draw of (cap_total - n_v1) from the sorted names when there are more."""
    pool = sorted(eligible)
    k = max(cap_total - n_v1, 0)
    return sorted(random.Random(seed).sample(pool, k)) if len(pool) > k else pool


def drop_duplicate_gold_v2(eligible: Sequence[str], gold: Dict[str, Set[str]], v1_sets: Dict[str, frozenset]
                           ) -> Tuple[List[str], Dict[str, str]]:
    """A new disease whose gold 1 set equals that of a retained v1 disease is dropped; among new ones the
    alphabetically first is kept."""
    keep, dropped = [], {}
    seen: Dict[frozenset, str] = {s: n for n, s in sorted(v1_sets.items())}
    for n in sorted(eligible):
        key = frozenset(gold[n])
        if key in seen:
            dropped[n] = seen[key]
        else:
            seen[key] = n
            keep.append(n)
    return keep, dropped


# ----------------------------------------------------------------------------- gold 3 (openFDA labels)
def _openfda_getter(url: str, params: dict) -> dict:
    """GET with exponential backoff (SSL/transient/429); 404 = no label; <= 3 requests per second; reduced response."""
    import requests
    last: Optional[Exception] = None
    for attempt in range(7):
        time.sleep(0.34)
        try:
            r = requests.get(url, params=params, timeout=60, headers={"User-Agent": "hselith-track2/0.3 (research)"})
            if r.status_code == 404:
                return {"total": 0, "labels": []}
            if r.status_code == 429:
                time.sleep(5 * (attempt + 1))
                continue
            r.raise_for_status()
            return reduce_openfda(r.json())
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"openFDA: no response for {params} ({last!r})")


def reduce_openfda(raw: dict) -> dict:
    """Keep only what gold 3 needs from an openFDA response (the cache would otherwise hold whole labels)."""
    labels = []
    for r in raw.get("results") or []:
        ind = r.get("indications_and_usage") or []
        labels.append({"set_id": r.get("set_id"), "effective_time": r.get("effective_time"),
                       "generic_name": (r.get("openfda") or {}).get("generic_name", []),
                       "indications": " ".join(ind) if isinstance(ind, list) else str(ind)})
    return {"total": (raw.get("meta") or {}).get("results", {}).get("total", len(labels)), "labels": labels}


def fetch_labels(drug: str, getter: Optional[Callable[[str, dict], Any]] = None,
                 cache_dir: Path = CACHE / "openfda") -> List[dict]:
    """Up to LABELS_PER_DRUG most recent current labels of a drug (generic or substance name phrase); cached."""
    q = drug.replace('"', " ").strip()
    params = {"search": f'openfda.generic_name:"{q}" OR openfda.substance_name:"{q}"', "limit": LABELS_PER_DRUG,
              "sort": "effective_time:desc"}
    return v1.http_get_json(OPENFDA, params, cache_dir, getter or _openfda_getter).get("labels", [])


def gold3_for_disease(phrases: Sequence[str], labels_by_drug: Dict[str, List[dict]]) -> Dict[str, List[dict]]:
    """lower-case drug name -> [{set_id, effective_time, matched_phrase}] for labels whose indications contain a phrase."""
    out: Dict[str, List[dict]] = {}
    rx = names_regex(phrases)
    if rx is None:
        return out
    for drug, labels in sorted(labels_by_drug.items()):
        hits = []
        for lab in labels:
            m = rx.search(lab.get("indications") or "")
            if m:
                hits.append({"set_id": lab.get("set_id"), "effective_time": lab.get("effective_time"),
                             "matched_phrase": m.group(0).lower()})
        if hits:
            out[drug] = hits
    return out


# ----------------------------------------------------------------------------- gold 2 (experimental arms only)
def parse_study_v2(study: dict) -> dict:
    ps = study.get("protocolSection", {})
    ai = ps.get("armsInterventionsModule") or {}
    p = v1.parse_study(study)
    p["arms"] = [{"label": a.get("label", ""), "type": a.get("type", ""),
                  "interventions": list(a.get("interventionNames") or [])} for a in ai.get("armGroups") or []]
    return p


def _arm_texts(arm_item: str, by_name: Dict[str, dict]) -> Tuple[List[str], bool]:
    """Name texts of an arm's intervention entry ('Drug: X' -> X plus the other names of intervention X); skip flag."""
    stripped = re.sub(r"^[^:]+:\s*", "", arm_item).strip()
    iv = by_name.get(stripped.lower())
    if iv is None:
        return [stripped], False
    return [iv["name"]] + list(iv["other"]), iv["type"] in GOLD2_SKIP_TYPES


def gold2_from_studies_v2(studies: Sequence[dict], dre: "re.Pattern[str]") -> Dict[str, Any]:
    """Gold 2 with arm information.

    Returns {'gold': drug -> [{nct, start_year, phases, arm_type, arm_label}] (EXPERIMENTAL arms only),
             'excluded': drug -> {'arm_types': [...], 'ncts': [...]} for drugs seen only in non-EXPERIMENTAL arms or only in
             studies without arm groups, 'n_qualifying_studies', 'n_qualifying_without_arms'}."""
    gold: Dict[str, Dict[str, dict]] = defaultdict(dict)
    other: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    n_q = n_noarm = 0
    for s in studies:
        p = parse_study_v2(s)
        if not p["nct"] or not v1.study_qualifies(p):
            continue
        n_q += 1
        by_name = {i["name"].lower(): i for i in p["interventions"]}
        if not p["arms"]:
            n_noarm += 1
            for iv in p["interventions"]:
                if iv["type"] in GOLD2_SKIP_TYPES:
                    continue
                for text in [iv["name"]] + iv["other"]:
                    for m in dre.finditer(text or ""):
                        other[m.group(0).lower()][p["nct"]].add("NO_ARM_GROUPS")
            continue
        for arm in p["arms"]:
            for item in arm["interventions"]:
                texts, skip = _arm_texts(item, by_name)
                if skip:
                    continue
                for text in texts:
                    for m in dre.finditer(text or ""):
                        d = m.group(0).lower()
                        if arm["type"] == EXPERIMENTAL:
                            gold[d].setdefault(p["nct"], {"nct": p["nct"], "start_year": p["start"].year,
                                                          "phases": sorted(p["phases"]), "arm_type": EXPERIMENTAL,
                                                          "arm_label": arm["label"]})
                        else:
                            other[d][p["nct"]].add(arm["type"] or "NO_ARM_TYPE")
    excluded = {d: {"arm_types": sorted({t for ts in v.values() for t in ts}), "ncts": sorted(v)}
                for d, v in sorted(other.items()) if d not in gold}
    return {"gold": {d: sorted(v.values(), key=lambda t: t["nct"]) for d, v in sorted(gold.items())},
            "excluded": excluded, "n_qualifying_studies": n_q, "n_qualifying_without_arms": n_noarm}


GOLD2_FIELDS = ("NCTId|StudyType|Phase|StartDate|InterventionName|InterventionOtherName|InterventionType|"
                "ArmGroupLabel|ArmGroupType|ArmGroupInterventionName")


def fetch_studies_v2(condition: str, getter: Optional[Callable[[str, dict], Any]] = None) -> List[dict]:
    """Interventional phase 2+ studies starting after 2018-12-31 for a condition search, with arm groups (paged, cached)."""
    studies: List[dict] = []
    token = None
    while True:
        params = {
            "query.cond": condition, "pageSize": 1000, "format": "json",
            "filter.advanced": "AREA[StudyType]INTERVENTIONAL AND AREA[Phase](PHASE2 OR PHASE3 OR PHASE4) AND "
                               f"AREA[StartDate]RANGE[{v1.GOLD2_AFTER.isoformat()},MAX]",
            "fields": GOLD2_FIELDS,
        }
        if token:
            params["pageToken"] = token
        data = v1.http_get_json(v1.CTGOV, params, CACHE / "ctgov_v2", getter)
        studies.extend(data.get("studies") or [])
        token = data.get("nextPageToken")
        if not token:
            return studies


# ----------------------------------------------------------------------------- genes
def _default_annotate(sym: str) -> dict:
    from track2 import disease_landscape as dl
    old = dl.CACHE
    dl.CACHE = CACHE / "landscape"  # same cache as v1; disease_landscape itself is not edited
    try:
        return dl.uniprot_annotation(sym)
    finally:
        dl.CACHE = old


def gene_record(symbol: str, annotate: Optional[Callable[[str], dict]] = None) -> dict:
    """Gene with the v1 narrow class (disease_landscape.classify_terms) and mechanism_class_general."""
    from track2 import disease_landscape as dl
    try:
        ann = (annotate or _default_annotate)(symbol)
    except Exception as e:  # noqa: BLE001 - keep the gene, class 'other', failure recorded
        ann = {"terms": [], "accession": None, "error": repr(e)[:200]}
    narrow, _ = dl.classify_terms(ann["terms"])
    general, matched = classify_general(ann["terms"])
    rec = {"symbol": symbol, "uniprot_accession": ann.get("accession"), "mechanism_class": narrow,
           "mechanism_class_general": general, "general_matched_terms": matched[:12], "n_terms": len(ann["terms"])}
    if ann.get("error"):
        rec["annotation_error"] = ann["error"]
    return rec


# ----------------------------------------------------------------------------- driver
def _sizes(g1: Set[str], g3: Set[str]) -> Dict[str, int]:
    return {"gold1": len(g1), "gold3": len(g3), "union": len(g1 | g3), "intersection": len(g1 & g3)}


def build_disease(name: str, origin: str, genes: Sequence[str], match: Dict[str, Any], gold1: Set[str],
                  uidx: Dict[str, dict], dre: "re.Pattern[str]", labels_by_drug: Dict[str, List[dict]],
                  studies_fn: Callable[[str], List[dict]], annotate: Optional[Callable[[str], dict]],
                  inheritance: Sequence[str]) -> dict:
    def drug(d: str) -> dict:
        return {"name": d, "drug_id": uidx[d]["drug_id"], "molecule_type": uidx[d].get("molecule_type"),
                "atc": uidx[d].get("atc", []), "atc_class": v1.atc_class(uidx[d].get("atc", []))}

    g3 = gold3_for_disease(sorted(set(match["name_set"]) | {clean_name(name)}), labels_by_drug)
    g3 = {d: h for d, h in g3.items() if d in uidx}
    studies = studies_fn(name)
    g2 = gold2_from_studies_v2(studies, dre)
    g2g = {d: t for d, t in g2["gold"].items() if d in uidx}
    g2x = {d: x for d, x in g2["excluded"].items() if d in uidx}
    g1s, g3s = set(gold1), set(g3)
    return {
        "name": name, "origin": origin,
        "identifiers": {"primekg_name": name, "mondo": match["mondo"], "mondo_missing": match["mondo_missing"],
                        "omim_xref": match["omim"], "omimps_xref": match["omimps"], "orphanet_xref": match["orphanet_xref"],
                        "hpo": match["hpo"], "orphanet": match["orphanet"], "id_matched": match["id_matched"],
                        "name_set": match["name_set"]},
        "inheritance": list(inheritance),
        "genes": [gene_record(s, annotate) for s in genes],
        "gold1": [dict(drug(d)) for d in sorted(g1s)],
        "gold3": [dict(drug(d), in_gold1=d in g1s, labels=h) for d, h in sorted(g3.items())],
        "gold2": [dict(drug(d), in_gold1=d in g1s, in_gold3=d in g3s, trials=t) for d, t in sorted(g2g.items())],
        "gold2_excluded_non_experimental": {d: x for d, x in sorted(g2x.items())},
        "gold_sizes": dict(_sizes(g1s, g3s), gold2=len(g2g), gold2_excluded=len(g2x)),
        "gold3_biologic_types": sorted(d for d in g3s if uidx[d].get("molecule_type") in v1.BIOLOGIC_THERAPY_TYPES),
        "gold2_query": {"endpoint": v1.CTGOV, "query.cond": name, "n_studies_returned": len(studies),
                        "n_qualifying_studies": g2["n_qualifying_studies"],
                        "n_qualifying_without_arm_groups": g2["n_qualifying_without_arms"]},
    }


def build(root: Path = ROOT, hpo: Optional[Dict[str, dict]] = None, orph: Optional[Dict[str, dict]] = None,
          orph_names: Optional[Dict[str, dict]] = None, mondo: Optional[Dict[str, dict]] = None,
          mondo_version: Optional[str] = None, primekg_ids: Optional[Dict[str, List[str]]] = None,
          universe: Optional[List[dict]] = None, gold1: Optional[Dict[str, Set[str]]] = None,
          v1_benchmark: Optional[dict] = None, curated: Optional[Set[str]] = None,
          labels_by_drug: Optional[Dict[str, List[dict]]] = None,
          studies_fn: Callable[[str], List[dict]] = fetch_studies_v2,
          annotate: Optional[Callable[[str], dict]] = None) -> Dict[str, Any]:
    drugs = universe if universe is not None else v1.load_universe(root / "data" / "universe.json")
    uidx = v1.universe_index(drugs)
    gold1 = gold1 if gold1 is not None else v1.primekg_gold1(root, set(uidx))
    v1b = v1_benchmark if v1_benchmark is not None else json.loads((root / "out" / "engine" / "benchmark.json").read_text(encoding="utf-8"))
    if hpo is None:
        src = root / "data" / "sources"
        hpo = v1.load_hpo(next(src.rglob("phenotype.hpoa")), next(src.rglob("genes_to_phenotype.txt")), next(src.rglob("hp.obo")))
    if orph is None:
        orph = v1.fetch_orphanet()
    if orph_names is None:
        orph_names = parse_orphanet_names(v1._cached_bytes(ORPHADATA_NAMES))
    if mondo is None:
        mondo, mondo_version = load_mondo(root / "data" / "sources" / "mondo" / "mondo.obo")
    if primekg_ids is None:
        primekg_ids = load_primekg_disease_ids(root / "data" / "sources" / "primekg" / "nodes.tab")
    if curated is None:
        from track2 import empirical_literature_control as elc
        curated = {c[0] for c in elc.CANDIDATES}
    index = DisorderIndex(hpo, orph)
    dre = v1.drug_name_regex(sorted(uidx))
    v1_names = [d["name"] for d in v1b["diseases"]]
    v1_by_name = {d["name"]: d for d in v1b["diseases"]}
    v1_reason = {e["disease"]: e["reason"] for e in v1b.get("selection_log", [])}
    v1_sets = {n: frozenset(g["name"] for g in v1_by_name[n]["gold1"]) for n in v1_names}

    pool = sorted({d for d, g in gold1.items() if len(g) >= MIN_GOLD} | (curated & set(gold1)))
    log, eligible, matches = [], [], {}
    for name in pool:
        if name in v1_by_name:
            log.append({"disease": name, "eligible": False, "retained_from_v1": True, "reason": "retained from v1",
                        "v1_reason": v1_reason.get(name)})
            continue
        m = match_disease_v2(name, primekg_ids.get(name, []), mondo, index, orph_names)
        ok, why = evaluate_candidate_v2(name, gold1[name], m, uidx, name in curated)
        matches[name] = m
        if ok:
            eligible.append(name)
        log.append({"disease": name, "eligible": ok, "retained_from_v1": False, "reason": why,
                    "v1_reason": v1_reason.get(name), "n_mondo_ids": len(m["mondo"])})
    eligible, dup = drop_duplicate_gold_v2(eligible, gold1, v1_sets)
    chosen = draw_new(eligible, len(v1_names))
    for e in log:
        if e["disease"] in dup:
            e["eligible"], e["reason"] = False, f"duplicate gold 1 set of {dup[e['disease']]} (not independent)"
        elif e["eligible"]:
            e["reason"] = "selected (new in v2)" if e["disease"] in chosen else "eligible, not drawn (cap)"
    # labels for gold 3: every universe drug (cached; only needed when gold 3 is computed)
    if labels_by_drug is None:
        labels_by_drug = {d: fetch_labels(d) for d in sorted(uidx)}
    diseases = []
    for n in v1_names:
        d1 = v1_by_name[n]
        m = match_disease_v2(n, primekg_ids.get(n, []), mondo, index, orph_names)
        diseases.append(build_disease(n, "v1", [g["symbol"] for g in d1["genes"]], m, set(v1_sets[n]), uidx, dre,
                                      labels_by_drug, studies_fn, annotate, d1.get("inheritance", [])))
    for n in chosen:
        m = matches[n]
        diseases.append(build_disease(n, "v2_new", m["genes"], m, set(gold1[n]), uidx, dre, labels_by_drug, studies_fn,
                                      annotate, m["inheritance"]))
    diseases.sort(key=lambda d: (d["origin"] != "v1", d["name"]))
    reasons: Dict[str, int] = defaultdict(int)
    for e in log:
        reasons[re.sub(r"[:(].*$", "", e["reason"]).strip()] += 1

    def agg(sel: Sequence[dict]) -> Dict[str, Any]:
        keys = ("gold1", "gold3", "union", "intersection", "gold2")
        sums = {k: sum(d["gold_sizes"][k] for d in sel) for k in keys}
        return {"n_diseases": len(sel), "sum_sizes": sums,
                "n_empty": {k: sum(1 for d in sel if d["gold_sizes"][k] == 0) for k in keys}}

    cls_pairs = Counter(g["mechanism_class_general"] for d in diseases for g in d["genes"])
    uniq = {g["symbol"]: g["mechanism_class_general"] for d in diseases for g in d["genes"]}
    return {
        "version": "v2", "seed": SEED, "cap_total": CAP_TOTAL, "target_total": TARGET_TOTAL, "cutoff_year_Y": v1.CUTOFF_YEAR,
        "n_benchmark": len(diseases), "n_v1_retained": len(v1_names), "n_v2_new": len(chosen),
        "n_candidates_evaluated": len(pool), "n_eligible_new": len(eligible), "target_met": len(diseases) >= TARGET_TOTAL,
        "n_universe_drugs": len(uidx),
        "sources": {"mondo_data_version": mondo_version, "mondo_file": "data/sources/mondo/mondo.obo",
                    "orphadata": [v1.ORPHADATA["gene_disease"], v1.ORPHADATA["inheritance"], ORPHADATA_NAMES],
                    "openfda": OPENFDA, "clinicaltrials": v1.CTGOV, "n_drugs_with_labels": sum(1 for v in labels_by_drug.values() if v),
                    "n_drugs_queried_openfda": len(labels_by_drug), "labels_per_drug_max": LABELS_PER_DRUG},
        "gold_rules": {
            "gold1": "universe drug with a PrimeKG `indication` edge to the disease (unchanged from v1)",
            "gold3": "universe drug whose current FDA label (openFDA indications_and_usage, up to 5 most recent labels) contains "
                     "the disease name or an exact synonym as a whole phrase (case-insensitive, no negation handling)",
            "gold2": "universe drug that is an intervention of an EXPERIMENTAL arm group of an interventional ClinicalTrials.gov "
                     "study of the disease, phase 2 or later, start date after 2018-12-31",
        },
        "summary": {"all": agg(diseases), "v1": agg([d for d in diseases if d["origin"] == "v1"]),
                    "v2_new": agg([d for d in diseases if d["origin"] == "v2_new"])},
        "mechanism_class_general_distribution": {"gene_disease_pairs": {c: cls_pairs.get(c, 0) for c in CLASSES},
                                                 "unique_genes": {c: sum(1 for v in uniq.values() if v == c) for c in CLASSES}},
        "diseases": diseases,
        "selection_log": sorted(log, key=lambda e: e["disease"]),
        "selection_reason_counts": dict(sorted(reasons.items())),
        "limitations": [
            "Gold 3 is a text match of the disease name or an exact synonym in the indications of up to 5 current labels per drug; negation and limitations of use are not handled and labels absent from openFDA are missed.",
            "Gold 2 uses the PrimeKG disease name as the ClinicalTrials.gov condition search only; a drug in an EXPERIMENTAL arm may still be background therapy of an add-on design.",
            "HPO (pinned files) has no synonym table: HPO disorders enter by identifier and by name containment only.",
            "Rare genetic for the new diseases = Orphanet-catalogued identity disorder (identifier mapped by MONDO or exact name/synonym) with Mendelian/mitochondrial/sex-linked inheritance; gene lists are the union over all matched subtypes.",
            "Gold 1 remains incomplete and partly circular; the v1-only result stays the primary pre-registered result.",
        ],
    }


def main() -> None:
    out = build()
    p = ROOT / "out" / "engine" / "benchmark_v2.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: out[k] for k in ("n_benchmark", "n_v1_retained", "n_v2_new", "n_eligible_new", "target_met",
                                          "summary", "mechanism_class_general_distribution", "selection_reason_counts")},
                     indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
