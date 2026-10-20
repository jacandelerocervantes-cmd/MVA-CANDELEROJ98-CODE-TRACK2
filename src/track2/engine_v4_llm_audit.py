r"""READ-ONLY, DEV-ONLY audit of the claim "at low literature volume there are almost no relevant sentences, so LLM-read evidence cannot be
tested" (protocol/addendum_engine_v4_llm_audit.md, written before this module ran).

Run (PYTHONPATH=src):
    python -m track2.engine_v4_llm_audit sample    prints the Q3a sentences to be read (no network)
    python -m track2.engine_v4_llm_audit run       full audit; Europe PMC only for Q2 (<= 4 req/s, cached under data/cache_engine/audit/)
    python -m track2.engine_v4_llm_audit md        re-renders the markdown from the JSON (reading and verdict prose edited in this module)

It never reads HELD-OUT labels, gold, shards or annotation files, runs no pipeline stage, edits no existing file, sends nothing to an LLM.
Multi-split files (candidate file, pair index, benchmark) are filtered to DEV diseases at load time.
"""
from __future__ import annotations

import hashlib
import json
import random
import re
import statistics
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from track2 import empirical_literature_control as elc
from track2 import engine_benchmark as v1
from track2 import engine_extractors as ex
from track2 import engine_score as es
from track2 import engine_v4 as v4
from track2 import engine_v4_dilution as dl
from track2 import engine_v4_llm as lm

ROOT = v4.ROOT
SEED = 20261015
OUT = v4.OUT
PROTOCOL = ROOT / "protocol" / "addendum_engine_v4_llm_audit.md"
AUDIT_CACHE = v4.CACHE / "audit"
OUT_JSON = OUT / "audit_llm_coverage_dev.json"
OUT_MD = OUT / "audit_llm_coverage_dev.md"
LABEL = "READ-ONLY, DEV-ONLY audit of the low-coverage claim (post hoc); no new annotation; opinion tables are audit readings, not ground truth"

N_SAMPLE_PAIRS = 150
N_READ = 60
N_LANG_DISEASES = 30
PAGE_FIRST = 20
R6_PAPERS = 10
MAX_DIS_TERMS = 6
PILOT_SHARDS = (2, 4, 6, 25, 26, 28, 57, 62)
K_PAPERS = lm.MAX_PAPERS_PER_PAIR
K_SENT = lm.MAX_SENT_PER_PAPER
SALTS = {"hydrochloride", "hcl", "sodium", "potassium", "calcium", "magnesium", "sulfate", "sulphate", "phosphate", "acetate", "maleate", "mesylate",
         "tartrate", "citrate", "succinate", "besylate", "fumarate", "bromide", "chloride", "hydrobromide", "nitrate", "carbonate", "lactate",
         "propionate", "dipropionate", "valerate", "pivalate", "benzoate", "diacetate", "decanoate", "cypionate", "palmitate", "furoate"}
FORMS = {"anhydrous", "free", "hydrate", "monohydrate", "dihydrate"}
STAGES = ["S0_candidate_pairs", "S1a_drug_named_in_a_kept_paper", "S1b_and_paper_has_disease_gene", "S2_sentence_extracted", "S3_verified_sentence",
          "S4_relevance_yes_either", "S5_relevance_yes_both", "S6_yes_both_and_benefit"]
SIDE = ["S6b_yes_either_and_benefit", "opus_yes", "sonnet_yes", "opus_yes_benefit", "sonnet_yes_benefit", "benefit_direction_only"]
VARIANTS = ["R0_baseline", "R1_sentence_cap_lifted", "R2_paper_cap_lifted", "R3_both_caps_lifted", "R4_R3_no_gene_gate", "R5_R3_drug_synonyms",
            "R6a_direct_query_top10_INN", "R6b_direct_query_top10_synonyms", "R7_union_R4_R5_R6b"]

UNVALIDATED: Dict[str, Dict[str, Any]] = {
    "sample_gold_positive_pairs": {"value": N_SAMPLE_PAIRS, "status": "UNVALIDATED", "conservative_default": N_SAMPLE_PAIRS, "note": "audit sample size fixed by the task"},
    "sentences_read": {"value": N_READ, "status": "UNVALIDATED", "conservative_default": N_READ, "note": "audit reading size fixed by the task"},
    "direct_query_papers": {"value": R6_PAPERS, "status": "UNVALIDATED", "conservative_default": R6_PAPERS, "note": "top 10 papers of the drug AND disease query"},
    "direct_query_first_records": {"value": PAGE_FIRST, "status": "UNVALIDATED", "conservative_default": PAGE_FIRST, "note": "recall measured on the first 20 records"},
    "disease_terms_in_query": {"value": MAX_DIS_TERMS, "status": "UNVALIDATED", "conservative_default": 1, "note": "name plus the shortest synonyms of the benchmark name set"},
    "language_sample_diseases": {"value": N_LANG_DISEASES, "status": "UNVALIDATED", "conservative_default": N_LANG_DISEASES, "note": "seeded sample"},
    "material_loss_rule": {"value": "relaxing one rule alone doubles drug+disease sentences per pair or adds >= 20 points of pairs with >= 1; language material if > 5 % non-English",
                           "status": "UNVALIDATED", "conservative_default": "same", "note": "audit convention"},
    "salt_words": {"value": sorted(SALTS | FORMS), "status": "UNVALIDATED", "conservative_default": sorted(SALTS | FORMS), "note": "synonym rule 3"},
    "inherited_constants": {"value": "top 30, 3 papers, 2 sentences, 300 records, year <= 2024", "status": "UNVALIDATED", "conservative_default": "as parent",
                            "note": "parent registry"},
}


# ============================================================================= small pure helpers
def salt_stripped(name: str) -> Optional[str]:
    """Lower-case name without trailing salt / form words; None if nothing is stripped or nothing remains."""
    toks = name.lower().split()
    n = len(toks)
    while len(toks) > 1 and toks[-1].strip(",") in SALTS | FORMS:
        toks.pop()
    if len(toks) == n or not toks:
        return None
    return " ".join(toks).strip(", ")


def is_plain_inn(name: str) -> bool:
    """A project name that is a bare INN: no salt / form word, no bracket, comma or digit."""
    low = name.lower()
    if any(c in low for c in "(),;") or any(c.isdigit() for c in low):
        return False
    return not any(t in SALTS | FORMS for t in low.split())


def name_variants(name: str) -> Set[str]:
    low = name.lower()
    out = set()
    if "-" in low:
        out.add(low.replace("-", " "))
    if " " in low:
        out.add(low.replace(" ", "-"))
    return {v for v in out if len(v) >= 4}


def build_synonyms(names: Sequence[str], ctgov_files: Iterable[Path], chembl_files: Iterable[Path]) -> Dict[str, Dict[str, Set[str]]]:
    """project name (lower) -> {term: {sources}}; sources: ctgov_other_names, chembl_synonyms, salt_stripped, hyphen_space."""
    known = {n.lower() for n in names}
    stripped = {}
    for n in known:
        s = salt_stripped(n)
        if s:
            stripped.setdefault(s, set()).add(n)
    syn: Dict[str, Dict[str, Set[str]]] = {}

    def add(nm: str, term: str, src: str) -> None:
        t = term.strip().lower()
        if len(t) >= 4 and t != nm:
            syn.setdefault(nm, {}).setdefault(t, set()).add(src)

    for n in known:
        s = salt_stripped(n)
        if s:
            add(n, s, "salt_stripped")
        for v in name_variants(n):
            add(n, v, "hyphen_space")
    for f in ctgov_files:
        try:
            doc = json.loads(Path(f).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        for st in doc.get("studies", []) or []:
            for iv in ((st.get("protocolSection") or {}).get("armsInterventionsModule") or {}).get("interventions", []) or []:
                nm = (iv.get("name") or "").lower().strip()
                targets = ({nm} if nm in known else set()) | stripped.get(nm, set())
                for tgt in targets:
                    if nm != tgt:
                        add(tgt, nm, "ctgov_other_names")
                    for o in iv.get("otherNames") or []:
                        add(tgt, o, "ctgov_other_names")
    for f in chembl_files:
        try:
            t = Path(f).read_text(encoding="utf-8", errors="replace")
            if '"molecule_synonyms"' not in t:
                continue
            doc = json.loads(t)
        except (OSError, ValueError):
            continue
        nm = (doc.get("pref_name") or "").lower()
        if nm in known:
            for s in doc.get("molecule_synonyms") or []:
                add(nm, s.get("molecule_synonym") or "", "chembl_synonyms")
    return syn


def crosstab(a: Sequence[str], b: Sequence[str], vocab: Sequence[str]) -> Dict[str, Any]:
    """Rows a, columns b, with n, raw agreement and Cohen's kappa."""
    tab = {x: {y: 0 for y in vocab} for x in vocab}
    for x, y in zip(a, b):
        tab[x][y] += 1
    n = len(a)
    agree = sum(tab[x][x] for x in vocab)
    return {"rows_first_columns_second": tab, "n": n, "agreement": (agree / n) if n else None,
            "kappa": float(ex.kappa(np.array(a), np.array(b))) if n else None}


def sample_pairs(pairs: Sequence[Tuple[str, str]], n: int, seed: int = SEED) -> List[Tuple[str, str]]:
    return sorted(random.Random(seed).sample(sorted(pairs), min(n, len(pairs))))


def disease_terms(name: str, name_set: Sequence[str]) -> Tuple[List[str], List[str]]:
    """(name-only terms, name + name-set terms) for sentence matching; terms of >= 4 characters."""
    only = [name] if len(name) >= 4 else []
    full = sorted({t.strip() for t in [name] + list(name_set) if t and len(t.strip()) >= 4}, key=lambda t: (-len(t), t))
    return only, full


def query_disease_terms(name: str, name_set: Sequence[str], max_terms: int = MAX_DIS_TERMS) -> List[str]:
    """Name plus the shortest other synonyms (no digits, commas or quotes, >= 5 characters)."""
    low = name.lower()
    extra = sorted({t.strip() for t in name_set if t and t.strip().lower() != low and len(t.strip()) >= 5 and not any(c.isdigit() for c in t) and "," not in t
                    and '"' not in t}, key=lambda t: (len(t), t))
    return [name.replace('"', " ").strip()] + extra[:max_terms - 1]


def clause(terms: Sequence[str]) -> str:
    return "(" + " OR ".join(f'TITLE_ABS:"{t.replace(chr(34), " ").strip()}"' for t in terms) + ")"


def direct_query(drug_terms: Sequence[str], dis_terms: Sequence[str], max_year: int = es.YEAR_PRIMARY, extra: str = "") -> str:
    return f"{clause(drug_terms)} AND {clause(dis_terms)} AND PUB_YEAR:[1800 TO {max_year}] AND HAS_ABSTRACT:y{extra}"


def recall(found: Sequence[str], retrieved: Set[str]) -> Optional[float]:
    return (sum(1 for p in found if p in retrieved) / len(found)) if found else None


def describe(vals: Sequence[float]) -> Dict[str, Any]:
    v = [float(x) for x in vals if x is not None]
    if not v:
        return {"n": 0}
    return {"n": len(v), "mean": float(np.mean(v)), "median": float(np.median(v)), "min": float(min(v)), "max": float(max(v)),
            "p25": float(np.percentile(v, 25)), "p75": float(np.percentile(v, 75))}


def language_summary(langs: Sequence[Optional[str]]) -> Dict[str, Any]:
    c = Counter((l or "unknown") for l in langs)
    n = sum(c.values())
    non = n - c.get("eng", 0)
    return {"n": n, "counts": dict(sorted(c.items(), key=lambda kv: -kv[1])), "non_english": non, "non_english_share": (non / n) if n else None}


# ============================================================================= Europe PMC (audit cache)
def epmc_search(query: str, page_size: int, cursor: str = "*", getter: Optional[Callable[[str, dict], dict]] = None,
                cache_dir: Path = AUDIT_CACHE) -> Dict[str, Any]:
    """One cached Europe PMC call reduced to {hitCount, next, records[{pid, year, language, title, abstract}]}; <= 4 req/s through the engine limiter."""
    params = {"query": query, "format": "json", "resultType": "core", "pageSize": page_size, "cursorMark": cursor}
    f = cache_dir / (hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()[:24] + ".json")
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    data = (getter or es.default_getter)(es.EPMC, params)
    recs = []
    for it in (data.get("resultList") or {}).get("result") or []:
        recs.append({"pid": f"{it.get('source', '')}:{it.get('id') or it.get('pmid') or ''}", "year": it.get("pubYear"), "language": it.get("language"),
                     "title": it.get("title") or "", "abstract": it.get("abstractText") or ""})
    out = {"query": query, "hitCount": int(data.get("hitCount") or 0), "next": data.get("nextCursorMark"), "records": recs}
    cache_dir.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out, ensure_ascii=False), encoding="utf-8")
    return out


def fetch_gene_records(symbols: Sequence[str], max_n: int = es.MAX_PAPERS, getter=None, cache_dir: Path = AUDIT_CACHE) -> Dict[str, Any]:
    """The project's gene query, paged as `engine_score.fetch_papers` does, keeping the language field."""
    query = es.build_gene_query(symbols, es.YEAR_PRIMARY)
    recs, seen, cursor, hits = [], set(), "*", None
    while len(recs) < max_n:
        page = epmc_search(query, 100, cursor, getter, cache_dir)
        hits = page["hitCount"]
        n0 = len(recs)
        for r in page["records"]:
            if not r["abstract"] or r["pid"] in seen:
                continue
            seen.add(r["pid"])
            recs.append(r)
        nxt = page["next"]
        if len(page["records"]) < 100 or not nxt or nxt == cursor or (len(recs) == n0 and not page["records"]):
            break
        cursor = nxt
    return {"query": query, "hitCount": hits, "records": recs[:max_n]}


def record_text(r: dict) -> str:
    return f"{r['title']}. {r['abstract']}"


# ============================================================================= text index of a disease
class DiseaseText:
    """Sentences of the retrieved papers of one disease with the production drug-name matches (cached per sentence)."""

    def __init__(self, papers: Sequence[dict], symbols: Sequence[str], dre: Optional["re.Pattern[str]"], gate: Optional[Set[str]],
                 dn_re: Optional["re.Pattern[str]"], ds_re: Optional["re.Pattern[str]"]):
        grx = es.gene_regex(symbols) if symbols else None
        self.pids = [p["pid"] for p in papers]
        self.texts = [p["text"] for p in papers]
        self.gene = [bool(grx.search(t)) if grx else True for t in self.texts]
        self.sents = [elc.split_sentences(t) for t in self.texts]
        self.prod: List[List[Set[str]]] = []
        for ss in self.sents:
            row = []
            for s in ss:
                if dre is None or (gate is not None and not (set(re.findall(r"\w+", s.lower())) & gate)):
                    row.append(set())
                else:
                    row.append({m.group(0).lower() for m in dre.finditer(s)})
            self.prod.append(row)
        self.dn = [[bool(dn_re and dn_re.search(s)) for s in ss] for ss in self.sents]
        self.ds = [[bool(ds_re and ds_re.search(s)) for s in ss] for ss in self.sents]
        self.dis_paper = [bool(ds_re and ds_re.search(t)) for t in self.texts]
        self.names_any = [set().union(*row) if row else set() for row in self.prod]


def extract_pair(dt: DiseaseText, kept: Iterable[int], match: Callable[[int, int], bool], k_papers: Optional[int] = None, k_sent: Optional[int] = None,
                 gene_gate: bool = True) -> List[Tuple[int, List[int]]]:
    """The production selection with switchable rules: papers in retrieval order, optional gene gate, the first k_sent sentences that match the
    drug, the first k_papers papers that have at least one. No label argument."""
    out: List[Tuple[int, List[int]]] = []
    for i in sorted(kept):
        if gene_gate and not dt.gene[i]:
            continue
        js = [j for j in range(len(dt.sents[i])) if match(i, j)]
        if not js:
            continue
        out.append((i, js[:k_sent] if k_sent is not None else js))
        if k_papers is not None and len(out) >= k_papers:
            break
    return out


def records_of(dt: DiseaseText, ext: Sequence[Tuple[int, List[int]]]) -> List[Dict[str, Any]]:
    return [{"pid": dt.pids[i], "sent": dt.sents[i][j], "dn": dt.dn[i][j], "ds": dt.ds[i][j], "paper_dis": dt.dis_paper[i]} for i, js in ext for j in js]


def records_of_papers(papers: Sequence[dict], drug_re: "re.Pattern[str]", dn_re, ds_re) -> List[Dict[str, Any]]:
    """Records of directly retrieved papers: every sentence that matches the drug."""
    out = []
    for p in papers:
        pdis = bool(ds_re.search(p["text"]))
        for s in elc.split_sentences(p["text"]):
            if drug_re.search(s):
                out.append({"pid": p["pid"], "sent": s, "dn": bool(dn_re.search(s)), "ds": bool(ds_re.search(s)), "paper_dis": pdis})
    return out


def dedupe_records(recs: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    seen, out = set(), []
    for r in recs:
        k = (r["pid"], r["sent"])
        if k not in seen:
            seen.add(k)
            out.append(r)
    return out


def pair_metrics(recs: Sequence[Dict[str, Any]], drug: str) -> Dict[str, int]:
    recs = dedupe_records(recs)
    dd = [r for r in recs if r["ds"]]
    return {"sent_drug": len(recs), "sent_drug_disease_name": sum(r["dn"] for r in recs), "sent_drug_disease_set": len(dd),
            "papers_drug": len({r["pid"] for r in recs}), "papers_drug_and_disease": len({r["pid"] for r in recs if r["paper_dis"]}),
            "sent_drug_disease_evaluative_rule": sum(ex.rule_direction(r["sent"], drug) in ("benefit", "harm") for r in dd)}


def summarise_variant(per_pair: Sequence[Dict[str, float]], thr: float = 1) -> Dict[str, Any]:
    out: Dict[str, Any] = {"n_pairs": len(per_pair)}
    for k in ("sent_drug", "sent_drug_disease_name", "sent_drug_disease_set", "papers_drug", "papers_drug_and_disease", "sent_drug_disease_evaluative_rule"):
        v = [p[k] for p in per_pair]
        out[k] = {"mean": float(np.mean(v)) if v else None, "median": float(np.median(v)) if v else None, "total": int(sum(v)),
                  "pairs_with_ge1": int(sum(1 for x in v if x >= thr)), "share_pairs_with_ge1": float(np.mean([x >= thr for x in v])) if v else None}
    return out


def match_prod(dt: DiseaseText, nm: str) -> Callable[[int, int], bool]:
    return lambda i, j: nm in dt.prod[i][j]


def match_syn(dt: DiseaseText, nm: str, terms: Sequence[str]) -> Callable[[int, int], bool]:
    rx = elc.term_regex(sorted(set(terms) | {nm}))
    return lambda i, j: nm in dt.prod[i][j] or bool(rx.search(dt.sents[i][j]))


def match_terms(dt: DiseaseText, terms: Sequence[str]) -> Callable[[int, int], bool]:
    rx = elc.term_regex(sorted(set(terms)))
    return lambda i, j: bool(rx.search(dt.sents[i][j]))


# ============================================================================= funnel
def label_flags(o_rel: str, s_rel: str, cons_dir: Optional[str]) -> Dict[str, bool]:
    ben = cons_dir == "benefit"
    either, both = o_rel == "yes" or s_rel == "yes", o_rel == "yes" and s_rel == "yes"
    return {"S3_verified_sentence": True, "S4_relevance_yes_either": either, "S5_relevance_yes_both": both, "S6_yes_both_and_benefit": both and ben,
            "S6b_yes_either_and_benefit": either and ben, "opus_yes": o_rel == "yes", "sonnet_yes": s_rel == "yes",
            "opus_yes_benefit": o_rel == "yes" and ben, "sonnet_yes_benefit": s_rel == "yes" and ben, "benefit_direction_only": ben}


def pair_stage_flags(entries: Sequence[dict], kept: Set[int], mention_any: bool, mention_gene: bool, disease: str, drug: str,
                     table: Dict[Tuple[str, str, str, str], Dict[str, bool]], k_papers: int = K_PAPERS) -> Dict[str, bool]:
    """Stage flags of one pair in one sample. entries = parent pair index of the drug (union over samples); the first k_papers kept ones count."""
    fl = {s: False for s in STAGES[1:] + SIDE}
    fl["S0_candidate_pairs"] = True
    fl["S1a_drug_named_in_a_kept_paper"] = mention_any
    fl["S1b_and_paper_has_disease_gene"] = mention_gene
    sel = [e for e in entries if e["i"] in kept][:k_papers]
    fl["S2_sentence_extracted"] = any(e["sents"] for e in sel)
    for e in sel:
        for s in e["sents"]:
            lab = table.get((disease, drug, lm.pmid_of(e["pid"]), s))
            if lab:
                for k, v in lab.items():
                    fl[k] = fl[k] or v
    return fl


def largest_losses(counts: Dict[str, float]) -> Dict[str, Any]:
    chain = [counts[s] for s in STAGES]
    trans = []
    for a, b, sa, sb in zip(chain, chain[1:], STAGES, STAGES[1:]):
        trans.append({"from": sa, "to": sb, "lost_pairs": a - b, "retention": (b / a) if a else None})
    big = max(trans, key=lambda t: t["lost_pairs"])
    low = min((t for t in trans if t["retention"] is not None), key=lambda t: t["retention"], default=None)
    return {"transitions": trans, "largest_absolute_loss": big, "lowest_retention": low}


# ============================================================================= context (DEV only)
def load_context() -> Dict[str, Any]:
    bench = v4.load_benchmark(ROOT)
    dev = list(bench["split"]["dev"])
    dev_set = set(dev)
    rank = {n: i for i, n in enumerate(sorted(d["name"] for d in bench["diseases"]))}
    dis = {}
    for d in bench["diseases"]:
        if d["name"] in dev_set:
            dis[d["name"]] = {"symbols": sorted({g["symbol"] for g in d["genes"]}), "name_set": list(d["identifiers"].get("name_set", []))}
    del bench
    cands_all = json.loads(lm.CAND_FILE.read_text(encoding="utf-8"))
    cands = {d: v for d, v in cands_all.items() if d in dev_set and v["split"] == "dev"}
    del cands_all
    pair_all = json.loads(lm.PAIR_INDEX_FILE.read_text(encoding="utf-8"))
    pair_idx = {d: v for d, v in pair_all.items() if d in dev_set}
    del pair_all
    items = lm._load_split_items("dev")
    ann = lm._annotations("dev", items)
    drugs = v1.load_universe(ROOT / "data" / "universe.json")
    names = sorted(v1.universe_index(drugs))
    return {"dev": dev, "rank": rank, "dis": dis, "cands": cands, "pair_idx": pair_idx, "items": items, "ann": ann, "names": names}


def item_table(items: Sequence[dict], ann) -> Tuple[Dict[Tuple[str, str, str, str], Dict[str, bool]], Dict[str, Dict[str, Any]]]:
    """(disease, drug, pmid, sentence) -> label flags; id -> per-item record with both annotators' labels and the consensus direction."""
    a, b = ann[lm.ANNOTATORS[0]], ann[lm.ANNOTATORS[1]]
    table, rec = {}, {}
    for it in items:
        ra, rb = a[it["id"]], b[it["id"]]
        c = lm.consensus(ra["relevance"], rb["relevance"], ra["direction"], rb["direction"], it["sentence"], it["drug"])
        table[(it["disease"], it["drug"], it["pmid"], it["sentence"])] = label_flags(ra["relevance"], rb["relevance"], c["direction"])
        rec[it["id"]] = {"disease": it["disease"], "drug": it["drug"], "sentence": it["sentence"], "pmid": it["pmid"], "o_rel": ra["relevance"], "s_rel": rb["relevance"],
                         "o_dir": ra["direction"], "s_dir": rb["direction"], "cons_dir": c["direction"]}
    return table, rec


def pair_gold(cands) -> Dict[Tuple[str, str], Dict[str, bool]]:
    return {(d, n): f for d, v in cands.items() for n, f in v["candidates"].items()}


# ============================================================================= Q1
def run_q1(ctx: Dict[str, Any], dts: Dict[str, DiseaseText], table) -> Dict[str, Any]:
    gold = pair_gold(ctx["cands"])
    keys = STAGES + SIDE
    per_frac: Dict[int, List[Dict[str, Dict[str, int]]]] = {}
    for fi, rep in dl.sample_keys():
        per_frac.setdefault(fi, []).append({g: {k: 0 for k in keys} for g in ("gold_positive", "other")})
    for d in ctx["dev"]:
        dt = dts[d]
        n = len(dt.pids)
        cand_names = sorted(ctx["cands"][d]["candidates"])
        for fi, rep in dl.sample_keys():
            kept = set(dl.subsample_indices(n, dl.FRACTIONS[fi], ctx["rank"][d], fi, rep).tolist())
            any_names: Set[str] = set()
            gene_names: Set[str] = set()
            for i in kept:
                any_names |= dt.names_any[i]
                if dt.gene[i]:
                    gene_names |= dt.names_any[i]
            acc = per_frac[fi][rep if fi else 0]
            for nm in cand_names:
                g = "gold_positive" if gold[(d, nm)]["g1"] else "other"
                fl = pair_stage_flags(ctx["pair_idx"][d]["pairs"].get(nm, []), kept, nm in any_names, nm in gene_names, d, nm, table)
                for k in keys:
                    acc[g][k] += bool(fl[k])
    out: Dict[str, Any] = {"definitions": {"gold_positive": "candidate pair with g1 (gold 1) true in the DEV candidate file", "other": "all other DEV candidate pairs"},
                           "fractions": {}}
    n_gp = sum(1 for f in gold.values() if f["g1"])
    n_all = len(gold)
    out["n_candidate_pairs"] = n_all
    out["n_gold_positive_pairs"] = n_gp
    out["n_other_pairs"] = n_all - n_gp
    out["n_other_with_gold2_or_gold3"] = sum(1 for f in gold.values() if not f["g1"] and (f["g2"] or f["g3"]))
    for fi, f in enumerate(dl.FRACTIONS):
        reps = per_frac[fi]
        fr: Dict[str, Any] = {}
        for g in ("gold_positive", "other"):
            tot = n_gp if g == "gold_positive" else n_all - n_gp
            stg = {}
            prev = None
            for k in keys:
                v = [r[g][k] for r in reps]
                m = float(np.mean(v))
                stg[k] = {"mean": m, "min": int(min(v)), "max": int(max(v)), "share_of_S0": m / tot if tot else None,
                          "retention_vs_previous_stage": (m / prev) if (prev and k in STAGES[1:]) else None}
                if k in STAGES:
                    prev = m
            fr[g] = {"n_pairs": tot, "stages": stg, "loss": largest_losses({s: stg[s]["mean"] for s in STAGES})}
        out["fractions"][str(f)] = fr
    out["n_replicates"] = {"1.0": 1, "other_fractions": dl.N_REPS}
    return out


# ============================================================================= Q2 (a): design facts
def design_facts(ctx, dts, example_disease: str) -> Dict[str, Any]:
    man = json.loads((lm.INPUT_DIR / "manifest.json").read_text(encoding="utf-8"))["splits"]["dev"]
    n_ret = [len(dts[d].pids) for d in ctx["dev"]]
    gene_share = [float(np.mean(dts[d].gene)) if dts[d].gene else 0.0 for d in ctx["dev"]]
    any_cand = []
    for d in ctx["dev"]:
        cand = set(ctx["cands"][d]["candidates"])
        any_cand.append(float(np.mean([bool(dts[d].names_any[i] & cand) for i in range(len(dts[d].pids))])) if dts[d].pids else 0.0)
    kept_mean = {str(f): float(np.mean([dl.n_sub(len(dts[d].pids), f) for d in ctx["dev"]])) for f in dl.FRACTIONS}
    return {
        "paper_retrieval": {
            "function": "engine_score.fetch_papers / build_gene_query (called by engine_v4_llm._pair_job, cache only)",
            "query_form": es.build_gene_query(ctx["dis"][example_disease]["symbols"][:3], es.YEAR_PRIMARY) + "   (example: first 3 gene symbols of one DEV disease)",
            "query_uses": "HGNC gene symbols of the disease only, TITLE_ABS phrase, OR-pooled; year <= 2024; abstract present; NO disease term, NO drug term, NO language filter",
            "record_cap": es.MAX_PAPERS, "ordering": "Europe PMC default order (no sort parameter), cursorMark paging, first 300 records with an abstract",
            "text_used": "title + '. ' + abstract (no full text)", "n_retrieved_per_dev_disease": describe(n_ret),
            "diseases_at_cap": int(sum(1 for x in n_ret if x >= es.MAX_PAPERS)),
            "share_of_retrieved_papers_with_a_gene_symbol_in_text": describe(gene_share),
            "share_of_retrieved_papers_naming_any_candidate_drug": describe(any_cand),
            "mean_papers_kept_by_fraction": kept_mean},
        "sentence_selection": {
            "function": "engine_v4_llm.paper_drug_sentences / pair_papers / union_pair_papers",
            "drug_rule": "elc.drug_regex over all 2705 universe names: exact whole-word, case-insensitive match of the single project name (ChEMBL preferred name); names under 4 characters are never matched; "
                         "longest alternative wins; no synonyms, brand names, salts, plurals or abbreviations; title+abstract sentences only",
            "paper_gate": "paper must contain a disease gene symbol (v1 gate)", "name_gate": "sentence must contain the first word of some drug name (speed only)",
            "caps": {"papers_per_pair": K_PAPERS, "sentences_per_paper": K_SENT}, "pair_universe": "top-30 B0 drugs union gold of the disease (candidate set P_d)",
            "verification": "sentence must be a verbatim substring of the Europe PMC title+abstract after one join-period repair; 67 DEV items dropped",
            "relevance_rule_for_features": "both annotators yes; direction 2-of-3 of Opus, Sonnet, rule"},
        "manifest_dev": {k: man[k] for k in ("n_candidate_pairs", "n_pairs_with_sentences", "n_items_before_verification", "n_items", "n_papers")}}


# ============================================================================= Q2 (b): relaxation, direct query, recall, language
def q2_relaxation(ctx, dts, sample, syn, dis_terms, getter=None) -> Dict[str, Any]:
    per_variant: Dict[str, List[Dict[str, int]]] = {v: [] for v in VARIANTS}
    per_variant_f03: Dict[str, List[Dict[str, int]]] = {v: [] for v in VARIANTS[:6]}
    r0_check = {"pairs": 0, "same_paper_list": 0, "same_sentence_count": 0}
    recall_rows, hit_rows, lang_rows = [], [], []
    r6_pid_overlap = []
    syn_cov = {"pairs_with_synonym": 0, "by_source": Counter()}
    for (d, nm) in sample:
        dt = dts[d]
        n = len(dt.pids)
        allp = range(n)
        terms = sorted(syn.get(nm, {}))
        if terms:
            syn_cov["pairs_with_synonym"] += 1
            for t in terms:
                for s in syn[nm][t]:
                    syn_cov["by_source"][s] += 1
        mp, ms = match_prod(dt, nm), match_syn(dt, nm, terms)
        r = {
            "R0_baseline": extract_pair(dt, allp, mp, K_PAPERS, K_SENT, True),
            "R1_sentence_cap_lifted": extract_pair(dt, allp, mp, K_PAPERS, None, True),
            "R2_paper_cap_lifted": extract_pair(dt, allp, mp, None, K_SENT, True),
            "R3_both_caps_lifted": extract_pair(dt, allp, mp, None, None, True),
            "R4_R3_no_gene_gate": extract_pair(dt, allp, mp, None, None, False),
            "R5_R3_drug_synonyms": extract_pair(dt, allp, ms, None, None, True)}
        # R0 reproduces the stored pair index?
        stored = ctx["pair_idx"][d]["pairs"].get(nm, [])[:K_PAPERS]
        r0_check["pairs"] += 1
        r0_check["same_paper_list"] += [i for i, _ in r["R0_baseline"]] == [e["i"] for e in stored]
        r0_check["same_sentence_count"] += sum(len(j) for _, j in r["R0_baseline"]) == sum(len(e["sents"]) for e in stored)
        recs = {k: records_of(dt, v) for k, v in r.items()}
        # direct queries
        dn_re, ds_re = elc.term_regex(dis_terms[d][0]), elc.term_regex(dis_terms[d][1])
        qd = query_disease_terms(d, ctx["dis"][d]["name_set"])
        pa = epmc_search(direct_query([nm], qd), PAGE_FIRST, "*", getter)
        pb = epmc_search(direct_query([nm] + terms, qd), PAGE_FIRST, "*", getter) if terms else pa
        genes_clause = clause(ctx["dis"][d]["symbols"])
        pg = epmc_search(direct_query([nm], qd, extra=f" AND {genes_clause}"), 1, "*", getter)
        pe = epmc_search(direct_query([nm], qd, extra=' AND LANG:"eng"'), 1, "*", getter)
        pn = epmc_search(direct_query([nm], qd, extra=' AND NOT LANG:"eng"'), 1, "*", getter)
        retrieved = set(dt.pids)
        first_a = [x["pid"] for x in pa["records"]]
        first_b = [x["pid"] for x in pb["records"]]
        recall_rows.append({"disease": d, "drug": nm, "hits": pa["hitCount"], "hits_with_gene_clause": pg["hitCount"], "n_first": len(first_a),
                            "recall_first20_INN": recall(first_a, retrieved), "recall_first20_synonyms": recall(first_b, retrieved),
                            "recall_first10_INN": recall(first_a[:R6_PAPERS], retrieved), "n_retrieved": n})
        hit_rows.append({"disease": d, "drug": nm, "hits": pa["hitCount"], "hits_eng": pe["hitCount"], "hits_not_eng": pn["hitCount"]})
        lang_rows.extend(x["language"] for x in pa["records"])
        mk = lambda page, rx: [{"pid": x["pid"], "text": record_text(x)} for x in page["records"][:R6_PAPERS]]
        inn_rx = elc.term_regex([nm])
        syn_rx = elc.term_regex(sorted({nm} | set(terms)))
        recs["R6a_direct_query_top10_INN"] = records_of_papers(mk(pa, inn_rx), inn_rx, dn_re, ds_re)
        recs["R6b_direct_query_top10_synonyms"] = records_of_papers(mk(pb, syn_rx), syn_rx, dn_re, ds_re)
        recs["R7_union_R4_R5_R6b"] = dedupe_records(recs["R4_R3_no_gene_gate"] + recs["R5_R3_drug_synonyms"] + recs["R6b_direct_query_top10_synonyms"])
        r6_pid_overlap.append(recall(first_a[:R6_PAPERS], retrieved))
        for v in VARIANTS:
            per_variant[v].append(pair_metrics(recs[v], nm))
        # f = 0.03, 5 replicates (mean per pair)
        fi = dl.FRACTIONS.index(0.03)
        reps = []
        for rep in range(dl.N_REPS):
            kept = dl.subsample_indices(n, 0.03, ctx["rank"][d], fi, rep).tolist()
            rr = {"R0_baseline": extract_pair(dt, kept, mp, K_PAPERS, K_SENT, True), "R1_sentence_cap_lifted": extract_pair(dt, kept, mp, K_PAPERS, None, True),
                  "R2_paper_cap_lifted": extract_pair(dt, kept, mp, None, K_SENT, True), "R3_both_caps_lifted": extract_pair(dt, kept, mp, None, None, True),
                  "R4_R3_no_gene_gate": extract_pair(dt, kept, mp, None, None, False), "R5_R3_drug_synonyms": extract_pair(dt, kept, ms, None, None, True)}
            reps.append({k: pair_metrics(records_of(dt, v), nm) for k, v in rr.items()})
        for k in VARIANTS[:6]:
            per_variant_f03[k].append({m: float(np.mean([x[k][m] for x in reps])) for m in reps[0][k]})
    summ = {v: summarise_variant(per_variant[v]) for v in VARIANTS}
    base = summ["R0_baseline"]
    for v in VARIANTS[1:]:
        s = summ[v]
        b_mean, v_mean = base["sent_drug_disease_set"]["mean"], s["sent_drug_disease_set"]["mean"]
        s["gain_vs_R0"] = {"mean_sent_drug_disease_set_ratio": (v_mean / b_mean) if b_mean else None,
                           "mean_sent_drug_disease_set_difference": v_mean - b_mean,
                           "share_pairs_ge1_difference_points": 100 * (s["sent_drug_disease_set"]["share_pairs_with_ge1"] - base["sent_drug_disease_set"]["share_pairs_with_ge1"]),
                           "major_artefact_rule": bool((b_mean and v_mean / b_mean >= 2) or
                                                       (100 * (s["sent_drug_disease_set"]["share_pairs_with_ge1"] - base["sent_drug_disease_set"]["share_pairs_with_ge1"]) >= 20))}
    summ_f03 = {v: summarise_variant(per_variant_f03[v], thr=1e-9) for v in per_variant_f03}
    rec_first = [x["recall_first20_INN"] for x in recall_rows if x["recall_first20_INN"] is not None]
    rec_syn = [x["recall_first20_synonyms"] for x in recall_rows if x["recall_first20_synonyms"] is not None]
    hits = [x["hits"] for x in recall_rows]
    gene_hits = [x["hits_with_gene_clause"] for x in recall_rows]
    elig = [x["hits_with_gene_clause"] / x["hits"] for x in recall_rows if x["hits"]]
    recall_out = {
        "pairs": len(recall_rows), "pairs_with_zero_direct_hits": int(sum(1 for h in hits if h == 0)), "hit_count_distribution": describe(hits),
        "recall_of_first20_records_in_the_300_retrieved_INN": describe(rec_first), "recall_first20_synonym_query": describe(rec_syn),
        "recall_first10_INN": describe([x["recall_first10_INN"] for x in recall_rows if x["recall_first10_INN"] is not None]),
        "pairs_with_recall_zero": int(sum(1 for x in rec_first if x == 0.0)), "pairs_with_recall_one": int(sum(1 for x in rec_first if x == 1.0)),
        "share_of_direct_hits_that_also_match_the_gene_clause": describe(elig),
        "hits_total": int(sum(hits)), "hits_with_gene_clause_total": int(sum(gene_hits)),
        "reading": "recall = share of the first 20 direct-query records that are among the 300 retrieved; the gene-clause share bounds what ANY gene-based retrieval can reach",
        "rows": recall_rows}
    hits_tot = sum(x["hits"] for x in hit_rows)
    lang_out = {"direct_query_hits_total": int(hits_tot), "direct_query_hits_english": int(sum(x["hits_eng"] for x in hit_rows)),
                "direct_query_hits_not_english": int(sum(x["hits_not_eng"] for x in hit_rows)),
                "direct_query_not_english_share_of_hits": (sum(x["hits_not_eng"] for x in hit_rows) / hits_tot) if hits_tot else None,
                "first20_records_language": language_summary(lang_rows)}
    return {"variants": summ, "variants_at_f0.03_mean_of_5_replicates": summ_f03, "R0_reproduces_stored_pair_index": r0_check, "synonym_coverage": {
        "pairs_with_at_least_one_synonym": syn_cov["pairs_with_synonym"], "n_pairs": len(sample), "synonym_terms_by_source": dict(syn_cov["by_source"])},
        "recall": recall_out, "language_direct_query": lang_out}


def q2_language_retrieved(ctx, dts, getter=None) -> Dict[str, Any]:
    gold_dis = sorted({d for d, v in ctx["cands"].items() if any(f["g1"] for f in v["candidates"].values())})
    pick = sorted(random.Random(SEED).sample(gold_dis, min(N_LANG_DISEASES, len(gold_dis))))
    rows, langs_all = [], []
    for d in pick:
        sy = ctx["dis"][d]["symbols"]
        g = fetch_gene_records(sy, es.MAX_PAPERS, getter)
        base = es.build_gene_query(sy, es.YEAR_PRIMARY)
        he = epmc_search(base + ' AND LANG:"eng"', 1, "*", getter)["hitCount"]
        hn = epmc_search(base + ' AND NOT LANG:"eng"', 1, "*", getter)["hitCount"]
        ls = language_summary([r["language"] for r in g["records"]])
        langs_all.extend(r["language"] for r in g["records"])
        cached = dts[d].pids
        fresh = [r["pid"] for r in g["records"]]
        rows.append({"disease": d, "n_fresh": len(fresh), "n_cached": len(cached), "same_ids_as_cache_share": (len(set(fresh) & set(cached)) / len(cached)) if cached else None,
                     "gene_query_hits": g["hitCount"], "hits_english": he, "hits_not_english": hn, "non_english_share_of_retrieved": ls["non_english_share"],
                     "retrieved_languages": ls["counts"]})
    return {"diseases": pick, "pooled_retrieved_records_language": language_summary(langs_all),
            "non_english_share_of_retrieved_by_disease": describe([r["non_english_share_of_retrieved"] for r in rows if r["non_english_share_of_retrieved"] is not None]),
            "gene_query_hits_not_english_share": describe([r["hits_not_english"] / r["gene_query_hits"] for r in rows if r["gene_query_hits"]]),
            "rows": rows}


def q2_cap(ctx, dts, getter=None) -> Dict[str, Any]:
    hits = {}
    for d in ctx["dev"]:
        hits[d] = epmc_search(es.build_gene_query(ctx["dis"][d]["symbols"], es.YEAR_PRIMARY), 1, "*", getter)["hitCount"]
    ratio = [len(dts[d].pids) / hits[d] for d in ctx["dev"] if hits[d]]
    return {"diseases": len(hits), "gene_query_hits": describe(list(hits.values())), "diseases_with_more_hits_than_300": int(sum(1 for h in hits.values() if h > es.MAX_PAPERS)),
            "share_of_gene_query_hits_retrieved": describe(ratio)}


def q2_naming(ctx, dts, syn, gold_pairs) -> Dict[str, Any]:
    """Over all gold-positive DEV pairs: plain-INN share, synonym availability, and sentences that name the drug only by a synonym."""
    drugs = sorted({nm for _, nm in gold_pairs})
    plain = {nm: is_plain_inn(nm) for nm in drugs}
    has_syn = {nm: bool(syn.get(nm)) for nm in drugs}
    prod_n = syn_only_n = 0
    drugs_syn_only = set()
    for d, nm in gold_pairs:
        terms = sorted(syn.get(nm, {}))
        if not terms:
            continue
        dt = dts[d]
        rx = elc.term_regex(terms)
        for i in range(len(dt.pids)):
            for j, s in enumerate(dt.sents[i]):
                p = nm in dt.prod[i][j]
                q = bool(rx.search(s))
                prod_n += p
                if q and not p:
                    syn_only_n += 1
                    drugs_syn_only.add(nm)
    by_src = Counter(s for nm in drugs for t in syn.get(nm, {}) for s in syn[nm][t])
    return {"gold_positive_pairs": len(gold_pairs), "distinct_gold_positive_drugs": len(drugs),
            "drugs_with_non_plain_project_name": int(sum(1 for v in plain.values() if not v)),
            "share_drugs_with_non_plain_project_name": float(np.mean([not v for v in plain.values()])),
            "drugs_with_at_least_one_project_synonym": int(sum(has_syn.values())), "share_drugs_with_synonym": float(np.mean(list(has_syn.values()))),
            "synonym_terms_by_source_over_gold_drugs": dict(by_src),
            "retrieved_300_sentences_naming_drug_by_project_name": int(prod_n), "retrieved_300_sentences_naming_drug_only_by_synonym": int(syn_only_n),
            "share_of_drug_sentences_only_by_synonym": (syn_only_n / (prod_n + syn_only_n)) if (prod_n + syn_only_n) else None,
            "drugs_with_at_least_one_synonym_only_sentence": len(drugs_syn_only),
            "share_of_gold_drugs_with_synonym_only_sentence": len(drugs_syn_only) / len(drugs) if drugs else None,
            "reading": "lower bound: only synonyms present in the project tables (ClinicalTrials.gov other names, partial ChEMBL records, salt stripping)"}


# ============================================================================= Q3
VOC_REL, VOC_DIR = ["yes", "no", "unclear"], ["benefit", "harm", "no_effect", "unclear"]

READING: Dict[str, Tuple[str, str]] = {
    'dev-02504': ('uninformative', 'drug only induces an osteonecrosis rat model'),
    'dev-04858': ('label_looks_right', 'treatments listed for a patient but no condition named in the sentence'),
    'dev-05568': ('uninformative', 'drug acts on osteosarcoma growth, not on the disease'),
    'dev-00682': ('uninformative', 'doxorubicin resistance in colorectal cancer, other disease'),
    'dev-07580': ('label_looks_right', 'donepezil is only the positive-control group of a herbal test'),
    'dev-00350': ('label_looks_right', 'drug is a treatment of the mice, but the sentence reports a miRNA screen'),
    'dev-03122': ('uninformative', 'apomorphine used as a tool to induce rotation'),
    'dev-00007': ('label_looks_right', 'cholesterol-lowering combination in a mouse, disease not stated'),
    'dev-05567': ('uninformative', 'patient with another movement disorder resistant to clonazepam'),
    'dev-05291': ('uninformative', 'levothyroxine for hypothyroidism, a different disease'),
    'dev-07273': ('uninformative', 'RET-fusion lung cancer trials, a different disease'),
    'dev-04223': ('label_looks_right', 'drug is a perturbation in a Western-blot assay on cell lines'),
    'dev-01790': ('uninformative', 'rat reproductive toxicity model'),
    'dev-00061': ('informative_label_too_strict', 'APL (an AML subtype) treated with mitoxantrone with remission reported; drug evaluated as leukaemia therapy'),
    'dev-01240': ('uninformative', 'liver cancer cells under glucose starvation'),
    'dev-06235': ('uninformative', 'title of a family report; vitamin E is the deficient nutrient, no treatment'),
    'dev-05360': ('uninformative', 'rheumatoid arthritis cohort, a different disease'),
    'dev-01879': ('label_looks_right', 'treatment effect stated but the condition is the abbreviation SAT, not the item disease'),
    'dev-05724': ('informative_label_too_strict', 'clobazam named among drugs that failed to control the seizures of a patient; a null treatment result in an epilepsy; Sonnet no looks strict'),
    'dev-04586': ('label_looks_right', 'drug list for an unspecified condition'),
    'dev-00961': ('uninformative', 'stallion testes gene expression'),
    'dev-03491': ('uninformative', 'title of a case report; vitamin E not given as treatment'),
    'dev-03596': ('uninformative', 'basic biology of vitamin E trafficking'),
    'dev-04830': ('uninformative', 'stallion testes sampling'),
    'dev-02889': ('uninformative', 'colorectal cancer cell stemness'),
    'dev-05386': ('uninformative', 'dexamethasone as the toxicant in rats'),
    'dev-03591': ('uninformative', 'cell-clock experiment, dexamethasone as synchroniser'),
    'dev-01686': ('label_looks_right', 'atorvastatin is the comparator in an assay of other compounds'),
    'dev-05371': ('label_looks_right', 'synergy with a statin in a cholesterol model, familial disease not tied'),
}   # id -> (class, one-line reason): my own reading of the `sample` output (audit opinion, not ground truth)
READING_CLASSES = ("uninformative", "informative_label_too_strict", "label_looks_right")


def q3a_pool(rec: Dict[str, Dict[str, Any]], gold) -> Tuple[List[str], Dict[Tuple[str, str], List[str]]]:
    by_pair: Dict[Tuple[str, str], List[str]] = {}
    for iid, r in rec.items():
        if gold.get((r["disease"], r["drug"]), {}).get("g1"):
            by_pair.setdefault((r["disease"], r["drug"]), []).append(iid)
    pool = sorted(p for p, ids in by_pair.items() if all(rec[i]["o_rel"] != "yes" and rec[i]["s_rel"] != "yes" for i in ids))
    return pool, by_pair


def q3a_sample(rec, gold, n: int = N_READ, seed: int = SEED) -> Dict[str, Any]:
    pool, by_pair = q3a_pool(rec, gold)
    rng = random.Random(seed)
    chosen = rng.sample(pool, min(n, len(pool)))
    ids = [rng.choice(sorted(by_pair[p])) for p in chosen]
    rest = sorted(i for p in pool for i in by_pair[p] if i not in set(ids))     # fewer than n pairs: fill with further sentences of the same pool pairs
    per_pair: Counter = Counter()
    for p in chosen:
        per_pair[p] += 1
    pair_of = {i: p for p in pool for i in by_pair[p]}
    rng.shuffle(rest)
    for i in rest:                                                                # at most 3 sentences per pair, for diversity
        if len(ids) >= n:
            break
        if per_pair[pair_of[i]] < 3:
            per_pair[pair_of[i]] += 1
            ids.append(i)
    return {"n_sentences_in_pool_pairs": sum(len(by_pair[p]) for p in pool), "n_one_per_pair": len(chosen), "n_gold_positive_pairs_with_items": len(by_pair), "n_pairs_no_yes_from_either_annotator": len(pool), "ids": ids}


def q3b(rec, gold) -> Dict[str, Any]:
    out = {}
    for name, sel in (("gold_positive_pairs", True), ("other_pairs", False)):
        ids = [i for i, r in rec.items() if bool(gold.get((r["disease"], r["drug"]), {}).get("g1")) == sel]
        rows = [rec[i] for i in ids]
        yes = {a: sum(r[k] == "yes" for r in rows) for a, k in (("opus", "o_rel"), ("sonnet", "s_rel"))}
        out[name] = {"n_items": len(rows), "relevance": crosstab([r["o_rel"] for r in rows], [r["s_rel"] for r in rows], VOC_REL),
                     "direction": crosstab([r["o_dir"] for r in rows], [r["s_dir"] for r in rows], VOC_DIR),
                     "share_yes": {a: (yes[a] / len(rows)) if rows else None for a in yes},
                     "share_yes_either": float(np.mean([r["o_rel"] == "yes" or r["s_rel"] == "yes" for r in rows])) if rows else None,
                     "share_yes_both": float(np.mean([r["o_rel"] == "yes" and r["s_rel"] == "yes" for r in rows])) if rows else None}
    return out


def read_pilot(root: Path = ROOT) -> Dict[str, Dict[str, Dict[str, Any]]]:
    d = lm.LABELS_DIR / "parts_v2"
    out: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for a in lm.ANNOTATORS:
        rows: Dict[str, Dict[str, Any]] = {}
        for sh in PILOT_SHARDS:
            doc = json.loads((d / f"{a}_dev_shard_{sh:03d}.json").read_text(encoding="utf-8"))
            rows.update({r["id"]: r for r in doc["labels"]})
        out[a] = rows
    return out


def q3d(rec, pilot) -> Dict[str, Any]:
    ids = sorted(pilot["opus"])
    assert ids == sorted(pilot["sonnet"]) and all(i in rec for i in ids), "pilot ids differ"
    res: Dict[str, Any] = {"n_items": len(ids)}
    v1r = {a: [rec[i]["o_rel" if a == "opus" else "s_rel"] for i in ids] for a in lm.ANNOTATORS}
    v2r = {a: [pilot[a][i]["relevance"] for i in ids] for a in lm.ANNOTATORS}
    v1d = {a: [rec[i]["o_dir" if a == "opus" else "s_dir"] for i in ids] for a in lm.ANNOTATORS}
    v2d = {a: [pilot[a][i]["direction"] for i in ids] for a in lm.ANNOTATORS}
    res["per_annotator_v1_vs_v2"] = {a: {"relevance": crosstab(v1r[a], v2r[a], VOC_REL), "direction": crosstab(v1d[a], v2d[a], VOC_DIR)} for a in lm.ANNOTATORS}
    res["opus_vs_sonnet_relevance"] = {"v1": crosstab(v1r["opus"], v1r["sonnet"], VOC_REL), "v2": crosstab(v2r["opus"], v2r["sonnet"], VOC_REL)}
    res["opus_vs_sonnet_direction"] = {"v1": crosstab(v1d["opus"], v1d["sonnet"], VOC_DIR), "v2": crosstab(v2d["opus"], v2d["sonnet"], VOC_DIR)}
    dis1 = [a != b for a, b in zip(v1r["opus"], v1r["sonnet"])]
    dis2 = [a != b for a, b in zip(v2r["opus"], v2r["sonnet"])]
    res["relevance_disagreement"] = {"v1_items": int(sum(dis1)), "v2_items": int(sum(dis2)),
                                     "removed_by_v2": int(sum(a and not b for a, b in zip(dis1, dis2))), "kept_in_v2": int(sum(a and b for a, b in zip(dis1, dis2))),
                                     "newly_created_by_v2": int(sum(b and not a for a, b in zip(dis1, dis2))),
                                     "share_of_v1_disagreement_removed_net": (1 - sum(dis2) / sum(dis1)) if sum(dis1) else None}
    cnt = lambda o, s: {"yes_either": int(sum(a == "yes" or b == "yes" for a, b in zip(o, s))), "yes_both": int(sum(a == "yes" and b == "yes" for a, b in zip(o, s))),
                        "opus_yes": int(sum(a == "yes" for a in o)), "sonnet_yes": int(sum(b == "yes" for b in s))}
    res["yes_counts"] = {"v1": cnt(v1r["opus"], v1r["sonnet"]), "v2": cnt(v2r["opus"], v2r["sonnet"])}
    return res


# ============================================================================= verdict
VERDICT_PROSE: List[str] = [   # hand-written interpretive paragraphs, written after the numbers above were read; every number is in the tables
    "**Short answer.** The claim is right about the outcome and wrong about the cause. With the current construction the LLM features are almost empty at 3 % of the literature "
    "(gold-positive pairs with a counted relevant benefit sentence: 3.4 of 465, 0.7 %), but they are also nearly empty at 100 % of the literature (37 of 465, 8.0 %; 74 of 465, 15.9 % if either annotator may say yes). "
    "So the problem is not 'low volume'; thinning multiplies an already poor base. The base is poor mainly because of how papers are retrieved (ii), partly because some gold pairs "
    "really have little drug-and-disease literature (i), and only marginally because of label strictness (iii). Language and drug naming are not material losses.",
    "**Where the pairs are lost (Q1, gold-positive, f = 1.0).** 465 candidate pairs -> 163 (35.1 %) have the drug named in any of the retrieved papers -> 159 have an extracted and verified sentence -> 138 have `yes` from "
    "either annotator -> 77 from both -> 37 with benefit direction. The largest single loss is the first step (302 pairs, retention 35 %); the extraction, the paper gate and the verification lose 4 pairs in total. "
    "At f = 0.03 the first step keeps 3.7 % of the pairs (17.4 of 465). Non-gold candidates behave the same way at the first step (57 % keep a mention at f = 1.0) but lose most pairs later, at the relevance step (545 of 2258 with a verified sentence have `yes` by either).",
    "**(ii) Retrieval is the artefact.** The paper query pools the gene symbols of the disease only, with no disease term, no drug term and no language filter, and keeps the first 300 records with an abstract; 89 of 103 DEV "
    "diseases are at that cap and the 300 are a median 8.7 % of the gene-query hits. A direct `drug AND disease` query shows how far this is from the drug-and-disease literature: of the first 20 direct records "
    "only 3.8 % (mean over the 93 pairs that have hits) are among the 300 retrieved, and only 11 % of all direct hits (3352 of 29871) even contain a disease gene symbol, so 89 % cannot be reached by any gene-symbol query whatever the cap. "
    "On the 150 seeded gold-positive pairs the baseline gives 0.11 sentences per pair naming drug and disease (11 pairs, 7 %); the top-10 direct-query papers give 3.9 (80 pairs, 53 %), and 0.83 vs 0.03 sentences per pair that the rule extractor reads as evaluative (benefit or harm). "
    "Rules relaxed one at a time inside the existing retrieval change little: lifting the sentence cap adds nothing to drug-and-disease sentences; lifting the paper cap gives 0.25 vs 0.11 (16 vs 11 pairs; flagged 'major' by the doubling convention only because the base is tiny); "
    "removing the gene gate adds nothing; drug synonyms add nothing to drug-and-disease sentences (0.25).",
    "**(i) Part of it is the literature.** 57 of the 150 gold-positive pairs (38 %) have zero direct hits for the drug AND the disease name or up to five of its synonyms in title or abstract (median 6.5 hits per pair). Some of these are rare diseases whose MONDO labels are not the words used in papers, "
    "so this is an upper bound on true absence, but for them no retrieval design will provide sentences; for those pairs the claim 'not testable' is correct. The 93 pairs with hits have a large literature (mean 199 hits over all 150 pairs, maximum 4177).",
    "**(iii) Label strictness is minor for coverage.** On DEV gold-positive items Sonnet says `yes` for 72.7 % and Opus for 33.0 % (kappa 0.35, versus 0.55 for other pairs); `either` vs `both` changes the gold-positive pairs with a counted `yes` sentence from 77 to 138 at f = 1.0 (x1.8) and from 8.4 to 13.8 at f = 0.03 (x1.6), "
    "below the doubling convention. In my reading of the gold-positive pairs where neither annotator said `yes` (only 21 pairs; 29 sentences read, audit opinion) 18 sentences were uninformative, 9 labels looked right and 2 looked too strict (7 %), below the quarter threshold. "
    "The disagreement is mostly Opus `unclear` against Sonnet `yes` (311 items on gold-positive pairs). The v2 guideline removes 102 of 125 v1 disagreements on the 800 pilot items (kappa 0.62 -> 0.85) but creates 28 and does so by making Sonnet stricter (`yes` by either 153 -> 70, by both 57 -> 41), so it lowers coverage further. "
    "The much higher `yes` rate on gold-positive than on other pairs (either 73 % vs 17 %) is compatible with real signal in the sentences and with annotator knowledge of indications (parent S4); this audit cannot separate the two.",
    "**Language and naming.** Explicit non-English records are 1.4 % of the 8423 retrieved records of 30 seeded diseases and 3.7 % of the first 20 direct records; 11 % of retrieved records have no language field and all are preprints, so they are not counted as non-English. "
    "Below the 5 % convention, not material. Drug naming: 2 of 202 gold-positive drugs have a non-plain project name; 104 have a project synonym (mostly ClinicalTrials.gov other names, so a lower bound); 16 (7.9 %) have sentences naming the drug only by a synonym, "
    "19.5 % of the drug sentences of pairs with synonyms, but they add nothing to drug-and-disease sentences, so naming is not the bottleneck either.",
    "**Recommendation (not implemented, nothing annotated).** The cheapest change with the largest gain is to retrieve per candidate pair instead of per disease: one Europe PMC query `drug (INN and synonyms) AND disease (name and a few synonyms)` for the first 10 to 20 records, and keep the sentences that name the drug "
    "(prefer those that also name the disease). Estimated gain from the 150-pair sample: pairs with at least one drug-and-disease sentence from 7 % to 53 %, sentences per pair from 0.11 to 3.9, i.e. about 7 times the pairs and 35 times the sentences. "
    "Cost: about one call per candidate pair (4442 DEV pairs, at most 4 requests per second), and about 4 more sentences per pair to annotate (about 17 thousand DEV items if the other pairs behave like gold-positive ones; fewer in practice). "
    "Low volume can then be simulated by keeping a fraction of each pair's own hits, which keeps the test meaningful for pairs with many hits; pairs with none stay empty, which is the true literature limit. "
    "Estimates come from gold-positive pairs only and count sentences, not relevant sentences; a second, much smaller lever is accepting `yes` by either annotator (x1.6 to x1.8 pairs), which does not address the 65 % of gold pairs never retrieved.",
]


def verdict_skeleton(res: Dict[str, Any]) -> Dict[str, Any]:
    q1 = res["Q1_funnel"]["fractions"]
    gp1, gp3 = q1["1.0"]["gold_positive"], q1["0.03"]["gold_positive"]
    v = res["Q2_relaxation"]["variants"]
    major = [k for k in VARIANTS[1:] if v[k]["gain_vs_R0"]["major_artefact_rule"]]
    lang = res["Q2_language_retrieved"]["pooled_retrieved_records_language"]
    first20 = res["Q2_relaxation"]["language_direct_query"]["first20_records_language"]
    explicit = lambda ls: (ls["non_english"] - ls["counts"].get("unknown", 0)) / ls["n"] if ls["n"] else None
    gps = lambda fr, k: q1[fr]["gold_positive"]["stages"][k]["share_of_S0"]
    r = res["Q2_relaxation"]["recall"]
    labs = res["Q3_labels"]
    st1 = lambda k, fr: q1[fr]["gold_positive"]["stages"][k]["mean"]
    cls = labs["a_reading"]["class_counts"]
    n_read = sum(cls.values())
    ratio = lambda fr: (st1("S4_relevance_yes_either", fr) / st1("S5_relevance_yes_both", fr)) if st1("S5_relevance_yes_both", fr) else None
    return {"gold_positive_at_f1.0": {"largest_absolute_loss": gp1["loss"]["largest_absolute_loss"], "lowest_retention": gp1["loss"]["lowest_retention"]},
            "gold_positive_at_f0.03": {"largest_absolute_loss": gp3["loss"]["largest_absolute_loss"], "lowest_retention": gp3["loss"]["lowest_retention"]},
            "gold_positive_share_of_S0": {fr: {k: gps(fr, k) for k in STAGES[1:] + ["S6b_yes_either_and_benefit"]} for fr in ("1.0", "0.1", "0.03")},
            "variants_meeting_major_artefact_rule": major,
            "relaxation_pairs_with_ge1_drug_disease_sentence": {k: v[k]["sent_drug_disease_set"]["pairs_with_ge1"] for k in VARIANTS},
            "relaxation_mean_drug_disease_sentences": {k: v[k]["sent_drug_disease_set"]["mean"] for k in VARIANTS},
            "recall_mean_first20": r["recall_of_first20_records_in_the_300_retrieved_INN"].get("mean"), "pairs_zero_direct_hits": r["pairs_with_zero_direct_hits"],
            "direct_hits_also_matching_gene_clause_share_of_total": (r["hits_with_gene_clause_total"] / r["hits_total"]) if r["hits_total"] else None,
            "language": {"explicit_non_english_share_retrieved": explicit(lang), "unknown_language_share_retrieved": lang["counts"].get("unknown", 0) / lang["n"],
                         "explicit_non_english_share_first20_direct": explicit(first20), "material_at_5_percent": bool((explicit(lang) or 0) > 0.05 or (explicit(first20) or 0) > 0.05),
                         "note": "records with no language field are all preprints (source PPR) in the retrieved set and are not counted as non-English"},
            "label_strictness": {"either_over_both_pairs_gold_positive_f1.0": ratio("1.0"), "either_over_both_pairs_gold_positive_f0.03": ratio("0.03"),
                                 "reading_too_strict_share": (cls.get("informative_label_too_strict", 0) / n_read) if n_read else None, "n_read": n_read}}


# ============================================================================= run
def _py(x):
    if isinstance(x, dict):
        return {str(k): _py(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_py(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def build_dts(ctx) -> Dict[str, DiseaseText]:
    names = ctx["names"]
    dre = elc.drug_regex(names)
    gate = es.name_gate(names)
    dts = {}
    t0 = time.time()
    for k, d in enumerate(ctx["dev"]):
        papers = es.fetch_papers(ctx["dis"][d]["symbols"], es.YEAR_PRIMARY, v4.CACHE / "epmc", dl.cache_only_getter)
        only, full = disease_terms(d, ctx["dis"][d]["name_set"])
        dts[d] = DiseaseText(papers, ctx["dis"][d]["symbols"], dre, gate, elc.term_regex(only or [d + "#"]), elc.term_regex(full or [d + "#"]))
        if (k + 1) % 20 == 0:
            print(json.dumps({"disease_text": k + 1, "of": len(ctx["dev"]), "seconds": round(time.time() - t0)}), flush=True)
    return dts


def cmd_sample() -> None:
    ctx = load_context()
    table, rec = item_table(ctx["items"], ctx["ann"])
    s = q3a_sample(rec, pair_gold(ctx["cands"]))
    print(json.dumps({k: v for k, v in s.items() if k != "ids"}))
    for k, i in enumerate(s["ids"], 1):
        r = rec[i]
        print(f"[{k}] {i} | {r['drug']} | {r['disease']} | O:{r['o_rel']}/{r['o_dir']} S:{r['s_rel']}/{r['s_dir']}\n    {r['sentence']}")


def run(getter=None) -> Dict[str, Any]:
    t0 = time.time()
    ctx = load_context()
    gold = pair_gold(ctx["cands"])
    table, rec = item_table(ctx["items"], ctx["ann"])
    dts = build_dts(ctx)
    res: Dict[str, Any] = {"label": LABEL, "protocol": str(PROTOCOL.relative_to(ROOT)), "protocol_sha256": v4.sha256_file(PROTOCOL), "seed": SEED,
                           "heldout_touched": False, "lock_file_created": False, "unvalidated_registry": UNVALIDATED,
                           "n_dev_diseases": len(ctx["dev"]), "n_items_dev": len(ctx["items"])}
    res["Q1_funnel"] = run_q1(ctx, dts, table)
    print(json.dumps({"q1_done": True, "seconds": round(time.time() - t0)}), flush=True)
    res["Q2_design"] = design_facts(ctx, dts, ctx["dev"][0])
    gold_pairs = sorted(p for p, f in gold.items() if f["g1"])
    sample = sample_pairs(gold_pairs, N_SAMPLE_PAIRS)
    syn = build_synonyms(ctx["names"], sorted((v4.CACHE / "ctgov").glob("*.json")), sorted((ROOT / "data" / "cache").glob("*.json")))
    dis_terms = {d: disease_terms(d, ctx["dis"][d]["name_set"]) for d in ctx["dev"]}
    res["Q2_relaxation"] = q2_relaxation(ctx, dts, sample, syn, dis_terms, getter)
    print(json.dumps({"q2_relaxation_done": True, "seconds": round(time.time() - t0)}), flush=True)
    res["Q2_language_retrieved"] = q2_language_retrieved(ctx, dts, getter)
    res["Q2_cap"] = q2_cap(ctx, dts, getter)
    res["Q2_drug_naming"] = q2_naming(ctx, dts, syn, gold_pairs)
    res["Q2_sample"] = {"n": len(sample), "seed": SEED, "pairs": [list(p) for p in sample]}
    # Q3
    s = q3a_sample(rec, gold)
    rows = []
    for i in s["ids"]:
        r = rec[i]
        cls, why = READING.get(i, (None, None))
        rows.append({"id": i, "drug": r["drug"], "disease": r["disease"], "sentence": r["sentence"], "opus": [r["o_rel"], r["o_dir"]], "sonnet": [r["s_rel"], r["s_dir"]],
                     "class": cls, "reason": why})
    cls_counts = Counter(r["class"] or "not_yet_read" for r in rows)
    ev = {}
    for fr in ("1.0", "0.1", "0.03"):
        ev[fr] = {g: {k: res["Q1_funnel"]["fractions"][fr][g]["stages"][k]["mean"] for k in ("S3_verified_sentence", "S4_relevance_yes_either", "S5_relevance_yes_both",
                                                                                              "S6_yes_both_and_benefit", "S6b_yes_either_and_benefit", "opus_yes", "sonnet_yes",
                                                                                              "opus_yes_benefit", "sonnet_yes_benefit", "benefit_direction_only")}
                  for g in ("gold_positive", "other")}
    res["Q3_labels"] = {"a_reading": {"pool": {k: v for k, v in s.items() if k != "ids"}, "reader": "the assistant (audit opinion, not ground truth); the pairs are known to be gold-positive, "
                                      "which may bias the reading toward 'informative'", "class_counts": dict(cls_counts), "rows": rows},
                        "b_crosstabs": q3b(rec, gold), "c_sensitivity": ev, "d_v2_pilot": q3d(rec, read_pilot())}
    res["Q4_verdict"] = {"computed": verdict_skeleton(res), "prose": VERDICT_PROSE}
    res["seconds"] = round(time.time() - t0, 1)
    v4.write_json(OUT_JSON, _py(res))
    OUT_MD.write_text(render_md(_py(res)), encoding="utf-8", newline="\n")
    print(json.dumps({"audit_done": True, "seconds": res["seconds"]}), flush=True)
    return res


# ============================================================================= markdown
def _fmt(x, nd=1):
    return "-" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def _pct(x):
    return "-" if x is None else f"{100 * x:.1f}%"


def render_md(res: Dict[str, Any]) -> str:
    L: List[str] = []
    A = L.append
    A("# Audit of the low-coverage claim for LLM-read evidence (DEV only, read-only)\n")
    A(f"{res['label']}. Protocol: `{res['protocol']}`. Seed {res['seed']}. HELD-OUT touched: {res['heldout_touched']}.\n")
    q1 = res["Q1_funnel"]
    A("## Q1. Funnel of DEV candidate pairs\n")
    A(f"{q1['n_candidate_pairs']} candidate pairs: {q1['n_gold_positive_pairs']} gold-positive (gold 1), {q1['n_other_pairs']} other ({q1['n_other_with_gold2_or_gold3']} of them gold 2 / gold 3 only). "
      "Counts are means over 5 replicates for f < 1 (min-max in the JSON). Each cell: pairs (share of the group).\n")
    for f in ("1.0", "0.3", "0.1", "0.03", "0.01"):
        A(f"### f = {f}\n")
        A("| stage | gold-positive | retention | other | retention |")
        A("|---|---|---|---|---|")
        for k in STAGES + SIDE:
            g, o = q1["fractions"][f]["gold_positive"]["stages"][k], q1["fractions"][f]["other"]["stages"][k]
            A(f"| {k} | {_fmt(g['mean'])} ({_pct(g['share_of_S0'])}) | {_pct(g['retention_vs_previous_stage'])} | {_fmt(o['mean'])} ({_pct(o['share_of_S0'])}) | {_pct(o['retention_vs_previous_stage'])} |")
        for grp in ("gold_positive", "other"):
            ls = q1["fractions"][f][grp]["loss"]
            A(f"\n{grp}: largest absolute loss {ls['largest_absolute_loss']['from']} -> {ls['largest_absolute_loss']['to']} ({_fmt(ls['largest_absolute_loss']['lost_pairs'])} pairs); "
              f"lowest retention {ls['lowest_retention']['from']} -> {ls['lowest_retention']['to']} ({_pct(ls['lowest_retention']['retention'])}).")
        A("")
    d = res["Q2_design"]
    A("## Q2. Retrieval and extraction design\n")
    for sec in ("paper_retrieval", "sentence_selection"):
        A(f"**{sec}**\n")
        for k, v in d[sec].items():
            A(f"- {k}: `{json.dumps(v, ensure_ascii=False)}`" if isinstance(v, (dict, list)) else f"- {k}: {v}")
        A("")
    A(f"Manifest (DEV): `{json.dumps(d['manifest_dev'])}`\n")
    q2 = res["Q2_relaxation"]
    A(f"### Relaxation experiment on {q2['synonym_coverage']['n_pairs']} seeded gold-positive pairs (no labelling)\n")
    A(f"R0 reproduces the stored pair index: paper list identical in {q2['R0_reproduces_stored_pair_index']['same_paper_list']}/{q2['R0_reproduces_stored_pair_index']['pairs']}, "
      f"sentence count identical in {q2['R0_reproduces_stored_pair_index']['same_sentence_count']}/{q2['R0_reproduces_stored_pair_index']['pairs']}. "
      f"Pairs with a project synonym beyond the INN: {q2['synonym_coverage']['pairs_with_at_least_one_synonym']} (terms by source: `{json.dumps(q2['synonym_coverage']['synonym_terms_by_source'])}`).\n")
    A("| variant | sentences naming drug (mean) | drug+disease sentence, name only (mean) | drug+disease sentence, name+synonyms (mean) | pairs with >=1 | papers drug+disease (mean) | evaluative (rule) drug+disease (mean) | ratio vs R0 | +points pairs>=1 | major |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for v in VARIANTS:
        s = q2["variants"][v]
        gn = s.get("gain_vs_R0", {})
        A(f"| {v} | {_fmt(s['sent_drug']['mean'], 2)} | {_fmt(s['sent_drug_disease_name']['mean'], 2)} | {_fmt(s['sent_drug_disease_set']['mean'], 2)} | "
          f"{s['sent_drug_disease_set']['pairs_with_ge1']} ({_pct(s['sent_drug_disease_set']['share_pairs_with_ge1'])}) | {_fmt(s['papers_drug_and_disease']['mean'], 2)} | "
          f"{_fmt(s['sent_drug_disease_evaluative_rule']['mean'], 2)} | {_fmt(gn.get('mean_sent_drug_disease_set_ratio'), 2)} | {_fmt(gn.get('share_pairs_ge1_difference_points'))} | {gn.get('major_artefact_rule', '-')} |")
    A("\nAt f = 0.03 (mean of 5 replicates, same pairs):\n")
    A("| variant | sentences naming drug | drug+disease (name+synonyms) | pairs with >=1 |")
    A("|---|---|---|---|")
    for v, s in q2["variants_at_f0.03_mean_of_5_replicates"].items():
        A(f"| {v} | {_fmt(s['sent_drug']['mean'], 2)} | {_fmt(s['sent_drug_disease_set']['mean'], 2)} | {s['sent_drug_disease_set']['pairs_with_ge1']} ({_pct(s['sent_drug_disease_set']['share_pairs_with_ge1'])}) |")
    r = q2["recall"]
    A("\n### Recall of the current retrieval against a direct `drug AND disease` query\n")
    A(f"Pairs {r['pairs']}; pairs with zero direct hits {r['pairs_with_zero_direct_hits']}; direct hits per pair `{json.dumps(r['hit_count_distribution'])}`.\n")
    A(f"Recall of the first 20 direct records among the 300 retrieved (INN query): `{json.dumps(r['recall_of_first20_records_in_the_300_retrieved_INN'])}`; pairs with recall 0: {r['pairs_with_recall_zero']}, recall 1: {r['pairs_with_recall_one']}; "
      f"synonym query: `{json.dumps(r['recall_first20_synonym_query'])}`.\n")
    A(f"Share of direct hits that also match the disease gene clause (what any gene-based retrieval can reach): `{json.dumps(r['share_of_direct_hits_that_also_match_the_gene_clause'])}`; "
      f"total hits {r['hits_total']}, with gene clause {r['hits_with_gene_clause_total']}.\n")
    A(f"Gene query cap: `{json.dumps(res['Q2_cap'])}`\n")
    A("### Language\n")
    A(f"Retrieved records (30 seeded diseases, language field of the core record): `{json.dumps(res['Q2_language_retrieved']['pooled_retrieved_records_language'])}`; "
      f"by disease `{json.dumps(res['Q2_language_retrieved']['non_english_share_of_retrieved_by_disease'])}`; share of gene-query hits not English `{json.dumps(res['Q2_language_retrieved']['gene_query_hits_not_english_share'])}`.\n")
    A(f"Direct query: `{json.dumps(q2['language_direct_query'])}`\n")
    A("### Drug naming\n")
    A(f"`{json.dumps(res['Q2_drug_naming'], ensure_ascii=False)}`\n")
    q3 = res["Q3_labels"]
    A("## Q3. Label audit\n")
    a = q3["a_reading"]
    A(f"### (a) Reading of gold-positive pairs with no `yes` from either annotator\n")
    A(f"Pool: {a['pool']['n_pairs_no_yes_from_either_annotator']} of {a['pool']['n_gold_positive_pairs_with_items']} gold-positive pairs with items. Reader: {a['reader']}. Class counts: `{json.dumps(a['class_counts'])}`.\n")
    A("| # | drug | disease | sentence | O (rel/dir) | S (rel/dir) | class | reason |")
    A("|---|---|---|---|---|---|---|---|")
    for k, row in enumerate(a["rows"], 1):
        sent = row["sentence"].replace("|", "/")
        A(f"| {k} | {row['drug']} | {row['disease']} | {sent} | {row['opus'][0]}/{row['opus'][1]} | {row['sonnet'][0]}/{row['sonnet'][1]} | {row['class'] or '-'} | {row['reason'] or '-'} |")
    A("\n### (b) Opus x Sonnet cross-tabs (rows Opus, columns Sonnet)\n")
    for g, v in q3["b_crosstabs"].items():
        A(f"**{g}** (n = {v['n_items']}): relevance agreement {_pct(v['relevance']['agreement'])}, kappa {_fmt(v['relevance']['kappa'], 3)}; direction agreement {_pct(v['direction']['agreement'])}, "
          f"kappa {_fmt(v['direction']['kappa'], 3)}; share yes `{json.dumps(v['share_yes'])}`, either {_pct(v['share_yes_either'])}, both {_pct(v['share_yes_both'])}\n")
        A(f"- relevance: `{json.dumps(v['relevance']['rows_first_columns_second'])}`")
        A(f"- direction: `{json.dumps(v['direction']['rows_first_columns_second'])}`\n")
    A("### (c) Coverage under different relevance rules (pairs, mean over replicates)\n")
    for fr, v in q3["c_sensitivity"].items():
        A(f"f = {fr}: `{json.dumps(v)}`\n")
    p = q3["d_v2_pilot"]
    A("### (d) v2 pilot on the same 800 items\n")
    A(f"Disagreement on relevance: `{json.dumps(p['relevance_disagreement'])}`; Opus vs Sonnet relevance kappa v1 {_fmt(p['opus_vs_sonnet_relevance']['v1']['kappa'], 3)} -> v2 {_fmt(p['opus_vs_sonnet_relevance']['v2']['kappa'], 3)} "
      f"(agreement {_pct(p['opus_vs_sonnet_relevance']['v1']['agreement'])} -> {_pct(p['opus_vs_sonnet_relevance']['v2']['agreement'])}); direction kappa {_fmt(p['opus_vs_sonnet_direction']['v1']['kappa'], 3)} -> {_fmt(p['opus_vs_sonnet_direction']['v2']['kappa'], 3)}; "
      f"yes counts `{json.dumps(p['yes_counts'])}`.\n")
    for an in ("opus", "sonnet"):
        A(f"- {an} relevance v1 (rows) vs v2 (columns): `{json.dumps(p['per_annotator_v1_vs_v2'][an]['relevance']['rows_first_columns_second'])}`; direction `{json.dumps(p['per_annotator_v1_vs_v2'][an]['direction']['rows_first_columns_second'])}`")
    A(f"- Opus x Sonnet relevance v1 `{json.dumps(p['opus_vs_sonnet_relevance']['v1']['rows_first_columns_second'])}` v2 `{json.dumps(p['opus_vs_sonnet_relevance']['v2']['rows_first_columns_second'])}`\n")
    A("## Q4. Verdict\n")
    A(f"Computed summary: `{json.dumps(res['Q4_verdict']['computed'], ensure_ascii=False)}`\n")
    for para in res["Q4_verdict"]["prose"] or VERDICT_PROSE:
        A(para + "\n")
    return "\n".join(L) + "\n"


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "sample":
        cmd_sample()
    elif cmd == "run":
        run()
    elif cmd == "md":
        res = json.loads(OUT_JSON.read_text(encoding="utf-8"))
        res["Q4_verdict"] = {"computed": verdict_skeleton(res), "prose": VERDICT_PROSE}
        v4.write_json(OUT_JSON, res)
        OUT_MD.write_text(render_md(res), encoding="utf-8", newline="\n")
    else:
        raise SystemExit("usage: python -m track2.engine_v4_llm_audit sample|run|md")


if __name__ == "__main__":
    main()
