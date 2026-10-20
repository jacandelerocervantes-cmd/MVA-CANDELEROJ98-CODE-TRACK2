"""Empirical positive and negative control for the literature layer of the repurposing report.

DETERMINISTIC, MENTION-BASED PROXY (no AI labelling). For each of ~26 rare genetic diseases that
have at least one approved drug indicated for them, the procedure of the real report is imitated
mechanically:

  1. query Europe PMC with 2-4 gene/pathway terms (never drug names), keep up to 200 abstracts;
  2. a drug "has a benefit sentence" when ONE sentence contains (a) the drug name (whole word,
     case-insensitive), (b) a benefit-verb stem from BENEFIT_STEMS and (c) a mechanism term of the
     disease query. The per-drug score is the number of such sentences (the proxy for rule R2:
     "at least one verified finding that may help", i.e. score >= 1);
  3. all approved drugs of data/universe.json are ranked by score; the gold standard is the set of
     universe drugs with an `indication` edge to the disease in PrimeKG (data/sources/primekg;
     DrugBank/DrugCentral derived). data/universe.json itself has no per-disease indications.

Measured: AUC (ties = 0.5), precision@10 (ties broken fractionally), hit rate (gold drugs with
score >= 1), enrichment of gold hits over 1000 random drug sets of the same size, and the
negative control (fraction of random NON-indicated drugs that pass score >= 1, and fraction of
random 9-drug sets with at least one passing drug). Confidence intervals: bootstrap over diseases
(macro statistics) and Wilson (pooled proportions).

WHAT THIS DOES NOT MEASURE: it does not read or judge sentences (no direction, no species, no
evidence level, no R1/R3/R4), it only counts co-mention; it does not test the quality of the nine
real proposals; the gold standard is incomplete (a non-indicated drug may be a true repurposing
success, so false positives are over-counted) and partly circular for drugs developed against the
disease target itself (easy positives, unlike repurposing). Drug names only: the universe has no
synonym table, so brand names and salts are missed (recall is underestimated). A random-drug
pass rate measures volume-of-literature bias of this proxy, not of the human/AI procedure.

Run: PYTHONPATH=src python -m track2.empirical_literature_control
"""
from __future__ import annotations

import hashlib
import html
import json
import math
import random
import re
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

ROOT = Path(__file__).resolve().parents[2]
SEED = 20260514
N_RANDOM_SETS = 1000
N_BOOT = 2000
MAX_ABSTRACTS = 200
N_DISEASES = 26
REPORT_SET_SIZE = 9  # the real report proposed nine drugs
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Fixed benefit-verb lexicon (stems, matched at a word start, case-insensitive).
BENEFIT_STEMS = (
    "improv", "rescu", "restor", "ameliorat", "allevi", "attenuat", "protect", "prevent",
    "benefit", "efficac", "effective", "reliev", "normali", "reduc", "suppress", "reverse",
    "stabili", "enhanc", "increas", "correct", "lower", "decreas", "inhibit", "ameliora",
)

# Candidate pool: (PrimeKG disease name, 2-4 gene/pathway query terms; no drug names).
# 26 are drawn with SEED among those that have >= 1 indicated universe drug in PrimeKG.
CANDIDATES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("cystic fibrosis", ("CFTR", "chloride channel")),
    ("Gaucher disease", ("GBA", "glucocerebrosidase", "glucosylceramide")),
    ("Fabry disease", ("GLA", "alpha-galactosidase", "globotriaosylceramide")),
    ("Duchenne muscular dystrophy", ("dystrophin", "exon skipping")),
    ("Huntington disease", ("huntingtin", "polyglutamine")),
    ("Wilson disease", ("ATP7B", "copper")),
    ("Dravet syndrome", ("SCN1A", "sodium channel")),
    ("spinal muscular atrophy", ("SMN1", "SMN2")),
    ("tuberous sclerosis", ("TSC1", "TSC2", "mTOR")),
    ("neurofibromatosis", ("NF1", "RAS", "MAPK")),
    ("phenylketonuria", ("PAH", "phenylalanine hydroxylase")),
    ("cystinosis", ("CTNS", "lysosomal cystine")),
    ("tyrosinemia", ("FAH", "tyrosine catabolism")),
    ("hereditary angioedema with C1Inh deficiency", ("SERPING1", "C1 inhibitor", "bradykinin")),
    ("homozygous familial hypercholesterolemia", ("LDLR", "PCSK9")),
    ("hypophosphatasia", ("ALPL", "alkaline phosphatase")),
    ("Leber hereditary optic neuropathy", ("MT-ND4", "mitochondrial complex I")),
    ("Marfan syndrome", ("FBN1", "TGF-beta")),
    ("Niemann-Pick disease", ("NPC1", "cholesterol trafficking")),
    ("Prader-Willi syndrome", ("SNRPN", "MAGEL2", "imprinting")),
    ("Berardinelli-Seip congenital lipodystrophy", ("AGPAT2", "BSCL2", "leptin")),
    ("mucopolysaccharidosis", ("glycosaminoglycan", "lysosomal storage", "iduronidase")),
    ("beta thalassemia", ("HBB", "globin", "ineffective erythropoiesis")),
    ("sickle cell anemia", ("HBB", "hemoglobin S", "polymerization")),
    ("hemophilia", ("F8", "F9", "coagulation factor")),
    ("paroxysmal nocturnal hemoglobinuria", ("PIGA", "complement")),
    ("hereditary hypophosphatemic rickets", ("PHEX", "FGF23")),
    ("urea cycle disorder", ("OTC", "CPS1", "ammonia")),
    ("cystinuria", ("SLC3A1", "SLC7A9")),
    ("polycystic kidney disease", ("PKD1", "PKD2")),
    ("myotonia congenita", ("CLCN1", "chloride channel")),
    ("hereditary hemochromatosis", ("HFE", "hepcidin")),
    ("acute intermittent porphyria", ("HMBS", "ALAS1")),
    ("erythropoietic protoporphyria", ("FECH", "protoporphyrin")),
    ("glycogen storage disease", ("GAA", "glycogen", "acid alpha-glucosidase")),
)


# ----------------------------------------------------------------------------- text scoring
def split_sentences(text: str) -> List[str]:
    text = re.sub(r"<[^>]+>", " ", html.unescape(text or ""))
    text = re.sub(r"\s+", " ", text).strip()
    return [s for s in re.split(r"(?<=[.!?])\s+(?=[A-Z0-9(\[])", text) if s]


def term_regex(terms: Sequence[str]) -> "re.Pattern[str]":
    parts = [r"(?<![\w])" + re.escape(t) + r"(?![\w])" for t in sorted(terms, key=len, reverse=True)]
    return re.compile("|".join(parts), re.IGNORECASE)


BENEFIT_RE = re.compile(r"(?<![\w])(?:" + "|".join(BENEFIT_STEMS) + r")\w*", re.IGNORECASE)


def drug_regex(names: Sequence[str]) -> "re.Pattern[str]":
    clean = sorted({n.lower() for n in names if n and len(n) >= 4}, key=len, reverse=True)
    return term_regex(clean)


def score_abstracts(abstracts: Sequence[Tuple[str, str]], mech_terms: Sequence[str],
                    drug_re: "re.Pattern[str]") -> Dict[str, int]:
    """Count, per drug (lower case), the sentences with drug + benefit stem + mechanism term.

    abstracts: sequence of (id, text). A sentence counts once per drug.
    """
    mech_re = term_regex(mech_terms)
    counts: Dict[str, int] = {}
    for _pid, text in abstracts:
        for sent in split_sentences(text):
            if not (BENEFIT_RE.search(sent) and mech_re.search(sent)):
                continue
            for name in {m.group(0).lower() for m in drug_re.finditer(sent)}:
                counts[name] = counts.get(name, 0) + 1
    return counts


# ----------------------------------------------------------------------------- metrics
def auc_ties(scores: Dict[str, float], pos: Set[str], universe: Sequence[str]) -> Optional[float]:
    """Mann-Whitney AUC, ties count 1/2. None when a class is empty."""
    p = [scores.get(d, 0.0) for d in universe if d in pos]
    n = [scores.get(d, 0.0) for d in universe if d not in pos]
    if not p or not n:
        return None
    ns = sorted(n)
    from bisect import bisect_left, bisect_right
    tot = 0.0
    for x in p:
        tot += bisect_left(ns, x) + 0.5 * (bisect_right(ns, x) - bisect_left(ns, x))
    return tot / (len(p) * len(n))


def precision_at_k(scores: Dict[str, float], pos: Set[str], universe: Sequence[str], k: int = 10) -> float:
    """Expected precision@k when ties are broken uniformly at random (exact, group-wise)."""
    groups: Dict[float, List[int]] = {}
    for d in universe:
        g = groups.setdefault(scores.get(d, 0.0), [0, 0])
        g[0] += 1
        g[1] += 1 if d in pos else 0
    left, exp_pos = k, 0.0
    for s in sorted(groups, reverse=True):
        size, npos = groups[s]
        take = min(left, size)
        exp_pos += take * npos / size
        left -= take
        if left <= 0:
            break
    return exp_pos / k


def wilson(k: int, n: int, z: float = 1.96) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if n == 0:
        return None, None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def boot_ci(values: Sequence[float], rng: random.Random, n_boot: int = N_BOOT) -> Dict[str, Optional[float]]:
    v = [x for x in values if x is not None]
    if not v:
        return {"mean": None, "lo": None, "hi": None, "n": 0}
    means = []
    for _ in range(n_boot):
        means.append(sum(rng.choice(v) for _ in v) / len(v))
    means.sort()
    return {"mean": sum(v) / len(v), "lo": means[int(0.025 * n_boot)], "hi": means[int(0.975 * n_boot) - 1],
            "n": len(v)}


def random_set_stats(scores: Dict[str, float], universe: Sequence[str], size: int, n_sets: int,
                     rng: random.Random, exclude: Set[str] = frozenset()) -> List[int]:
    """Number of drugs with score >= 1 in each random set of `size` drugs."""
    pool = [d for d in universe if d not in exclude]
    return [sum(1 for d in rng.sample(pool, size) if scores.get(d, 0) >= 1) for _ in range(n_sets)]


def evaluate_disease(scores: Dict[str, float], gold: Set[str], universe: Sequence[str],
                     rng: random.Random, n_sets: int = N_RANDOM_SETS) -> Dict[str, object]:
    gold = {g for g in gold if g in set(universe)}
    hits = sum(1 for g in gold if scores.get(g, 0) >= 1)
    null_hits = random_set_stats(scores, universe, len(gold), n_sets, rng)
    mean_null = sum(null_hits) / n_sets
    p_emp = (1 + sum(1 for x in null_hits if x >= hits)) / (n_sets + 1)
    # negative control: random NON-indicated drugs
    non_gold = [d for d in universe if d not in gold]
    n_pass = sum(1 for d in non_gold if scores.get(d, 0) >= 1)
    sets9 = random_set_stats(scores, universe, REPORT_SET_SIZE, n_sets, rng, exclude=gold)
    return {
        "n_gold": len(gold), "gold_hits": hits,
        "auc": auc_ties(scores, gold, universe),
        "precision_at_10": precision_at_k(scores, gold, universe, 10),
        "null_mean_hits": mean_null, "null_p_ge_observed": p_emp,
        "enrichment_hits": (hits / mean_null) if mean_null > 0 else None,
        "n_nongold": len(non_gold), "nongold_pass": n_pass,
        "frac_sets9_with_pass": sum(1 for x in sets9 if x >= 1) / n_sets,
    }


def aggregate(per: Dict[str, Dict[str, object]], rng: random.Random) -> Dict[str, object]:
    ds = sorted(per)
    col = lambda k: [per[d][k] for d in ds if per[d].get(k) is not None]
    g_hits, g_n = sum(per[d]["gold_hits"] for d in ds), sum(per[d]["n_gold"] for d in ds)
    ng_pass, ng_n = sum(per[d]["nongold_pass"] for d in ds), sum(per[d]["n_nongold"] for d in ds)
    gp = wilson(g_hits, g_n)
    npp = wilson(ng_pass, ng_n)
    # pooled enrichment: bootstrap over diseases of sum(hits)/sum(null mean hits)
    ratios = []
    for _ in range(N_BOOT):
        pick = [rng.choice(ds) for _ in ds]
        num = sum(per[d]["gold_hits"] for d in pick)
        den = sum(per[d]["null_mean_hits"] for d in pick)
        if den > 0:
            ratios.append(num / den)
    ratios.sort()
    tot_null = sum(per[d]["null_mean_hits"] for d in ds)
    return {
        "n_diseases": len(ds),
        "auc_macro": boot_ci(col("auc"), rng),
        "precision_at_10_macro": boot_ci(col("precision_at_10"), rng),
        "frac_diseases_auc_gt_0.5": boot_ci([1.0 if a > 0.5 else 0.0 for a in col("auc")], rng),
        "frac_diseases_null_p_le_0.05": boot_ci([1.0 if p <= 0.05 else 0.0 for p in col("null_p_ge_observed")], rng),
        "gold_hit_rate_pooled": {"k": g_hits, "n": g_n, "rate": gp[0], "lo": gp[1], "hi": gp[2]},
        "pooled_enrichment_hits": {
            "value": (g_hits / tot_null) if tot_null > 0 else None,
            "lo": ratios[int(0.025 * len(ratios))] if ratios else None,
            "hi": ratios[int(0.975 * len(ratios)) - 1] if ratios else None},
        "negative_control_pass_rate_random_drug": {
            "k": ng_pass, "n": ng_n, "rate": npp[0], "lo": npp[1], "hi": npp[2],
            "threshold": "score >= 1 benefit sentence (proxy of R2)"},
        "negative_control_frac_random_sets9_with_any_pass": boot_ci(col("frac_sets9_with_pass"), rng),
    }


# ----------------------------------------------------------------------------- data access
def select_diseases(available: Dict[str, Set[str]], seed: int = SEED, n: int = N_DISEASES,
                    candidates=CANDIDATES) -> List[Tuple[str, Tuple[str, ...]]]:
    """Seeded draw of n candidates that have >= 1 indicated universe drug; sorted by name."""
    ok = [c for c in candidates if available.get(c[0])]
    rnd = random.Random(seed)
    pick = rnd.sample(ok, min(n, len(ok)))
    return sorted(pick, key=lambda c: c[0])


def load_universe(path: Path) -> List[dict]:
    return json.loads(path.read_text(encoding="utf-8"))["drugs"]


def load_primekg_indications(primekg_dir: Path, universe_names: Set[str]) -> Dict[str, Set[str]]:
    """disease name -> set of lower-case universe drug names with an `indication` edge."""
    import pandas as pd
    nodes = pd.read_csv(primekg_dir / "nodes.tab")
    ntype, nname = nodes["node_type"].to_dict(), nodes["node_name"].to_dict()
    # node_index is the row label in PrimeKG; map explicitly
    idx = nodes["node_index"].tolist()
    ntype, nname = dict(zip(idx, nodes["node_type"])), dict(zip(idx, nodes["node_name"]))
    out: Dict[str, Set[str]] = {}
    for ch in pd.read_csv(primekg_dir / "edges.csv", chunksize=2_000_000,
                          usecols=["relation", "x_index", "y_index"]):
        ch = ch[ch["relation"] == "indication"]
        for x, y in zip(ch["x_index"], ch["y_index"]):
            if ntype[x] == "drug" and ntype[y] == "disease" and nname[x].lower() in universe_names:
                out.setdefault(nname[y], set()).add(nname[x].lower())
    return out


def epmc_fetch(query: str, cache_dir: Path, max_n: int = MAX_ABSTRACTS,
               getter: Optional[Callable[[str, dict], dict]] = None) -> List[Tuple[str, str]]:
    """Up to max_n (id, abstract) pairs for `query`, cached on disk (polite: <= 5 requests/s)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(f"{query}|{max_n}|core".encode()).hexdigest()[:24]
    f = cache_dir / f"{key}.json"
    if f.exists():
        return [tuple(x) for x in json.loads(f.read_text(encoding="utf-8"))]
    if getter is None:
        import requests

        def getter(url, params):  # type: ignore[no-redef]
            for attempt in range(6):
                try:
                    time.sleep(0.25)
                    r = requests.get(url, params=params, timeout=60)
                    r.raise_for_status()
                    return r.json()
                except Exception:  # SSL / transient errors: exponential backoff
                    if attempt == 5:
                        raise
                    time.sleep(2 ** attempt)
    res: List[Tuple[str, str]] = []
    page = 1
    while len(res) < max_n:
        data = getter(EPMC, {"query": f"({query}) AND HAS_ABSTRACT:y", "format": "json",
                             "resultType": "core", "pageSize": 100, "page": page})
        items = (data.get("resultList") or {}).get("result") or []
        for it in items:
            if it.get("abstractText"):
                res.append((str(it.get("id") or it.get("pmid") or ""), it["abstractText"]))
        if len(items) < 100:
            break
        page += 1
    res = res[:max_n]
    f.write_text(json.dumps(res), encoding="utf-8")
    return res


def build_query(terms: Sequence[str]) -> str:
    return " OR ".join(f'"{t}"' for t in terms)


# ----------------------------------------------------------------------------- driver
def run(root: Path = ROOT, getter=None) -> Dict[str, object]:
    drugs = load_universe(root / "data" / "universe.json")
    names = sorted({d["name"].lower() for d in drugs})
    small = {d["name"].lower() for d in drugs if d.get("molecule_type") == "Small molecule"}
    avail = load_primekg_indications(root / "data" / "sources" / "primekg", set(names))
    chosen = select_diseases(avail)
    dre = drug_regex(names)
    rng = random.Random(SEED)
    cache = root / "data" / "cache_empirical"
    per, per_small, detail = {}, {}, {}
    for dname, terms in chosen:
        abstracts = epmc_fetch(build_query(terms), cache, getter=getter)
        counts = score_abstracts(abstracts, terms, dre)
        gold = avail[dname]
        per[dname] = evaluate_disease(counts, gold, names, rng)
        gs = gold & small
        if gs:
            per_small[dname] = evaluate_disease(counts, gs, names, rng)
        top = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:10]
        detail[dname] = {"query_terms": list(terms), "n_abstracts": len(abstracts),
                         "gold_drugs": sorted(gold), "top10_by_score": top,
                         "gold_scores": {g: counts.get(g, 0) for g in sorted(gold)}}
        per[dname].update(detail[dname])
    out = {
        "method": "mention-based deterministic proxy; see module docstring",
        "seed": SEED, "n_random_sets": N_RANDOM_SETS, "n_bootstrap": N_BOOT,
        "max_abstracts": MAX_ABSTRACTS, "n_universe_drugs": len(names),
        "benefit_stems": list(BENEFIT_STEMS),
        "gold_standard": "PrimeKG `indication` edges restricted to data/universe.json drugs "
                         "(data/universe.json has no per-disease indications)",
        "aggregate_all_gold": aggregate({d: {k: v for k, v in r.items() if not isinstance(v, (list, dict))}
                                         for d, r in per.items()}, rng),
        "aggregate_small_molecule_gold": aggregate(per_small, rng) if per_small else None,
        "per_disease": per,
        "limitations": [
            "Counts co-mention only; no AI/human reading, direction or evidence level; rules R1,R3,R4 not modelled.",
            "Gold standard incomplete: unlisted drugs may be true repurposing successes (false positives over-counted).",
            "Many gold drugs were developed against the disease mechanism itself, so recovery is easier than for repurposing.",
            "Drug names only (no synonyms/brands): recall underestimated.",
            "Random-drug pass rate measures this proxy's volume-of-literature bias, not that of the real procedure.",
            "Abstracts only (max 200 by Europe PMC relevance); full texts not used.",
        ],
    }
    outp = root / "out" / "empirical" / "literature_control.json"
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(out, indent=1, sort_keys=True), encoding="utf-8")
    return out


if __name__ == "__main__":
    res = run()
    print(json.dumps({k: res[k] for k in ("aggregate_all_gold", "aggregate_small_molecule_gold")},
                     indent=1, sort_keys=True))
