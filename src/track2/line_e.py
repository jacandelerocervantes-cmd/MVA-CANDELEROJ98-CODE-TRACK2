"""Line E, signature reversal (protocol/addendum_line_e_signature_reversal.md, written before any score).

Run from the project root, reading the downloaded files with an isolated interpreter:
    python -I -c "import sys; sys.path.insert(0, 'src'); from track2 import line_e; line_e.main()"
Exploratory; hypotheses only; never evidence of efficacy.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import os

import numpy as np

from track2.fda_approved import candidate_names

SEED, N_NULL, N_PRIMARY, FDR_CUT, FLOOR, HIGH_SCORE = 20261018, 1000, 150, 0.1, 0.05, 0.10
P_SIZE, P_MIN_D, P_SEED = 100, 0.5, 20261019
GPL_FILE = "data/geo/GPL6104/GPL6104.annot.gz"
KO_FILE = "data/lincs/crispr_ko/LINCS_L1000_CRISPR_KO_Consensus_Sigs.txt"
CHEM_FILE = "data/lincs/chem_pert/LINCS_L1000_Chem_Pert_Consensus_Sigs.txt"
GEO_FILE = "data/geo/GSE22206/GSE22206_series_matrix.txt.gz"
TRACK1_CSV = "../HSelith_MVA-Track1/submissions/track1/candeleroj98_genomewide_ranking.csv"  # relative to the project root; override with the TRACK1_RANKING_CSV environment variable
TRACK1_COLUMNS = ["rank", "gene", "type", "n_high_impact", "log10_joint_rarity", "phenotype_score"]  # the 'variants' column is never read
EXTRA_MECHANISM = {"CENPE"}
ALIASES = {"PHENYLBUTYRATE": "PHENYLBUTANOIC ACID", "SAHA": "VORINOSTAT"}
NINE = ["sirolimus", "everolimus", "metformin", "sodium phenylbutyrate", "trametinib", "selumetinib", "bortezomib", "chloroquine", "bosutinib"]
NINE_LIB = {"sodium phenylbutyrate": "phenylbutyrate"}
CONTROL_PAIRS = {"MTOR": ["sirolimus", "everolimus", "temsirolimus"], "MAP2K1": ["trametinib", "selumetinib", "cobimetinib", "binimetinib"],
                 "EGFR": ["erlotinib", "gefitinib", "afatinib", "lapatinib", "osimertinib"], "BRAF": ["vemurafenib", "dabrafenib"],
                 "HMGCR": ["atorvastatin", "simvastatin", "lovastatin", "fluvastatin", "rosuvastatin", "pitavastatin", "pravastatin"],
                 "ABL1": ["imatinib", "dasatinib", "nilotinib", "bosutinib"]}


# ---------------------------------------------------------------- parsing
def parse_library(text: str) -> dict[str, dict[str, list[str]]]:
    """Enrichr text library -> {name: {'Up': [...], 'Down': [...]}} (gene order kept)."""
    lib: dict[str, dict[str, list[str]]] = {}
    for line in text.split("\n"):
        if not line.strip():
            continue
        f = line.rstrip("\r").split("\t")
        term = f[0].strip()
        if term.endswith(" Up"):
            name, d = term[:-3], "Up"
        elif term.endswith(" Down"):
            name, d = term[:-5], "Down"
        else:
            continue
        lib.setdefault(name, {})[d] = [g for g in f[2:] if g]
    return {k: v for k, v in lib.items() if "Up" in v and "Down" in v}


def signed_vector(up: Iterable[str], down: Iterable[str], n: int | None, index: dict[str, int]) -> np.ndarray:
    """+1 for the first n up genes, -1 for the first n down genes (a gene in both is set to 0)."""
    v = np.zeros(len(index), dtype=np.float32)
    u = [g for g in list(up)[:n] if g in index]
    d = [g for g in list(down)[:n] if g in index]
    for g in u:
        v[index[g]] = 1.0
    for g in d:
        v[index[g]] = 0.0 if v[index[g]] == 1.0 else -1.0
    return v


def reversal(drugs: np.ndarray, sig: np.ndarray) -> np.ndarray:
    """R = |Dup&Sdown| + |Ddown&Sup| - |Dup&Sup| - |Ddown&Sdown| = (Dup - Ddown) . (Sdown - Sup); drugs (n_drugs x genes), sig (genes) or (genes x k)."""
    return drugs @ (-sig)


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def test_signature(drugs: np.ndarray, sig: np.ndarray, null_r: np.ndarray) -> dict[str, np.ndarray]:
    """null_r: (n_drugs x n_null) reversal scores of the control signatures."""
    r = reversal(drugs, sig)
    n = null_r.shape[1]
    p_rev = (1 + (null_r >= r[:, None]).sum(1)) / (n + 1)
    p_mim = (1 + (null_r <= r[:, None]).sum(1)) / (n + 1)
    q_rev, q_mim = bh(p_rev), bh(p_mim)
    mean = null_r.mean(1)
    return {"R": r, "p_rev": p_rev, "p_mim": p_mim, "q_rev": q_rev, "q_mim": q_mim, "null_mean": mean,
            "usable": (q_rev < FDR_CUT) & (r > 0) & (r > mean), "hazard": (q_mim < FDR_CUT) & (r < 0) & (r < mean)}


def combined_vote(lib_ko: dict[str, dict[str, list[str]]], weights: dict[str, float], n: int | None, index: dict[str, int]) -> np.ndarray:
    """Weighted signed vote over the knockout signatures of the weighted genes; positive net -> Up, negative -> Down; truncated to n each side."""
    score = np.zeros(len(index), dtype=np.float64)
    order = np.zeros(len(index), dtype=np.float64)
    for g, w in weights.items():
        if g not in lib_ko:
            continue
        for sign, key in ((1.0, "Up"), (-1.0, "Down")):
            genes = lib_ko[g][key]
            for rank, x in enumerate(genes):
                if x in index:
                    score[index[x]] += sign * w
                    order[index[x]] += sign * w * (1.0 - rank / (len(genes) + 1)) * 1e-6  # tie-break by within-signature rank
    score = score + order
    v = np.zeros(len(index), dtype=np.float32)
    pos = np.where(score > 0)[0]
    neg = np.where(score < 0)[0]
    v[pos[np.argsort(-score[pos], kind="stable")][:n]] = 1.0
    v[neg[np.argsort(score[neg], kind="stable")][:n]] = -1.0
    return v


def track1_weights(rows: list[dict[str, Any]], extra_genes: Iterable[str], min_score: float = 0.0) -> dict[str, float]:
    """w = phenotype_score with a floor of FLOOR; structural-list genes with no score get FLOOR. min_score drops genes below the stated threshold."""
    w: dict[str, float] = {}
    for r in rows:
        s = float(r.get("phenotype_score") or 0.0)
        if s >= min_score:
            w[r["gene"]] = max(s, FLOOR)
    if min_score <= 0:
        for g in extra_genes:
            w.setdefault(g, FLOOR)
    return w


# ---------------------------------------------------------------- drugs
def match_drugs(lib_names: Iterable[str], universe: list[dict[str, Any]]) -> tuple[dict[str, dict[str, Any]], list[str]]:
    """library name -> universe drug (small molecules only); exact upper name, then salt-stripped names, then the alias table."""
    small = [d for d in universe if d.get("molecule_type") == "Small molecule"]
    by_name: dict[str, dict[str, Any]] = {}
    for d in small:
        for nm in candidate_names(d["name"])[:2]:
            by_name.setdefault(nm, d)
        by_name[d["name"].upper()] = d
    matched, unmatched = {}, []
    for n in lib_names:
        u = ALIASES.get(n.upper(), n.upper())
        hit = next((by_name[c] for c in candidate_names(u)[:2] if c in by_name), None)
        if hit is None:
            unmatched.append(n)
        else:
            matched[n] = hit
    # one library entry per drug: keep the first
    seen, out = set(), {}
    for n, d in matched.items():
        if d["drug_id"] not in seen:
            seen.add(d["drug_id"])
            out[n] = d
    return out, unmatched


def drug_matrix(lib_chem: dict[str, dict[str, list[str]]], names: list[str], n: int | None, index: dict[str, int]) -> np.ndarray:
    return np.stack([signed_vector(lib_chem[x]["Up"], lib_chem[x]["Down"], n, index) for x in names]).astype(np.float32)


def candidate_table(rows: list[dict[str, Any]], lines_usable: dict[str, set[str]], safety: dict[str, Any]) -> list[dict[str, Any]]:
    """Combination rule of section 7 over line-E-usable drugs. rows: dicts with name, e_usable, atc, fda_label."""
    out = []
    for r in rows:
        if not r["e_usable"]:
            continue
        name = r["name"].upper()
        existing = sorted(k for k, s in lines_usable.items() if name in s)
        cytotoxic = any(a.startswith("L01") for a in r.get("atc", [])) or name in safety["caution_names"]
        ok = bool(existing) and not cytotoxic and bool(r["fda_label"])
        out.append({"name": name, "existing_lines": existing, "oncology_cytotoxic_or_caution_class": cytotoxic, "fda_label": bool(r["fda_label"]),
                    "proposal_safety_entries": safety["proposal"].get(name, []), "candidate": ok, "age_support": r.get("age_support"), "label_level": r.get("label_level")})
    return out


def read_annotation(path: Path) -> dict[str, str]:
    """GPL6104.annot(.gz) -> {probe: gene symbol}; probes without a symbol or with several symbols ('///') are dropped."""
    op = gzip.open if str(path).endswith(".gz") else open
    out: dict[str, str] = {}
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("!platform_table_begin"):
                header = next(fh).rstrip("\n").split("\t")
                continue
            if header is None or line.startswith("!platform_table_end"):
                continue
            f = line.rstrip("\n").split("\t")
            row = dict(zip(header, f))
            sym = (row.get("Gene symbol") or "").strip()
            if sym and "///" not in sym:
                out[row["ID"]] = sym
    return out


def geo_signature(matrix_gz: Path, annotation: dict[str, str] | None, n: int = P_SIZE, min_abs_d: float = P_MIN_D) -> dict[str, Any]:
    """Sections 6 and 10: LCL proband (2 arrays) versus unrelated-control LCL (4 arrays); without annotation nothing is computed."""
    with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().split("\n")
    titles = [t.strip('"') for t in next(l for l in lines if l.startswith("!Sample_title")).split("\t")[1:]]
    meta = {"n_arrays": len(titles), "array_titles": titles}
    if not annotation:
        return {"computed": False, "reason": "probe-to-gene annotation (GPL6104) not available; signature not computed", **meta}
    i = next(k for k, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    case = [k for k, t in enumerate(titles) if t.startswith("LCL-proband")]
    ctrl = [k for k, t in enumerate(titles) if t.startswith("LCL-control")]
    use = case + ctrl
    ids, mat = [], []
    for l in lines[i + 2:]:
        if not l or l.startswith("!"):
            continue
        f = l.split("\t")
        ids.append(f[0].strip('"'))
        mat.append([float(f[1 + k]) if f[1 + k] not in ("", "null") else np.nan for k in use])
    m = np.array(mat)
    nc = len(case)
    best: dict[str, int] = {}
    mean = np.nanmean(m, axis=1)
    for r, pid in enumerate(ids):
        g = annotation.get(pid)
        if g and np.isfinite(mean[r]) and (g not in best or mean[r] > mean[best[g]]):
            best[g] = r
    genes = sorted(best)
    g_m = m[[best[g] for g in genes]]
    a, b = g_m[:, :nc], g_m[:, nc:]
    d = a.mean(1) - b.mean(1)
    pooled = np.sqrt((a.var(1, ddof=1) * (a.shape[1] - 1) + b.var(1, ddof=1) * (b.shape[1] - 1)) / (a.shape[1] + b.shape[1] - 2) * (1 / a.shape[1] + 1 / b.shape[1]))
    t = d / (pooled + np.nanmedian(pooled))
    gmean = g_m.mean(1)
    expressed = gmean >= np.median(gmean)
    ok = [k for k in range(len(genes)) if expressed[k] and abs(d[k]) >= min_abs_d and np.isfinite(t[k])]
    up = [genes[k] for k in sorted((k for k in ok if d[k] > 0), key=lambda k: -t[k])[:n]]
    down = [genes[k] for k in sorted((k for k in ok if d[k] < 0), key=lambda k: t[k])[:n]]
    return {"computed": True, "up": up, "down": down, "n_genes_mapped": len(genes), "n_expressed": int(expressed.sum()), "n_passing_effect": len(ok),
            "expressed_genes": [genes[k] for k in range(len(genes)) if expressed[k]], "n_cases": nc, "n_controls": len(ctrl), "array_titles": titles, "n_arrays": len(titles)}


def random_signature_null(expressed: list[str], index: dict[str, int], n_up: int, n_down: int, n_null: int, seed: int) -> np.ndarray:
    """(genes x n_null) signed vectors of random disjoint n_up + n_down gene sets from the expressed genes, truncated to the library index."""
    rng = np.random.default_rng(seed)
    cols = []
    for _ in range(n_null):
        pick = rng.choice(len(expressed), size=n_up + n_down, replace=False)
        cols.append(signed_vector([expressed[k] for k in pick[:n_up]], [expressed[k] for k in pick[n_up:]], None, index))
    return np.stack(cols, axis=1)


# ---------------------------------------------------------------- main
def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _summ(t: dict[str, np.ndarray], names: list[str], i: int) -> dict[str, Any]:
    return {"R": int(t["R"][i]), "null_mean": round(float(t["null_mean"][i]), 2), "p_rev": float(t["p_rev"][i]), "q_rev": float(t["q_rev"][i]),
            "p_mimic": float(t["p_mim"][i]), "q_mimic": float(t["q_mim"][i]), "usable": bool(t["usable"][i]), "hazard": bool(t["hazard"][i])}


def main(root: str = ".") -> None:
    import pandas as pd

    root_p = Path(root)
    ko = parse_library((root_p / KO_FILE).read_text(encoding="utf-8"))
    chem = parse_library((root_p / CHEM_FILE).read_text(encoding="utf-8"))
    universe = json.loads((root_p / "data/universe.json").read_text(encoding="utf-8"))["drugs"]
    landscape = json.loads((root_p / "out/mva_landscape.json").read_text(encoding="utf-8"))["genes"]
    mva_genes = sorted(g["symbol"] for g in landscape if g["tier"] == "curated")
    mva_in_lib = [g for g in mva_genes if g in ko]
    t1df = pd.read_csv(os.environ.get("TRACK1_RANKING_CSV") or (root_p / TRACK1_CSV), usecols=TRACK1_COLUMNS)
    rows = t1df.to_dict("records")
    ext = json.loads((root_p / "config/track1_gene_level_extract.json").read_text(encoding="utf-8"))
    extra = sorted(set(ext["differential"]["structural"]) | set(ext["differential"]["ranking"]))
    t1_genes = sorted({r["gene"] for r in rows} | set(extra))
    t1_in_lib = [g for g in t1_genes if g in ko]
    weights = track1_weights(rows, extra)
    weights_high = track1_weights(rows, extra, HIGH_SCORE)
    mechanism = set(mva_genes) | EXTRA_MECHANISM

    # genes universe: every gene occurring in any set
    genes = sorted({g for lib in (ko, chem) for v in lib.values() for s in v.values() for g in s})
    index = {g: i for i, g in enumerate(genes)}
    matched, unmatched = match_drugs(list(chem), universe)
    names = list(matched)
    dmeta = [matched[n] for n in names]
    pool = sorted(g for g in ko if g not in set(mva_genes) | set(t1_genes))
    rng = np.random.default_rng(SEED)
    draw = list(rng.choice(pool, size=N_NULL, replace=False))

    # safety / label layers
    caution = json.loads((root_p / "out/caution_list.json").read_text(encoding="utf-8"))
    caution_names = {f["name"].upper() for c in caution["clases"].values() for f in c["farmacos"]}
    psafe = json.loads((root_p / "config/proposal_safety.json").read_text(encoding="utf-8"))["entries"]
    safety = {"caution_names": caution_names, "proposal": {k: [{"section": e["section"], "quote": e["quote"], "flag": e.get("flag")} for e in v] for k, v in psafe.items()}}
    fda = json.loads((root_p / "out/fda_approval_check.json").read_text(encoding="utf-8"))["drugs"]
    apply_mva = json.loads((root_p / "out/engine/apply_mva.json").read_text(encoding="utf-8"))
    age_tab = {f["name"].upper(): f for f in json.loads((root_p / "out/age_label_table.json").read_text(encoding="utf-8"))["farmacos"]}
    from track2.label_levels import LEVEL

    def fda_label(name: str) -> bool:
        u = name.upper()
        return bool(fda.get(u, {}).get("confirmed")) or any(k.upper() == u and v.get("fda_label") for k, v in apply_mva["label_support"].items()) or \
            any(k.upper() == u and v.get("eligible") for k, v in apply_mva["drugs_at_fda_checked"].items())

    def age_support(name: str) -> dict[str, Any]:
        u = name.upper()
        for k, v in apply_mva["label_support"].items():
            if k.upper() == u:
                return {"source": "out/engine/apply_mva.json", "bands": v["bands"], "unknown_age": "never assumed"}
        a = age_tab.get(u)
        if a:
            return {"source": "out/age_label_table.json (three groups: 0-<3, 3-<12, >=12)", "bands": a["estado_por_grupo"], "unknown_age": "never assumed"}
        return {"source": None, "bands": "not_assessed", "unknown_age": "never assumed"}

    # existing lines
    fl = json.loads((root_p / "out/four_lines_review2.json").read_text(encoding="utf-8"))
    sd = fl["support_detail"]
    lines_usable: dict[str, set[str]] = {}
    if fl["lines"]["A"]["usable"]:
        lines_usable["A"] = {k.upper() for k in sd["random_query_fraction_of_A_specific"]}
    if fl["lines"]["B1_phenotype_module"]["usable"]:
        lines_usable["B1"] = {k.upper() for k in sd["B1_robust"]}
    if fl["lines"]["B2_mitosis_module"]["usable"]:
        lines_usable["B2"] = {k.upper() for k in sd["B2_W3_two_record"]}
    if fl["lines"]["C"]["usable"]:
        lines_usable["C"] = set()

    # nulls per truncation
    result: dict[str, Any] = {"label": "exploratory; hypotheses for follow-up; never evidence of efficacy", "protocol": "protocol/addendum_line_e_signature_reversal.md",
                              "seed": SEED, "n_null": N_NULL, "fdr_cut": FDR_CUT}
    result["data"] = {"ko_library": {"n_signatures": len(ko), "sha256": _sha(root_p / KO_FILE)}, "chem_library": {"n_signatures": len(chem), "sha256": _sha(root_p / CHEM_FILE)},
                      "cell_line_specific_signatures": "not available (consensus libraries)"}
    n_universe_small = sum(1 for d in universe if d.get("molecule_type") == "Small molecule")
    result["drug_matching"] = {"n_library_drugs": len(chem), "n_matched_to_approved_small_molecules": len(names), "n_project_small_molecules": n_universe_small,
                               "n_unmatched_library_names": len(unmatched)}
    result["signature_availability"] = {"mva_curated_genes": mva_genes, "mva_in_library": mva_in_lib, "mva_absent_from_library": sorted(set(mva_genes) - set(mva_in_lib)),
                                        "track1_genes": t1_genes, "n_track1_genes": len(t1_genes), "track1_in_library": t1_in_lib, "n_track1_in_library": len(t1_in_lib),
                                        "track1_weights": {g: round(w, 4) for g, w in weights.items()}, "track1_high_score_weights": {g: round(w, 4) for g, w in weights_high.items()}}
    result["analysis"] = {}
    store: dict[int | None, Any] = {}
    for n in (N_PRIMARY, 50, None):
        D = drug_matrix(chem, names, n, index)
        null_sigs = np.stack([signed_vector(ko[g]["Up"], ko[g]["Down"], n, index) for g in draw], axis=1)
        null_r = reversal(D, null_sigs)
        sigs: dict[str, np.ndarray] = {}
        for g in sorted(set(mva_in_lib) | set(t1_in_lib) | {"BUB1B"}):
            sigs[f"gene:{g}"] = signed_vector(ko[g]["Up"], ko[g]["Down"], n, index)
        if mva_in_lib:
            sigs["mva_set"] = combined_vote(ko, {g: 1.0 for g in mva_in_lib}, n, index)
        sigs["track1_combined"] = combined_vote(ko, weights, n, index)
        sigs["track1_high_score_only"] = combined_vote(ko, weights_high, n, index)
        res = {k: test_signature(D, v, null_r) for k, v in sigs.items()}
        # null for combined signatures: random genes carrying the same weights
        for key, w in (("track1_combined", weights), ("track1_high_score_only", weights_high), ("mva_set", {g: 1.0 for g in mva_in_lib})):
            wl = [x for g, x in w.items() if g in ko]
            if not wl:
                continue
            cols = []
            for _ in range(N_NULL):
                gs = rng.choice(pool, size=len(wl), replace=False)
                cols.append(combined_vote(ko, {g: x for g, x in zip(gs, wl)}, n, index))
            res[key] = test_signature(D, sigs[key], reversal(D, np.stack(cols, axis=1)))
        store[n] = (D, res, null_r)
        result["analysis"][f"N={n if n else 'full'}"] = {k: {"n_usable": int(t["usable"].sum()), "n_hazard": int(t["hazard"].sum()),
                                                            "usable": sorted(names[i] for i in np.where(t["usable"])[0])} for k, t in res.items()}
    D, res, null_r = store[N_PRIMARY]
    geo = geo_signature(root_p / GEO_FILE, read_annotation(root_p / GPL_FILE))
    if geo["computed"]:
        p_sig = signed_vector(geo["up"], geo["down"], None, index)
        p_null = reversal(D, random_signature_null(geo["expressed_genes"], index, len(geo["up"]), len(geo["down"]), N_NULL, P_SEED))
        res["P"] = test_signature(D, p_sig, p_null)
        result["analysis"]["N=P(100+100)"] = {"P": {"n_usable": int(res["P"]["usable"].sum()), "n_hazard": int(res["P"]["hazard"].sum()),
                                                    "usable": sorted(names[i] for i in np.where(res["P"]["usable"])[0])}}
    # primary table
    support_keys = ["gene:BUB1B", "track1_combined", "P"]
    table = []
    for i, nm in enumerate(names):
        sup = [k for k, t in res.items() if t["usable"][i]]
        e_usable = any(k in sup for k in support_keys)
        d = dmeta[i]
        table.append({"library_name": nm, "name": d["name"], "drug_id": d["drug_id"], "atc": d.get("atc") or [], "e_usable": e_usable,
                      "supported_by_signatures": sorted(sup), "reversal_not_tied_to_mechanism": bool(sup) and not e_usable and not any(s.split(":")[-1] in mechanism or s == "mva_set" for s in sup),
                      "bub1b": _summ(res["gene:BUB1B"], names, i), "track1_combined": _summ(res["track1_combined"], names, i), "P": _summ(res["P"], names, i) if "P" in res else None,
                      "hazard_signatures": sorted(k for k, t in res.items() if t["hazard"][i]),
                      "fda_label": fda_label(d["name"]), "age_support": age_support(d["name"]), "label_level": LEVEL.get(d["name"].upper())})
    table_by_lib = {r["library_name"]: r for r in table}
    table.sort(key=lambda r: (-r["bub1b"]["R"], r["name"]))
    result["top20_reversal_bub1b"] = [{k: r[k] for k in ("name", "bub1b", "e_usable", "supported_by_signatures")} for r in table[:20]]
    result["hazards_bub1b_or_combined"] = [{"name": r["name"], "signatures": r["hazard_signatures"], "bub1b": r["bub1b"]} for r in table if r["hazard_signatures"]]
    result["line_e_usable_drugs"] = [r["name"] for r in table if r["e_usable"]]
    result["supported_by_track1_combined"] = [r["name"] for r in table if "track1_combined" in r["supported_by_signatures"]]
    result["not_tied_to_mechanism"] = [r["name"] for r in table if r["reversal_not_tied_to_mechanism"]]
    result["nine_report_drugs"] = {}
    by_lib = {n.lower(): i for i, n in enumerate(names)}
    for nm in NINE:
        key = NINE_LIB.get(nm, nm)
        i = by_lib.get(key)
        if i is None:
            result["nine_report_drugs"][nm] = {"line_e": "not matched / not in the chemical library", "usable": False}
            continue
        r = table[[x["library_name"] for x in table].index(names[i])]
        result["nine_report_drugs"][nm] = {"bub1b": r["bub1b"], "track1_combined": r["track1_combined"], "P": r["P"], "usable": r["e_usable"],
                                           "supported_by_signatures": r["supported_by_signatures"], "hazard_signatures": r["hazard_signatures"],
                                           "rank_by_reversal_bub1b": 1 + int((res["gene:BUB1B"]["R"] > res["gene:BUB1B"]["R"][i]).sum())}
    # candidates (strict rule; and the reading with B2 held exploratory)
    cand_rows = [{"name": r["name"], "e_usable": r["e_usable"], "atc": r["atc"], "fda_label": r["fda_label"], "age_support": r["age_support"], "label_level": r["label_level"]} for r in table]
    result["candidates_b2_prespecified"] = candidate_table(cand_rows, lines_usable, safety)
    result["candidates_b2_exploratory"] = candidate_table(cand_rows, {k: v for k, v in lines_usable.items() if k != "B2"}, safety)
    result["n_candidates"] = {"b2_prespecified": sum(c["candidate"] for c in result["candidates_b2_prespecified"]),
                              "b2_exploratory": sum(c["candidate"] for c in result["candidates_b2_exploratory"])}
    result["existing_lines_usable_sets"] = {k: sorted(v) for k, v in lines_usable.items()}
    result["existing_line_drugs_in_line_e"] = {n: ("matched" if n.lower() in by_lib else "no library entry") for n in sorted(set().union(*lines_usable.values()))}
    # sensitivity: agreement of usable sets
    sens = {}
    for nn, lab in ((50, "N=50"), (None, "N=full")):
        for k in ("gene:BUB1B", "track1_combined", "mva_set"):
            if k in store[nn][1] and k in res:
                a, b = set(np.where(res[k]["usable"])[0]), set(np.where(store[nn][1][k]["usable"])[0])
                sens[f"{k} {lab} vs N=150"] = {"n_150": len(a), "n_other": len(b), "overlap": len(a & b)}
    result["sensitivity"] = sens
    # controls
    ctrl = {}
    for g, drugs in CONTROL_PAIRS.items():
        if g not in ko:
            ctrl[g] = {"in_library": False}
            continue
        s = signed_vector(ko[g]["Up"], ko[g]["Down"], N_PRIMARY, index)
        t = test_signature(D, s, null_r)
        present = [x for x in drugs if x in by_lib]
        order_mim = np.argsort(t["R"], kind="stable")
        rank_mim = {names[i]: 1 + int(np.where(order_mim == i)[0][0]) for i in map(by_lib.get, present)}
        rank_rev = {names[i]: 1 + int((t["R"] > t["R"][i]).sum()) for i in map(by_lib.get, present)}
        n_d = len(names)
        best = min(rank_mim.values()) if rank_mim else None
        ctrl[g] = {"in_library": True, "drugs_present": present, "rank_by_mimicry": rank_mim, "rank_by_reversal": rank_rev, "n_drugs": n_d,
                   "best_mimicry_rank": best, "recovered_top5pct": bool(best is not None and best <= 0.05 * n_d),
                   "q_mimic": {names[i]: float(t["q_mim"][i]) for i in map(by_lib.get, present)}}
    result["positive_controls"] = {"primary": "MTOR", "results": ctrl}
    rng2 = np.random.default_rng(SEED + 1)
    neg = rng2.choice(pool, size=20, replace=False)
    counts = []
    for g in neg:
        s = signed_vector(ko[g]["Up"], ko[g]["Down"], N_PRIMARY, index)
        t = test_signature(D, s, null_r)
        counts.append(int(t["usable"].sum()))
    result["negative_control"] = {"n_pseudo_signatures": 20, "usable_drugs_per_signature": counts, "mean_fraction_usable": round(float(np.mean(counts)) / len(names), 4),
                                  "max": int(max(counts))}
    result["geo"] = {"accession": "GSE22206", "file_sha256": _sha(root_p / GEO_FILE), "platform_sha256": _sha(root_p / GPL_FILE),
                     **{k: v for k, v in geo.items() if k not in ("array_titles", "expressed_genes", "up", "down")}, "n_up": len(geo.get("up", [])), "n_down": len(geo.get("down", []))}
    if "P" in res:
        po = sorted(range(len(names)), key=lambda i: -res["P"]["R"][i])
        result["P_signature"] = {"top20_reversal": [{"name": table_by_lib[names[i]]["name"], **_summ(res["P"], names, i)} for i in po[:20]],
                                 "n_usable": int(res["P"]["usable"].sum()), "n_hazard": int(res["P"]["hazard"].sum()),
                                 "hazards": [table_by_lib[names[i]]["name"] for i in np.where(res["P"]["hazard"])[0]],
                                 "usable": [table_by_lib[names[i]]["name"] for i in np.where(res["P"]["usable"])[0]],
                                 "limitation": "one adult patient's lymphoblastoid cell line (n = 1 patient, 2 arrays) versus 2 controls; cannot represent the child"}
    result["table_all_drugs_n150"] = table
    (root_p / "out").mkdir(exist_ok=True)
    (root_p / "out/line_e_signature_reversal.json").write_text(json.dumps(result, indent=1, sort_keys=True, default=float), encoding="utf-8")
    (root_p / "out/line_e_signature_reversal.md").write_text(render_md(result), encoding="utf-8")


def render_md(r: dict[str, Any]) -> str:
    L = ["# Line E: signature reversal (fifth line)", "", f"Label: {r['label']}. Protocol: `{r['protocol']}`. The rule of at least two usable lines is unchanged; a candidate here is a hypothesis for follow-up only.", ""]
    dm, av = r["drug_matching"], r["signature_availability"]
    L += ["## Data and matching", f"* Knockout signatures in library: {r['data']['ko_library']['n_signatures']}; chemical signatures: {r['data']['chem_library']['n_signatures']}.",
          f"* Library drugs matched to approved small molecules of the project: {dm['n_matched_to_approved_small_molecules']} of {dm['n_library_drugs']} library names ({dm['n_project_small_molecules']} project small molecules).",
          f"* Curated MVA genes with a knockout signature: {', '.join(av['mva_in_library'])}; absent: {', '.join(av['mva_absent_from_library'])}.",
          f"* Track 1 genes: {av['n_track1_genes']}; with a knockout signature: {av['n_track1_in_library']} ({', '.join(av['track1_in_library'])}).",
          f"* GEO {r['geo']['accession']}: {r['geo'].get('reason', 'computed, see signature P')}.", "", "## Top 20 reversal drugs (BUB1B signature, N=150)", "",
          "| # | Drug | R | q (reversal) | usable | supported by |", "|---|---|---|---|---|---|"]
    for i, x in enumerate(r["top20_reversal_bub1b"], 1):
        b = x["bub1b"]
        L.append(f"| {i} | {x['name']} | {b['R']} | {b['q_rev']:.3g} | {x['e_usable']} | {', '.join(x['supported_by_signatures']) or '-'} |")
    if "P_signature" in r:
        ps = r["P_signature"]
        L += ["", f"## Signature P (patient-derived, GSE22206): {r['geo']['n_up']} Up + {r['geo']['n_down']} Down genes; {r['geo']['n_genes_mapped']} genes mapped; usable {ps['n_usable']}, hazards {ps['n_hazard']}", "",
              f"Limitation: {ps['limitation']}.", "", "| # | Drug | R | q (reversal) | usable |", "|---|---|---|---|---|"]
        for i, x in enumerate(ps["top20_reversal"][:20], 1):
            L.append(f"| {i} | {x['name']} | {x['R']} | {x['q_rev']:.3g} | {x['usable']} |")
    L += ["", "## The nine report drugs", "", "| Drug | BUB1B R | q_rev | Track1 combined R | q_rev | P R | P q_rev | line E usable | rank (BUB1B) |", "|---|---|---|---|---|---|---|---|---|"]
    for k, v in r["nine_report_drugs"].items():
        if "bub1b" not in v:
            L.append(f"| {k} | {v['line_e']} | | | | | | False | |")
        else:
            L.append(f"| {k} | {v['bub1b']['R']} | {v['bub1b']['q_rev']:.3g} | {v['track1_combined']['R']} | {v['track1_combined']['q_rev']:.3g} | {(v['P'] or {}).get('R', '')} | {('%.3g' % v['P']['q_rev']) if v.get('P') else ''} | {v['usable']} | {v['rank_by_reversal_bub1b']} |")
    L += ["", "## Controls", "", "Positive controls (inhibitors should mimic their gene knockout; rank by mimicry among matched drugs, top 5% = recovered):", ""]
    for g, c in r["positive_controls"]["results"].items():
        L.append(f"* {g}: " + ("not in library" if not c["in_library"] else f"drugs present {c['drugs_present']}; mimicry ranks {c['rank_by_mimicry']}; reversal ranks {c['rank_by_reversal']}; recovered (top 5%): {c['recovered_top5pct']}"))
    nc = r["negative_control"]
    L += [f"* Negative control: usable drugs per random pseudo-signature (20): mean fraction {nc['mean_fraction_usable']}, max {nc['max']} drugs.", "",
          "## Signatures and sensitivity", "", "| Setting | Signature | usable | hazards |", "|---|---|---|---|"]
    for lab, d in r["analysis"].items():
        for k, v in d.items():
            L.append(f"| {lab} | {k} | {v['n_usable']} | {v['n_hazard']} |")
    L += ["", "Agreement of usable sets: " + "; ".join(f"{k}: {v['overlap']} shared of {v['n_150']}/{v['n_other']}" for k, v in r["sensitivity"].items()), "",
          "## Candidates by the combination rule", f"Usable sets of the existing lines: {r['existing_lines_usable_sets']}.",
          f"Candidates (B2 as prespecified): {r['n_candidates']['b2_prespecified']}; with B2 held exploratory: {r['n_candidates']['b2_exploratory']}.",
          f"Drugs supported by line E (BUB1B and/or Track 1 combined): {len(r['line_e_usable_drugs'])}; by the Track 1 combined signature: {r['supported_by_track1_combined']}.", ""]
    if r["candidates_b2_prespecified"]:
        L += ["| Drug | existing lines | cytotoxic/caution | FDA label | candidate | age support (per band) |", "|---|---|---|---|---|---|"]
        for c in r["candidates_b2_prespecified"]:
            L.append(f"| {c['name']} | {', '.join(c['existing_lines']) or 'none'} | {c['oncology_cytotoxic_or_caution_class']} | {c['fda_label']} | {c['candidate']} | {c['age_support']['bands']} |")
    else:
        L.append("No drug is usable for line E, so the candidate table is empty. This negative result is reported as is.")
    L += ["", "## Hazards (concordant, mimic the knockout; FDR < 0.1)", ""]
    for h in r["hazards_bub1b_or_combined"][:30]:
        L.append(f"* {h['name']}: {', '.join(h['signatures'])} (BUB1B R = {h['bub1b']['R']})")
    L += ["", "Reversal not tied to the mechanism: " + (", ".join(r["not_tied_to_mechanism"]) or "none"), "",
          "## Limits", "A knockout signature is not the partial loss of function of the child; the consensus signatures come from cancer cell lines; the nulls use unrelated gene knockouts; hypotheses only, no dosing, no treatment recommendation."]
    return "\n".join(L) + "\n"
