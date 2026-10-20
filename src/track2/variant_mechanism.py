"""Protein-level mechanism characterisation of two BUB1B variants from public sources.

Variants (protein level only, MANE transcript NM_001211.6 / ENST00000287598, UniProt O60566):
  * p.Leu737*      nonsense: NMD rule (50-55 nt rule) from the Ensembl exon structure, and the
                   UniProt/InterPro features that the truncated product would lack.
  * p.Asn1002Lys   missense: Ensembl Compara conservation, AlphaFold model geometry
                   (pLDDT, solvent accessibility, neighbours, distance to known pathogenic
                   missense), ClinVar regional pathogenicity (+-20 residues, Fisher exact test),
                   percentile of the calibrated predictors among ClinVar BUB1B missense variants
                   (scores from Ensembl VEP REST), and the published MVA1 statement (Europe PMC).

No patient-level data are used or written: no genomic coordinates, cDNA positions, alleles,
genotypes or read metrics.  Public ClinVar records are reduced to protein change + class.
Every number is queried from a source at run time (responses cached under
data/cache_variant_mechanism/, access date stored with each cached response).  Method constants
(van der Waals radii, probe radius, Fisher window, NMD 50-55 nt rule, geometry thresholds for the
secondary-structure heuristic) are declared in ``METHOD_CONSTANTS`` and echoed in the output.

Run: ``python -m track2.variant_mechanism`` (PYTHONPATH=src; network needed on first run).
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = ROOT / "data" / "cache_variant_mechanism"
OUT_JSON = ROOT / "out" / "variant_mechanism.json"
CALIBRATION_JSON = ROOT / "out" / "empirical" / "variant_calibration.json"

UNIPROT = "O60566"
ENST = "ENST00000287598"
ENSG = "ENSG00000156970"
MANE_REFSEQ = "NM_001211.6"
STOP_POS = 737
MIS_POS = 1002
MIS_REF, MIS_ALT = "Asn", "Lys"
EXPECTED_STOP_AA = "L"
EXPECTED_MIS_AA = "N"
PMID_REQUESTED = "20530692"      # as given in the task; resolved by code (see literature)
PMID_SUIJKERBUIJK = "20516114"   # found by title/author search in Europe PMC

METHOD_CONSTANTS = {
    "nmd_rule_nt": [50, 55],
    "nmd_rule_text": "premature stop codon more than 50-55 nt upstream of the last exon-exon junction "
                     "is expected to trigger NMD (standard rule)",
    "sasa_probe_radius_A": 1.4,
    "sasa_sphere_points": 250,
    "vdw_radii_A": {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80},
    "neighbour_radius_A": 8.0,
    "contact_number_radius_A_cbeta": 10.0,
    "fisher_window_half_width": 20,
    "conservation_window_half_width": 10,
    "plddt_confident_min": 70.0,
    "buried_sasa_percentile_max": 25.0,
    "conserved_fraction_min": 0.75,
    "ss_heuristic": "CA-trace geometry: helix-like if d(i,i+3)<=5.8 A and d(i,i+4)<=6.8 A; "
                    "strand-like if d(i-1,i+1)>=6.2 A and d(i,i+2)>=6.2 A; else coil. "
                    "Heuristic (not DSSP); reported only where pLDDT>=70.",
    "relative_accessibility": "no tabulated maximum-ASA values are used; accessibility is expressed as the "
                              "percentile of the residue's SASA among all residues of the same amino acid "
                              "type in the same model",
}

_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
        "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
        "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
_AA1_TO_3 = {v: k for k, v in _AA3.items()}
_CODONS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_BASES = "TCAG"
CODON_TABLE = {a + b + c: _CODONS[16 * i + 4 * j + k]
               for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


# --------------------------------------------------------------------------- HTTP + cache
_LAST_CALL = [0.0]


def _today() -> str:
    return _dt.date.today().isoformat()


def cached_request(url: str, params: dict | None = None, json_body=None, headers: dict | None = None,
                   binary: bool = False, min_interval: float = 0.25, retries: int = 6) -> dict:
    """GET (or POST with json_body) with on-disk cache and exponential backoff.

    Returns {"url", "accessed", "text"}; "accessed" is the date of the original (cached) download.
    """
    full = url + ("?" + urllib.parse.urlencode(params, doseq=True) if params else "")
    key = hashlib.sha1((full + "|" + json.dumps(json_body, sort_keys=True)).encode()).hexdigest()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = CACHE_DIR / (key + ".json")
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    hdrs = {"User-Agent": "track2-variant-mechanism/1.0", "Accept": "application/json"}
    hdrs.update(headers or {})
    data = None
    if json_body is not None:
        data = json.dumps(json_body).encode()
        hdrs["Content-Type"] = "application/json"
    last = None
    for attempt in range(retries):
        wait = min_interval - (time.time() - _LAST_CALL[0])
        if wait > 0:
            time.sleep(wait)
        _LAST_CALL[0] = time.time()
        try:
            req = urllib.request.Request(full, data=data, headers=hdrs)
            with urllib.request.urlopen(req, timeout=120) as r:
                text = r.read().decode("utf-8", errors="replace")
            rec = {"url": full, "accessed": _today(), "text": text}
            path.write_text(json.dumps(rec), encoding="utf-8")
            return rec
        except urllib.error.HTTPError as e:
            last = e
            if e.code not in (429, 500, 502, 503, 504):
                raise
            time.sleep(min(2 ** attempt, 30))
        except (urllib.error.URLError, OSError) as e:  # includes SSL errors
            last = e
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"request failed after {retries} attempts: {url}: {last}")


def get_json(url, params=None, **kw):
    rec = cached_request(url, params, **kw)
    return json.loads(rec["text"]), rec["accessed"], rec["url"]


# --------------------------------------------------------------------------- generic helpers
def translate(cds: str) -> str:
    return "".join(CODON_TABLE.get(cds[i:i + 3].upper(), "X") for i in range(0, len(cds) - 2, 3))


def parse_protein_change(title: str):
    """'NM_001211.6(BUB1B):c.742G>T (p.Ala248Ser)' -> (248, 'A', 'S') or None (BUB1B MANE only)."""
    if not re.search(r"NM_001211\.\d+\(BUB1B\)", title or ""):
        return None
    m = re.search(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)", title)
    if not m or m.group(1) not in _AA3 or m.group(3) not in _AA3:
        return None
    return int(m.group(2)), _AA3[m.group(1)], _AA3[m.group(3)]


def classify_clinvar(description: str) -> str:
    d = (description or "").strip().lower()
    if d in ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"):
        return "P_LP"
    if d in ("benign", "likely benign", "benign/likely benign"):
        return "B_LB"
    return "other"


def fisher_window(variants: list[dict], centre: int, half: int) -> dict:
    """variants: [{'pos', 'cls'}] with cls in P_LP / B_LB (others ignored).  Two-sided Fisher exact test."""
    from scipy.stats import fisher_exact
    lo, hi = centre - half, centre + half
    a = b = c = d = 0  # a: P in window, b: B in window, c: P outside, d: B outside
    for v in variants:
        inside = lo <= v["pos"] <= hi
        if v["cls"] == "P_LP":
            a += inside
            c += not inside
        elif v["cls"] == "B_LB":
            b += inside
            d += not inside
    n_in, n_out = a + b, c + d
    res = {"window": [lo, hi], "in_window": {"pathogenic": a, "benign": b, "n": n_in},
           "outside_window": {"pathogenic": c, "benign": d, "n": n_out},
           "fraction_pathogenic_in_window": (a / n_in) if n_in else None,
           "fraction_pathogenic_outside": (c / n_out) if n_out else None}
    if n_in and n_out and (a + c) and (b + d):
        odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        res["fisher_odds_ratio"] = None if math.isinf(odds) else float(odds)
        res["fisher_p_two_sided"] = float(p)
    else:
        res["fisher_odds_ratio"] = None
        res["fisher_p_two_sided"] = None
        res["fisher_note"] = "test not computable: an empty row or column in the 2x2 table"
    return res


def fraction_at_least_as_damaging(values: list[float], x: float, higher_is_damaging: bool) -> dict:
    vals = [v for v in values if v is not None]
    if not vals:
        return {"n": 0, "fraction_at_least_as_damaging": None}
    k = sum(1 for v in vals if (v >= x if higher_is_damaging else v <= x))
    return {"n": len(vals), "n_at_least_as_damaging": k, "fraction_at_least_as_damaging": k / len(vals)}


def percentile_rank(values: list[float], x: float) -> float | None:
    """Percentage of values strictly below x plus half of ties (mid-rank)."""
    vals = list(values)
    if not vals:
        return None
    lt = sum(1 for v in vals if v < x)
    eq = sum(1 for v in vals if v == x)
    return 100.0 * (lt + 0.5 * eq) / len(vals)


# --------------------------------------------------------------------------- (1) NMD
def transcript_order_exons(exons: list[dict], strand: int) -> list[dict]:
    return sorted(exons, key=lambda e: e["start"], reverse=(strand < 0))


def exon_lengths(exons_tx: list[dict]) -> list[int]:
    return [e["end"] - e["start"] + 1 for e in exons_tx]


def genomic_to_tx(exons_tx: list[dict], g: int, strand: int) -> int | None:
    """1-based transcript position of genomic coordinate g (None if intronic)."""
    cum = 0
    for e in exons_tx:
        if e["start"] <= g <= e["end"]:
            return cum + ((g - e["start"] + 1) if strand > 0 else (e["end"] - g + 1))
        cum += e["end"] - e["start"] + 1
    return None


def nmd_assessment(lengths: list[int], cds_start_tx: int, codon: int) -> dict:
    """Standard 50-55 nt rule.  lengths: exon lengths in transcript order; cds_start_tx: 1-based tx
    position of CDS nt 1; codon: protein position of the premature stop."""
    first_nt = cds_start_tx + (codon - 1) * 3        # tx position of the first nt of the stop codon
    last_nt = first_nt + 2
    cum, bounds = 0, []
    for L in lengths:
        cum += L
        bounds.append(cum)                           # junction after tx position bounds[k]
    n = len(lengths)
    exon_no = next(i + 1 for i, b in enumerate(bounds) if first_nt <= b)
    split = next(i + 1 for i, b in enumerate(bounds) if last_nt <= b) != exon_no
    last_junction = bounds[-2] if n > 1 else None
    out = {"n_exons": n, "stop_exon_number": exon_no, "stop_in_last_exon": exon_no == n,
           "stop_codon_split_by_junction": split, "stop_exon_length_nt": lengths[exon_no - 1],
           "last_exon_length_nt": lengths[-1], "junctions_downstream_of_stop": n - exon_no,
           "rule_nt": list(METHOD_CONSTANTS["nmd_rule_nt"])}
    if last_junction is None or exon_no == n:
        out.update(distance_to_last_junction_nt=None if last_junction is None else last_junction - last_nt,
                   distance_class="NMD_escape_predicted: stop in the last exon (no downstream junction)",
                   nmd_predicted=False)
        return out
    dist = last_junction - last_nt                   # nt between end of stop codon and last junction
    lo, hi = METHOD_CONSTANTS["nmd_rule_nt"]
    if dist > hi:
        cls, pred = f"NMD_predicted: >{hi} nt upstream of the last exon-exon junction", True
    elif dist >= lo:
        cls, pred = f"borderline: {lo}-{hi} nt upstream of the last exon-exon junction", None
    else:
        cls, pred = f"NMD_escape_predicted: <{lo} nt upstream of the last exon-exon junction", False
    out.update(distance_to_last_junction_nt=dist, distance_class=cls, nmd_predicted=pred)
    return out


def lost_features(features: list[dict], first_lost: int, length: int) -> dict:
    """Split functional features (start,end,type,description) into lost / partly lost / retained
    for a protein truncated before residue ``first_lost``."""
    lost, partial, kept = [], [], []
    for f in features:
        s, e = f["start"], f["end"]
        item = {k: f[k] for k in ("type", "description", "start", "end", "source") if k in f}
        if s >= first_lost:
            lost.append(item)
        elif e >= first_lost:
            partial.append(item)
        else:
            kept.append(item)
    key = lambda x: (x["start"], x["end"], x["type"])
    return {"truncated_product": {"retained_residues": [1, first_lost - 1],
                                  "lost_residues": [first_lost, length],
                                  "n_lost": length - first_lost + 1,
                                  "fraction_lost": (length - first_lost + 1) / length},
            "lost_entirely": sorted(lost, key=key), "partially_lost": sorted(partial, key=key),
            "retained": sorted(kept, key=key)}


def scan_motif(seq: str, motif: str) -> list[int]:
    return [m.start() + 1 for m in re.finditer(f"(?={motif})", seq)]


# --------------------------------------------------------------------------- (2a) conservation
def window_identity(human_aln: str, other_aln: str, pos: int, half: int) -> dict:
    """Aligned residue at human position ``pos`` and identity over human residues pos-half..pos+half."""
    hp = 0
    col = None
    window_cols = []
    for i, ch in enumerate(human_aln):
        if ch not in "-.":
            hp += 1
            if hp == pos:
                col = i
            if pos - half <= hp <= pos + half:
                window_cols.append(i)
    if col is None:
        raise ValueError("position beyond the human sequence")
    res = other_aln[col]
    ident = sum(1 for i in window_cols if other_aln[i] == human_aln[i])
    gaps = sum(1 for i in window_cols if other_aln[i] in "-.")
    return {"aligned_residue": "-" if res in "-." else res,
            "identical_at_position": res == human_aln[col],
            "window_residues": len(window_cols), "window_identical": ident, "window_gaps": gaps,
            "window_identity_fraction": ident / len(window_cols) if window_cols else None}


# --------------------------------------------------------------------------- (2b) structure
def parse_pdb(text: str) -> list[dict]:
    atoms = []
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        el = line[76:78].strip() or line[12:16].strip()[0]
        atoms.append({"name": line[12:16].strip(), "res": line[17:20].strip(), "chain": line[21],
                      "resi": int(line[22:26]), "xyz": (float(line[30:38]), float(line[38:46]),
                                                        float(line[46:54])),
                      "b": float(line[60:66]), "el": el.upper()})
    return atoms


def shrake_rupley(coords, radii, probe: float, n_points: int):
    """Per-atom solvent accessible surface area (A^2), numpy + KD-tree implementation."""
    import numpy as np
    from scipy.spatial import cKDTree
    coords = np.asarray(coords, float)
    radii = np.asarray(radii, float)
    i = np.arange(n_points) + 0.5
    phi = math.pi * (3.0 - math.sqrt(5.0)) * np.arange(n_points)
    z = 1 - 2 * i / n_points
    r = np.sqrt(np.maximum(0.0, 1 - z * z))
    sphere = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)
    ext = radii + probe
    tree = cKDTree(coords)
    out = np.zeros(len(coords))
    maxr = ext.max()
    for a in range(len(coords)):
        nbrs = [j for j in tree.query_ball_point(coords[a], ext[a] + maxr) if j != a]
        pts = coords[a] + sphere * ext[a]
        free = np.ones(n_points, bool)
        for j in nbrs:
            if np.linalg.norm(coords[j] - coords[a]) < ext[a] + ext[j]:
                free &= np.sum((pts - coords[j]) ** 2, axis=1) >= ext[j] ** 2
                if not free.any():
                    break
        out[a] = 4.0 * math.pi * ext[a] ** 2 * free.sum() / n_points
    return out


def residue_table(atoms: list[dict]) -> dict:
    res: dict[int, dict] = {}
    for k, a in enumerate(atoms):
        r = res.setdefault(a["resi"], {"name": a["res"], "idx": [], "plddt": a["b"]})
        r["idx"].append(k)
    return res


def ss_heuristic(ca: dict, i: int) -> str:
    def d(a, b):
        if a not in ca or b not in ca:
            return None
        return math.dist(ca[a], ca[b])
    h3, h4 = d(i, i + 3), d(i, i + 4)
    h3b, h4b = d(i - 3, i), d(i - 4, i)
    if (h3 is not None and h4 is not None and h3 <= 5.8 and h4 <= 6.8) or \
       (h3b is not None and h4b is not None and h3b <= 5.8 and h4b <= 6.8):
        return "helix-like"
    s1, s2 = d(i - 1, i + 1), d(i, i + 2)
    if s1 is not None and s2 is not None and s1 >= 6.2 and s2 >= 6.2:
        return "strand-like"
    return "coil"


def structure_features(pdb_text: str, pos: int, expected_aa3: str, known_path: list[dict]) -> dict:
    import numpy as np
    atoms = parse_pdb(pdb_text)
    radii_tab = METHOD_CONSTANTS["vdw_radii_A"]
    heavy = [a for a in atoms if a["el"] in radii_tab]
    coords = np.array([a["xyz"] for a in heavy])
    radii = [radii_tab[a["el"]] for a in heavy]
    sasa = shrake_rupley(coords, radii, METHOD_CONSTANTS["sasa_probe_radius_A"],
                         METHOD_CONSTANTS["sasa_sphere_points"])
    res = residue_table(heavy)
    if pos not in res or res[pos]["name"].capitalize() != expected_aa3:
        raise ValueError(f"model residue {pos} is {res.get(pos, {}).get('name')}, expected {expected_aa3}")
    res_sasa = {r: float(sum(sasa[k] for k in v["idx"])) for r, v in res.items()}
    same = [res_sasa[r] for r, v in res.items() if v["name"] == res[pos]["name"]]
    plddt = {r: v["plddt"] for r, v in res.items()}
    # C-beta (CA for Gly) contact number
    cb = {}
    ca = {}
    for r, v in res.items():
        by = {heavy[k]["name"]: heavy[k]["xyz"] for k in v["idx"]}
        ca[r] = by.get("CA")
        cb[r] = by.get("CB", by.get("CA"))
    keys = [r for r in res if cb[r] is not None]
    cbarr = np.array([cb[r] for r in keys])
    from scipy.spatial import cKDTree
    ctree = cKDTree(cbarr)
    rad = METHOD_CONSTANTS["contact_number_radius_A_cbeta"]
    cn = {r: len(ctree.query_ball_point(cb[r], rad)) - 1 for r in keys}
    # residues within 8 A (any heavy atom)
    tree = cKDTree(coords)
    cut = METHOD_CONSTANTS["neighbour_radius_A"]
    near = set()
    for k in res[pos]["idx"]:
        for j in tree.query_ball_point(coords[k], cut):
            near.add(heavy[j]["resi"])
    near.discard(pos)
    near_sorted = sorted(near)

    def min_dist(r1, r2):
        a = coords[res[r1]["idx"]]
        b = coords[res[r2]["idx"]]
        return float(np.min(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)))

    pl = plddt[pos]
    win = [plddt[r] for r in range(pos - 10, pos + 11) if r in plddt]
    ca_ok = {r: v for r, v in ca.items() if v is not None}
    ss = ss_heuristic(ca_ok, pos)
    pathogenic = []
    for kp in known_path:
        p = kp["pos"]
        if p not in res:
            continue
        pathogenic.append({
            "variant": kp["label"], "source": kp["source"], "position": p,
            "sequence_distance_residues": abs(p - pos),
            "ca_ca_distance_A": round(math.dist(ca[pos], ca[p]), 2),
            "min_heavy_atom_distance_A": round(min_dist(pos, p), 2),
            "within_8A_any_atom": p in near, "plddt_at_site": round(plddt[p], 1)})
    pathogenic.sort(key=lambda x: (x["ca_ca_distance_A"], x["position"]))
    sasa_pct = percentile_rank(same, res_sasa[pos])
    cn_pct = percentile_rank(list(cn.values()), cn[pos])
    return {"plddt": round(pl, 2), "plddt_mean_pm10": round(sum(win) / len(win), 2),
            "plddt_confident": pl >= METHOD_CONSTANTS["plddt_confident_min"],
            "sasa_A2": round(res_sasa[pos], 2),
            "n_same_aa_residues_in_model": len(same),
            "sasa_percentile_among_same_aa": round(sasa_pct, 1),
            "sasa_same_aa_median_A2": float(np.median(same)),
            "sasa_same_aa_p90_A2": float(np.percentile(same, 90)),
            "sasa_relative_to_same_aa_p90": round(res_sasa[pos] / float(np.percentile(same, 90)), 3)
            if np.percentile(same, 90) > 0 else None,
            "cbeta_contact_number_10A": cn[pos], "contact_number_percentile": round(cn_pct, 1),
            "secondary_structure_heuristic": ss if pl >= METHOD_CONSTANTS["plddt_confident_min"] else
            "not assigned (pLDDT below confident threshold)",
            "secondary_structure_heuristic_raw": ss,
            "residues_within_8A": [f"{res[r]['name'].capitalize()}{r}" for r in near_sorted],
            "n_residues_within_8A": len(near_sorted),
            "known_pathogenic_missense_distances": pathogenic,
            "known_pathogenic_within_8A": [x["variant"] for x in pathogenic if x["within_8A_any_atom"]],
            "n_residues_modelled": len(res),
            "model_mean_plddt": round(sum(plddt.values()) / len(plddt), 2)}


# --------------------------------------------------------------------------- (2d) literature
def normalise_text(t: str) -> str:
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", t).strip()


def verify_sentences(abstract: str, sentences: dict[str, str]) -> dict:
    norm = normalise_text(abstract)
    return {k: {"text": s, "verified_verbatim_in_europepmc_abstract": normalise_text(s) in norm}
            for k, s in sentences.items()}


# --------------------------------------------------------------------------- mechanism statements
def nonsense_statement(nmd: dict | None, lost: dict | None, lit_ok: bool) -> dict:
    clauses = []
    if nmd:
        dist, exon, n = nmd.get("distance_to_last_junction_nt"), nmd["stop_exon_number"], nmd["n_exons"]
        if nmd.get("nmd_predicted") is True:
            clauses.append({"clause": f"predicted loss of function by nonsense-mediated decay (stop in exon "
                                      f"{exon} of {n}, {dist} nt upstream of the last exon-exon junction, "
                                      f"beyond the 50-55 nt rule)",
                            "feature": "nmd.distance_class", "value": nmd["distance_class"]})
        elif nmd.get("nmd_predicted") is None:
            clauses.append({"clause": f"NMD borderline (stop in exon {exon} of {n}, {dist} nt from the last "
                                      f"junction, inside the 50-55 nt window); decay not predictable by rule",
                            "feature": "nmd.distance_class", "value": nmd["distance_class"]})
        else:
            clauses.append({"clause": f"NMD escape predicted (stop in exon {exon} of {n}); a truncated "
                                      f"protein could be made", "feature": "nmd.distance_class",
                            "value": nmd["distance_class"]})
    else:
        clauses.append({"clause": "NMD status could not be assessed (Ensembl exon structure unavailable)",
                        "feature": "nmd", "value": None})
    if lost:
        tp = lost["truncated_product"]
        names = [f"{f['description']} ({f['start']}-{f['end']})" for f in lost["lost_entirely"]
                 + lost["partially_lost"] if f["type"] in ("Domain", "Active site", "Region", "Motif")
                 or f["source"].startswith("InterPro")]
        dom = [f for f in lost["lost_entirely"] + lost["partially_lost"]
               if f["type"] == "Domain" and "kinase" in (f["description"] or "").lower()]
        txt = (f"any truncated product would lack residues {tp['lost_residues'][0]}-{tp['lost_residues'][1]} "
               f"({100 * tp['fraction_lost']:.0f}% of the protein)")
        if dom:
            txt += f", including the kinase-like domain ({dom[0]['start']}-{dom[0]['end']}, UniProt)"
        clauses.append({"clause": txt, "feature": "lost_features", "value": names})
    if lit_ok:
        clauses.append({"clause": "published MVA1 analysis attributes low BubR1 abundance in truncating "
                                  "mutants to absence of transcripts", "feature": "literature.suijkerbuijk2010",
                        "value": "verified abstract sentence"})
    return {"statement": "; ".join(c["clause"] for c in clauses) + ".", "clauses": clauses,
            "what_would_test_it": [
                "RT-qPCR / allele-specific expression of BUB1B transcript in patient cells with and without "
                "an NMD inhibitor (cycloheximide or UPF1 knockdown): decay predicts recovery of the "
                "truncating-allele transcript on inhibition",
                "Western blot of BubR1 (C-terminal and N-terminal epitope antibodies) in patient cells versus "
                "controls: loss of the full-length band and absence of a ~736-residue product support decay "
                "without a stable truncated protein"]}


def missense_statement(feat: dict) -> dict:
    clauses = []
    cons, st, reg = feat.get("conservation"), feat.get("structure"), feat.get("clinvar_window")
    dom = feat.get("domain")
    if cons and cons.get("all_one2one_orthologs"):
        c = cons["all_one2one_orthologs"]
        frac = c["fraction_identical_at_position_of_aligned"]
        word = "conserved" if frac is not None and frac >= METHOD_CONSTANTS["conserved_fraction_min"] \
            else "not strongly conserved"
        clauses.append({"clause": f"missense at a {word} position (Asn in {c['n_identical']}/{c['n_aligned']} "
                                  f"aligned one-to-one orthologs; mean +-10 window identity "
                                  f"{c['mean_window_identity']:.2f})", "feature": "conservation",
                        "value": frac})
    else:
        clauses.append({"clause": "conservation could not be assessed (Ensembl Compara unavailable)",
                        "feature": "conservation", "value": None})
    if dom:
        clauses.append({"clause": f"lies in the protein kinase-like domain ({dom[0]}-{dom[1]}, UniProt)",
                        "feature": "uniprot.domain", "value": dom})
    stab = None
    if st:
        buried = st["sasa_percentile_among_same_aa"] <= METHOD_CONSTANTS["buried_sasa_percentile_max"]
        if not st["plddt_confident"]:
            clauses.append({"clause": f"AlphaFold local confidence is low (pLDDT {st['plddt']}), so burial "
                                      f"and stability cannot be inferred", "feature": "structure.plddt",
                            "value": st["plddt"]})
            stab = "not determinable (low-confidence model)"
        else:
            clauses.append({"clause": f"{'buried' if buried else 'not buried'} in the AlphaFold model "
                                      f"(pLDDT {st['plddt']}, SASA {st['sasa_A2']} A2 = "
                                      f"{st['sasa_percentile_among_same_aa']}th percentile of Asn residues; "
                                      f"{st['n_residues_within_8A']} residues within 8 A)",
                            "feature": "structure.sasa_percentile_among_same_aa",
                            "value": st["sasa_percentile_among_same_aa"]})
            conserved = cons and cons.get("all_one2one_orthologs") and \
                cons["all_one2one_orthologs"]["fraction_identical_at_position_of_aligned"] is not None and \
                cons["all_one2one_orthologs"]["fraction_identical_at_position_of_aligned"] >= \
                METHOD_CONSTANTS["conserved_fraction_min"]
            stab = "expected (buried + conserved + confident model)" if buried and conserved else \
                "not expected from structure alone (criteria for buried and conserved not both met)"
        near = st["known_pathogenic_missense_distances"]
        if near:
            n0 = near[0]
            clauses.append({"clause": f"nearest known pathogenic missense is {n0['variant']} "
                                      f"({n0['sequence_distance_residues']} residues away, "
                                      f"{n0['ca_ca_distance_A']} A CA-CA); "
                                      f"{len(st['known_pathogenic_within_8A'])} within 8 A",
                            "feature": "structure.known_pathogenic_missense_distances",
                            "value": n0["variant"]})
    clauses.append({"clause": f"predicted effect on protein stability/abundance: {stab or 'not determinable'}",
                    "feature": "structure+conservation", "value": stab})
    if reg and reg.get("fisher_p_two_sided") is not None:
        w = reg["in_window"]
        clauses.append({"clause": f"ClinVar missense in +-20 residues: {w['pathogenic']} P/LP vs "
                                  f"{w['benign']} B/LB (Fisher p={reg['fisher_p_two_sided']:.3g}); "
                                  f"regional signal "
                                  f"{'present' if reg['fisher_p_two_sided'] < 0.05 else 'not demonstrable'}",
                        "feature": "clinvar_window.fisher_p_two_sided", "value": reg["fisher_p_two_sided"]})
    elif reg:
        clauses.append({"clause": "ClinVar regional test not computable (too few classified missense)",
                        "feature": "clinvar_window", "value": None})
    return {"statement": "; ".join(c["clause"] for c in clauses) + ".", "clauses": clauses,
            "predicted_stability_effect": stab,
            "what_would_test_it": [
                "Western blot of BubR1 level in patient cells versus controls; in a compound-heterozygous "
                "state the truncating allele gives no protein, so residual BubR1 reflects the missense "
                "product",
                "cycloheximide chase of BubR1 to test increased turnover (the mechanism reported for MVA1 "
                "missense mutants)",
                "re-expression of Asn1002Lys versus wild-type BubR1 in BUB1B-depleted cells with a "
                "nocodazole mitotic-arrest (checkpoint) and chromosome alignment assay"]}


# --------------------------------------------------------------------------- data collection
def sec_uniprot(errors):
    try:
        d, acc, url = get_json(f"https://rest.uniprot.org/uniprotkb/{UNIPROT}.json")
    except Exception as e:  # noqa: BLE001
        errors["uniprot"] = str(e)
        return None
    seq = d["sequence"]["value"]
    feats = []
    for f in d.get("features", []):
        loc = f["location"]
        feats.append({"type": f["type"], "description": f.get("description", ""),
                      "start": loc["start"]["value"], "end": loc["end"]["value"],
                      "alt": f.get("alternativeSequence"), "source": "UniProt"})
    return {"sequence": seq, "features": feats, "accessed": acc, "url": url,
            "length": d["sequence"]["length"], "entry_version": d.get("entryAudit", {}).get("entryVersion")}


def uniprot_known_pathogenic(feats):
    out = []
    for f in feats:
        if f["type"] != "Natural variant" or f["start"] != f["end"] or not f.get("alt"):
            continue
        if not re.search(r"MVA1|PCS|mosaic variegated|premature chromatid", f["description"], re.I):
            continue
        alts = f["alt"].get("alternativeSequences") or []
        orig = f["alt"].get("originalSequence")
        if orig and len(orig) == 1 and alts and all(len(a) == 1 for a in alts):
            out.append({"pos": f["start"], "label": f"p.{_AA1_TO_3[orig]}{f['start']}{_AA1_TO_3[alts[0]]}",
                        "source": "UniProt natural variant: " + f["description"].split(";")[0]})
    return out


def sec_interpro(errors):
    try:
        d, acc, url = get_json(f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{UNIPROT}")
    except Exception as e:  # noqa: BLE001
        errors["interpro"] = str(e)
        return None
    feats = []
    for r in d["results"]:
        m = r["metadata"]
        if m["type"] not in ("domain", "repeat"):
            continue
        for loc in r["proteins"][0]["entry_protein_locations"]:
            for fr in loc["fragments"]:
                feats.append({"type": "Domain" if m["type"] == "domain" else "Repeat",
                              "description": f"{m['name']} [{m['accession']}, {m['source_database']}]",
                              "start": fr["start"], "end": fr["end"], "source": "InterPro"})
    return {"features": feats, "accessed": acc, "url": url}


def sec_ensembl_transcript(errors, uni_seq):
    try:
        tr, acc, url = get_json(f"https://rest.ensembl.org/lookup/id/{ENST}", {"expand": 1})
        xr, acc2, _ = get_json(f"https://rest.ensembl.org/xrefs/id/{ENST}", {"external_db": "RefSeq_mRNA"})
        cds, acc3, _ = get_json(f"https://rest.ensembl.org/sequence/id/{ENST}", {"type": "cds"})
    except Exception as e:  # noqa: BLE001
        errors["ensembl_transcript"] = str(e)
        return None
    strand = tr["strand"]
    exons_tx = transcript_order_exons(tr["Exon"], strand)
    lengths = exon_lengths(exons_tx)
    tl = tr["Translation"]
    g_cds_start = tl["start"] if strand > 0 else tl["end"]
    g_cds_end = tl["end"] if strand > 0 else tl["start"]
    cs = genomic_to_tx(exons_tx, g_cds_start, strand)
    ce = genomic_to_tx(exons_tx, g_cds_end, strand)
    cds_seq = cds["seq"]
    prot = translate(cds_seq)
    nmd = nmd_assessment(lengths, cs, STOP_POS)
    checks = {"cds_length_from_exons_nt": ce - cs + 1 if cs and ce else None,
              "cds_length_from_sequence_nt": len(cds_seq),
              "translation_matches_uniprot": prot.rstrip("*") == uni_seq,
              "residue_737_is_Leu": prot[STOP_POS - 1] == EXPECTED_STOP_AA,
              "residue_1002_is_Asn": prot[MIS_POS - 1] == EXPECTED_MIS_AA,
              "refseq_xref": [x.get("display_id") for x in xr],
              "mane_refseq_confirmed": MANE_REFSEQ in [x.get("display_id") for x in xr],
              "ensembl_transcript_version": tr.get("version"), "n_exons": len(lengths)}
    # the stop codon created by a 737 truncation sits at CDS nt (737-1)*3+1..+3; the codon must be intact
    return {"nmd": nmd, "checks": checks, "accessed": acc,
            "url": f"https://rest.ensembl.org/lookup/id/{ENST}?expand=1", "mane_xref_accessed": acc2,
            "cds_accessed": acc3}


def sec_conservation(errors, uni_seq):
    try:
        d, acc, url = get_json(f"https://rest.ensembl.org/homology/id/human/{ENSG}",
                               {"type": "orthologues", "format": "full", "sequence": "protein", "aligned": 1})
    except Exception as e:  # noqa: BLE001
        errors["ensembl_compara"] = str(e)
        return None
    homs = d["data"][0]["homologies"]
    half = METHOD_CONSTANTS["conservation_window_half_width"]
    mismatched = sum(1 for h in homs if h["source"]["align_seq"].replace("-", "") != uni_seq)
    if mismatched:
        errors["ensembl_compara_sequence"] = (f"{mismatched} pairwise alignments have a human sequence "
                                              f"different from UniProt (excluded)")
    named_species = ["mus_musculus", "rattus_norvegicus", "gallus_gallus", "xenopus_tropicalis", "danio_rerio",
                     "drosophila_melanogaster", "caenorhabditis_elegans", "saccharomyces_cerevisiae",
                     "pan_troglodytes", "bos_taurus"]
    rows, all_rows = [], []
    for h in homs:
        t = h["target"]
        src = h["source"]["align_seq"]          # each homology has its own pairwise alignment
        if src.replace("-", "") != uni_seq:
            continue
        w = window_identity(src, t["align_seq"], MIS_POS, half)
        row = {"species": t["species"], "target_id": t["id"], "homology_type": h["type"],
               "percent_identity_whole_protein": round(t["perc_id"], 1), **w}
        all_rows.append(row)
    for sp in named_species:
        for row in all_rows:
            if row["species"] == sp:
                rows.append(dict(row))
    for row in rows:  # gene symbols for the named subset (informative for Mad3 / BUB1 paralogues)
        try:
            ld, _, _ = get_json(f"https://rest.ensembl.org/lookup/id/{row['target_id']}")
            row["gene_symbol"] = ld.get("display_name")
        except Exception:  # noqa: BLE001
            row["gene_symbol"] = None
    o2o = [r for r in all_rows if r["homology_type"] == "ortholog_one2one"]
    aligned = [r for r in o2o if r["aligned_residue"] != "-"]
    from collections import Counter
    cnt = Counter(r["aligned_residue"] for r in aligned)
    summ = {"n_one2one_orthologs": len(o2o), "n_aligned": len(aligned),
            "n_identical": sum(1 for r in aligned if r["identical_at_position"]),
            "fraction_identical_at_position_of_aligned":
                (sum(1 for r in aligned if r["identical_at_position"]) / len(aligned)) if aligned else None,
            "residue_counts_at_position": dict(sorted(cnt.items())),
            "mean_window_identity": sum(r["window_identity_fraction"] for r in aligned) / len(aligned)
            if aligned else None}
    return {"named_orthologs": rows, "all_one2one_orthologs": summ, "n_homologies_total": len(homs),
            "window": [MIS_POS - half, MIS_POS + half], "accessed": acc, "url": url}


def sec_structure(errors, known_path, uni_seq):
    try:
        meta, acc, url = get_json(f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT}")
        m = meta[0]
        rec = cached_request(m["pdbUrl"], headers={"Accept": "*/*"})
    except Exception as e:  # noqa: BLE001
        errors["alphafold"] = str(e)
        return None
    atoms = parse_pdb(rec["text"])
    seq = "".join(_AA3.get(r["name"].capitalize(), "X") for r in
                  [dict(name=v["name"]) for _, v in sorted(residue_table(atoms).items())])
    feat = structure_features(rec["text"], MIS_POS, MIS_REF, known_path)
    feat.update({"model_id": m["modelEntityId"], "model_version": m["latestVersion"],
                 "model_file": m["pdbUrl"].rsplit("/", 1)[-1], "model_sequence_matches_uniprot": seq == uni_seq,
                 "accessed": rec["accessed"], "metadata_accessed": acc, "url": m["pdbUrl"]})
    return feat


def sec_clinvar(errors):
    term = 'BUB1B[gene] AND "missense variant"[molecular consequence]'
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    try:
        ids, start, count = [], 0, None
        while count is None or start < count:
            d, acc, _ = get_json(base + "esearch.fcgi", {"db": "clinvar", "term": term, "retmax": 1000,
                                                         "retstart": start, "retmode": "json"},
                                 min_interval=0.4)
            r = d["esearchresult"]
            count = int(r["count"])
            ids += r["idlist"]
            start += 1000
        summ = []
        for i in range(0, len(ids), 200):
            rec = cached_request(base + "esummary.fcgi", {"db": "clinvar", "id": ",".join(ids[i:i + 200]),
                                                          "retmode": "json"}, min_interval=0.4)
            res = json.loads(rec["text"])["result"]
            summ += [res[u] for u in res["uids"]]
    except Exception as e:  # noqa: BLE001
        errors["clinvar"] = str(e)
        return None
    recs = []
    for v in summ:
        pc = parse_protein_change(v.get("title", ""))
        if not pc:
            continue
        spdi = (v.get("variation_set") or [{}])[0].get("canonical_spdi", "")
        recs.append({"pos": pc[0], "ref": pc[1], "alt": pc[2],
                     "cls": classify_clinvar(v.get("germline_classification", {}).get("description", "")),
                     "description": v.get("germline_classification", {}).get("description", ""),
                     "review": v.get("germline_classification", {}).get("review_status", ""),
                     "_spdi": spdi})  # coordinates are used only transiently for VEP, never written out
    return {"records": recs, "n_esearch": count, "accessed": acc}


def sec_vep_scores(errors, recs):
    """Annotate ClinVar missense records with SIFT/PolyPhen/CADD/REVEL/AlphaMissense (Ensembl VEP REST)."""
    todo = []
    for r in recs:
        p = r["_spdi"].split(":")
        if len(p) == 4 and p[2] and p[3]:
            todo.append((r, f"{p[0].split('.')[0][-2:].lstrip('0')} {int(p[1]) + 1} . {p[2]} {p[3]} . . ."))
    accessed = None
    try:
        for i in range(0, len(todo), 100):
            chunk = todo[i:i + 100]
            rec = cached_request("https://rest.ensembl.org/vep/human/region?" + urllib.parse.urlencode(
                dict(AlphaMissense=1, CADD=1, REVEL=1, SIFT="b", PolyPhen="b", mane=1, hgvs=1)),
                json_body={"variants": [c[1] for c in chunk]}, min_interval=0.5)
            accessed = rec["accessed"]
            out = {x["input"]: x for x in json.loads(rec["text"])}
            for r, inp in chunk:
                x = out.get(inp)
                if not x:
                    continue
                for t in x.get("transcript_consequences", []):
                    if t.get("transcript_id") == ENST and t.get("protein_start") == r["pos"]:
                        am = t.get("alphamissense") or {}
                        r["scores"] = {"CADD": t.get("cadd_phred"), "SIFT": t.get("sift_score"),
                                       "PolyPhen2": t.get("polyphen_score"), "REVEL": t.get("revel"),
                                       "AlphaMissense": am.get("am_pathogenicity")}
                        break
    except Exception as e:  # noqa: BLE001
        errors["ensembl_vep"] = str(e)
        return None
    return {"accessed": accessed, "url": "https://rest.ensembl.org/vep/human/region"}


def sec_literature(errors):
    out = {"pmid_requested_in_task": PMID_REQUESTED, "pmid_used": PMID_SUIJKERBUIJK}
    try:
        base = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        chk, acc, _ = get_json(base, {"query": f"EXT_ID:{PMID_REQUESTED} AND SRC:MED", "resultType": "core",
                                      "format": "json"})
        r0 = chk["resultList"]["result"]
        out["requested_pmid_resolves_to"] = r0[0]["title"] if r0 else None
        d, acc, url = get_json(base, {"query": f"EXT_ID:{PMID_SUIJKERBUIJK} AND SRC:MED",
                                      "resultType": "core", "format": "json"})
        r = d["resultList"]["result"][0]
    except Exception as e:  # noqa: BLE001
        errors["europepmc"] = str(e)
        return None
    sentences = {
        "abundance_mechanism":
            "Low protein abundance is the direct result of the absence of transcripts from truncating "
            "mutants combined with high protein turnover of missense mutants.",
        "kinase_domain_location":
            "In this group of missense mutants, the amino acid change consistently occurs in or near the "
            "BUBR1 kinase domain."}
    out.update({"title": r.get("title"), "authors": r.get("authorString"), "journal": r.get("journalTitle"),
                "year": r.get("pubYear"), "pmcid": r.get("pmcid"), "accessed": acc, "url": url,
                "title_matches_expected": "BUBR1 dysfunction" in (r.get("title") or ""),
                "sentences": verify_sentences(r.get("abstractText", ""), sentences)})
    out["quotes_verified_by_code"] = all(v["verified_verbatim_in_europepmc_abstract"]
                                         for v in out["sentences"].values())
    return out


# --------------------------------------------------------------------------- main
def main() -> int:
    errors: dict = {}
    result: dict = {"generated": _today(), "gene": "BUB1B", "transcript": f"{MANE_REFSEQ} / {ENST}",
                    "uniprot": UNIPROT, "scope": "protein-level; no patient-level data used",
                    "method_constants": METHOD_CONSTANTS}
    uni = sec_uniprot(errors)
    seq = uni["sequence"] if uni else ""
    ipr = sec_interpro(errors)
    ens = sec_ensembl_transcript(errors, seq) if uni else None
    lit = sec_literature(errors)

    # ---- variant 1
    v1: dict = {"variant": f"p.Leu{STOP_POS}* ({MANE_REFSEQ})"}
    if ens:
        v1["nmd"] = {**ens["nmd"], "source": "Ensembl REST lookup/id (exon structure) + sequence/id (CDS)",
                     "accessed": ens["accessed"], "url": ens["url"], "consistency_checks": ens["checks"]}
    lost = None
    if uni:
        func_types = ("Domain", "Region", "Motif", "Active site", "Binding site", "Site", "Chain")
        feats = [f for f in uni["features"] if f["type"] in func_types and f["type"] != "Chain"]
        if ipr:
            feats += ipr["features"]
        mods = [f for f in uni["features"] if f["type"] in ("Modified residue", "Cross-link")]
        lost = lost_features(feats, STOP_POS, uni["length"])
        lost["lost_modified_residues"] = sorted(f["start"] for f in mods if f["start"] >= STOP_POS)
        lost["natural_variants_in_lost_region"] = sorted(
            {f["description"].split(";")[0] + f" @{f['start']}" for f in uni["features"]
             if f["type"] == "Natural variant" and f["start"] >= STOP_POS})
        lost["motif_scans_of_sequence"] = {
            "note": "KEN, Gle2-binding sequence (GLEBS) and ABBA motifs are not annotated as features in the "
                    "queried UniProt/InterPro records; KEN tripeptide occurrences below come from a plain "
                    "sequence scan (positions of K), not from a functional annotation",
            "KEN_tripeptide_start_positions": scan_motif(seq, "KEN"),
            "KEN_positions_lost_by_truncation": [p for p in scan_motif(seq, "KEN") if p >= STOP_POS]}
        lost["sources"] = {"uniprot": {"accessed": uni["accessed"], "url": uni["url"]},
                           "interpro": {"accessed": ipr["accessed"], "url": ipr["url"]} if ipr else
                           "not obtained"}
        v1["truncated_product_features"] = lost
    v1["mechanism"] = nonsense_statement(v1.get("nmd"), lost,
                                         bool(lit and lit.get("quotes_verified_by_code")))
    result["variant_Leu737Ter"] = v1

    # ---- variant 2
    known = uniprot_known_pathogenic(uni["features"]) if uni else []
    v2: dict = {"variant": f"p.{MIS_REF}{MIS_POS}{MIS_ALT} ({MANE_REFSEQ})"}
    if uni:
        dom = [f for f in uni["features"] if f["type"] == "Domain" and "kinase" in f["description"].lower()
               and f["start"] <= MIS_POS <= f["end"]]
        v2["domain"] = [dom[0]["start"], dom[0]["end"]] if dom else None
        v2["domain_source"] = {"description": dom[0]["description"] if dom else None, "source": "UniProt",
                               "accessed": uni["accessed"]}
        v2["reference_residue_is_Asn_in_UniProt"] = seq[MIS_POS - 1] == EXPECTED_MIS_AA
        v2["uniprot_secondary_structure_at_position"] = [
            f["type"] for f in uni["features"] if f["type"] in ("Helix", "Beta strand", "Turn")
            and f["start"] <= MIS_POS <= f["end"]]
    v2["conservation"] = sec_conservation(errors, seq) if uni else None
    clin = sec_clinvar(errors)
    v2["known_pathogenic_missense_sources"] = {"uniprot_MVA1_PCS": known}
    clin_pl = []
    if clin:
        clin_pl = [{"pos": r["pos"], "label": f"p.{_AA1_TO_3[r['ref']]}{r['pos']}{_AA1_TO_3[r['alt']]}",
                    "source": f"ClinVar {r['description']} ({r['review']})"} for r in clin["records"]
                   if r["cls"] == "P_LP"]
        v2["known_pathogenic_missense_sources"]["clinvar_P_LP"] = clin_pl
    allp, seen = [], set()
    for k in known + clin_pl:
        if k["label"] not in seen:
            seen.add(k["label"])
            allp.append(k)
    v2["structure"] = sec_structure(errors, allp, seq) if uni else None

    if clin:
        recs = clin["records"]
        same = [r for r in recs if r["pos"] == MIS_POS and r["ref"] == "N" and r["alt"] == "K"]
        v2["clinvar_variant_record_present"] = [{"classification": r["description"], "review_status": r["review"]}
                                                for r in same]
        used = [r for r in recs if r not in same]
        from collections import Counter
        cls_counts = Counter(r["description"] or "(none)" for r in recs)
        cw = fisher_window(used, MIS_POS, METHOD_CONSTANTS["fisher_window_half_width"])
        cw.update({"source": "NCBI E-utilities ClinVar (esearch+esummary), BUB1B missense on NM_001211.x",
                   "accessed": clin["accessed"], "n_missense_records_parsed": len(recs),
                   "classification_counts_all_missense": dict(sorted(cls_counts.items())),
                   "classes_used": "P_LP = Pathogenic/Likely pathogenic/Pathogenic-LP; B_LB = Benign/Likely "
                                   "benign; VUS, conflicting and others excluded; the variant itself excluded",
                   "p_lp_protein_changes": sorted(k["label"] for k in clin_pl),
                   "n_other_classes_in_window": sum(1 for r in used if r["cls"] == "other"
                                                    and abs(r["pos"] - MIS_POS)
                                                    <= METHOD_CONSTANTS["fisher_window_half_width"])})
        # sensitivity: add UniProt MVA1/PCS missense to the pathogenic class
        extra = [{"pos": k["pos"], "cls": "P_LP"} for k in known
                 if k["label"] not in {c["label"] for c in clin_pl}]
        cw["sensitivity_adding_uniprot_MVA1_PCS_as_pathogenic"] = fisher_window(
            [r for r in used if r["cls"] in ("P_LP", "B_LB")] + extra, MIS_POS,
            METHOD_CONSTANTS["fisher_window_half_width"])
        cw["interpretation_caveat"] = ("ClinVar BUB1B missense records are dominated by VUS; the number of "
                                       "classified variants is small, so the test has little power and "
                                       "ascertainment of classified variants is not random")
        v2["clinvar_window"] = cw
        # predictor percentiles
        vep = sec_vep_scores(errors, recs)
        try:
            cal = json.loads(CALIBRATION_JSON.read_text(encoding="utf-8"))
            patient_scores = cal["inputs"]["scores"]
        except Exception as e:  # noqa: BLE001
            errors["variant_calibration_json"] = str(e)
            patient_scores = {}
        if vep and patient_scores:
            pct = {}
            for name, x in patient_scores.items():
                hi = name != "SIFT"
                allv = [r["scores"][name] for r in used if r.get("scores") and r["scores"].get(name) is not None]
                bl = [r["scores"][name] for r in used if r["cls"] == "B_LB" and r.get("scores")
                      and r["scores"].get(name) is not None]
                pct[name] = {"variant_score": x, "direction": "higher = more damaging" if hi
                             else "lower = more damaging",
                             "among_all_clinvar_missense": {
                                 **fraction_at_least_as_damaging(allv, x, hi),
                                 "percentile_rank_of_variant_score": None if not allv else round(
                                     percentile_rank(allv, x) if hi else 100 - percentile_rank(allv, x), 1)},
                             "among_benign_likely_benign": {
                                 **fraction_at_least_as_damaging(bl, x, hi),
                                 "percentile_rank_of_variant_score": None if not bl else round(
                                     percentile_rank(bl, x) if hi else 100 - percentile_rank(bl, x), 1)},
                             "among_pathogenic_likely_pathogenic": [
                                 r["scores"][name] for r in used if r["cls"] == "P_LP" and r.get("scores")
                                 and r["scores"].get(name) is not None]}
            v2["predictor_percentiles"] = {
                "predictor_scores_source": "out/empirical/variant_calibration.json (inputs.scores)",
                "clinvar_scores_source": "Ensembl VEP REST (SIFT, PolyPhen-2, CADD PHRED, REVEL, AlphaMissense "
                                         "for ClinVar missense SNVs on ENST00000287598)",
                "accessed": vep["accessed"], "by_predictor": pct,
                "n_clinvar_missense_with_any_score": sum(1 for r in used if r.get("scores")),
                "percentile_definition": "percentile_rank = percentage of ClinVar missense variants with a "
                                         "less damaging score (mid-rank for ties)"}
        else:
            v2["predictor_percentiles"] = "could not be obtained (see errors)"
    else:
        v2["clinvar_window"] = "could not be obtained (see errors)"
        v2["predictor_percentiles"] = "could not be obtained (ClinVar unavailable)"
    v2["literature"] = lit
    v2["mechanism"] = missense_statement({"conservation": v2["conservation"], "structure": v2["structure"],
                                          "clinvar_window": v2["clinvar_window"]
                                          if isinstance(v2["clinvar_window"], dict) else None,
                                          "domain": v2.get("domain")})
    result["variant_Asn1002Lys"] = v2
    result["literature"] = lit
    result["errors_sources_not_obtained"] = errors
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {OUT_JSON}; sources not obtained: {sorted(errors) or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
