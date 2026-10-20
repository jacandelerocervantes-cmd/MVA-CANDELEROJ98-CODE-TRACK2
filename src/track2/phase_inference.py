"""Phase (cis versus trans) of the two BUB1B alleles as explicit Bayesian inference (limit B of the differential).

Run:  PYTHONPATH=src python -m track2.phase_inference      (offline; writes out/phase_inference.json)

This module quantifies a question; it does not assert an answer. No genomic coordinate, cDNA position, allele,
genotype or read metric of the patient is used: the model needs only the fact that two rare BUB1B alleles are present,
the gene-level verdict of out/differential_diagnosis.json, one public gene length, and one published incidence.

MODEL
-----
Hypotheses: H_trans (the two alleles lie on different chromosomes) and H_cis (same chromosome).

Prior. If both rare alleles were inherited, each came from the mother or the father with probability 1/2,
independently (independent transmission of two variants, recombination between them neglected). They lie on different
parental chromosomes (trans) when they come from different parents: P = 1/2; from the same parent they travel on the
same transmitted haplotype (cis): P = 1/2. If exactly one allele is de novo, it arose on the maternal or the paternal
chromosome; it is trans to the inherited allele iff it arose on the other parent's chromosome, so
P(trans | one de novo) = q*(1/2) + (1-q)*(1/2) = 1/2 for any paternal-origin fraction q, because the inherited allele
is equally likely to be maternal or paternal. Hence prior P(H_trans) = 1/2 whatever the de novo probability pi_dn
(`prior_trans`). A grid of other priors is reported as a robustness check.

Likelihood of "affected" (recessive disorder). Let c = P(BUB1B is the causal gene) and let "W0" be the alternative world
in which another cause explains the phenotype (the phenotype is then certain and phase-independent, b = 1).
    P(affected | trans) = c*f + (1-c)*1
    P(affected | cis)   = c*a + (1-c)*1
f = penetrance of the biallelic genotype (sensitivity range, an assumption: no sourced value is claimed).
a = probability that a cis carrier (effectively a single-allele carrier at the gene) is nevertheless affected through
another cause; it is bounded below by the population incidence of PCS/MVA syndrome, 1:1,000,000, read from the
Europe PMC full text PMC11102876 (PMID 38770464); the grid is widened up to 1e-3 to allow for undiagnosed other causes.
Posterior: P(trans | affected) = p*Lt / (p*Lt + (1-p)*Lc), with p the prior of H_trans. With c = 1 this reduces to
f / (f + a) at p = 1/2. The posterior is conditional on BUB1B being the causal gene; `posterior_trans` with c < 1 shows
the dilution.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DIFF = ROOT / "out" / "differential_diagnosis.json"
OUT = ROOT / "out" / "phase_inference.json"

# ---- sourced constants
MVA_INCIDENCE = 1e-6
MVA_INCIDENCE_SOURCE = {
    "pmid": "38770464",
    "pmcid": "PMC11102876",
    "clause": "the incidence of PCS/MVA syndrome is only 1:1,000,000",
    "retrieved_from": "Europe PMC full text REST (fullTextXML), searched 2026-10-15",
}
BUB1B_GENE_LENGTH_BP = 60154  # Ensembl REST lookup/symbol/homo_sapiens/BUB1B, release 116, GRCh38, ENSG00000156970 (end-start+1)
ENSEMBL_SOURCE = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/BUB1B (release 116)"
READ_LENGTHS_BP = (10_000, 15_000, 20_000, 25_000)  # standard long-read range stated in the task (10-25 kb)

# ---- assumption grids (not sourced; documented as sensitivity ranges)
PENETRANCE_GRID = (0.5, 0.7, 0.9, 1.0)
ALT_CAUSE_GRID = (1e-6, 1e-5, 1e-4, 1e-3)
CAUSAL_GRID = (1.0, 0.99, 0.9, 0.75, 0.5, 0.25, 0.1, 0.0)
PRIOR_GRID = (0.3, 0.5, 0.7)
HEADLINE = {"penetrance": 0.9, "alt_cause": MVA_INCIDENCE, "p_causal": 1.0, "prior": 0.5}


# ---------------------------------------------------------------- gene-level fact
def bub1b_sole_causal(diff: dict) -> bool:
    """True iff BUB1B is the only gene with verdict 'causal candidate' in the differential diagnosis."""
    genes = diff["genes"]
    causal = sorted(g for g, r in genes.items() if r.get("verdict") == "causal candidate")
    return causal == ["BUB1B"]


# ---------------------------------------------------------------- prior and posterior
def prior_trans(pi_dn: float = 0.0, q_paternal: float = 0.5) -> float:
    """P(H_trans) before any data. Both inherited: 1/2. One de novo: q*(1/2)+(1-q)*(1/2) = 1/2. Mixture is 1/2."""
    if not (0 <= pi_dn <= 1 and 0 <= q_paternal <= 1):
        raise ValueError("probabilities must lie in [0, 1]")
    inherited = 0.5
    one_de_novo = q_paternal * 0.5 + (1 - q_paternal) * 0.5
    return (1 - pi_dn) * inherited + pi_dn * one_de_novo


def likelihoods(f: float, a: float, c: float = 1.0) -> tuple:
    """(P(affected | trans), P(affected | cis)); see module docstring."""
    for x in (f, a, c):
        if not 0 <= x <= 1:
            raise ValueError("probabilities must lie in [0, 1]")
    return c * f + (1 - c), c * a + (1 - c)


def posterior_trans(f: float, a: float, c: float = 1.0, prior: float = 0.5) -> float:
    lt, lc = likelihoods(f, a, c)
    num = prior * lt
    den = num + (1 - prior) * lc
    return num / den


def sensitivity_table(c: float = 1.0, prior: float = 0.5) -> list:
    return [
        {"penetrance": f, "alt_cause": a, "posterior_trans": posterior_trans(f, a, c, prior)}
        for f in PENETRANCE_GRID
        for a in ALT_CAUSE_GRID
    ]


def causal_table(f: float, a: float, prior: float = 0.5) -> list:
    return [{"p_causal": c, "posterior_trans": posterior_trans(f, a, c, prior)} for c in CAUSAL_GRID]


def causal_threshold(target: float, f: float, a: float, prior: float = 0.5) -> float | None:
    """Smallest P(BUB1B causal) at which the posterior of trans is still >= target (bisection; posterior rises with c)."""
    if posterior_trans(f, a, 1.0, prior) < target:
        return None
    lo, hi = 0.0, 1.0
    if posterior_trans(f, a, 0.0, prior) >= target:
        return 0.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if posterior_trans(f, a, mid, prior) >= target:
            hi = mid
        else:
            lo = mid
    return hi


# ---------------------------------------------------------------- parental Sanger test
PATTERNS = {
    "T": "variant 1 in one parent only and variant 2 in the other parent only",
    "C": "both variants in one parent, neither in the other",
    "DN_M": "exactly one of the two variants found, in the mother only (the other variant absent in both parents: de novo)",
    "DN_F": "exactly one of the two variants found, in the father only (the other variant absent in both parents: de novo)",
}


def pattern_likelihoods(pi_dn: float = 0.0, q_paternal: float = 0.5, eps: float = 0.0) -> dict:
    """P(pattern | hypothesis). eps = recombination between the variants, mosaicism or genotyping error (assumption).

    Both inherited (1-pi_dn): trans -> T (1-eps) / C (eps); cis -> C (1-eps) / T (eps).
    One de novo (pi_dn): inherited allele maternal or paternal with 1/2 each; trans iff the de novo allele arose on the
    other parent's chromosome (paternal with probability q): P(DN_M | trans) = q, P(DN_F | trans) = 1-q, mirrored for cis.
    Both de novo (probability pi_dn squared) is neglected.
    """
    ok = 1 - pi_dn
    return {
        "trans": {"T": ok * (1 - eps), "C": ok * eps, "DN_M": pi_dn * q_paternal, "DN_F": pi_dn * (1 - q_paternal)},
        "cis": {"T": ok * eps, "C": ok * (1 - eps), "DN_M": pi_dn * (1 - q_paternal), "DN_F": pi_dn * q_paternal},
    }


def posterior_after_parental_test(pattern: str, prior: float = 0.5, pi_dn: float = 0.0, q_paternal: float = 0.5, eps: float = 0.0) -> float:
    pl = pattern_likelihoods(pi_dn, q_paternal, eps)
    num = prior * pl["trans"][pattern]
    den = num + (1 - prior) * pl["cis"][pattern]
    if den == 0:
        raise ValueError("pattern impossible under both hypotheses")
    return num / den


def parental_test_summary() -> dict:
    rows = []
    for pi in (0.0, 0.05, 0.10):
        for eps in (0.0, 0.01):
            row = {"pi_dn": pi, "eps": eps, "p_decisive": 1 - pi}
            for pat in PATTERNS:
                try:
                    row["posterior_trans_" + pat] = posterior_after_parental_test(pat, 0.5, pi, 0.5, eps)
                except ValueError:
                    row["posterior_trans_" + pat] = None
            rows.append(row)
    ext = {
        "DN_M_if_q_is_0": posterior_after_parental_test("DN_M", 0.5, 0.05, 0.0),
        "DN_M_if_q_is_1": posterior_after_parental_test("DN_M", 0.5, 0.05, 1.0),
    }
    return {
        "expected_patterns": PATTERNS,
        "expected_result_if_trans": "pattern T (each parent carries a different one of the two variants)",
        "expected_result_if_cis": "pattern C (one parent carries both variants, the other none)",
        "grid": rows,
        "de_novo_pattern_extremes_pi_dn_0.05": ext,
        "assumptions": "eps (recombination between variants, mosaicism, test error) is an assumption; pattern 'both de novo' is neglected; parentage is assumed biological.",
    }


# ---------------------------------------------------------------- long-read phasing
def p_span_given_cover(read_len: float, dist: float) -> float:
    """Probability that a read of length L covering variant 1 also covers variant 2 at distance d: max(0, L-d)/L.

    A read of fixed length L placed uniformly covers a given base from L start positions; it covers both bases when its
    start lies in a window of L-d positions, so the conditional probability is (L-d)/L for d < L, else 0.
    """
    if read_len <= 0 or dist < 0:
        raise ValueError("invalid lengths")
    return max(0.0, read_len - dist) / read_len


def poisson_cdf_below(m: int, lam: float) -> float:
    return sum(math.exp(-lam) * lam**k / math.factorial(k) for k in range(m))


def p_phase_decisive(depth: float, read_len: float, dist: float, min_reads: int = 1) -> float:
    """P(at least `min_reads` reads span both variants) = 1 - sum_{k<m} e^-lam lam^k / k!, lam = D*(L-d)/L.

    D is the mean depth of aligned reads over the site (both haplotypes pooled). Every spanning read is informative:
    in trans each read shows one alternate and one reference allele; in cis a read shows two alternates or two
    references. Read length is treated as fixed (a length distribution replaces (L-d)/L by its base-weighted mean).
    """
    lam = depth * p_span_given_cover(read_len, dist)
    return 1 - poisson_cdf_below(min_reads, lam)


def min_depth(read_len: float, dist: float, min_reads: int = 1, target: float = 0.95) -> float | None:
    """Smallest depth with P(decisive) >= target; None if no single read can span (d >= L)."""
    p = p_span_given_cover(read_len, dist)
    if p == 0:
        return None
    lo, hi = 0.0, 1.0
    while p_phase_decisive(hi, read_len, dist, min_reads) < target:
        hi *= 2
    for _ in range(80):
        mid = (lo + hi) / 2
        if p_phase_decisive(mid, read_len, dist, min_reads) >= target:
            hi = mid
        else:
            lo = mid
    return hi


def long_read_table(distances: tuple, min_reads_options=(1, 3), target: float = 0.95) -> list:
    rows = []
    for dist in distances:
        for L in READ_LENGTHS_BP:
            row = {"distance_bp": dist, "read_length_bp": L, "p_span_given_cover": p_span_given_cover(L, dist)}
            for m in min_reads_options:
                d = min_depth(L, dist, m, target)
                row[f"min_depth_for_{int(target*100)}pct_with_{m}_reads"] = None if d is None else round(d, 2)
            for depth in (10, 20, 30):
                row[f"p_decisive_depth_{depth}x_min_reads_1"] = p_phase_decisive(depth, L, dist, 1)
            rows.append(row)
    return rows


# ---------------------------------------------------------------- assembly
def build(diff: dict) -> dict:
    sole = bub1b_sole_causal(diff)
    f0, a0 = HEADLINE["penetrance"], HEADLINE["alt_cause"]
    grid = sensitivity_table(1.0, 0.5)
    vals = [r["posterior_trans"] for r in grid]
    distances = (1_000, 5_000, 10_000, BUB1B_GENE_LENGTH_BP)
    return {
        "assumes": "BUB1B is the causal gene; the gene-level fact used is only that BUB1B is the sole 'causal candidate' in out/differential_diagnosis.json",
        "bub1b_sole_causal_candidate_in_differential": sole,
        "hypotheses": {"H_trans": "the two BUB1B alleles are on different chromosomes", "H_cis": "the two BUB1B alleles are on the same chromosome"},
        "prior": {
            "p_trans": prior_trans(),
            "p_trans_if_de_novo_probability_0.1": prior_trans(0.1),
            "derivation": "both inherited: different parents with probability 1/2 each under independent transmission; one de novo: trans iff it arose on the other parent's chromosome, 1/2 for any paternal fraction q; hence 1/2",
            "prior_grid_posteriors_at_headline": [{"prior": p, "posterior_trans": posterior_trans(f0, a0, 1.0, p)} for p in PRIOR_GRID],
        },
        "likelihood": {
            "model": "P(affected|trans)=c*f+(1-c); P(affected|cis)=c*a+(1-c)",
            "alt_cause_source": MVA_INCIDENCE_SOURCE,
            "penetrance_note": "f is an assumption varied over the grid, not a sourced value",
        },
        "headline": {
            **HEADLINE,
            "posterior_trans": posterior_trans(f0, a0, 1.0, 0.5),
            "range_over_penetrance_and_alt_cause_grid": [min(vals), max(vals)],
        },
        "sensitivity_penetrance_by_alt_cause": grid,
        "sensitivity_p_causal_at_headline": causal_table(f0, a0),
        "sensitivity_p_causal_at_worst_grid_cell": causal_table(0.5, 1e-3),
        "p_causal_needed_for_posterior_trans": {
            "0.95": causal_threshold(0.95, f0, a0),
            "0.90": causal_threshold(0.90, f0, a0),
            "0.75": causal_threshold(0.75, f0, a0),
        },
        "parental_sanger_test": parental_test_summary(),
        "long_read_phasing": {
            "formula": "p_span|cover = max(0, L-d)/L; lambda = D*p_span|cover; P(decisive) = 1 - sum_{k<m} exp(-lambda) lambda^k / k!",
            "gene_body_length_bp": BUB1B_GENE_LENGTH_BP,
            "gene_body_length_source": ENSEMBL_SOURCE,
            "distances_are_scenarios": "the real distance between the two alleles is not used here; the largest scenario is the whole gene body",
            "read_length_range_bp": [READ_LENGTHS_BP[0], READ_LENGTHS_BP[-1]],
            "table": long_read_table(distances),
            "beyond_read_length": "if d >= L no single read spans both alleles; phase then needs parental testing or chaining through intermediate heterozygous markers",
        },
    }


def main() -> None:
    diff = json.loads(DIFF.read_text(encoding="utf8"))
    doc = build(diff)
    OUT.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf8")
    h = doc["headline"]
    print("headline posterior P(trans):", round(h["posterior_trans"], 7), "range:", [round(x, 5) for x in h["range_over_penetrance_and_alt_cause_grid"]])


if __name__ == "__main__":
    main()
