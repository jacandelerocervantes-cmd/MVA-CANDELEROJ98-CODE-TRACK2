"""Text of sections 3.8 to 3.11 and Supplement S13 of the report (gene landscape, differential diagnosis, mechanism of the variants, phase, case profile),
written only from the JSON files of `out/` produced by `disease_landscape`, `differential`, `variant_level_consequence`, `variant_mechanism`, `phase_inference`,
`case_profile` and `case_contrast`. Nothing but wording is typed here."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

REVIEW_N = 31
MAP = {"causal candidate": "mechanism_candidate", "contextual": "pathway_context_only", "insufficient evidence": "insufficient_evidence"}
GENOTYPE = {"compound_het": "compound heterozygous", "homozygous": "homozygous", "hemizygous": "hemizygous", None: "unpaired structural call"}


def _load(out: Path, name: str) -> Any | None:
    p = out / f"{name}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def _classes(cc: dict[str, int]) -> str:
    return ", ".join(f"{n} {k.replace('structural call not', 'structural calls not') if n > 1 else k}" for k, n in cc.items())


def _clean(text: str) -> str:
    """Statements keep their measured values; units are written with symbols and the final full stop is left to the sentence that embeds them."""
    text = re.sub(r"(\d) A2\b", r"\1 Å²", text)
    text = re.sub(r"(\d) A\b", r"\1 Å", text)
    return text.rstrip(". ")


def _agreement(diff: dict[str, Any], screen_path: Path) -> tuple[int, int, list[str]]:
    if not screen_path.exists():
        return 0, 0, []
    screen = json.loads(screen_path.read_text(encoding="utf-8")).get("genes", [])
    status = {g["gene"]: g.get("status") for g in screen if isinstance(g, dict) and "gene" in g}
    agree, total, diff_genes = 0, 0, []
    for gene, row in diff["genes"].items():
        if gene not in status:
            continue
        total += 1
        verdict = row["verdict"]
        mapped = next((v for k, v in MAP.items() if verdict.startswith(k.split()[0])), None)
        if mapped == status[gene]:
            agree += 1
        else:
            diff_genes.append(gene)
    return agree, total, diff_genes


def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / (1 + z * z / n)
    return max(0.0, centre - half), min(1.0, centre + half)


def _interval(x: dict[str, Any]) -> str:
    """95 % Wilson interval from the denominator stated in the source (k/n, or a percentage of 'the group of N')."""
    v = x.get("value") or ""
    m = re.fullmatch(r"(\d+)/(\d+)", v)
    if m:
        k, n = int(m.group(1)), int(m.group(2))
    else:
        mp = re.fullmatch(r"([\d.]+)%", v)
        if not (mp and re.search(r"group of (\d+)", x.get("population", ""))):
            return ""
        n = REVIEW_N  # the percentages of the 2022 review are exact fractions of 31, although its text says 'group of 30'
        k = round(float(mp.group(1)) / 100 * n)
    lo, hi = _wilson(k, n)
    return f"; 95 % interval {100 * lo:.0f} to {100 * hi:.0f} %, n = {n}"


def main_sections(out: Path) -> str:
    land, diff, cons = _load(out, "mva_landscape"), _load(out, "differential_diagnosis"), _load(out, "variant_consequence_classes")
    mech, phase, contrast = _load(out, "variant_mechanism"), _load(out, "phase_inference"), _load(out, "case_contrast")
    t: list[str] = []
    if land and diff:
        by_tier, by_class = land["n_genes_by_tier"], land["n_genes_by_mechanism_class"]
        in_t1 = sorted(g["symbol"] for g in land["genes"] if g.get("in_track1_output"))
        agree, total, other = _agreement(diff, Path("config/mechanism_screen.json"))
        t.append("### 3.8 Genes of mosaic variegated aneuploidy and differential diagnosis of the Track 1 genes\n")
        t.append(f"The genes associated with MVA were derived by code from the Human Phenotype Ontology annotations and Orphanet (release and access date in `out/mva_landscape.json`), with the mechanism class of each gene assigned by a fixed rule from UniProt keywords and Gene Ontology terms and one verbatim Europe PMC sentence per gene verified by code. "
                 f"The result is {by_tier['curated']} curated genes and {by_tier.get('proposed, not curated', 0)} proposed in recent reports but not yet curated; {by_class.get('spindle-assembly checkpoint', 0)} act on the spindle-assembly checkpoint, {by_class.get('centrosome/centriole', 0)} on the centrosome and {by_class.get('DNA repair/replication', 0)} on DNA repair or replication (Supplement S13). "
                 f"Only {' and '.join('*' + g + '*' for g in in_t1)} of these genes appear in the Track 1 output. The same procedure, without changes to its rules, was run for Fanconi anaemia (Supplement S13).\n")
        t.append("The genes of the Track 1 output were then classified by rules written before they were applied (`protocol/addendum_differential.md`): whether the gene is a curated MVA gene, whether the genotype fits the annotated inheritance, the consequence class, the phenotype fit and ClinVar.\n")
        t.append("| Gene | In the MVA landscape | Genotype | Variant classes | Phenotype score (rank) | ClinVar | Verdict |\n|---|---|---|---|---|---|---|")
        for gene, r in diff["genes"].items():
            cls = _classes(cons["genes"][gene]["class_counts"]) if cons and gene in cons["genes"] else "n/a"
            c4 = r["c4"]
            ph = f"{c4['phenotype_score']:.3f} ({c4['phenotype_rank']})" if c4.get("phenotype_rank") is not None and c4.get("phenotype_score") is not None else ("no score" if c4.get("has_hpo") else "no HPO annotation")
            t.append(f"| *{gene}* | {r['c1']} | {GENOTYPE.get(r['genotype_type'], r['genotype_type'])} | {cls} | {ph} | {r['c5']} | {r['verdict']} |")
        t.append("")
        t.append(f"*BUB1B* is the only causal candidate under these rules. The rule-based verdicts agree with the hand-made mechanism screen for {agree} of {total} genes" + (f"; the difference is {', '.join('*' + g + '*' for g in other)}, which the screen keeps as pathway context and the rules, which depend on the landscape, rate as insufficient evidence" if other else "") + ". The phase of the two variants is not part of these rules and is treated in Section 3.10.\n")
    if mech:
        a, b = mech["variant_Leu737Ter"], mech["variant_Asn1002Lys"]
        t.append("### 3.9 Mechanism of the two *BUB1B* variants\n")
        t.append(f"Both statements below are written by rule from features measured in public sources (`out/variant_mechanism.json`), and each clause names the feature that triggered it. **p.Leu737\\*:** {_clean(a['mechanism']['statement'])}. "
                 f"**p.Asn1002Lys:** {_clean(b['mechanism']['statement'])}.\n")
        cal = _load(out / "empirical", "variant_calibration")
        if cal:
            pp = cal["per_predictor"]
            t.append(f"With calibrated thresholds the computational evidence for p.Asn1002Lys is {pp['AlphaMissense']['strength'].replace('PP3_', '').lower()} for AlphaMissense, {pp['PolyPhen2']['strength'].replace('PP3_', '').lower()} for PolyPhen-2 and none for REVEL, CADD and SIFT (Section 3.7), so the variant remains of uncertain significance by formal criteria. "
                     "Its predicted effect is on the abundance of BubR1, consistent with the published analysis of MVA1 missense alleles (Suijkerbuijk 2010).\n")
        tests = a["mechanism"].get("what_would_test_it", []) + b["mechanism"].get("what_would_test_it", [])
        if tests:
            t.append("What would test these statements: " + "; ".join(_clean(x) for x in tests) + "." + "\n")
    if phase:
        h = phase["headline"]
        lo, hi = h["range_over_penetrance_and_alt_cause_grid"]
        rows = phase.get("sensitivity_p_causal_at_headline", [])
        pick = {round(r["p_causal"], 2): r["posterior_trans"] for r in rows}
        t.append("### 3.10 Phase of the two variants\n")
        t.append("Phase cannot be resolved from short reads. We therefore state what the data allow as an explicit inference and not as a finding. Under the assumption that *BUB1B* is the causal gene, the prior probability that the two alleles are on different chromosomes (in trans) is "
                 f"{h['prior']:.1f}, because two inherited rare alleles come from different parents with that probability; the observation of an affected child with a recessive disorder, with penetrance {h['penetrance']} and a probability of {h['alt_cause']:.6f} that a carrier of both alleles in cis is affected by another cause (incidence of PCS/MVA of 1 in 1,000,000, from a verified sentence of PMID {phase['likelihood']['alt_cause_source']['pmid']}), gives a posterior probability of trans of {h['posterior_trans']:.6f} (range {lo:.3f} to {hi:.6f} over the grid of penetrance and alternative cause). "
                 + ("If the probability that *BUB1B* is the causal gene is lowered, the posterior falls to " + ", ".join(f"{pick[k]:.3f} at {k}" for k in sorted(pick, reverse=True) if k in (0.99, 0.9, 0.5, 0.0)) + ". " if pick else "")
                 + "The inference is conditional on the penetrance and the other-cause probability, which are assumptions. Phase would be decided by testing the two variants in the parents (one in each parent indicates trans; both in one parent indicates cis) or by long reads: a single read can phase two alleles only if it is longer than their distance, and the gene body is "
                 f"{phase['long_read_phasing']['gene_body_length_bp']:,} bp, so alleles further apart than the read length need parental testing.\n")
    if contrast:
        t.append("### 3.11 Minimal case profile and comparison with published cases\n")
        t.append("Each individual is described by the same minimal record (`src/track2/case_profile.py`): age band, sex, phenotypic features as HPO terms with present, absent or unknown, tumour history and status, and genomic findings in categorical terms only, with an explicit \"unknown\" state for every field and a validator that rejects coordinates, alleles, genotype calls and read metrics. "
                 "The profile of the child was compared with the frequencies reported in published *BUB1B* (MVA1) cases, taken from verified sentences (`out/case_contrast.json`). The comparison lists what to check; it is not a diagnosis.\n")
        t.append("| Feature | Published frequency (source) | Child | Reading |\n|---|---|---|---|")
        for r in contrast["rows"]:
            pf = r["published_frequency"]
            if isinstance(pf, list) and pf:
                freq = "; ".join((f"{x['value']}" if x.get("value") else "named, no frequency") + f" (PMID {x['pmid']}" + ("; 2008 series" if x["pmid"] == "18548531" else "") + _interval(x) + ")" for x in pf)
            else:
                freq = "not reported in the sources"
            t.append(f"| {r['feature'][0].upper() + r['feature'][1:]} | {freq} | {r['child']} | {r['reading']} |")
        t.append("")
        t.append(f"Intervals are 95 % Wilson intervals from the denominators of the sources. The eight percentages of the 2022 review (PMID 35804254) are exact fractions of {REVIEW_N} cases, although its text says a group of 30, so {REVIEW_N} is used.\n")
        g = contrast["genotype"]
        t.append(f"Genotype: {_clean(g['reading'])}. " + f"The frequency of two truncating variants and the counts of truncating plus missense against two truncating are not reported in the sources. {contrast['limits']}\n")
    return "\n".join(t)


def supplement(out: Path) -> str:
    t = ["## S13. Gene landscape of the diseases analysed\n"]
    for name, title in (("mva_landscape", "S13.1 Mosaic variegated aneuploidy"), ("landscape_fanconi_anaemia", "S13.2 Fanconi anaemia (same procedure, no change of rules)")):
        d = _load(out, name)
        if not d:
            continue
        t.append(f"### {title}\n")
        t.append(f"Derived by code from public sources (access {d['header']['access_date']}); association evidence is that of the cited sources. Genes by mechanism class: " + ", ".join(f"{v} {k}" for k, v in d["n_genes_by_mechanism_class"].items()) + f"; {d['n_genes_with_verified_sentence']} of {d['n_genes']} genes have a sentence verified against the source.\n")
        t.append("| Gene | Subtype | Mechanism class | Inheritance (as annotated) | Tier | Sources | Sentence |\n|---|---|---|---|---|:-:|---|")
        for g in d["genes"]:
            ev = g["association_evidence"]
            t.append(f"| *{g['symbol']}* | {g['disease_subtype']} | {g['mechanism_class']} | {'; '.join(g['inheritance_as_annotated'])} | {g['tier']} | {g['n_sources']} | {ev.get('evidence', '')} (PMID {ev.get('pmid', 'n/a')}) |")
        t.append("")
    return "\n".join(t)


def abstract_sentence(out: Path) -> str:
    diff = _load(out, "differential_diagnosis")
    land = _load(out, "mva_landscape")
    if not (diff and land):
        return ""
    causal = [g for g, r in diff["genes"].items() if r["verdict"] == "causal candidate"]
    return (f"A rule-based differential diagnosis of the Track 1 genes against {land['n_genes_by_tier']['curated']} curated MVA genes leaves " + (", ".join(f"*{g}*" for g in causal) if causal else "no gene") + " as the only causal candidate (Section 3.8).")


def scalability_sentence(out: Path) -> str:
    d = _load(out, "landscape_fanconi_anaemia")
    if not d:
        return ""
    return (f"The same landscape procedure, with unchanged rules, gave for Fanconi anaemia {d['n_genes']} genes in {len(d['n_genes_by_mechanism_class'])} mechanism classes, {d['n_genes_with_verified_sentence']} of them with a verified sentence (Supplement S13).")
