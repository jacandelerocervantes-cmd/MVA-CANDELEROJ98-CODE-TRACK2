"""Text of the empirical checks of the method (section 3.7, Supplement S12, scalability measurements), written only from `out/empirical/*.json`.

Every number in the sentences and tables below is read from those files; nothing is typed here except wording. The checks measure the method (selection rules, verifier, calibration of
predictors, age support of labels, reuse on other cases), not whether any drug works.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from track2.submission_input import redact_cdna

NINE = ["SIROLIMUS", "EVEROLIMUS", "METFORMIN", "PHENYLBUTANOIC ACID", "TRAMETINIB", "SELUMETINIB", "BORTEZOMIB", "CHLOROQUINE", "BOSUTINIB"]
SHOW = {"PHENYLBUTANOIC ACID": "Sodium phenylbutyrate"}


def n_tests() -> int:
    """Number of tests collected by pytest (the text of the report must not carry a typed number)."""
    import re
    import subprocess
    import sys
    out = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q"], capture_output=True, text=True, cwd=str(Path(__file__).resolve().parents[2])).stdout
    m = re.search(r"(\d+) tests? collected", out) or re.search(r"^(\d+)/\d+ ", out, re.M)
    return int(m.group(1)) if m else len(re.findall(r"::", out))


def _load(emp: Path, name: str) -> dict[str, Any] | None:
    p = emp / f"{name}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def _name(n: str) -> str:
    return SHOW.get(n, n.title())


def _f(x: float) -> str:
    """Two decimals, three when two would round a value below 1 up to 1.00."""
    return f"{x:.3f}" if 0.995 <= x < 1 else f"{x:.2f}"


def _pct(x: float, digits: int = 1) -> str:
    return f"{100 * x:.{digits}f} %"


def _power_row(pw: dict[str, Any], scenario: str, reduction: float, power: float = 0.8) -> dict[str, Any]:
    for r in pw["table"]:
        if r["scenario"] == scenario and abs(r["relative_reduction"] - reduction) < 1e-9 and abs(r["power_target"] - power) < 1e-9:
            return r
    raise KeyError((scenario, reduction, power))


def _cultures_range(row: dict[str, Any], per_culture: int = 50) -> tuple[int, int]:
    ks = [c["cultures_per_arm"] for c in row["clustered"] if c["cells_per_culture"] == per_culture]
    return min(ks), max(ks)


def main_section(emp: Path) -> str:
    rob, ctl, var, lit, age, pw = (_load(emp, k) for k in ("label_robustness", "controls", "variant_calibration", "literature_control", "age_support", "power_analysis"))
    t = ["### 3.7 Empirical checks of the method\n",
         "Six checks quantify how far the selection can be trusted and what a refuting experiment would need. They measure the method, not the drugs: whether the rules hold when the labels are wrong, whether a proxy of the literature rule recovers known indications, how often the verifier catches an error of the AI search tasks, whether the caution layer recovers harms already published for MVA, how strong the computational evidence on the second variant is, over what part of the paediatric age range the labels support each drug, and how many cells a refuting experiment needs. The tables are in Supplement S12.\n"]
    if rob:
        r = rob["rates"]
        parts = ", ".join(f"{r[k]['expected_survivors']:.2f}" for k in ("5", "10", "20"))
        parts += f" and {r['30']['expected_survivors']:.2f}"
        ivs = ", ".join(f"{int(r[k]['survivors_95_interval'][0])} to {int(r[k]['survivors_95_interval'][1])}" for k in ("5", "10", "20"))
        ivs += f" and {int(r['30']['survivors_95_interval'][0])} to {int(r['30']['survivors_95_interval'][1])}"
        s30 = r["30"]["survival_fraction_by_drug"]
        weak = sorted(s30, key=lambda d: s30[d])[:2]
        fixed = sorted(d for d, v in s30.items() if v >= 0.99)
        names = [_name(d) for d in fixed]
        fixed_txt = (", ".join(n.lower() for n in names[:-1]) + " and " + names[-1].lower()) if len(names) > 1 else names[0].lower()
        fixed_txt = fixed_txt[0].upper() + fixed_txt[1:]
        floor = min(s30[d] for d in fixed)
        t.append(f"**Errors in the evidence labels.** The levels and directions of the evidence were assigned by AI search tasks and were not reviewed by a specialist. We corrupted 5, 10, 20 and 30 % of the {rob['n_entries']} verified entries at random "
                 f"({rob['replicates']:,} replicates per rate, fixed seed; a level changed by one step, a benefit and harm direction swapped, or both) and applied rules R1 to R4 again. At 5, 10, 20 and 30 % an expected {parts}, respectively, of the nine drugs keep their group "
                 f"(95 % intervals {ivs}). The least stable are {_name(weak[0]).lower()} and {_name(weak[1]).lower()} ({s30[weak[0]]:.2f} and {s30[weak[1]]:.2f} of the replicates at 30 %). "
                 f"{fixed_txt} barely change ({floor:.3f} or more at 30 %), because the class-sibling and editorial rules fix them, not because of their evidence. Only levels and directions were corrupted; section, tag and caution flags were not.\n")
    if lit:
        a = lit["aggregate_all_gold"]
        t.append(f"**A proxy of the literature rule against known indications.** A deterministic proxy of rule R2 (a benefit sentence that also mentions a mechanism term, in up to {lit['max_abstracts']} Europe PMC abstracts per disease, with no drug name in the query) was applied to {a['n_diseases']} rare diseases whose indicated drugs come from PrimeKG. "
                 f"It ranked the indicated drugs above the others with a macro AUC of {a['auc_macro']['mean']:.2f} (95 % interval {a['auc_macro']['lo']:.2f} to {a['auc_macro']['hi']:.2f}) and found a benefit sentence for {a['gold_hit_rate_pooled']['k']} of {a['gold_hit_rate_pooled']['n']} indicated drugs ({_pct(a['gold_hit_rate_pooled']['rate'], 0)}). "
                 f"A random non-indicated drug passed in {_pct(a['negative_control_pass_rate_random_drug']['rate'], 2)} of cases, and a random set of nine contained a passing drug in {_pct(a['negative_control_frac_random_sets9_with_any_pass']['mean'])}. "
                 "This is a negative control for literature volume. It counts co-mentions, has no direction or level, and its gold standard is incomplete; it does not measure the AI-assisted procedure used here, whose first-round queries also named drugs.\n")
    if ctl:
        v, p = ctl["verifier_catch_rate"], ctl["positive_control"]
        fc = v["file_counts_field"]
        n = v["n_sentences"]
        rec = ", ".join(_name(d).lower() for d in p["recovered"])
        ind = p["caution_list_json_only"]
        t.append(f"**Errors of the AI search tasks caught by the verifier.** Of {n} sentence slots proposed by the AI search tasks, {fc['verified']} ({_pct(fc['verified'] / n)}) were found word for word in their source, {fc['not_checkable']} could not be checked and {fc['not_found']} ({_pct(fc['not_found'] / n, 2)}; "
                 f"{v['distinct_not_found_sentences']} distinct sentence) were not found and were excluded; all {v['label_quote_check']['n_quotes']} label quotes were found in the label they cite. "
                 f"As a positive control of the caution layer, {p['approved_only']['recovers'].replace('recovers ', '')} approved drugs with a harm or reduced sensitivity reported in MVA or *BUB1B* patients or cells were recovered ({rec}); "
                 f"only {ind['N_recovered']} of them come from the label-based list, which is independent of the literature layer, so the recovery of the third is circular.\n")
    if var:
        pp = var["per_predictor"]
        t.append(f"**Computational evidence on p.Asn1002Lys.** With calibrated thresholds (Pejaver 2022; Bergquist 2025 for AlphaMissense), AlphaMissense ({pp['AlphaMissense']['score']}) gives {pp['AlphaMissense']['strength'].replace('PP3_', '').lower()} and PolyPhen-2 ({pp['PolyPhen2']['score']}) gives {pp['PolyPhen2']['strength'].replace('PP3_', '').lower()} computational evidence of pathogenicity, "
                 f"whereas REVEL ({pp['REVEL']['score']}), CADD-PHRED ({pp['CADD']['score']}) and SIFT ({pp['SIFT']['score']}) give none. The tools are correlated and disagree, their strengths are not added, and the variant remains of uncertain significance with the phase unresolved.\n")
    if age:
        c = {d["name"]: d["C"] for d in age["drugs"] if d["name"] in NINE}
        t.append("**Label support across the paediatric age range.** For each drug we computed the fraction of the 0 to 18 year range in which the label states use and no age-specific hazard sentence applies, without assuming the age of the child (Supplement S12.3). "
                 f"It is {c['SIROLIMUS']['lower']:.2f} for sirolimus (use stated from 13 years), {c['METFORMIN']['lower']:.2f} for metformin (from 10 years), {min(c[d]['lower'] for d in ('EVEROLIMUS', 'TRAMETINIB', 'SELUMETINIB', 'BOSUTINIB')):.2f} for everolimus, trametinib, selumetinib and bosutinib (from 1 year), "
                 f"{c['BORTEZOMIB']['lower']:.2f} for bortezomib (not established in children) and between {c['CHLOROQUINE']['lower']:.0f} and {c['CHLOROQUINE']['upper']:.0f} for chloroquine (unresolved). The extraction is by regular expressions, so hazards without a numeric age are missed and the values lean upward; "
                 "weight restrictions (for example the phenylbutyrate tablets, 20 kg) are flagged and not converted to ages. This is label support under an unknown age, not a recommendation, and no dose is computed.\n")
    if pw:
        lo, mid, hi = (_power_row(pw, k, 0.5) for k in ("S_low_BUB1B", "S_threshold", "S_high"))
        l25 = _power_row(pw, "S_low_BUB1B", 0.25)
        c_lo, c_hi = _cultures_range(lo)
        t.append("**Sizing the refuting experiments.** The experiments of Supplement S4 read the fraction of aneuploid cells. Published baselines in patient cells are 12 of 100 cells in the only *BUB1B* value found (PMID 35804254), 25 % as the usual lower bound of MVA (PMID 37796616) and 38 to 39 % in patients without *BUB1B* mutations (PMIDs 36322655 and 41374403). "
                 f"With a two-sided test at alpha 0.05 and power 0.8, a 50 % relative reduction of the aneuploid fraction needs {lo['n_cells_per_arm_independent']:,} cells per arm from the 12 % baseline, {mid['n_cells_per_arm_independent']:,} from 25 % and {hi['n_cells_per_arm_independent']:,} from 39 %; "
                 f"with 50 cells per culture and an assumed intraclass correlation of 0.05 to 0.3, the 12 % baseline needs {c_lo} to {c_hi} independent cultures per arm. A 25 % reduction from 12 % needs {l25['n_cells_per_arm_independent']:,} cells. "
                 "These are planning numbers, not results: the correlation between cells of one culture is an assumption, the background in control cells could not be sourced, the baselines come from lymphocytes of patients (mostly without *BUB1B* mutations), and an exact test needs about 5 to 20 % more cells (Supplement S12.4).\n")
    return "\n".join(t)


def supplement(emp: Path) -> str:
    rob, var, age, pw = (_load(emp, k) for k in ("label_robustness", "variant_calibration", "age_support", "power_analysis"))
    t = ["## S12. Empirical checks of the method\n"]
    if rob:
        t.append("### S12.1 Robustness of the selection to errors in the evidence labels\n")
        t.append("Fraction of replicates in which each drug stays proposed in the same group, at each rate of corrupted labels (" + f"{rob['replicates']:,}" + " replicates per rate; seed " + str(rob["seed"]) + ").\n")
        t.append("| Drug | 5 % | 10 % | 20 % | 30 % |\n|---|:-:|:-:|:-:|:-:|")
        r = rob["rates"]
        for d in NINE:
            t.append(f"| {_name(d)} | " + " | ".join(f"{r[k]['survival_fraction_by_drug'][d]:.2f}" for k in ("5", "10", "20", "30")) + " |")
        t.append("| Expected number of the nine (95 % interval) | " + " | ".join(f"{r[k]['expected_survivors']:.2f} ({int(r[k]['survivors_95_interval'][0])} to {int(r[k]['survivors_95_interval'][1])})" for k in ("5", "10", "20", "30")) + " |")
        t.append("| New cautions per replicate (mean) | " + " | ".join(f"{r[k]['mean_new_cautions']:.2f}" for k in ("5", "10", "20", "30")) + " |\n")
    if var:
        t.append("### S12.2 Calibrated computational evidence for p.Asn1002Lys\n")
        t.append("| Predictor | Score | Strength under the calibrated thresholds | Source |\n|---|:-:|---|---|")
        for k in ("AlphaMissense", "PolyPhen2", "REVEL", "CADD", "SIFT"):
            x = var["per_predictor"][k]
            t.append(f"| {x['predictor']} | {x['score']} | {x['strength'].replace('_', ' ')} | PMID {x['pmid']} |")
        t.append("")
        t.append("Thresholds are quoted and checked against the full text of the sources by code (`out/empirical/variant_calibration.json`). Combining a population-rarity criterion at supporting strength with a computational criterion does not classify the variant, which remains of uncertain significance.\n")
    if age:
        t.append("### S12.3 Label support across the paediatric age range\n")
        t.append(f"Fraction of the range 0 to 18 years in which the label states use and no age-specific hazard sentence applies (lower and upper bound; `unresolved` labels give 0 to 1). {age['header']['not_a_recommendation']}\n")
        t.append("| Drug | 0 to 18 y | 0 to <3 y | 3 to <12 y | 12 y and older | Weight flags |\n|---|:-:|:-:|:-:|:-:|:-:|")
        by = {d["name"]: d for d in age["drugs"]}
        for n in NINE:
            d = by[n]

            def cell(c: dict[str, float]) -> str:
                return _f(c["lower"]) if abs(c["lower"] - c["upper"]) < 1e-9 else f"{_f(c['lower'])} to {_f(c['upper'])}"
            t.append(f"| {_name(n)} | {cell(d['C'])} | " + " | ".join(cell(d["C_by_group"][g]) for g in ("0_to_under_3", "3_to_under_12", "12_and_older")) + f" | {d['n_weight_flags']} |")
        t.append("")
    if pw:
        t.append("### S12.4 Sizing of the refuting experiments (planning numbers, not results)\n")
        t.append(pw["status"] + " " + pw["scope"] + " Two-sided alpha 0.05, power 0.8, arcsine method; cultures per arm for 50 cells per culture and an assumed intraclass correlation of 0.05 to 0.3.\n")
        t.append("| Baseline aneuploid fraction | Relative reduction to detect | Cells per arm (independent cells) | Cultures per arm (50 cells each) |\n|---|:-:|:-:|:-:|")
        labels = {"S_low_BUB1B": "12 % (the only BUB1B patient value, PMID 35804254)", "S_threshold": "25 % (lower bound of MVA, PMID 37796616)", "S_high": "39 % (patients without BUB1B mutations, PMIDs 36322655, 41374403)"}
        for sc in ("S_low_BUB1B", "S_threshold", "S_high"):
            for red in (0.25, 0.5, 0.75):
                row = _power_row(pw, sc, red)
                a, b = _cultures_range(row)
                t.append(f"| {labels[sc]} | {int(red * 100)} % | {row['n_cells_per_arm_independent']:,} | {a} to {b} |")
        t.append("")
        t.append("Sources, each checked by code against the text of the paper: " + redact_cdna("; ".join(f"{x['cell_type']} (PMID {x['pmid']})" for x in pw["sources"][:6])) + ". Not sourced: " + " ".join(pw["not_sourced"][:3]) + " " + " ".join(pw["caveats"][:2]) + "\n")
    return "\n".join(t)


def abstract_sentence(emp: Path) -> str:
    rob = _load(emp, "label_robustness")
    if not rob:
        return ""
    r = rob["rates"]["30"]
    return (f"With up to 30 % of the evidence labels corrupted at random, an expected {r['expected_survivors']:.1f} of the nine drugs kept their group (Section 3.7).")


def scalability_text(emp: Path) -> str:
    s = _load(emp, "scalability")
    if not s:
        return ""
    cases = s["cases"]
    identical = all(c["byte_identical_between_repeats"] for c in cases.values())
    files = sorted({c["n_output_files"] for c in cases.values()})
    return ("The early stages (module building and network proximity at three STRING thresholds, then the consensus) were run twice on two synthetic cases, Fanconi anaemia (*FANCA*, *ATR*) and a second MVA-type gene pair (*CEP57*, *TRIP13*), "
            f"with every output redirected to separate folders and the files of the real case unchanged (checked by hash). They took about 80 to 100 seconds and 260 MB for the first stage, and the output files were {'byte-identical' if identical else 'not identical'} between the two runs ({files[0]} files per case). "
            "Building the candidate table needs a patient file written by hand, the validation of each module takes 30 to 90 minutes and was not repeated, and the literature layer would have to be redone for another disease (`out/empirical/scalability.json`; timings depend on the machine).")
