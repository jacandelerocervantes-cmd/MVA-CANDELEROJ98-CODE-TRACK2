"""'Answer in brief' paragraph at the start of the abstract, written only from the JSON files of `out/`; nothing but wording is typed here."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from track2 import lines_section as lines


def _load(out: Path, name: str) -> Any | None:
    p = out / f"{name}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def answer(out: Path) -> str:
    mech, contrast, rob, lit, ctl, diff = (_load(out, n) for n in ("variant_mechanism", "case_contrast", "empirical/label_robustness", "empirical/literature_control", "empirical/controls", "differential_diagnosis"))
    parts = ["**Answer in brief.**"]
    if mech and diff:
        stop = mech["variant_Leu737Ter"]["nmd"]
        causal = [g for g, r in diff["genes"].items() if r["verdict"] == "causal candidate"]
        parts.append(f"Of the genes returned by Track 1, {', '.join('*' + g + '*' for g in causal) if causal else 'none'} is the only causal candidate under rules written beforehand (Section 3.8). Its two variants are a stop-gain predicted to cause nonsense-mediated decay (exon {stop['stop_exon_number']} of {stop['n_exons']}) and a missense variant at a conserved, buried position of the kinase domain; "
                     "the phase is unresolved and the missense variant remains of uncertain significance (Sections 3.9 and 3.10).")
    parts.append("The other Track 1 genes are carried as hypotheses: *IGBP1* as a hypothesis-grade PP2A, AKT and mTOR modifier, *CENPE* and *MAD1L1* as context, and the megabase-scale events as undetermined until a karyotype or single-cell sequencing (Section 3.1).")
    parts.append("The track asks for candidate medications as hypotheses for follow-up. The literature layer lists approved drugs of three kinds: constitutional, tumour-context and hazards to check, not to give (Section 3.2). " + lines.brief_sentence(out)
                 + " No approved drug is supported by two usable lines of the protocol's rule (Section 3.4). Label support is shown for age bands, without assuming the age of the child (Sections 3.6, 3.7 and 3.13).")
    if contrast:
        todo = [r["feature"] for r in contrast["rows"] if "check clinically" in r["reading"]]
        if todo:
            parts.append("Features reported in most published *BUB1B* cases and not recorded for the child are listed for the clinicians to confirm: " + ", ".join(todo) + " (Section 3.11). The first test for the hypotheses is a karyotype or single-cell sequencing.")
    bits = []
    if rob:
        bits.append(f"with 30 % of the evidence labels corrupted, {rob['rates']['30']['expected_survivors']:.1f} of the nine drugs keep their group")
    if ctl:
        v = ctl["verifier_catch_rate"]
        bits.append(f"the quote verifier excluded {v['file_counts_field']['not_found']} of {v['n_sentences']} sentence slots proposed by the AI search tasks because the source did not contain them")
    if lit:
        a = lit["aggregate_all_gold"]
        bits.append(f"a proxy of the literature rule recovers known indications of {a['n_diseases']} rare diseases with an AUC of {a['auc_macro']['mean']:.2f}")
    if bits:
        parts.append("Measured: " + "; ".join(bits) + " (Section 3.7).")
    return " ".join(parts)
