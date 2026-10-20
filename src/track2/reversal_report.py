"""Markdown summary of the transcriptomic-reversal layer (every number comes from the computed outputs)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

CAVEATS = [
    "Signatures come from **cancer and immortalised cell lines** (L1000), not from the proband's tissues.",
    "A CRISPR knockout of BUB1B in a cell line is only a **proxy** for the constitutional *hypomorphic* loss in the proband "
    "(partial loss, all tissues, developmental context). The other checkpoint genes are pathway-level proxies only.",
    "A drug whose signature moves against the knockout signature is **not** shown to correct mosaic variegated aneuploidy; "
    "the score is hypothesis-generating, carries no direction of clinical benefit, and is not a dose or treatment statement.",
    "Consensus chemical signatures merge cell lines, doses and times; cell-line detail comes from a second, descriptive "
    "query of the individual-signature library restricted to its top reversers (absence there is not absence of effect).",
    "Drugs are matched to LINCS perturbagens by cleaned name only; ambiguous matches are flagged `name_match: ambiguous`.",
]


def write_summary(path: Path, scores: dict[str, Any], meta: dict[str, Any]) -> None:
    avail = {d: v for d, v in scores.items() if v["available"]}
    top = sorted(avail.items(), key=lambda kv: (-kv[1]["score"], kv[0]))[:25]
    nul = meta["null"]
    lines = ["# Transcriptomic reversal layer (SigCom LINCS, L1000)", "", "> Hypothesis-generating evidence. Not medical advice.", "",
             "## Caveats", ""] + [f"- {c}" for c in CAVEATS] + ["", "## Method", "",
             f"- Query: top {meta['top_n']} up and {meta['top_n']} down genes of the LINCS CRISPR-KO consensus signature of "
             "BUB1B (primary) and of " + ", ".join(k for k in meta["query_sets"] if k != "BUB1B") + " (proxies).",
             "- Ranking: SigCom LINCS `POST /data-api/api/v1/enrich/ranktwosided` against `l1000_mean_cp` (chemical consensus "
             "signatures). Score = best reversal z-sum (-z-sum, larger = stronger reversal) among the signatures of a drug.",
             f"- Null: {nul['n_permutations']} random up/down gene sets of the same size (seed {meta['seed']}) through the same "
             f"pipeline; empirical p pools {nul['n_null_values']} drug-level null scores; FDR is Benjamini-Hochberg across scored drugs.",
             "- Direction: `reverses` if FDR <= 0.05 for reversal, `mimics` if FDR <= 0.05 for mimicry, else `none`.", "",
             f"- Structured control: the same pipeline on {len(meta['ko_control_genes'])} unrelated CRISPR-KO consensus signatures "
             f"({', '.join(meta['ko_control_genes'])}); `ko_control.p` = (1 + controls scoring >= BUB1B) / (1 + {len(meta['ko_control_genes'])}). "
             "A drug that reverses most knockout signatures is a non-specific reverser.", "",
             "## Results", "",
             f"- Approved drugs in the universe: {meta['n_universe_drugs']}; with a LINCS chemical consensus signature (scored): "
             f"{meta['n_drugs_with_score']}; the rest are `available: false`.",
             f"- `reverses` at FDR <= 0.05: {meta['n_direction_reverses']}; `mimics`: {meta['n_direction_mimics']}.",
             f"- Null reversal score quantiles (random sets): {nul['null_reversal_score_quantiles']}; maximum {nul['null_max_reversal_score']}.",
             f"- Drugs `reverses` AND ko_control p <= 0.05 (specific to BUB1B-like signature): {meta['n_reverses_and_specific']}.",
             f"- Random sets called `reverses` per set (leave-one-out, same rule): mean {nul['mean_false_calls_per_random_set']} "
             f"(per set: {nul['false_reverses_calls_per_random_set_leave_one_out']}).", "",
             "## Top 25 reversers", "",
             "| Drug | Score | p | FDR | Direction | Signatures | Cell lines (top-reverser layer) | In candidates | Proxy sets reversing | p vs unrelated KO controls |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for _, v in top:
        px = sum(p["direction"] == "reverses" for p in v["proxy_sac"].values())
        lines.append(f"| {v['name']} | {v['score']} | {v['p']} | {v['fdr']} | {v['direction']} | {v['n_signatures']} | "
                     f"{', '.join(v['cell_lines']) or '-'} | {v['in_candidates']} | {px}/{len(v['proxy_sac'])} | {v['ko_control']['p']} |")
    lines += ["", "Full table: `reversal_scores.json`; gene lists and parameters: `reversal_meta.json`.", ""]
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")
