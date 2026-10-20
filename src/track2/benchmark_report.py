"""Output of the multi-disease benchmark: `benchmark_results.json`, `benchmark_scores.json`, `benchmark_summary.md`."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from track2.http_cache import sha256_file

ADDENDUM = Path("protocol/addendum_benchmark.md")
ADDENDUM_HASH = Path("protocol/addendum_benchmark.sha256")
LABELS = {"700": "STRING 700 (main)", "400": "STRING 400", "900": "STRING 900", "consensus": "Consensus (mean z-rank of 400/700/900)"}


def _public(result: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in result.items() if k not in ("z", "data")}


def _fmt(value: Any, digits: int = 3) -> str:
    return "-" if value is None else f"{value:.{digits}f}"


def _summary_table(summaries: dict[str, Any]) -> list[str]:
    lines = ["| Analysis | Diseases | Mean AUC | 95% CI | Median AUC | AUC > 0.5 | Permutation null mean | Permutation p |", "|---|---|---|---|---|---|---|---|"]
    for key, label in LABELS.items():
        s = summaries.get(key)
        if s is None:
            lines.append(f"| {label} | <2 evaluable | - | - | - | - | - | - |")
            continue
        p = s["permutation"]
        lines.append(
            f"| {label} | {s['n_diseases']} | {_fmt(s['mean_auc'])} | {_fmt(s['ci95_low'])} to {_fmt(s['ci95_high'])} | "
            f"{_fmt(s['median_auc'])} | {s['n_auc_above_0_5']}/{s['n_diseases']} | {_fmt(p['null_mean'])} | {_fmt(p['p_one_sided'], 4)} |"
        )
    return lines


def _disease_table(diseases: list[dict[str, Any]], results: dict[Any, dict[str, dict[str, Any]]]) -> list[str]:
    lines = ["| Disease id | Name | OT targets (score >= 0.5) | Mapped module nodes | Known drugs | AUC 700 | Median pct 700 | AUC 400 | AUC 900 | AUC consensus |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for d in diseases:
        main = results[700][d["efo_id"]]
        cells = [d["efo_id"], d["name"], str(len(d["symbols"])), str(main.get("module_nodes", "-")),
                 str(main["metrics"]["n_known"]) if main["evaluable"] else "-",
                 _fmt(main["metrics"]["auc"]) if main["evaluable"] else "n/e",
                 _fmt(main["metrics"]["median_percentile"]) if main["evaluable"] else "-"]
        for key in (400, 900, "consensus"):
            r = results[key][d["efo_id"]]
            cells.append(_fmt(r["metrics"]["auc"]) if r["evaluable"] else "n/e")
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def write_reports(out_dir: Path, diseases: list[dict[str, Any]], results: dict[Any, dict[str, dict[str, Any]]],
                  summaries: dict[str, Any], decision: dict[str, Any], universe: dict[str, Any], datasources: list[dict[str, Any]]) -> None:
    table = []
    for d in diseases:
        row: dict[str, Any] = {"efo_id": d["efo_id"], "name": d["name"], "n_ot_targets_score_ge_0_5": len(d["symbols"])}
        row["analyses"] = {str(k): _public(results[k][d["efo_id"]]) for k in results}
        table.append(row)
    payload = {
        "addendum_sha256": sha256_file(ADDENDUM), "addendum_hash_file": ADDENDUM_HASH.read_text(encoding="utf-8").split()[0],
        "chembl_version": universe["chembl_version"], "association_datasource_settings": datasources,
        "n_diseases_selected": len(diseases), "diseases": table, "summary": summaries, "decision": decision,
    }
    (out_dir / "benchmark_results.json").write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    scores = {str(k): {eid: {i: round(v, 6) for i, v in r["z"].items()} for eid, r in results[k].items() if r.get("evaluable")}
              for k in results if k != "consensus"}
    (out_dir / "benchmark_scores.json").write_text(json.dumps(scores, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    main = summaries["700"]
    lines = [
        "# Multi-disease benchmark (pre-registered, `protocol/addendum_benchmark.md`)", "",
        f"Addendum SHA-256: `{payload['addendum_sha256']}`. {len(diseases)} diseases selected by the pre-registered rule "
        "(Open Targets approved indications of ChEMBL approved drugs; modules from associated targets with score >= 0.5, "
        "`clinical_precedence` excluded; first diseases in ascending id order).", "",
        "## Decision", "",
        f"Rule: {decision['rule']}.", "",
        f"**Verdict: {decision['verdict']}.**" + ("" if main is None else
            f" Mean AUC {main['mean_auc']:.3f}, 95% CI {main['ci95_low']:.3f} to {main['ci95_high']:.3f}."), "",
        "## Summary across diseases", "", *_summary_table(summaries), "",
        "AUC: known approved-indication drugs vs all other scorable drugs, score = -z. 0.5 = chance. "
        "n/e = not evaluable at that threshold (module < 10 nodes or too few known drugs).", "",
        "## Per disease", "", *_disease_table(diseases, results), "",
        "## Caveats", "",
        "- Diseases share targets and drugs, so they are not independent; the bootstrap treats them as exchangeable.",
        "- Unlabeled drugs may be true indications; this pulls AUC towards 0.5.",
        "- The id-order rule selects mostly `EFO_` terms; the sample is not representative of all diseases.",
        "- Open Targets literature evidence (`europepmc`) can mention known drug targets; this leak is not removed.",
        "- Proximity does not give the direction of effect. These are hypotheses, not medical advice.", "",
    ]
    (out_dir / "benchmark_summary.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
