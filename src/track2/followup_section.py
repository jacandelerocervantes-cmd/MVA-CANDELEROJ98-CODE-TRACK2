"""Section 3.15 (which hypotheses can be refuted by computation and which need an experiment) and supplement S15 (follow-up experiments).

Everything numeric or listed is read from `out/experiment_plan.json` (and `out/line_c_adjusted.json`); only wording is typed here.
The classification rule is fixed in `protocol/addendum_refutability.md`; the designs are in `protocol/follow_up_experiments_v2.md`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_FILES = ("experiment_plan", "line_c_adjusted")


def _load(out: Path, name: str) -> Any:
    p = Path(out) / f"{name}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def _have(out: Path) -> bool:
    return all((Path(out) / f"{n}.json").exists() for n in _FILES)


def _cell(x: Any, geo: dict[str, Any] | None = None) -> str:
    """Table cell; the plan quotes the patient-signature arrays as '2 vs 2' while the analysis has 2 case and 4 control arrays, so the counts are re-read from the analysis."""
    text = str(x).replace("|", "/").replace("\n", " ")
    if geo:
        text = text.replace("2 vs 2 arrays", f"{geo['n_cases']} vs {geo['n_controls']} arrays").replace("2 versus 2 arrays", f"{geo['n_cases']} versus {geo['n_controls']} arrays")
    return text


def abstract_sentence(out: Path) -> str:
    out = Path(out)
    if not _have(out):
        return ""
    e = _load(out, "experiment_plan")
    n_cf = len(e["computational_followups"])
    n_run = sum(bool(x["run"]) for x in e["computational_followups"])
    return (f"One computational refutation was done (the adjustment of line C) and {n_cf - n_run} more checks are defined and not run; the hypotheses that computation cannot reach need the experiments of Section 3.15 and Supplement S15, "
            "beginning with a karyotype or single-cell sequencing (Step 0).")


def main_section(out: Path) -> str:
    out = Path(out)
    if not _have(out):
        return ""
    e = _load(out, "experiment_plan")
    cnt = e["counts"]
    t = ["### 3.15 Which hypotheses can be refuted by computation and which need an experiment\n"]
    t.append("The rule was written before the table (`protocol/addendum_refutability.md`). A hypothesis is computationally refutable if public data already in the project can measure the gene, pathway or cell-line effect it states; it needs an experiment if it requires cells with partial BUB1B loss, a model of the case or a clinical outcome. "
             "'Refuted' is kept for a line that shows an effect in the opposite direction; a failed check closes a row, and a passed check only means not yet refuted.\n")
    t.append("**Done.** " + " ".join(f"{x['what']}: {x['result']}." for x in e["computational_done"]) + "\n")
    t.append("**Defined and not run (computational follow-ups).**\n")
    t.append("| ID | Check | Closes the row if | Run |\n|---|---|---|---|")
    for x in e["computational_followups"]:
        t.append(f"| {x['id']} | {_cell(x['name'])} | {_cell(x['refutation'])} | {'yes' if x['run'] else 'no'} |")
    t.append("")
    t.append("**Needs a laboratory experiment (ordered by information value).**\n")
    t.append("| ID | Candidates | Refuted if |\n|---|---|---|")
    for x in sorted(e["experiments"], key=lambda z: z["order"]):
        t.append(f"| {x['id']} | {_cell(', '.join(x['candidates']))} | {_cell(x['refuted_if'])} |")
    t.append("")
    t.append(f"Of the {cnt['triage_rows']} candidates from all lines, {cnt['comp_first']} go to a computational check first, {cnt['lab_design_benefit']} to a laboratory design for benefit, {cnt['hazard_class']} to a laboratory design for harm (cytotoxic classes) and {cnt['closed']} are closed by the adjustment (the full triage is in Supplement S15). "
             "The computational checks can say that a drug is not favoured by public data and can raise a hazard flag; they cannot say anything about cells with partial BUB1B loss, about a child of unknown age, or about benefit and harm in a person. "
             "The experiments can test that, but only after Step 0 (E0) confirms aneuploidy or an impaired missense variant, and a null result with low power would not refute a hypothesis. Neither a computational nor an experimental pass would justify use in a person.\n")
    return "\n".join(t) + "\n"


def supplement(out: Path) -> str:
    out = Path(out)
    if not _have(out):
        return ""
    e = _load(out, "experiment_plan")
    geo = (_load(out, "line_e_signature_reversal") or {}).get("geo")
    t = ["## S15. Follow-up experiments\n"]
    t.append("Condensed from `protocol/follow_up_experiments_v2.md` and `out/experiment_plan.json`. In vitro only; no dose, schedule or recommendation for any person; no experiment has been run. "
             "Model for every drug: an isogenic human line with partial BUB1B reduction (two clones), its euploid parent and a rescue line; patient-derived cells only with consent.\n")
    t.append("**Triage** (age readings are in Section 3.13).\n")
    t.append("| Candidate | Lines | Computational status | Deciding uncertainty | Next step |\n|---|---|---|---|---|")
    for r in e["triage"]:
        t.append(f"| {_cell(r['candidate'])} | {_cell(', '.join(r['lines']))} | {_cell(r['computational_status'], geo)} | {_cell(r['deciding_uncertainty'], geo)} | {_cell(r['next_step'], geo)} |")
    t.append("")
    reps: dict[str, list[str]] = {}
    for x in sorted(e["experiments"], key=lambda z: z["order"]):
        reps.setdefault(_cell(x["replicates"], geo), []).append(x["id"])
    t.append("**Replicates.** " + " ".join(f"{', '.join(ids)}: {r}." for r, ids in reps.items()) + "\n")
    t.append("**Stop-go order.** " + "; ".join(f"{s['stage']}: go if {s['go'].lower()}, stop if {s['stop'].lower()}" for s in e["stop_go"]) + ". A drug that survives a stage is not accepted, only not yet refuted; a drug refuted in two clones is recorded as a documented dead end.\n")
    return "\n".join(t) + "\n"
