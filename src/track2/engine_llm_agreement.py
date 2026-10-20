"""Agreement between independent annotators of the direction of the 139 verified sentences.

    python -m track2.engine_llm_agreement

Sources compared: the original AI labels (`ai`), two blind LLM annotators that saw only the sentence, the drug and the guideline in `protocol/annotation_guideline.md`
(`out/engine/annotations_<model>.json`), and the rule-based extractor of `track2.engine_extractors`. For every pair: accuracy and Cohen's kappa with bootstrap 95 % intervals over papers.
Each sentence also receives a reliability flag from the three semantic annotators (the original AI label and the two blind annotators): `robust` if all three agree, `majority` if two agree,
`unresolved` if all differ. The flag is a property of the label, not a correction of it. Writes `out/engine/llm_agreement.json`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from track2 import engine_extractors as ex

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "out" / "engine"
ANNOTATORS = ("opus", "sonnet")


def load_annotations(engine: Path, n: int) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for name in ANNOTATORS:
        p = engine / f"annotations_{name}.json"
        if not p.exists():
            continue
        labels = {row["id"]: row["label"] for row in json.loads(p.read_text(encoding="utf-8"))["labels"]}
        if sorted(labels) != list(range(n)):
            raise ValueError(f"{p.name}: expected ids 0..{n - 1}")
        bad = set(labels.values()) - set(ex.LABELS)
        if bad:
            raise ValueError(f"{p.name}: unknown labels {sorted(bad)}")
        out[name] = [labels[i] for i in range(n)]
    return out


def reliability(labels: list[str]) -> str:
    counts = sorted((labels.count(c) for c in set(labels)), reverse=True)
    if counts[0] == len(labels):
        return "robust"
    return "majority" if counts[0] >= 2 else "unresolved"


def run(out: Path, engine: Path) -> dict[str, Any]:
    rows = ex.load_entries(out)
    n = len(rows)
    groups = [r["pmid"] for r in rows]
    sources: dict[str, list[str]] = {"ai": [r["ai"] for r in rows], "rule": [ex.rule_direction(r["sentence"], r["drug"]) for r in rows]}
    sources.update(load_annotations(engine, n))
    names = list(sources)
    pairs = {f"{a}_vs_{b}": ex.agreement(sources[a], sources[b], groups) for i, a in enumerate(names) for b in names[i + 1:]}
    semantic = [s for s in ("ai", *ANNOTATORS) if s in sources]
    flags = [reliability([sources[s][i] for s in semantic]) for i in range(n)] if len(semantic) >= 3 else None
    result: dict[str, Any] = {
        "n_sentences": n, "n_papers": len(set(groups)), "annotators": names, "semantic_annotators": semantic,
        "guideline": "protocol/annotation_guideline.md", "blind_annotators_saw": "drug and sentence only, plus the guideline",
        "label_counts": {s: {c: sources[s].count(c) for c in ex.LABELS} for s in names},
        "pairwise_agreement": {k: {m: v for m, v in d.items() if m in ("accuracy", "kappa")} for k, d in pairs.items()},
    }
    if flags:
        result["reliability_counts"] = {f: flags.count(f) for f in ("robust", "majority", "unresolved")}
        by_dir: dict[str, dict[str, int]] = {}
        for r, f in zip(rows, flags):
            by_dir.setdefault(r["ai"], {}).setdefault(f, 0)
            by_dir[r["ai"]][f] += 1
        result["reliability_by_original_label"] = {k: dict(sorted(v.items())) for k, v in sorted(by_dir.items())}
        result["items"] = [{"id": i, "pmid": r["pmid"], "drug": r["drug"], "reliability": f, **{s: sources[s][i] for s in semantic}} for i, (r, f) in enumerate(zip(rows, flags))]
    return result


def main() -> None:
    doc = run(ROOT / "out", ENGINE)
    (ENGINE / "llm_agreement.json").write_text(json.dumps(doc, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    print(doc["annotators"], doc.get("reliability_counts"))
    for k, v in doc["pairwise_agreement"].items():
        print(k, v["accuracy"]["estimate"], v["kappa"]["estimate"], v["kappa"]["ci95"])


if __name__ == "__main__":
    main()
