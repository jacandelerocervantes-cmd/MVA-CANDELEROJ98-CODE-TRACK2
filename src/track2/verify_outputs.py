"""Record or verify SHA-256 checksums of the deterministic outputs.

    python -m track2.verify_outputs --write      # record the current hashes in out/checksums.sha256
    python -m track2.verify_outputs              # compare the files with the record; exit 1 on any difference

Only outputs that depend on the code, the fixed seeds and the cached inputs are listed. Outputs that embed live
service text (DailyMed labels, openFDA counts, PubMed) are reproducible only with the same cache, whose source hashes are in
`out/data_manifest.json`; they are listed under `NETWORK_DEPENDENT` and reported separately, never silently accepted.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

CHECKSUMS = Path("out/checksums.sha256")
DETERMINISTIC = [
    "out/consensus_ranking.json", "out/dti_e2_validation.json", "out/dti_e2_edges.json", "out/dti_e1_module_edges.json",
    "out/validation_results.json",
]
GLOBS = ["out/threshold_scores_sub_*.json", "out/scores_sub_*.json", "out/dti_threshold_scores_*.json", "out/validation_sub_*.json", "out/submission_genes.json",
         # addendum 11: results of the four lines, the class analysis and the ledger (seeded; they depend on data/sources/record.json hashes)
         "out/four_lines_review2.json", "out/class_analysis.json", "out/forking_ledger.json", "out/kg_query_e5*.json", "out/kg_nullquery_e5.json", "out/kg_validation.json", "out/kg_extra.json",
         "out/phenotype_validation*.json", "out/phenotype_proximity_e5.json", "out/phenotype_module_e5.json", "out/module_rule_validation.json", "out/vulnerability_e6.json",
         "out/proposals_review2.json", "out/overlap_check.json", "out/literature_hypotheses/growth_muscle.json", "out/literature_hypotheses/proteotoxic_autophagy.json", "out/literature_hypotheses/senescence_p53.json", "out/literature_hypotheses/tumour_cin_safety.json", "out/literature_hypotheses/context.json",
         "out/signed_effect.json", "out/track1_scenarios.json", "out/differential_diagnosis.json", "out/variant_consequence_classes.json", "out/empirical/label_robustness.json", "out/empirical/controls.json", "out/empirical/variant_calibration.json", "out/empirical/age_support.json", "out/empirical/power_analysis.json", "out/core_n2/*.json", "submissions/*.md", "data/sources/record.json"]
NETWORK_DEPENDENT = ["out/candidates.json", "out/novel_hypotheses_e0_e1_e2.json", "out/plausibility_e0_e1_e2.json", "out/proposals.json", "out/pediatric_arm.json", "out/age_views.json", "out/prism_l4.json", "out/docking_validation.json", "out/candidate_register.json",
                     "out/mva_landscape.json", "out/landscape_fanconi_anaemia.json", "out/phase_inference.json", "out/variant_mechanism.json", "out/case_profile_child.json", "out/case_contrast.json", "out/case_reference_bub1b.json", "out/empirical/literature_control.json", "out/empirical/scalability.json", "out/caution_list.json", "out/age_label_table.json", "out/quote_check.json", "out/label_provenance.json", "out/literature_hypotheses/verification.json", "out/fda_approval_check.json",
                     "out/engine/*.json", "out/engine/llm_labels/*.json"]   # the engine outputs read caches outside the repository; the two label files before them embed live label text (DailyMed, openFDA) read on the date recorded in each file


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def files(root: Path = Path(".")) -> tuple[list[Path], list[Path]]:
    det = {root / p for p in DETERMINISTIC if (root / p).exists()}
    det |= set(root.glob("protocol/*.md")) | set(root.glob("protocol/*.sha256")) | set(root.glob("config/*.json"))  # the rules and the review files, in their final state
    for g in GLOBS:
        det |= set(root.glob(g))
    det = {p for p in det if not p.name.endswith("_track2_report.md")}  # the report is edited by hand until submission and is not a pipeline output
    net = sorted({q for p in NETWORK_DEPENDENT for q in (root.glob(p) if any(c in p for c in "*?[") else [root / p]) if q.exists()})
    return sorted(det), sorted(net)


def write(root: Path = Path(".")) -> int:
    det, net = files(root)
    lines = ["# deterministic"] + [f"{sha256(p)}  {p.relative_to(root).as_posix()}" for p in det] + ["# network-dependent (same cache needed)"] + \
            [f"{sha256(p)}  {p.relative_to(root).as_posix()}" for p in net]
    (root / CHECKSUMS).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"recorded {len(det)} deterministic and {len(net)} network-dependent files")
    return 0


def verify(root: Path = Path(".")) -> int:
    record = (root / CHECKSUMS).read_text(encoding="utf-8").splitlines()
    section, bad, missing, checked = "deterministic", [], [], 0
    for line in record:
        if line.startswith("#"):
            section = "network" if "network" in line else "deterministic"
            continue
        digest, name = line.split("  ", 1)
        path = root / name
        if not path.exists():
            missing.append(name)
        elif sha256(path) != digest:
            bad.append((section, name))
        checked += 1
    det_bad = [n for s, n in bad if s == "deterministic"]
    net_bad = [n for s, n in bad if s == "network"]
    print(f"checked {checked} files; deterministic differences: {len(det_bad)}; network-dependent differences: {len(net_bad)}; missing: {len(missing)}")
    for n in det_bad:
        print("  DIFFERENT (deterministic):", n)
    for n in net_bad:
        print("  different (network-dependent; expected if the cache differs):", n)
    for n in missing:
        print("  MISSING:", n)
    return 1 if det_bad or missing else 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    sys.exit(write() if args.write else verify())


if __name__ == "__main__":
    main()
