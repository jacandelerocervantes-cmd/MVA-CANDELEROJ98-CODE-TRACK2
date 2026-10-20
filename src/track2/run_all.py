"""One command to reproduce the whole analysis, stage by stage, on Windows, Linux and macOS.

    python -m track2.run_all --list                      # what the stages are and roughly how long they take
    python -m track2.run_all --dry-run                   # print the commands without running them
    python -m track2.run_all                             # run every stage in order
    python -m track2.run_all --from-stage 8              # resume from a stage (stages are independent once their inputs exist)
    python -m track2.run_all --only 12 13                # run just these stages
    python -m track2.run_all --skip 3 4                  # skip the slow validation stages

Run it from the repository root with `PYTHONPATH=src` (or an editable install). The first run downloads ChEMBL, STRING, Reactome, DailyMed,
openFDA and PubMed responses into `data/cache/`; later runs reuse them. Every stage is a module with a fixed seed, so the same inputs give
the same outputs; `track2.verify_outputs` checks that against `out/checksums.sha256`. The patient VCF is not distributed, so the genome
covariate and pharmacogenomic genotyping stages are not part of this script (see the README).
"""

from __future__ import annotations

import os
import argparse
import subprocess
import sys
import time
from typing import NamedTuple

PY = [sys.executable, "-m"]


class Stage(NamedTuple):
    number: int
    name: str
    minutes: str
    commands: list[list[str]]


def stages() -> list[Stage]:
    m = "config/mechanisms_from_submission"
    return [
        Stage(1, "Gene modules and proximity from the Track 1 submission (STRING 400/700/900)", "10-40 (first run: +15 for ChEMBL)",
              [PY + ["track2.from_submission", "--submission", os.environ.get("TRACK1_SUBMISSION", "inputs/track1_submission.csv"), "--out", "out"]]),
        Stage(2, "Hand-defined comparison modules and the Fanconi control (needed by the frozen validation)", "5-15",
              [PY + ["track2.run_pipeline", "--mechanism", f"config/mechanisms/{x}.json", "--out", "out/comparison_manual"] for x in
               ("bub1b_sac_lof", "bub1b_down_proteotoxic", "bub1b_down_p53_senescence", "bub1b_down_growth_mtor_ampk")]
              + [PY + ["track2.run_pipeline", "--mechanism", "config/mechanisms/fanconi_anemia_demo.json", "--out", "out"]]),
        Stage(3, "Frozen validation protocol P1-P4 (hand-defined checkpoint module)", "30-60", [PY + ["track2.validation", "--out", "out"]]),
        Stage(4, "Validation P1/P3/P4 on each submission-derived module", "30-90 each (run in parallel by hand if you wish)",
              [PY + ["track2.validate_submission_modules", "--all", "--mechanism-dir", m, "--out", "out"]]),
        Stage(5, "Candidate records, schema check, PMID verification", "3-10", [PY + ["track2.build_candidates", "--out", "out"]]),
        Stage(6, "Evidence layers: literature and trials, pharmacovigilance, direction of effect", "15-30",
              [PY + ["track2.evidence"], PY + ["track2.faers_disproportionality"],
               PY + ["track2.direction", "--universe", "data/universe.json", "--candidates", "out/candidates.json", "--out", "out"]]),
        Stage(7, "Measured properties, affinity and ADME from ChEMBL", "5-10", [PY + ["track2.physchem"]]),
        Stage(8, "Consensus across thresholds and weighting schemes", "<1", [PY + ["track2.consensus", "--out", "out"]]),
        Stage(9, "Curated literature context and the extrapolation arm", "5-10",
              [PY + ["track2.literature_context", "--candidates", "out/candidates.json", "--out", "out"],
               PY + ["track2.build_candidates", "--out", "out/comparison_manual", "--mechanism-dir", "config/mechanisms_extrapolation", "--dest", "out/arm_extrapolation"],
               PY + ["track2.literature_context", "--candidates", "out/arm_extrapolation/candidates.json", "--out", "out", "--name", "literature_context_extrapolation.json"]]),
        Stage(10, "Paediatric-labelled search", "5-15", [PY + ["track2.pediatric_arm", "--out", "out"]]),
        Stage(11, "Expanded drug-target links: measured (E1) and predicted (E2) binding", "40-90 (ChEMBL calls)",
              [PY + ["track2.dti_edges"], PY + ["track2.dti_similarity", "--out", "out"]]),
        Stage(12, "Proximity with the expanded links, and the sensitivity analysis without generic complexes", "30-60",
              [PY + ["track2.augmented_proximity", "--variant", "e0_e1", "--out", "out"],
               PY + ["track2.augmented_proximity", "--variant", "e0_e1_e2", "--out", "out"],
               PY + ["track2.core_sensitivity", "--out", "out"]]),
        Stage(13, "New hypotheses, FDA scope, plausibility tiers, proposals and age groups", "5-15",
              [PY + ["track2.novel_hypotheses", "--variant", v, "--out", "out"] for v in ("e0_e1", "e0_e1_e2")]
              + [PY + ["track2.plausibility", "--variant", v, "--out", "out"] for v in ("e0_e1", "e0_e1_e2")]
              + [PY + ["track2.novel_hypotheses", "--variant", "e0_e1_e2", "--out", "out", "--top-n", "50", "--pool-min", "7", "--tag", "_wide"],
                 PY + ["track2.plausibility", "--variant", "e0_e1_e2", "--out", "out", "--tag", "_wide"]]
              + [PY + ["track2.proposals", "--variant", "e0_e1_e2", "--out", "out"], PY + ["track2.age_view", "--out", "out"]]),
        Stage(14, "Independent functional evidence L4: DepMap expression and PRISM drug sensitivity (downloads about 540 MB, CC BY 4.0)", "5-15",
              [PY + ["track2.depmap_download"], PY + ["track2.prism_l4", "--out", "out"], PY + ["track2.proposals", "--variant", "e0_e1_e2", "--out", "out"]]),
        Stage(15, "Docking (E3): Vina binary, PDB structures, validation per target, screen and top candidates", "120-400 (8 cores)",
              [PY + ["track2.get_vina"], PY + ["track2.docking_structures"], PY + ["track2.docking", "--stage", "validate"], PY + ["track2.docking", "--stage", "screen"],
               PY + ["track2.docking", "--stage", "top"]]),
        Stage(16, "Integrate, review page, data snapshot (the addendum-10 report generator `track2.build_report` is no longer run: the report comes from stage 21)", "2-5",
              [PY + ["track2.integrate", "--out", "out"], PY + ["track2.review_page", "--out", "out"], PY + ["track2.snapshot", "--out", "out"],
               PY + ["track2.candidate_register", "--out", "out"]]),
        Stage(17, "Sources of the four lines (HPO, MONDO, PrimeKG, DepMap CRISPR and covariates), SHA-256 recorded in data/sources/record.json (about 1.1 GB, one source at a time)", "10-30",
              [PY + ["track2.fetch_sources", "--only", s] for s in ("hpo_annotations", "hpo_ontology", "mondo", "primekg", "depmap_crispr", "depmap_covariates")]),
        Stage(18, "Line A (knowledge graph): validation, extra checks, patient query without the parental term (E5), specificity against 500 random queries (E4)", "30-90",
              [PY + ["track2.kg_validation", "--out", "out"], PY + ["track2.kg_extra", "--out", "out"], PY + ["track2.kg_query", "--out", "out"], PY + ["track2.kg_query", "--out", "out", "--e5"],
               PY + ["track2.kg_nullquery", "--out", "out", "--e5", "--queries", "500"]]
              + [PY + ["track2.kg_query", "--out", "out", "--e5", "--leave-out", t] for t in ("HP:0000121", "HP:0001508", "HP:0001518", "HP:0001622", "HP:0003202", "HP:0004322")]
              + [PY + ["track2.kg_query", "--out", "out", "--e5", "--only", o] for o in ("core", "nonspecific")]),
        Stage(19, "Line B (proximity): phenotype module (E5), validation with and without the causal genes (E2/E3), proximity with two-record links (E8), BUB1B module rule (E1)", "60-180",
              [PY + ["track2.phenotype_module"], PY + ["track2.phenotype_validation", "--out", "out"], PY + ["track2.phenotype_validation", "--out", "out", "--leave-causal"],
               PY + ["track2.phenotype_proximity", "--module", "out/phenotype_module_e5.json", "--out", "out", "--strict-links", "--tag", "_e5"],
               PY + ["track2.core_n2", "--out", "out/core_n2"], PY + ["track2.track1_scenarios", "--out", "out"], PY + ["track2.module_rule_validation", "--out", "out", "--workers", "2"]]),
        Stage(20, "Lines C and D: vulnerability with TP53, MSI and CIN and 1,000 null sets (E6); signed direction (exploratory)", "30-90",
              [PY + ["track2.vulnerability_scan", "--out", "out", "--null-sets", "1000", "--workers", "4", "--covariates", "data/sources/depmap_covariates"], PY + ["track2.signed_effect", "--out", "out"]]),
        Stage(21, "Combine the lines (E1, E4, E7, E9), class analysis (E10), caution list (E11), age tables (E15), forking ledger (E12), report, methods, form and PDF", "20-40",
              [PY + ["track2.four_lines_review2", "--out", "out"], PY + ["track2.class_analysis", "--out", "out", "--permutations", "10000"], PY + ["track2.caution_list"],
               PY + ["track2.fda_approval_check"], PY + ["track2.verify_literature"], PY + ["track2.literature_summary"], PY + ["track2.overlap_check"],
               PY + ["track2.age_label_table"], PY + ["track2.verify_quotes"], PY + ["track2.label_provenance"], PY + ["track2.forking_ledger", "--out", "out"], PY + ["track2.build_article", "--out", "out"],
               PY + ["track2.build_methods_review2", "--out", "out"], PY + ["track2.fill_methods_form"],
               PY + ["track2.md_to_pdf", "submissions/candeleroj98_track2_report.md", "submissions/candeleroj98_track2_report.pdf"]]),
        Stage(22, "Verify the outputs against the recorded checksums", "<1", [PY + ["track2.verify_outputs"]]),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--list", action="store_true", help="list the stages and exit")
    parser.add_argument("--dry-run", action="store_true", help="print the commands without running them")
    parser.add_argument("--from-stage", type=int, default=1)
    parser.add_argument("--only", type=int, nargs="*")
    parser.add_argument("--skip", type=int, nargs="*", default=[])
    args = parser.parse_args()
    table = stages()
    if args.list:
        for s in table:
            print(f"{s.number:>2}. {s.name}  [~{s.minutes} min]")
        return
    chosen = [s for s in table if (s.number in args.only if args.only else s.number >= args.from_stage) and s.number not in args.skip]
    numbers = {s.number for s in chosen}
    if not args.dry_run:
        if 1 in numbers:
            csv_path = os.environ.get("TRACK1_SUBMISSION", "inputs/track1_submission.csv")
            if not os.path.isfile(csv_path):
                sys.exit(f"Stage 1 needs the Track 1 submission CSV, which this repository does not include ({csv_path} not found).\n"
                         "Take it from the Track 1 submission and set the environment variable TRACK1_SUBMISSION to its path, or start from a later stage with --from-stage.")
        if numbers & {18, 19, 20, 21} and 17 not in numbers and not os.path.isfile("data/sources/record.json"):
            sys.exit("Stages 18 to 21 need the public source files downloaded by stage 17 (data/sources/ is empty or missing).\n"
                     "Run stage 17 first: python -m track2.run_all --only 17")
    for s in chosen:
        print(f"\n=== Stage {s.number}: {s.name} (~{s.minutes} min)", flush=True)
        for cmd in s.commands:
            print("$ python", " ".join(cmd[1:]), flush=True)  # the interpreter path is not printed (it would expose the user's home directory)
            if args.dry_run:
                continue
            start = time.time()
            if subprocess.run(cmd, check=False).returncode != 0:
                sys.exit(f"Stage {s.number} failed on: {' '.join(cmd)}\nFix the cause and resume with --from-stage {s.number}.")
            print(f"  done in {time.time() - start:.0f} s", flush=True)


if __name__ == "__main__":
    main()
