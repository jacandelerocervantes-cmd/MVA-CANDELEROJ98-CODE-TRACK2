"""Empirical scalability harness: rerun the reusable early stages on synthetic Track 1-format cases and measure them.

    python -m track2.empirical_scalability --list
    python -m track2.empirical_scalability --dry-run
    python -m track2.empirical_scalability --snapshot out/empirical/snapshot_before.json     # hashes of out/ and config/
    python -m track2.empirical_scalability                                                  # run all cases twice, write out/empirical/scalability.json

Synthetic inputs go to data/empirical_cases/<case>.csv (positions are made up and belong to no patient). Every run writes only to
out_case_<case>_run<k>/ and config_case_<case>_run<k>/; the real case's out/ and config/ are never written (a guard refuses those
paths). Reusable stages measured: stage 1 (track2.from_submission: gene weights, module building, proximity at STRING 400/700/900)
and stage 8 (track2.consensus). Timing and memory fields are machine-dependent.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any

try:
    import psutil
except ImportError:  # timing still works; memory is then reported as null
    psutil = None

CASES: dict[str, dict[str, Any]] = {
    "fanconi_anemia": {"description": "Fanconi anaemia (public knowledge: FANCA, ATR); synthetic positions",
                       "rows": [("FANCA", "16", 89000001, "A", "G", 0.91, "primary"), ("ATR", "3", 142000001, "C", "T", 0.62, "secondary")]},
    "mva_cep57_trip13": {"description": "Mosaic variegated aneuploidy-type gene set (CEP57, TRIP13); synthetic positions",
                         "rows": [("CEP57", "11", 95000001, "G", "A", 0.88, "primary"), ("TRIP13", "5", 900001, "T", "C", 0.57, "secondary")]},
}
SEED_NOTE = "ITERATIONS and the permutation seed are fixed inside track2.run_pipeline; the repeat run uses the same code, so it tests determinism"
PROTECTED = ("out", "config")


def write_case_csv(path: Path, rows: list[tuple]) -> None:
    """Track 1-format CSV (same columns as submission_input.read_submission); the gene is the first word of notes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        w = csv.writer(handle, lineterminator="\n")
        w.writerow(["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"])
        for gene, chrom, pos, ref, alt, epcr, ftype in rows:
            w.writerow(["SYNTHETIC_" + gene, "chr" + chrom, pos, ref, alt, "", "", "", "", epcr, ftype, f"{gene} synthetic_position_not_a_patient"])


def guard(path: Path, root: Path) -> None:
    resolved = path.resolve()
    for name in PROTECTED:
        protected = (root / name).resolve()
        if resolved == protected or protected in resolved.parents:
            raise SystemExit(f"Refusing to write inside the real case's {name}/: {path}")


def hash_tree(root: Path, exclude: tuple[str, ...] = ()) -> dict[str, str]:
    out = {}
    for p in sorted(root.rglob("*")):
        if p.is_file() and "__pycache__" not in p.parts and not any(e in p.relative_to(root).parts[:2] for e in exclude):
            out[p.relative_to(root).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def snapshot(root: Path) -> dict[str, str]:
    snap: dict[str, str] = {}
    for name in PROTECTED:
        if (root / name).is_dir():
            snap.update({f"{name}/{k}": v for k, v in hash_tree(root / name, exclude=("empirical",)).items()})
    return snap


def run_timed(cmd: list[str], env: dict[str, str] | None = None, poll: float = 0.05) -> dict[str, Any]:
    """Run a command; return wall seconds, peak resident memory (MB, process tree; null without psutil), return code and output tail."""
    start = time.perf_counter()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, encoding="utf-8", errors="replace")
    peak = 0
    if psutil is not None:
        try:
            ps = psutil.Process(proc.pid)
        except psutil.Error:
            ps = None
        chunks: list[str] = []
        reader = threading.Thread(target=lambda: chunks.append(proc.stdout.read()), daemon=True)
        reader.start()
        while proc.poll() is None:
            try:
                procs = [ps] + ps.children(recursive=True) if ps else []
                peak = max(peak, sum(p.memory_info().rss for p in procs))
            except psutil.Error:
                pass
            time.sleep(poll)
        reader.join()
        output = "".join(chunks)
    else:
        output = proc.communicate()[0]
    wall = time.perf_counter() - start
    return {"returncode": proc.returncode, "wall_seconds": round(wall, 3), "peak_rss_mb": round(peak / 1e6, 1) if psutil is not None else None,
            "output_tail": output.strip().splitlines()[-5:]}


def stage_commands(csv_path: Path, out: Path, cfg: Path) -> list[tuple[str, list[str]]]:
    py = [sys.executable, "-m"]
    return [("stage1_from_submission", py + ["track2.from_submission", "--submission", str(csv_path), "--out", str(out), "--config-dir", str(cfg)]),
            ("stage8_consensus", py + ["track2.consensus", "--out", str(out), "--mechanism-dir", str(cfg)])]


def read_counts(out: Path, cfg: Path) -> dict[str, Any]:
    genes = json.loads((out / "submission_genes.json").read_text(encoding="utf-8"))
    modules = sorted(p.stem for p in cfg.glob("*.json"))
    scored = {}
    for p in sorted(out.glob("manifest_*.json")):
        m = json.loads(p.read_text(encoding="utf-8"))
        scored[m["mechanism"]] = {"n_scored": {k: v["n_scored"] for k, v in sorted(m["sensitivity"].items())}, "n_pass": {k: v["n_pass"] for k, v in sorted(m["sensitivity"].items())}}
    consensus = out / "consensus_ranking.json"
    return {"n_genes": len(genes["genes"]), "n_modellable_modules": genes["n_modellable"], "modules": modules,
            "module_size_in_network": {g["gene"]: g["module_size_in_network"] for g in genes["genes"]}, "drugs_scored_per_module": scored,
            "n_drugs_in_consensus": len(json.loads(consensus.read_text(encoding="utf-8"))["drugs"]) if consensus.exists() else None}


def run_case(name: str, root: Path, env: dict[str, str], repeats: int = 2) -> dict[str, Any]:
    csv_path = root / "data" / "empirical_cases" / f"{name}.csv"
    write_case_csv(csv_path, CASES[name]["rows"])
    runs, hashes = [], []
    for k in range(1, repeats + 1):
        out, cfg = root / f"out_case_{name}_run{k}", root / f"config_case_{name}_run{k}"
        for d in (out, cfg):
            guard(d, root)
        stages = {}
        for label, cmd in stage_commands(csv_path, out, cfg):
            stages[label] = run_timed(cmd, env)
            if stages[label]["returncode"] != 0:
                break
        runs.append(stages)
        hashes.append({**hash_tree(out), **{f"cfg/{k2}": v for k2, v in hash_tree(cfg).items()}} if out.is_dir() else {})
    ok = all(s["returncode"] == 0 for st in runs for s in st.values())
    # submission_genes.json embeds the input hash only, so it is comparable across runs
    counts = read_counts(root / f"out_case_{name}_run1", root / f"config_case_{name}_run1") if ok else None
    return {"description": CASES[name]["description"], "input_csv": csv_path.relative_to(root).as_posix(),
            "input_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(), "counts": counts, "runs_machine_dependent": runs,
            "n_output_files": len(hashes[0]), "byte_identical_between_repeats": bool(ok and hashes[0] and all(h == hashes[0] for h in hashes[1:])),
            "differing_files": sorted(f for f in set(hashes[0]) | set(hashes[-1]) if hashes[0].get(f) != hashes[-1].get(f)) if ok else None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--snapshot", type=Path, help="write the SHA-256 snapshot of out/ and config/ to this file and exit")
    parser.add_argument("--cases", nargs="*", default=list(CASES))
    parser.add_argument("--result", type=Path, default=Path("out/empirical/scalability.json"))
    args = parser.parse_args()
    root = Path.cwd()
    env = {**os.environ, "PYTHONPATH": "src"}
    if args.list:
        for n, c in CASES.items():
            print(n, "-", c["description"])
        return
    if args.snapshot:
        args.snapshot.parent.mkdir(parents=True, exist_ok=True)
        snap = snapshot(root)
        args.snapshot.write_text(json.dumps(snap, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        print(len(snap), "files hashed")
        return
    if args.dry_run:
        for n in args.cases:
            for k in (1, 2):
                for label, cmd in stage_commands(Path(f"data/empirical_cases/{n}.csv"), Path(f"out_case_{n}_run{k}"), Path(f"config_case_{n}_run{k}")):
                    print(f"[{n} run{k}] {label}: python", " ".join(cmd[1:]))
        return
    before = snapshot(root)
    result = {"cases": {n: run_case(n, root, env) for n in args.cases}, "machine": {"platform": platform.platform(), "python": platform.python_version(),
              "psutil": psutil.__version__ if psutil else None, "cpu_count": os.cpu_count()},
              "machine_dependent_fields": ["runs_machine_dependent.*.wall_seconds", "runs_machine_dependent.*.peak_rss_mb", "machine"],
              "not_run_stages": {"4": "validate_submission_modules: 30-90 min per module", "5": "build_candidates reads config/patient.json (patient-specific gates): needs a manual patient file",
                                 "18-21": "rewrite the real case's outputs; never run here"},
              "determinism_note": SEED_NOTE, "real_case_out_config_unchanged": snapshot(root) == before}
    guard(args.result, Path("/nonexistent")) if False else None
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.result.write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({n: {"identical": c["byte_identical_between_repeats"], "counts": c["counts"]} for n, c in result["cases"].items()}, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
