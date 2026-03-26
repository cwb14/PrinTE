#!/usr/bin/env python3
import argparse
import itertools
import math
import os
import random
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional


def dir_name_from_combo(ir, dr, sr, lb):
    return f"insertion_rates_{ir}_deletion_rates_{dr}_solo_ratio_{sr}_length_bias_{lb}"


def format_sci(x: float) -> str:
    """
    Format floats like 1e-07 -> '1e-7', 2.5e-12 -> '2.5e-12'
    so directory names stay readable/stable.
    """
    s = f"{x:.12g}"
    return s.replace("E", "e")


def build_logspace_values(start: float, end: float, count: int) -> List[str]:
    """
    Build `count` values log10-spaced between start and end (inclusive).
    Handles start > end or start < end.
    """
    if count < 2:
        raise ValueError("--*-count must be >= 2 to define a range.")
    if start <= 0 or end <= 0:
        raise ValueError("--*-start and --*-end must be positive for log spacing.")

    log_s = math.log10(start)
    log_e = math.log10(end)

    values = []
    for i in range(count):
        t = i / (count - 1)
        log_v = log_s + (log_e - log_s) * t
        v = 10 ** log_v
        values.append(format_sci(v))
    return values


def build_param_grids(args):
    """
    insertion_rates: log-spaced between ins_start and ins_end (count = ins_count)
    deletion_rates:  log-spaced between del_start and del_end (count = del_count)
    solo_ratio:      int range
    length_bias:     int range
    """
    insertion_rates = build_logspace_values(args.ins_start, args.ins_end, args.ins_count)
    deletion_rates = build_logspace_values(args.del_start, args.del_end, args.del_count)

    solo_ratios = list(range(args.sr_start, args.sr_end + 1, args.sr_step))
    length_biases = list(range(args.k_start, args.k_end + 1, args.k_step))

    return insertion_rates, deletion_rates, solo_ratios, length_biases


def build_all_combinations(insertion_rates, deletion_rates, solo_ratios, length_biases):
    return list(itertools.product(insertion_rates, deletion_rates, solo_ratios, length_biases))


def validate_inputs(args):
    """
    Light validation: ensure required script/files exist (as typed),
    so tab-completion helps users catch mistakes early.
    """
    for label, p, must_exist in [
        ("--printe-script", args.printe_script, True),
        ("--bed", args.bed, True),
        ("--fasta", args.fasta, True),
        ("--te-lib", args.te_lib, True),
        ("--ratios", args.ratios, True),
    ]:
        pp = Path(p).expanduser()
        if must_exist and not pp.exists():
            raise FileNotFoundError(f"{label} path does not exist: {p}")

    if args.ge <= 0 or args.st <= 0:
        raise ValueError("--ge and --st must be > 0")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.tstv <= 0:
        raise ValueError("--tstv must be > 0")


def relpath_from(workdir: Path, user_path: str, base_cwd: Path) -> str:
    """
    Convert a user-supplied path (relative to base_cwd) into a path usable from workdir.
    - If user_path is absolute, return it.
    - Else treat it as relative to base_cwd (where script was launched),
      then compute relative path from workdir.
    """
    p = Path(user_path).expanduser()
    if p.is_absolute():
        return str(p)

    abs_target = (base_cwd / p).resolve()
    return os.path.relpath(str(abs_target), str(workdir.resolve()))


def build_printe_cmd(
    combo: Tuple[str, str, int, int],
    args,
    workdir: Path,
    base_cwd: Path
) -> List[str]:
    ir, dr, sr, lb = combo

    printe_script = relpath_from(workdir, args.printe_script, base_cwd)
    te_lib = relpath_from(workdir, args.te_lib, base_cwd)
    bed = relpath_from(workdir, args.bed, base_cwd)
    fasta = relpath_from(workdir, args.fasta, base_cwd)
    ratios = relpath_from(workdir, args.ratios, base_cwd)

    cmd = [
        "bash", printe_script,
        "-ge", str(args.ge),
        "-st", str(args.st),
        "-F", f"{ir},{dr}",
        "-m", str(args.mut),
        "-mxgs", str(args.mxgs),
        "-mngs", str(args.mngs),
        "-sr", str(sr),
        "-k", str(lb),
        "-TsTv", str(args.tstv),
        "--ex_LTR",
        "--no_postproc",
        "--TE_lib", str(te_lib),
        "-t", str(args.threads),
        "--bed", str(bed),
        "--fasta", str(fasta),
        "-r", str(ratios),
    ]
    return cmd


def launch_local_job(combo, args, base_cwd: Path):
    ir, dr, sr, lb = combo
    dir_name = dir_name_from_combo(ir, dr, sr, lb)
    workdir = Path(dir_name)
    workdir.mkdir(exist_ok=True)

    cmd = build_printe_cmd(combo, args, workdir, base_cwd)

    print(f"Starting combo in {dir_name}:")
    print("  " + " ".join(cmd))

    proc = subprocess.Popen(cmd, cwd=str(workdir))
    return proc, dir_name


def wait_for_slot(active, max_procs):
    while len(active) >= max_procs:
        for proc, name in list(active):
            ret = proc.poll()
            if ret is not None:
                print(f"[DONE] {name} (exit code {ret})")
                active.remove((proc, name))
                break
        else:
            proc, name = active[0]
            print(f"[WAIT] Waiting for {name} to finish...")
            ret = proc.wait()
            print(f"[DONE] {name} (exit code {ret})")
            active.remove((proc, name))


def wait_for_all(active):
    for proc, name in active:
        print(f"[FINAL WAIT] Waiting for {name} to finish...")
        ret = proc.wait()
        print(f"[DONE] {name} (exit code {ret})")
    active.clear()


def prepare_combos(args) -> Tuple[List[Tuple[str, str, int, int]], int, int]:
    """
    Returns:
      available_combos (dirs not already present),
      total_grid_count,
      skipped_existing_count
    """
    insertion_rates, deletion_rates, solo_ratios, length_biases = build_param_grids(args)
    combos = build_all_combinations(insertion_rates, deletion_rates, solo_ratios, length_biases)

    total = len(combos)
    available = []
    skipped = 0
    for ir, dr, sr, lb in combos:
        dname = dir_name_from_combo(ir, dr, sr, lb)
        if Path(dname).exists():
            skipped += 1
        else:
            available.append((ir, dr, sr, lb))
    return available, total, skipped


def write_combos_tsv(path: Path, combos: List[Tuple[str, str, int, int]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("ir\tdr\tsr\tlb\n")
        for ir, dr, sr, lb in combos:
            out.write(f"{ir}\t{dr}\t{sr}\t{lb}\n")


def write_slurm_array_script(
    sbatch_path: Path,
    args,
    combos_tsv: Path,
    n_tasks: int,
):
    """
    Writes an sbatch script that runs this same python file in --run-one mode.
    """
    sbatch_path.parent.mkdir(parents=True, exist_ok=True)

    job_name = args.slurm_job_name or "printe_grid"
    array_spec = f"0-{n_tasks-1}"
    script_abs = Path(__file__).resolve()

    # Slurm uses 1-based arrays sometimes in user configs, but 0-based is fine.
    # We'll map SLURM_ARRAY_TASK_ID -> line index + 2 (header line + 1).
    contents = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --array={array_spec}
#SBATCH --cpus-per-task={args.slurm_cpus}
#SBATCH --mem={args.slurm_mem}
#SBATCH --time={args.slurm_time}
#SBATCH --output={args.slurm_outdir}/%x_%A_%a.out
#SBATCH --error={args.slurm_outdir}/%x_%A_%a.err
{f"#SBATCH --partition={args.slurm_partition}" if args.slurm_partition else ""}
{f"#SBATCH --account={args.slurm_account}" if args.slurm_account else ""}

set -euo pipefail

echo "[INFO] Host: $(hostname)"
echo "[INFO] Date: $(date)"
echo "[INFO] PWD:  $(pwd)"
echo "[INFO] Task: ${{SLURM_ARRAY_TASK_ID}}"

python "{script_abs}" \\
  --run-one \\
  --combo-file "{combos_tsv}" \\
  --combo-index "${{SLURM_ARRAY_TASK_ID}}" \\
  --printe-script "{args.printe_script}" \\
  --ge "{args.ge}" \\
  --st "{args.st}" \\
  --mut "{args.mut}" \\
  --mxgs "{args.mxgs}" \\
  --mngs "{args.mngs}" \\
  --tstv "{args.tstv}" \\
  --threads "{args.threads}" \\
  --bed "{args.bed}" \\
  --fasta "{args.fasta}" \\
  --te-lib "{args.te_lib}" \\
  --ratios "{args.ratios}"
"""

    with sbatch_path.open("w") as f:
        f.write(contents)


def read_combo_from_tsv(combo_file: Path, combo_index: int) -> Tuple[str, str, int, int]:
    """
    combo_index is 0-based index over data lines (excluding header).
    """
    with combo_file.open("r") as f:
        header = f.readline()
        if not header:
            raise ValueError(f"Empty combo file: {combo_file}")
        # Read until desired line
        for i, line in enumerate(f):
            if i == combo_index:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 4:
                    raise ValueError(f"Bad combo line at index {combo_index}: {line}")
                ir, dr, sr, lb = parts
                return ir, dr, int(sr), int(lb)
    raise IndexError(f"combo_index {combo_index} out of range for {combo_file}")


def run_one_combo_from_file(args):
    base_cwd = Path.cwd().resolve()
    combo = read_combo_from_tsv(Path(args.combo_file), args.combo_index)

    ir, dr, sr, lb = combo
    dname = dir_name_from_combo(ir, dr, sr, lb)
    workdir = Path(dname)

    if workdir.exists():
        print(f"[SKIP] {dname} already exists; presumed complete.")
        raise SystemExit(0)

    workdir.mkdir(exist_ok=True)

    cmd = build_printe_cmd(combo, args, workdir, base_cwd)

    print(f"[RUN-ONE] combo_index={args.combo_index} -> {dname}")
    print("  " + " ".join(cmd))

    # Run and propagate return code to Slurm
    completed = subprocess.run(cmd, cwd=str(workdir))
    raise SystemExit(completed.returncode)


def parse_args():
    p = argparse.ArgumentParser(
        description="Grid-search wrapper for PrinTE.sh with random sampling. "
                    "Supports local parallelism OR Slurm job arrays."
    )

    # Modes
    p.add_argument("--mode", choices=["local", "slurm"], default="local",
                   help="Run mode: local (default) or slurm (generate array jobs).")
    p.add_argument("--run-one", action="store_true",
                   help=argparse.SUPPRESS)  # internal for array tasks
    p.add_argument("--combo-file", default=None,
                   help=argparse.SUPPRESS)
    p.add_argument("--combo-index", type=int, default=None,
                   help=argparse.SUPPRESS)

    # Sampling / parallelism (local)
    p.add_argument("--samples", "-n", type=int, default=200,
                   help="Number of random parameter combinations to run (default: 200).")
    p.add_argument("--max-procs", "-p", type=int, default=2,
                   help="Maximum number of parallel PrinTE.sh runs (default: 2).")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed (optional, for reproducibility).")

    # Parameter ranges – insertion_rates (log-spaced)
    p.add_argument("--ins-start", type=float, default=1e-12)
    p.add_argument("--ins-end", type=float, default=5e-12)
    p.add_argument("--ins-count", type=int, default=25)

    # Parameter ranges – deletion_rates (log-spaced)
    p.add_argument("--del-start", type=float, default=1e-11)
    p.add_argument("--del-end", type=float, default=6e-11)
    p.add_argument("--del-count", type=int, default=30)

    # Parameter ranges – solo_ratio
    p.add_argument("--sr-start", type=int, default=5)
    p.add_argument("--sr-end", type=int, default=95)
    p.add_argument("--sr-step", type=int, default=5)

    # Parameter ranges – length_bias
    p.add_argument("--k-start", type=int, default=0)
    p.add_argument("--k-end", type=int, default=20)
    p.add_argument("--k-step", type=int, default=2)

    # PrinTE settings / inputs
    p.add_argument("--printe-script", required=True)
    p.add_argument("--ge", type=int, required=True)
    p.add_argument("--st", type=int, required=True)
    p.add_argument("--mut", required=True)
    p.add_argument("--mxgs", required=True)
    p.add_argument("--mngs", required=True)
    p.add_argument("--tstv", type=float, required=True)
    p.add_argument("--threads", "-t", type=int, required=True)
    p.add_argument("--bed", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--te-lib", required=True)
    p.add_argument("--ratios", required=True)

    # Slurm options (only used in --mode slurm)
    p.add_argument("--slurm-dir", default="slurm_grid",
                   help="Directory to write combos.tsv and sbatch script (default: slurm_grid).")
    p.add_argument("--slurm-submit", action="store_true",
                   help="If set, automatically runs sbatch on the generated script.")
    p.add_argument("--slurm-job-name", default=None,
                   help="Slurm job name (default: printe_grid).")
    p.add_argument("--slurm-partition", default=None,
                   help="Slurm partition/queue (optional).")
    p.add_argument("--slurm-account", default=None,
                   help="Slurm account (optional).")
    p.add_argument("--slurm-cpus", type=int, default=4,
                   help="cpus-per-task (default: 4). NOTE: this should usually match --threads for PrinTE "
                        "only if PrinTE is truly thread-scaled on-node; otherwise keep smaller.")
    p.add_argument("--slurm-mem", default="8G",
                   help="Memory per task, e.g. 8G, 32G (default: 8G).")
    p.add_argument("--slurm-time", default="24:00:00",
                   help="Walltime, e.g. 02:00:00, 1-00:00:00 (default: 24:00:00).")
    p.add_argument("--slurm-outdir", default="slurm_logs",
                   help="Directory for Slurm stdout/err logs (default: slurm_logs).")

    return p.parse_args()


def main():
    args = parse_args()

    # Internal array task runner
    if args.run_one:
        if args.combo_file is None or args.combo_index is None:
            raise SystemExit("--run-one requires --combo-file and --combo-index")
        validate_inputs(args)
        run_one_combo_from_file(args)
        return

    if args.seed is not None:
        random.seed(args.seed)

    validate_inputs(args)

    available_combos, total, skipped = prepare_combos(args)

    print(f"Total possible combinations (full grid): {total}")
    print(f"Existing parameter dirs (skipped): {skipped}")
    print(f"Available new combinations: {len(available_combos)}")

    if not available_combos:
        print("No new combinations to run: all parameter directories already exist.")
        return

    n_samples = min(args.samples, len(available_combos))
    if n_samples < args.samples:
        print(f"Requested {args.samples} samples, but only {len(available_combos)} "
              f"new combos available. Using {n_samples}.")

    sampled = random.sample(available_combos, n_samples)
    print(f"Randomly selected {n_samples} NEW combinations.\n")

    if args.mode == "local":
        print(f"[MODE] local with up to {args.max_procs} parallel jobs.\n")
        base_cwd = Path.cwd().resolve()
        active = []
        for combo in sampled:
            wait_for_slot(active, args.max_procs)
            proc, name = launch_local_job(combo, args, base_cwd)
            active.append((proc, name))
        wait_for_all(active)
        print("All grid-search jobs completed (local).")
        return

    # Slurm mode
    slurm_dir = Path(args.slurm_dir)
    slurm_logs = Path(args.slurm_outdir)
    slurm_logs.mkdir(parents=True, exist_ok=True)

    combos_tsv = slurm_dir / "combos.tsv"
    sbatch_path = slurm_dir / "submit_array.sbatch"

    write_combos_tsv(combos_tsv, sampled)
    write_slurm_array_script(sbatch_path, args, combos_tsv, n_tasks=len(sampled))

    print(f"[MODE] slurm")
    print(f"Wrote combos: {combos_tsv}")
    print(f"Wrote sbatch: {sbatch_path}")
    print(f"Array size:   {len(sampled)} tasks (0..{len(sampled)-1})")
    print()
    print("To submit manually:")
    print(f"  sbatch {sbatch_path}")
    print()

    if args.slurm_submit:
        print("[SUBMIT] Running sbatch...")
        res = subprocess.run(["sbatch", str(sbatch_path)], check=False, capture_output=True, text=True)
        print(res.stdout.strip())
        if res.returncode != 0:
            print(res.stderr.strip())
            raise SystemExit(res.returncode)


if __name__ == "__main__":
    main()
