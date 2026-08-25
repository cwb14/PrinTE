"""Shared utilities for PrinTE grid and guided search."""

import math
import os
import stat
from pathlib import Path


def format_sci(x: float) -> str:
    """
    Format floats like 1e-07 -> '1e-7', 2.5e-12 -> '2.5e-12'
    so directory names stay readable/stable.
    """
    s = f"{x:.12g}"
    s = s.replace("E", "e")
    # Strip leading zeros from exponent: e-07 -> e-7, e+03 -> e+3
    if "e-" in s:
        base, exp = s.split("e-")
        s = f"{base}e-{exp.lstrip('0') or '0'}"
    elif "e+" in s:
        base, exp = s.split("e+")
        s = f"{base}e+{exp.lstrip('0') or '0'}"
    elif "e" in s:
        base, exp = s.split("e")
        s = f"{base}e{exp.lstrip('0') or '0'}"
    return s


def build_logspace_values(start: float, end: float, count: int) -> list[str]:
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


def dir_name_from_combo(ir, dr, sr, lb):
    return f"insertion_rates_{ir}_deletion_rates_{dr}_solo_ratio_{sr}_length_bias_{lb}"


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


def validate_inputs(args):
    """
    Light validation: ensure required script/files exist (as typed),
    so tab-completion helps users catch mistakes early.
    """
    te_lib = getattr(args, "te_lib", None)
    clean_lib = getattr(args, "clean_lib", None)
    lib_label = "--clean-lib" if clean_lib else "--te-lib"
    lib_path = clean_lib or te_lib
    if not lib_path:
        raise FileNotFoundError("Provide either --te-lib or --clean-lib")

    for label, p, must_exist in [
        ("--printe-script", args.printe_script, True),
        ("--bed", args.bed, True),
        ("--fasta", args.fasta, True),
        (lib_label, lib_path, True),
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


def build_printe_cmd(
    combo: tuple[str, str, int, int],
    args,
    workdir: Path,
    base_cwd: Path
) -> list[str]:
    ir, dr, sr, lb = combo

    printe_script = relpath_from(workdir, args.printe_script, base_cwd)
    bed = relpath_from(workdir, args.bed, base_cwd)
    fasta = relpath_from(workdir, args.fasta, base_cwd)
    ratios = relpath_from(workdir, args.ratios, base_cwd)

    clean_lib = getattr(args, "clean_lib", None)
    if clean_lib:
        lib_flag, lib_path = "--clean_lib", relpath_from(workdir, clean_lib, base_cwd)
    else:
        lib_flag, lib_path = "--TE_lib", relpath_from(workdir, args.te_lib, base_cwd)

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
        lib_flag, str(lib_path),
        "-t", str(args.threads),
        "--bed", str(bed),
        "--fasta", str(fasta),
        "-r", str(ratios),
    ]
    return cmd


def write_combos_tsv(path: Path, combos: list[tuple[str, str, int, int]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("ir\tdr\tsr\tlb\n")
        for ir, dr, sr, lb in combos:
            out.write(f"{ir}\t{dr}\t{sr}\t{lb}\n")


def read_combo_from_tsv(combo_file: Path, combo_index: int) -> tuple[str, str, int, int]:
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


def write_slurm_array_script(
    sbatch_path: Path,
    args,
    combos_tsv: Path,
    n_tasks: int,
    runner_script: Path,
):
    """Probe the cluster, plan an efficient/portable sbatch array, and write it.

    Thin orchestrator over slurm_introspect (discover) + slurm_sbatch (plan +
    render). `n_tasks` is the number of combos; the emitted array size may be
    smaller when combos are packed onto exclusive whole nodes. Signature kept
    stable for guided_search.py and the Nextflow sweep.
    """
    from .slurm_introspect import probe_cluster
    from .slurm_sbatch import (
        build_request_from_args,
        plan_sbatch,
        render_chunk_submitter,
        render_sbatch_script,
        summarize_plan,
    )

    sbatch_path = Path(sbatch_path)
    sbatch_path.parent.mkdir(parents=True, exist_ok=True)

    req = build_request_from_args(args, n_combos=int(n_tasks))
    profile = probe_cluster(req.partition, timeout_s=req.probe_timeout,
                            do_probe=not req.no_probe)
    plan = plan_sbatch(profile, req)

    script = render_sbatch_script(plan, args, Path(runner_script), Path(combos_tsv))
    sbatch_path.write_text(script)

    submitter = render_chunk_submitter(plan, sbatch_path)
    if submitter:
        sub_path = sbatch_path.parent / "submit_chunks.sh"
        sub_path.write_text(submitter)
        sub_path.chmod(sub_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    print(summarize_plan(plan))
    if submitter:
        print(f"[SLURM] Array exceeds MaxArraySize -> submit with: bash "
              f"{sbatch_path.parent / 'submit_chunks.sh'}")
    return sbatch_path


def write_local_array_script(
    sh_path: Path,
    args,
    combos_tsv: Path,
    n_tasks: int,
    runner_script: Path,
    threads_per_sample: int,
):
    """
    Write a non-SLURM launcher: a bash script that backgrounds all combos via
    `nohup ... &`, splitting `args.threads` across them so each task gets
    `threads_per_sample` threads.

    Per-task stdout/stderr go to <args.slurm_outdir>/task_<i>.{out,err}.
    The script `wait`s on all backgrounded children at the end so a single
    wrapper PID indicates overall completion.

    Recommended invocation (so processes survive logout):
        nohup bash <sh_path> > submit.log 2>&1 &
    """
    sh_path.parent.mkdir(parents=True, exist_ok=True)

    last_idx = n_tasks - 1
    logdir = args.slurm_outdir or "slurm_logs"
    total_threads = args.threads

    clean_lib = getattr(args, "clean_lib", None)
    lib_arg = f'--clean-lib "{clean_lib}"' if clean_lib else f'--te-lib "{args.te_lib}"'

    contents = f"""#!/bin/bash
# Non-SLURM launcher: runs all {n_tasks} combos in parallel.
# Each combo gets {threads_per_sample} thread(s)  ({total_threads} total / {n_tasks} tasks).
#
# Recommended (survives logout):
#   nohup bash "$0" > submit.log 2>&1 &

set -u

LOGDIR="{logdir}"
mkdir -p "$LOGDIR"

echo "[INFO] Host: $(hostname)"
echo "[INFO] Date: $(date)"
echo "[INFO] PWD:  $(pwd)"
echo "[INFO] Launching {n_tasks} task(s), {threads_per_sample} thread(s) each"

PIDS=()
for i in $(seq 0 {last_idx}); do
  nohup python "{runner_script}" \\
    --run-one \\
    --combo-file "{combos_tsv}" \\
    --combo-index "$i" \\
    --printe-script "{args.printe_script}" \\
    --ge "{args.ge}" \\
    --st "{args.st}" \\
    --mut "{args.mut}" \\
    --mxgs "{args.mxgs}" \\
    --mngs "{args.mngs}" \\
    --tstv "{args.tstv}" \\
    --threads "{threads_per_sample}" \\
    --bed "{args.bed}" \\
    --fasta "{args.fasta}" \\
    {lib_arg} \\
    --ratios "{args.ratios}" \\
    > "$LOGDIR/task_${{i}}.out" 2> "$LOGDIR/task_${{i}}.err" &
  PIDS+=($!)
done

echo "[INFO] Launched ${{#PIDS[@]}} background task(s). PIDs: ${{PIDS[*]}}"
echo "[INFO] Per-task logs: $LOGDIR/task_<i>.{{out,err}}"
echo "[INFO] Waiting for all tasks to finish..."
wait
echo "[INFO] All tasks complete: $(date)"
"""

    with sh_path.open("w") as f:
        f.write(contents)
    sh_path.chmod(sh_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
