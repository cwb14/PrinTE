#!/usr/bin/env python3
import random
import subprocess
import itertools
import argparse
import math
from pathlib import Path


def dir_name_from_combo(ir, dr, sr, lb):
    return f"insertion_rates_{ir}_deletion_rates_{dr}_solo_ratio_{sr}_length_bias_{lb}"


def format_sci(x: float) -> str:
    """
    Format floats like 1e-07 -> '1e-7', 2.5e-12 -> '2.5e-12'
    so directory names stay readable/stable.
    """
    s = f"{x:.12g}"
    # Python sometimes emits 'E'; normalize to 'e'
    s = s.replace("E", "e")
    return s


def build_logspace_values(start: float, end: float, count: int):
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
    Use CLI arguments to build parameter grids.
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


def launch_job(combo, settings):
    ir, dr, sr, lb = combo
    dir_name = dir_name_from_combo(ir, dr, sr, lb)
    workdir = Path(dir_name)
    workdir.mkdir(exist_ok=True)

    cmd = [
        "bash", settings["PRINTE_SCRIPT"],
        "-ge", str(settings["GE"]),
        "-st", str(settings["ST"]),
        "-F", f"{ir},{dr}",
        "-m", str(settings["MUT"]),
        "-mxgs", str(settings["MXGS"]),
        "-mngs", str(settings["MNGS"]),
        "-sr", str(sr),
        "-k", str(lb),
        "-TsTv", str(settings["TSTV"]),
        "--ex_LTR",
        "--TE_lib", str(Path("..") / settings["TE_LIB"]),
        "-t", str(settings["THREADS"]),
        "--bed", str(Path("..") / settings["BED"]),
        "--fasta", str(Path("..") / settings["FASTA"]),
        "-r", str(settings["RATIOS"]),
    ]

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


def main():
    parser = argparse.ArgumentParser(
        description="Grid-search wrapper for PrinTE.sh with random sampling and parallel jobs."
    )

    # Sampling / parallelism
    parser.add_argument("--samples", "-n", type=int, default=200,
                        help="Number of random parameter combinations to run (default: 200).")
    parser.add_argument("--max-procs", "-p", type=int, default=2,
                        help="Maximum number of parallel PrinTE.sh runs (default: 2).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed (optional, for reproducibility).")

    # Parameter ranges – insertion_rates (direct numeric range)
    parser.add_argument("--ins-start", type=float, default=1e-12,
                        help="Insertion-rate range start (e.g., 1e-7). Default: 1e-12")
    parser.add_argument("--ins-end", type=float, default=5e-12,
                        help="Insertion-rate range end (e.g., 1e-14). Default: 5e-12")
    parser.add_argument("--ins-count", type=int, default=25,
                        help="Number of insertion-rate values to generate (log-spaced, inclusive). Default: 25")

    # Parameter ranges – deletion_rates (direct numeric range)
    parser.add_argument("--del-start", type=float, default=1e-11,
                        help="Deletion-rate range start (e.g., 2e-8). Default: 1e-11")
    parser.add_argument("--del-end", type=float, default=6e-11,
                        help="Deletion-rate range end (e.g., 6e-15). Default: 6e-11")
    parser.add_argument("--del-count", type=int, default=30,
                        help="Number of deletion-rate values to generate (log-spaced, inclusive). Default: 30")

    # Parameter ranges – solo_ratio
    parser.add_argument("--sr-start", type=int, default=5,
                        help="Solo ratio start (default: 5)")
    parser.add_argument("--sr-end", type=int, default=95,
                        help="Solo ratio end (default: 95)")
    parser.add_argument("--sr-step", type=int, default=5,
                        help="Solo ratio step (default: 5)")

    # Parameter ranges – length_bias
    parser.add_argument("--k-start", type=int, default=0,
                        help="Length bias start (default: 0)")
    parser.add_argument("--k-end", type=int, default=20,
                        help="Length bias end (default: 20)")
    parser.add_argument("--k-step", type=int, default=2,
                        help="Length bias step (default: 2)")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    # Fixed arguments (edit here if needed)
    settings = {
        "PRINTE_SCRIPT": "../PrinTE/PrinTE.sh",
        "GE": 5100000,
        "ST": 1020000,
        "MUT": "2e-9",
        "MXGS": "175M",
        "MNGS": "80M",
        "TSTV": 2.0,
        "THREADS": 100,
        "BED": "./burnin.bed",
        "FASTA": "./burnin.fasta",
        "TE_LIB": "./maize_rice_arab_curated_TE.lib.fa",
        "RATIOS": "../PrinTE/ratios_ltr_only.tsv",
    }

    insertion_rates, deletion_rates, solo_ratios, length_biases = build_param_grids(args)
    combos = build_all_combinations(insertion_rates, deletion_rates, solo_ratios, length_biases)

    total = len(combos)
    print(f"Total possible combinations (full grid): {total}")

    # Filter out combos whose directories already exist
    available_combos = []
    skipped = 0
    for ir, dr, sr, lb in combos:
        dname = dir_name_from_combo(ir, dr, sr, lb)
        if Path(dname).exists():
            skipped += 1
        else:
            available_combos.append((ir, dr, sr, lb))

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
    print(f"Randomly selected {n_samples} NEW combinations.")
    print(f"Running with up to {args.max_procs} parallel jobs.\n")

    active = []
    for combo in sampled:
        wait_for_slot(active, args.max_procs)
        proc, name = launch_job(combo, settings)
        active.append((proc, name))

    wait_for_all(active)
    print("All grid-search jobs completed.")


if __name__ == "__main__":
    main()
