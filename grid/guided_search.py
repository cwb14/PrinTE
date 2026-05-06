#!/usr/bin/env python3
"""
Surrogate-model-guided adaptive parameter search for PrinTE simulations.

Subcommands:
  init    -- Round 0 exploration using Latin Hypercube Sampling
  harvest -- Parse simulation results into training data
  next    -- Train surrogate model and select guided candidates
"""

import argparse
import csv
import json
import math
import pickle
import random as stdlib_random
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats.qmc import LatinHypercube
from sklearn.ensemble import RandomForestRegressor

from grid_utils import (
    build_printe_cmd,
    dir_name_from_combo,
    format_sci,
    read_combo_from_tsv,
    validate_inputs,
    write_combos_tsv,
    write_local_array_script,
    write_slurm_array_script,
)


# ---------------------------------------------------------------------------
# Regex patterns for parsing directory names and projection log lines
# ---------------------------------------------------------------------------
_DIR_RE = re.compile(
    r"^insertion_rates_(?P<ins>[^_]+(?:e[+-]?\d+)?)"
    r"_deletion_rates_(?P<del>[^_]+(?:e[+-]?\d+)?)"
    r"_solo_ratio_(?P<sr>\d+)"
    r"_length_bias_(?P<k>\d+)$"
)

_PROJ_RE = re.compile(
    r"Projection\s+\(iter\s+\d+,\s*n=\d+,\s*\d+%\s*observed,\s*R²=(?P<r2>-?[0-9.]+)\):\s*"
    r"point=(?P<point>\d+),\s*empirical\s+95%\s+PI=\[(?P<lo>\d+)\.\.(?P<hi>\d+)\]"
)

_TERM_TOO_LARGE_RE = re.compile(r"95% prediction interval lower bound .* exceeds maximum")
_TERM_TOO_SMALL_RE = re.compile(r"95% prediction interval upper bound .* is below minimum")
_TERM_CONTINUING_RE = re.compile(r"Projection within bounds .* Continuing")

# Training data TSV column order
_TSV_COLUMNS = [
    "round", "ins", "del", "sr", "k", "genome_size",
    "source", "direction", "r2", "pi_lo", "pi_hi", "weight",
]


def lhs_sample(
    n: int,
    ins_start: float, ins_end: float,
    del_start: float, del_end: float,
    sr_start: int, sr_end: int, sr_step: int,
    k_start: int, k_end: int, k_step: int,
    seed: int = None,
) -> List[Tuple[str, str, int, int]]:
    """
    Latin Hypercube Sample over the 4D parameter space.
    ins/del are sampled in log10 space (continuous), then formatted via format_sci.
    sr and k are sampled continuously, then snapped to the nearest valid integer step.
    """
    sr_values = list(range(sr_start, sr_end + 1, sr_step))
    k_values = list(range(k_start, k_end + 1, k_step))

    sampler = LatinHypercube(d=4, seed=seed)
    unit_samples = sampler.random(n=n)

    log_ins_lo = math.log10(min(ins_start, ins_end))
    log_ins_hi = math.log10(max(ins_start, ins_end))
    log_del_lo = math.log10(min(del_start, del_end))
    log_del_hi = math.log10(max(del_start, del_end))

    combos = []
    for row in unit_samples:
        log_ins = log_ins_lo + row[0] * (log_ins_hi - log_ins_lo)
        ins_val = 10 ** log_ins

        log_del = log_del_lo + row[1] * (log_del_hi - log_del_lo)
        del_val = 10 ** log_del

        sr_idx = int(round(row[2] * (len(sr_values) - 1)))
        sr_val = sr_values[sr_idx]

        k_idx = int(round(row[3] * (len(k_values) - 1)))
        k_val = k_values[k_idx]

        combos.append((format_sci(ins_val), format_sci(del_val), sr_val, k_val))

    return combos


# ---------------------------------------------------------------------------
# Harvest: parse simulation results into training data
# ---------------------------------------------------------------------------

def parse_dir_name(name: str) -> Optional[Tuple[str, str, int, int]]:
    """
    Extract (ins, del, sr, k) from a parameter directory name.

    Expected format:
        insertion_rates_{ir}_deletion_rates_{dr}_solo_ratio_{sr}_length_bias_{k}

    Returns None if the name does not match.
    """
    m = _DIR_RE.match(name)
    if not m:
        return None
    return (m.group("ins"), m.group("del"), int(m.group("sr")), int(m.group("k")))


def parse_last_projection(log_path: Path) -> Optional[Dict]:
    """
    Parse the last Projection line and its termination direction from a
    pipeline.log file.

    Returns a dict with keys: point, r2, pi_lo, pi_hi, direction
    or None if the log does not exist or contains no projection lines.

    Direction is one of: "too_large", "too_small", "continuing".
    """
    log_path = Path(log_path)
    if not log_path.exists():
        return None

    text = log_path.read_text()
    lines = text.splitlines()

    # Find the last projection line and its index
    last_proj = None
    last_proj_idx = -1
    for i, line in enumerate(lines):
        m = _PROJ_RE.search(line)
        if m:
            last_proj = m
            last_proj_idx = i

    if last_proj is None:
        return None

    # Determine direction from lines after the last projection
    direction = "continuing"  # default if no termination line follows
    for line in lines[last_proj_idx + 1:]:
        if _TERM_TOO_LARGE_RE.search(line):
            direction = "too_large"
            break
        elif _TERM_TOO_SMALL_RE.search(line):
            direction = "too_small"
            break
        elif _TERM_CONTINUING_RE.search(line):
            direction = "continuing"
            # Don't break: a later projection might override this

    return {
        "point": int(last_proj.group("point")),
        "r2": float(last_proj.group("r2")),
        "pi_lo": int(last_proj.group("lo")),
        "pi_hi": int(last_proj.group("hi")),
        "direction": direction,
    }


def harvest_results(
    run_dir: Path,
    ge: int,
    mngs: int,
    mxgs: int,
    current_round: int,
) -> List[Dict]:
    """
    Scan simulation directories under run_dir, build training rows.

    Parameters
    ----------
    run_dir : Path
        Parent directory containing parameter subdirectories.
    ge : int
        Total generations (used to find the terminal FASTA file).
    mngs, mxgs : int
        Min/max genome size bounds for the target range.
    current_round : int
        Round number to tag in training data.

    Returns a list of row dicts (one per parseable simulation).
    Skips directories with no parseable output.
    """
    run_dir = Path(run_dir)
    terminal_name = f"gen{ge}_final.fasta"
    rows = []
    skipped = 0

    for entry in sorted(run_dir.iterdir()):
        if not entry.is_dir():
            continue

        params = parse_dir_name(entry.name)
        if params is None:
            continue

        ins, del_, sr, k = params
        terminal_file = entry / terminal_name

        if terminal_file.exists():
            # Completed simulation: use actual file size
            genome_size = terminal_file.stat().st_size
            rows.append({
                "round": current_round,
                "ins": ins,
                "del": del_,
                "sr": sr,
                "k": k,
                "genome_size": genome_size,
                "source": "completed",
                "direction": "hit",
                "r2": "",
                "pi_lo": "",
                "pi_hi": "",
                "weight": 1.0,
            })
        else:
            # Not completed: try to parse projection from pipeline.log
            log_path = entry / "pipeline.log"
            proj = parse_last_projection(log_path)
            if proj is None:
                skipped += 1
                continue

            rows.append({
                "round": current_round,
                "ins": ins,
                "del": del_,
                "sr": sr,
                "k": k,
                "genome_size": proj["point"],
                "source": "projected",
                "direction": proj["direction"],
                "r2": proj["r2"],
                "pi_lo": proj["pi_lo"],
                "pi_hi": proj["pi_hi"],
                "weight": proj["r2"],
            })

    if skipped:
        import sys
        print(f"harvest: skipped {skipped} dir(s) with no parseable output",
              file=sys.stderr)

    return rows


def write_training_data(
    path: Path,
    rows: List[Dict],
    append: bool = False,
) -> None:
    """
    Write (or append) training data rows to a TSV file.

    When appending, deduplicates by (ins, del, sr, k) key, keeping the
    newer row (from the incoming rows list) over existing data.
    """
    path = Path(path)

    if append and path.exists():
        existing = read_training_data(path)
    else:
        existing = []

    # Build a keyed dict: existing first, then new rows overwrite
    seen = {}
    for row in existing:
        key = (str(row["ins"]), str(row["del"]), str(row["sr"]), str(row["k"]))
        seen[key] = row
    for row in rows:
        key = (str(row["ins"]), str(row["del"]), str(row["sr"]), str(row["k"]))
        seen[key] = row

    merged = list(seen.values())

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in merged:
            writer.writerow({col: row.get(col, "") for col in _TSV_COLUMNS})


def read_training_data(path: Path) -> List[Dict]:
    """
    Read training data TSV back with proper types.

    Type conversions:
    - round, sr, k, genome_size -> int
    - weight -> float
    - r2, pi_lo, pi_hi -> float (if non-empty, else "")
    """
    path = Path(path)
    rows = []

    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for raw in reader:
            row = dict(raw)

            # Integer fields
            row["round"] = int(row["round"])
            row["sr"] = int(row["sr"])
            row["k"] = int(row["k"])
            row["genome_size"] = int(float(row["genome_size"]))

            # Float field (always present)
            row["weight"] = float(row["weight"])

            # Optional float fields (may be empty string)
            for col in ("r2", "pi_lo", "pi_hi"):
                val = row[col]
                if val and val.strip():
                    row[col] = float(val)
                else:
                    row[col] = ""

            rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# Surrogate model: train, score, and select candidates
# ---------------------------------------------------------------------------

def train_surrogate(
    rows: List[Dict],
) -> Tuple[RandomForestRegressor, float]:
    """
    Train a Random Forest surrogate on training data rows.

    Features: [log10(ins), log10(del), sr, k]
    Target:   log10(genome_size)

    Returns (fitted_model, oob_r2_score).
    Raises ValueError if fewer than 30 rows are provided.
    """
    n = len(rows)
    if n < 30:
        raise ValueError(
            f"Need at least 30 training observations, got {n}"
        )

    X = np.array([
        [math.log10(float(r["ins"])),
         math.log10(float(r["del"])),
         r["sr"],
         r["k"]]
        for r in rows
    ])
    y = np.array([math.log10(float(r["genome_size"])) for r in rows])
    weights = np.array([float(r["weight"]) for r in rows])

    model = RandomForestRegressor(
        n_estimators=200,
        oob_score=True,
        random_state=0,
        n_jobs=-1,
    )
    model.fit(X, y, sample_weight=weights)

    return model, model.oob_score_


def score_candidates(
    model: RandomForestRegressor,
    candidates: List[Tuple[str, str, int, int]],
    target_midpoint: float,
    alpha: float,
) -> List[Tuple[Tuple[str, str, int, int], float, float, float]]:
    """
    Score candidate parameter combos using an acquisition function.

    Score = -|predicted_log_size - log10(target_midpoint)| + alpha * uncertainty

    Uncertainty is the standard deviation of per-tree predictions.

    Returns list of (combo, score, predicted_log_gs, uncertainty).
    """
    log_target = math.log10(target_midpoint)

    X = np.array([
        [math.log10(float(c[0])),
         math.log10(float(c[1])),
         c[2],
         c[3]]
        for c in candidates
    ])

    # Collect per-tree predictions to compute mean and std
    tree_preds = np.array([tree.predict(X) for tree in model.estimators_])
    means = tree_preds.mean(axis=0)
    stds = tree_preds.std(axis=0)

    results = []
    for i, combo in enumerate(candidates):
        pred = float(means[i])
        unc = float(stds[i])
        score = -abs(pred - log_target) + alpha * unc
        results.append((combo, score, pred, unc))

    return results


def select_candidates(
    model: RandomForestRegressor,
    target_midpoint: float,
    alpha: float,
    n_select: int,
    explore_frac: float,
    candidate_pool_size: int,
    existing_keys: set,
    ins_start: float,
    ins_end: float,
    del_start: float,
    del_end: float,
    sr_start: int,
    sr_end: int,
    sr_step: int,
    k_start: int,
    k_end: int,
    k_step: int,
    seed: int = None,
) -> List[Tuple[str, str, int, int]]:
    """
    Generate a pool of random candidates, score with the surrogate model,
    and select the top-N via guided exploitation + random exploration.

    ins/del are sampled log-uniformly; sr/k are sampled uniformly from
    their valid step grids.

    Returns a list of (ins_str, del_str, sr, k) tuples.
    """
    rng = stdlib_random.Random(seed)

    sr_values = list(range(sr_start, sr_end + 1, sr_step))
    k_values = list(range(k_start, k_end + 1, k_step))

    log_ins_lo = math.log10(min(ins_start, ins_end))
    log_ins_hi = math.log10(max(ins_start, ins_end))
    log_del_lo = math.log10(min(del_start, del_end))
    log_del_hi = math.log10(max(del_start, del_end))

    # Generate random candidate pool, filtering out existing keys
    pool = []
    attempts = 0
    max_attempts = candidate_pool_size * 10  # safety valve
    while len(pool) < candidate_pool_size and attempts < max_attempts:
        attempts += 1
        log_ins = rng.uniform(log_ins_lo, log_ins_hi)
        ins_str = format_sci(10 ** log_ins)

        log_del = rng.uniform(log_del_lo, log_del_hi)
        del_str = format_sci(10 ** log_del)

        sr = rng.choice(sr_values)
        k = rng.choice(k_values)

        key = (ins_str, del_str, sr, k)
        if key not in existing_keys:
            pool.append(key)

    # Score all candidates
    scored = score_candidates(model, pool, target_midpoint, alpha)

    # Sort by score descending
    scored.sort(key=lambda x: x[1], reverse=True)

    # Split into guided (exploitation) and exploration slots
    n_explore = max(1, int(n_select * explore_frac))
    n_guided = n_select - n_explore

    # Take top-N guided candidates
    guided = [item[0] for item in scored[:n_guided]]

    # Fill exploration slots with random picks from the remaining pool
    remaining = [item[0] for item in scored[n_guided:]]
    if remaining:
        explore = rng.sample(remaining, min(n_explore, len(remaining)))
    else:
        explore = []

    return guided + explore


# ---------------------------------------------------------------------------
# Search state management: persist round-over-round bookkeeping as JSON
# ---------------------------------------------------------------------------

def create_search_state(
    path,
    ge: int,
    mngs: int,
    mxgs: int,
    alpha: float,
    printe_args: Dict,
    slurm_args: Dict,
    param_ranges: Optional[Dict] = None,
) -> None:
    """
    Initialise a search_state.json file for a new guided search run.

    Parameters
    ----------
    path : str or Path
        Where to write the JSON file.
    ge : int
        Total generations for the simulation.
    mngs, mxgs : int
        Min / max genome size bounds defining the target range.
    alpha : float
        Initial exploration weight for the acquisition function.
    printe_args, slurm_args : dict
        Passthrough arguments for PrinTE and SLURM job submission.
    param_ranges : dict, optional
        Extra parameter range info to merge into the state dict.
    """
    state = {
        "round": 0,
        "ge": ge,
        "mngs": mngs,
        "mxgs": mxgs,
        "target_midpoint": (mngs + mxgs) // 2,
        "alpha": alpha,
        "total_sims_launched": 0,
        "total_completions": 0,
        "printe_args": printe_args,
        "slurm_args": slurm_args,
        "rounds": [],
    }
    if param_ranges is not None:
        state.update(param_ranges)

    path = Path(path)
    path.write_text(json.dumps(state, indent=2))


def load_search_state(path) -> Dict:
    """Read and return the search state dict from a JSON file."""
    path = Path(path)
    return json.loads(path.read_text())


def update_search_state(
    path,
    round_num: int,
    samples: int,
    round_type: str,
    explore_frac: Optional[float],
    completions: int,
    alpha_decay: float = 0.7,
) -> None:
    """
    Record one round's results and advance the search state.

    Appends round info, increments counters, and decays alpha so that
    later rounds shift from exploration toward exploitation.

    Parameters
    ----------
    path : str or Path
        Path to search_state.json.
    round_num : int
        Zero-based round number just completed.
    samples : int
        Number of simulations launched this round.
    round_type : str
        Label for the round (e.g. "exploration", "guided").
    explore_frac : float or None
        Fraction of candidates chosen randomly (None for pure exploration rounds).
    completions : int
        Number of simulations that completed this round.
    alpha_decay : float
        Multiplicative decay factor applied to alpha each round (default 0.7).
    """
    state = load_search_state(path)

    round_info: Dict = {"round": round_num, "samples": samples, "type": round_type}
    if explore_frac is not None:
        round_info["explore_frac"] = explore_frac

    state["rounds"].append(round_info)
    state["round"] = round_num + 1
    state["total_sims_launched"] += samples
    state["total_completions"] += completions
    state["alpha"] *= alpha_decay

    path = Path(path)
    path.write_text(json.dumps(state, indent=2))


# ---------------------------------------------------------------------------
# CLI: argument parsing
# ---------------------------------------------------------------------------

def add_printe_args(parser):
    """Add PrinTE pass-through arguments."""
    parser.add_argument("--printe-script", required=True)
    parser.add_argument("--ge", type=int, required=True)
    parser.add_argument("--st", type=int, required=True)
    parser.add_argument("--mut", required=True)
    parser.add_argument("--mxgs", type=int, required=True)
    parser.add_argument("--mngs", type=int, required=True)
    parser.add_argument("--tstv", type=float, required=True)
    parser.add_argument("--threads", "-t", type=int, required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--fasta", required=True)
    lib_grp = parser.add_mutually_exclusive_group(required=True)
    lib_grp.add_argument("--te-lib")
    lib_grp.add_argument("--clean-lib",
                         help="Pre-built cleaned TE library; skips library processing in PrinTE.sh")
    parser.add_argument("--ratios", required=True)


def add_slurm_args(parser):
    """Add scheduler / SLURM job submission arguments."""
    parser.add_argument(
        "--scheduler",
        choices=["local", "slurm"],
        default="local",
        help="local: emit submit_array.sh that launches all combos in parallel via "
             "nohup, splitting --threads across them. slurm: emit submit_array.sbatch "
             "(an array job; each task gets the full --threads value). Default: local.",
    )
    parser.add_argument("--slurm-dir", default="slurm_grid")
    parser.add_argument("--slurm-submit", action="store_true",
        help="Auto-launch after writing the script: 'sbatch <path>' for slurm; "
             "'nohup bash <path> > submit.log 2>&1 &' for local.")
    parser.add_argument("--slurm-job-name", default="printe_grid")
    parser.add_argument("--slurm-partition", default=None)
    parser.add_argument("--slurm-account", default=None)
    parser.add_argument("--slurm-cpus", type=int, default=4)
    parser.add_argument("--slurm-mem", default="8G")
    parser.add_argument("--slurm-time", default="24:00:00")
    parser.add_argument("--slurm-outdir", default="slurm_logs")


def add_param_range_args(parser):
    """Add parameter range arguments for the 4D search space."""
    parser.add_argument("--ins-start", type=float, default=1e-12)
    parser.add_argument("--ins-end", type=float, default=5e-12)
    parser.add_argument("--del-start", type=float, default=1e-11)
    parser.add_argument("--del-end", type=float, default=6e-11)
    parser.add_argument("--sr-start", type=int, default=5)
    parser.add_argument("--sr-end", type=int, default=95)
    parser.add_argument("--sr-step", type=int, default=5)
    parser.add_argument("--k-start", type=int, default=0)
    parser.add_argument("--k-end", type=int, default=20)
    parser.add_argument("--k-step", type=int, default=2)


def parse_args():
    """Build the main argument parser with init/harvest/next subcommands."""
    p = argparse.ArgumentParser(
        description="Surrogate-model-guided adaptive parameter search for PrinTE."
    )
    p.add_argument("--verbose", "-v", action="store_true")

    sub = p.add_subparsers(dest="command")

    # init
    init_p = sub.add_parser("init", help="Round 0: LHS exploration")
    init_p.add_argument("--samples", "-n", type=int, default=200)
    init_p.add_argument("--seed", type=int, default=None)
    init_p.add_argument("--alpha", type=float, default=0.5)
    add_param_range_args(init_p)
    add_printe_args(init_p)
    add_slurm_args(init_p)

    # harvest
    harv_p = sub.add_parser("harvest", help="Parse simulation results")
    harv_p.add_argument("--run-dir", default=".")
    harv_p.add_argument("--ge", type=int, default=None)
    harv_p.add_argument("--mngs", type=int, default=None)
    harv_p.add_argument("--mxgs", type=int, default=None)

    # next
    next_p = sub.add_parser("next", help="Guided round using surrogate model")
    next_p.add_argument("--samples", "-n", type=int, default=150)
    next_p.add_argument("--explore-frac", type=float, default=0.15)
    next_p.add_argument("--alpha", type=float, default=None)
    next_p.add_argument("--alpha-decay", type=float, default=0.7)
    next_p.add_argument("--candidate-pool", type=int, default=10000)
    next_p.add_argument("--seed", type=int, default=None)

    return p.parse_args()


def parse_run_one_args():
    """Separate parser for --run-one mode (called from SLURM array tasks)."""
    p = argparse.ArgumentParser(description=argparse.SUPPRESS)
    p.add_argument("--run-one", action="store_true", required=True)
    p.add_argument("--combo-file", required=True)
    p.add_argument("--combo-index", type=int, required=True)
    add_printe_args(p)
    return p.parse_args()


# ---------------------------------------------------------------------------
# CLI: subcommand handlers
# ---------------------------------------------------------------------------

def _write_launch_script(scheduler, *, slurm_dir, args, combos_tsv, n_tasks, runner_script):
    """
    Dispatch to the scheduler-specific launcher writer.

    Returns (submit_path, per_task_threads). per_task_threads is an int for
    local mode and None for slurm mode. Prints budget warnings inline.
    """
    if scheduler == "slurm":
        sbatch_path = slurm_dir / "submit_array.sbatch"
        write_slurm_array_script(sbatch_path, args, combos_tsv,
                                 n_tasks=n_tasks, runner_script=runner_script)
        return sbatch_path, None

    if scheduler == "local":
        total = int(args.threads)
        per = max(1, total // n_tasks)
        if total < n_tasks:
            print(f"WARNING: --threads {total} < {n_tasks} tasks; "
                  f"each task gets 1 thread (CPUs oversubscribed by "
                  f"{n_tasks - total} slots).")
        elif total % n_tasks != 0:
            idle = total - per * n_tasks
            print(f"NOTE: {per} thread(s) * {n_tasks} tasks = {per * n_tasks}; "
                  f"{idle} thread(s) unused.")
        sh_path = slurm_dir / "submit_array.sh"
        write_local_array_script(sh_path, args, combos_tsv, n_tasks=n_tasks,
                                 runner_script=runner_script,
                                 threads_per_sample=per)
        return sh_path, per

    raise SystemExit(f"Unknown --scheduler value: {scheduler!r}")


def _auto_submit(scheduler, submit_path):
    """Auto-launch the just-written script. Called when --slurm-submit is set."""
    if scheduler == "slurm":
        print("[SUBMIT] Running sbatch...")
        res = subprocess.run(["sbatch", str(submit_path)],
                             check=False, capture_output=True, text=True)
        print(res.stdout.strip())
        if res.returncode != 0:
            print(res.stderr.strip())
            raise SystemExit(res.returncode)
    else:
        print("[SUBMIT] Launching local job in background via nohup...")
        with open("submit.log", "ab") as logf:
            subprocess.Popen(
                ["nohup", "bash", str(submit_path)],
                stdout=logf, stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        print("[SUBMIT] Detached. Tail submit.log and slurm_logs/task_*.{out,err} for progress.")


def cmd_init(args):
    """Round 0: LHS exploration."""
    validate_inputs(args)

    combos = lhs_sample(
        n=args.samples,
        ins_start=args.ins_start, ins_end=args.ins_end,
        del_start=args.del_start, del_end=args.del_end,
        sr_start=args.sr_start, sr_end=args.sr_end, sr_step=args.sr_step,
        k_start=args.k_start, k_end=args.k_end, k_step=args.k_step,
        seed=args.seed,
    )

    # Filter existing directories
    original = len(combos)
    combos = [(ir, dr, sr, k) for ir, dr, sr, k in combos
              if not Path(dir_name_from_combo(ir, dr, sr, k)).exists()]
    skipped = original - len(combos)

    print(f"LHS generated {original} samples, {skipped} already exist, {len(combos)} new.")

    if not combos:
        print("No new combinations to run.")
        return

    slurm_dir = Path(args.slurm_dir)
    Path(args.slurm_outdir).mkdir(parents=True, exist_ok=True)
    combos_tsv = slurm_dir / "combos.tsv"

    write_combos_tsv(combos_tsv, combos)
    submit_path, per_threads = _write_launch_script(
        scheduler=args.scheduler, slurm_dir=slurm_dir, args=args,
        combos_tsv=combos_tsv, n_tasks=len(combos),
        runner_script=Path(__file__).resolve(),
    )

    # Persist state
    printe_args = {k: getattr(args, k) for k in
        ["printe_script", "ge", "st", "mut", "mxgs", "mngs",
         "tstv", "threads", "bed", "fasta", "te_lib", "clean_lib", "ratios"]}
    slurm_args_dict = {k: getattr(args, k) for k in
        ["slurm_dir", "slurm_job_name", "slurm_partition", "slurm_account",
         "slurm_cpus", "slurm_mem", "slurm_time", "slurm_outdir"]}
    slurm_args_dict["scheduler"] = args.scheduler
    param_ranges = {
        "ins_start": args.ins_start, "ins_end": args.ins_end,
        "del_start": args.del_start, "del_end": args.del_end,
        "sr_start": args.sr_start, "sr_end": args.sr_end, "sr_step": args.sr_step,
        "k_start": args.k_start, "k_end": args.k_end, "k_step": args.k_step,
    }

    state_path = Path("search_state.json")
    create_search_state(state_path, ge=args.ge, mngs=args.mngs, mxgs=args.mxgs,
                        alpha=args.alpha, printe_args=printe_args,
                        slurm_args=slurm_args_dict, param_ranges=param_ranges)
    update_search_state(state_path, round_num=0, samples=len(combos),
                        round_type="exploration", explore_frac=None, completions=0,
                        alpha_decay=1.0)  # no decay on init

    print(f"Wrote combos:  {combos_tsv}")
    if args.scheduler == "slurm":
        print(f"Wrote sbatch:  {submit_path}")
    else:
        print(f"Wrote launcher: {submit_path}")
        print(f"  Each task: {per_threads} thread(s)  ({args.threads} total / {len(combos)} tasks)")
    print(f"Wrote state:   {state_path}")
    print(f"Array size:    {len(combos)} tasks")
    if args.scheduler == "slurm":
        print(f"\nTo submit:  sbatch {submit_path}")
    else:
        print(f"\nTo run:  nohup bash {submit_path} > submit.log 2>&1 &")

    if args.slurm_submit:
        _auto_submit(args.scheduler, submit_path)


def cmd_harvest(args):
    """Parse simulation results into training data."""
    state_path = Path("search_state.json")
    run_dir = Path(args.run_dir)

    ge = args.ge
    mngs = args.mngs
    mxgs = args.mxgs
    current_round = 0

    if state_path.exists():
        state = load_search_state(state_path)
        if ge is None: ge = state["ge"]
        if mngs is None: mngs = state["mngs"]
        if mxgs is None: mxgs = state["mxgs"]
        current_round = state["round"] - 1

    if ge is None or mngs is None or mxgs is None:
        raise SystemExit("Cannot determine --ge, --mngs, --mxgs. Provide them or run 'init' first.")

    print(f"Harvesting from {run_dir}")
    print(f"  Terminal gen: gen{ge}_final.fasta")
    print(f"  Target range: [{mngs:,} .. {mxgs:,}]")

    rows = harvest_results(run_dir, ge, mngs, mxgs, current_round)

    completed = sum(1 for r in rows if r["source"] == "completed")
    projected = sum(1 for r in rows if r["source"] == "projected")
    too_large = sum(1 for r in rows if r["direction"] == "too_large")
    too_small = sum(1 for r in rows if r["direction"] == "too_small")

    tsv_path = Path("training_data.tsv")
    write_training_data(tsv_path, rows, append=tsv_path.exists())
    all_rows = read_training_data(tsv_path)

    print(f"\nHarvested: {completed} completed, {projected} projected")
    print(f"  Direction: {too_large} too_large, {too_small} too_small, {completed} hit")
    print(f"  Total training data: {len(all_rows)} rows")
    print(f"  Written to: {tsv_path}")

    directions = {r["direction"] for r in all_rows}
    if len(directions) == 1 and "hit" not in directions:
        d = directions.pop()
        print(f"\nWARNING: All {len(all_rows)} observations are {d}. The viable region may be "
              f"outside your parameter bounds. Consider shifting --ins-start/--ins-end "
              f"or --del-start/--del-end.")


def cmd_next(args):
    """Train surrogate model and select guided candidates."""
    state_path = Path("search_state.json")
    tsv_path = Path("training_data.tsv")

    if not state_path.exists():
        raise SystemExit("search_state.json not found. Run 'init' first.")
    if not tsv_path.exists():
        raise SystemExit("training_data.tsv not found. Run 'harvest' first.")

    state = load_search_state(state_path)
    rows = read_training_data(tsv_path)

    alpha = args.alpha if args.alpha is not None else state["alpha"]
    target_midpoint = state["target_midpoint"]

    print(f"Training surrogate on {len(rows)} observations...")
    model, oob_r2 = train_surrogate(rows)
    importances = model.feature_importances_

    if oob_r2 < 0.3:
        print(f"\nWARNING: Model OOB R2 = {oob_r2:.3f} (low). Consider increasing --explore-frac.")

    existing_keys = {(r["ins"], r["del"], r["sr"], r["k"]) for r in rows}

    selected = select_candidates(
        model=model, target_midpoint=target_midpoint, alpha=alpha,
        n_select=args.samples, explore_frac=args.explore_frac,
        candidate_pool_size=args.candidate_pool,
        existing_keys=existing_keys,
        ins_start=state.get("ins_start", 1e-17), ins_end=state.get("ins_end", 1e-10),
        del_start=state.get("del_start", 1e-17), del_end=state.get("del_end", 1e-10),
        sr_start=state.get("sr_start", 5), sr_end=state.get("sr_end", 95),
        sr_step=state.get("sr_step", 5),
        k_start=state.get("k_start", 0), k_end=state.get("k_end", 10),
        k_step=state.get("k_step", 1),
        seed=args.seed,
    )

    if not selected:
        print("No novel candidates remaining.")
        return

    n_guided = args.samples - max(1, int(args.samples * args.explore_frac))
    n_explore = len(selected) - n_guided
    completed = sum(1 for r in rows if r["source"] == "completed")
    projected = sum(1 for r in rows if r["source"] == "projected")

    print(f"\nRound {state['round']} summary:")
    print(f"  Training data: {len(rows)} sims ({completed} completed, {projected} projected)")
    print(f"  Model R2 (OOB): {oob_r2:.3f}")
    print(f"  Top features: log10(ins) {importances[0]:.2f}, log10(del) {importances[1]:.2f}, "
          f"sr {importances[2]:.2f}, k {importances[3]:.2f}")
    print(f"  Alpha: {alpha:.3f}")
    print(f"  Candidates selected: {n_guided} guided + {n_explore} random = {len(selected)}")

    sa = state.get("slurm_args", {})
    pa = state.get("printe_args", {})
    scheduler = sa.get("scheduler", "slurm")  # state files predating --scheduler default to slurm

    slurm_dir = Path(sa.get("slurm_dir", "slurm_grid"))
    Path(sa.get("slurm_outdir", "slurm_logs")).mkdir(parents=True, exist_ok=True)
    combos_tsv = slurm_dir / "combos.tsv"

    write_combos_tsv(combos_tsv, selected)

    slurm_ns = argparse.Namespace(**{**pa, **sa})
    submit_path, per_threads = _write_launch_script(
        scheduler=scheduler, slurm_dir=slurm_dir, args=slurm_ns,
        combos_tsv=combos_tsv, n_tasks=len(selected),
        runner_script=Path(__file__).resolve(),
    )

    # Save model
    round_num = state["round"]
    model_path = Path(f"model_round_{round_num}.pkl")
    with model_path.open("wb") as f:
        pickle.dump(model, f)

    update_search_state(state_path, round_num=round_num, samples=len(selected),
                        round_type="guided", explore_frac=args.explore_frac,
                        completions=0, alpha_decay=args.alpha_decay)

    print(f"\n  Wrote combos: {combos_tsv}")
    if scheduler == "slurm":
        print(f"  Wrote sbatch: {submit_path}")
    else:
        print(f"  Wrote launcher: {submit_path}")
        print(f"  Each task: {per_threads} thread(s)  "
              f"({pa.get('threads')} total / {len(selected)} tasks)")
    print(f"  Saved model:  {model_path}")
    if scheduler == "slurm":
        print(f"\n  To submit:  sbatch {submit_path}")
    else:
        print(f"\n  To run:  nohup bash {submit_path} > submit.log 2>&1 &")


# ---------------------------------------------------------------------------
# CLI: run-one mode (for SLURM array tasks) and main entry point
# ---------------------------------------------------------------------------

def run_one_combo(args):
    """Internal: run a single combo from a SLURM array task."""
    validate_inputs(args)
    combo = read_combo_from_tsv(Path(args.combo_file), args.combo_index)
    ir, dr, sr, lb = combo
    dname = dir_name_from_combo(ir, dr, sr, lb)
    workdir = Path(dname)

    if workdir.exists():
        print(f"[SKIP] {dname} already exists; presumed complete.")
        raise SystemExit(0)

    workdir.mkdir(exist_ok=True)
    base_cwd = Path.cwd().resolve()
    cmd = build_printe_cmd(combo, args, workdir, base_cwd)

    print(f"[RUN-ONE] combo_index={args.combo_index} -> {dname}")
    print("  " + " ".join(cmd))

    completed = subprocess.run(cmd, cwd=str(workdir))
    raise SystemExit(completed.returncode)


def main():
    """Entry point: dispatch to subcommand or --run-one mode."""
    import sys
    if "--run-one" in sys.argv:
        args = parse_run_one_args()
        run_one_combo(args)
        return

    args = parse_args()

    if args.command == "init":
        cmd_init(args)
    elif args.command == "harvest":
        cmd_harvest(args)
    elif args.command == "next":
        cmd_next(args)
    else:
        raise SystemExit("Usage: guided_search.py {init|harvest|next}. Use --help for details.")


if __name__ == "__main__":
    main()
