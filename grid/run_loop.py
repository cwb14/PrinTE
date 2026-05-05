#!/usr/bin/env python3
"""
Unattended driver for the guided-search loop.

Each round runs:
  1. guided_search.py harvest --run-dir .
  2. (count hits in training_data.tsv; stop if >= --target-hits)
  3. purge_intermediates.py --execute     (skipped with --no-purge)
  4. guided_search.py next --samples N --explore-frac F
  5. bash slurm_grid/submit_array.sh      (FOREGROUND — blocks until done)

The script is safe to start while a previous round is still running: at
startup it polls pipeline.log files and waits for any in-progress sims to
complete before harvesting, so CPU is never oversubscribed.

Run it under nohup so it survives ssh logout:

    nohup python ../PrinTE/grid/run_loop.py \\
        --target-hits 200 \\
        --max-rounds 50 \\
        --samples 100 --explore-frac 0.1 \\
        > orchestrator.log 2>&1 &

Tail orchestrator.log to monitor. To stop early: kill the orchestrator PID;
the round currently in `bash submit_array.sh` will continue (its tasks were
started with nohup) but no new round will be launched.

Resumable: re-launching the script just picks up at the next harvest. State
lives in search_state.json, training_data.tsv, and the on-disk param dirs.
"""

import argparse
import csv
import os
import re
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

_HERE = Path(__file__).resolve().parent
GUIDED_SEARCH = _HERE / "guided_search.py"
PURGE = _HERE / "purge_intermediates.py"

_DIR_RE = re.compile(
    r"^insertion_rates_[^_]+(?:e[+-]?\d+)?"
    r"_deletion_rates_[^_]+(?:e[+-]?\d+)?"
    r"_solo_ratio_\d+_length_bias_\d+$"
)
_FINISHED_MARKER = "Pipeline completed at"


def log(msg: str) -> None:
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {msg}", flush=True)


def is_finished(log_path: Path) -> bool:
    if not log_path.exists():
        return False
    try:
        with open(log_path, "rb") as fh:
            try:
                fh.seek(-4096, 2)
            except OSError:
                fh.seek(0)
            return _FINISHED_MARKER in fh.read().decode(errors="ignore")
    except OSError:
        return False


def count_inflight(run_dir: Path) -> int:
    """Count parameter dirs whose pipeline.log lacks the completion marker."""
    n = 0
    for d in run_dir.iterdir():
        if not d.is_dir() or not _DIR_RE.match(d.name):
            continue
        log_path = d / "pipeline.log"
        if log_path.exists() and not is_finished(log_path):
            n += 1
    return n


def count_hits(tsv_path: Path) -> int:
    if not tsv_path.exists():
        return 0
    n = 0
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("direction") == "hit":
                n += 1
    return n


def wait_for_inflight(run_dir: Path, poll_seconds: int) -> None:
    n = count_inflight(run_dir)
    if n == 0:
        return
    log(f"Found {n} in-progress sim(s) at startup. Waiting before first harvest...")
    last_n = n
    same_count = 0
    while True:
        n = count_inflight(run_dir)
        if n == 0:
            log("All in-progress sims complete; proceeding.")
            return
        if n == last_n:
            same_count += 1
        else:
            same_count = 0
            last_n = n
        if same_count and same_count % 30 == 0:
            log(f"  WARNING: in-progress count stuck at {n} for "
                f"{same_count * poll_seconds}s — sims may be hung. "
                "Investigate or kill orchestrator + clean up.")
        else:
            log(f"  {n} sim(s) still in-progress; sleeping {poll_seconds}s")
        time.sleep(poll_seconds)


def run_cmd(label: str, cmd, cwd: Path, capture_to: Path = None) -> int:
    # start_new_session=True puts the child (and its descendants) in their
    # own session/process group, detached from any controlling tty. Combined
    # with the orchestrator itself running under nohup, this means an ssh
    # logout cannot deliver SIGHUP to any leg of the pipeline.
    log(f"[{label}] $ {' '.join(str(x) for x in cmd)}  (cwd={cwd})")
    if capture_to is not None:
        with open(capture_to, "ab") as fh:
            res = subprocess.run(cmd, cwd=str(cwd),
                                 stdout=fh, stderr=subprocess.STDOUT,
                                 start_new_session=True)
    else:
        res = subprocess.run(cmd, cwd=str(cwd), start_new_session=True)
    log(f"[{label}] exit={res.returncode}")
    return res.returncode


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--run-dir", default=".",
                    help="Working dir containing search_state.json (default: cwd)")
    ap.add_argument("--target-hits", type=int, default=200,
                    help="Stop when training_data.tsv has this many 'hit' rows (default: 200)")
    ap.add_argument("--max-rounds", type=int, default=100,
                    help="Hard cap on rounds (safety; default: 100)")
    ap.add_argument("--samples", type=int, default=100,
                    help="--samples passed to 'next' (default: 100)")
    ap.add_argument("--explore-frac", type=float, default=0.1,
                    help="--explore-frac passed to 'next' (default: 0.1)")
    ap.add_argument("--poll-seconds", type=int, default=60,
                    help="Polling interval while waiting for in-progress sims (default: 60)")
    ap.add_argument("--no-purge", action="store_true",
                    help="Skip the purge_intermediates step")
    ap.add_argument("--start-round", type=int, default=1,
                    help="Round number to label the first iteration (default: 1)")
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.is_dir():
        ap.error(f"--run-dir does not exist: {run_dir}")
    if not (run_dir / "search_state.json").exists():
        ap.error(f"search_state.json not found in {run_dir}; run "
                 f"`guided_search.py init` first.")
    submit_script = run_dir / "slurm_grid" / "submit_array.sh"
    tsv_path = run_dir / "training_data.tsv"

    # Quiet SIGHUP so backgrounded launches survive logout cleanly.
    signal.signal(signal.SIGHUP, signal.SIG_IGN)

    log(f"=== Orchestrator starting (PID {os.getpid()}) ===")
    log(f"  run_dir:       {run_dir}")
    log(f"  target_hits:   {args.target_hits}")
    log(f"  max_rounds:    {args.max_rounds}")
    log(f"  samples/round: {args.samples}")
    log(f"  explore_frac:  {args.explore_frac}")
    log(f"  purge:         {'OFF (--no-purge)' if args.no_purge else 'ON'}")

    wait_for_inflight(run_dir, args.poll_seconds)

    for r in range(args.start_round, args.start_round + args.max_rounds):
        log(f"=== Round {r} ===")

        if run_cmd("harvest",
                   [sys.executable, str(GUIDED_SEARCH),
                    "harvest", "--run-dir", "."],
                   cwd=run_dir) != 0:
            log("Harvest failed; aborting.")
            sys.exit(2)

        hits = count_hits(tsv_path)
        log(f"  Hits so far: {hits} / target {args.target_hits}")
        if hits >= args.target_hits:
            log("Target hit count reached. Done.")
            return

        if not args.no_purge:
            # Don't abort on purge failure — it's housekeeping, not pipeline.
            run_cmd("purge",
                    [sys.executable, str(PURGE),
                     "--run-dir", ".", "--execute"],
                    cwd=run_dir)

        if run_cmd("next",
                   [sys.executable, str(GUIDED_SEARCH), "next",
                    "--samples", str(args.samples),
                    "--explore-frac", str(args.explore_frac)],
                   cwd=run_dir) != 0:
            log("'next' failed; aborting.")
            sys.exit(2)

        if not submit_script.exists():
            log(f"Expected {submit_script} after 'next', but it is missing. "
                "Aborting.")
            sys.exit(2)

        round_log = run_dir / f"round_{r:03d}.log"
        log(f"  launching {submit_script.name} (foreground); "
            f"per-round log: {round_log.name}")
        rc = run_cmd("round",
                     ["bash", str(submit_script)],
                     cwd=run_dir, capture_to=round_log)
        if rc != 0:
            log(f"  WARNING: submit_array.sh exited non-zero ({rc}). "
                "Continuing — harvest will pick up whatever finished.")

    log(f"Reached --max-rounds {args.max_rounds} without hitting target. Done.")


if __name__ == "__main__":
    main()
