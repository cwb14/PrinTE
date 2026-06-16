#!/usr/bin/env python3
"""
Unattended driver for the guided-search loop.

Each round runs:
  1. guided_search.py harvest --run-dir .
  2. (count hits in training_data.tsv; stop if >= --target-hits)
  3. purge_intermediates.py --execute     (skipped with --no-purge)
  4. guided_search.py next --samples N --explore-frac F
  5. Launch sims — mode depends on scheduler detected from search_state.json:
       local  : bash slurm_grid/submit_array.sh  (FOREGROUND — blocks until done)
       slurm  : sbatch submit_array.sbatch (or submit_chunks.sh), then polls
                squeue until no active tasks remain for the job name.

Scheduler detection
-------------------
The scheduler is read from search_state.json → slurm_args.scheduler.  If the
key is absent the scheduler defaults to 'local' (preserving behaviour for any
run directory created before SLURM support was added).  In slurm mode the
script verifies that sbatch and squeue are on PATH at startup and aborts with
a clear message if they are not.

Re-launching / resumability
----------------------------
Re-launching the script is safe: on startup it reconciles the current round
before the first harvest.  In slurm mode, reconcile_round_slurm attaches to an
in-flight array job, skips a round that already completed, submits one that was
never launched, or proceeds directly to harvest when tasks ran with partial
failures.  In local mode the script waits for any in-progress pipeline.log
files to reach their completion marker before harvesting.

State lives in search_state.json, training_data.tsv, and the on-disk param dirs.

Run it under nohup so it survives ssh logout:

    nohup python ../PrinTE/grid/run_loop.py \\
        --target-hits 200 \\
        --max-rounds 50 \\
        --samples 100 --explore-frac 0.1 \\
        > orchestrator.log 2>&1 &

Tail orchestrator.log to monitor. To stop early: kill the orchestrator PID.
In slurm mode the already-submitted array job will run to completion, but no
new round will be launched.
"""

import argparse
import csv
import getpass
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

_HERE = Path(__file__).resolve().parent
GUIDED_SEARCH = _HERE / "guided_search.py"
PURGE = _HERE / "purge_intermediates.py"

# grid_utils is import-light (stdlib only); reuse its dir-name format so run_loop
# and guided_search agree on parameter directory names.
sys.path.insert(0, str(_HERE))
from grid_utils import dir_name_from_combo  # noqa: E402

_DIR_RE = re.compile(
    r"^insertion_rates_[^_]+(?:e[+-]?\d+)?"
    r"_deletion_rates_[^_]+(?:e[+-]?\d+)?"
    r"_solo_ratio_\d+_length_bias_\d+$"
)
_FINISHED_MARKER = "Pipeline completed at"


def read_search_state(run_dir: Path) -> dict:
    """Parse search_state.json from the run directory."""
    with open(run_dir / "search_state.json") as fh:
        return json.load(fh)


def detect_scheduler(state: dict) -> str:
    """Return the scheduler ('slurm' or 'local') from persisted state.

    Defaults to 'local' when the key is absent. This is intentional and differs
    from guided_search.py's legacy 'slurm' default: run_loop.py was local-only
    before SLURM support existed, so any state file lacking the 'scheduler' key
    was necessarily a local run (SLURM rounds were driven manually, never by
    run_loop). Defaulting legacy files to 'local' preserves run_loop's history;
    do not change this to 'slurm'.
    """
    return state.get("slurm_args", {}).get("scheduler", "local")


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


def read_combos(combos_tsv: Path):
    """Parse slurm_grid/combos.tsv -> list of (ir, dr, sr, lb) (ir/dr str, sr/lb int)."""
    combos = []
    with open(combos_tsv) as fh:
        next(fh, None)  # skip header
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            ir, dr, sr, lb = line.split("\t")
            combos.append((ir, dr, int(sr), int(lb)))
    return combos


def round_progress(run_dir: Path, combos: list) -> tuple:
    """Return (n_expected, n_started, n_done) for the given round's combos.

    n_started = param dirs that exist; n_done = those whose pipeline.log is finished.
    """
    n_expected = len(combos)
    n_started = 0
    n_done = 0
    for ir, dr, sr, lb in combos:
        d = run_dir / dir_name_from_combo(ir, dr, sr, lb)
        if d.is_dir():
            n_started += 1
            if is_finished(d / "pipeline.log"):
                n_done += 1
    return n_expected, n_started, n_done


def squeue_active(job_name: str) -> int:
    """Count this user's active array tasks for job_name.

    Returns the count (>=0), or -1 on a transient squeue error. Raises SystemExit
    if squeue is not installed.
    """
    user = getpass.getuser()
    try:
        res = subprocess.run(
            ["squeue", "-h", "-r", "-u", user, "-n", job_name, "-o", "%i"],
            capture_output=True, text=True, check=False,
        )
    except FileNotFoundError:
        raise SystemExit(
            "scheduler is 'slurm' but `squeue` was not found on PATH. Run run_loop.py "
            "where Slurm client commands are available (e.g. a login node)."
        )
    if res.returncode != 0:
        log(f"  WARNING: squeue exited {res.returncode}: {res.stderr.strip()}")
        return -1
    return sum(1 for ln in res.stdout.splitlines() if ln.strip())


def sbatch_submit(run_dir: Path, sbatch_rel: str) -> str:
    """sbatch --parsable the array script from run_dir; return the job id."""
    try:
        res = subprocess.run(
            ["sbatch", "--parsable", sbatch_rel],
            cwd=str(run_dir), capture_output=True, text=True, check=False,
        )
    except FileNotFoundError:
        raise SystemExit(
            "scheduler is 'slurm' but `sbatch` was not found on PATH. Run run_loop.py "
            "where Slurm client commands are available (e.g. a login node)."
        )
    if res.returncode != 0:
        raise SystemExit(f"sbatch failed (exit {res.returncode}): {res.stderr.strip()}")
    job_id = res.stdout.strip().split(";")[0]
    log(f"  submitted array job {job_id} ({sbatch_rel})")
    return job_id


def wait_for_slurm_round(job_name: str, poll_seconds: int) -> None:
    """Block until no active Slurm tasks remain for job_name."""
    # Initial sleep so a just-submitted job is registered before the first check.
    time.sleep(poll_seconds)
    unknown_streak = 0
    waited = poll_seconds
    while True:
        n = squeue_active(job_name)
        if n == 0:
            log(f"  slurm: no active '{job_name}' tasks remain; round complete "
                f"(waited ~{waited}s).")
            return
        if n < 0:
            unknown_streak += 1
            if unknown_streak >= 5:
                raise SystemExit(
                    "squeue failed 5 times in a row; cannot track job state. "
                    "Investigate the scheduler, then re-launch run_loop.py to resume."
                )
            log(f"  slurm: squeue unavailable (streak {unknown_streak}); "
                f"retrying in {poll_seconds}s")
        else:
            unknown_streak = 0
            # A stable active count is normal for SLURM (tasks queue/run for hours),
            # so we report cumulative elapsed time instead of warning on no change.
            log(f"  slurm: {n} active '{job_name}' task(s) after ~{waited}s; "
                f"sleeping {poll_seconds}s")
        time.sleep(poll_seconds)
        waited += poll_seconds


def launch_round_slurm(run_dir: Path, slurm_dir: str, job_name: str,
                       poll_seconds: int) -> None:
    """Submit the current round's array job (single or chunked) and wait for it."""
    chunks = run_dir / slurm_dir / "submit_chunks.sh"
    if chunks.exists():
        log(f"  slurm: launching chunked submitter {slurm_dir}/submit_chunks.sh")
        rc = run_cmd("chunks", ["bash", f"{slurm_dir}/submit_chunks.sh"],
                     cwd=run_dir)
        if rc != 0:
            raise SystemExit(f"submit_chunks.sh failed (exit {rc})")
    else:
        sbatch_submit(run_dir, f"{slurm_dir}/submit_array.sbatch")
    wait_for_slurm_round(job_name, poll_seconds)


def reconcile_round_slurm(run_dir: Path, slurm_dir: str, job_name: str,
                          poll_seconds: int) -> None:
    """Ensure the current round is complete before the loop's first harvest.

    Idempotent across controller restarts: attaches to an in-flight array, leaves a
    finished one alone, or submits one that was never launched.
    """
    combos_tsv = run_dir / slurm_dir / "combos.tsv"
    if not combos_tsv.exists():
        raise SystemExit(
            f"{combos_tsv} not found; run guided_search.py init (or next) first."
        )
    combos = read_combos(combos_tsv)
    n_exp, n_started, n_done = round_progress(run_dir, combos)

    active = squeue_active(job_name)
    tries = 0
    while active < 0 and tries < 5:
        time.sleep(poll_seconds)
        active = squeue_active(job_name)
        tries += 1
    if active < 0:
        raise SystemExit("squeue repeatedly failed at startup; cannot reconcile state.")

    if active > 0:
        log(f"reconcile: round combos={n_exp}, dirs={n_started}, finished={n_done}, "
            f"squeue active={active} -> attaching to in-flight job")
        wait_for_slurm_round(job_name, poll_seconds)
    elif n_done >= n_exp:
        log(f"reconcile: round combos={n_exp}, finished={n_done}, squeue empty "
            f"-> already complete; proceeding to harvest")
    elif n_started == 0:
        log(f"reconcile: round combos={n_exp}, dirs=0, squeue empty "
            f"-> not yet submitted; submitting now")
        launch_round_slurm(run_dir, slurm_dir, job_name, poll_seconds)
    else:
        log(f"reconcile: round combos={n_exp}, dirs={n_started}, finished={n_done}, "
            f"squeue empty -> ran with failures; proceeding to harvest (partial)")


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
    state = read_search_state(run_dir)
    scheduler = detect_scheduler(state)
    slurm_args = state.get("slurm_args", {})
    slurm_dir = slurm_args.get("slurm_dir", "slurm_grid")
    job_name = slurm_args.get("slurm_job_name", "printe_grid")

    if scheduler == "slurm":
        for tool in ("sbatch", "squeue"):
            if shutil.which(tool) is None:
                ap.error(
                    f"scheduler is 'slurm' but `{tool}` is not on PATH. Launch "
                    f"run_loop.py where Slurm client commands exist (e.g. a login node)."
                )

    submit_script = run_dir / slurm_dir / "submit_array.sh"  # local launcher only
    tsv_path = run_dir / "training_data.tsv"

    # Quiet SIGHUP so backgrounded launches survive logout cleanly.
    signal.signal(signal.SIGHUP, signal.SIG_IGN)

    log(f"=== Orchestrator starting (PID {os.getpid()}) ===")
    log(f"  run_dir:       {run_dir}")
    log(f"  scheduler:     {scheduler}")
    log(f"  target_hits:   {args.target_hits}")
    log(f"  max_rounds:    {args.max_rounds}")
    log(f"  samples/round: {args.samples}")
    log(f"  explore_frac:  {args.explore_frac}")
    log(f"  purge:         {'OFF (--no-purge)' if args.no_purge else 'ON'}")

    if scheduler == "slurm":
        reconcile_round_slurm(run_dir, slurm_dir, job_name, args.poll_seconds)
    else:
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

        # Run purge before the target-reached check so the final round's
        # intermediates get cleaned up too — otherwise the last 100 sims keep
        # ~14 GB each of gen*_final.{fasta,bed,lib,mut.txt} on disk forever.
        if not args.no_purge:
            # Don't abort on purge failure — it's housekeeping, not pipeline.
            run_cmd("purge",
                    [sys.executable, str(PURGE),
                     "--run-dir", ".", "--execute"],
                    cwd=run_dir)

        if hits >= args.target_hits:
            log("Target hit count reached. Done.")
            return

        if run_cmd("next",
                   [sys.executable, str(GUIDED_SEARCH), "next",
                    "--samples", str(args.samples),
                    "--explore-frac", str(args.explore_frac)],
                   cwd=run_dir) != 0:
            log("'next' failed; aborting.")
            sys.exit(2)

        if scheduler == "slurm":
            sbatch_path = run_dir / slurm_dir / "submit_array.sbatch"
            chunks_path = run_dir / slurm_dir / "submit_chunks.sh"
            if not sbatch_path.exists() and not chunks_path.exists():
                log(f"Expected {sbatch_path} (or {chunks_path.name}) after 'next', "
                    f"but neither exists. Aborting.")
                sys.exit(2)
            log(f"  launching slurm round {r} and waiting for completion")
            launch_round_slurm(run_dir, slurm_dir, job_name, args.poll_seconds)
            # n_done counts sims whose pipeline.log shows completion. A task that
            # failed/timed out leaves the squeue queue without that marker, so a
            # shortfall here correctly flags genuinely-incomplete sims.
            combos = read_combos(run_dir / slurm_dir / "combos.tsv")
            n_exp, _, n_done = round_progress(run_dir, combos)
            if n_done < n_exp:
                log(f"  WARNING: round {r}: {n_done}/{n_exp} sims completed; "
                    f"harvesting partial results next round.")
        else:
            if not submit_script.exists():
                log(f"Expected {submit_script} after 'next', but it is missing. "
                    "Aborting.")
                sys.exit(2)
            round_log = run_dir / f"round_{r:03d}.log"
            log(f"  launching {submit_script.name} (foreground); "
                f"per-round log: {round_log.name}")
            rc = run_cmd("round", ["bash", str(submit_script)],
                         cwd=run_dir, capture_to=round_log)
            if rc != 0:
                log(f"  WARNING: submit_array.sh exited non-zero ({rc}). "
                    "Continuing — harvest will pick up whatever finished.")

    log(f"Reached --max-rounds {args.max_rounds} without hitting target. Done.")


if __name__ == "__main__":
    main()
