#!/usr/bin/env python3
"""
Purge intermediate files from finished PrinTE guided-search simulations.

Keeps only what `guided_search.py harvest` needs to rebuild training_data.tsv:
  - The parameter directory itself (its name encodes ins/del/sr/k; the
    `[SKIP] {dname} already exists` check in run_one_combo prevents re-runs)
  - pipeline.log (always — the regression projection is parsed from here for
    failed sims)
  - pipeline.report, pipeline.error (small; retained for debugging)
  - gen{ge}_final.fasta and gen{ge}_final.bed/.lib for HITS (the terminal
    fasta's file size IS the genome_size measurement; bed/lib are the data)

Skips dirs whose pipeline.log does NOT contain "Pipeline completed at" — those
are still running or were killed mid-flight.

Default mode is dry-run; pass --execute to actually delete.

Usage
-----
    # See what would be freed (safe):
    python purge_intermediates.py --run-dir .

    # Actually delete:
    python purge_intermediates.py --run-dir . --execute

    # After every harvest in your loop:
    python guided_search.py harvest --run-dir ./
    python ../PrinTE/grid/purge_intermediates.py --run-dir . --execute
    python guided_search.py next --samples 100 --explore-frac 0.1
"""

import argparse
import json
import re
import sys
from pathlib import Path

# Same regex used by guided_search.parse_dir_name
_DIR_RE = re.compile(
    r"^insertion_rates_(?P<ins>[^_]+(?:e[+-]?\d+)?)"
    r"_deletion_rates_(?P<del>[^_]+(?:e[+-]?\d+)?)"
    r"_solo_ratio_(?P<sr>\d+)"
    r"_length_bias_(?P<k>\d+)$"
)

# Files always kept regardless of hit/miss
_KEEP_ALWAYS = {"pipeline.log", "pipeline.report", "pipeline.error"}

# fnmatch-style globs for files safe to delete from FAILED sims
_PURGE_GLOBS = (
    "gen*_final.fasta",
    "gen*_final.bed",
    "gen*_final.lib",
    "gen*_mut.txt",
    "TAIR10.cds.fa",
    "weighted_candidate_selection.pdf",
)

_FINISHED_MARKER = "Pipeline completed at"


def is_finished(log_path: Path) -> bool:
    """True if pipeline.log contains the completion marker (success OR failure)."""
    if not log_path.exists():
        return False
    try:
        with open(log_path, "rb") as fh:
            try:
                fh.seek(-4096, 2)
            except OSError:
                fh.seek(0)
            tail = fh.read().decode(errors="ignore")
    except OSError:
        return False
    return _FINISHED_MARKER in tail


def purge_dir(d: Path, terminal_fasta: str, dry_run: bool, verbose: bool):
    """Delete intermediate files in a finished sim dir.

    For hits (terminal fasta exists), all gen{ge}_final.* files are preserved
    in addition to pipeline.{log,report,error}. For failed sims, all gen*
    files are deleted.

    Returns (bytes_freed, files_deleted, status) where status is "hit" or "failed".
    """
    is_hit = (d / terminal_fasta).exists()
    terminal_prefix = terminal_fasta.rsplit(".", 1)[0]  # e.g. "gen5400000_final"

    bytes_freed = 0
    files_deleted = 0

    for f in d.iterdir():
        if not f.is_file():
            continue
        name = f.name
        if name in _KEEP_ALWAYS:
            continue
        # For hits: preserve every terminal-generation file (.fasta, .bed, .lib, etc.)
        if is_hit and name.startswith(terminal_prefix + "."):
            continue
        if not any(f.match(g) for g in _PURGE_GLOBS):
            continue

        sz = f.stat().st_size
        if not dry_run:
            try:
                f.unlink()
            except OSError as e:
                print(f"  ERROR: could not delete {f}: {e}", file=sys.stderr)
                continue
        bytes_freed += sz
        files_deleted += 1
        if verbose:
            verb = "would delete" if dry_run else "deleted"
            print(f"  {verb}: {name} ({sz / 1e6:.1f} MB)")

    return bytes_freed, files_deleted, ("hit" if is_hit else "failed")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--run-dir", default=".",
                    help="Parent dir containing parameter subdirs (default: cwd)")
    ap.add_argument("--ge", type=int, default=None,
                    help="Terminal generation. Auto-detected from search_state.json if omitted.")
    ap.add_argument("--execute", action="store_true",
                    help="Actually delete files. Default is dry-run.")
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="List every file action.")
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.is_dir():
        ap.error(f"--run-dir does not exist: {run_dir}")

    ge = args.ge
    if ge is None:
        state_path = run_dir / "search_state.json"
        if not state_path.exists():
            ap.error("Cannot determine --ge. Pass it explicitly or run from a "
                     "directory containing search_state.json.")
        ge = json.loads(state_path.read_text())["ge"]

    terminal_fasta = f"gen{ge}_final.fasta"
    print(f"Run dir:        {run_dir}")
    print(f"Terminal fasta: {terminal_fasta}")
    print(f"Mode:           {'EXECUTE (deletions live)' if args.execute else 'DRY-RUN (no changes)'}")
    print()

    total_bytes = 0
    total_files = 0
    n_hit = 0
    n_failed = 0
    n_skip_running = 0
    n_skip_unparseable = 0

    for entry in sorted(run_dir.iterdir()):
        if not entry.is_dir():
            continue
        if not _DIR_RE.match(entry.name):
            continue

        log_path = entry / "pipeline.log"
        if not is_finished(log_path):
            n_skip_running += 1
            if args.verbose:
                print(f"  SKIP (in-progress / no completion marker): {entry.name}")
            continue

        try:
            bytes_freed, files_deleted, status = purge_dir(
                entry, terminal_fasta,
                dry_run=not args.execute,
                verbose=args.verbose,
            )
        except OSError as e:
            print(f"  ERROR scanning {entry.name}: {e}", file=sys.stderr)
            n_skip_unparseable += 1
            continue

        total_bytes += bytes_freed
        total_files += files_deleted
        if status == "hit":
            n_hit += 1
        else:
            n_failed += 1

    print()
    print(f"Scanned:   {n_hit} hit, {n_failed} failed, "
          f"{n_skip_running} in-progress (skipped)"
          + (f", {n_skip_unparseable} errors" if n_skip_unparseable else ""))
    verb = "Freed" if args.execute else "Would free"
    print(f"{verb}:     {total_bytes / 1e9:.2f} GB across {total_files} files")
    if not args.execute and total_files:
        print()
        print("Re-run with --execute to actually delete.")


if __name__ == "__main__":
    main()
