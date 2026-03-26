#!/usr/bin/env python3
"""
scan_results.py  –  Summarise completed PrinTE grid-search runs.

For every directory matching insertion_rates_*_deletion_rates_* in the
working directory, the script:
  1. Finds the highest-numbered gen<N>_final.fasta inside it.
  2. Checks whether that file exists (→ "complete").
  3. Records the file size and compares it to --target-size (±10%).

Outputs three stanzas:
  ALL        – every matching directory
  COMPLETE   – directories that contain the highest gen*_final.fasta
  BALLPARK   – complete dirs whose final fasta is within 10 % of --target-size

Usage:
  python scan_results.py --target-size 211630673
  python scan_results.py --target-size 211630673 --dir /path/to/runs
  python scan_results.py --target-size 211630673 --ball-pct 15
"""

import argparse
import os
import re
import sys
from pathlib import Path


DIR_RE = re.compile(
    r"insertion_rates_([^_]+)_deletion_rates_([^_]+)"
)
GEN_RE = re.compile(r"gen(\d+)_final\.fasta$")


def find_highest_gen_fasta(d: Path):
    """Return (gen_number, Path) for the highest gen*_final.fasta, or None."""
    best_n = -1
    best_path = None
    for f in d.iterdir():
        m = GEN_RE.match(f.name)
        if m:
            n = int(m.group(1))
            if n > best_n:
                best_n = n
                best_path = f
    if best_n == -1:
        return None
    return best_n, best_path


def file_size(p: Path) -> int:
    try:
        return p.stat().st_size
    except OSError:
        return 0


def parse_rates(name: str):
    """Extract (insertion_rate, deletion_rate) floats from directory name."""
    m = DIR_RE.search(name)
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2))


def fmt(v):
    """Format a float for display."""
    if v is None:
        return "NA"
    return f"{v:g}"


class Stats:
    def __init__(self):
        self.count = 0
        self.min_I = None
        self.max_I = None
        self.min_D = None
        self.max_D = None

    def add(self, I, D):
        self.count += 1
        self.min_I = I if self.min_I is None else min(self.min_I, I)
        self.max_I = I if self.max_I is None else max(self.max_I, I)
        self.min_D = D if self.min_D is None else min(self.min_D, D)
        self.max_D = D if self.max_D is None else max(self.max_D, D)

    def print_block(self, label):
        print(f"{label}:")
        print(f"  count = {self.count}")
        print(f"  insertion_rate min = {fmt(self.min_I)}  max = {fmt(self.max_I)}")
        print(f"  deletion_rate  min = {fmt(self.min_D)}  max = {fmt(self.max_D)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target-size", "-ts", type=float, required=True,
                    help="Target genome size in bytes (fasta file size to aim for).")
    ap.add_argument("--ball-pct", "-bp", type=float, default=10.0,
                    help="Percentage tolerance around --target-size for BALLPARK (default: 10).")
    ap.add_argument("--suggest-pct", "-sp", type=float, default=50.0,
                    help="How far to expand BALLPARK bounds for SUGGESTED NEW PARAMETERS, "
                         "as a percentage (default: 50, i.e. min*0.5 / max*1.5).")
    ap.add_argument("--dir", "-d", default=".",
                    help="Directory to scan (default: current directory).")
    args = ap.parse_args()

    scan_root = Path(args.dir).resolve()
    ball_min = args.target_size * (1 - args.ball_pct / 100)
    ball_max = args.target_size * (1 + args.ball_pct / 100)

    all_stats  = Stats()
    ok_stats   = Stats()
    ball_stats = Stats()

    for entry in sorted(scan_root.iterdir()):
        if not entry.is_dir():
            continue
        if not DIR_RE.search(entry.name):
            continue

        I, D = parse_rates(entry.name)
        if I is None:
            continue

        all_stats.add(I, D)

        result = find_highest_gen_fasta(entry)
        if result is None:
            continue  # no gen*_final.fasta at all

        _, fasta_path = result
        fsize = file_size(fasta_path)

        ok_stats.add(I, D)

        if ball_min <= fsize <= ball_max:
            ball_stats.add(I, D)

    # ── Output ────────────────────────────────────────────────────────────────
    all_stats.print_block("ALL directories")
    print()
    ok_stats.print_block("COMPLETE directories (contain highest gen*_final.fasta)")
    print()
    ball_stats.print_block(
        f"BALLPARK (highest gen*_final.fasta size within "
        f"{args.target_size:g} \u00b1 {args.ball_pct:g}%)"
    )

    if ball_stats.count > 0:
        lo = 1 - args.suggest_pct / 100
        hi = 1 + args.suggest_pct / 100
        print()
        print(f"SUGGESTED NEW PARAMETERS (ballpark bounds expanded by ±{args.suggest_pct:g}%):")
        print(f"  insertion_rate min = {fmt(ball_stats.min_I * lo)}  "
              f"max = {fmt(ball_stats.max_I * hi)}")
        print(f"  deletion_rate  min = {fmt(ball_stats.min_D * lo)}  "
              f"max = {fmt(ball_stats.max_D * hi)}")


if __name__ == "__main__":
    main()
