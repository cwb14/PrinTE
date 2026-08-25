#!/usr/bin/env python3
"""Aggregate per-simulation napkin_check.tsv (or napkin_estimate.tsv) files
from a grid-search run directory into one calibration table for offline
refinement of the back-of-napkin model. Read-only over the sim dirs."""
import argparse
from pathlib import Path


def aggregate(run_dir: Path, fname: str, out_path: Path) -> int:
    run_dir = Path(run_dir)
    header_written = False
    n = 0
    with open(out_path, "w") as out:
        for sub in sorted(run_dir.iterdir()):
            f = sub / fname
            if not (sub.is_dir() and f.exists()):
                continue
            lines = f.read_text().splitlines()
            if not lines:
                continue
            if not header_written:
                out.write("dir\t" + lines[0] + "\n")
                header_written = True
            for line in lines[1:]:
                if line.strip():
                    out.write(sub.name + "\t" + line + "\n")
                    n += 1
    return n


def main():
    p = argparse.ArgumentParser(description="Aggregate napkin calibration files across a grid run")
    p.add_argument("--run-dir", default=".")
    p.add_argument("--file", default="napkin_check.tsv",
                   help="per-sim file to collect (napkin_check.tsv or napkin_estimate.tsv)")
    p.add_argument("--out", default="napkin_calibration.tsv")
    args = p.parse_args()
    n = aggregate(Path(args.run_dir), args.file, Path(args.out))
    print(f"Aggregated {n} row(s) from {args.file} into {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
