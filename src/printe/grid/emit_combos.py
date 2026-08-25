#!/usr/bin/env python3
"""Write the parameter grid as a CSV, one row per combination.

This exists so the Nextflow sweep can scatter over the same parameter space the
Python front-ends search, without either side reimplementing the spacing. The
log-spacing and the directory-name format both come from grid_utils.
"""

import argparse
import csv
import sys

from .grid_utils import build_logspace_values, dir_name_from_combo

FIELDS = ["ins", "del", "sr", "k", "dirname"]


def int_range(start, end, step):
    if step <= 0:
        raise ValueError("step must be positive")
    lo, hi = (start, end) if start <= end else (end, start)
    values = list(range(lo, hi + 1, step))
    if not values:
        raise ValueError(f"empty range: {start}..{end} by {step}")
    return values


def build_combos(args):
    ins = build_logspace_values(args.ins_start, args.ins_end, args.ins_count)
    dele = build_logspace_values(args.del_start, args.del_end, args.del_count)
    srs = int_range(args.sr_start, args.sr_end, args.sr_step)
    ks = int_range(args.k_start, args.k_end, args.k_step)
    return [
        (i, d, s, k, dir_name_from_combo(i, d, s, k))
        for i in ins
        for d in dele
        for s in srs
        for k in ks
    ]


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ins-start", type=float, required=True,
                   help="Insertion rate at one end of the range")
    p.add_argument("--ins-end", type=float, required=True,
                   help="Insertion rate at the other end")
    p.add_argument("--ins-count", type=int, default=3,
                   help="How many log-spaced insertion rates (default: %(default)s)")
    p.add_argument("--del-start", type=float, required=True,
                   help="Deletion rate at one end of the range")
    p.add_argument("--del-end", type=float, required=True,
                   help="Deletion rate at the other end")
    p.add_argument("--del-count", type=int, default=3,
                   help="How many log-spaced deletion rates (default: %(default)s)")
    p.add_argument("--sr-start", type=int, default=90, help="Solo-LTR percent, low end")
    p.add_argument("--sr-end", type=int, default=95, help="Solo-LTR percent, high end")
    p.add_argument("--sr-step", type=int, default=5, help="Solo-LTR percent step")
    p.add_argument("--k-start", type=int, default=0, help="Length-bias slope, low end")
    p.add_argument("--k-end", type=int, default=10, help="Length-bias slope, high end")
    p.add_argument("--k-step", type=int, default=5, help="Length-bias slope step")
    p.add_argument("--out", default="-", help="Output CSV, or - for stdout (default: -)")
    return p


def main():
    args = build_parser().parse_args()
    combos = build_combos(args)
    handle = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    try:
        writer = csv.writer(handle)
        writer.writerow(FIELDS)
        writer.writerows(combos)
    finally:
        if handle is not sys.stdout:
            handle.close()
    print(f"{len(combos)} parameter combinations", file=sys.stderr)


if __name__ == "__main__":
    main()
