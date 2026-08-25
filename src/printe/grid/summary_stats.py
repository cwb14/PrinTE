#!/usr/bin/env python3

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Bin values and report total count per bin."
    )
    parser.add_argument(
        "--bins",
        type=int,
        required=True,
        help="Total number of bins (>= 2). Bin 0 is always 0.0000000–0.0000000."
    )
    parser.add_argument(
        "--bin_max",
        type=float,
        required=True,
        help="Maximum value to consider; values larger than this are ignored."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input file with one value per line, or '-' for stdin."
    )
    return parser.parse_args()


def read_values(path, bin_max):
    if path == "-":
        fh = sys.stdin
    else:
        fh = open(path)

    vals = []
    total_read = 0

    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            v = float(line)
        except ValueError:
            continue

        total_read += 1
        if v <= bin_max:
            vals.append(v)

    if fh is not sys.stdin:
        fh.close()

    return vals, total_read


def main():
    args = parse_args()

    if args.bins < 2:
        sys.stderr.write("ERROR: --bins must be at least 2.\n")
        sys.exit(1)

    bin_max = args.bin_max
    zero_bin_start = 0.0
    zero_bin_end = 0.0
    nonzero_start = 0.0000001

    values, total_read = read_values(args.input, bin_max)
    n_used = len(values)

    # Prepare bins
    n_nonzero_bins = args.bins - 1
    if bin_max <= nonzero_start:
        width = 0.0
    else:
        width = (bin_max - nonzero_start) / n_nonzero_bins

    bin_counts = [0] * args.bins

    for v in values:
        if v == 0.0:
            bin_counts[0] += 1
        else:
            if bin_max <= nonzero_start or width <= 0:
                idx = args.bins - 1
            else:
                if v < nonzero_start:
                    idx = 1
                else:
                    pos = (v - nonzero_start) / width
                    idx = int(pos) + 1
                    if idx >= args.bins:
                        idx = args.bins - 1
            bin_counts[idx] += 1

    # Output
    print(f"# total_values_read\t{total_read}")
    print(f"# total_values_used_(<=bin_max)\t{n_used}")
    print("# bin_index\tbin_start\tbin_end\tcount")

    # Bin 0
    print(f"0\t{zero_bin_start:.7f}\t{zero_bin_end:.7f}\t{bin_counts[0]}")

    # Other bins
    for i in range(1, args.bins):
        if width > 0:
            start = nonzero_start + (i - 1) * width
            end = nonzero_start + i * width
        else:
            start = nonzero_start
            end = nonzero_start

        print(
            f"{i}\t{start:.7f}\t{end:.7f}\t{bin_counts[i]}"
        )


if __name__ == "__main__":
    main()
