#!/usr/bin/env python3
import sys
import re

USAGE = "Usage: python script.py infile.fa outfile.tsv"

def main():
    if len(sys.argv) != 3:
        print(USAGE, file=sys.stderr)
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    # Compile once for speed; capture the FIRST (leftmost) occurrence only.
    # Pattern tolerates optional spaces after the colon (e.g., LTRlen:123 or LTRlen: 123).
    ltrlen_re = re.compile(r"LTRlen\s*:\s*(\d+)")

    with open(in_path, "r", encoding="utf-8") as fin, open(out_path, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.startswith(">"):
                continue  # sequence lines are irrelevant
            header = line[1:].rstrip("\n").rstrip("\r")  # full header, minus leading ">"
            m = ltrlen_re.search(header)
            if not m:
                # Skip headers with no LTRlen
                continue
            ltrlen = m.group(1)  # first (leftmost) match
            fout.write(f"{header}\t{ltrlen}\n")

if __name__ == "__main__":
    main()
