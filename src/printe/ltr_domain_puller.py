#!/usr/bin/env python3
import argparse
import re
import sys

LTRLEN_RE = re.compile(r"LTRlen\s*:\s*(\d+)")

def detect_input_type(path):
    """
    Return 'fasta' if any non-empty line starts with '>',
    otherwise 'tsv'. This is robust to leading blank/comment lines.
    """
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith(">"):
                    return "fasta"
                # First meaningful line isn't a FASTA header -> treat as TSV
                return "tsv"
    except FileNotFoundError:
        print(f"Error: cannot open '{path}'", file=sys.stderr)
        sys.exit(1)
    # Default to TSV if file is empty
    return "tsv"

def process_fasta(in_path, fout):
    with open(in_path, encoding="utf-8") as fin:
        for line in fin:
            if not line.startswith(">"):
                continue  # sequence lines are irrelevant
            header = line[1:].rstrip("\n").rstrip("\r")  # minus leading ">"
            m = LTRLEN_RE.search(header)
            if not m:
                continue  # skip headers with no LTRlen
            ltrlen = m.group(1)
            # For FASTA input, output the entire header (as in the original)
            fout.write(f"{header}\t{ltrlen}\n")

def process_tsv(in_path, fout):
    with open(in_path, encoding="utf-8") as fin:
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                # Not enough columns; skip
                continue
            col1 = parts[0]  # transpose column1 exactly
            col2 = parts[1]  # search only in column2
            m = LTRLEN_RE.search(col2)
            if not m:
                # Skip rows with no LTRlen in column2
                continue
            ltrlen = m.group(1)  # leftmost occurrence due to regex .search
            fout.write(f"{col1}\t{ltrlen}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Pull LTRlen values out of a FASTA or key TSV into a two-column domain file."
    )
    parser.add_argument("infile", help="Input .fa or .tsv carrying LTRlen: tags")
    parser.add_argument("outfile", help="Output two-column TSV: name, LTR length")
    args = parser.parse_args()

    in_path, out_path = args.infile, args.outfile
    mode = detect_input_type(in_path)

    with open(out_path, "w", encoding="utf-8") as fout:
        if mode == "fasta":
            process_fasta(in_path, fout)
        else:
            process_tsv(in_path, fout)

if __name__ == "__main__":
    main()
