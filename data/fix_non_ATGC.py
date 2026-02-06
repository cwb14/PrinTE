#!/usr/bin/env python3

import sys
import gzip

ALLOWED = set("ATGCatgc")


def sanitize_sequence(seq: str):
    """
    Replace any non-ATGC character with 'T'.
    Return the cleaned sequence and the count of replaced bases.
    """
    cleaned = []
    replaced = 0
    for ch in seq:
        if ch in ALLOWED:
            cleaned.append(ch)
        else:
            # Ignore whitespace just in case, but sequences shouldn't have it
            if ch.strip() == "":
                continue
            cleaned.append("T")
            replaced += 1
    return "".join(cleaned), replaced


def open_maybe_gzip(path: str):
    """
    Open plain FASTA or gzipped FASTA as text.
    '-' means stdin (assumed plain text).
    """
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def process_fasta_stream(fin, fout):
    """
    Read FASTA from fin, write cleaned FASTA to fout.
    Writes sequences unwrapped (single line per record), like your original script.
    Prints stats to stderr.
    """
    total_seqs = 0
    modified_seqs = 0
    total_replaced_bases = 0

    header = None
    seq_chunks = []

    def flush_record():
        nonlocal header, seq_chunks, total_seqs, modified_seqs, total_replaced_bases
        if header is None:
            return
        total_seqs += 1
        raw_seq = "".join(seq_chunks)
        cleaned_seq, replaced = sanitize_sequence(raw_seq)
        if replaced > 0:
            modified_seqs += 1
            total_replaced_bases += replaced

        fout.write(header + "\n")
        fout.write(cleaned_seq + "\n")

    for line in fin:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush_record()
            header = line
            seq_chunks = []
        else:
            # Strip spaces/tabs; keep no wrapping
            seq_chunks.append(line.strip())

    flush_record()

    print("[STATS]", file=sys.stderr)
    print(f"  Total sequences:       {total_seqs}", file=sys.stderr)
    print(f"  Sequences modified:    {modified_seqs}", file=sys.stderr)
    print(f"  Bases replaced with T: {total_replaced_bases}", file=sys.stderr)


def usage():
    msg = (
        "Usage:\n"
        "  python fix_non_ATGC.py input.fa > output.fa\n"
        "  python fix_non_ATGC.py input.fa.gz > output.fa\n"
        "  cat input.fa | python fix_non_ATGC.py - > output.fa\n"
    )
    print(msg, file=sys.stderr)


def main():
    if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
        usage()
        sys.exit(0 if (len(sys.argv) == 2 and sys.argv[1] in ("-h", "--help")) else 1)

    in_path = sys.argv[1]

    try:
        with open_maybe_gzip(in_path) as fin:
            process_fasta_stream(fin, sys.stdout)
    except BrokenPipeError:
        # Allows piping into `head` etc. without stacktrace
        sys.exit(0)
    except FileNotFoundError:
        print(f"[ERROR] File not found: {in_path}", file=sys.stderr)
        sys.exit(2)
    except gzip.BadGzipFile:
        print(f"[ERROR] Input ends with .gz but is not a valid gzip file: {in_path}", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()
