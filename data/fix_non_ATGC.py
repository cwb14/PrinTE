#!/usr/bin/env python3

import os
import glob
import gzip
import subprocess

ALLOWED = set("ATGCatgc")
SKIP_FILE = "TAIR10.pep.fa.gz"


def sanitize_sequence(seq):
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


def process_fasta_gz(in_gz):
    base_name = os.path.basename(in_gz)
    if base_name == SKIP_FILE:
        print(f"[SKIP] {base_name} is a protein file; leaving it unchanged.")
        return

    # Strip .gz, then add a suffix so we don't overwrite original
    out_plain = in_gz[:-3] + ".ATGCfixed.fa"
    out_gz = out_plain + ".gz"

    if os.path.exists(out_gz):
        print(f"[WARN] Output {out_gz} already exists; skipping to avoid overwrite.")
        return

    print(f"[INFO] Processing {in_gz} -> {out_gz}")

    total_seqs = 0
    modified_seqs = 0
    total_replaced_bases = 0

    with gzip.open(in_gz, "rt") as fin, open(out_plain, "w") as fout:
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
            # Write header and *unwrapped* sequence (single line)
            fout.write(header + "\n")
            fout.write(cleaned_seq + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # New record
                flush_record()
                header = line
                seq_chunks = []
            else:
                # Sequence line; strip spaces/tabs, no wrapping kept
                seq_chunks.append(line.strip())

        # Flush last record
        flush_record()

    print(f"[STATS] {in_gz}")
    print(f"        Total sequences:       {total_seqs}")
    print(f"        Sequences modified:    {modified_seqs}")
    print(f"        Bases replaced with T: {total_replaced_bases}")
    print(f"[INFO] Compressing {out_plain} -> {out_gz} with gzip -9")

    # Compress with gzip -9
    subprocess.run(["gzip", "-9", out_plain], check=True)

    print(f"[DONE] Wrote cleaned file: {out_gz}\n")


def main():
    gz_files = sorted(glob.glob("*.gz"))
    if not gz_files:
        print("[INFO] No .gz files found in current directory.")
        return

    print("[INFO] Found .gz files:")
    for f in gz_files:
        print(f"   - {f}")
    print("")

    for f in gz_files:
        process_fasta_gz(f)


if __name__ == "__main__":
    main()
