#!/usr/bin/env python3
"""
bedtools.py: Compute overlap stats between SCN/PASS and BED files with reciprocal overlap threshold.

Usage:
    python bedtools.py -scn infile.scn -bed infile.bed -r 0.9 [-print overlapping]
    python bedtools.py -pass infile.pass.list -bed infile.bed -r 0.9 [-print unique-scn/pass]
"""
import argparse
from collections import defaultdict

def parse_scn(path):
    entries = []
    with open(path) as f:
        for line in f:
            raw = line.rstrip("\n")
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            start = int(parts[0])
            end = int(parts[1])
            chrom = parts[11]
            entries.append((chrom, start, end, raw))
    return entries

def parse_pass(path):
    entries = []
    with open(path) as f:
        for line in f:
            raw = line.rstrip("\n")
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            loc = parts[0]
            try:
                chrom, coords = loc.split(':')
                start_str, end_str = coords.split('..')
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue
            entries.append((chrom, start, end, raw))
    return entries

def parse_bed(path):
    entries = []
    with open(path) as f:
        for line in f:
            raw = line.rstrip("\n")
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            entries.append((chrom, start, end, raw))
    return entries

def reciprocal_overlap(a_start, a_end, b_start, b_end):
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    overlap_len = max(0, overlap_end - overlap_start)
    if overlap_len == 0:
        return 0.0, 0.0
    len_a = a_end - a_start
    len_b = b_end - b_start
    return overlap_len / len_a, overlap_len / len_b

def main():
    parser = argparse.ArgumentParser(description="Compute overlap stats between SCN/PASS and BED files.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-scn', help='Input SCN-formatted file')
    group.add_argument('-pass', dest='passfile', help='Input PASS-formatted file')
    parser.add_argument('-bed', required=True, help='Input BED-formatted file')
    parser.add_argument('-r', type=float, default=0.0, help='Reciprocal overlap threshold (0-1)')
    parser.add_argument('-print', dest='print_mode', choices=['overlapping', 'unique-scn/pass', 'unique-bed'], help='Lines to print')
    args = parser.parse_args()

    # Load entries
    scn_entries = parse_scn(args.scn) if args.scn else parse_pass(args.passfile)
    bed_entries = parse_bed(args.bed)

    # Index by chromosome
    scn_by_chrom = defaultdict(list)
    for idx, (chrom, start, end, _) in enumerate(scn_entries):
        scn_by_chrom[chrom].append((start, end, idx))

    bed_by_chrom = defaultdict(list)
    for idx, (chrom, start, end, _) in enumerate(bed_entries):
        bed_by_chrom[chrom].append((start, end, idx))

    for chrom in scn_by_chrom:
        scn_by_chrom[chrom].sort(key=lambda x: x[0])
    for chrom in bed_by_chrom:
        bed_by_chrom[chrom].sort(key=lambda x: x[0])

    matched_scn = set()
    matched_bed = set()
    bed_to_scn = defaultdict(list)

    # find reciprocally overlapping entries
    for chrom, scn_list in scn_by_chrom.items():
        if chrom not in bed_by_chrom:
            continue
        bed_list = bed_by_chrom[chrom]
        for start_s, end_s, idx_scn in scn_list:
            for start_b, end_b, idx_bed in bed_list:
                ro_scn, ro_bed = reciprocal_overlap(start_s, end_s, start_b, end_b)
                if ro_scn >= args.r and ro_bed >= args.r:
                    matched_scn.add(idx_scn)
                    matched_bed.add(idx_bed)
                    bed_to_scn[idx_bed].append(idx_scn)
                    break

    total_scn = len(scn_entries)
    total_bed = len(bed_entries)
    overlapped_scn = len(matched_scn)
    overlapped_bed = len(matched_bed)
    unique_scn = total_scn - overlapped_scn
    unique_bed = total_bed - overlapped_bed

    # Summary stats
    print(f"Overlapping entries (SCN matched): {overlapped_scn} ({overlapped_bed} unique)")
    print(f"Entries unique to SCN/PASS file: {unique_scn}")
    print(f"Entries unique to BED file: {unique_bed}")

    # Optional detailed printing
    if args.print_mode == 'overlapping':
        for idx_bed in sorted(bed_to_scn):
            print('--')
            print(f"bed: {bed_entries[idx_bed][3]}")
            for idx_scn in bed_to_scn[idx_bed]:
                print(f"scn: {scn_entries[idx_scn][3]}")
    elif args.print_mode == 'unique-scn/pass':
        for idx in sorted(set(range(total_scn)) - matched_scn):
            print(f"scn: {scn_entries[idx][3]}")
    elif args.print_mode == 'unique-bed':
        for idx in sorted(set(range(total_bed)) - matched_bed):
            print(f"bed: {bed_entries[idx][3]}")

if __name__ == '__main__':
    main()
