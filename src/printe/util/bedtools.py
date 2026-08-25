#!/usr/bin/env python3
"""
bedtools.py: Compute overlap stats between SCN/PASS and BED files with reciprocal overlap threshold.

Usage:
    python bedtools.py -pass_scn infile.[pass.list|scn] -bed infile.bed -r 0.9 [-print overlapping]
"""
import argparse
import signal
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
            try:
                start = int(parts[0])
                end = int(parts[1])
            except ValueError:
                continue
            chrom = parts[-1]
            # normalize coordinates
            if start > end:
                start, end = end, start
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
            # Expect something like:
            #   chr:start..end   OR   chr:start-end
            if ':' not in loc:
                continue
            chrom, coords = loc.split(':', 1)
            # Strip annotation suffix (e.g. #LTR/Gypsy/Selgy)
            if '#' in coords:
                coords = coords[:coords.index('#')]

            # Support both ".." and "-" as range separators
            if '..' in coords:
                start_str, end_str = coords.split('..', 1)
            elif '-' in coords:
                start_str, end_str = coords.split('-', 1)
            else:
                # Unknown format
                continue

            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue

            # normalize coordinates
            if start > end:
                start, end = end, start
            entries.append((chrom, start, end, raw))
    return entries


def parse_pass_scn(path):
    entries = []
    with open(path) as f:
        # find first non-comment line
        first_raw = None
        first_line = None
        for line in f:
            raw = line.rstrip("\n")
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            first_raw = raw
            first_line = s
            break

        if first_line is None:
            return []

        first_tok = first_line.split()[0]
        is_pass = (':' in first_tok) and ('..' in first_tok or '-' in first_tok)

        # helper: parse one PASS line (raw + stripped)
        def parse_pass_line(raw, s):
            parts = s.split()
            loc = parts[0]
            if ':' not in loc:
                return None
            chrom, coords = loc.split(':', 1)
            # Strip annotation suffix (e.g. #LTR/Gypsy/Selgy)
            if '#' in coords:
                coords = coords[:coords.index('#')]
            if '..' in coords:
                start_str, end_str = coords.split('..', 1)
            elif '-' in coords:
                start_str, end_str = coords.split('-', 1)
            else:
                return None
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                return None
            if start > end:
                start, end = end, start
            return (chrom, start, end, raw)

        # helper: parse one SCN line
        def parse_scn_line(raw, s):
            parts = s.split()
            if len(parts) < 12:
                return None
            try:
                start = int(parts[0])
                end = int(parts[1])
            except ValueError:
                return None
            chrom = parts[-1]
            if start > end:
                start, end = end, start
            return (chrom, start, end, raw)

        # parse the first line + the rest from the same stream
        parser = parse_pass_line if is_pass else parse_scn_line

        first_parsed = parser(first_raw, first_line)
        if first_parsed:
            entries.append(first_parsed)

        for line in f:
            raw = line.rstrip("\n")
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parsed = parser(raw, s)
            if parsed:
                entries.append(parsed)

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
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            # normalize coordinates
            if start > end:
                start, end = end, start
            entries.append((chrom, start, end, raw))
    return entries


# ---------------------------------------------------------------------------
# Intact LTR-RT subsetting (--printe-intact).
# Classification logic is ported from ../bin/extract_intact_LTR.py, which is
# the canonical source. If the category rules change there, update them here
# too: same first-pass categories, the same "rule 3" nearby-fragment check,
# and the same final filter on intact TEs whose TE_class is "LTR".
# ---------------------------------------------------------------------------
def _extract_te_class(feature_id):
    """TE_class from a feature_id of form TE_name#TE_class/TE_superfamily."""
    try:
        _te_name, rest = feature_id.split("#", 1)
        te_class, _te_superfamily = rest.split("/", 1)
        return te_class
    except ValueError:
        return None


def parse_bed_intact(path):
    """Parse a BED file, returning only intact LTR-RT entries.

    Same (chrom, start, end, raw) tuples as parse_bed(), but restricted to
    records classified as an intact TE whose TE_class is "LTR", following the
    rules in ../bin/extract_intact_LTR.py.
    """
    records = []
    with open(path) as f:
        for line in f:
            raw = line.rstrip("\n")
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) < 6:
                continue  # need name/tsd/strand columns to classify
            chrom, start_s, end_s, name, tsd, strand = parts[:6]
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            records.append({
                "chrom": chrom, "start": start, "end": end, "raw": raw,
                "name": name, "tsd": tsd, "strand": strand,
                "feature_id": name.split(";")[0],
                "additional": name.split(";")[1:],
                "category": None,
            })

    # First pass: initial category per record.
    for rec in records:
        fid = rec["feature_id"]
        if fid.startswith("gene"):
            rec["category"] = "gene"
        elif "_FRAG" in fid:
            rec["category"] = "Fragmented TE"
        elif "_SOLO" in fid:
            rec["category"] = "SoloLTR"
        elif rec["additional"] and "CUT_BY" in rec["additional"][0]:
            rec["category"] = "Fragmented TE"
        else:
            rec["category"] = "Potential intact TE"

    # Rule 3: a potential-intact record carrying extra attributes is really a
    # fragment if a nearby (±100 in file order) non-gene record shares its TSD
    # and strand and a prefix-related name.
    for i, rec in enumerate(records):
        if rec["category"] == "Potential intact TE" and rec["additional"]:
            for j in range(max(0, i - 100), min(len(records), i + 101)):
                if i == j:
                    continue
                other = records[j]
                if other["feature_id"].startswith("gene"):
                    continue
                if other["tsd"] == rec["tsd"] and other["strand"] == rec["strand"]:
                    if rec["name"].startswith(other["name"]) or other["name"].startswith(rec["name"]):
                        rec["category"] = "Fragmented TE"
                        break

    # Surviving potential-intact records are intact; keep LTR-class ones.
    for rec in records:
        if rec["category"] == "Potential intact TE":
            rec["category"] = "Intact TE"

    intact = []
    for rec in records:
        if rec["category"] != "Intact TE":
            continue
        if _extract_te_class(rec["feature_id"]) != "LTR":
            continue
        intact.append((rec["chrom"], rec["start"], rec["end"], rec["raw"]))
    return intact


def compute_metrics(overlapped_scn, overlapped_bed, total_scn, total_bed):
    """Precision/recall/F1/FDR from overlap counts.

    Overlaps are not strictly 1:1, so true positives have two values:
    predicted-side (overlapped_scn) and truth-side (overlapped_bed). Precision
    is computed from the predicted side, recall from the truth side. Ratios
    that are undefined (zero denominator) are returned as None.
    """
    fp = total_scn - overlapped_scn
    fn = total_bed - overlapped_bed
    precision = overlapped_scn / total_scn if total_scn > 0 else None
    recall = overlapped_bed / total_bed if total_bed > 0 else None
    if precision is None or recall is None or (precision + recall) == 0:
        f1 = None
    else:
        f1 = 2 * precision * recall / (precision + recall)
    fdr = (1 - precision) if precision is not None else None
    return {
        "tp_pred": overlapped_scn,
        "tp_truth": overlapped_bed,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "fdr": fdr,
    }


def reciprocal_overlap(a_start, a_end, b_start, b_end):
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    overlap_len = max(0, overlap_end - overlap_start)
    if overlap_len <= 0:
        return 0.0, 0.0
    len_a = a_end - a_start
    len_b = b_end - b_start
    return overlap_len / len_a, overlap_len / len_b


def main():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    parser = argparse.ArgumentParser(description="Compute overlap stats between SCN/PASS and BED files.")
    parser.add_argument('-pass_scn', required=True, help='Input SCN- or PASS-formatted file')
    parser.add_argument('-bed', required=True, help='Input BED-formatted file')
    parser.add_argument('-r', type=float, default=0.0, help='Reciprocal overlap threshold (0-1)')
    parser.add_argument('-print', dest='print_mode', choices=['overlapping', 'unique-input', 'unique-bed'], help='Lines to print')
    parser.add_argument('--printe-intact', dest='printe_intact', action='store_true',
                        help='Subset the BED to intact LTR-RT entries only before '
                             'computing overlap, using PrinTE classification rules '
                             '(see bin/extract_intact_LTR.py). Use when the BED is a '
                             'full genome annotation but the SCN/PASS file contains '
                             'only intact LTR-RTs, so false negatives stay meaningful.')
    args = parser.parse_args()

    in_entries = parse_pass_scn(args.pass_scn)
    if args.printe_intact:
        bed_entries = parse_bed_intact(args.bed)
        print(f"Subset BED via --printe-intact: {len(bed_entries)} intact LTR-RT entries")
    else:
        bed_entries = parse_bed(args.bed)
    scn_entries = in_entries

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

    # Summary stats. The two overlap counts differ because matches are not
    # strictly 1:1; they serve as predicted-side and truth-side TP respectively.
    m = compute_metrics(overlapped_scn, overlapped_bed, total_scn, total_bed)

    def fmt(x):
        return f"{x:.4f}" if x is not None else "NA"

    print(f"Overlapping entries: {overlapped_scn} ({overlapped_bed} unique)        [TP: {overlapped_scn} predicted-side / {overlapped_bed} truth-side]")
    print(f"Entries unique to SCN/PASS file: {unique_scn}        [FP]")
    print(f"Entries unique to BED file: {unique_bed}        [FN]")
    print()
    print(f"Precision = {fmt(m['precision'])}   (TP_pred {overlapped_scn} / predictions {total_scn})")
    print(f"Recall    = {fmt(m['recall'])}   (TP_truth {overlapped_bed} / truths {total_bed})   [= sensitivity]")
    print(f"F1        = {fmt(m['f1'])}")
    print(f"FDR       = {fmt(m['fdr'])}   (1 - precision)")

    # Optional detailed printing
    if args.print_mode == 'overlapping':
        for idx_bed in sorted(bed_to_scn):
            print('--')
            print(f"bed: {bed_entries[idx_bed][3]}")
            for idx_scn in bed_to_scn[idx_bed]:
                print(f"in: {scn_entries[idx_scn][3]}")
    elif args.print_mode == 'unique-input':
        for idx in sorted(set(range(total_scn)) - matched_scn):
            print(f"in: {scn_entries[idx][3]}")
    elif args.print_mode == 'unique-bed':
        for idx in sorted(set(range(total_bed)) - matched_bed):
            print(f"bed: {bed_entries[idx][3]}")

if __name__ == '__main__':
    main()
