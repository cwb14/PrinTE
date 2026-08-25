#!/usr/bin/env python3
"""
Merge LTR retrotransposon parts in a RepeatMasker-style FASTA library.

Supersedes TE_lib_stitcher.py, which keys on the first "LTR" or "I" anywhere in
the header and so silently skips names containing an I (ATCOPIA*, for one).

- Headers look like: >[name]#[class]/[superfamily]
- For class == "LTR", sequences may be split into:
    * LTR parts: name contains "LTR" (e.g. ATCOPIA67_LTR)
    * Internal parts: name contains "I" or "INT" (e.g. ATCOPIA67_I or ATCOPIA67_INT)
  The rest of the name (left/right of the marker) is identical between a pair.
- Non-LTR classes pass through unchanged.
- LTRs without a matching internal are assumed intact and pass through unchanged.
- Internals without a matching LTR are, by default, passed through unchanged
  (use --drop-unmatched-internal to drop them).
- When an internal pairs with multiple LTRs, produce one intact element per LTR,
  appending _1, _2, ... to the base name.

Messages about merges and pass-throughs are printed to stderr.
FASTA output is written to stdout (or to --out).

Usage:
    merge_ltr_parts.py --in mylib.fasta > merged.fasta
    merge_ltr_parts.py -i mylib.fasta -o merged.fasta --drop-unmatched-internal
"""
import argparse
import re
import sys
from collections import OrderedDict


def parse_fasta(handle):
    header = None
    seq_lines = []
    for raw in handle:
        line = raw.rstrip('\n')
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield header[1:], ''.join(seq_lines)
            header = line
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        yield header[1:], ''.join(seq_lines)

def parse_header_id(id_str):
    tok = id_str.strip().split()[0]
    if '#' in tok:
        name, rest = tok.split('#',1)
        if '/' in rest:
            cl, superfam = rest.split('/',1)
        else:
            cl, superfam = rest, None
    else:
        name, cl, superfam = tok, None, None
    return tok, name, cl, superfam

def clean_name(n):
    n = re.sub(r'__+', '_', n)
    return n.strip('_')

def classify_name(name):
    # priority: LTR first, then INT, then _I, then standalone I between characters
    if 'LTR' in name:
        base = name.replace('LTR','',1)
        return 'ltr', clean_name(base)
    if 'INT' in name:
        base = name.replace('INT','',1)
        return 'internal', clean_name(base)
    if '_I' in name:
        base = name.replace('_I','',1)
        return 'internal', clean_name(base)
    m = re.search(r'(?<=\w)I(?=\w)', name)
    if m:
        base = name[:m.start()] + name[m.end():]
        return 'internal', clean_name(base)
    return 'unknown', name

def wrap_seq(seq, width=60):
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))

def merge_ltr_parts(in_handle, out_handle, drop_unmatched_internal=False):
    # First pass: collect entries and group candidates by (base_name, class, superfamily)
    entries = []  # list of dicts
    grouping = OrderedDict()  # key -> {'ltr':[idx], 'internal':[idx]}
    key_order = []
    for header, seq in parse_fasta(in_handle):
        tok, name, cl, superfam = parse_header_id(header)
        seg_type = None
        base = None
        if cl == 'LTR':
            seg_type, base = classify_name(name)
        entry = {
            'orig_header': header,
            'tok': tok,
            'name': name,
            'class': cl,
            'superfam': superfam,
            'segment_type': seg_type,
            'base_name': base,
            'seq': seq
        }
        entries.append(entry)
        if cl == 'LTR' and seg_type in ('ltr','internal'):
            key = (base, cl, superfam)
            if key not in grouping:
                grouping[key] = {'ltr': [], 'internal': []}
                key_order.append(key)
            grouping[key][seg_type].append(len(entries)-1)

    consumed = set()
    # Merge pairs
    for key in key_order:
        base, cl, superfam = key
        ltr_idxs = grouping[key]['ltr']
        int_idxs = grouping[key]['internal']
        if not ltr_idxs or not int_idxs:
            continue
        total_pairs = len(ltr_idxs) * len(int_idxs)
        pair_count = 0
        for i_idx in int_idxs:
            for l_idx in ltr_idxs:
                pair_count += 1
                newname = base if total_pairs == 1 else f"{base}_{pair_count}"
                new_header = f"{newname}#{cl}/{superfam}"
                ltr = entries[l_idx]
                inte = entries[i_idx]
                merged_seq = ltr['seq'] + inte['seq'] + ltr['seq']
                # log
                print(f">{inte['tok']} and >{ltr['tok']} appear paired. Merging to >{new_header}.",
                      file=sys.stderr)
                # write
                out_handle.write(f">{new_header}\n{wrap_seq(merged_seq)}\n")
        # mark parts as consumed so we don't pass them through again
        for idx in set(ltr_idxs + int_idxs):
            consumed.add(idx)

    # Pass-through for everything else
    for idx, e in enumerate(entries):
        if idx in consumed:
            continue
        if e['class'] != 'LTR':
            # Non-LTR classes pass through unchanged (silent)
            out_handle.write(f">{e['orig_header']}\n{wrap_seq(e['seq'])}\n")
            continue
        # LTR class
        if e['segment_type'] == 'ltr':
            print(f">{e['tok']} has no matching internal. Passing through unmodified.", file=sys.stderr)
            out_handle.write(f">{e['orig_header']}\n{wrap_seq(e['seq'])}\n")
        elif e['segment_type'] == 'unknown':
            print(f">{e['tok']} appears intact or not partitioned. Passing through unmodified.", file=sys.stderr)
            out_handle.write(f">{e['orig_header']}\n{wrap_seq(e['seq'])}\n")
        elif e['segment_type'] == 'internal':
            if drop_unmatched_internal:
                print(f">{e['tok']} is an internal with no matching LTR. Dropping (use --drop-unmatched-internal off to retain).",
                      file=sys.stderr)
            else:
                print(f">{e['tok']} is an internal with no matching LTR. Passing through unmodified.",
                      file=sys.stderr)
                out_handle.write(f">{e['orig_header']}\n{wrap_seq(e['seq'])}\n")
        else:
            out_handle.write(f">{e['orig_header']}\n{wrap_seq(e['seq'])}\n")

def main():
    ap = argparse.ArgumentParser(description="Merge LTR internal/LTR parts into intact elements in a RepeatMasker-style FASTA.")
    ap.add_argument('-i','--in', dest='inp', required=True, help="Input FASTA (use '-' for stdin)")
    ap.add_argument('-o','--out', dest='out', default='-', help="Output FASTA (default: stdout)")
    ap.add_argument('--drop-unmatched-internal', action='store_true',
                    help="Drop internal-only sequences that have no matching LTR (default: keep them).")
    args = ap.parse_args()

    in_handle = sys.stdin if args.inp == '-' else open(args.inp)
    out_handle = sys.stdout if args.out == '-' else open(args.out, 'w')

    try:
        merge_ltr_parts(in_handle, out_handle, drop_unmatched_internal=args.drop_unmatched_internal)
    finally:
        if in_handle is not sys.stdin:
            in_handle.close()
        if out_handle is not sys.stdout:
            out_handle.close()

if __name__ == '__main__':
    main()
