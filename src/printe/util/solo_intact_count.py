#!/usr/bin/env python3
"""Count solo-LTRs vs intact LTR-RTs in PrinTE BED files and report the ratio.

Definitions are taken directly from PrinTE's own code so the numbers match
what PrinTE itself would report:

  * Intact LTR-RT  -- as in printe/extract_intact_LTR.py:
      non-gene, feature_ID has NO "_FRAG", NO "_SOLO",
      TE_class == "LTR", and the element is NOT *disrupted* (see below).

  * SoloLTR        -- as in printe/util/plot_solo_intact.py, refined:
      feature_ID has "_SOLO", NO "_FRAG", TE_class == "LTR",
      and the solo remnant is NOT *disrupted* by a later TE insertion.

  * "disrupted"    -- exactly printe/TE_exciser.py
                      qualifies_as_intact_ltr():
      the first ';' attribute contains "CUT_BY", OR another record
      within +/-100 BED records has identical TSD and strand and a
      NAME that is an exact-or-prefix match (the element was split in
      two by a nested insertion).  Disruption is only tested for
      records that carry at least one ';' supplemental attribute,
      mirroring extract_intact_LTR.py.

  * solo:intact ratio -- solo / intact (as in plot_solo_intact.py;
      defined as 0.0 when both are 0; reported as +inf with a warning
      if solo > 0 while intact == 0).

The +/-100 neighbour test is evaluated over a bounded sliding window, so
arbitrarily large BED files stream in O(window) memory -- the BED is never
loaded whole.

Stdlib only (argparse, re, os, sys, collections).
"""
import argparse
import os
import re
import sys
from collections import defaultdict, deque

NEIGHBORHOOD = 100  # records before/after, inclusive -- matches PrinTE


# --------------------------------------------------------------------------- #
# Parsing helpers (kept byte-identical in spirit to PrinTE's parsers)
# --------------------------------------------------------------------------- #
def parse_attributes(name):
    """Split BED NAME into (feature_id, [supplemental attributes])."""
    if ";" in name:
        parts = name.split(";")
        return parts[0], parts[1:]
    return name, []


def extract_te_class(feature_id):
    """feature_ID 'te_name#TE_class/TE_superfamily[~junk]' -> (class, super).

    Returns (None, None) when there is no '#': genes and malformed IDs.
    """
    try:
        _, rest = feature_id.split("#", 1)
        te_class, te_super = rest.split("/", 1)
        return te_class, te_super.split("~")[0]
    except ValueError:
        return None, None


def extract_generation(path):
    """Generation label for a PrinTE BED, mirroring plot_solo_intact.py.

    'burnin.bed' -> 0 ; 'gen<N>_final*' -> N ; otherwise None.
    """
    base = os.path.basename(path)
    if base.lower() == "burnin.bed":
        return 0
    m = re.search(r"gen(\d+)_final", base)
    return int(m.group(1)) if m else None


class Rec:
    """Minimal per-line record; only the fields the rules touch."""
    __slots__ = ("name", "tsd", "strand", "feature_id", "supp")

    def __init__(self, name, tsd, strand):
        self.name = name
        self.tsd = tsd
        self.strand = strand
        self.feature_id, self.supp = parse_attributes(name)


def is_disrupted(rec, window, center):
    """PrinTE's qualifies_as_intact_ltr(), negated.

    `window` is a list of Rec; `center` is rec's index in it.  Only tested
    for records carrying a supplemental attribute (extract_intact_LTR.py
    rule 3 / TE_exciser.py condition (b) both gate on `e.supp`).
    """
    if not rec.supp:
        return False
    # (a) first supplemental attribute contains CUT_BY  -> disrupted
    if "CUT_BY" in rec.supp[0]:
        return True
    # (b) split into pieces by a nested insertion: a neighbour within +/-100
    #     records with identical TSD+strand and an exact/prefix NAME match.
    lo = max(0, center - NEIGHBORHOOD)
    hi = min(len(window), center + NEIGHBORHOOD + 1)
    for j in range(lo, hi):
        if j == center:
            continue
        other = window[j]
        if other.feature_id.startswith("gene"):
            continue
        if other.tsd == rec.tsd and other.strand == rec.strand:
            if rec.name.startswith(other.name) or other.name.startswith(rec.name):
                return True
    return False


# --------------------------------------------------------------------------- #
# Core: stream one BED, classify each record, tally
# --------------------------------------------------------------------------- #
def count_bed(path, verbose=False):
    """Return a dict of counts for one PrinTE BED file.

    Streams the file through a bounded deque so the +/-100 neighbourhood is
    exact while memory stays O(2*NEIGHBORHOOD+1) regardless of file size.
    """
    counts = {
        "solo": 0,            # clean solo-LTR
        "intact": 0,          # intact LTR-RT
        "disrupted_solo": 0,  # _SOLO but split/CUT_BY by a later insertion
        "fragmented_ltr": 0,  # LTR candidate killed by CUT_BY / split
        "frag": 0,            # _FRAG (pre-fragmented) LTR or any _FRAG TE
        "gene": 0,
        "non_ltr": 0,         # non-LTR TE (not solo/intact relevant)
        "malformed": 0,       # < 6 columns
        "total": 0,
    }
    solo_by_super = defaultdict(int)
    intact_by_super = defaultdict(int)

    buf = deque()      # Rec sliding window
    cur = 0            # index in buf of the record being classified

    def classify(rec, window, center):
        counts["total"] += 1
        fid = rec.feature_id
        if fid.startswith("gene"):
            counts["gene"] += 1
            return
        if "_FRAG" in fid:
            counts["frag"] += 1
            return
        te_class, te_super = extract_te_class(fid)
        if te_class != "LTR":
            counts["non_ltr"] += 1
            return
        disrupted = is_disrupted(rec, window, center)
        if "_SOLO" in fid:
            if disrupted:
                counts["disrupted_solo"] += 1
            else:
                counts["solo"] += 1
                solo_by_super[te_super] += 1
        else:
            if disrupted:
                counts["fragmented_ltr"] += 1
            else:
                counts["intact"] += 1
                intact_by_super[te_super] += 1

    src = open(path)
    try:
        line_iter = iter(src)
        primed = False
        while True:
            # Pull one record ahead so the window can extend +NEIGHBORHOOD.
            rec = None
            for raw in line_iter:
                if not raw.strip() or raw.startswith("#"):
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 6:
                    counts["malformed"] += 1
                    if verbose and counts["malformed"] <= 5:
                        sys.stderr.write(
                            f"  warn: skipping malformed line "
                            f"({len(parts)} cols): {raw.strip()[:80]}\n")
                    continue
                rec = Rec(parts[3], parts[4], parts[5])
                break
            if rec is not None:
                buf.append(rec)
            # Prime: fill the look-ahead half before classifying anything.
            if not primed:
                if rec is not None and len(buf) <= NEIGHBORHOOD:
                    continue
                primed = True
            # Classify the current record if one is waiting.
            if cur < len(buf):
                classify(buf[cur], buf, cur)
                cur += 1
            # Bound memory: drop fully-consumed left context (kept >=100 back).
            while cur > NEIGHBORHOOD and len(buf) > 2 * NEIGHBORHOOD + 1:
                buf.popleft()
                cur -= 1
            if rec is None and cur >= len(buf):
                break
    finally:
        src.close()

    counts["solo_by_super"] = dict(solo_by_super)
    counts["intact_by_super"] = dict(intact_by_super)
    return counts


def ratio(solo, intact):
    """solo:intact, with PrinTE's 0/0 := 0.0 convention."""
    if intact > 0:
        return solo / intact
    return 0.0 if solo == 0 else float("inf")


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def fmt_ratio(solo, intact):
    if intact == 0:
        return "0.0000 (no intact LTR-RTs)" if solo == 0 else "inf (intact=0)"
    return f"{solo / intact:.4f}"


def report_file(path, c, verbose):
    gen = extract_generation(path)
    gtag = f"  [generation {gen}]" if gen is not None else ""
    print(f"{path}{gtag}")
    print(f"  SoloLTR (clean) : {c['solo']}")
    print(f"  Intact LTR-RT   : {c['intact']}")
    print(f"  solo:intact     : {fmt_ratio(c['solo'], c['intact'])}")
    if verbose:
        print("  -- detail --")
        print(f"  disrupted solo (excluded)      : {c['disrupted_solo']}")
        print(f"  fragmented LTR, CUT_BY/split    : {c['fragmented_ltr']}")
        print(f"  _FRAG features (excluded)       : {c['frag']}")
        print(f"  genes (skipped)                 : {c['gene']}")
        print(f"  non-LTR TEs (skipped)           : {c['non_ltr']}")
        print(f"  malformed lines (skipped)       : {c['malformed']}")
        print(f"  records examined                : {c['total']}")
        if c["intact_by_super"]:
            print("  intact LTR-RT by superfamily:")
            for k in sorted(c["intact_by_super"]):
                print(f"      {k:<24} {c['intact_by_super'][k]}")
        if c["solo_by_super"]:
            print("  solo-LTR by superfamily:")
            for k in sorted(c["solo_by_super"]):
                print(f"      {k:<24} {c['solo_by_super'][k]}")
        # sanity checks
        if c["solo"] == 0 and c["intact"] == 0:
            sys.stderr.write(
                f"  note: {path} has 0 solo and 0 intact LTR-RTs "
                f"(no LTR-RTs, or none survived as intact).\n")
    if verbose and c["malformed"] > 5:
        sys.stderr.write(
            f"  warn: {c['malformed']} malformed lines skipped in {path}\n")
    print()


def main():
    p = argparse.ArgumentParser(
        description="Count solo-LTRs and intact LTR-RTs in PrinTE BED files "
                    "and report the solo:intact ratio. Definitions follow "
                    "PrinTE (extract_intact_LTR.py / plot_solo_intact.py / "
                    "TE_exciser.py).")
    p.add_argument("bed", nargs="+",
                   help="PrinTE BED file(s) (e.g. gen260000_final.bed, "
                        "burnin.bed). Each is counted independently.")
    p.add_argument("-o", "--tsv", metavar="PATH",
                   help="Also write a machine-readable TSV "
                        "(one row per BED) for cross-genome plotting.")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Per-file breakdown, per-superfamily tables, and "
                        "sanity checks (default output is minimal).")
    args = p.parse_args()

    # Fail fast on bad/foundational input.
    missing = [b for b in args.bed if not os.path.isfile(b)]
    if missing:
        sys.stderr.write("ERROR: BED file(s) not found:\n  " +
                         "\n  ".join(missing) + "\n")
        sys.exit(2)

    rows = []
    tot_solo = tot_intact = 0
    for b in args.bed:
        if args.verbose:
            sys.stderr.write(f"[solo_intact_count] scanning {b} ...\n")
        c = count_bed(b, verbose=args.verbose)
        report_file(b, c, args.verbose)
        rows.append((b, extract_generation(b), c))
        tot_solo += c["solo"]
        tot_intact += c["intact"]

    if len(args.bed) > 1:
        print("TOTAL (all BED files)")
        print(f"  SoloLTR (clean) : {tot_solo}")
        print(f"  Intact LTR-RT   : {tot_intact}")
        print(f"  solo:intact     : {fmt_ratio(tot_solo, tot_intact)}")
        print()

    if args.tsv:
        with open(args.tsv, "w") as fh:
            fh.write("file\tgeneration\tsolo\tintact\tsolo_intact_ratio"
                     "\tdisrupted_solo\tfragmented_ltr\n")
            for b, gen, c in rows:
                fh.write(
                    f"{b}\t{'' if gen is None else gen}\t{c['solo']}"
                    f"\t{c['intact']}\t{ratio(c['solo'], c['intact'])}"
                    f"\t{c['disrupted_solo']}\t{c['fragmented_ltr']}\n")
        if args.verbose:
            sys.stderr.write(f"[solo_intact_count] wrote {args.tsv}\n")


if __name__ == "__main__":
    main()
