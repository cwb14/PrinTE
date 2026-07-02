#!/usr/bin/env python3
"""
Filter LTR_RT_lib.fa to detector-fair intact LTR-RTs.

The custom library harvested by de-novo annotation carries a large `unknown`
superfamily tail of small / short-LTR / domain-less copies that no LTR detector
(minlenltr=100, +0.9 reciprocal overlap) can re-find. Inserting them as intact
"truth" deflates recall. This drops the elements that fail the detector's
geometry envelope. Named clades (Copia/Gypsy/...) are kept regardless of LTR
length by default, since a recognizable protein-coding internal region lets the
detector anchor them even with a short LTR.

Header format expected:
  >acc:start-end#LTR/superfamily/clade~LTRlen:NNN

Usage:
  python filter_ltr_lib.py -i LTR_RT_lib.fa -o LTR_RT_lib.filt.fa \
      --min-ltrlen 100 --max-ltrlen 7000 --min-total 1000
  # add --require-named to also drop ALL unknown-superfamily elements
  # add --named-bypass-ltrlen/--no-named-bypass-ltrlen to toggle the clade rule
"""
import argparse, re, sys

HDR = re.compile(r':(\d+)-(\d+)#([^/]+)/([^/]+)/([^~]+)~LTRlen:(\d+)')


def parse_header(h):
    m = HDR.search(h)
    if not m:
        return None
    total = abs(int(m.group(2)) - int(m.group(1)))
    return dict(total=total, superfamily=m.group(4), clade=m.group(5),
               ltrlen=int(m.group(6)), internal=total - 2*int(m.group(6)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--min-ltrlen", type=int, default=100,
                    help="drop if LTR shorter than this (detector minlenltr; default 100)")
    ap.add_argument("--max-ltrlen", type=int, default=7000,
                    help="drop if LTR longer than this (detector maxlenltr; default 7000)")
    ap.add_argument("--min-internal", type=int, default=100,
                    help="drop if internal region shorter than this (mindistltr; default 100)")
    ap.add_argument("--min-total", type=int, default=1000,
                    help="drop if total element shorter than this (default 1000)")
    ap.add_argument("--require-named", action="store_true",
                    help="also drop every unknown-superfamily element")
    ap.add_argument("--named-bypass-ltrlen", dest="named_bypass", action="store_true",
                    default=True, help="named clades keep short LTRs (default: on)")
    ap.add_argument("--no-named-bypass-ltrlen", dest="named_bypass", action="store_false")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    kept = dropped = 0
    reasons = {"short_ltr": 0, "long_ltr": 0, "short_internal": 0,
               "small_total": 0, "unknown": 0, "unparsed": 0}

    def keep(rec):
        named = rec["superfamily"].lower() not in ("unknown", "na")
        if args.require_named and not named:
            reasons["unknown"] += 1; return False
        if rec["ltrlen"] > args.max_ltrlen:
            reasons["long_ltr"] += 1; return False
        if rec["internal"] < args.min_internal:
            reasons["short_internal"] += 1; return False
        if rec["total"] < args.min_total:
            reasons["small_total"] += 1; return False
        if rec["ltrlen"] < args.min_ltrlen and not (named and args.named_bypass):
            reasons["short_ltr"] += 1; return False
        return True

    with open(args.input) as fi, open(args.output, "w") as fo:
        write = False
        for ln in fi:
            if ln.startswith(">"):
                rec = parse_header(ln)
                if rec is None:
                    write = False; dropped += 1; reasons["unparsed"] += 1
                elif keep(rec):
                    write = True; kept += 1
                else:
                    write = False; dropped += 1
            if write:
                fo.write(ln)

    tot = kept + dropped
    print(f"kept {kept}/{tot} ({100*kept/tot:.1f}%)   dropped {dropped}", file=sys.stderr)
    for k, v in reasons.items():
        if v:
            print(f"  dropped[{k}] = {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
