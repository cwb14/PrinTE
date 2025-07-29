#!/usr/bin/env python3
import sys
import re
import argparse

### LTR-RTs are fragmented into "LTR" and "Internal" in many DB downloads. This script is to stitch LTR pieces into intact LTR-RTs.

### Download Repeat library.
### Merge repeats from all species
# cat RepBase*fasta/*ref > all.fa

### There are a few duplicates (same header & same sequence). Remove duplicates and call it 'all.remove_dups.fa'.
### From header, remove spaces, tabs, hyphens, slash, and underscore.
# cat all.remove_dups.fa | awk '/^>/{$0=gensub(/[ _\t\/-]/, "", "g"); print; next} {print}' > all.remove_dups.squeezed_headers.fa
### I think its a good move since those seperating charcters can sometimes be inconsistent between LTR-RT pieces. Better to remove than build regex to parse variability.

### This script produces 'all.remove_dups.squeezed_headers.stitched.ltr.fa', where LTR-RT fragments (internal and LTR) are stitched to intact LTR-RTs.
### Non-LTRs and unpaired LTR fragments pass through the script unmodified.
### Extract LTR-RTs and remaining fragments.
# awk '/^>/ {p = tolower($0) ~ /copia|gypsy|ltr|ty3|ty1/} p' all.remove_dups.squeezed_headers.stitched.fa > all.remove_dups.squeezed_headers.stitched.ltr.fa

### Make multi-line fasta single-line | add "_LTR" suffix to header.
# cat all.remove_dups.squeezed_headers.stitched.ltr.fa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk '/^>/ {$0=$0"_LTR"} {print}' > all.remove_dups.squeezed_headers.stitched.ltr.clean.fa

### There are still many LTR-RT fragments in the db that were not paired using the regex below.
### Use 'seq_divergence.py' to find them.
# python /home/chris/data/TEGE/benchmarking/mutations/TESS/PrinTE/bin/seq_divergence.py -i all.remove_dups.squeezed_headers.stitched.ltr.clean.fa -o all.remove_dups.squeezed_headers.stitched.ltr.fa.clean.tsv -t 100 --min_align 50 --max_off 20 --miu 3e-8 --blast_outfmt '6 qseqid sseqid sstart send slen qstart qend qlen length nident btop'
### The output TSV contains IDs for the varifiably intact ones. I can use these IDs to purge fragments (IDs absent from TSV).


def parse_fasta(stream):
    """Yield (header, sequence) pairs (header without the leading '>')."""
    header = None
    seq_parts = []
    for line in stream:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_parts)
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)

def main():
    parser = argparse.ArgumentParser(
        description="Merge LTR + internal into intact LTR-RTs using header tags"
    )
    parser.add_argument("infile", help="Input multi-FASTA")
    parser.add_argument("-o", "--output", help="Write merged FASTA here (default: stdout)")
    args = parser.parse_args()

    # Open output
    outfh = open(args.output, "w") if args.output else sys.stdout

    # Read all records and bucket by key = header with LTR/I removed
    buckets = {}
    for hdr, seq in parse_fasta(open(args.infile)):
        m = re.match(r"^(.*?)(LTR|I)(.*)$", hdr)
        if m:
            pre, tag, post = m.groups()
            key = pre + post
            buckets.setdefault(key, {})[tag] = (hdr, seq)
        else:
            # no LTR/I tag → treat as “other”
            buckets.setdefault(hdr, {})["OTHER"] = (hdr, seq)

    # Process each bucket
    for key, parts in buckets.items():
        if "LTR" in parts and "I" in parts:
            ltr_hdr,   ltr_seq   = parts["LTR"]
            int_hdr,   int_seq   = parts["I"]
            intact_hdr = key
            intact_seq = ltr_seq + int_seq + ltr_seq

            # Write FASTA
            print(f">{intact_hdr}", file=outfh)
            print(intact_seq, file=outfh)

            # Log merge
            print(
                f"LTR {ltr_hdr} merge with internal {int_hdr} to form {intact_hdr}",
                file=sys.stderr
            )

        else:
            # Unpaired (either only LTR, only I, or OTHER) → pass through
            for hdr, seq in parts.values():
                print(f">{hdr}", file=outfh)
                print(seq, file=outfh)

    if args.output:
        outfh.close()

if __name__ == "__main__":
    main()
