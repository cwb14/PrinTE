#!/usr/bin/env python3


# Extracts an updated TE library from bed and fasta. 
# Non-LTRs are extracted directly.
# LTR-RTs are modified so that the 3' LTR matches the 5' LTR.
# Script works well, but some minor confusion when seq_divergence.py has a tough time identifying exact LTR lengths (probably doesnt matter).
# To resolve, I can consider pre-processing TE library such that the 5' and 3' LTRs are forced to be exactly identical and exact specified length. # I added library mode to do this.
# Our length based deletions leads to pereferential loss of long elements. I added '--weight_by' to offset this. The fasta provided here is the base TE library. 

"""
extract_te_with_weighted_sampling.py

Extracts TE sequences from a genome (--bed + --genome) or processes a TE library (--lib),
with optional length-weighted resampling to match a guide distribution.

Usage:
  Genome mode:
    extract_te_with_weighted_sampling.py --bed file1.bed file2.bed --genome genome.fa \
        --out_fasta out.fa [--weight_by guide.fa] [--plot_kde_comparison]
  Library mode:
    extract_te_with_weighted_sampling.py --lib te_library.fa --out_fasta out.fa

Options:
  -h, --help            Show this help message and exit
  --h, -help            Show this help message and exit (aliases)
  --bed BED [BED ...]   One or more BED files specifying TE coordinates (genome mode)
  --genome FASTA        Reference genome FASTA (required with --bed)
  --lib FASTA           TE library FASTA to correct LTR ends (library mode)
  --out_fasta FILE      Output FASTA path (required)
  --weight_by FASTA     Guide FASTA whose length distribution is targeted
  --duplication_mode    Retain all original sequences and duplicate additional copies according to importance weights (only with --weight_by)
  --plot_kde_comparison Save KDE comparison plot as <out_basename>_kde_comparison.pdf

Functions:
  parse_line           Parse a BED line into fields
  parse_attributes     Split NAME field into feature_id and attributes
  extract_TE_info      Decode feature_id into (name, class, superfamily, ltr_len)
  process_bed_file     Classify BED records into intact/fragmented TEs & genes
  load_genome          Load genome FASTA into a dict
  extract_intact_TEs   Extract and LTR-fix intact TE sequences
  process_library_fasta  Fix LTR ends in an existing library FASTA
  write_fasta          Write sequences to FASTA with line-wrapping
  weighted_resample    KDE-based importance sampling & optional plotting
  main                 Entry point: parse args and dispatch modes
"""
import argparse
import os
import sys
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq

def parse_line(line):
    """
    Parse a BED line into its 6 mandatory columns.
    Returns a dict or None if invalid.
    """
    parts = line.strip().split("\t")
    if len(parts) < 6:
        return None
    chrom, start, end, name, tsd, strand = parts[:6]
    return {'chrom': chrom, 'start': start, 'end': end, 'name': name,
            'tsd': tsd, 'strand': strand}


def parse_attributes(name):
    """
    Split the NAME field into feature_id before ';' and list of additional attrs.
    """
    if ";" in name:
        parts = name.split(";")
        return parts[0], parts[1:]
    return name, []


def extract_TE_info(feature_id):
    """
    From a feature_id like "TE#LTR/Super~LTRlen:100", extract:
      - te_name (before '#')
      - te_class (e.g., 'LTR')
      - te_superfamily
      - ltr_len (int) or None
    Returns (te_name, te_class, te_superfamily, ltr_len).
    """
    try:
        te_name, rest = feature_id.split("#", 1)
        te_class, remainder = rest.split("/", 1)
        ltr_len = None
        if te_class == "LTR" and "~" in remainder:
            te_superfamily, ltr_info = remainder.split("~", 1)
            if ltr_info.startswith("LTRlen:"):
                try:
                    ltr_len = int(ltr_info.split(':', 1)[1])
                except ValueError:
                    ltr_len = None
        else:
            te_superfamily = remainder.split("~")[0]
        return te_name, te_class, te_superfamily, ltr_len
    except Exception:
        return None, None, None, None


def process_bed_file(bed_file):
    """
    Load and classify records from a BED file (genome mode).
    Returns a list of records (each is a dict) with additional fields:
      - feature_id: parsed from the NAME column
      - additional: additional attributes (if any)
      - category: classification (e.g. Intact TE, Fragmented TE, etc.)
    """
    records = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            rec = parse_line(line)
            if not rec:
                continue
            fid, attrs = parse_attributes(rec['name'])
            rec.update({'feature_id': fid, 'additional': attrs,
                        'category': None})
            records.append(rec)

    # Initial classification
    for rec in records:
        fid, add = rec['feature_id'], rec['additional']
        if fid.startswith("gene"):
            rec['category'] = "Intact gene" if not add else "Fragmented gene"
        else:
            if "_FRAG" in fid:
                rec['category'] = "Fragmented TE"
            elif "_SOLO" in fid:
                rec['category'] = "SoloLTR"
            elif add and "CUT_BY" in add[0]:
                rec['category'] = "Fragmented TE"
            else:
                rec['category'] = "Potential intact TE"

    # Pair potential intact TEs to reclassify as fragmented if needed
    for i, rec in enumerate(records):
        if rec['category'] != "Potential intact TE":
            continue
        for j in range(max(0, i - 100), min(len(records), i + 101)):
            if i == j:
                continue
            other = records[j]
            if other['feature_id'].startswith("gene"):
                continue
            if (other['tsd'] == rec['tsd'] and
                other['strand'] == rec['strand'] and
                (rec['name'].startswith(other['name']) or
                 other['name'].startswith(rec['name']))):
                rec['category'] = "Fragmented TE"
                break

    # Finalize intact TEs
    for rec in records:
        if rec['category'] == "Potential intact TE":
            rec['category'] = "Intact TE"

    return records


def load_genome(genome_fasta):
    """Load genome FASTA into dict {seqid: Seq}."""
    genome = {}
    for rec in SeqIO.parse(genome_fasta, "fasta"):
        genome[rec.id] = rec.seq
    return genome


def extract_intact_TEs(records, genome):
    """
    From genome‐mode records, extract intact TE sequences (with LTR fix
    From the list of BED records, extract sequences of intact TEs (all types) 
    from the genome FASTA.
    
    For each intact TE:
      - The genomic region is extracted based on the BED coordinates.
      - If the BED indicates the minus strand, the sequence is reverse complemented.
    
    Special handling for TEs where TE_class is LTR:
      - The feature_id contains an LTR length field (e.g., ~LTRlen:4).
      - The extracted sequence is assumed to consist of a 5' LTR, an internal region,
        and a 3' LTR. The script will replace the 3' LTR with a copy of the 5' LTR.
        For example, if the sequence is "GCTAGCGGCACG" with LTRlen 4, then the 5' LTR is "GCTA"
        and the internal region is "GCGG", and the 3' LTR ("CACG") is replaced by "GCTA"
        yielding a final sequence of "GCTAGCGGGCTA".
    
    Returns a list of tuples: (fasta_header, sequence)
    """
    out = []
    for rec in records:
        if rec['category'] != "Intact TE":
            continue
        seq = str(genome.get(rec['chrom'], "")[int(rec['start']):int(rec['end'])])
        if not seq:
            continue
        if rec['strand'] == "-":
            seq = str(Seq(seq).reverse_complement())

        _, te_class, _, ltr_len = extract_TE_info(rec['feature_id'])
        if te_class == "LTR" and ltr_len and len(seq) >= 2 * ltr_len:
            five = seq[:ltr_len]
            internal = seq[ltr_len:-ltr_len]
            seq = five + internal + five

#       out.append((rec['name'], seq))
        out.append((rec['feature_id'], seq))
    return out


def process_library_fasta(lib_fasta):
    """
    Library mode: read an existing TE library FASTA,
    fix LTR entries by copying 5' LTR to 3' end.
    """
    out = []
#   for rec in SeqIO.parse(lib_fasta, "fasta"):
#       header = rec.description
    for rec in SeqIO.parse(lib_fasta, "fasta"):
        # pull only the first token of the header, then strip any ";"-attrs
        name = rec.description.split()[0]
        feature_id, _ = parse_attributes(name)
        seq = str(rec.seq)
#       _, te_class, _, ltr_len = extract_TE_info(header.split()[0])
        _, te_class, _, ltr_len = extract_TE_info(feature_id)
        if te_class == "LTR" and ltr_len and len(seq) >= 2 * ltr_len:
            five = seq[:ltr_len]
            internal = seq[ltr_len:-ltr_len]
            seq = five + internal + five
#       out.append((header, seq))
        out.append((feature_id, seq))
    return out

def write_fasta(entries, out_file):
    """Write (header, seq) pairs to FASTA, wrapping at 60 bp."""
    with open(out_file, "w") as f:
        for header, seq in entries:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def weighted_resample(entries, guide_fasta, out_base, *,
                      plot=False, seed=42,
                      duplication_mode=False):
    """
    KDE-based importance resampling that minimises the number of unique
    records that are excluded when --duplication_mode is *not* requested.

    * If --duplication_mode is given we keep the original behaviour
      (retain all originals + add extra copies).
    * Otherwise we:
        1. Convert KDE weights -> expected copy counts.
        2. Give each record floor(expected) copies.
        3. Hand out the still-unassigned copies by a single weighted
           draw **without replacement** using the fractional parts.
    """
    import math
    rng = np.random.default_rng(seed)

    # ----------  KDE densities & importance weights  ----------
    te_lengths   = np.fromiter((len(seq) for _, seq in entries), dtype=float)
    guide_lengths = np.fromiter((len(r.seq) for r in SeqIO.parse(guide_fasta, "fasta")),
                                dtype=float)

    kde_te    = gaussian_kde(te_lengths,    bw_method="scott")
    kde_guide = gaussian_kde(guide_lengths, bw_method="scott")

    dens_te     = kde_te(te_lengths)
    dens_guide  = kde_guide(te_lengths)
    weights     = dens_guide / (dens_te + 1e-8)          # importance weights

    # ------------------------------------------------------------------
    #  --duplication_mode : keep the original behaviour (all originals +
    #                       extra copies proportional to weight)
    # ------------------------------------------------------------------
    if duplication_mode:
        mean_w      = weights.mean()
        dup_counts  = np.maximum(1, np.rint(weights / mean_w).astype(int))
        resampled   = [entry
                       for idx, cnt in enumerate(dup_counts)
                       for entry in [entries[idx]] * cnt]

        print(f"Original sequences : {len(entries)}")
        print(f"Duplicated copies  : {resampled.__len__() - len(entries)}")
        print(f"Total sequences    : {len(resampled)}", flush=True)

    # ------------------------------------------------------------------
    #  default : minimise exclusions, allow duplicates only when floor(ci) ≥ 2
    # ------------------------------------------------------------------
    else:
        n = len(entries)
        probs        = weights / weights.sum()
        exp_counts   = probs * n                   # expected copies (may be <1 or >1)
        base_counts  = np.floor(exp_counts).astype(int)
        remaining    = n - base_counts.sum()      # slots that still need to be filled
        if remaining:                             # distribute the residual slots
            residual = exp_counts - base_counts
            # If all residuals are 0 (rare edge-case), fall back to uniform draw
            residual_prob = residual / residual.sum() if residual.sum() else \
                            np.full_like(residual, 1 / len(residual))
            extra_idx = rng.choice(len(entries), size=remaining,
                                   replace=False, p=residual_prob)
            base_counts[extra_idx] += 1

        # Build the final sample list
        resampled = [entry
                     for idx, cnt in enumerate(base_counts)
                     for entry in [entries[idx]] * cnt]

        lost      = (base_counts == 0).sum()
        dups      = (base_counts > 1).sum()
        print(f"Total sequences   : {n}")
        print(f"Unique lost       : {lost}")
        print(f"Entries duplicated: {dups}", flush=True)

    # ----------  Optional KDE comparison plot  ----------
    if plot:
        sampled_lengths = np.fromiter((len(seq) for _, seq in resampled), dtype=float)
        x = np.linspace(min(te_lengths.min(), guide_lengths.min()),
                        max(te_lengths.max(), guide_lengths.max()), 1000)

        plt.figure()
        plt.plot(x, kde_te(x),                       label="Original TE")
        plt.plot(x, kde_guide(x),                    label="Guide")
        plt.plot(x, gaussian_kde(sampled_lengths)(x), label="Resampled")
        plt.xlabel("Sequence length")
        plt.ylabel("Density")
        plt.legend()
        pdf = f"{out_base}_kde_comparison.pdf"
        plt.savefig(pdf)
        plt.close()
        print(f"KDE plot saved to {pdf}", flush=True)

    return resampled

def main():
    parser = argparse.ArgumentParser(
        description="Extract TE sequences (genome or library mode) with optional length-weighted sampling",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--h', '-help', action='help', help='Show this help message and exit')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--lib', help='Library-mode input FASTA')
    group.add_argument('--bed', nargs='+', help='Genome-mode input BED file(s)')

    parser.add_argument('--genome', help='Genome FASTA (required with --bed)')
    parser.add_argument('--out_fasta', required=True, help='Output FASTA path')
    parser.add_argument('--weight_by', help='Guide FASTA for length-weighted sampling')
    parser.add_argument('--duplication_mode', action='store_true',
                        help='Retain all original sequences and duplicate additional copies according to importance weights (only with --weight_by)')
    parser.add_argument('--plot_kde_comparison', action='store_true',
                        help='Save KDE comparison plot to PDF')

    args = parser.parse_args()

    if args.lib:
        entries = process_library_fasta(args.lib)
    else:
        if not args.genome:
            parser.error('--genome is required when using --bed')
        genome = load_genome(args.genome)
        recs = []
        for b in args.bed:
            recs.extend(process_bed_file(b))
        entries = extract_intact_TEs(recs, genome)

        if not entries:
            print('No intact TE entries found.', file=sys.stderr)
            sys.exit(1)

        if args.weight_by:
            base = os.path.splitext(os.path.basename(args.out_fasta))[0]
            entries = weighted_resample(
                entries,
                args.weight_by,
                out_base=base,
                plot=args.plot_kde_comparison,
                duplication_mode=args.duplication_mode
            )

    write_fasta(entries, args.out_fasta)
    print(f"Processed {len(entries)} entries → {args.out_fasta}", flush=True)

if __name__ == '__main__':
    main()
