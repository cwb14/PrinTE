#!/usr/bin/env python3

# Extracts an updated TE library from bed and fasta. 
# Non-LTRs are extracted directly.
# LTR-RTs are modified so that the 3' LTR matches the 5' LTR.
# Script works well, but some minor confusion when seq_divergence.py has a tough time identifying exact LTR lengths (probably doesnt matter).
# To resolve, I can consider pre-processing TE library such that the 5' and 3' LTRs are forced to be exactly identical and exact specified length. # I added library mode to do this.
# Need to incorporate this into the pipeline so that TEs are not drawn from the ancestral library every generation. 


import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq

def parse_line(line):
    """Parse a BED line into its columns."""
    parts = line.strip().split("\t")
    if len(parts) < 6:
        return None
    chrom, start, end, name, tsd, strand = parts[:6]
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'name': name,
        'tsd': tsd,
        'strand': strand
    }

def parse_attributes(name):
    """
    Split the NAME into feature_ID and additional attributes.
    If no semicolon is present, feature_ID is the entire name.
    """
    if ";" in name:
        parts = name.split(";")
        feature_id = parts[0]
        additional = parts[1:]
    else:
        feature_id = name
        additional = []
    return feature_id, additional

def extract_TE_info(feature_id):
    """
    Extract TE information from feature_id.
    Expected format for LTRs: TE_name#LTR/TE_superfamily~LTRlen:X 
    For non-LTRs, any content after "~" is ignored.
    Returns:
      te_name, te_class, te_superfamily, ltr_len (int or None)
    """
    try:
        te_name, rest = feature_id.split("#", 1)
        te_class, remainder = rest.split("/", 1)
        ltr_len = None
        if te_class == "LTR" and "~" in remainder:
            te_superfamily, ltr_info = remainder.split("~", 1)
            if ltr_info.startswith("LTRlen:"):
                try:
                    ltr_len = int(ltr_info[len("LTRlen:"):])
                except ValueError:
                    ltr_len = None
        else:
            te_superfamily = remainder.split("~")[0] if "~" in remainder else remainder
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
            feature_id, additional = parse_attributes(rec['name'])
            rec.update({
                'feature_id': feature_id,
                'additional': additional,
                'category': None
            })
            records.append(rec)

    # 1st pass: gene vs TE candidates
    for rec in records:
        fid, add = rec['feature_id'], rec['additional']
        if fid.startswith("gene"):
            rec['category'] = "Intact gene" if not add else "Fragmented gene"
        else:
            if "_SOLO" in fid:
                rec['category'] = "SoloLTR"
            elif add and "CUT_BY" in add[0]:
                rec['category'] = "Fragmented TE"
            else:
                rec['category'] = "Potential intact TE"

    # rule (3): pair up potential intact TEs by TSD & strand & name‐prefix
    for i, rec in enumerate(records):
        if rec['category']=="Potential intact TE" and rec['additional']:
            for j in range(max(0, i-100), min(len(records), i+101)):
                if i==j: continue
                other = records[j]
                if other['feature_id'].startswith("gene"): continue
                if (other['tsd']==rec['tsd'] and
                    other['strand']==rec['strand'] and
                    (rec['name'].startswith(other['name']) or
                     other['name'].startswith(rec['name']))):
                    rec['category']="Fragmented TE"
                    break

    # finalize intact TEs
    for rec in records:
        if rec['category']=="Potential intact TE":
            rec['category']="Intact TE"

    return records

def load_genome(genome_fasta):
    """Load genome FASTA into a dict of SeqRecords."""
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
        if rec['category']!="Intact TE":
            continue
        chrom = rec['chrom']
        start, end = int(rec['start']), int(rec['end'])
        if chrom not in genome:
            continue
        seq = str(genome[chrom][start:end])
        if rec['strand']=="-":
            seq = str(Seq(seq).reverse_complement())

        _, te_class, _, ltr_len = extract_TE_info(rec['feature_id'])
        if te_class=="LTR" and ltr_len and len(seq)>=2*ltr_len:
            five = seq[:ltr_len]
            internal = seq[ltr_len:len(seq)-ltr_len]
            seq = five + internal + five

        header = rec['name']
        out.append((header, seq))
    return out

def process_library_fasta(lib_fasta):
    """
    Library mode: read an existing TE library FASTA,
    fix LTR entries by copying 5' LTR to 3' end.
    """
    out = []
    for rec in SeqIO.parse(lib_fasta, "fasta"):
        header = rec.description    # full header line, no leading '>'
        seq = str(rec.seq)
        feature_id, _ = parse_attributes(header)
        _, te_class, _, ltr_len = extract_TE_info(feature_id)
        if te_class=="LTR" and ltr_len and len(seq)>=2*ltr_len:
            five = seq[:ltr_len]
            internal = seq[ltr_len:len(seq)-ltr_len]
            seq = five + internal + five
        out.append((header, seq))
    return out

def write_fasta(entries, out_file):
    """Write (header, seq) pairs to FASTA, wrapping at 60 bp."""
    with open(out_file, "w") as f:
        for header, seq in entries:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="Two modes: genome mode (--bed + --genome) extracts and fixes LTRs from genome+BED; "
                    "library mode (--lib) fixes LTRs in a supplied TE library FASTA."
                    "For TEs with TE_class = LTR, special extraction is applied based on LTR length. "
                    "Genome mode extracts intact TEs from BED files using a genome FASTA. "
                    "Sequences on the minus strand are reverse complemented. "
                    "The FASTA header is taken from the BED NAME field."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--lib", help="Library‐mode: input TE library FASTA")
    mode.add_argument("--bed", nargs="+", help="Genome‐mode: input BED file(s)")
    parser.add_argument("--genome", help="Genome FASTA (required with --bed)")
    parser.add_argument("--out_fasta", required=True, help="Output FASTA")
    args = parser.parse_args()

    if args.lib:
        entries = process_library_fasta(args.lib)
        write_fasta(entries, args.out_fasta)
        print(f"Processed {len(entries)} library entries → {args.out_fasta}", flush=True)

    else:
        # genome mode
        if not args.genome:
            parser.error("--genome is required when using --bed")
        genome = load_genome(args.genome)
        all_entries = []
        for b in args.bed:
            recs = process_bed_file(b)
            all_entries.extend(extract_intact_TEs(recs, genome))
        if not all_entries:
            print("No intact TE entries found.", flush=True)
        else:
            write_fasta(all_entries, args.out_fasta)
            print(f"Extracted {len(all_entries)} intact TE entries → {args.out_fasta}", flush=True)

if __name__ == "__main__":
    main()
