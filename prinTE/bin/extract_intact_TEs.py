#!/usr/bin/env python3

# Extracts an updated TE library from bed and fasta. 
# Non-LTRs are extracted directly.
# LTR-RTs are modified so that the 3' LTR matches the 5' LTR.
# Script works well, but some minor confusion when seq_divergence.py has a tough time identifying exact LTR lengths (probably doesnt matter).
# To resolve, I can consider pre-processing TE library such that the 5' and 3' LTRs are forced to be exactly identical and exact specified length. 
# Need to incorporate this into the pipeline so that TEs are not drawn from the ancestral library every generation. 

import argparse
import os
from Bio import SeqIO

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
      te_name, te_class, te_superfamily, ltr_len (an integer if available, else None)
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
    Process one BED file.
    Returns a list of records (each is a dict) with additional fields:
      - feature_id: parsed from the NAME column
      - additional: additional attributes (if any)
      - category: classification (e.g. Intact TE, Fragmented TE, etc.)
    """
    records = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            record = parse_line(line)
            if record:
                feature_id, additional = parse_attributes(record['name'])
                record['feature_id'] = feature_id
                record['additional'] = additional
                record['category'] = None  # will be assigned below
                records.append(record)

    # First pass: classify genes and TE candidates (before applying rule 3)
    for record in records:
        if record['feature_id'].startswith("gene"):
            # Gene entries: intact if no additional attributes; fragmented otherwise.
            if len(record['additional']) == 0:
                record['category'] = "Intact gene"
            else:
                record['category'] = "Fragmented gene"
        else:
            # TE candidate
            if "_SOLO" in record['feature_id']:
                record['category'] = "SoloLTR"
            else:
                if record['additional'] and "CUT_BY" in record['additional'][0]:
                    record['category'] = "Fragmented TE"
                else:
                    record['category'] = "Potential intact TE"

    # Rule (3): For potential intact TEs with at least one additional attribute,
    # check within +/- 100 lines for another TE with the same TSD and strand and
    # where one NAME is a prefix of the other.
    for i, record in enumerate(records):
        if record['category'] == "Potential intact TE" and record['additional']:
            for j in range(max(0, i - 100), min(len(records), i + 101)):
                if i == j:
                    continue
                other = records[j]
                if other['feature_id'].startswith("gene"):
                    continue
                if other['tsd'] == record['tsd'] and other['strand'] == record['strand']:
                    if record['name'].startswith(other['name']) or other['name'].startswith(record['name']):
                        record['category'] = "Fragmented TE"
                        break

    # Final classification: Convert remaining "Potential intact TE" to "Intact TE"
    for record in records:
        if record['category'] == "Potential intact TE":
            record['category'] = "Intact TE"

    return records

def load_genome(genome_fasta):
    """
    Load genome FASTA into a dictionary mapping chromosome names to sequences.
    """
    genome = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        genome[record.id] = record.seq
    return genome

def extract_intact_TEs(records, genome):
    """
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
    te_entries = []
    for record in records:
        if record['category'] != "Intact TE":
            continue

        # Get the genomic coordinates (0-based integers)
        try:
            chrom = record['chrom']
            start = int(record['start'])
            end = int(record['end'])
        except Exception:
            continue
        if chrom not in genome:
            continue

        # Extract the sequence from the genome.
        seq = str(genome[chrom][start:end])
        
        # If on minus strand, compute reverse complement.
        if record['strand'] == "-":
            from Bio.Seq import Seq
            seq = str(Seq(seq).reverse_complement())

        # Extract TE info from feature_id.
        te_name, te_class, te_superfamily, ltr_len = extract_TE_info(record['feature_id'])
        
        # Special handling for LTRs.
        if te_class == "LTR" and ltr_len is not None:
            if len(seq) >= 2 * ltr_len:
                five_prime = seq[:ltr_len]
                internal = seq[ltr_len:len(seq)-ltr_len]
                # Replace the 3' LTR with a copy of the 5' LTR.
                seq = five_prime + internal + five_prime
            # If the sequence is too short, leave it unchanged.
        
        # Use the full BED NAME as the FASTA header.
        header = record['name']
        te_entries.append((header, seq))
    return te_entries

def write_fasta(entries, out_file):
    """
    Write a list of (header, sequence) tuples to a FASTA file.
    """
    with open(out_file, "w") as f:
        for header, seq in entries:
            f.write(f">{header}\n")
            # Wrap sequence every 60 characters
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="Extract intact TEs from BED files using a genome FASTA. "
                    "For TEs with TE_class = LTR, special extraction is applied based on LTR length. "
                    "Sequences on the minus strand are reverse complemented. "
                    "The FASTA header is taken from the BED NAME field."
    )
    parser.add_argument("--bed", nargs="+", required=True, help="Input BED file(s)")
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--out_fasta", required=True, help="Output multi-sequence FASTA file")
    args = parser.parse_args()

    # Load genome FASTA
    genome = load_genome(args.genome)

    # Process each BED file and extract intact TEs.
    all_entries = []
    for bed_file in args.bed:
        records = process_bed_file(bed_file)
        entries = extract_intact_TEs(records, genome)
        all_entries.extend(entries)

    if not all_entries:
        print("No intact TE entries found.", flush=True)
    else:
        write_fasta(all_entries, args.out_fasta)
        print(f"Extracted {len(all_entries)} intact TE entries to {args.out_fasta}", flush=True)

if __name__ == "__main__":
    main()
