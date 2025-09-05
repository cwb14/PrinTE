#!/usr/bin/env python3
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
    Extract TE information from feature_ID.
    Expected format: TE_name#TE_class/TE_superfamily[~junk]
    The function removes junk information starting with '~' from TE_superfamily.
    Returns: te_name, te_class, te_superfamily
    """
    try:
        te_name, rest = feature_id.split("#", 1)
        te_class, te_superfamily = rest.split("/", 1)
        # Remove any junk after a '~' in te_superfamily if present
        te_superfamily = te_superfamily.split("~")[0]
        return te_name, te_class, te_superfamily
    except Exception:
        return None, None, None

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
            if "_FRAG" in record['feature_id']:
                record['category'] = "Fragmented TE"
            elif "_SOLO" in record['feature_id']:
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

def extract_intact_LTRs(records, genome):
    """
    From the list of BED records, extract sequences of intact TEs with TE_class == LTR.
    Returns a list of tuples: (fasta_header, sequence)
    """
    ltr_entries = []
    for record in records:
        if record['category'] != "Intact TE":
            continue
        # Extract TE info from feature_id
        _, te_class, _ = extract_TE_info(record['feature_id'])
        if te_class != "LTR":
            continue

        # Get the genomic coordinates (convert to 0-based integers)
        try:
            chrom = record['chrom']
            start = int(record['start'])
            end = int(record['end'])
        except Exception:
            continue  # skip if conversion fails

        if chrom not in genome:
            continue  # skip if chromosome not found in genome

        seq = genome[chrom][start:end]
        # Use the full NAME from the BED file as FASTA header.
        header = record['name']
        ltr_entries.append((header, str(seq)))
    return ltr_entries

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
        description="Extract intact LTRs from BED files using a genome FASTA. "
                    "Intact LTRs are intact TEs with TE_class = LTR. "
                    "The FASTA header is taken from the BED NAME field."
    )
    parser.add_argument("--bed", nargs="+", required=True, help="Input BED file(s)")
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--out_fasta", required=True, help="Output multi-sequence FASTA file")
    args = parser.parse_args()

    # Load genome FASTA
    genome = load_genome(args.genome)

    # Process each BED file and extract intact LTRs
    all_entries = []
    for bed_file in args.bed:
        records = process_bed_file(bed_file)
        entries = extract_intact_LTRs(records, genome)
        all_entries.extend(entries)

    if not all_entries:
        print("No intact LTR entries found.", flush=True)
    else:
        write_fasta(all_entries, args.out_fasta)
        print(f"Extracted {len(all_entries)} intact LTR entries to {args.out_fasta}", flush=True)

if __name__ == "__main__":
    main()
