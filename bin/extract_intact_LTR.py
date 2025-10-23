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

    # Rule (3): for potential intact TEs with at least one additional attribute, look nearby
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
    Returns a list of tuples: (original_header, sequence_string)
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
        # Save the original BED NAME as the "original_header"
        original_header = record['name']
        ltr_entries.append((original_header, str(seq)))
    return ltr_entries

def wrap_and_write_sequence(handle, header, seq, width=60):
    handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i+width] + "\n")

def write_fasta_numeric(entries, out_fasta, out_key, start_index=1):
    """
    Write sequences to FASTA with numeric headers and produce a key TSV.
    entries: list of (original_header, sequence_string)
    out_fasta: path to FASTA with numeric headers
    out_key: path to TSV mapping numeric_id -> original_header
    """
    count = 0
    with open(out_fasta, "w") as f_fa, open(out_key, "w") as f_key:
        for idx, (orig_header, seq) in enumerate(entries, start=start_index):
            wrap_and_write_sequence(f_fa, str(idx), seq, width=60)
            f_key.write(f"{idx}\t{orig_header}\n")
            count += 1
    return count

def main():
    parser = argparse.ArgumentParser(
        description="Extract intact LTRs from BED files using a genome FASTA. "
                    "Outputs a multi-sequence FASTA with numeric headers (>1, >2, ...) "
                    "and a key mapping numeric IDs to the original BED NAME."
    )
    parser.add_argument("--bed", nargs="+", required=True, help="Input BED file(s)")
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--out_fasta", required=True, help="Output multi-sequence FASTA file (numeric headers)")
    parser.add_argument("--out_key", default=None,
                        help="Output TSV mapping numeric_id to original header (default: <out_fasta>.key.tsv)")
    parser.add_argument("--start_index", type=int, default=1,
                        help="Starting numeric index for FASTA headers (default: 1)")
    args = parser.parse_args()

    out_key = args.out_key if args.out_key else f"{args.out_fasta}.key.tsv"

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
        return

    # Write numeric FASTA and key
    n = write_fasta_numeric(all_entries, args.out_fasta, out_key, start_index=args.start_index)
    print(f"Extracted {n} intact LTR entries to {args.out_fasta}", flush=True)
    print(f"Key (numeric_id -> original header) written to {out_key}", flush=True)

if __name__ == "__main__":
    main()
