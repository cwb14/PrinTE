#!/usr/bin/env python3
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    Process one BED file and return overall counts, TE breakdowns, and the list of records.
    Returns:
      overall_counts: dict mapping category to count
      intact_TE_counts: dict mapping TE_class/TE_superfamily to count (for Intact TE)
      frag_TE_counts: dict mapping TE_class/TE_superfamily to count (for Fragmented TE)
      records: list of record dictionaries (each augmented with classification info)
    """
    overall_counts = {
        "Intact gene": 0,
        "Fragmented gene": 0,
        "Intact TE": 0,
        "SoloLTR": 0,
        "Fragmented TE": 0
    }
    intact_TE_counts = {}
    frag_TE_counts = {}

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
                record['category'] = None  # to be assigned later
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
            elif "_FRAG" in record['feature_id']:
                record['category'] = "Fragmented TE"
            elif record['additional'] and "CUT_BY" in record['additional'][0]:
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

    # Final classification step: Convert remaining "Potential intact TE" to "Intact TE"
    for record in records:
        if record['category'] == "Potential intact TE":
            record['category'] = "Intact TE"

        overall_counts[record['category']] += 1

        # For TE breakdown counts, only consider Intact TE and Fragmented TE with expected format.
        if record['category'] in ["Intact TE", "Fragmented TE"]:
            _, te_class, te_superfamily = extract_TE_info(record['feature_id'])
            if te_class is None or te_superfamily is None:
                continue
            key = f"{te_class}/{te_superfamily}"
            if record['category'] == "Intact TE":
                intact_TE_counts[key] = intact_TE_counts.get(key, 0) + 1
            else:
                frag_TE_counts[key] = frag_TE_counts.get(key, 0) + 1

    return overall_counts, intact_TE_counts, frag_TE_counts, records

def merge_dicts(dict_list, keys):
    """
    Given a list of dictionaries and a set of all keys, return a dictionary mapping each key
    to a list of counts from each dictionary in sorted order.
    """
    merged = {}
    for key in keys:
        merged[key] = [d.get(key, 0) for d in dict_list]
    return merged

def extract_fasta_for_LTR(records, genome_dict):
    """
    Given a list of BED records and a genome (as a dict from SeqIO), extract sequences
    for records that are Intact TE and have TE_class == "LTR". Returns a list of SeqRecord.
    """
    fasta_records = []
    for record in records:
        if record['category'] != "Intact TE":
            continue
        # Extract TE information and check for LTR TE_class.
        _, te_class, _ = extract_TE_info(record['feature_id'])
        if te_class != "LTR":
            continue
        chrom = record['chrom']
        try:
            start = int(record['start'])
            end = int(record['end'])
        except ValueError:
            continue
        if chrom not in genome_dict:
            continue
        # Extract the sequence (BED coordinates are assumed to be 0-based, half-open)
        seq = genome_dict[chrom].seq[start:end]
        # Build a header that includes the feature_id and coordinates
        header = f"{record['feature_id']}|{chrom}:{start}-{end}"
        fasta_record = SeqRecord(Seq(str(seq)), id=header, description="")
        fasta_records.append(fasta_record)
    return fasta_records

def main():
    parser = argparse.ArgumentParser(
        description="Classify BED entries into gene and TE categories, produce summary stats, and optionally extract LTR TE sequences."
    )
    parser.add_argument("--bed", nargs="+", required=True, help="Input BED file(s)")
    parser.add_argument("--genome", help="Genome FASTA file for sequence extraction")
    parser.add_argument("--out_prefix", required=True, help="Output prefix for stats files and extracted FASTA")
    args = parser.parse_args()

    overall_list = []
    intact_TE_list = []
    frag_TE_list = []
    all_records = []
    file_bases = []

    # Process each BED file
    for bed_file in args.bed:
        file_base = os.path.splitext(os.path.basename(bed_file))[0]
        file_bases.append(file_base + "_Count")
        overall, intact_TE, frag_TE, records = process_bed_file(bed_file)
        overall_list.append(overall)
        intact_TE_list.append(intact_TE)
        frag_TE_list.append(frag_TE)
        all_records.extend(records)

    # Categories to output for overall stats
    categories = ["Intact gene", "Fragmented gene", "Intact TE", "SoloLTR", "Fragmented TE"]
    # Merge overall counts
    overall_merged = merge_dicts(overall_list, categories)
    with open(args.out_prefix + "_overall.tsv", "w") as out:
        header = "\t".join(["Category"] + file_bases)
        out.write(header + "\n")
        for cat in categories:
            counts = "\t".join(str(x) for x in overall_merged[cat])
            out.write(f"{cat}\t{counts}\n")

    # Merge intact TE counts (keys are TE_class/TE_superfamily)
    intact_keys = set()
    for d in intact_TE_list:
        intact_keys.update(d.keys())
    intact_merged = merge_dicts(intact_TE_list, intact_keys)
    with open(args.out_prefix + "_intact.tsv", "w") as out:
        header = "\t".join(["TE_class/TE_superfamily"] + file_bases)
        out.write(header + "\n")
        for key in sorted(intact_keys):
            counts = "\t".join(str(x) for x in intact_merged[key])
            out.write(f"{key}\t{counts}\n")

    # Merge fragmented TE counts
    frag_keys = set()
    for d in frag_TE_list:
        frag_keys.update(d.keys())
    frag_merged = merge_dicts(frag_TE_list, frag_keys)
    with open(args.out_prefix + "_frag.tsv", "w") as out:
        header = "\t".join(["TE_class/TE_superfamily"] + file_bases)
        out.write(header + "\n")
        for key in sorted(frag_keys):
            counts = "\t".join(str(x) for x in frag_merged[key])
            out.write(f"{key}\t{counts}\n")

    # If a genome FASTA is provided, extract sequences for intact TEs with TE_class "LTR"
    if args.genome:
        # Load the genome FASTA into a dictionary for quick access.
        genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
        ltr_fasta_records = extract_fasta_for_LTR(all_records, genome_dict)
        fasta_outfile = args.out_prefix + ".fa"
        with open(fasta_outfile, "w") as out_fa:
            SeqIO.write(ltr_fasta_records, out_fa, "fasta")
        print(f"Extracted {len(ltr_fasta_records)} LTR TE sequences to {fasta_outfile}")

if __name__ == "__main__":
    main()
