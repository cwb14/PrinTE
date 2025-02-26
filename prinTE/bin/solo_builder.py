#!/usr/bin/env python3
"""
Script to convert a fraction of intact LTR-RTs to soloLTRs,
updating both the genome FASTA (by deleting internal regions) and
the corresponding BED file (adjusting coordinates accordingly).

Usage:
  python convert_soloLTR.py --genome input.fasta --bed input.bed --output out_prefix --soloLTR_freq 5 --seed 42

The BED file is assumed to have 6 columns:
  chrom, start, end, NAME, TSD, strand
where:
  - NAME is a semicolon-delimited string whose first field is the feature_ID.
  - TSD is a short sequence (e.g. "caaag").
  - strand is the typical plus/minus orientation.
An intact LTR-RT is one whose feature_ID (first attribute) contains “LTRlen” 
and does not already contain “_SOLO”. In addition, the following rules are applied:
  1) Exclude LTR-RTs where the first supplemental field (the field after the first semicolon)
     contains "CUT_BY".
  2) If supplemental info is present and a TSD exists (from column 5), then if within ±100 lines 
     a record is found that has supplemental info (i.e. its NAME contains a semicolon) and whose NAME 
     is a partial match (one name is a substring of the other) and that has identical TSD and strand, 
     then the record is considered non-intact.
For a conversion, the script parses the LTRlen (e.g. “LTRlen:300”) from the feature_ID, and then resets 
the feature’s coordinates so that only the 5′ LTR (of length LTRlen) remains. The genome FASTA is modified 
by “deleting” the region from start+LTRlen to end, and all coordinates in the BED file are updated accordingly.
"""

import argparse
import random
import sys
import re
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a fraction of intact LTR-RTs to soloLTRs, updating both FASTA and BED files."
    )
    parser.add_argument("--genome", required=True, help="Input genome FASTA file.")
    parser.add_argument("--soloLTR_freq", type=float, required=True,
                        help="Percentage (e.g. 5 for 5%%) of intact LTR-RTs to convert to soloLTR.")
    parser.add_argument("--bed", required=True, help="Existing BED file with TE/gene coordinates.")
    parser.add_argument("--output", required=True, help="Output prefix (for .bed and .fasta).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    return parser.parse_args()

def read_fasta(fasta_file):
    """
    Reads a FASTA file and returns a dictionary:
      {chrom: (description, sequence_string)}
    """
    fasta = {}
    header = None
    seq_lines = []
    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    fasta[header.split()[0]] = (" ".join(header.split()[1:]) if len(header.split()) > 1 else "", "".join(seq_lines))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            fasta[header.split()[0]] = (" ".join(header.split()[1:]) if len(header.split()) > 1 else "", "".join(seq_lines))
    return fasta

def write_fasta(fasta_dict, output_file):
    with open(output_file, "w") as outfh:
        for chrom in fasta_dict:
            desc, seq = fasta_dict[chrom]
            header = f">{chrom}"
            if desc:
                header += " " + desc
            outfh.write(header + "\n")
            # Wrap sequence to 60 characters per line.
            for i in range(0, len(seq), 60):
                outfh.write(seq[i:i+60] + "\n")

def read_bed(bed_file):
    """
    Reads a BED file (assumed 6 columns) and returns a list of dictionaries:
      { chrom, start, end, name, tsd, strand, orig_line }
    Here we assume columns:
      1: chrom, 2: start, 3: end, 4: NAME, 5: TSD, 6: strand.
    """
    records = []
    with open(bed_file, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                sys.exit("Error: BED file must have at least 6 columns.")
            rec = {
                "chrom": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "name": parts[3],
                "tsd": parts[4],
                "strand": parts[5],
                "line": line  # original line for reference if needed.
            }
            records.append(rec)
    return records

def parse_feature_ID(name_field):
    """
    Given the NAME field (a semicolon-delimited list of attributes), return:
      feature_id, other_attributes_list
    where feature_id is the first attribute.
    """
    attrs = name_field.split(";")
    feature_id = attrs[0]
    other_attrs = attrs[1:] if len(attrs) > 1 else []
    return feature_id, other_attrs

def update_feature_ID(feature_id):
    """
    Append '_SOLO' to feature_id.
    """
    return feature_id + "_SOLO"

def parse_LTRlen(feature_id):
    """
    Given a feature_id string that contains "LTRlen:XXX", extract and return the integer XXX.
    Returns None if not found.
    """
    m = re.search(r"LTRlen:(\d+)", feature_id)
    if m:
        return int(m.group(1))
    return None

def main():
    args = parse_args()
    random.seed(args.seed)

    # Read the genome FASTA.
    genome_dict = read_fasta(args.genome)  # {chrom: (desc, sequence)}

    # Read the BED file records.
    bed_records = read_bed(args.bed)

    # Identify intact LTR-RTs eligible for conversion.
    intact_indices = []
    for i, rec in enumerate(bed_records):
        feature_id, supplemental = parse_feature_ID(rec["name"])
        # Basic intact criteria.
        if "LTRlen" not in feature_id or "_SOLO" in feature_id:
            continue
        # New Rule 1: If supplemental info is present, check the first supplemental field.
        if supplemental and "CUT_BY" in supplemental[0]:
            continue
        # New Rule 2: For records that have supplemental info (i.e. at least one semicolon in NAME)
        # and that have a TSD (from column 5) and a strand (column 6), scan ±100 lines.
        # Only consider duplicates if both the candidate record and the compared record have supplemental info.
        if supplemental and rec["tsd"] and rec["strand"]:
            duplicate_found = False
            start_window = max(0, i - 100)
            end_window = min(len(bed_records), i + 101)
            for j in range(start_window, end_window):
                if j == i:
                    continue
                other = bed_records[j]
                # Only consider other records that also have supplemental info.
                if ";" not in other["name"]:
                    continue
                # If one record's NAME is contained in the other's (i.e. a partial match)
                # and both TSD and strand are identical, then mark as duplicate.
                if ((rec["name"] in other["name"]) or (other["name"] in rec["name"])) \
                   and (other["tsd"] == rec["tsd"]) \
                   and (other["strand"] == rec["strand"]):
                    duplicate_found = True
                    break
            if duplicate_found:
                continue

        intact_indices.append(i)

    # Determine how many to convert.
    n_intact = len(intact_indices)
    n_convert = int(round(n_intact * (args.soloLTR_freq / 100.0)))
    convert_indices = set(random.sample(intact_indices, n_convert)) if n_convert > 0 else set()

    # Process conversion events: compute deletion info and update BED records.
    deletions_by_chrom = defaultdict(list)
    for i, rec in enumerate(bed_records):
        if i in convert_indices:
            S = rec["start"]
            E = rec["end"]
            feature_id, other_attrs = parse_feature_ID(rec["name"])
            L = parse_LTRlen(feature_id)
            if L is None:
                sys.exit(f"Error: Could not parse LTRlen in feature_ID: {feature_id}")
            deletion_start = S + L
            deletion_end = E
            deletion_length = E - (S + L)
            if deletion_length < 0:
                sys.exit(f"Error: Negative deletion length for record: {rec}")
            deletions_by_chrom[rec["chrom"]].append({
                "dstart": deletion_start,
                "dend": deletion_end,
                "length": deletion_length
            })
            # Update the BED record: new end is start + LTRlen.
            rec["end"] = S + L
            # Append _SOLO to the feature_ID and reassemble the NAME.
            new_feature_id = update_feature_ID(feature_id)
            rec["name"] = ";".join([new_feature_id] + other_attrs)

    # For each chromosome, sort deletion events by starting coordinate.
    for chrom in deletions_by_chrom:
        deletions_by_chrom[chrom].sort(key=lambda x: x["dstart"])

    # Function to compute total deletion offset before a given coordinate.
    def get_offset(chrom, x):
        offset = 0
        if chrom not in deletions_by_chrom:
            return 0
        for d in deletions_by_chrom[chrom]:
            if d["dstart"] < x:
                offset += d["length"]
            else:
                break
        return offset

    # Update BED record coordinates for all entries.
    for rec in bed_records:
        chrom = rec["chrom"]
        offset_start = get_offset(chrom, rec["start"])
        offset_end = get_offset(chrom, rec["end"])
        rec["start"] = rec["start"] - offset_start
        rec["end"] = rec["end"] - offset_end

    # Update genome FASTA sequences per chromosome by removing deletion segments.
    updated_genome = {}
    for chrom in genome_dict:
        desc, seq = genome_dict[chrom]
        if chrom in deletions_by_chrom:
            new_seq_parts = []
            last_index = 0
            for d in deletions_by_chrom[chrom]:
                dstart = d["dstart"]
                dend = d["dend"]
                new_seq_parts.append(seq[last_index:dstart])
                last_index = dend
            new_seq_parts.append(seq[last_index:])
            new_seq = "".join(new_seq_parts)
        else:
            new_seq = seq
        updated_genome[chrom] = (desc, new_seq)

    # Write updated FASTA and BED files.
    out_fasta = args.output + ".fasta"
    out_bed = args.output + ".bed"
    write_fasta(updated_genome, out_fasta)
    with open(out_bed, "w") as bedout:
        for rec in bed_records:
            # Write in the same column order: chrom, start, end, NAME, TSD, strand.
            bedout.write("\t".join([
                rec["chrom"],
                str(rec["start"]),
                str(rec["end"]),
                rec["name"],
                rec["tsd"],
                rec["strand"]
            ]) + "\n")
    
    sys.stderr.write(f"Processed {len(bed_records)} BED entries; converted {len(convert_indices)} intact LTR-RTs to soloLTRs.\n")
    sys.stderr.write(f"Output FASTA: {out_fasta}\n")
    sys.stderr.write(f"Output BED: {out_bed}\n")

if __name__ == "__main__":
    main()
