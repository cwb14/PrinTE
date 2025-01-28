#!/usr/bin/env python3

import argparse
import sys
import csv
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description='Append alignment features to GFF attributes.')
    parser.add_argument('-aln', '--alignment', required=True, help='Path to alignmentfile.tsv')
    parser.add_argument('-gff', '--gff', required=True, help='Path to annotation.gff')
    parser.add_argument('-model', '--model', choices=['K2P', 'raw', 'JC69'], default='K2P',
                        help='Mode to select Div and Time columns (default: K2P)')
    return parser.parse_args()

def extract_te_id(gff_attributes):
    """Extract TE_ID from the GFF attributes using the Name= field."""
    match = re.search(r'Name=([^;]+)', gff_attributes)
    if match:
        return match.group(1)
    return None

def load_alignment(alignment_file):
    """Load alignment TSV and return a dict mapping TE_ID to its data."""
    alignment_dict = {}
    with open(alignment_file, 'r') as aln_fh:
        reader = csv.DictReader(aln_fh, delimiter='\t')
        for row in reader:
            # Extract TE_ID from sseqid or qseqid
            # Assuming TE_ID is after the last '_' in sseqid or qseqid
            # Example: ATLANTYS2#LTR/Ty3_TE0000039#LTR_retrotransposon
            sseqid = row['sseqid']
            te_id_match = re.search(r'_TE(\d+)', sseqid)
            if te_id_match:
                te_id = f"TE{te_id_match.group(1)}"
                alignment_dict[te_id] = row
            else:
                # Alternatively, extract after 'Name=' in GFF, which is simpler
                # But since alignmentfile.tsv may have different TE_ID format,
                # Ensure the matching TE_ID extraction logic is correct
                # For the example, TE_ID is TE0000039
                te_id_search = re.search(r'TE(\d+)', sseqid)
                if te_id_search:
                    te_id = f"TE{te_id_search.group(1)}"
                    alignment_dict[te_id] = row
    return alignment_dict

def load_alignment_correct(alignment_file):
    """Load alignment TSV and return a dict mapping TE_ID to its data."""
    alignment_dict = {}
    with open(alignment_file, 'r') as aln_fh:
        reader = csv.DictReader(aln_fh, delimiter='\t')
        for row in reader:
            # Extract TE_ID from sseqid or qseqid
            # Assuming TE_ID is after 'Name=' in GFF, which corresponds to part after '#' and before '_'
            # Given sseqid: ATLANTYS2#LTR/Ty3_TE0000039#LTR_retrotransposon
            # TE_ID: TE0000039
            # Extract 'TE0000039' from sseqid
            sseqid = row['sseqid']
            te_id_match = re.search(r'(TE\d+)', sseqid)
            if te_id_match:
                te_id = te_id_match.group(1)
                alignment_dict[te_id] = row
    return alignment_dict

def append_features_to_gff(gff_file, alignment_dict, model, output_fh):
    """Process GFF file and append features based on alignment_dict and model."""
    with open(gff_file, 'r') as gff_fh:
        for line in gff_fh:
            line = line.rstrip('\n')
            if line.startswith('#') or not line.strip():
                # Comment or empty line, write as is
                print(line, file=output_fh)
                continue
            fields = line.split('\t')
            if len(fields) != 9:
                # Not a valid GFF line, write as is
                print(line, file=output_fh)
                continue
            attributes = fields[8]
            te_id = extract_te_id(attributes)
            if te_id and te_id in alignment_dict:
                aln_data = alignment_dict[te_id]
                # Add CalcLTRlen
                calc_ltrlen = f"CalcLTRlen:{aln_data['tot_len']}"
                # Select Div and Time based on model
                div_key = f"{model}_d"
                time_key = f"{model}_T"
                div_value = aln_data.get(div_key, 'NA')
                time_value = aln_data.get(time_key, 'NA')
                div = f"Div:{div_value}"
                time = f"Time:{time_value}"
                # Append to attributes, separated by semicolon
                new_attributes = f"{attributes};{calc_ltrlen};{div};{time}"
                fields[8] = new_attributes
                new_line = '\t'.join(fields)
                print(new_line, file=output_fh)
            else:
                # No matching TE_ID, write line as is
                print(line, file=output_fh)

def main():
    args = parse_arguments()
    # Load alignment data
    alignment_dict = load_alignment_correct(args.alignment)
    # Process GFF and append features
    append_features_to_gff(args.gff, alignment_dict, args.model, sys.stdout)

if __name__ == "__main__":
    main()
