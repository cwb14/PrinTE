#!/usr/bin/env python3
import sys
import re

def print_help():
    """
    Print help message.
    """
    help_text = """
    Usage: gff_to_bed.py input.gff output.bed

    This script converts TEgenomeSimulator GFF files to BED format for input into prinTE.
    
    Columns in the output BED file:
    1. Chromosome (from column 1 of GFF)
    2. Start position (from column 4 of GFF)
    3. End position (from column 5 of GFF)
    4. ID with modifications:
       - Drops everything after and including '_TE0'
       - Appends 'CUT_BY:' information if 'Cut_by=' is present
       - Appends 'NESTED_IN:' information if 'Nest_in=' is present
    5. TSD (Target Site Duplication):
       - Prefers 'TSD_5' over 'TSD_3' if available
       - Defaults to 'NA' if neither is present or is 'None'
    6. Strand (from column 7 of GFF)

    Example:
        gff_to_bed.py input.gff output.bed
    """
    print(help_text)
    sys.exit(1)

def extract_id_parts(id_val):
    """ Extracts base ID before '_TE0' and the TE number. """
    if '_TE' in id_val:
        parts = id_val.split('_TE')
        base_id = parts[0]
        match = re.search(r'^(0\d+)', parts[1])
        if match:
            te_id = "TE" + match.group(1)
        else:
            te_id = "TE" + parts[1].split('_')[0]
        return base_id, te_id
    else:
        return id_val, None

def parse_attributes(attr_str):
    """ Parses GFF attribute column into a dictionary. """
    attrs = {}
    for attr in attr_str.strip().split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = value
        elif ':' in attr:
            key, value = attr.split(':', 1)
            attrs[key] = value
    return attrs

def main():
    # Handle incorrect usage and display help message
    if len(sys.argv) != 3 or sys.argv[1] in ('-h', '--help'):
        print_help()
        
    gff_file = sys.argv[1]
    bed_file = sys.argv[2]
    
    # First pass: build a mapping of TE IDs
    te_mapping = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == "":
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            attr_dict = parse_attributes(parts[8])
            if "ID" in attr_dict:
                id_full = attr_dict["ID"]
                base_id, te_id = extract_id_parts(id_full)
                if te_id:
                    te_mapping[te_id] = base_id
                    # print(f"\tMapping created: {te_id} -> {base_id}")

    # Second pass: process GFF lines and write to BED
    with open(gff_file) as fin, open(bed_file, 'w') as fout:
        for line in fin:
            if line.startswith('#') or line.strip() == "":
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom = parts[0]
            start = parts[3]
            end = parts[4]
            strand = parts[6]
            attr_dict = parse_attributes(parts[8])
            
            # Extract main ID
            if "ID" not in attr_dict:
                # print("No ID attribute found in line:", line.strip())
                continue
            id_full = attr_dict["ID"]
            base_id, te_id = extract_id_parts(id_full)
            bed_id = base_id
            # print(f"\tExtracted base ID: {bed_id} from {id_full}")
            
            # Append Cut_by information
            if "Cut_by" in attr_dict:
                cut_by_full = attr_dict["Cut_by"]
                cut_base, _ = extract_id_parts(cut_by_full)
                bed_id += ";" + "CUT_BY:" + cut_base
                # print(f"\tAppended Cut_by: {cut_base}")
            
            # Append Nest_in information
            if "Nest_in" in attr_dict:
                nest_te = attr_dict["Nest_in"]
                if nest_te in te_mapping:
                    nest_base = te_mapping[nest_te]
                    bed_id += ";" + "NESTED_IN:" + nest_base
                    # print(f"\tAppended Nest_in: {nest_base} (mapped from {nest_te})")
                else:
                    # print(f"\tWarning: Nest_in TE id {nest_te} not found in mapping. Skipping Nest_in.")
                    pass
            
            # Extract TSD (Target Site Duplication)
            tsd = "NA"
            if "TSD_5" in attr_dict and attr_dict["TSD_5"] != "None":
                tsd = attr_dict["TSD_5"]
                # print(f"\tUsing TSD_5: {tsd}")
            elif "TSD_3" in attr_dict and attr_dict["TSD_3"] != "None":
                tsd = attr_dict["TSD_3"]
                # print(f"\tUsing TSD_3: {tsd}")
            else:
                # print("\tNo valid TSD found; using NA")
                pass
            
            # Write to BED file
            bed_line = "\t".join([chrom, start, end, bed_id, tsd, strand])
            fout.write(bed_line + "\n")
    
    print("Conversion complete. Output written to", bed_file)

if __name__ == '__main__':
    main()
