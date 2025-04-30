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
       - Uses the GFF 'Name' attribute (up to first '_')
       - Adds '_FRAG' if Integrity != 1.0, directly before the '#'
       - Appends '#<Classification>' from the GFF
       - Appends ';CUT_BY:<base>#<Classification>' if 'Cut_by=' is present
       - Appends ';NESTED_IN:<base>#<Classification>' if 'Nest_in=' is present
    5. TSD (Target Site Duplication):
       - Prefers 'TSD_5' over 'TSD_3' if available and not 'None'
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
    
    # First pass: build a mapping of TE IDs to Name-based classification
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
                _, te_id = extract_id_parts(id_full)
                if te_id:
                    # Use the GFF 'Name' attribute for mapping
                    name_val = attr_dict.get("Name")
                    if name_val:
                        name_base = name_val.split('_')[0]
                    else:
                        name_base = id_full.split('#')[0]
                    classification = attr_dict.get("Classification", "")
                    # Map TE ID to Name-based base and classification
                    te_mapping[te_id] = f"{name_base}#{classification}"

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
            
            # Extract NAME attribute for bed ID
            name_val = attr_dict.get("Name")
            if name_val:
                base_name = name_val.split('_')[0]
            else:
                # Fallback to ID parsing if Name missing
                id_full = attr_dict.get("ID")
                if not id_full:
                    continue
                base_name, _ = extract_id_parts(id_full)
                base_name = base_name.split('#')[0]

            # Determine integrity and fragmentation
            integrity = attr_dict.get("Integrity")
            frag_tag = "_FRAG" if integrity != "1.0" else ""

            # Classification string for this TE
            classification = attr_dict.get("Classification", "")

            # Construct the initial bed ID with fragmentation and classification
            bed_id = f"{base_name}{frag_tag}#{classification}"

            # Append Cut_by information
            if "Cut_by" in attr_dict:
                cut_by_full = attr_dict["Cut_by"]
                cut_base, _ = extract_id_parts(cut_by_full)
                bed_id += ";CUT_BY:" + cut_base

            # Append Nest_in information
            if "Nest_in" in attr_dict:
                nest_te = attr_dict["Nest_in"]
                if nest_te in te_mapping:
                    nest_base = te_mapping[nest_te]
                    bed_id += ";NESTED_IN:" + nest_base
                else:
                    pass

            # Extract TSD (Target Site Duplication)
            tsd = "NA"
            if "TSD_5" in attr_dict and attr_dict["TSD_5"] != "None":
                tsd = attr_dict["TSD_5"]
            elif "TSD_3" in attr_dict and attr_dict["TSD_3"] != "None":
                tsd = attr_dict["TSD_3"]

            # Write to BED file
            bed_line = "\t".join([chrom, start, end, bed_id, tsd, strand])
            fout.write(bed_line + "\n")
    
    print("Conversion complete. Output written to", bed_file)

if __name__ == '__main__':
    main()
