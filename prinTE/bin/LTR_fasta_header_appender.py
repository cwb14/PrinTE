import argparse
import re
from Bio import SeqIO
import sys

def main():
    parser = argparse.ArgumentParser(description='Append LTR information to FASTA headers using semicolon separators.')
    parser.add_argument('-fasta', required=True, help='Input FASTA file')
    parser.add_argument('-domains', required=True, help='Domains TSV file')
    parser.add_argument('-div_type', default='K2P', choices=['raw', 'K2P', 'JC69', 'none'],
                        help='Divergence type (raw, K2P, JC69, none)')
    args = parser.parse_args()

    # Read domains data into a dictionary using original headers as keys
    domains_data = {}
    with open(args.domains, 'r') as domains_file:
        header = domains_file.readline().strip().split('\t')
        try:
            tot_len_idx = header.index('tot_len')
            div_col = f"{args.div_type}_d" if args.div_type != 'none' else None
            div_idx = header.index(div_col) if div_col else None
        except ValueError as e:
            sys.stderr.write(f"Error: Required column not found in domains file: {e}\n")
            sys.exit(1)

        for line in domains_file:
            fields = line.strip().split('\t')
            if len(fields) < tot_len_idx + 1 or (div_idx is not None and len(fields) < div_idx + 1):
                continue  # Skip lines that don't have enough columns
            qseqid = fields[0]
            tot_len = fields[tot_len_idx]
            div_value = fields[div_idx] if div_idx is not None else None
            domains_data[qseqid] = {'tot_len': tot_len, 'div': div_value}

    # Process each FASTA record to correct headers and ensure unique names
    corrected_records = []
    original_headers = []
    name_counts = {}

    for record in SeqIO.parse(args.fasta, 'fasta'):
        original_header = record.description
        original_headers.append(original_header)

        # Split the header into name, class, and superfamily parts
        name_part = class_part = superfamily_part = 'Unknown'
        parts = original_header.split('#', 1)
        name_part = parts[0] if parts else 'Unknown'
        remaining_after_name = parts[1] if len(parts) > 1 else ''

        if remaining_after_name:
            class_super = remaining_after_name.split('/', 1)
            class_part = class_super[0] if class_super else 'Unknown'
            remaining_after_class = class_super[1] if len(class_super) > 1 else ''
            if remaining_after_class:
                super_match = re.match(r'^([^\s]*)', remaining_after_class)
                superfamily_part = super_match.group(1) if super_match else 'Unknown'
            else:
                superfamily_part = 'Unknown'
        else:
            class_part = 'Unknown'
            superfamily_part = 'Unknown'

        # Ensure name_part is not empty
        if not name_part.strip():
            name_part = 'Unknown'

        # Ensure class_part and superfamily_part are not empty
        class_part = class_part.strip() if class_part.strip() else 'Unknown'
        superfamily_part = superfamily_part.strip() if superfamily_part.strip() else 'Unknown'

        # Generate corrected header (without uniqueness)
        corrected_header = f"{name_part}#{class_part}/{superfamily_part}"

        # Ensure name_part is unique
        base_name = name_part  # The original name from header processing
        if base_name in name_counts:
            name_counts[base_name] += 1
            new_name = f"{base_name}_{name_counts[base_name] - 1}"
        else:
            name_counts[base_name] = 1
            new_name = base_name

        # Create new_header with unique name
        new_header = f"{new_name}#{class_part}/{superfamily_part}"

        # Update the record's id and description
        record.id = new_header
        record.description = new_header
        corrected_records.append(record)

    # Append LTR information based on original headers
    for i, record in enumerate(corrected_records):
        original_header = original_headers[i]
        if original_header in domains_data:
            data = domains_data[original_header]
            if args.div_type == 'none':
                new_desc = f"{record.description}~LTRlen:{data['tot_len']}"
            else:
                new_desc = f"{record.description}~LTRlen:{data['tot_len']}~LTRdiv:{data['div']}"
            record.id = new_desc
            record.description = new_desc

    # Write the modified records to standard output
    SeqIO.write(corrected_records, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
