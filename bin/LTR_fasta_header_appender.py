import argparse
import re
from Bio import SeqIO
import sys

def main():
    parser = argparse.ArgumentParser(
        description='Append LTR information to FASTA headers using semicolon separators.')
    parser.add_argument('-fasta', required=True, help='Input FASTA file')
    parser.add_argument('-domains', required=True, help='Domains TSV file')
    parser.add_argument('-div_type', default='K2P', choices=['raw', 'K2P', 'JC69', 'none'],
                        help='Divergence type (raw, K2P, JC69, none)')
    # Only applies to LTR elements
    parser.add_argument(
        '-exclude_no_hits',
        action='store_true',
        help='Exclude LTR-RT FASTA sequences (class_part=LTR) that do not occur in domains file'
    )
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
    class_parts = []  # track class_part for each record
    name_counts = {}

    for record in SeqIO.parse(args.fasta, 'fasta'):
        original_header = record.description
        original_headers.append(original_header)

        # Split the header into name, class, and superfamily parts
        parts = original_header.split('#', 1)
        name_part = parts[0] if parts[0].strip() else 'Unknown'
        class_part = superfamily_part = 'Unknown'
        remaining_after_name = parts[1] if len(parts) > 1 else ''

        if remaining_after_name:
            class_super = remaining_after_name.split('/', 1)
            class_part = class_super[0].strip() if class_super[0].strip() else 'Unknown'
            remaining_after_class = class_super[1] if len(class_super) > 1 else ''
            if remaining_after_class:
                super_match = re.match(r'^([^\s]*)', remaining_after_class)
                superfamily_part = super_match.group(1).strip() if super_match and super_match.group(1).strip() else 'Unknown'

        class_parts.append(class_part)

        # Generate corrected header (without uniqueness yet)
        corrected_header = f"{name_part}#{class_part}/{superfamily_part}"

        # Ensure name_part is unique
        base_name = name_part
        if base_name in name_counts:
            name_counts[base_name] += 1
            new_name = f"{base_name}_{name_counts[base_name] - 1}"
        else:
            name_counts[base_name] = 1
            new_name = base_name

        # Create new header with unique name
        new_header = f"{new_name}#{class_part}/{superfamily_part}"
        record.id = new_header
        record.description = new_header
        corrected_records.append(record)

    # Prepare final records list based on the -exclude_no_hits flag
    final_records = []
    for original_header, record, class_part in zip(original_headers, corrected_records, class_parts):

        # Only exclude if:
        #   - -exclude_no_hits is set
        #   - this is an LTR element
        #   - and the original header has no domain hit
        if args.exclude_no_hits and class_part == 'LTR' and original_header not in domains_data:
            continue

        # Append LTR information if available (for any class where there *is* a hit)
        if original_header in domains_data:
            data = domains_data[original_header]
            if args.div_type == 'none':
                new_desc = f"{record.description}~LTRlen:{data['tot_len']}"
            else:
                new_desc = f"{record.description}~LTRlen:{data['tot_len']}~LTRdiv:{data['div']}"
            record.id = new_desc
            record.description = new_desc

        final_records.append(record)

    # Write the modified records to standard output
    SeqIO.write(final_records, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
