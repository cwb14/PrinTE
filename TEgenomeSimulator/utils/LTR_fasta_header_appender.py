import argparse
from Bio import SeqIO
import sys

def main():
    parser = argparse.ArgumentParser(description='Append LTR information to FASTA headers using semicolon separators.')
    parser.add_argument('-fasta', required=True, help='Input FASTA file')
    parser.add_argument('-domains', required=True, help='Domains TSV file')
    parser.add_argument('-div_type', default='K2P', choices=['raw', 'K2P', 'JC69', 'none'],
                        help='Divergence type (raw, K2P, JC69, none)')
    args = parser.parse_args()

    # Read domains data into a dictionary
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

    # Process each FASTA record
    records = []
    for record in SeqIO.parse(args.fasta, 'fasta'):
        qseqid_fasta = record.description
        if qseqid_fasta in domains_data:
            data = domains_data[qseqid_fasta]
            # Format with semicolons; exclude LTRdiv if div_type is 'none'
            if args.div_type == 'none':
                new_desc = f"{qseqid_fasta};LTRlen:{data['tot_len']}"
            else:
                new_desc = f"{qseqid_fasta};LTRlen:{data['tot_len']};LTRdiv:{data['div']}"
            record.id = new_desc
            record.description = new_desc
        records.append(record)

    # Write the modified records to standard output
    SeqIO.write(records, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
