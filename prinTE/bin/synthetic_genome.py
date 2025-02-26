#!/usr/bin/env python3

import argparse
import sys
import random
import math
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Simulate a random genome FASTA file.')
    parser.add_argument('-size', type=str, default='300Mb',
                        help='Total genome size (e.g., 300Mb, 500kb, 2Gb). Default is 300Mb.')
    parser.add_argument('-cds', type=str, default=None,
                        help='Optional CDS FASTA file to insert into the genome. These are inserted at roughly evenly spaced intervals.')
    parser.add_argument('-chr_number', type=int, default=5,
                        help='Number of chromosomes. Default is 5.')
    parser.add_argument('-seed', type=int, default=None,
                        help='Optional random seed for reproducibility.')
    parser.add_argument('-out_prefix', type=str, default='output',
                        help="Output file prefix. Example: '-out_prefix genome' generates 'genome.fa', 'genome.cds', and 'genome.bed'. Default is 'output'.")
    # New arguments:
    parser.add_argument('-cds_num', type=int, default=None,
                        help='Specifies the number of CDS sequences to insert.')
    parser.add_argument('-cds_percent', type=float, default=None,
                        help='Specifies the percent of the genome that should be CDS.')
    return parser.parse_args()

def parse_size(size_str):
    size_str = size_str.strip().upper()
    if size_str.endswith('KB'):
        return int(float(size_str[:-2]) * 1e3)
    elif size_str.endswith('MB'):
        return int(float(size_str[:-2]) * 1e6)
    elif size_str.endswith('GB'):
        return int(float(size_str[:-2]) * 1e9)
    else:
        raise ValueError("Size must end with 'kb', 'Mb', or 'Gb'.")

def read_fasta(filepath):
    sequences = []
    headers = []
    current_seq = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(''.join(current_seq).upper())
                        current_seq = []
                    headers.append(line[1:])  # Store header without '>'
                else:
                    current_seq.append(line)
            if current_seq:
                sequences.append(''.join(current_seq).upper())
    except IOError as e:
        sys.exit(f"Error reading CDS file: {e}")
    return headers, sequences

def distribute_cds_to_chromosomes(cds_list, chr_number):
    chr_cds = [[] for _ in range(chr_number)]
    for idx, cds in enumerate(cds_list):
        chr_index = idx % chr_number
        gene_name = f"gene{idx+1}"
        chr_cds[chr_index].append((gene_name, cds))
    return chr_cds

def generate_random_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

def generate_chromosome_sequence(chrom_size, cds_list):
    chromosome = []
    gene_positions = []
    current_pos = 0  # 0-based position

    if not cds_list:
        chromosome_seq = generate_random_sequence(chrom_size)
        return chromosome_seq, gene_positions

    num_cds = len(cds_list)
    total_cds_length = sum(len(cds) for _, cds in cds_list)

    if total_cds_length > chrom_size:
        sys.exit("Error: Total length of CDS sequences exceeds chromosome size.")

    # Calculate the number of random regions (num_cds + 1)
    num_regions = num_cds + 1
    remaining_length = chrom_size - total_cds_length

    # Distribute remaining_length as evenly as possible
    base_region_length = remaining_length // num_regions
    extra = remaining_length % num_regions
    region_lengths = [base_region_length + (1 if i < extra else 0) for i in range(num_regions)]

    for i in range(num_cds):
        # Add random sequence before CDS
        if num_regions > 0:
            rand_length = region_lengths[i]
            rand_seq = generate_random_sequence(rand_length)
            chromosome.append(rand_seq)
            current_pos += rand_length

        # Add CDS
        gene_name, cds_seq = cds_list[i]
        chromosome.append(cds_seq)
        cds_length = len(cds_seq)
        start = current_pos
        end = current_pos + cds_length
        gene_positions.append((gene_name, start, end))
        current_pos += cds_length

    # Add random sequence after last CDS
    if num_regions > 0:
        rand_length = region_lengths[-1]
        rand_seq = generate_random_sequence(rand_length)
        chromosome.append(rand_seq)
        current_pos += rand_length

    chromosome_seq = ''.join(chromosome)
    return chromosome_seq, gene_positions

def write_fasta(filepath, sequences):
    with open(filepath, 'w') as f:
        for header, seq in sequences:
            f.write(f"{header}\n")
            # Write sequence in lines of up to 60 characters
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

def write_cds(filepath, gene_sequences):
    with open(filepath, 'w') as f:
        for gene_name, seq in gene_sequences:
            f.write(f">{gene_name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

def write_bed(filepath, bed_entries):
    with open(filepath, 'w') as f:
        for entry in bed_entries:
            chrom, start, end, gene = entry
            f.write(f"{chrom}\t{start}\t{end}\t{gene}\n")

def main():
    args = parse_arguments()

    # Check that -cds_num and -cds_percent are not both provided.
    if args.cds_num is not None and args.cds_percent is not None:
        sys.exit("Error: Cannot specify both -cds_num and -cds_percent.")

    # If either -cds_num or -cds_percent is provided, a CDS file must be provided.
    if (args.cds_num is not None or args.cds_percent is not None) and not args.cds:
        sys.exit("Error: -cds_num and -cds_percent require a CDS file provided with -cds.")

    # Set the random seed if provided
    if args.seed is not None:
        random.seed(args.seed)
    else:
        # Optional: Set a default seed for reproducible behavior if seed is not set.
        random.seed(42)

    try:
        total_size = parse_size(args.size)
    except ValueError as e:
        sys.exit(f"Error parsing size: {e}")

    cds_list = None
    if args.cds:
        original_headers, cds_seqs = read_fasta(args.cds)
        if not cds_seqs:
            sys.exit("No CDS sequences found in the provided CDS file.")

        # Process -cds_num
        if args.cds_num is not None:
            if args.cds_num > len(cds_seqs):
                print(f"Warning: Requested {args.cds_num} CDS but only {len(cds_seqs)} available. Using all available CDS.")
                cds_list = cds_seqs
            else:
                cds_list = random.sample(cds_seqs, args.cds_num)
        # Process -cds_percent
        elif args.cds_percent is not None:
            desired_percent = float(args.cds_percent)
            target_cds_total = (desired_percent / 100) * total_size
            available_cds_total = sum(len(seq) for seq in cds_seqs)
            if available_cds_total < target_cds_total:
                print(f"Warning: Available CDS total length ({available_cds_total} bases) is less than needed for {desired_percent}% of genome size {total_size} bases ({target_cds_total} bases). Using all available CDS.")
                cds_list = cds_seqs
            else:
                # Select a random subset (shuffling then accumulating until the target is met)
                random.shuffle(cds_seqs)
                chosen_cds = []
                cumulative = 0
                for seq in cds_seqs:
                    chosen_cds.append(seq)
                    cumulative += len(seq)
                    if cumulative >= target_cds_total:
                        break
                cds_list = chosen_cds
                # Adjust the genome size so that CDS equals exactly the desired percentage.
                new_total_size = int(round(cumulative * 100 / desired_percent))
                print(f"Adjusting genome size from {total_size} to {new_total_size} to achieve desired CDS percent of {desired_percent}%.")
                total_size = new_total_size
        else:
            # Neither -cds_num nor -cds_percent was provided; use all CDS sequences.
            cds_list = cds_seqs

        # Distribute CDS across chromosomes with gene names.
        chr_cds = distribute_cds_to_chromosomes(cds_list, args.chr_number)
    else:
        chr_cds = [[] for _ in range(args.chr_number)]

    # Determine chromosome sizes.
    chr_number = args.chr_number
    chr_size = total_size // chr_number
    chr_sizes = [chr_size] * chr_number
    remaining = total_size - chr_size * chr_number
    if remaining > 0:
        chr_sizes[-1] += remaining

    # Prepare output file paths.
    prefix = args.out_prefix
    fa_path = f"{prefix}.fa"
    cds_path = f"{prefix}.cds" if args.cds else None
    bed_path = f"{prefix}.bed" if args.cds else None

    fasta_entries = []
    all_gene_sequences = []  # List of tuples (gene_name, seq)
    bed_entries = []  # List of tuples (chrom, start, end, gene)

    # Generate and collect each chromosome.
    for i in range(chr_number):
        chr_name = f"chr{i+1}"
        chr_size_current = chr_sizes[i]
        cds_current = chr_cds[i]
        chr_sequence, gene_positions = generate_chromosome_sequence(chr_size_current, cds_current)

        # Add to FASTA entries.
        fasta_entries.append((f">{chr_name}", chr_sequence))

        if args.cds:
            # Collect gene sequences and BED entries.
            for gene_name, start, end in gene_positions:
                # Using gene_name (e.g. "gene3") to index into the overall chosen CDS.
                idx = int(gene_name.replace("gene", "")) - 1
                if idx < len(cds_list):
                    seq = cds_list[idx]
                    all_gene_sequences.append((gene_name, seq))
                else:
                    sys.exit(f"Gene name {gene_name} exceeds the number of CDS sequences.")
                bed_entries.append((chr_name, start, end, gene_name))

    # Write the FASTA file.
    write_fasta(fa_path, fasta_entries)
    print(f"Genome FASTA written to '{fa_path}'.")

    if args.cds:
        # Write the CDS file with renamed headers.
        write_cds(cds_path, all_gene_sequences)
        print(f"CDS file written to '{cds_path}'.")
        # Write the BED file.
        write_bed(bed_path, bed_entries)
        print(f"BED file written to '{bed_path}'.")

        # Compute and print the percent of the genome that is CDS.
        total_cds_bases = sum(end - start for _, start, end in (pos for entry in [generate_chromosome_sequence(chr_sizes[i], chr_cds[i])[1] for i in range(chr_number)] for pos in entry))
        # Alternatively, since each BED entry covers one CDS, you could also do:
        # total_cds_bases = sum(end - start for _, start, end, _ in bed_entries)
        percent_cds = (total_cds_bases / total_size) * 100
        print(f"CDS accounts for {percent_cds:.2f}% of the genomic landscape.")

if __name__ == "__main__":
    main()
