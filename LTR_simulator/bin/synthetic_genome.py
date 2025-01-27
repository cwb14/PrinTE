#!/usr/bin/env python3

import argparse
import sys
import random
import math

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
                else:
                    current_seq.append(line)
            if current_seq:
                sequences.append(''.join(current_seq).upper())
    except IOError as e:
        sys.exit(f"Error reading CDS file: {e}")
    return sequences

def distribute_cds_to_chromosomes(cds_list, chr_number):
    chr_cds = [[] for _ in range(chr_number)]
    for idx, cds in enumerate(cds_list):
        chr_index = idx % chr_number
        chr_cds[chr_index].append(cds)
    return chr_cds

def generate_random_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

def generate_chromosome_sequence(chrom_size, cds_list):
    if not cds_list:
        return generate_random_sequence(chrom_size)
    
    num_cds = len(cds_list)
    total_cds_length = sum(len(cds) for cds in cds_list)
    
    if total_cds_length > chrom_size:
        sys.exit("Error: Total length of CDS sequences exceeds chromosome size.")
    
    # Calculate the number of random regions (num_cds + 1)
    num_regions = num_cds + 1
    remaining_length = chrom_size - total_cds_length
    if num_regions == 0:
        region_length = 0
    else:
        # Distribute remaining_length as evenly as possible
        base_region_length = remaining_length // num_regions
        extra = remaining_length % num_regions
        region_lengths = [base_region_length + (1 if i < extra else 0) for i in range(num_regions)]
    
    chromosome = []
    for i in range(num_cds):
        # Add random sequence before CDS
        if num_regions > 0:
            chr_seq = generate_random_sequence(region_lengths[i])
            chromosome.append(chr_seq)
        # Add CDS
        chromosome.append(cds_list[i])
    # Add random sequence after last CDS
    if num_regions > 0:
        chr_seq = generate_random_sequence(region_lengths[-1])
        chromosome.append(chr_seq)
    
    return ''.join(chromosome)

def main():
    args = parse_arguments()
    
    # Set the random seed if provided
    if args.seed is not None:
        random.seed(args.seed)
    else:
        # Optional: Set a default seed for consistent behavior when seed is not provided
        # Comment out the next line if you prefer non-reproducible random genomes when seed is not set
        random.seed(42)
    
    try:
        total_size = parse_size(args.size)
    except ValueError as e:
        sys.exit(f"Error parsing size: {e}")
    
    chr_number = args.chr_number
    if chr_number <= 0:
        sys.exit("Number of chromosomes must be a positive integer.")
    
    chr_size = total_size // chr_number
    # Handle any remaining bases by adding to the last chromosome
    chr_sizes = [chr_size] * chr_number
    remaining = total_size - chr_size * chr_number
    if remaining > 0:
        chr_sizes[-1] += remaining
    
    # Read CDS sequences if provided
    if args.cds:
        cds_list = read_fasta(args.cds)
        if not cds_list:
            sys.exit("No CDS sequences found in the provided CDS file.")
        # Distribute CDS across chromosomes
        chr_cds = distribute_cds_to_chromosomes(cds_list, chr_number)
    else:
        chr_cds = [[] for _ in range(chr_number)]
    
    # Generate and output each chromosome
    for i in range(chr_number):
        chr_name = f">chr{i+1}"
        chr_size_current = chr_sizes[i]
        cds_current = chr_cds[i]
        chr_sequence = generate_chromosome_sequence(chr_size_current, cds_current)
        # Format the sequence to have lines of up to 60 characters
        formatted_sequence = '\n'.join([chr_sequence[j:j+60] for j in range(0, len(chr_sequence), 60)])
        print(chr_name)
        print(formatted_sequence)

if __name__ == "__main__":
        main()
