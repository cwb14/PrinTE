#!/usr/bin/env python3

import argparse
import random
import sys
import multiprocessing
from functools import partial

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary with headers as keys and sequences as values.
    """
    fasta_dict = {}
    header = None
    seq_chunks = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # skip empty lines
                if line.startswith('>'):
                    if header:
                        fasta_dict[header] = ''.join(seq_chunks).upper()
                    header = line[1:].split()[0]  # Take the first word after '>'
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
            if header:
                fasta_dict[header] = ''.join(seq_chunks).upper()
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        sys.exit(1)
    return fasta_dict

def write_fasta(file_path, fasta_dict):
    """
    Writes a dictionary of sequences to a FASTA file.
    """
    with open(file_path, 'w') as f:
        for header, sequence in fasta_dict.items():
            f.write(f">{header}\n")
            # Write sequence in lines of 60 characters
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')

def process_chromosome(args):
    """
    Worker function to perform insertions on a single chromosome.
    """
    chrom_name, chrom_seq, num_insertions, ltr_headers, ltr_dict = args
    inserted_ltrs = []

    if num_insertions == 0:
        return chrom_name, ''.join(chrom_seq), inserted_ltrs

    seq = chrom_seq.copy()
    chr_length = len(seq)

    if chr_length < 5:
        # Cannot perform insertions
        return chrom_name, ''.join(seq), inserted_ltrs

    # Select insertion positions
    insertion_positions = [random.randint(5, len(seq)) for _ in range(num_insertions)]
    # Sort positions in descending order to avoid index shifting
    insertion_positions_sorted = sorted(insertion_positions, reverse=True)

    for insertion_pos in insertion_positions_sorted:
        # Extract TSD
        tsd_start = insertion_pos - 5
        tsd_end = insertion_pos
        tsd = ''.join(seq[tsd_start:tsd_end])

        # Select LTR
        selected_ltr = random.choice(ltr_headers)
        ltr_seq = ltr_dict[selected_ltr]

        # Insert LTR after TSD
        seq[insertion_pos:insertion_pos] = list(ltr_seq)

        # Insert duplicated TSD after LTR
        duplicated_tsd = list(tsd)
        seq[insertion_pos + len(ltr_seq):insertion_pos + len(ltr_seq)] = duplicated_tsd

        # Log inserted LTR
        inserted_ltrs.append({
            'chromosome': chrom_name,
            'insertion_position': insertion_pos,
            'tsd': tsd,
            'inserted_ltr': selected_ltr
        })

    return chrom_name, ''.join(seq), inserted_ltrs

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Randomly insert LTR sequences into a genome with TSD duplication.")
    parser.add_argument('-genome', required=True, help='Path to the genome FASTA file (genome.fa)')
    parser.add_argument('-LTR', required=True, help='Path to the LTR FASTA file (LTR.fa)')
    parser.add_argument('-n', type=int, required=True, help='Number of LTR insertions to perform')
    parser.add_argument('-output', required=True, help='Path to output the modified genome FASTA file')

    args = parser.parse_args()

    genome_path = args.genome
    ltr_path = args.LTR
    num_insertions = args.n
    output_genome_path = args.output

    # Read genome and LTR sequences
    print("Reading genome FASTA file...")
    genome_dict = parse_fasta(genome_path)
    if not genome_dict:
        print("Error: Genome FASTA file is empty or improperly formatted.")
        sys.exit(1)
    print(f"Loaded {len(genome_dict)} chromosome(s) from genome.")

    print("Reading LTR FASTA file...")
    ltr_dict = parse_fasta(ltr_path)
    if not ltr_dict:
        print("Error: LTR FASTA file is empty or improperly formatted.")
        sys.exit(1)
    ltr_headers = list(ltr_dict.keys())
    print(f"Loaded {len(ltr_headers)} LTR sequence(s).")

    chromosomes = list(genome_dict.keys())
    if not chromosomes:
        print("Error: No chromosomes found in the genome FASTA file.")
        sys.exit(1)

    # Assign insertions to chromosomes randomly
    print(f"Assigning {num_insertions} insertions to chromosomes...")
    insertions_per_chrom = {chrom: 0 for chrom in chromosomes}
    for _ in range(num_insertions):
        chrom = random.choice(chromosomes)
        insertions_per_chrom[chrom] += 1

    # Prepare arguments for multiprocessing
    pool_args = []
    for chrom, count in insertions_per_chrom.items():
        chrom_seq = list(genome_dict[chrom])
        pool_args.append((chrom, chrom_seq, count, ltr_headers, ltr_dict))

    # Initialize multiprocessing pool
    num_processes = multiprocessing.cpu_count()
    print(f"Starting multiprocessing with {num_processes} processes...")
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(process_chromosome, pool_args)

    # Collect modified genome and insertion logs
    modified_genome = {}
    all_inserted_ltrs = []
    for chrom_name, modified_seq, inserted_ltrs in results:
        modified_genome[chrom_name] = modified_seq
        all_inserted_ltrs.extend(inserted_ltrs)

    # Write the modified genome to the specified output FASTA file
    print(f"Writing modified genome to {output_genome_path}...")
    write_fasta(output_genome_path, modified_genome)

    # Logging insertion details
    total_insertions_done = len(all_inserted_ltrs)
    print(f"Total insertions performed: {total_insertions_done}")

    if total_insertions_done > 0:
        # Log details of the first insertion globally
        first_insertion = all_inserted_ltrs[0]
        print("First Insertion Details:")
        print(f"  Chromosome: {first_insertion['chromosome']}")
        print(f"  Insertion Position: {first_insertion['insertion_position']}")
        print(f"  TSD: {first_insertion['tsd']}")
        print(f"  Inserted LTR: {first_insertion['inserted_ltr']}")
        print("-" * 50)

    print("Insertion process completed successfully.")

if __name__ == "__main__":
    main()
