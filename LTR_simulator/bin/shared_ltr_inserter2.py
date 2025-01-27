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

def process_chromosome_single(args):
    """
    Worker function to perform insertions on a single chromosome in single-genome mode.
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

def process_chromosome_double(args):
    """
    Worker function to perform insertions on a pair of chromosomes in two-genome mode.
    Inserts LTRs at the same positions in both genomes.
    """
    chrom_name, chrom_seq1, chrom_seq2, num_insertions, ltr_headers, ltr_dict = args
    inserted_ltrs = []

    if num_insertions == 0:
        return (chrom_name, ''.join(chrom_seq1), ''.join(chrom_seq2), inserted_ltrs)

    seq1 = chrom_seq1.copy()
    seq2 = chrom_seq2.copy()
    chr_length1 = len(seq1)
    chr_length2 = len(seq2)

    if chr_length1 < 5 or chr_length2 < 5:
        # Cannot perform insertions
        return (chrom_name, ''.join(seq1), ''.join(seq2), inserted_ltrs)

    # Select insertion positions based on genome1
    insertion_positions = [random.randint(5, len(seq1)) for _ in range(num_insertions)]
    # Sort positions in descending order to avoid index shifting
    insertion_positions_sorted = sorted(insertion_positions, reverse=True)

    for insertion_pos in insertion_positions_sorted:
        # Extract TSD from genome1
        tsd_start1 = insertion_pos - 5
        tsd_end1 = insertion_pos
        tsd1 = ''.join(seq1[tsd_start1:tsd_end1])

        # Extract TSD from genome2
        tsd_start2 = insertion_pos - 5
        tsd_end2 = insertion_pos
        tsd2 = ''.join(seq2[tsd_start2:tsd_end2])

        # Select LTR
        selected_ltr = random.choice(ltr_headers)
        ltr_seq = ltr_dict[selected_ltr]

        # Insert LTR into genome1 after TSD
        seq1[insertion_pos:insertion_pos] = list(ltr_seq)
        # Insert duplicated TSD into genome1 after LTR
        duplicated_tsd1 = list(tsd1)
        seq1[insertion_pos + len(ltr_seq):insertion_pos + len(ltr_seq)] = duplicated_tsd1

        # Insert LTR into genome2 after TSD
        seq2[insertion_pos:insertion_pos] = list(ltr_seq)
        # Insert duplicated TSD into genome2 after LTR
        duplicated_tsd2 = list(tsd2)
        seq2[insertion_pos + len(ltr_seq):insertion_pos + len(ltr_seq)] = duplicated_tsd2

        # Log inserted LTR
        inserted_ltrs.append({
            'chromosome': chrom_name,
            'insertion_position': insertion_pos,
            'tsd_genome1': tsd1,
            'tsd_genome2': tsd2,
            'inserted_ltr': selected_ltr
        })

    return (chrom_name, ''.join(seq1), ''.join(seq2), inserted_ltrs)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Randomly insert LTR sequences into one or two genomes with TSD duplication.")
    
    # Define mutually exclusive groups for single and two-genome modes
    group = parser.add_mutually_exclusive_group(required=True)
    
    # Single-genome mode arguments
    group.add_argument('-genome', help='Path to the genome FASTA file (genome.fa)')
    
    # Two-genome mode arguments
    group.add_argument('-genome1', help='Path to the first genome FASTA file (genome1.fa)')
    
    # Add remaining arguments
    parser.add_argument('-genome2', help='Path to the second genome FASTA file (genome2.fa)')
    parser.add_argument('-LTR', required=True, help='Path to the LTR FASTA file (LTR.fa)')
    parser.add_argument('-n', type=int, required=True, help='Number of LTR insertions to perform')
    parser.add_argument('-output', help='Path to output the modified genome FASTA file (single-genome mode)')
    parser.add_argument('-output1', help='Path to output the modified first genome FASTA file (two-genome mode)')
    parser.add_argument('-output2', help='Path to output the modified second genome FASTA file (two-genome mode)')

    args = parser.parse_args()

    # Determine mode based on provided arguments
    single_genome_mode = False
    two_genome_mode = False

    if args.genome:
        single_genome_mode = True
        # Ensure output is provided
        if not args.output:
            print("Error: -output is required in single-genome mode.")
            sys.exit(1)
    elif args.genome1:
        two_genome_mode = True
        # Ensure genome2 and outputs are provided
        if not args.genome2 or not args.output1 or not args.output2:
            print("Error: -genome2, -output1, and -output2 are required in two-genome mode.")
            sys.exit(1)
    else:
        print("Error: Either -genome or -genome1 must be specified.")
        sys.exit(1)

    # Common arguments
    ltr_path = args.LTR
    num_insertions = args.n

    # Read LTR sequences
    print("Reading LTR FASTA file...")
    ltr_dict = parse_fasta(ltr_path)
    if not ltr_dict:
        print("Error: LTR FASTA file is empty or improperly formatted.")
        sys.exit(1)
    ltr_headers = list(ltr_dict.keys())
    print(f"Loaded {len(ltr_headers)} LTR sequence(s).")

    if single_genome_mode:
        genome_path = args.genome
        output_genome_path = args.output

        # Read genome sequences
        print("Reading genome FASTA file...")
        genome_dict = parse_fasta(genome_path)
        if not genome_dict:
            print("Error: Genome FASTA file is empty or improperly formatted.")
            sys.exit(1)
        print(f"Loaded {len(genome_dict)} chromosome(s) from genome.")

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
            results = pool.map(process_chromosome_single, pool_args)

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

    elif two_genome_mode:
        genome1_path = args.genome1
        genome2_path = args.genome2
        output1_genome_path = args.output1
        output2_genome_path = args.output2

        # Read both genome sequences
        print("Reading first genome FASTA file...")
        genome1_dict = parse_fasta(genome1_path)
        if not genome1_dict:
            print("Error: First genome FASTA file is empty or improperly formatted.")
            sys.exit(1)
        print(f"Loaded {len(genome1_dict)} chromosome(s) from first genome.")

        print("Reading second genome FASTA file...")
        genome2_dict = parse_fasta(genome2_path)
        if not genome2_dict:
            print("Error: Second genome FASTA file is empty or improperly formatted.")
            sys.exit(1)
        print(f"Loaded {len(genome2_dict)} chromosome(s) from second genome.")

        # Verify that both genomes have the same chromosomes and lengths
        chromosomes1 = list(genome1_dict.keys())
        chromosomes2 = list(genome2_dict.keys())

        if len(chromosomes1) != len(chromosomes2):
            print("Error: The two genomes have different numbers of chromosomes.")
            sys.exit(1)

        chromosomes1_sorted = sorted(chromosomes1)
        chromosomes2_sorted = sorted(chromosomes2)

        for chrom1, chrom2 in zip(chromosomes1_sorted, chromosomes2_sorted):
            if chrom1 != chrom2:
                print(f"Error: Chromosome names do not match: {chrom1} vs {chrom2}.")
                sys.exit(1)
            if len(genome1_dict[chrom1]) != len(genome2_dict[chrom2]):
                print(f"Error: Chromosome lengths do not match for {chrom1}.")
                sys.exit(1)

        chromosomes = chromosomes1_sorted
        print("Verification passed: Both genomes have the same chromosomes with identical lengths.")

        # Assign insertions to chromosomes randomly
        print(f"Assigning {num_insertions} insertions to chromosomes...")
        insertions_per_chrom = {chrom: 0 for chrom in chromosomes}
        for _ in range(num_insertions):
            chrom = random.choice(chromosomes)
            insertions_per_chrom[chrom] += 1

        # Prepare arguments for multiprocessing
        pool_args = []
        for chrom in chromosomes:
            count = insertions_per_chrom[chrom]
            chrom_seq1 = list(genome1_dict[chrom])
            chrom_seq2 = list(genome2_dict[chrom])
            pool_args.append((chrom, chrom_seq1, chrom_seq2, count, ltr_headers, ltr_dict))

        # Initialize multiprocessing pool
        num_processes = multiprocessing.cpu_count()
        print(f"Starting multiprocessing with {num_processes} processes...")
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.map(process_chromosome_double, pool_args)

        # Collect modified genomes and insertion logs
        modified_genome1 = {}
        modified_genome2 = {}
        all_inserted_ltrs = []
        for result in results:
            chrom_name, modified_seq1, modified_seq2, inserted_ltrs = result
            modified_genome1[chrom_name] = modified_seq1
            modified_genome2[chrom_name] = modified_seq2
            all_inserted_ltrs.extend(inserted_ltrs)

        # Write the modified genomes to the specified output FASTA files
        print(f"Writing modified first genome to {output1_genome_path}...")
        write_fasta(output1_genome_path, modified_genome1)

        print(f"Writing modified second genome to {output2_genome_path}...")
        write_fasta(output2_genome_path, modified_genome2)

        # Logging insertion details
        total_insertions_done = len(all_inserted_ltrs)
        print(f"Total insertions performed: {total_insertions_done}")

        if total_insertions_done > 0:
            # Log details of the first insertion globally
            first_insertion = all_inserted_ltrs[0]
            print("First Insertion Details:")
            print(f"  Chromosome: {first_insertion['chromosome']}")
            print(f"  Insertion Position: {first_insertion['insertion_position']}")
            print(f"  TSD Genome1: {first_insertion['tsd_genome1']}")
            print(f"  TSD Genome2: {first_insertion['tsd_genome2']}")
            print(f"  Inserted LTR: {first_insertion['inserted_ltr']}")
            print("-" * 50)

        print("Insertion process completed successfully.")

if __name__ == "__main__":
    main()
