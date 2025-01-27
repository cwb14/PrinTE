#!/usr/bin/env python3

import argparse
import sys
import random
import numpy as np
from scipy.stats import beta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""
        This script simulates mutations on sequences from a multi-FASTA file.
        For each sequence, it calculates the minimum and maximum number of expected mutations based on
        the provided mutation rate, number of generations, and maximum divergence percentage.
        It then applies a random number of mutations within this range, sampled from a Beta distribution
        with a mean near the minimum and a long tail extending to the maximum. The steepness of the
        drop-off can be controlled via the shape parameter.

        Arguments:
          -fasta        Path to the input multi-FASTA file.
          -rate         Mutation rate per generation (e.g., 1.3e-8).
          -generations  Number of generations to simulate (e.g., 10000).
          -max_perc_div  The maximum divergence percentage (e.g., 10).
          -shape        Shape parameter for the distribution (0 to 1, default: 0.5).
                        Lower values create a steeper drop-off (more mutations near minimum).
                        Higher values create a more gradual decline towards the maximum.
          -out          Path to the output multi-FASTA file.
          --multiplier  Number of times to sample and mutate each input sequence (default: 1).
                        For example, '--multiplier 4' will generate four mutated versions of each input sequence.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-fasta', required=True, help='Path to the input multi-FASTA file.')
    parser.add_argument('-rate', type=float, required=True, help='Mutation rate per generation (e.g., 1.3e-8).')
    parser.add_argument('-generations', type=int, required=True, help='Number of generations to simulate (e.g., 10000).')
    parser.add_argument('-max_perc_div', type=float, required=True, help='The maximum divergence percentage (e.g., 10).')
    parser.add_argument('-shape', type=float, default=0.5, help='Shape parameter for the distribution (0 to 1, default: 0.5).')
    parser.add_argument('-out', required=True, help='Path to the output multi-FASTA file.')
    parser.add_argument('--multiplier', type=int, default=1, help='Number of times to sample and mutate each input sequence (default: 1).')

    return parser.parse_args()


def calculate_mutation_bounds(generations, mutation_rate, seq_length, max_perc_div):
    min_mut = generations * mutation_rate * seq_length
    max_mut = (seq_length * (max_perc_div / 100)) / 2

    print(f"Calculating mutation bounds:")
    print(f"  Generations          : {generations}")
    print(f"  Mutation rate       : {mutation_rate}")
    print(f"  Sequence length     : {seq_length}")
    print(f"  Max divergence (%)  : {max_perc_div}")
    print(f"  Minimum expected muts: {min_mut:.2f}")
    print(f"  Maximum expected muts: {max_mut:.2f}")

    if max_mut < min_mut:
        print("WARNING: Maximum divergence is lower than the minimum expected mutations.")
        print("         Adjusting maximum expected mutations to be slightly higher than minimum.")
        max_mut = min_mut * 1.1  # Ensure max_mut > min_mut
        print(f"  Adjusted maximum expected mutations: {max_mut:.2f}")

    return int(np.ceil(min_mut)), int(np.ceil(max_mut))


def sample_mutation_count(min_mut, max_mut, shape):
    """
    Samples the number of mutations using a Beta distribution scaled to [min_mut, max_mut].
    
    Parameters:
    - min_mut (int): Minimum number of mutations.
    - max_mut (int): Maximum number of mutations.
    - shape (float): Shape parameter between 0 and 1 controlling the steepness.

    Returns:
    - int: Number of mutations to introduce.
    """
    if not 0 <= shape <= 1:
        print("Error: Shape parameter must be between 0 and 1.")
        sys.exit(1)

    # Define Beta distribution parameters
    # alpha = 1 ensures the distribution is skewed towards 0 as beta increases
    alpha = 1
    beta_param = 1 + 9 * (1 - shape)  # shape=1 -> beta=1 (uniform), shape=0 -> beta=10 (steep)

    # Sample from Beta distribution
    sampled_value = beta.rvs(a=alpha, b=beta_param)

    # Scale sampled value to [min_mut, max_mut]
    mutation_count = min_mut + sampled_value * (max_mut - min_mut)

    return int(np.round(mutation_count))


def mutate_sequence(seq_record, min_mut, max_mut, shape, mutation_suffix):
    """
    Introduces mutations into a sequence record.

    Parameters:
    - seq_record (SeqRecord): The original sequence record.
    - min_mut (int): Minimum number of mutations.
    - max_mut (int): Maximum number of mutations.
    - shape (float): Shape parameter for mutation count distribution.
    - mutation_suffix (str): Suffix to append to the sequence ID for uniqueness.

    Returns:
    - SeqRecord: The mutated sequence record.
    """
    seq = list(str(seq_record.seq).upper())
    seq_length = len(seq)
    num_mutations = sample_mutation_count(min_mut, max_mut, shape)

    print(f"Mutating sequence '{seq_record.id}': {num_mutations} mutations will be introduced.")

    if num_mutations > seq_length:
        print(f"  Number of mutations ({num_mutations}) exceeds sequence length ({seq_length}). Adjusting to mutate every position.")
        num_mutations = seq_length

    # Select unique positions to mutate
    mutation_positions = random.sample(range(seq_length), num_mutations)
    nucleotides = ['A', 'C', 'G', 'T']

    for pos in mutation_positions:
        original_nt = seq[pos]
        if original_nt not in nucleotides:
            # Skip non-standard nucleotides
            continue
        possible_mutations = [nt for nt in nucleotides if nt != original_nt]
        new_nt = random.choice(possible_mutations)
        seq[pos] = new_nt

    mutated_seq = Seq(''.join(seq))
    new_id = f"{seq_record.id}{mutation_suffix}"
    return SeqRecord(mutated_seq, id=new_id, description=seq_record.description)


def main():
    args = parse_arguments()

    # Validate shape parameter
    if not 0 <= args.shape <= 1:
        print("Error: The shape parameter must be between 0 and 1.")
        sys.exit(1)

    # Validate multiplier
    if args.multiplier < 1:
        print("Error: The multiplier must be an integer greater than or equal to 1.")
        sys.exit(1)

    # Read input FASTA file
    try:
        fasta_sequences = list(SeqIO.parse(args.fasta, "fasta"))
        if not fasta_sequences:
            print(f"No sequences found in the FASTA file: {args.fasta}")
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: The file '{args.fasta}' does not exist.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)

    mutated_sequences = []

    # Process each sequence
    for seq_record in fasta_sequences:
        seq_length = len(seq_record.seq)
        min_mut, max_mut = calculate_mutation_bounds(
            args.generations,
            args.rate,
            seq_length,
            args.max_perc_div
        )

        for i in range(1, args.multiplier + 1):
            mutation_suffix = f"_mut{i}" if args.multiplier > 1 else ""
            mutated_seq_record = mutate_sequence(seq_record, min_mut, max_mut, args.shape, mutation_suffix)
            mutated_sequences.append(mutated_seq_record)

    # Write mutated sequences to output file
    try:
        SeqIO.write(mutated_sequences, args.out, "fasta")
        print(f"\nMutated sequences have been written to '{args.out}'.")
        print(f"Total input sequences : {len(fasta_sequences)}")
        print(f"Total mutated sequences: {len(mutated_sequences)}")
    except Exception as e:
        print(f"Error writing mutated FASTA file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
