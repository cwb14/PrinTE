#!/usr/bin/env python3
"""
Random DNA Sequence Generator

This program generates random DNA sequences or combines reference sequences with random bases
to create new sequences. It supports multiprocessing for efficient sequence generation.

Author: Tianyu Lu (tlu83@wisc.edu)
Date: 2024-11-27
"""

import os
import argparse
import random
from typing import List, Optional
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

VERSION = "v0.2.0"
INFO = ("by Tianyu (Sky) Lu (tlu83@wisc.edu)\n"
        "released under LGPL-2.1 license")

# PRE-DEFINED PARAMETERS
NUCLEOTIDES = ["A", "T", "G", "C"]
PROCESS_CORE_RATIO = 1  # Number of processes = CPU cores * PROCESS_CORE_RATIO
DEFAULT_RANDOM_BASE_RATIO = 0.2
DEFAULT_ALLOCATED_CPU_CORES = os.cpu_count() - 2 if os.cpu_count() > 2 else 1
DEFAULT_REF_LIB_DIR_REL_PATH = "RefLib"
DEFAULT_OUTPUT_DIR_REL_PATH = "RandSeqGen-Result"

PROGRAM_ROOT_DIR_ABS_PATH = os.path.dirname(__file__)
DEFAULT_REF_LIB_ABS_PATH = os.path.join(PROGRAM_ROOT_DIR_ABS_PATH, DEFAULT_REF_LIB_DIR_REL_PATH)
DEFAULT_OUTPUT_DIR_ABS_PATH = os.path.join(os.getcwd(), DEFAULT_OUTPUT_DIR_REL_PATH)


# ======================================================================================================================


def process_length(length_str: str) -> int:
    """
    Process a length string to return the number of bases.

    This function converts various length formats (plain number, kb, mb) into
    the actual number of bases for sequence generation.

    Args:
        length_str (str): A string specifying sequence length. Can be:
            - A plain number (e.g., "100")
            - Kilobases with 'kb' suffix (e.g., "1kb")
            - Megabases with 'mb' suffix (e.g., "1mb")

    Returns:
        int: The number of bases represented by the input string.

    Raises:
        SystemExit: If the length string format is invalid.

    Examples:
        >>> process_length("100")
        100
        >>> process_length("1kb")
        1000
        >>> process_length("1mb")
        1000000
    """
    length_str = length_str.lower()

    # Check for kb/mb suffix
    if length_str.endswith("kb"):
        return int(float(length_str[:-2]) * 1000)
    elif length_str.endswith("mb"):
        return int(float(length_str[:-2]) * 1000000)

    # Check if pure number or b suffix
    if length_str.isdigit() or length_str.endswith("b"):
        return int(length_str.replace('b', ''))

    raise SystemExit("ERROR: Invalid length format. Use number, kb or mb (e.g. 100, 1kb, 1mb)")


def load_reference_sequences(ref_lib_abs_path: Optional[str], len_limit: Optional[int] = None) -> Optional[List[str]]:
    """
    Load reference sequences from all FASTA files in the reference directory.

    Searches for files with .fa or .fasta extensions in the reference directory
    and loads sequences that meet the length criteria into memory. Sequences are 
    stored as strings for efficient processing.

    Args:
        ref_lib_abs_path: Path to the reference sequence directory
        len_limit: Optional maximum length limit for sequences to load

    Raises:
        SystemExit: If no FASTA files are found in the reference directory
    """
    if ref_lib_abs_path is None:
        return None

    if not os.path.exists(ref_lib_abs_path):
        raise SystemExit(f"ERROR: {ref_lib_abs_path} does not exist!")

    if not os.path.isdir(ref_lib_abs_path):
        raise SystemExit(f"ERROR: {ref_lib_abs_path} is not a directory!")

    ref_sequences: List[str] = []
    for file in os.listdir(ref_lib_abs_path):
        if file.endswith((".fa", ".fasta")):
            file_path = os.path.join(ref_lib_abs_path, file)
            if len_limit is not None:
                ref_sequences.extend([str(record.seq) for record in SeqIO.parse(file_path, "fasta") if
                                      len(record.seq) <= len_limit])
            else:
                ref_sequences.extend([str(record.seq) for record in SeqIO.parse(file_path, "fasta")])

    ref_sequences.sort(key=len)
    return ref_sequences


def create_sequence_record(seq: str, id: str) -> SeqRecord:
    """
    Create a BioPython SeqRecord object for a sequence.

    Args:
        seq (str): The sequence string
        id (str): The id of the sequence

    Returns:
        SeqRecord: BioPython SeqRecord object for the sequence
    """
    return SeqRecord(Seq(seq), id=id, description="")


def generate_random_sequence(length: int) -> str:
    """
    Generate a random DNA sequence of specified length.

    Args:
        length (int): The desired length of the sequence

    Returns:
        str: A random DNA sequence composed of nucleotides from NUCLEOTIDES
    """
    return "".join(random.choices(NUCLEOTIDES, k=length))


def generate_random_segment_lengths(num_segments: int, total_length: int) -> List[int]:
    """
    Generate random segment lengths that sum to total_length.

    Args:
        num_segments (int): Number of segments to generate
        total_length (int): Total length to divide into segments

    Returns:
        List[int]: List of segment lengths that sum to total_length

    Note:
        This is used to determine the distribution of random nucleotides
        throughout the sequence in reference mode.
    """
    if num_segments <= 1:
        return [total_length]

    # Generate random points to split the sequence
    split_points: List[int] = sorted(random.sample(range(1, total_length), num_segments - 1))

    # Calculate segment lengths
    segments: List[int] = [split_points[0]]
    for i in range(1, len(split_points)):
        segments.append(split_points[i] - split_points[i - 1])
    segments.append(total_length - split_points[-1])

    return segments


def binary_search_suitable_refs(sorted_ref_sequences: List[str], len_limit: int) -> List[str]:
    """
    Find reference sequences that fit within length limit using binary search.

    Args:
        sorted_ref_sequences (List[str]): Pre-sorted list of reference sequences
        len_limit (int): Maximum length of reference sequences to find

    Returns:
        List[str]: List of reference sequences that fit within target_length

    Note:
        Uses binary search because reference sequences are pre-sorted by length.
    """
    left, right = 0, len(sorted_ref_sequences) - 1
    while left <= right:
        mid = (left + right) // 2
        if len(sorted_ref_sequences[mid]) <= len_limit:
            left = mid + 1
        else:
            right = mid - 1

    return sorted_ref_sequences[:left]


def save_single_batch(records: List[SeqRecord], batch: int, output_dir_path: str):
    """
    Save a single batch of sequences to a FASTA file.

    Args:
        records (List[SeqRecord]): List of BioPython SeqRecord objects to save
        batch (int): Current batch number
        output_dir_path (str): Directory path where the output files will be saved

    Note:
        Files are saved in FASTA format with names like 'sequences_batch_1.fa'
    """
    output_file = os.path.join(output_dir_path, f"sequences_batch_{batch}.fa")
    SeqIO.write(records, output_file, "fasta")
    print(f"Generated batch {batch}")


def _find_ref_lib_dir_abs_path(path: str) -> Optional[str]:
    """
    Locate the reference library directory.

    This function attempts to find the reference library directory in two locations:
    1. The directly specified path
    2. Under the default RefLib directory in the program root

    Args:
        path (str): The path to search for reference library. Can be:
            - An absolute path
            - A relative path
            - A directory name to look for in the default RefLib location

    Returns:
        str: The absolute path to the found reference library directory.

    Raises:
        SystemExit: If no valid reference library directory is found.

    Examples:
        >>> _find_ref_lib_dir_abs_path("/abs/path/to/refs")  # Returns path if valid
        '/abs/path/to/refs'
        >>> _find_ref_lib_dir_abs_path("rice_refs")  # Checks in default RefLib directory
        '/path/to/program/RefLib/rice_refs'
    """
    if path is None:
        return None

    # Check direct path
    if os.path.exists(path) and os.path.isdir(path):
        return os.path.abspath(path)

    # Check in built-in RefLib directory
    default_ref_lib_abs_path = os.path.join(DEFAULT_REF_LIB_ABS_PATH, path)
    if os.path.exists(default_ref_lib_abs_path) and os.path.isdir(default_ref_lib_abs_path):
        return default_ref_lib_abs_path

    raise SystemExit(f"ERROR: Specified reference library directory not found at {path} or built-in reference library!")


class SeqGenerator:
    """
    Class for generating random DNA sequences with parallelization.

    This class handles both purely random sequence generation and reference-based
    sequence generation with random base insertion.
    """

    def __init__(self, length: str, number: int, batch: int, cpu_cores: int, output_dir: str,
                 ref_lib: Optional[str] = None, random_ratio: float = DEFAULT_RANDOM_BASE_RATIO,
                 flag_verbose: bool = False):
        self.seq_length: int = process_length(length)
        self.seq_number: int = number
        self.batch_number: int = batch
        self.cpu_cores: int = cpu_cores
        self.output_dir_abs_path: str = os.path.abspath(output_dir)
        self.ref_lib_abs_path: str = _find_ref_lib_dir_abs_path(ref_lib)
        self.random_ratio: float = random_ratio
        self.flag_verbose: bool = flag_verbose

        self.ref_sequences: List[str] = load_reference_sequences(self.ref_lib_abs_path,
                                                                 self.seq_length // self.random_ratio)

        self._pre_check()

    def _pre_check(self):
        """
        Perform pre-execution checks and setup.

        This method:
        1. Creates the output directory if it doesn't exist
        2. Validates the random base ratio is between 0 and 1
        3. Checks if reference sequences are available when in reference mode
        4. Ensures sequence length and number are positive integers

        Raises:
            SystemExit: If any validation checks fail
        """
        if self.random_ratio < 0 or self.random_ratio > 1:
            raise ValueError("ERROR: Ratio of random bases must be between 0 and 1!")

        if self.cpu_cores < 1:
            raise ValueError("ERROR: Number of processors must be at least 1!")

        if self.ref_lib_abs_path is not None and len(self.ref_sequences) == 0:
            print("WARN: Length of sequence does not meet the requirement "
                  "(length of sequence must >= random_ratio * shortest reference sequence's length). "
                  "Referenced mode has been disabled and program will run in random mode.")
            self.ref_sequences = None

    def _generate_single_referenced_sequence(self) -> str:
        """
        Generate a sequence combining reference segments with random nucleotides.

        This method:
        1. Calculates random segment distribution
        2. Interleaves random segments with reference sequences
        3. Maintains intact reference sequences
        4. Falls back to random nucleotides when no suitable reference is found

        Returns:
            str: Generated sequence combining reference and random segments
        """
        random_length = int(self.seq_length * self.random_ratio)
        num_random_segments = max(2, int(self.seq_length / 1000))

        random_segment_lengths = generate_random_segment_lengths(num_random_segments, random_length)
        # ref_length = self.seq_length - random_length

        sequence: List[str] = []
        current_length = 0
        random_idx: int = 0

        while current_length < self.seq_length:
            # Add random segment
            if random_idx < len(random_segment_lengths):
                sequence.append(generate_random_sequence(random_segment_lengths[random_idx]))
                current_length += random_segment_lengths[random_idx]
                random_idx += 1

            # Try to add reference segment
            if current_length < self.seq_length:
                remaining = self.seq_length - current_length
                suitable_refs = binary_search_suitable_refs(self.ref_sequences, remaining)

                if suitable_refs:
                    ref_seq = random.choice(suitable_refs)
                    sequence.append(ref_seq)
                    current_length += len(ref_seq)
                else:
                    sequence.append(generate_random_sequence(remaining))
                    current_length += remaining

        return "".join(sequence)

    def _generate_single_sequence(self, seq_id: str, verbose: bool = False) -> SeqRecord:
        """
        Generate a single sequence with appropriate mode.

        Args:
            seq_id (str): ID for the sequence to be generated

        Returns:
            SeqRecord: Generated sequence as a BioPython SeqRecord object
        """

        if self.ref_sequences:
            seq: str = self._generate_single_referenced_sequence()
        else:
            seq: str = generate_random_sequence(self.seq_length)

        if verbose:
            print(f"  Generated: {seq_id}")
        else:
            print('*', end="")
        return create_sequence_record(seq, seq_id)

    def _generate_batch(self, batch: int, verbose: bool = False) -> List[SeqRecord]:
        """
        Generate a batch of sequences using multiprocessing.

        Args:
            batch (int): Batch number to generate

        Returns:
            List[SeqRecord]: List of generated sequences for this batch
        """
        print(f"\nGenerating batch {batch}...")

        mp_args_list = [(f"seq_{i}_batch_{batch}_len_{self.seq_length}", verbose) for i in range(self.seq_number)]

        with mp.Pool(self.cpu_cores * PROCESS_CORE_RATIO) as pool:
            sequences = pool.starmap(self._generate_single_sequence, mp_args_list)

        return sequences

    def execute(self):
        """
        Execute the entire sequence generation process.

        This method:
        1. Sets up multiprocessing pool based on CPU cores
        2. Generates sequences in batches using parallel processing
        3. Saves each batch to a separate FASTA file
        4. Handles both random and reference-based sequence generation
        """
        batch_results: List[List[SeqRecord]] = []
        for batch in range(self.batch_number):
            batch_results.append(self._generate_batch(batch, self.flag_verbose))

        print()
        os.makedirs(self.output_dir_abs_path, exist_ok=True)

        with mp.Pool(self.cpu_cores * PROCESS_CORE_RATIO) as pool:
            mp_args_list = [(records, batch, self.output_dir_abs_path) for batch, records in enumerate(batch_results)]
            pool.starmap(save_single_batch, mp_args_list)


# ======================================================================================================================


if __name__ == "__main__":
    """
    Main entry point for the sequence generator.
    """
    parser = argparse.ArgumentParser(prog="RandomSequenceGenerator")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}\n{INFO}")

    parser.add_argument("-l", "--length", help="Length of each sequence (e.g. 100, 1kb, 1mb)", required=True)
    parser.add_argument("-n", "--number", help="Number of sequences in each file (Optional, default: 1)",
                        type=int, default=1)
    parser.add_argument("-b", "--batch", help="Number of generated fasta files (Optional, default: 1)",
                        type=int, default=1)
    parser.add_argument("-p", "--processor", help="Number of processor cores allowed to be used by RandSeqGen "
                                                  f"(Optional, default: {DEFAULT_ALLOCATED_CPU_CORES}). "
                                                  "Please use the number of logical processor cores.",
                        type=int, default=DEFAULT_ALLOCATED_CPU_CORES)
    parser.add_argument("-o", "--output", help="Output directory path "
                                               f"(Optional, default: {DEFAULT_OUTPUT_DIR_ABS_PATH})",
                        default=DEFAULT_OUTPUT_DIR_ABS_PATH)

    parser.add_argument("-r", "--reference", help="Path to reference library directory (Optional)."
                                                  "If not specified, the program will generate complete "
                                                  "random sequence.", default=None)
    parser.add_argument("--ratio", help="Ratio of random bases in reference mode "
                                        f"(Optional, default: {DEFAULT_RANDOM_BASE_RATIO})",
                        type=float, default=DEFAULT_RANDOM_BASE_RATIO)

    parser.add_argument("--verbose", help="Verbose mode (Optional). Will show more execution details.",
                        action="store_true")

    parsed_args = parser.parse_args()

    length = parsed_args.length
    number = parsed_args.number
    batch = parsed_args.batch
    cpu_cores = parsed_args.processor
    output_dir_path = parsed_args.output
    ref_lib = parsed_args.reference
    random_ratio = parsed_args.ratio
    flag_verbose = parsed_args.verbose

    SeqGenerator_instance = SeqGenerator(length=length, number=number, batch=batch, cpu_cores=cpu_cores,
                                         output_dir=output_dir_path, ref_lib=ref_lib, random_ratio=random_ratio,
                                         flag_verbose=flag_verbose)

    SeqGenerator_instance.execute()
