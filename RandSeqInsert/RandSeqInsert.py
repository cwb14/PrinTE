#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)

"""
Random Sequence Insertion Generator

This program takes acceptor sequences and performs random insertions using sequences
from a donor library. It supports multiprocessing for efficient sequence generation.

Note: This program uses 1-based indexing for all sequence positions (both input and output).

Author: Tianyu Lu (tianyu@lu.fm)
Date: 2024-11-27
"""

import os
import argparse
import random
from typing import Tuple, List, Union, Optional, Iterable
import multiprocessing as mp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time

from utils import convert_humanized_int, create_sequence_record, sort_multiple_lists, save_multi_fasta_from_dict
from core import SequenceTree

VERSION = "v1.1.0"
INFO = ("by Tianyu (Sky) Lu (tianyu@lu.fm) "
        "released under GPLv3")

PROGRAM_ROOT_DIR_ABS_PATH = os.path.dirname(__file__)

# PRE-DEFINED PARAMETERS
VISUAL_DIR_NAME = "visualization"
DEFAULT_ALLOCATED_CPU_CORES = os.cpu_count() - 2 if os.cpu_count() > 2 else 1

DEFAULT_OUTPUT_DIR_REL_PATH = "RandSeqInsert-Result"
DEFAULT_OUTPUT_DIR_ABS_PATH = os.path.join(os.getcwd(), DEFAULT_OUTPUT_DIR_REL_PATH)


def load_sequences(path_list: Optional[List[str]], len_limit: Optional[int] = None,
                   flag_filter_n: bool = False) -> Optional[Tuple[List[str], List[str]]]:
    """
    Load donor sequences from donor library files.
    Only loads sequences that meet the length criteria into memory.
    Sequences are stored as strings for efficient processing.

    Args:
        path_list: List of paths to donor sequence files
        len_limit: Optional maximum length limit for sequences to load
        flag_filter_n: Flag to filter out sequences containing N

    Returns:
        Optional[Tuple[List[str], List[str]]]: Tuple of (sequences, sequence_ids) as strings, or None if no donor paths provided
    """
    if not path_list:
        return None

    donor_sequences: List[str] = []
    donor_ids: List[str] = []
    
    for file_path in path_list:
        for record in SeqIO.parse(file_path, "fasta"):
            # Check length and N filter criteria
            if (len(record) <= len_limit if len_limit else True) and \
               ("N" not in record.seq.upper() if flag_filter_n else True):
                donor_sequences.append(str(record.seq.upper()))
                donor_ids.append(record.id)

    # Sort by length, keeping sequences and IDs synchronized
    sorted_lengths = [len(seq) for seq in donor_sequences]
    sorted_lengths, donor_sequences, donor_ids = sort_multiple_lists(
        sorted_lengths, donor_sequences, donor_ids)
    
    return donor_sequences, donor_ids

def _find_donor_lib_abs_path_list(path: Optional[str]) -> Optional[List[str]]:
    """
    Locate the donor library files.

    Args:
        path: The path to search for donor library

    Returns:
        Optional[List[str]]: A list of absolute paths to the found donor library files,
                           or None if path is None
    """
    if not path:
        return None

    # Check if path exists and is a file
    if os.path.exists(path) and os.path.isfile(path):
        return [os.path.abspath(path)]

    # Check if path exists and is a directory
    if os.path.exists(path) and os.path.isdir(path):
        donor_lib_files: List[str] = []
        for file in os.listdir(path):
            if file.endswith((".fa", ".fasta")):
                donor_lib_files.append(os.path.join(path, file))
        if donor_lib_files:
            return list(map(os.path.abspath, donor_lib_files))

    # Check in built-in DonorLib directory
    default_donor_lib_abs_path = os.path.join(PROGRAM_ROOT_DIR_ABS_PATH, path)
    if os.path.exists(default_donor_lib_abs_path):
        if os.path.isfile(default_donor_lib_abs_path):
            return [default_donor_lib_abs_path]
        elif os.path.isdir(default_donor_lib_abs_path):
            donor_lib_files: List[str] = []
            for file in os.listdir(default_donor_lib_abs_path):
                if file.endswith((".fa", ".fasta")):
                    donor_lib_files.append(os.path.join(default_donor_lib_abs_path, file))
            if donor_lib_files:
                return donor_lib_files

    raise SystemExit(f"[ERROR] Specified donor library not found at {path} or built-in donor library!")

def _load_multiple_donor_libs(path_list: List[str], weight_list: Optional[List[float]] = None,
                            len_limit: Optional[int] = None,
                            flag_filter_n: bool = False) -> Optional[Tuple[List[str], List[str], List[int], List[float]]]:
    """
    Load multiple donor libraries by combining their sequences.

    Args:
        path_list: List of paths to donor libraries
        weight_list: List of weights for each donor library
        len_limit: Maximum length limit for sequences to load
        flag_filter_n: Flag to filter out sequences containing N

    Returns:
        Tuple containing:
            - Combined list of donor sequences from all libraries
            - Combined list of donor sequence IDs from all libraries
            - List of lengths of donor sequences
            - List of weights for each donor sequence
    """
    if not path_list:
        return None

    all_donor_sequence_list: List[str] = []
    all_donor_id_list: List[str] = []
    all_donor_len_list: List[int] = []
    all_donor_weight_list: List[float] = []

    for lib_path, weight in zip(path_list, weight_list or [None] * len(path_list)):
        single_donor_lib_abs_path_list = _find_donor_lib_abs_path_list(lib_path)
        if single_donor_lib_abs_path_list:
            single_donor_lib_result = load_sequences(single_donor_lib_abs_path_list, len_limit, flag_filter_n)
            if single_donor_lib_result:
                single_donor_lib_sequences, single_donor_lib_ids = single_donor_lib_result
                all_donor_sequence_list.extend(single_donor_lib_sequences)
                all_donor_id_list.extend(single_donor_lib_ids)
                all_donor_len_list.extend([len(seq) for seq in single_donor_lib_sequences])
                if weight:
                    len_single_donor_lib_sequences = len(single_donor_lib_sequences)
                    all_donor_weight_list.extend([weight / len_single_donor_lib_sequences] * len_single_donor_lib_sequences)

    # Sort by length, keeping all lists synchronized
    all_donor_len_list, all_donor_sequence_list, all_donor_id_list, all_donor_weight_list = sort_multiple_lists(
        all_donor_len_list, all_donor_sequence_list, all_donor_id_list, all_donor_weight_list)

    # Uniform weights by default: If no weight provided, set all weights to 1
    if not weight_list:
        all_donor_weight_list = [1] * len(all_donor_sequence_list)

    return all_donor_sequence_list, all_donor_id_list, all_donor_len_list, all_donor_weight_list

# ======================================================================================================================
# Main Sequence Generator Class

class SeqGenerator:
    # noinspection PyShadowingBuiltins
    def __init__(self, input: str, insertion: Optional[Union[str, int]] = None, batch: int = 1,
                 processors: int = DEFAULT_ALLOCATED_CPU_CORES, output_dir_path: str = None,
                 donor_lib: Optional[List[str]] = None, donor_lib_weight: Optional[List[float]] = None,
                 donor_len_limit: Optional[int] = None, flag_filter_n: bool = False, flag_track: bool = False,
                 tsd_length: Optional[int] = None, flag_visual: bool = False, flag_recursive: bool = False,
                 iteration: int = 1, flag_debug: bool = False, seed: Optional[int] = None):
        """
        Initialize the sequence generator.

        Args:
            input: Path to the input sequence file
            insertion: Number of insertions per sequence
            batch: Number of independent result files to generate
            processors: Number of processors to use
            output_dir_path: Directory to save output files
            donor_lib: List of donor library file paths
            donor_lib_weight: Weights for donor libraries
            donor_len_limit: Maximum length limit for donor sequences to load
            flag_filter_n: Whether to filter out sequences containing N
            flag_track: Whether to track donor sequences used
            tsd_length: Length of Target Site Duplication (TSD) to generate
            flag_visual: Whether to generate graphviz visualization of the sequence tree
            flag_recursive: Whether to use recursive insertion method
            iteration: Number of insertion iterations to perform on each sequence
            flag_debug: Whether to enable debug mode
            seed: Random seed for reproducibility
        """
        self.input = input
        self.insertion: int = convert_humanized_int(insertion)
        self.batch: int = batch
        self.processors: int = processors
        self.output_dir_path: str = output_dir_path
        self.donor_len_limit: Optional[int] = donor_len_limit
        self.flag_filter_n: bool = flag_filter_n
        self.flag_track: bool = flag_track
        self.tsd_length: Optional[int] = tsd_length
        self.flag_visual: bool = flag_visual
        self.flag_recursive: bool = flag_recursive
        self.iteration: int = iteration
        self.flag_debug: bool = flag_debug
        self.seed: Optional[int] = seed

        # Load input sequence
        try:
            self.input = list(SeqIO.parse(input, "fasta"))
            if not self.input:
                raise ValueError(f"No sequences found in input file {input}")
        except (FileNotFoundError, IOError) as e:
            raise ValueError(f"Error loading input file {input}: {str(e)}")

        # Load donor sequences
        if not donor_lib:
            raise ValueError("Donor library is required for insertion. Please provide at least one donor library path using the --donor parameter.")

        try:
            donor_lib_data = _load_multiple_donor_libs(donor_lib, donor_lib_weight, donor_len_limit, flag_filter_n)
            if not donor_lib_data or not donor_lib_data[0]:
                raise ValueError("No valid donor sequences loaded. Check if donor files exist and contain valid sequences.")
            # Assign donor sequences and weights, ignore lengths as they're not needed
            self.donor_sequences, self.donor_ids, _, self.donor_weights = donor_lib_data
            # Normalize weights if they were provided
            if donor_lib_weight and sum(self.donor_weights) > 0:
                weight_sum = sum(self.donor_weights)
                self.donor_weights = [w / weight_sum for w in self.donor_weights]
        except Exception as e:
            # Convert any exception during donor loading to a more informative error
            raise ValueError(f"Error loading donor libraries: {str(e)}")

    def __print_header(self):
        """
        Print the program header with basic information.
        """
        print(f"=== RandSeqInsert {VERSION} ===")
        print(f"Number of input sequences: {len(self.input)}")
        print(f"Insertion settings: {self.insertion} insertions per sequence")
        if self.iteration > 1:
            print(f"Performing {self.iteration} iterations of insertion")
        if self.tsd_length:
            print(f"TSD settings: Generating TSD of length {self.tsd_length} at insertion sites")
        print(f"Donor library: {len(self.donor_sequences)} sequences loaded")
        print(f"Generating {self.batch} independent result file(s)")
        if self.seed is not None:
            print(f"Seed: {self.seed}")
        if self.flag_recursive:
            print(f"Using recursive insertion method")
        if self.flag_visual:
            print(f"Tree and Graph visualization enabled")
        if self.flag_debug:
            print(f"Debug mode enabled: detailed nesting graph information will be printed")

    def __pre_check(self):
        """
        Perform pre-execution checks.
        """
        if bool(self.insertion) ^ bool(self.donor_sequences):
            raise ValueError("\"insertion\" and \"donor_sequences\" must be specified or not specified at the same time")

    def __process_single_batch_multiprocessing(self, batch_num, seed: int = None) -> Optional[float]:
        """
        Process a single batch of sequences using multiprocessing.

        Args:
            batch_num: Batch number (starting from 1)
            seed: Optional seed for this batch

        Returns:
            float: Time taken to process this batch (seconds)
        """
        if not self.input:
            return None

        batch_start_time = time.time()

        total_sequences = len(self.input)
        if seed:
            mp_args_gen = ((i, seed + i) for i in range(total_sequences))
        else:
            mp_args_gen = list(range(total_sequences))

        print(f"Using multiprocessing with {self.processors} worker processes to process {total_sequences} sequences")

        all_sequences = []
        all_donors = [] if self.flag_track else None
        all_reconstructed_donors = [] if self.flag_track else None
        all_seq_trees = []

        with mp.Pool(self.processors) as pool:
            results = pool.imap_unordered(self._imap_worker_process_single_sequence, mp_args_gen, chunksize=1)

            for i, (processed_seq, donors, reconstructed, seq_tree) in enumerate(results, 1):
                all_sequences.append(processed_seq)
                all_seq_trees.append((processed_seq.id, seq_tree))
                if self.flag_track:
                    if donors:
                        all_donors.extend(donors)
                    if reconstructed:
                        all_reconstructed_donors.extend(reconstructed)

                if i % 10 == 0 or i == total_sequences:
                    print(f"Completed {i}/{total_sequences} sequences")

        os.makedirs(self.output_dir_path, exist_ok=True)
        print(f"Output directory: {self.output_dir_path}")
        self.__save_batch_results(self.output_dir_path, all_sequences, all_donors, all_reconstructed_donors, batch_num, all_seq_trees)

        return time.time() - batch_start_time

    def _imap_worker_process_single_sequence(self, args):
        if isinstance(args, Iterable):
            return self.__process_single_sequence(*args)
        return self.__process_single_sequence(args)

    def __process_single_sequence(self, idx_or_record, seed = None) -> Tuple[SeqRecord, Optional[List[SeqRecord]], Optional[List[SeqRecord]], Optional['SequenceTree']]:
        """
        Process a single sequence with random insertions, potentially over multiple iterations.

        Args:
            idx_or_record: Either an index to look up in self.input, or a SeqRecord to process directly

        Returns:
            Tuple containing:
                - The processed sequence
                - Used donors (if tracking enabled)
                - Reconstructed donors (if tracking enabled)
                - The SequenceTree object
        """
        # If idx_or_record is an integer, retrieve the sequence record
        # If it's already a SeqRecord, use it directly
        if isinstance(idx_or_record, int):
            seq_record = self.input[idx_or_record]
        else:
            seq_record = idx_or_record

        # Check if there are donor sequences to insert
        if not self.insertion or not self.donor_sequences:
            # If no insertions requested or no donor sequences available, return the original sequence
            return seq_record, None, None, None

        # Create a new sequence tree for this sequence
        seq_tree = SequenceTree(str(seq_record.seq), 0)
        # Track all used donors across iterations if tracking is enabled
        all_used_donors = [] if self.flag_track else None
        all_reconstructed_donors = [] if self.flag_track else None

        # Set the seed for random
        if not seed:
            seed = time.perf_counter_ns()
        random.seed(seed)
        if self.flag_debug:
            print(f"[DEBUG] Seed for {seq_record.id}: [{seed}]")

        # Perform multiple iterations of insertion
        for iteration in range(1, self.iteration + 1):
            # Get the current total sequence length (updated after each iteration)
            current_seq = str(seq_tree)
            total_length = len(current_seq)

            # Generate insertion positions and donor sequences
            # Use 1-based positions (1 to total_length+1)
            insert_positions = random.choices(range(1, total_length + 2), k=self.insertion)
            
            # Select donor indices to ensure sequence and ID correspondence
            selected_donor_indices = random.choices(range(len(self.donor_sequences)), weights=self.donor_weights, k=self.insertion)
            selected_donors = [self.donor_sequences[i] for i in selected_donor_indices]
            selected_donor_ids = [self.donor_ids[i] for i in selected_donor_indices]

            # Sort positions in descending order to maintain position integrity during insertion
            # Also sort donor sequences and IDs accordingly
            insert_positions, selected_donors, selected_donor_ids = sort_multiple_lists(
                insert_positions, selected_donors, selected_donor_ids, reverse=True)

            # Insert donors into the tree with their real IDs
            for pos, donor_seq, donor_id in zip(insert_positions, selected_donors, selected_donor_ids):
                seq_tree.insert(pos, donor_seq, donor_id, self.tsd_length, self.flag_recursive, self.flag_debug)

        # Collect donor sequences once after all iterations if tracking is enabled
        if self.flag_track:
            used_donors, reconstructed_donors = seq_tree.donors(seq_record.id)
            all_used_donors = used_donors if used_donors else []
            all_reconstructed_donors = reconstructed_donors if reconstructed_donors else []

        # Generate final sequence
        new_id = f"{seq_record.id}_ins{self.insertion}"
        if self.iteration > 1:
            new_id += f"_iter{self.iteration}"
        if self.tsd_length:
            new_id += f"_tsd{self.tsd_length}"
        new_seq_record = create_sequence_record(str(seq_tree), new_id)

        # Generate visualization if requested
        if self.flag_visual:
            # Generate tree visualization
            graphviz_str = seq_tree.to_graphviz_dot()
            tree_visual_dir_path = os.path.join(self.output_dir_path, VISUAL_DIR_NAME)
            os.makedirs(tree_visual_dir_path, exist_ok=True)
            with open(os.path.join(tree_visual_dir_path, f"{seq_record.id}_tree_visual.dot"), "w") as f:
                f.write(graphviz_str)

            # Generate nesting graph visualization
            graph_visual_dir_path = os.path.join(self.output_dir_path, VISUAL_DIR_NAME)
            os.makedirs(graph_visual_dir_path, exist_ok=True)
            nesting_graphviz_str = seq_tree.event_journal.to_graphviz_dot()
            with open(os.path.join(graph_visual_dir_path, f"{seq_record.id}_event_visual.dot"), "w") as f:
                f.write(nesting_graphviz_str)

            print(f"Generated tree and graph visualizations for {seq_record.id}")

        # If debug mode is enabled, print the tree structure details
        if self.flag_debug:
            print(f"\n[DEBUG] Tree structure information for {seq_record.id}:")
            print(f"Total length: {len(str(seq_tree))}")
            print(f"Node count: {len(seq_tree.collect_active_nodes())}")
            print("[DEBUG] End of tree structure information\n")

        return new_seq_record, all_used_donors, all_reconstructed_donors, seq_tree

    def execute(self):
        """
        Execute the sequence generation process.
        """
        self.__print_header()

        start_time = time.time()
        self.__pre_check()

        # Process each batch
        for batch_num in range(1, self.batch + 1):
            print(f"Processing batch {batch_num}/{self.batch}")
            seed = (self.seed + batch_num) if self.seed else None
            batch_elapsed_time = self.__process_single_batch_multiprocessing(batch_num, seed)
            print(f"Batch {batch_num} completed in {batch_elapsed_time:.2g} seconds")

        # Print summary
        self.__print_summary(time.time() - start_time)

    def __print_summary(self, total_elapsed_time: float):
        """
        Print a summary of the processing results.

        Args:
            total_elapsed_time: Total time taken for all batches
        """
        print(f"\nAll batches completed in {total_elapsed_time:.2g} seconds")
        print(f"Results saved to \"{os.path.abspath(self.output_dir_path)}\"")

    def _generate_bed_records_for_all_donors(self, chrom: str, seq_tree: 'SequenceTree') -> List[str]:
        """
        Generate BED format records for all donor insertions.
        Uses tree structure directly to calculate positions and get TSD from nodes.
        
        Args:
            chrom: Chromosome/sequence ID
            seq_tree: SequenceTree object containing insertion information
            
        Returns:
            List[str]: List of BED format records
        """
        bed_records = []
        
        def _collect_donors_with_positions(node, current_pos=0):
            """Traverse tree and collect donor nodes with their positions."""
            if not node:
                return current_pos
            
            # Process left subtree
            left_end_pos = _collect_donors_with_positions(node.left, current_pos)
            
            # Process current node
            node_start_pos = left_end_pos
            if node.is_donor:
                # Calculate clean sequence position (excluding TSD)
                clean_start = node_start_pos + len(node.tsd_5)
                clean_length = node.length
                
                # Skip empty donors
                if clean_length > 0:
                    clean_end = clean_start + clean_length
                    
                    # Get TSD sequence for BED (use 5' TSD, or "NA" if none)
                    tsd_seq = node.tsd_5 or "NA"
                    
                    # Build name field
                    name_parts = []
                    if node.donor_id:
                        name_parts.append(node.donor_id)
                    else:
                        name_parts.append(f"donor_{node.uid}")
                    
                    # Join name parts with semicolons
                    name = ";".join(name_parts)
                    
                    # Create BED record: chrom start end name tsd strand
                    bed_record = f"{chrom}\t{clean_start}\t{clean_end}\t{name}\t{tsd_seq}\t+"
                    bed_records.append(bed_record)
            
            # Process right subtree
            right_end_pos = _collect_donors_with_positions(node.right, node_start_pos + node.length + node.get_tsd_length())
            
            return right_end_pos
        
        _collect_donors_with_positions(seq_tree.root)
        return bed_records

    def __save_batch_results(self, output_dir: str, sequences: List[SeqRecord], 
                             donors: Optional[List[SeqRecord]], 
                             reconstructed_donors: Optional[List[SeqRecord]],
                             batch_num: int = 1, seq_trees: List[Tuple[str, 'SequenceTree']] = None):
        """
        Save the batch results to output files.

        Args:
            output_dir: Directory to save the results
            sequences: List of sequence records to save
            donors: List of donor records to save if tracking is enabled
            reconstructed_donors: List of reconstructed donor records
            batch_num: The batch number (1-based index)
            seq_trees: List of tuples containing sequence ID and SequenceTree object
        """
        if not sequences:
            print("No sequences to save")
            return

        # Add batch suffix to filenames if multiple batches
        suffix = f"_batch{batch_num}" if self.batch > 1 else ""

        print(f"Saving {len(sequences)} processed sequences")
        output_dict = {f"sequences{suffix}.fa": sequences}

        if self.flag_track:
            if donors:
                print(f"Saving {len(donors)} donor sequence records")
                output_dict[f"used_donors{suffix}.fa"] = donors

            if reconstructed_donors:
                # Group by reconstruction type
                full_recs = [d for d in reconstructed_donors if "reconstruction_type" in d.annotations and d.annotations["reconstruction_type"] == "full"]
                clean_recs = [d for d in reconstructed_donors if "reconstruction_type" in d.annotations and d.annotations["reconstruction_type"] == "clean"]

                if full_recs:
                    print(f"Saving {len(full_recs)} full reconstructed donor records")
                    output_dict[f"reconstructed_donors{suffix}.fa"] = full_recs

                if clean_recs:
                    print(f"Saving {len(clean_recs)} clean reconstructed donor records")
                    output_dict[f"clean_reconstructed_donors{suffix}.fa"] = clean_recs

        # Generate BED file if tracking is enabled and we have seq_trees
        if self.flag_track and seq_trees:
            bed_records = []
            for seq_id, seq_tree in seq_trees:
                if seq_tree:
                    bed_records.extend(self._generate_bed_records_for_all_donors(seq_id, seq_tree))
            
            if bed_records:
                bed_filename = f"donors{suffix}.bed"
                bed_filepath = os.path.join(output_dir, bed_filename)
                print(f"Saving {len(bed_records)} BED records to {bed_filename}")
                with open(bed_filepath, "w") as bed_file:
                    # Write track header
                    bed_file.write(f'track name="RandSeqInsert_Donors" description="Transposon insertions generated by RandSeqInsert" itemRgb="On"\n')
                    for record in bed_records:
                        bed_file.write(record + "\n")

        save_multi_fasta_from_dict(output_dict, output_dir)

def main():
    """Main function to handle command line arguments and execute the program."""
    parser = argparse.ArgumentParser(prog="RandSeqInsert",
                                     description="RandomSequenceInsertion: A high-performance Python tool for randomly inserting "
                                               "genomic fragments from donor libraries into acceptor sequences. "
                                               "It supports multiprocessing for efficient processing of large datasets.",
                                     epilog="")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}\n{INFO}")

    # Core Arguments
    core_group = parser.add_argument_group("Core Arguments")
    core_group.add_argument("-i", "--input",
                       help="Input sequence file in FASTA format. Contains the acceptor sequences to be inserted into.",
                       type=str, required=True, metavar="FILE")
    core_group.add_argument("-o", "--output", default=DEFAULT_OUTPUT_DIR_ABS_PATH, metavar="DIR",
                       help=f"Output directory path. Generated sequences and related files will be saved here. Default: \"{DEFAULT_OUTPUT_DIR_ABS_PATH}\"")

    core_group.add_argument("-is", "--insert", metavar="INT/STR",
                       help="Number of insertions per sequence. Accepts either a plain number (e.g., 100) or a string with k/m suffix (e.g., 1k, 1m).",
                       type=str)
    core_group.add_argument("-it", "--iteration", type=int, default=1, metavar="INT",
                            help="Number of insertion iterations to perform on each sequence. Each iteration will use the sequence from the previous iteration as input. Default: 1")

    # Donor Library Arguments
    donor_group = parser.add_argument_group("Donor Library Arguments")
    donor_group.add_argument("-d", "--donor", nargs="+", metavar="FILE/DIR",
                            help="Donor sequence library file or directory paths. Multiple FASTA format donor files can be specified. Sequences from these files will be selected and inserted into the acceptor sequences.")
    donor_group.add_argument("-w", "--weight", type=float, nargs="+", metavar="FLOAT",
                       help="Weights for donor libraries. Controls the probability of selecting sequences from different donor libraries. The number of weights should match the number of donor libraries.")
    donor_group.add_argument("-l", "--limit", default=None, type=int, metavar="INT",
                       help="Donor sequence length limit. Only loads donor sequences with length less than or equal to this value. Default: no limit.")

    # Control Arguments
    ctrl_group = parser.add_argument_group("Control Arguments")
    ctrl_group.add_argument("-b", "--batch", default=1, type=int, metavar="INT",
                            help="Number of independent result files to generate. Runs the entire process multiple times with different random seeds to generate multiple output sets. Default: 1")
    ctrl_group.add_argument("-p", "--processors", default=DEFAULT_ALLOCATED_CPU_CORES, type=int, metavar="INT",
                       help=f"Number of processors to use for parallel processing. Default: {DEFAULT_ALLOCATED_CPU_CORES}")
    ctrl_group.add_argument("--seed", type=int, metavar="INT",
                       help="Random seed for reproducibility. Setting this will make insertions deterministic. Each process in multiprocessing will derive its seed from this base seed.")
    ctrl_group.add_argument("--tsd", default=0, type=int, metavar="LENGTH",
                       help="Enable Target Site Duplication (TSD) with specified length. When a donor sequence is inserted, TSD of this length will be generated at the insertion site.")

    # Flags
    flag_group = parser.add_argument_group("Flags")
    flag_group.add_argument("--filter_n", action="store_true",
                                help="Filter out donor sequences containing N.")
    flag_group.add_argument("--track", action="store_true",
                       help="Track and save used donor sequences. Enable this option to generate an additional FASTA file in the output directory recording all used donor sequences and their insertion positions.")
    flag_group.add_argument("--visual", action="store_true",
                       help="Generate Graphviz DOT files visualizing the tree structure of each sequence. Files will be named {seqid}_tree_visual.dot and saved in the output directory.")
    flag_group.add_argument("--recursive", action="store_true",
                       help="Use recursive insertion method instead of iterative one.")
    flag_group.add_argument("--debug", action="store_true",
                       help="Enable debug mode to print detailed information about the nesting graph.")

    parsed_args = parser.parse_args()

    generator = SeqGenerator(
        input=parsed_args.input,
        insertion=parsed_args.insert,
        batch=parsed_args.batch,
        processors=parsed_args.processors,
        output_dir_path=parsed_args.output,
        donor_lib=parsed_args.donor,
        donor_lib_weight=parsed_args.weight,
        donor_len_limit=parsed_args.limit,
        flag_filter_n=parsed_args.filter_n,
        flag_track=parsed_args.track,
        tsd_length=parsed_args.tsd,
        flag_visual=parsed_args.visual,
        flag_recursive=parsed_args.recursive,
        iteration=parsed_args.iteration,
        flag_debug=parsed_args.debug,
        seed=parsed_args.seed
    )
    generator.execute()

if __name__ == "__main__":
    main()
