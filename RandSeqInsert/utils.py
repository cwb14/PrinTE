#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)

import os
import random
import re
from typing import Tuple, List, Dict, Iterable, Sequence, Union, Optional, Callable
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def convert_humanized_number(text: str, base: Union[float, int], units: Sequence[Union[str, Iterable[str]]]) -> float:
    """
    Convert human-readable numbers with unit suffixes to numerical values.

    Parses strings containing numbers with optional unit suffixes and converts them
    to raw numerical values based on the provided base and unit definitions.

    Args:
        text: Input string containing the number and optional unit suffix
        base: Numerical base for unit conversion (e.g., 1000 or 1024)
        units: Sequence of unit suffixes ordered by increasing scale

    Returns:
        float: The converted numerical value

    Examples:
        >>> convert_humanized_number("1.14 kbp", 1000, ["bp", ("kbp", "kb"), ("Mbp", "Mb")])
        1140.0
        >>> convert_humanized_number("5.14", 1024, ('B', "KiB", "MiB", "GiB"))
        5.14
        >>> convert_humanized_number("1.919e-3M", 1000, ("", 'K', 'M', 'B', 'T'))
        1919.0
    """
    match = re.fullmatch(
        r"\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*([a-zA-Z]*)",
        text.strip()
    )
    if not match:
        raise ValueError(f"Invalid format: \"{text}\"")

    num_str, raw_unit = match.groups()

    try:
        number = float(num_str)
    except ValueError:
        raise ValueError(f"Invalid number: \"{num_str}\"")

    if not raw_unit:
        return number

    target_unit = raw_unit.lower()
    try:
        unit_index = next(index for index, entry in enumerate(units)
                          if any(unit.lower() == target_unit for unit in
                                 (entry if isinstance(entry, Iterable) and
                                           not isinstance(entry, str) else (entry,))))
    except StopIteration:
        raise ValueError(f"Unknown unit \"{raw_unit}\"")

    return number * (base ** unit_index)


def convert_humanized_int(number: Union[str, int]) -> int:
    """
    Convert a human-readable number string (e.g., '1k', '2M') to an integer.

    Args:
        number: Either an integer or a string with optional k/m/b/t suffix

    Returns:
        int: The converted integer value
    """
    if isinstance(number, int):
        return number
    return int(convert_humanized_number(number, 1000, ("", 'K', 'M', 'B', 'T')))


def sort_multiple_lists(base: list, *lists: list,
                        key: Optional[Callable] = None, reverse: bool = False) -> Union[list, tuple]:
    """
    Synchronously sorts multiple lists based on the sorting order of the base list.

    Args:
        base: The base list used as the reference for sorting
        *lists: Other lists to be sorted synchronously
        key: A function to execute to decide the order
        reverse: Whether to sort in descending order

    Returns:
        tuple: A tuple containing all sorted lists
    """
    if not base:
        raise ValueError("[ERROR] Base list cannot be empty.")
    if not lists:
        return sorted(base, key=key, reverse=reverse)

    # Combine all lists and validate length consistency
    all_lists = [base] + list(lists)
    base_len = len(base)
    len_set = set(len(l) for l in all_lists)
    if len(len_set) != 1 and not (len(len_set) == 2 and 0 in len_set):
        raise ValueError("[ERROR] All non-empty lists must have the same length as the base list.")

    # Generate sorting indices
    if key is None:
        sorted_idx = sorted(range(base_len), key=lambda i: base[i], reverse=reverse)
    else:
        sorted_idx = sorted(range(base_len), key=lambda i: key(base[i]), reverse=reverse)

    # Reorganize all lists
    sorted_lists = []
    for l in all_lists:
        if l:
            sorted_l = [l[i] for i in sorted_idx]
            sorted_lists.append(sorted_l)
        else:
            sorted_lists.append([])

    return tuple(sorted_lists)


def create_sequence_record(seq: str, id: str) -> SeqRecord:
    """
    Create a BioPython SeqRecord object from a sequence string and ID.

    Args:
        seq: The biological sequence as a string
        id: The unique identifier for the sequence

    Returns:
        SeqRecord: A SeqRecord object containing the sequence and its metadata
    """
    return SeqRecord(Seq(seq), id=id, description="")


def save_multi_fasta_from_dict(records_dict: Dict[str, List[SeqRecord]], output_dir_path: str):
    """
    Save multiple FASTA files from a dictionary of sequence records.

    Args:
        records_dict: Dictionary mapping file names to lists of SeqRecord objects
        output_dir_path: Directory path where the output files will be saved
    """
    if not records_dict:
        return

    os.makedirs(output_dir_path, exist_ok=True)
    for file_name, records in records_dict.items():
        if records:  # Skip empty record lists
            output_file_path: str = os.path.join(output_dir_path, file_name)
            SeqIO.write(records, output_file_path, "fasta")


# ----------------------------------------------------------------------------------------------------------------------
# TSD (Target Site Duplication) Functions


DEFAULT_TSD_SNP_MUTATION_RATE = 0.05
DEFAULT_TSD_INDEL_MUTATION_RATE = 0.05


def add_snp_mutation(seq: str) -> str:
    """Add a single SNP mutation at a random position in the sequence.

    Args:
        seq (str): Input DNA sequence

    Returns:
        str: Mutated sequence (maintaining the same length)
    """
    if not seq:
        return seq

    # Randomly choose mutation position
    pos = random.randrange(len(seq))
    curr_base = seq[pos]
    # Choose a different base
    new_base = random.choice(list(set("ATGC") - {curr_base}))

    # Construct mutated sequence
    return seq[:pos] + new_base + seq[pos+1:]


def add_indel_mutation(seq: str) -> str:
    """Add a single indel mutation at a random position in the sequence.

    Args:
        seq (str): Input DNA sequence

    Returns:
        str: Mutated sequence (maintaining the same length)
    """
    if len(seq) < 2:
        return seq

    # 50% probability of insertion or deletion
    is_insertion = random.random() < 0.5
    new_base = random.choice("ATGC")

    if is_insertion:
        # Choose insertion position (n + 1 positions)
        pos = random.randrange(len(seq) + 1)
        return seq[:pos] + new_base + seq[pos:]
    else:
        # Choose deletion position (n positions)
        pos = random.randrange(len(seq))
        return seq[:pos] + seq[pos+1:]


def generate_TSD(seq_slice_right: str,
                 length: Optional[int] = None,
                 snp_mutation_rate: float = DEFAULT_TSD_SNP_MUTATION_RATE,
                 indel_mutation_rate: float = DEFAULT_TSD_INDEL_MUTATION_RATE) -> Tuple[str, str]:
    """Generate TSD sequences with potential mutations.

    Args:
        seq_slice_right (str): sequence fragment on right side (3' direction)
        length (int): TSD length (defaults to seq_slice_right length)
        snp_mutation_rate (float): SNP mutation probability
        indel_mutation_rate (float): InDel mutation probability

    Returns:
        Tuple[str, str]: (5' TSD sequence, 3' TSD sequence)
    """
    tsd_length = min(len(seq_slice_right), length or len(seq_slice_right))
    # if tsd_length < length:
    #     print(f"[Warning] Requested TSD length {length} exceeds available sequence length {len(seq_slice_right)}. "
    #           f"Using {tsd_length} instead.")

    # Initialize both ends of TSD
    tsd_5 = tsd_3 = seq_slice_right[:tsd_length]

    # Apply SNP mutation
    if random.random() < snp_mutation_rate:
        if random.random() < 0.5:
            tsd_5 = add_snp_mutation(tsd_5)
        else:
            tsd_3 = add_snp_mutation(tsd_3)

    # Apply InDel mutation
    if random.random() < indel_mutation_rate:
        if random.random() < 0.5:
            tsd_5 = add_indel_mutation(tsd_5)
        else:
            tsd_3 = add_indel_mutation(tsd_3)

    return tsd_5, tsd_3

# ----------------------------------------------------------------------------------------------------------------------
