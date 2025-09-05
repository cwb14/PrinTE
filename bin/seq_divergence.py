#!/usr/bin/env python3

import argparse
import subprocess
import sys
import math
import tempfile
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""
        Estimate sequence divergence with JC69 and K2P models for LTR retrotransposons.

        The script processes a multi-sequence FASTA file, performing self-alignment
        on sequences containing 'LTR' in their headers to estimate divergence of terminal repeats.
        """
    )
    parser.add_argument(
        "-i", "--input", required=True, type=str,
        help="Input multi-sequence FASTA file."
    )
    parser.add_argument(
        "-o", "--output", required=True, type=str,
        help="Output file to store divergence results."
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of parallel threads to use. Default is 1."
    )
    parser.add_argument(
        "-f", "--fasta_path", type=str, default="",
        help="Path to blastn executable if not in PATH. Default is to use blastn from PATH."
    )
    parser.add_argument(
        "--min_align", type=int, default=100,
        help="Minimum alignment length. Default is 100 bp."
    )
    parser.add_argument(
        "--max_off", type=int, default=20,
        help="Maximum terminal offset for LTR-LTR alignments. Default is 20 bp."
    )
    parser.add_argument(
        "--miu", type=float, default=3e-8,
        help="Mutation rate per bp per year. Default is 3e-8."
    )
    parser.add_argument(
        "--blast_outfmt", type=str, default="6 qseqid sseqid sstart send slen qstart qend qlen length nident btop",
        help="BLAST output format. Default is '6 qseqid sseqid sstart send slen qstart qend qlen length nident btop'."
    )
    return parser.parse_args()

def run_blastn(blast_path: str, query_seq: SeqRecord, subject_seq: SeqRecord, outfmt: str) -> List[str]:
    """
    Run blastn with given query and subject sequences.
    Returns the BLAST output as a list of lines.
    """
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as query_fasta, \
         tempfile.NamedTemporaryFile(mode='w+', delete=False) as subject_fasta:
        SeqIO.write(query_seq, query_fasta, "fasta")
        SeqIO.write(subject_seq, subject_fasta, "fasta")
        query_fasta_path = query_fasta.name
        subject_fasta_path = subject_fasta.name

    blast_cmd = [
        blast_path,
        "-query", query_fasta_path,
        "-subject", subject_fasta_path,
        "-dust", "no",
        "-outfmt", outfmt
    ]

    try:
        result = subprocess.run(
            blast_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        blast_output = result.stdout.strip().split('\n') if result.stdout.strip() else []
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST for {query_seq.id}: {e.stderr}", file=sys.stderr)
        blast_output = []
    finally:
        Path(query_fasta_path).unlink(missing_ok=True)
        Path(subject_fasta_path).unlink(missing_ok=True)

    return blast_output

def parse_btop(btop: str) -> str:
    """
    Remove digits from btop string.
    """
    import re
    return re.sub(r'\d+', '', btop)

def calculate_divergence(n_transition: int, n_transversion: int, nident: int, qlen: int, miu: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate divergence metrics.
    Returns raw_d, K2P_d, JC69_d, raw_T, K2P_T, JC69_T
    """
    tot_len = n_transition + n_transversion + nident
    if tot_len == 0:
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    raw_d = (n_transition + n_transversion) / tot_len
    if raw_d < 0.66:
        JC69_d = -3/4 * math.log(1 - 4 * raw_d / 3)
    else:
        JC69_d = raw_d
    P = n_transition / tot_len
    Q = n_transversion / tot_len
    try:
        K2P_d = -0.5 * math.log((1 - 2*P - Q) * math.sqrt(1 - 2*Q))
    except ValueError:
        K2P_d = raw_d  # Fallback if math domain error
    raw_T = raw_d / (2 * miu)
    JC69_T = JC69_d / (2 * miu)
    K2P_T = K2P_d / (2 * miu)
    return (raw_d, K2P_d, JC69_d, raw_T, K2P_T, JC69_T)

def process_sequence(
    record: SeqRecord,
    blast_path: str,
    min_align: int,
    max_off: int,
    miu: float,
    blast_outfmt: str
) -> str:
    """
    Process a single LTR sequence:
    - Perform self-BLAST
    - Parse BLAST results
    - Calculate divergence
    Returns a formatted result string.
    """
    qseqid = record.id
    seq_str = str(record.seq).upper()

    # Create a SeqRecord for subject (self)
    subject_record = SeqRecord(Seq.Seq(seq_str), id=qseqid, description="")

    # Run BLAST
    blast_results = run_blastn(blast_path, record, subject_record, blast_outfmt)

    if not blast_results:
        return ""

    # Initialize accumulators
    sum_tot_len = 0
    sum_raw_d = 0.0
    sum_K2P_d = 0.0
    sum_JC69_d = 0.0
    sum_raw_T = 0.0
    sum_K2P_T = 0.0
    sum_JC69_T = 0.0
    total_n_transition = 0
    total_n_transversion = 0
    hit_count = 0

    for line in blast_results:
        fields = line.strip().split('\t')
        if len(fields) < 11:
            continue  # Skip incomplete lines
        (
            qseqid_b, sseqid, sstart, send, slen,
            qstart, qend, qlen, length, nident, btop
        ) = fields[:11]

        if int(length) < min_align:
            continue

        # Ensure the hit is within 100bp from start or end
        if not (int(sstart) < 100 or (int(slen) - int(send)) < 100):
            continue

        # Self-alignment specific filters
        if qseqid_b == sseqid:
            if int(nident) == int(qlen):
                continue  # Skip perfect self-alignment
            # Check terminal offsets
            if not (
                (int(sstart) < max_off and (int(qlen) - int(qend)) < max_off) or
                (int(qstart) < max_off and (int(slen) - int(send)) < max_off)
            ):
                continue

        # Passed all filters
        hit_count += 1

        # Process btop
        btop_clean = parse_btop(btop)
        if len(btop_clean) % 2 != 0:
            continue  # Expect even length

        n_transition = 0
        n_transversion = 0
        n_indel = 0

        for i in range(0, len(btop_clean), 2):
            snp = btop_clean[i:i+2]
            if '-' in snp:
                n_indel += 1
            elif snp in ['AG', 'GA', 'CT', 'TC']:
                n_transition += 1
            elif snp in ['AC', 'CA', 'AT', 'TA', 'GC', 'CG', 'GT', 'TG']:
                n_transversion += 1

        # Calculate divergence
        raw_d, K2P_d, JC69_d, raw_T, K2P_T, JC69_T = calculate_divergence(
            n_transition, n_transversion, int(nident), int(qlen), miu
        )

        # Accumulate
        sum_tot_len += (n_transition + n_transversion + int(nident))
        sum_raw_d += raw_d
        sum_K2P_d += K2P_d
        sum_JC69_d += JC69_d
        sum_raw_T += raw_T
        sum_K2P_T += K2P_T
        sum_JC69_T += JC69_T
        total_n_transition += n_transition
        total_n_transversion += n_transversion

    if hit_count == 0:
        return ""

    # Compute averages
    avg_tot_len = sum_tot_len / hit_count
    avg_raw_d = sum_raw_d / hit_count
    avg_K2P_d = sum_K2P_d / hit_count
    avg_JC69_d = sum_JC69_d / hit_count
    avg_raw_T = sum_raw_T / hit_count
    avg_K2P_T = sum_K2P_T / hit_count
    avg_JC69_T = sum_JC69_T / hit_count

    # Format the result
    result_line = (
        f"{qseqid}\t{sseqid}\t{qlen}\t{slen}\t"
        f"{int(avg_tot_len)}\t{avg_raw_d:.4f}\t{avg_K2P_d:.4f}\t{avg_JC69_d:.4f}\t"
        f"{int(avg_raw_T)}\t{int(avg_K2P_T)}\t{int(avg_JC69_T)}\t"
        f"{total_n_transition}\t{total_n_transversion}\t{hit_count}\n"
    )

    return result_line

def main():
    args = parse_arguments()

    # Determine blastn path
    if args.fasta_path:
        blastn_path = Path(args.fasta_path) / "blastn"
    else:
        blastn_path = "blastn"  # Assume blastn is in PATH

    # Verify blastn is accessible
    try:
        subprocess.run([blastn_path, "-version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except Exception as e:
        print(f"Error: blastn not found or not executable. {e}", file=sys.stderr)
        sys.exit(1)

    # Read and filter sequences
    try:
        records = list(SeqIO.parse(args.input, "fasta"))
    except Exception as e:
        print(f"Error reading input FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    ltr_records = [record for record in records if 'LTR' in record.description]
    if not ltr_records:
        print("No sequences with 'LTR' found in the input FASTA.", file=sys.stderr)
        sys.exit(1)

    # Prepare for multiprocessing
    process_func = partial(
        process_sequence,
        blast_path=str(blastn_path),
        min_align=args.min_align,
        max_off=args.max_off,
        miu=args.miu,
        blast_outfmt=args.blast_outfmt
    )

    # Run multiprocessing pool
    with Pool(processes=args.threads) as pool:
        results = pool.map(process_func, ltr_records)

    # Write output
    header = (
        "qseqid\tsseqid\tqlen\tslen\ttot_len\traw_d\tK2P_d\tJC69_d\t"
        "raw_T\tK2P_T\tJC69_T\tn_transition\tn_transversion\tblast_count\n"
    )
    with open(args.output, 'w') as out_f:
        out_f.write(header)
        for res in results:
            if res:
                out_f.write(res)

    print(f"Divergence estimation completed. Results saved to {args.output}")

if __name__ == "__main__":
    main()
