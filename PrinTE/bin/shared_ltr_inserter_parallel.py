#!/usr/bin/env python3

import argparse
import random
import sys
import bisect
import re
import numpy as np
import math
import matplotlib.pyplot as plt
import multiprocessing

def parse_fasta(file_path):
    fasta_dict = {}
    header = None
    seq_chunks = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header:
                        fasta_dict[header] = ''.join(seq_chunks).upper()
                    header = line[1:].split()[0]
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
    with open(file_path, 'w') as f:
        for header, sequence in fasta_dict.items():
            f.write(f">{header}\n")
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')

def write_bed(file_path, bed_entries):
    with open(file_path, 'w') as f:
        for entry in bed_entries:
            chrom = entry['chromosome']
            start = entry['start']
            end = entry['end']
            name = entry['name']
            strand = entry['strand']
            tsd = entry['tsd']
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{tsd}\t{strand}\n")

def merge_intervals(intervals):
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def parse_bed(file_path):
    bed_dict = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) < 4:
                    continue
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if chrom not in bed_dict:
                    bed_dict[chrom] = []
                bed_dict[chrom].append((start, end))
        for chrom in bed_dict:
            merged = merge_intervals(bed_dict[chrom])
            bed_dict[chrom] = merged
    except FileNotFoundError:
        print(f"Error: Gene BED file {file_path} not found.")
        sys.exit(1)
    return bed_dict

def parse_gene_bed(file_path):
    genes = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) < 4:
                    continue
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                genes.append({'chromosome': chrom, 'start': start, 'end': end, 'name': name})
    except FileNotFoundError:
        print(f"Error: Gene BED file {file_path} not found.")
        sys.exit(1)
    return genes

def extract_te_info(header):
    """
    Extract TE_class and TE_superfamily from the TE header.
    Expected formats:
    '>TE_name#TE_class/TE_superfamily'
    '>TE_name#TE_class/TE_superfamily;TE_supplemental_info'
    """
    # Use regex to extract class and superfamily
    match = re.match(r"[^#]+#([^/]+)/([^~]+)", header)
    if match:
        te_class = match.group(1).strip()
        te_superfamily = match.group(2).strip()
        return te_class, te_superfamily
    else:
        # Default values if parsing fails
        return "unknown", "unknown"

def get_tsd_length(te_class, te_superfamily):
    """
    Determine the TSD length based on TE class and superfamily.
    Returns an integer representing the TSD length.
    If TSD length is variable, selects a random length within the specified range.
    """
    tsd_mapping = {
        # Fixed TSD lengths
        ('LTR', 'Copia'): lambda: 5,
        ('LTR', 'Gypsy'): lambda: 5,
        ('LTR', 'Solo'): lambda: 5,
        ('LTR', 'Ty3'): lambda: 5,
        ('LTR', 'unknown'): lambda: 5,
        ('DNA', 'Harbinger'): lambda: 3,
        ('DNA', 'Mariner'): lambda: 2,
        ('DNAauto', 'Helitron'): lambda: 0,
        ('DNA', 'Helitron'): lambda: 0,
        ('DNAnona', 'Helitron'): lambda: 0,

        # Variable TSD lengths (range inclusive)
        ('MITE', 'Stow'): lambda: random.randint(2, 10),
        ('MITE', 'Tourist'): lambda: random.randint(2, 10),
        ('SINE', 'tRNA'): lambda: random.randint(5, 20),
        ('SINE', 'unknown'): lambda: random.randint(5, 20),
        ('LINE', 'L1'): lambda: random.randint(5, 20),
        ('LINE', 'unknown'): lambda: random.randint(5, 20),
        ('DNAauto', 'CACTA'): lambda: random.randint(2, 4),
        ('DNA', 'CACTA'): lambda: random.randint(2, 4),
        ('DNAnona', 'CACTA'): lambda: random.randint(2, 4),
        ('DNA', 'hAT'): lambda: random.randint(5, 8),
        ('DNAnona', 'hAT'): lambda: random.randint(5, 8),
        ('DNAauto', 'hAT'): lambda: random.randint(5, 8),
        ('DNAauto', 'MuDR'): lambda: random.randint(8, 9),
        ('DNAnona', 'MULE'): lambda: random.randint(8, 9),
        ('DNA', 'MuDR'): lambda: random.randint(8, 9),
        ('DNAauto', 'MULE'): lambda: random.randint(8, 9),
    }

    key = (te_class, te_superfamily)
    if key in tsd_mapping:
        return tsd_mapping[key]()
    else:
        print(f"Warning: TSD length not defined for TE class '{te_class}' and superfamily '{te_superfamily}'. Using default TSD length of 5.")
        return 5

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(complement)[::-1]

def compute_allowed_intervals(seq_length, exclusion_intervals):
    """
    Given the full length of a sequence and a list of exclusion intervals (start, end),
    return a list of allowed intervals (complement of exclusions) within [0, seq_length].
    """
    allowed = []
    if not exclusion_intervals:
        return [(0, seq_length)]
    merged = merge_intervals(exclusion_intervals)
    current = 0
    for (ex_start, ex_end) in merged:
        if ex_start > current:
            allowed.append((current, ex_start))
        current = max(current, ex_end)
    if current < seq_length:
        allowed.append((current, seq_length))
    return allowed

#def mutate_sequence(seq, mutation_rate):
#    """
#    Mutate the given sequence by randomly substituting bases with probability mutation_rate.
#    mutation_rate should be a fraction (e.g., 0.03 for 3%).
#    """
#    mutated = []
#    for base in seq:
#        if random.random() < mutation_rate:
#            bases = ['A', 'C', 'G', 'T']
#            if base.upper() in bases:
#                bases.remove(base.upper())
#            mutated.append(random.choice(bases))
#        else:
#            mutated.append(base)
#    return ''.join(mutated)

def mutate_sequence(seq, mutation_rate, ts_tv_ratio):
    """
    Mutate seq by substitution at rate `mutation_rate`.  When substituting,
    choose a transition vs. transversion event with ratio ts_tv_ratio.
    Prints counts of Ts, Tv, and the observed Ts/Tv ratio.
    """
    ts_map = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
    tv_map = {
        'A': ['C', 'T'], 'G': ['C', 'T'],
        'C': ['A', 'G'], 'T': ['A', 'G']
    }

    transitioned = 0
    transverted = 0
    mutated = []

    for base in seq:
        if random.random() < mutation_rate:
            b = base.upper()
            if b in ts_map:
                # decide between Ts and Tv
                p_ts = ts_tv_ratio / (ts_tv_ratio + 1.0)
                if random.random() < p_ts:
                    new = ts_map[b]
                    transitioned += 1
                else:
                    new = random.choice(tv_map[b])
                    transverted += 1
                # preserve case
                mutated.append(new if base.isupper() else new.lower())
            else:
                # non-ACGT characters unchanged
                mutated.append(base)
        else:
            mutated.append(base)

    # After processing, report counts
    obs_ratio = (transitioned / transverted) if transverted > 0 else float('inf')
#    print(f"[mutate_sequence] Transitions: {transitioned}; Transversions: {transverted}; Observed Ts/Tv = {obs_ratio:.2f}")

    return ''.join(mutated)

def plot_decay_function(k, Mmax, pdf_out):
    """
    Generate a publication-quality PDF of the exponential decay function used for sampling mutation percentages.
    The probability density function (PDF) is:
      f(M) = (k / (Mmax * (1 - exp(-k)))) * exp(-k * (M / Mmax))
    for M in [0, Mmax].
    """
    M_values = np.linspace(0, Mmax, 500)
    # Compute normalized PDF
    norm_const = 1 - math.exp(-k)
    f_values = (k / (Mmax * norm_const)) * np.exp(-k * (M_values / Mmax))

    plt.figure(figsize=(6,4))
    plt.plot(M_values, f_values, lw=2)
    plt.xlabel("Mutation Percent (%)", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.title("Exponential Decay Function for Mutation Sampling", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(pdf_out, format='pdf')
    plt.close()
    print(f"Decay function plot saved to {pdf_out}.")

def process_chromosome(task):
    """
    Process TE insertions on a single chromosome.
    task is a tuple containing:
      (chrom, chrom_data, target_mode, target_value, te_headers, te_info_dict,
       te_dict, te_ratio_weights, te_type_to_headers, k, Mmax, seed_offset, ts_tv_ratio)
    target_mode: 'n' (number of insertions) or 'p' (target TE bp)
    target_value: integer target insertions or target bp to insert
    """
    (chrom, data, target_mode, target_value, te_headers, te_info_dict,
     te_dict, te_ratio_weights, te_type_to_headers, k, Mmax, seed_offset, ts_tv_ratio) = task

    # Optionally seed random generators per chromosome
    if seed_offset is not None:
        random.seed(seed_offset)
        np.random.seed(seed_offset)

    seq = list(data['seq'])
    chr_length = data['length']
    # Copy the exclusion intervals (gene plus already inserted TEs)
    exclusion_intervals = list(data['exclusion'])
    local_te_bed_entries = []
    total_TE_bp_inserted = 0
    total_insertions_done = 0
    attempts = 0
    max_global_attempts = 100000

    def try_insertion():
        nonlocal exclusion_intervals, seq, total_TE_bp_inserted, total_insertions_done
        max_attempts = 10
        for attempt in range(max_attempts):
            # Select a TE.
            if te_ratio_weights:
                keys = list(te_ratio_weights.keys())
                weights = [te_ratio_weights[kkey] for kkey in keys]
                chosen_key = random.choices(keys, weights=weights, k=1)[0]
                if chosen_key in te_type_to_headers and te_type_to_headers[chosen_key]:
                    selected_te = random.choice(te_type_to_headers[chosen_key])
                else:
                    # Fallback to uniform sampling
                    selected_te = random.choice(te_headers)
            else:
                selected_te = random.choice(te_headers)

            te_class = te_info_dict[selected_te]['class']
            te_superfamily = te_info_dict[selected_te]['superfamily']
            tsd_length = get_tsd_length(te_class, te_superfamily)
            te_sequence = te_dict[selected_te]

            # Sample a mutation percent using the exponential decay function.
            u = random.random()
            mutation_percent = - (Mmax / k) * math.log(1 - u * (1 - math.exp(-k)))
            mutation_rate = mutation_percent / 100.0
            te_sequence = mutate_sequence(te_sequence, mutation_rate, ts_tv_ratio)
            # Uncomment the next line to print detailed mutation info per chromosome.
            # print(f"[{chrom}] Mutating TE {selected_te} at {mutation_percent:.2f}% rate.")

            TE_length = len(te_sequence)
            deletion_length = TE_length + tsd_length  # Total backbone sequence to delete.

            # Compute allowed intervals.
            allowed = compute_allowed_intervals(chr_length, exclusion_intervals)
            candidate_intervals = []
            total_candidates = 0
            for (a, b) in allowed:
                lower_bound = a + tsd_length
                upper_bound = b - deletion_length
                if lower_bound <= upper_bound:
                    candidate_intervals.append((lower_bound, upper_bound))
                    total_candidates += (upper_bound - lower_bound + 1)
            if total_candidates <= 0:
                continue

            # Randomly choose a deletion start position.
            r = random.randint(0, total_candidates - 1)
            chosen_deletion_start = None
            for (low, high) in candidate_intervals:
                interval_length = high - low + 1
                if r < interval_length:
                    chosen_deletion_start = low + r
                    break
                r -= interval_length
            if chosen_deletion_start is None:
                print(f"[{chrom}] Error: Failed to select a deletion position.")
                sys.exit(1)
            deletion_start = chosen_deletion_start

            # Extract TSD from immediately upstream (if applicable).
            if tsd_length > 0:
                tsd_seq = ''.join(seq[deletion_start - tsd_length: deletion_start])
            else:
                tsd_seq = ''

            strand = random.choice(['+', '-'])
            if strand == '-':
                te_seq_final = reverse_complement(te_sequence)
            else:
                te_seq_final = te_sequence

            replacement = te_seq_final + tsd_seq
            seq[deletion_start : deletion_start + deletion_length] = list(replacement)

            local_te_bed_entries.append({
                'chromosome': chrom,
                'start': deletion_start,
                'end': deletion_start + len(te_seq_final),
                'name': selected_te,
                'strand': strand,
                'tsd': tsd_seq if tsd_length > 0 else 'NA'
            })

            # Update exclusion intervals: 20bp on each side of the deletion.
            new_excl_start = max(0, deletion_start - 20)
            new_excl_end = min(chr_length, deletion_start + deletion_length + 20)
            exclusion_intervals.append((new_excl_start, new_excl_end))
            exclusion_intervals = merge_intervals(exclusion_intervals)
            return TE_length
        return 0

    # Loop until the target is met or maximum attempts are reached.
    while attempts < max_global_attempts:
        if target_mode == 'n' and total_insertions_done >= target_value:
            break
        if target_mode == 'p' and total_TE_bp_inserted >= target_value:
            break
        inserted = try_insertion()
        if inserted > 0:
            total_TE_bp_inserted += inserted
            total_insertions_done += 1
            # Uncomment the next line to see per-insertion logs.
#            print(f"[{chrom}] Insertion {total_insertions_done} inserted TE of length {inserted}bp.")
        attempts += 1

    modified_seq = ''.join(seq)
    return (chrom, modified_seq, local_te_bed_entries, total_insertions_done, total_TE_bp_inserted)

def main():
    parser = argparse.ArgumentParser(
        description="Randomly replace backbone sequence with TE sequences (with TSD duplication) while deleting TE_length+TSD_length bp from the genome. The replacement is done in regions at least 20bp away from genes or previously inserted TEs. Mutation percentages for TEs are sampled from an exponential decay distribution."
    )
    parser.add_argument('-genome', required=True, help='Path to the genome FASTA file (genome.fa)')
    parser.add_argument('-TE', required=True, help='Path to the TE FASTA file (TE.fa)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-n', type=int, help='Total number of TE insertions to perform (will be distributed among chromosomes)')
    group.add_argument('-p', type=float, help='Target percent of the genome that should be TE (e.g. -p 20 for 20%% TE; applied per chromosome)')
    parser.add_argument('-bed', required=True, help='Path to the gene BED file')
    parser.add_argument('-output', required=True, help='Prefix for output files (will create prefix.fa and prefix.bed)')
    parser.add_argument('-seed', type=int, default=None, help='Seed for the random number generator (optional)')
    parser.add_argument('-TE_ratio', help='Path to TE ratios file (optional)', default=None)
    parser.add_argument('-stat_out', help='File to output statistics (optional)', default=None)
    # New parameters for exponential decay mutation sampling:
    parser.add_argument('-k', type=float, default=10.0, help='Exponential decay parameter (default 10)')
    parser.add_argument('-Mmax', type=float, default=10.0, help='Maximum mutation percent (default 10)')
    parser.add_argument('-pdf_out', help='Output PDF file for mutation decay function plot (default mutation_decay.pdf)', default="mutation_decay.pdf")
    # New multiprocessing parameter: number of chromosomes to process concurrently.
    parser.add_argument('-m', type=int, default=1, help='Number of chromosomes to process concurrently (multiprocessing; default 1)')
    parser.add_argument('-TsTv', type=float, default=1.0,
                        help='Transition/transversion ratio (default 1.0: Ts and Tv equally likely)')
    args = parser.parse_args()

    # Set global random seed if provided.
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
        print(f"Global random seed set to {args.seed}.")
    else:
        print("No random seed provided. Results will be non-reproducible.")

    # Plot the decay function for mutation percentages.
    # But not if Mmax zero (ie, no mutations added).
#    plot_decay_function(args.k, args.Mmax, args.pdf_out)
    if args.Mmax > 0:
        plot_decay_function(args.k, args.Mmax, args.pdf_out)
    else:
        print("Mmax is 0. Skipping mutation decay plot since mutation is disabled.")

    genome_path = args.genome
    te_path = args.TE
    bed_path = args.bed
    output_prefix = args.output
    k = args.k
    Mmax = args.Mmax

    print("Reading genome FASTA file...")
    genome_dict = parse_fasta(genome_path)
    if not genome_dict:
        print("Error: Genome FASTA file is empty or improperly formatted.")
        sys.exit(1)
    print(f"Loaded {len(genome_dict)} chromosome(s) from genome.")

    print("Reading TE FASTA file...")
    te_dict = parse_fasta(te_path)
    if not te_dict:
        print("Error: TE FASTA file is empty or improperly formatted.")
        sys.exit(1)
    te_headers = list(te_dict.keys())
    print(f"Loaded {len(te_headers)} TE sequence(s).")

    # Map each TE header to its class and superfamily.
    te_info_dict = {}
    for header in te_headers:
        te_class, te_superfamily = extract_te_info(header)
        te_info_dict[header] = {'class': te_class, 'superfamily': te_superfamily}

    # If a TE_ratio file is provided, parse it and build a weighted distribution.
    te_ratio_weights = None
    if args.TE_ratio:
        print("Reading TE ratios file...")
        te_ratio_dict = {}
        try:
            with open(args.TE_ratio, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    te_class, te_superfamily, freq = parts[0], parts[1], float(parts[2])
                    te_ratio_dict[(te_class, te_superfamily)] = freq
        except FileNotFoundError:
            print(f"Error: TE ratios file {args.TE_ratio} not found.")
            sys.exit(1)
        total_ratio = sum(te_ratio_dict.values())
        te_ratio_weights = {}
        for key, freq in te_ratio_dict.items():
            te_ratio_weights[key] = freq / total_ratio
        print("TE ratio weights (normalized):")
        for key, weight in te_ratio_weights.items():
            print(f"  {key[0]} {key[1]}: {weight:.4f}")

    # Build a lookup from (TE_class, TE_superfamily) to list of TE headers.
    te_type_to_headers = {}
    for header, info in te_info_dict.items():
        key_tuple = (info['class'], info['superfamily'])
        if key_tuple not in te_type_to_headers:
            te_type_to_headers[key_tuple] = []
        te_type_to_headers[key_tuple].append(header)

    print("Reading gene BED file...")
    gene_bed_dict = parse_bed(bed_path)
    original_genes = parse_gene_bed(bed_path)
    # Prepare gene BED entries for output.
    gene_bed_entries = []
    for gene in original_genes:
        gene_bed_entries.append({
            'chromosome': gene['chromosome'],
            'start': gene['start'],
            'end': gene['end'],
            'name': gene['name'],
            'strand': '+',
            'tsd': 'NA'
        })

    # Prepare per-chromosome data.
    chromosomes = list(genome_dict.keys())
    chrom_data = {}
    total_genome_bp = 0
    for chrom in chromosomes:
        seq = genome_dict[chrom]
        chr_length = len(seq)
        total_genome_bp += chr_length
        exclusion_intervals = []
        for gene in [g for g in original_genes if g['chromosome'] == chrom]:
            start = max(0, gene['start'] - 20)
            end = min(chr_length, gene['end'] + 20)
            exclusion_intervals.append((start, end))
        exclusion_intervals = merge_intervals(exclusion_intervals)
        chrom_data[chrom] = {'seq': seq, 'length': chr_length, 'exclusion': exclusion_intervals}

    # Determine per-chromosome targets based on -n or -p.
    tasks = []
    if args.n is not None:
        total_insertions = args.n
        # Distribute insertions proportional to chromosome length.
        for idx, chrom in enumerate(chromosomes):
            chr_length = chrom_data[chrom]['length']
            target_n = int(round(total_insertions * (chr_length / total_genome_bp)))
            # Ensure that if rounding gives 0 but chromosome is non-empty, we try at least once.
            if chr_length > 0 and target_n == 0:
                target_n = 1
            tasks.append((chrom, chrom_data[chrom], 'n', target_n,
                          te_headers, te_info_dict, te_dict, te_ratio_weights, te_type_to_headers, k, Mmax,
                          (args.seed + idx) if args.seed is not None else None,
                          args.TsTv))

        print(f"Performing a total of {total_insertions} TE insertions distributed across chromosomes.")
    else:
        target_percent = args.p
        for idx, chrom in enumerate(chromosomes):
            chr_length = chrom_data[chrom]['length']
            target_bp = chr_length * (target_percent / 100)
            tasks.append((chrom, chrom_data[chrom], 'p', target_bp,
                          te_headers, te_info_dict, te_dict, te_ratio_weights, te_type_to_headers, k, Mmax,
                          (args.seed + idx) if args.seed is not None else None,
                          args.TsTv))
        print(f"Performing TE insertions on each chromosome until approximately {target_percent}% TE content is reached.")

    # Process chromosomes concurrently using multiprocessing.
    pool = multiprocessing.Pool(processes=args.m)
    results = pool.map(process_chromosome, tasks)
    pool.close()
    pool.join()

    # Collect modified sequences and TE BED entries.
    modified_genome = {}
    te_bed_entries = []
    total_insertions_done_global = 0
    total_TE_bp_inserted_global = 0
    for res in results:
        chrom, mod_seq, chrom_te_entries, ins_done, te_bp = res
        modified_genome[chrom] = mod_seq
        te_bed_entries.extend(chrom_te_entries)
        total_insertions_done_global += ins_done
        total_TE_bp_inserted_global += te_bp

    # Combine gene BED entries and TE BED entries.
    output_bed_entries = gene_bed_entries + te_bed_entries
    output_bed_entries_sorted = sorted(output_bed_entries, key=lambda x: (x['chromosome'], x['start']))

    output_genome_path = f"{output_prefix}.fa"
    output_bed_path = f"{output_prefix}.bed"

    print(f"Writing modified genome to {output_genome_path}...")
    write_fasta(output_genome_path, modified_genome)

    print(f"Writing BED file to {output_bed_path}...")
    write_bed(output_bed_path, output_bed_entries_sorted)

    genome_length = sum(len(seq) for seq in modified_genome.values())
    gene_total_bp = sum(gene['end'] - gene['start'] for gene in original_genes)
    TE_total_bp = sum(entry['end'] - entry['start'] for entry in te_bed_entries)
    gene_percent = (gene_total_bp / genome_length) * 100 if genome_length > 0 else 0
    TE_percent = (TE_total_bp / genome_length) * 100 if genome_length > 0 else 0
    gene_count = len(original_genes)
    TE_count = len(te_bed_entries)

    stat_message = (f"The burn-in genome is {genome_length}bp in length with {gene_count} genes "
                    f"({gene_percent:.2f}%) and {TE_count} TEs ({TE_percent:.2f}%).\n"
                    f"Total insertions done: {total_insertions_done_global} (TE bp inserted: {total_TE_bp_inserted_global}).")
    print("Replacement process completed successfully.")
    print(stat_message)

    if args.stat_out:
        try:
            with open(args.stat_out, 'w') as stat_file:
                stat_file.write(stat_message + "\n")
            print(f"Statistics also written to {args.stat_out}.")
        except Exception as e:
            print(f"Error writing statistics to {args.stat_out}: {e}")

if __name__ == "__main__":
    main()
