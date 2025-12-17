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
    match = re.match(r"[^#]+#([^/]+)/([^~]+)", header)
    if match:
        te_class = match.group(1).strip()
        te_superfamily = match.group(2).strip()
        return te_class, te_superfamily
    else:
        return "unknown", "unknown"

def get_tsd_length(te_class, te_superfamily):
    """
    Determine the TSD length based on TE class and superfamily.
    Case-insensitive. Returns an integer. Variable TSDs are sampled within ranges.
    """
    c = (te_class or "").strip().lower()
    s = (te_superfamily or "").strip().lower()

    _const5 = lambda: 5
    _hat_rng = lambda: random.randint(5, 8)
    _cacta_rng = lambda: random.randint(2, 4)
    _tourist_rng = lambda: random.randint(2, 10)
    _mule_rng = lambda: random.randint(8, 9)
    _l1_rng = lambda: random.randint(5, 20)
    _unknown_tir_rng = lambda: random.choice([2, 10])

    tsd_mapping = {
        ('ltr', 'copia'): _const5,
        ('ltr', 'gypsy'): _const5,
        ('ltr', 'solo'):  _const5,
        ('ltr', 'ty3'):   _const5,
        ('ltr', 'unknown'): _const5,
        ('ltr', 'crm'): _const5,
        ('ltr', 'trim'): _const5,

        ('dna', 'harbinger'): lambda: 3,
        ('dna', 'mariner'):   lambda: 2,

        ('dnaauto', 'helitron'): lambda: 0,
        ('dna',     'helitron'): lambda: 0,
        ('dnanona', 'helitron'): lambda: 0,

        ('mite', 'stow'):     _tourist_rng,
        ('mite', 'tourist'):  _tourist_rng,
        ('sine', 'trna'):     _l1_rng,
        ('sine', 'unknown'):  _l1_rng,
        ('line', 'l1'):       _l1_rng,
        ('line', 'unknown'):  _l1_rng,
        ('line', 'rte'):      _l1_rng,

        ('dnaauto', 'cacta'): _cacta_rng,
        ('dna',     'cacta'): _cacta_rng,
        ('dnanona', 'cacta'): _cacta_rng,

        ('dna',     'hat'):   _hat_rng,
        ('dnanona', 'hat'):   _hat_rng,
        ('dnaauto', 'hat'):   _hat_rng,

        ('dnaauto', 'mudr'):  _mule_rng,
        ('dnanona', 'mule'):  _mule_rng,
        ('dna',     'mudr'):  _mule_rng,
        ('dnaauto', 'mule'):  _mule_rng,
        ('dna', 'mutator'): _mule_rng,

        ('dna',  'dta'): _hat_rng,
        ('mite', 'dta'): _hat_rng,

        ('dna',  'dtc'): _cacta_rng,
        ('mite', 'dtc'): _cacta_rng,

        ('dna',  'dth'): _tourist_rng,
        ('mite', 'dth'): _tourist_rng,

        ('dna',  'dtm'): _mule_rng,
        ('mite', 'dtm'): _mule_rng,

        ('dna',  'dtt'): _tourist_rng,
        ('mite', 'dtt'): _tourist_rng,

        ('dnaauto', 'cactg'): _cacta_rng,
        ('dnanona', 'cactg'): _cacta_rng,

        ('dnaauto', 'mle'): _tourist_rng,
        ('dnanona', 'mle'): _tourist_rng, 

        ('dnanona', 'muletir'): _mule_rng,

        ('dnanona', 'tourist'): _tourist_rng,
        ('dna', 'tc1_mariner'): _tourist_rng,

        ('dnaauto', 'pile'): _unknown_tir_rng,
        ('dnaauto', 'pole'): _unknown_tir_rng,
        ('dnanona', 'pile'): _unknown_tir_rng,
        ('dnanona', 'pole'): _unknown_tir_rng,
        ('dnanona', 'unknown'): _unknown_tir_rng,
        ('dna', 'unknown'): _unknown_tir_rng,
    }

    func = tsd_mapping.get((c, s))
    if func:
        return func()
    else:
        print(f"Warning: TSD length not defined for TE class '{te_class}' and superfamily '{te_superfamily}'. Using default TSD length of 5.")
        return 5

def reverse_complement(seq):
    complement = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(complement)[::-1]

def compute_allowed_intervals(seq_length, exclusion_intervals):
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

def mutate_sequence(seq, mutation_rate, ts_tv_ratio):
    ts_map = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
    tv_map = {
        'A': ['C', 'T'], 'G': ['C', 'T'],
        'C': ['A', 'G'], 'T': ['A', 'G']
    }

    transitioned = 0
    transverted = 0
    mutated_bases = 0
    mutated = []

    for base in seq:
        if random.random() < mutation_rate:
            b = base.upper()
            if b in ts_map:
                p_ts = ts_tv_ratio / (ts_tv_ratio + 1.0)
                if random.random() < p_ts:
                    new = ts_map[b]
                    transitioned += 1
                else:
                    new = random.choice(tv_map[b])
                    transverted += 1
                mutated_bases += 1
                mutated.append(new if base.isupper() else new.lower())
            else:
                mutated.append(base)
        else:
            mutated.append(base)

    return ''.join(mutated), mutated_bases, transitioned, transverted

def plot_decay_function(k, Mmax, pdf_out):
    M_values = np.linspace(0, Mmax, 500)
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

def plot_decay_with_empirical(k, Mmax, mutation_bins, intact_percents, frag_percents, pdf_out):
    # Decide x-range (in percent units)
    if mutation_bins is not None:
        xmax = max(1.0, max(b['end'] for b in mutation_bins) * 100.0)  # rates -> %
    else:
        xmax = max(1.0, float(Mmax))

    plt.figure(figsize=(6,4))

    # Use a common plotting grid so "target" and "empirical" are comparable
    plot_edges = np.linspace(0.0, xmax, 60)
    plot_widths = np.diff(plot_edges)
    plot_mids = (plot_edges[:-1] + plot_edges[1:]) / 2.0

    if mutation_bins is not None:
        # Convert original bins (in 0–1 rate space) into a smooth density on plot_edges
        # by distributing each bin's probability mass into overlapping plot bins.
        target_mass = np.zeros(len(plot_widths), dtype=float)

        prev = 0.0
        for b in mutation_bins:
            w = b['cum'] - prev   # probability mass in this original bin
            prev = b['cum']
            if w <= 0:
                continue

            s = b['start'] * 100.0  # to %
            e = b['end']   * 100.0
            if e <= s:
                continue

            # distribute mass proportional to overlap length
            # find overlap with each plot bin
            for i in range(len(plot_widths)):
                a = plot_edges[i]
                c = plot_edges[i+1]
                overlap = max(0.0, min(e, c) - max(s, a))
                if overlap > 0:
                    target_mass[i] += w * (overlap / (e - s))

        # Convert mass per plot bin to density per percent
        target_density = np.divide(
            target_mass, plot_widths,
            out=np.zeros_like(target_mass),
            where=plot_widths > 0
        )

        plt.plot(plot_mids, target_density, lw=2, label="Target (bins)")

    else:
        # Original exponential mode (as before)
        M_values = np.linspace(0, xmax, 500)
        norm_const = 1 - math.exp(-k)
        f_values = (k / (xmax * norm_const)) * np.exp(-k * (M_values / xmax))
        plt.plot(M_values, f_values, lw=2, label="Target (exponential)")

    # Helper for empirical curves (histogram density -> line) using the SAME plot_edges
    def plot_empirical(data, label):
        if not data:
            return
        data = np.array(data, dtype=float)
        data = data[(data >= 0) & (data <= xmax)]
        if data.size == 0:
            return
        hist, _ = np.histogram(data, bins=plot_edges, density=True)
        plt.plot(plot_mids, hist, lw=2, label=label)

    plot_empirical(intact_percents, "Empirical intact")
    plot_empirical(frag_percents, "Empirical fragmented")

    plt.xlabel("Mutation Percent (%)", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.title("Target vs Empirical Mutation Distributions", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(pdf_out, format='pdf')
    plt.close()
    print(f"Target + empirical mutation plot saved to {pdf_out}.")


def parse_mutation_bins(file_path):
    """
    Parse a mutation bin file:
        <bin_start> <bin_end> <frequency>

    bin_start / bin_end are mutation *rates* (0–1, e.g. 0.006 = 0.6% of bases mutated).
    frequency is a weight; the script normalizes them so they sum to 1.
    """
    bins = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                start = float(parts[0])
                end = float(parts[1])
                freq = float(parts[2])
                bins.append({'start': start, 'end': end, 'freq': freq})
    except FileNotFoundError:
        print(f"Error: Mutation bin file {file_path} not found.")
        sys.exit(1)

    if not bins:
        print(f"Error: Mutation bin file {file_path} contained no valid rows.")
        sys.exit(1)

    total = sum(b['freq'] for b in bins)
    if total <= 0:
        print(f"Error: Mutation bin frequencies in {file_path} sum to zero or negative.")
        sys.exit(1)

    # Normalize and build cumulative distribution
    cum = 0.0
    for b in bins:
        cum += b['freq'] / total
        b['cum'] = cum
    bins[-1]['cum'] = 1.0  # ensure final bin closes at 1.0

    return bins


def sample_mutation_rate(mutation_bins, k, Mmax):
    """
    Return a mutation RATE (0–1).

    If mutation_bins is not None:
        1) choose a bin according to its frequency
        2) draw uniformly between [bin_start, bin_end] for that bin
    Else:
        use the existing exponential approach based on k and Mmax.
    """
    # Bin-based mode
    if mutation_bins is not None:
        r = random.random()
        for b in mutation_bins:
            if r <= b['cum']:
                if b['start'] == b['end']:
                    return b['start']
                return random.uniform(b['start'], b['end'])
        # Fallback, should never really happen
        last = mutation_bins[-1]
        return (last['start'] + last['end']) / 2.0

    # Original exponential mode (unchanged logic)
    if Mmax <= 0:
        return 0.0
    u = random.random()
    mutation_percent = - (Mmax / k) * math.log(1 - u * (1 - math.exp(-k)))
    return mutation_percent / 100.0

def extract_ltr_length(header):
    """
    Extract LTR length from TE header metadata like:
      ...#LTR/Copia~LTRlen:1000~LTRlen:1000
    Returns an int or None if not found / invalid.
    If multiple LTRlen values exist, uses the first but warns if they disagree.
    """
    matches = re.findall(r"LTRlen:(\d+)", header)
    if not matches:
        return None
    vals = [int(x) for x in matches if x.isdigit()]
    if not vals:
        return None
    if len(set(vals)) > 1:
        print(f"Warning: Multiple different LTRlen values in header '{header}'. Using first: {vals[0]}")
    return vals[0]


def mutate_ltr3_from_ltr5(ltr5, target_divergence_percent, ts_tv_ratio):
    """
    Create a mutated 3' LTR by copying the 5' LTR and mutating EXACTLY
    round(LTRlen * divergence%) positions.

    Returns: (ltr3_mut, n_mut, n_ts, n_tv)
    """
    ltr5 = ltr5.upper()
    L = len(ltr5)
    if L == 0:
        return ltr5, 0, 0, 0

    # exact number of differences
    n_mut = int(round(L * (target_divergence_percent / 100.0)))
    n_mut = max(0, min(L, n_mut))
    if n_mut == 0:
        return ltr5, 0, 0, 0

    ts_map = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
    tv_map = {
        'A': ['C', 'T'], 'G': ['C', 'T'],
        'C': ['A', 'G'], 'T': ['A', 'G']
    }

    transitioned = 0
    transverted = 0

    ltr3 = list(ltr5)
    positions = random.sample(range(L), k=n_mut)

    p_ts = ts_tv_ratio / (ts_tv_ratio + 1.0) if (ts_tv_ratio is not None) else 0.5

    for i in positions:
        b = ltr5[i]
        if b not in ts_map:
            # For non-ACGT, leave it unchanged (and reduce actual divergence by 1)
            # Alternatively, you could force-change N to A/C/G/T, but this keeps it conservative.
            continue

        if random.random() < p_ts:
            new = ts_map[b]
            transitioned += 1
        else:
            new = random.choice(tv_map[b])
            transverted += 1

        # Ensure we actually change the base (should always be true here)
        if new == b:
            # ultra-defensive fallback
            choices = [x for x in "ACGT" if x != b]
            new = random.choice(choices)

        ltr3[i] = new

    # NOTE: If some chosen positions were non-ACGT, true divergence may be slightly < requested.
    # If you want to guarantee EXACT divergence even with Ns, we can re-sample positions until n_mut changes occur.
    return ''.join(ltr3), n_mut, transitioned, transverted

def add_frag_tag(header):
    """
    Insert _FRAG before the first '#' in header; if no '#', append.
    Example:
      TE_00004120#DNA/Helitron -> TE_00004120_FRAG#DNA/Helitron
    """
    if '#' in header:
        i = header.index('#')
        return header[:i] + "_FRAG" + header[i:]
    else:
        return header + "_FRAG"

def parse_te_ratio_file(path, need_fragmented=False):
    """
    Returns two dicts:
      - intact_weights: {(class, superfamily): weight}
      - frag_weights: {(class, superfamily): weight} or None if not needed
    If need_fragmented=True, require a 4-column file with fragmented_frequency.
    If need_fragmented=False, accept 3-column file (intact only).
    """
    try:
        with open(path, 'r') as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    except FileNotFoundError:
        print(f"Error: TE ratios file {path} not found.")
        sys.exit(1)

    intact = {}
    frag = {}

    for line in lines:
        parts = line.split()
        if need_fragmented:
            if len(parts) < 4:
                print("Error: Fragmented insertions requested but TE_ratio file lacks 'fragmented_frequency' column.")
                sys.exit(1)
            te_class, te_superfamily, intact_freq, frag_freq = parts[0], parts[1], float(parts[2]), float(parts[3])
            intact[(te_class, te_superfamily)] = intact_freq
            frag[(te_class, te_superfamily)] = frag_freq
        else:
            if len(parts) < 3:
                # ignore malformed
                continue
            te_class, te_superfamily, freq = parts[0], parts[1], float(parts[2])
            intact[(te_class, te_superfamily)] = freq

    # Normalize
    def normalize(d):
        total = sum(d.values())
        if total <= 0:
            return {k: 1.0/len(d) for k in d} if d else {}
        return {k: v/total for k, v in d.items()}

    intact_weights = normalize(intact)
    frag_weights = normalize(frag) if frag else ({} if need_fragmented else None)

    return intact_weights, frag_weights

def process_chromosome(task):
    """
    Process TE insertions on a single chromosome.

    task is a tuple containing:
      (chrom, data,
       mode_intact, target_intact, mode_frag, target_frag,
       te_headers, te_info_dict, te_dict,
       te_type_to_headers, intact_weights, frag_weights,
       k, Mmax, seed_offset, ts_tv_ratio)
    mode_*: 'n' (count) or 'p' (bp) or None
    target_*: numeric target (int for n, float/bp for p) or 0
    """
    (chrom, data,
     mode_intact, target_intact, mode_frag, target_frag,
     te_headers, te_info_dict, te_dict,
     te_type_to_headers, intact_weights, frag_weights,
     mutation_bins,
     k, Mmax, seed_offset, ts_tv_ratio) = task

    if seed_offset is not None:
        random.seed(seed_offset)
        np.random.seed(seed_offset)

    seq = list(data['seq'])
    chr_length = data['length']
    exclusion_intervals = list(data['exclusion'])
    local_te_bed_entries = []
    empirical_intact = []
    empirical_frag = []

    total_TE_bp_inserted_intact = 0
    total_insertions_done_intact = 0

    total_TE_bp_inserted_frag = 0
    total_insertions_done_frag = 0

    attempts = 0
    max_global_attempts = 100000

    # helpers
    def pick_header(weights):
        if weights:
            keys = list(weights.keys())
            wts = [weights[k] for k in keys]
            chosen_key = random.choices(keys, weights=wts, k=1)[0]
            # chosen_key from TE_ratio is (class, superfamily).
            if chosen_key in te_type_to_headers and te_type_to_headers[chosen_key]:
                return random.choice(te_type_to_headers[chosen_key])
        # fallback uniform
        return random.choice(te_headers)

    def remaining_for(mode, target, done_count, done_bp):
        if mode == 'n':
            return max(0, int(target) - done_count)
        elif mode == 'p':
            return max(0.0, float(target) - done_bp)
        else:
            return 0

    def choose_type():
        """Interleave by remaining requirement."""
        rem_i = remaining_for(mode_intact, target_intact, total_insertions_done_intact, total_TE_bp_inserted_intact)
        rem_f = remaining_for(mode_frag, target_frag, total_insertions_done_frag, total_TE_bp_inserted_frag)
        both_zero = (rem_i == 0 and rem_f == 0)
        if both_zero:
            return None
        if rem_i == 0:
            return 'frag'
        if rem_f == 0:
            return 'intact'
        # Weight by remaining requirement
        wi = rem_i
        wf = rem_f
        # If percent mode, they are bp (floats); if count, ints. Both ok for weights.
        return random.choices(['intact', 'frag'], weights=[wi, wf], k=1)[0]

    def compute_allowed(deletion_length, tsd_length):
        allowed = compute_allowed_intervals(chr_length, exclusion_intervals)
        candidate_intervals = []
        total_candidates = 0
        for (a, b) in allowed:
            lower_bound = a + tsd_length
            upper_bound = b - deletion_length
            if lower_bound <= upper_bound:
                candidate_intervals.append((lower_bound, upper_bound))
                total_candidates += (upper_bound - lower_bound + 1)
        return candidate_intervals, total_candidates

    def try_insertion(which):
        nonlocal exclusion_intervals, seq
        max_attempts_local = 10
        for _ in range(max_attempts_local):
            if which == 'intact':
                selected_te = pick_header(intact_weights)
                te_class = te_info_dict[selected_te]['class']
                te_superfamily = te_info_dict[selected_te]['superfamily']
                tsd_length = get_tsd_length(te_class, te_superfamily)
                te_sequence = te_dict[selected_te]
                # mutation
                u = random.random()
                # mutation_percent = - (Mmax / k) * math.log(1 - u * (1 - math.exp(-k))) if Mmax > 0 else 0.0
                # mutation (rate is 0–1)
                mutation_rate = sample_mutation_rate(mutation_bins, k, Mmax)

                te_sequence_before = te_sequence

                # NEW: LTRlen-aware mutation for intact LTRs
                ltr_len = None
                if (te_class or "").strip().lower() == "ltr":
                    ltr_len = extract_ltr_length(selected_te)

                if ltr_len is not None and ltr_len > 0 and len(te_sequence_before) >= 2 * ltr_len:
                    five_ltr = te_sequence_before[:ltr_len]
                    internal = te_sequence_before[ltr_len:len(te_sequence_before) - ltr_len]

                    # target divergence percent comes from the SAME distribution you already use
                    target_div_pct = mutation_rate * 100.0

                    three_ltr_mut, n_mut, n_ts, n_tv = mutate_ltr3_from_ltr5(
                        five_ltr, target_div_pct, ts_tv_ratio
                    )

                    te_sequence = five_ltr + internal + three_ltr_mut

                    # empirical percent: divergence within the 3' LTR relative to 5' LTR
                    emp_pct = (n_mut / float(ltr_len)) * 100.0 if ltr_len > 0 else 0.0
                    empirical_intact.append(emp_pct)

                else:
                    # fallback: original behavior (mutate entire TE)
                    te_sequence, n_mut, n_ts, n_tv = mutate_sequence(te_sequence_before, mutation_rate, ts_tv_ratio)

                    emp_pct = (n_mut / len(te_sequence_before) * 100.0) if len(te_sequence_before) > 0 else 0.0
                    empirical_intact.append(emp_pct)

                TE_length = len(te_sequence)
                deletion_length = TE_length + tsd_length

                candidate_intervals, total_candidates = compute_allowed(deletion_length, tsd_length)
                if total_candidates <= 0:
                    continue

                r = random.randint(0, total_candidates - 1)
                chosen_deletion_start = None
                for (low, high) in candidate_intervals:
                    interval_length = high - low + 1
                    if r < interval_length:
                        chosen_deletion_start = low + r
                        break
                    r -= interval_length
                if chosen_deletion_start is None:
                    continue
                deletion_start = chosen_deletion_start

                tsd_seq = ''.join(seq[deletion_start - tsd_length: deletion_start]) if tsd_length > 0 else ''

                strand = random.choice(['+', '-'])
                te_seq_final = reverse_complement(te_sequence) if strand == '-' else te_sequence

                replacement = te_seq_final + tsd_seq
                seq[deletion_start : deletion_start + deletion_length] = list(replacement)

                name_for_bed = selected_te
                local_te_bed_entries.append({
                    'chromosome': chrom,
                    'start': deletion_start,
                    'end': deletion_start + len(te_seq_final),
                    'name': name_for_bed,
                    'strand': strand,
                    'tsd': tsd_seq if tsd_length > 0 else 'NA'
                })

                # expand exclusion ±20 bp around the full deletion footprint
                new_excl_start = max(0, deletion_start - 20)
                new_excl_end = min(chr_length, deletion_start + deletion_length + 20)
                exclusion_intervals.append((new_excl_start, new_excl_end))
                exclusion_intervals = merge_intervals(exclusion_intervals)

                return TE_length, len(te_seq_final), 'intact'  # (deleted TE len, inserted len, type)

            else:  # 'frag'
                selected_te = pick_header(frag_weights if frag_weights is not None else None)
                te_sequence_full = te_dict[selected_te]
                # fragmenting: chop between 20% and 80% from one end
                pct_chop = random.randint(20, 80) / 100.0
                chop_left = random.choice([True, False])
                cut = int(round(len(te_sequence_full) * pct_chop))
                if chop_left:
                    te_sequence = te_sequence_full[cut:]
                else:
                    te_sequence = te_sequence_full[:len(te_sequence_full) - cut]
                if len(te_sequence) <= 0:
                    continue  # extremely unlikely with 20–80 range, but guard anyway

                # mutation
                u = random.random()
                # mutation_percent = - (Mmax / k) * math.log(1 - u * (1 - math.exp(-k))) if Mmax > 0 else 0.0
                # mutation (rate is 0–1)
                mutation_rate = sample_mutation_rate(mutation_bins, k, Mmax)
                te_sequence_before = te_sequence
                te_sequence, n_mut, n_ts, n_tv = mutate_sequence(te_sequence_before, mutation_rate, ts_tv_ratio)

                emp_pct = (n_mut / len(te_sequence_before) * 100.0) if len(te_sequence_before) > 0 else 0.0
                empirical_frag.append(emp_pct)

                tsd_length = 0
                TE_length = len(te_sequence)  # no TSD
                deletion_length = TE_length

                candidate_intervals, total_candidates = compute_allowed(deletion_length, tsd_length)
                if total_candidates <= 0:
                    continue

                r = random.randint(0, total_candidates - 1)
                chosen_deletion_start = None
                for (low, high) in candidate_intervals:
                    interval_length = high - low + 1
                    if r < interval_length:
                        chosen_deletion_start = low + r
                        break
                    r -= interval_length
                if chosen_deletion_start is None:
                    continue
                deletion_start = chosen_deletion_start

                strand = random.choice(['+', '-'])
                te_seq_final = reverse_complement(te_sequence) if strand == '-' else te_sequence

                replacement = te_seq_final  # no TSD
                seq[deletion_start : deletion_start + deletion_length] = list(replacement)

                name_for_bed = add_frag_tag(selected_te)
                local_te_bed_entries.append({
                    'chromosome': chrom,
                    'start': deletion_start,
                    'end': deletion_start + len(te_seq_final),
                    'name': name_for_bed,
                    'strand': strand,
                    'tsd': 'NA'
                })

                new_excl_start = max(0, deletion_start - 20)
                new_excl_end = min(chr_length, deletion_start + deletion_length + 20)
                exclusion_intervals.append((new_excl_start, new_excl_end))
                exclusion_intervals = merge_intervals(exclusion_intervals)

                return TE_length, len(te_seq_final), 'frag'

        return 0, 0, None

    def target_met():
        ri = remaining_for(mode_intact, target_intact, total_insertions_done_intact, total_TE_bp_inserted_intact)
        rf = remaining_for(mode_frag, target_frag, total_insertions_done_frag, total_TE_bp_inserted_frag)
        return (ri == 0 and rf == 0)

    # main loop
    while attempts < max_global_attempts:
        if target_met():
            break
        which = choose_type()
        if which is None:
            break
        deleted_len, ins_len, typ = try_insertion(which)
        if ins_len > 0 and typ == 'intact':
            total_TE_bp_inserted_intact += ins_len
            total_insertions_done_intact += 1
        elif ins_len > 0 and typ == 'frag':
            total_TE_bp_inserted_frag += ins_len
            total_insertions_done_frag += 1
        attempts += 1

    modified_seq = ''.join(seq)
    return (chrom, modified_seq, local_te_bed_entries,
            total_insertions_done_intact, total_TE_bp_inserted_intact,
            total_insertions_done_frag, total_TE_bp_inserted_frag,
            empirical_intact, empirical_frag)

def main():
    parser = argparse.ArgumentParser(
        description=("Randomly replace backbone sequence with TE sequences (with or without TSDs). "
                     "Supports intact insertions (with TSD) and fragmented insertions (no TSD). "
                     "Mutation percents drawn from an exponential decay distribution.")
    )
    parser.add_argument('-genome', required=True, help='Path to the genome FASTA file (genome.fa)')
    parser.add_argument('-TE', required=True, help='Path to the TE FASTA file (TE.fa)')

    # Intact target: count OR percent
    group_intact = parser.add_mutually_exclusive_group()
    group_intact.add_argument('-n_intact', type=int, help='Total number of INTACT TE insertions to perform (distributed across chromosomes)')
    group_intact.add_argument('-p_intact', type=float, help='Target percent of genome to be INTACT TE bp per chromosome (e.g. 20 for 20%%)')

    # Fragmented target: count OR percent
    group_frag = parser.add_mutually_exclusive_group()
    group_frag.add_argument('-n_frag', type=int, help='Total number of FRAGMENTED TE insertions to perform (distributed across chromosomes)')
    group_frag.add_argument('-p_frag', type=float, help='Target percent of genome to be FRAGMENTED TE bp per chromosome (e.g. 5 for 5%%)')

    parser.add_argument('-bed', required=True, help='Path to the gene BED file')
    parser.add_argument('-output', required=True, help='Prefix for output files (will create prefix.fa and prefix.bed)')
    parser.add_argument('-seed', type=int, default=None, help='Seed for the random number generator (optional)')
    parser.add_argument('-TE_ratio', help='Path to TE ratios file (see docs below)', default=None)
    parser.add_argument('-stat_out', help='File to output statistics (optional)', default=None)

    # Mutation parameters
    parser.add_argument('-k', type=float, default=10.0, help='Exponential decay parameter (default 10)')
    parser.add_argument('-Mmax', type=float, default=10.0, help='Maximum mutation percent (default 10)')
    parser.add_argument('-pdf_out', help='Output PDF file for mutation decay function plot (default mutation_decay.pdf)', default="mutation_decay.pdf")

    parser.add_argument(
        '-mutation_bins',
        help='Optional file of mutation bins: <bin_start> <bin_end> <frequency> per line. '
             'If provided, overrides -k and -Mmax and disables the decay plot.',
        default=None
    )

    # Multiprocessing
    parser.add_argument('-m', type=int, default=1, help='Number of chromosomes to process concurrently (default 1)')

    parser.add_argument('-TsTv', type=float, default=1.0,
                        help='Transition/transversion ratio (default 1.0: Ts and Tv equally likely)')

    args = parser.parse_args()

    # Sanity: require at least one target among the four to be provided and >0
    provided = [
        ('n_intact', args.n_intact),
        ('p_intact', args.p_intact),
        ('n_frag', args.n_frag),
        ('p_frag', args.p_frag),
    ]
    if not any(v for _, v in provided):
        print("Error: You must provide at least one of -n_intact, -p_intact, -n_frag, or -p_frag.")
        sys.exit(1)

    # Set seed
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
        print(f"Global random seed set to {args.seed}.")
    else:
        print("No random seed provided. Results will be non-reproducible.")

    # Mutation bins: override k/Mmax if provided
    mutation_bins = None
    if args.mutation_bins:
        print(f"Reading mutation bin file: {args.mutation_bins}")
        mutation_bins = parse_mutation_bins(args.mutation_bins)

    # Plot decay function only if we're using the exponential mode
    # Do NOT plot yet; empirical lists are only available after insertions.

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

    # Build a lookup from (TE_class, TE_superfamily) to list of TE headers.
    te_type_to_headers = {}
    for header, info in te_info_dict.items():
        key_tuple = (info['class'], info['superfamily'])
        te_type_to_headers.setdefault(key_tuple, []).append(header)

    # TE ratios (weights)
    need_fragmented = (args.n_frag and args.n_frag > 0) or (args.p_frag and args.p_frag > 0)
    intact_weights = None
    frag_weights = None
    if args.TE_ratio:
        intact_weights, frag_weights = parse_te_ratio_file(args.TE_ratio, need_fragmented=bool(need_fragmented))
        print("TE ratio weights (normalized):")
        for key, weight in (intact_weights or {}).items():
            print(f"  intact {key[0]} {key[1]}: {weight:.4f}")
        if need_fragmented:
            for key, weight in (frag_weights or {}).items():
                print(f"  frag   {key[0]} {key[1]}: {weight:.4f}")
    else:
        # If no TE_ratio provided, sampling is uniform.
        intact_weights = None
        frag_weights = None
        if need_fragmented:
            print("Warning: Fragmented insertions requested but no TE_ratio provided. Sampling uniformly for both intact and fragmented.")

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

    # Determine per-chromosome targets
    tasks = []

    # helpers for distributing targets proportional to chromosome length
    def dist_counts(total_count):
        # proportional, ensure sum equals requested; guarantee at least 1 for non-empty chromosomes
        alloc = []
        running = 0
        for chrom in chromosomes:
            chr_len = chrom_data[chrom]['length']
            n = int(round(total_count * (chr_len / total_genome_bp))) if total_genome_bp > 0 else 0
            if chr_len > 0 and total_count > 0 and n == 0:
                n = 1
            alloc.append(n)
            running += n
        # Adjust if rounding drifted
        delta = running - (total_count if total_count is not None else 0)
        i = 0
        while delta > 0 and any(a > 0 for a in alloc):
            if alloc[i % len(alloc)] > 0:
                alloc[i % len(alloc)] -= 1
                delta -= 1
            i += 1
        while delta < 0:
            alloc[i % len(alloc)] += 1
            delta += 1
            i += 1
        return {chromosomes[i]: alloc[i] for i in range(len(chromosomes))}

    def dist_bp(percent):
        # percent is e.g. 20 meaning "20% of chr bp"
        return {chrom: chrom_data[chrom]['length'] * (percent / 100.0) for chrom in chromosomes}

    # Build per-chromosome targets for intact
    if args.n_intact is not None and args.n_intact > 0:
        perchr_n_intact = dist_counts(args.n_intact)
        mode_intact = 'n'
        perchr_t_intact = perchr_n_intact
        print(f"Performing a total of {args.n_intact} INTACT insertions distributed across chromosomes.")
    elif args.p_intact is not None and args.p_intact > 0:
        perchr_t_intact = dist_bp(args.p_intact)
        mode_intact = 'p'
        print(f"Performing INTACT insertions on each chromosome until ~{args.p_intact}% TE content is reached.")
    else:
        perchr_t_intact = {chrom: 0 for chrom in chromosomes}
        mode_intact = None

    # Build per-chromosome targets for fragmented
    if args.n_frag is not None and args.n_frag > 0:
        perchr_n_frag = dist_counts(args.n_frag)
        mode_frag = 'n'
        perchr_t_frag = perchr_n_frag
        print(f"Performing a total of {args.n_frag} FRAGMENTED insertions distributed across chromosomes.")
    elif args.p_frag is not None and args.p_frag > 0:
        perchr_t_frag = dist_bp(args.p_frag)
        mode_frag = 'p'
        print(f"Performing FRAGMENTED insertions on each chromosome until ~{args.p_frag}% TE content is reached.")
    else:
        perchr_t_frag = {chrom: 0 for chrom in chromosomes}
        mode_frag = None

    # Compose tasks
    for idx, chrom in enumerate(chromosomes):
        t_intact = perchr_t_intact[chrom]
        t_frag = perchr_t_frag[chrom]
        tasks.append(
            (chrom, chrom_data[chrom],
             mode_intact, t_intact,
             mode_frag, t_frag,
             te_headers, te_info_dict, te_dict,
             te_type_to_headers, intact_weights, frag_weights,
             mutation_bins,
             k, Mmax,
             (args.seed + idx) if args.seed is not None else None,
             args.TsTv)
        )

    # Process chromosomes concurrently
    pool = multiprocessing.Pool(processes=args.m)
    results = pool.map(process_chromosome, tasks)
    pool.close()
    pool.join()

    # Collect modified sequences and TE BED entries.
    modified_genome = {}
    te_bed_entries = []
    total_insertions_done_intact_global = 0
    total_TE_bp_inserted_intact_global = 0
    total_insertions_done_frag_global = 0
    total_TE_bp_inserted_frag_global = 0
    
    empirical_intact_all = []
    empirical_frag_all = []

    for res in results:
        (chrom, mod_seq, chrom_te_entries,
         ins_done_intact, te_bp_intact,
         ins_done_frag, te_bp_frag,
         emp_i, emp_f) = res

        modified_genome[chrom] = mod_seq
        te_bed_entries.extend(chrom_te_entries)
        total_insertions_done_intact_global += ins_done_intact
        total_TE_bp_inserted_intact_global += te_bp_intact
        total_insertions_done_frag_global += ins_done_frag
        total_TE_bp_inserted_frag_global += te_bp_frag

        empirical_intact_all.extend(emp_i)
        empirical_frag_all.extend(emp_f)

    print(f"Empirical counts: intact={len(empirical_intact_all)} frag={len(empirical_frag_all)}")

    if args.Mmax > 0 or mutation_bins is not None:
        plot_decay_with_empirical(args.k, args.Mmax, mutation_bins,
                                  empirical_intact_all, empirical_frag_all,
                                  args.pdf_out)
    else:
        print("Mmax is 0 and no mutation bins provided. Skipping mutation plot.")

    # Combine gene BED entries and TE BED entries.
    output_bed_entries = gene_bed_entries + te_bed_entries
    output_bed_entries_sorted = sorted(output_bed_entries, key=lambda x: (x['chromosome'], x['start']))

    output_genome_path = f"{output_prefix}.fasta"
    output_bed_path = f"{output_prefix}.bed"

    print(f"Writing modified genome to {output_genome_path}...")
    write_fasta(output_genome_path, modified_genome)

    print(f"Writing BED file to {output_bed_path}...")
    write_bed(output_bed_path, output_bed_entries_sorted)

    genome_length = sum(len(seq) for seq in modified_genome.values())
    gene_total_bp = sum(gene['end'] - gene['start'] for gene in original_genes)

    # Compute TE-only stats from entries (intact+frag already mixed in entries)
    TE_total_bp = sum(entry['end'] - entry['start'] for entry in te_bed_entries)

    # Percentages
    gene_percent = (gene_total_bp / genome_length) * 100 if genome_length > 0 else 0
    TE_percent_total = (TE_total_bp / genome_length) * 100 if genome_length > 0 else 0

    # Split intact vs frag percentages using the tracked bp totals
    TE_percent_intact = (total_TE_bp_inserted_intact_global / genome_length) * 100 if genome_length > 0 else 0
    TE_percent_frag = (total_TE_bp_inserted_frag_global / genome_length) * 100 if genome_length > 0 else 0

    gene_count = len(original_genes)
    TE_count_intact = total_insertions_done_intact_global
    TE_count_frag = total_insertions_done_frag_global
    TE_count_total = TE_count_intact + TE_count_frag

    stat_message = (
        f"The burn-in genome is {genome_length}bp in length with {gene_count} genes "
        f"({gene_percent:.2f}%) and {TE_count_total} TEs ({TE_percent_total:.2f}%).\n"
        f"  INTACT:  insertions={TE_count_intact}, TE bp inserted={total_TE_bp_inserted_intact_global}, "
        f"({TE_percent_intact:.2f}% of genome)\n"
        f"  FRAG   : insertions={TE_count_frag}, TE bp inserted={total_TE_bp_inserted_frag_global}, "
        f"({TE_percent_frag:.2f}% of genome)"
    )

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
