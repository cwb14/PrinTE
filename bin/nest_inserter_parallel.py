#!/usr/bin/env python3
"""
Insert TEs into a genome that already contains genes and TEs.

Example usage:
  python insert_tes.py \
      --genome genome.fasta \
      --TE new_TEs.fasta \
      --rate 1e-8 \
      --generations 10 \
      --bed existing_features.bed \
      --output my_output \
      --seed 42 \
      --fix_in 1e-6 \
      -b 1e-2 \
      -bf burn_in.txt \
      --TE_ratio TE_ratio.txt \
      --disable_genes
"""

import argparse
import bisect
import random
import re
import sys
import multiprocessing
import numpy as np  # New import for Poisson simulation
import os
import cutpaste_common as _cc

# Resolved relative to the working directory: the exciser writes this file and
# the next generation's inserter consumes it. PrinTE.sh always invokes both
# from the run directory, so the relative path is the cross-step contract.
_DEBT_FILE = "cutpaste_debt.tsv"

# Read-only TE library shared with worker processes via fork copy-on-write.
# Set in main() before the Pool is created; workers READ these (never mutate),
# so on Linux (fork) no private per-worker copy is pickled/made.
_TE_BY_CATEGORY = None
_TE_ALL_LIST = None

# Pre-computed translation table for reverse complement (C-speed)
_RC_TABLE = bytes.maketrans(b'ACGTacgtNn', b'TGCAtgcaNn')

def parse_args():
    parser = argparse.ArgumentParser(
        description="Randomly insert TE sequences into a genome, possibly nesting inside existing TEs. Output updated BED and FASTA."
    )
    parser.add_argument("--genome", required=True, help="Input genome FASTA file.")
    parser.add_argument("--TE", required=True, help="FASTA file of TEs to insert.")
    parser.add_argument("--rate", type=float, default=1e-8,
                        help="Rate of TE insertions per intact TE per generation. Default=1e-8")
    parser.add_argument("--generations", type=int, default=1,
                        help="Number of generations to simulate. Default=1")
    parser.add_argument("--bed", required=True,
                        help="Existing BED file with TE/gene coordinates.")
    parser.add_argument("--output", required=True, help="Output prefix (for .bed and .fasta).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility.")
    parser.add_argument("--fix_in", type=float, default=None,
                        help="Fixed rate of TE insertions per base per generation (overrides rate and birth_rate, if provided).")
    parser.add_argument("-b", "--birth_rate", type=float, default=0.0,
                        help="Birth rate of new TEs. Supports scientific (e.g. 1e-2) and numeric (e.g. 10) formats. Default=0.0")
    parser.add_argument("-bf", "--birth_file",
                        help="File with burn-in genome statistics. Must contain a line like: '... 49 TEs ...'")
    parser.add_argument("--TE_ratio", dest="TE_ratio_file",
                        help="File with TE category ratios. Format: <te_class> <te_superfamily> <non-normalized ratio> per line.")
    parser.add_argument("--euch_het_buffer", type=int, default=None,
                        help="Buffer (in bp) around gene features to be considered euchromatin. Only interpreted in rate mode.")
    parser.add_argument("--euch_het_bias", type=float, default=None,
                        help="Bias factor to increase probability of TE insertion in euchromatin. Only interpreted in rate mode.")
    parser.add_argument("-m", "--max_processes", type=int, default=1,
                        help="Number of chromosomes to process simultaneously. Default=1")
    parser.add_argument("--disable_genes", action="store_true",
                        help="Disable insertion into genes (only effective if --fix_in is provided).")
    parser.add_argument("--cutpaste_reinsertion", type=float, default=1.0,
                        help="Expected successful re-insertions per excised cut-and-paste "
                             "element. 1.0 conserves copy number (strict one-to-one "
                             "excision-reinsertion); <1.0 models failed re-insertion "
                             "(net loss); >1.0 models replication/repair-coupled "
                             "amplification (net gain). Default=1.0")
    args = parser.parse_args()
    if args.cutpaste_reinsertion < 0:
        parser.error("--cutpaste_reinsertion must be >= 0.")
    return args

def read_fasta(fasta_path):
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences[current_header] = "".join(current_seq)
    return sequences

def convert_genome_to_dict_of_lists(genome_dict):
    """Convert genome sequences to bytearrays for memory-efficient editing.

    Using bytearray instead of list-of-single-char-strings reduces memory
    ~50x (1 byte/base vs ~50 bytes/Python-str-object) and speeds up slice
    insertions ~8x (shifting bytes vs 8-byte pointers).
    """
    out = {}
    for chrom, seq in genome_dict.items():
        out[chrom] = bytearray(seq.encode('ascii'))
    return out

def convert_genome_back_to_fasta(genome_dict_of_lists):
    out = {}
    for chrom, seq_ba in genome_dict_of_lists.items():
        out[chrom] = seq_ba.decode('ascii')
    return out

def extract_te_info(header):
    match = re.match(r"[^#]+#([^/]+)/([^/~;]+)", header)
    if match:
        te_class = match.group(1).strip()
        te_superfamily = match.group(2).strip()
        return te_class, te_superfamily
    else:
        return "unknown", "unknown"

def build_te_index(te_raw):
    """Pre-build a {(class, superfamily): [(header, seq), ...]} lookup dict.

    Eliminates the O(num_TEs) regex scan that previously ran on every
    insertion call in pick_random_TE_by_category.
    """
    te_by_category = {}
    te_all_list = list(te_raw.items())
    for header, seq in te_all_list:
        key = extract_te_info(header)
        te_by_category.setdefault(key, []).append((header, seq))
    return te_by_category, te_all_list

def get_tsd_length(te_class, te_superfamily):
    """
    Determine the TSD length based on TE class and superfamily.
    Case-insensitive. Variable TSDs are sampled within ranges.
    """
    # Normalize for case-insensitive matching
    c = (te_class or "").strip().lower()
    s = (te_superfamily or "").strip().lower()

    # Helper lambdas for shared behaviors
    _const5        = lambda: 5
    _l1_rng        = lambda: random.randint(5, 20)   # LINE/L1
    _cacta_rng     = lambda: random.randint(2, 4)    # CACTA
    _hat_rng       = lambda: random.randint(5, 8)    # hAT
    _mule_rng      = lambda: random.randint(8, 9)    # MULE / MuDR
    _tourist_rng   = lambda: random.randint(2, 10)   # Tourist / Stowaway
    _unk_tir_rng   = lambda: random.randint(2, 10)   # Unknown TIR (explicit alias)

    tsd_mapping = {
        # Fixed TSDs (original)
        ('ltr', 'copia'):   _const5,
        ('ltr', 'gypsy'):   _const5,
        ('ltr', 'solo'):    _const5,
        ('ltr', 'ty3'):     _const5,
        ('ltr', 'unknown'): _const5,
        ('ltr', 'crm'):     _const5,
        ('ltr', 'trim'):    _const5,     # NEW alias

        ('dna', 'harbinger'): lambda: 3,
        ('dna', 'mariner'):   lambda: 2,

        ('dnaauto', 'helitron'): lambda: 0,
        ('dna',     'helitron'): lambda: 0,
        ('dnanona', 'helitron'): lambda: 0,

        # Variable TSDs (original)
        ('mite', 'stow'):     _tourist_rng,  # Stowaway 2–10
        ('mite', 'tourist'):  _tourist_rng,
        ('sine', 'trna'):     _l1_rng,
        ('sine', 'unknown'):  _l1_rng,
        ('line', 'l1'):       _l1_rng,
        ('line', 'unknown'):  _l1_rng,
        ('line', 'rte'):      _l1_rng,      # LINE/RTE

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

        # DTA uses same as hAT
        ('dna',  'dta'): _hat_rng,
        ('mite', 'dta'): _hat_rng,

        # DTC uses same as CACTA
        ('dna',  'dtc'): _cacta_rng,
        ('mite', 'dtc'): _cacta_rng,

        # DTH uses same as Tourist
        ('dna',  'dth'): _tourist_rng,
        ('mite', 'dth'): _tourist_rng,

        # DTM uses same as MULE
        ('dna',  'dtm'): _mule_rng,
        ('mite', 'dtm'): _mule_rng,

        # DTT uses same as Stowaway
        ('dna',  'dtt'): _tourist_rng,
        ('mite', 'dtt'): _tourist_rng,

        # DNAauto
        ('dnaauto', 'cactg'): _cacta_rng,   # CACTG  CACTA
        ('dnaauto', 'mle'):   _tourist_rng, # MLE  Stowaway
        ('dnaauto', 'pile'):  _unk_tir_rng, # unknown TIR (2–10)
        ('dnaauto', 'pole'):  _unk_tir_rng, # unknown TIR (2–10)

        # DNAnona
        ('dnanona', 'cactg'):   _cacta_rng,   # CACTG  CACTA
        ('dnanona', 'mle'):     _tourist_rng, # MLE  Stowaway
        ('dnanona', 'muletir'): _mule_rng,    # MULEtir  MULE
        ('dnanona', 'pile'):    _unk_tir_rng, # unknown TIR (2–10)
        ('dnanona', 'pole'):    _unk_tir_rng, # unknown TIR (2–10)
        ('dnanona', 'tourist'): _tourist_rng, # Tourist  Stowaway behavior
        ('dnanona', 'unknown'): _unk_tir_rng, # unknown TIR (2–10)
        ('dna', 'mutator'):     _mule_rng,     # Mutator/MuDR behaves like MULE (8–9)
        ('dna', 'tc1_mariner'): _tourist_rng,  # Tc1/Mariner behaves like Tourist/Stowaway (2–10)
        ('dna', 'unknown'):     _unk_tir_rng   # Unknown DNA → 2–10 bp (if desired, instead of default 5)

    }

    func = tsd_mapping.get((c, s))
    if func:
        return func()
    else:
        # Keep a helpful warning with the original casing
        print(f"Warning: TSD length not defined for TE class '{te_class}' and superfamily '{te_superfamily}'. Using default of 5.")
        return 5

def parse_bed(bed_path):
    features = []
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            chrom, start, end, name, tsd, strand = parts[:6]
            start = int(start)
            end = int(end)
            features.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': name,
                'strand': strand,
                'tsd': tsd
            })
    return features

def write_bed(features, out_bed):
    with open(out_bed, 'w') as f:
        for feat in features:
            # if this feature was just inserted, append a 7th column "new_TE"
            new_tag = '\tnew_TE' if feat.get('new', False) else ''
            f.write(
                f"{feat['chrom']}\t"
                f"{feat['start']}\t"
                f"{feat['end']}\t"
                f"{feat['name']}\t"
                f"{feat['tsd']}\t"
                f"{feat['strand']}"
                f"{new_tag}\n"
            )

def _write_fasta_record(fout, chrom, seq_list):
    """Write one FASTA record to a binary file object, no str-decode copy.

    `seq_list` is a bytearray/bytes written directly, so the bytes produced are
    identical to '>{chrom}\\n{seq}\\n' written in text mode on POSIX. The header
    is utf-8 encoded to match the original text-mode behavior on Linux (a
    default-locale text file), so non-ASCII contig names don't crash.
    """
    fout.write(b">")
    fout.write(chrom.encode("utf-8"))
    fout.write(b"\n")
    fout.write(seq_list)
    fout.write(b"\n")


def pick_random_TE(te_all_list):
    return random.choice(te_all_list)

def pick_random_TE_by_category(te_by_category, te_all_list, te_class, te_superfamily):
    matching = te_by_category.get((te_class, te_superfamily))
    if not matching:
        print(f"Warning: No TE found in TE library for category {te_class}/{te_superfamily}. Picking random from all TEs.")
        return random.choice(te_all_list)
    return random.choice(matching)

def reverse_complement(seq):
    return seq.encode('ascii')[::-1].translate(_RC_TABLE).decode('ascii')

def partial_name_match(name1, name2):
    if ';' not in name1 or ';' not in name2:
        return False
    parts1 = name1.split(';')
    parts2 = name2.split(';')
    if len(parts1) <= len(parts2) and parts1 == parts2[:len(parts1)]:
        return True
    elif len(parts2) < len(parts1) and parts2 == parts1[:len(parts2)]:
        return True
    return False

def get_intact_te_stats(features):
    """Return (intact_count, distribution_dict) in a single pass.

    Replaces the former count_intact_TE_count + get_intact_TE_distribution
    pair, which duplicated the same expensive logic.
    """
    n = len(features)
    intact_flags = [True] * n
    for i, feat in enumerate(features):
        feature_id = feat['name'].split(';')[0]
#       if 'gene' in feature_id.lower():
        if 'gene' in feature_id.lower() or "_FRAG" in feature_id:
            intact_flags[i] = False
    for i, feat in enumerate(features):
        feature_id = feat['name'].split(';')[0]
        if "_SOLO" in feature_id:
            intact_flags[i] = False
    for i, feat in enumerate(features):
        if intact_flags[i]:
            parts = feat['name'].split(';')
            if len(parts) > 1:
                if 'CUT_BY' in parts[1]:
                    intact_flags[i] = False
    for i, feat in enumerate(features):
        if not intact_flags[i]:
            continue
        parts = feat['name'].split(';')
        if len(parts) <= 1:
            continue
        window_start = max(0, i - 100)
        window_end = min(n, i + 101)
        for j in range(window_start, window_end):
            if j == i:
                continue
            other = features[j]
            if (feat['tsd'] == other['tsd'] and
                feat['strand'] == other['strand'] and
                partial_name_match(feat['name'], other['name'])):
                intact_flags[i] = False
                break
    count = 0
    distribution = {}
    for i, feat in enumerate(features):
        if intact_flags[i]:
            count += 1
            feature_id = feat['name'].split(';')[0]
            te_class, te_superfamily = extract_te_info(feature_id)
            key = (te_class, te_superfamily)
            distribution[key] = distribution.get(key, 0) + 1
    return count, distribution

# --- New helper functions for euchromatin/heterochromatin intervals ---
def merge_intervals(intervals):
    """Merge overlapping intervals. Each interval is a tuple (start, end)."""
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def compute_bias_intervals(genome, features, buffer_size, bias_factor):
    """
    For each chromosome, compute euchromatin intervals based on gene features with a buffer,
    merge overlapping intervals, then compute heterochromatin as the complement.
    Return a dict mapping chrom -> list of tuples (start, end, effective_weight),
    where effective_weight = length * (bias_factor for euchromatin, 1 for heterochromatin).
    """
    intervals = {}
    for chrom in genome:
        chrom_length = len(genome[chrom])
        euch_intervals = []
        # Pick gene features (using the same logic as intact TE count: feature id contains "gene")
        for feat in features:
            if feat['chrom'] != chrom:
                continue
            feature_id = feat['name'].split(';')[0]
            if "gene" in feature_id.lower():
                start = max(0, feat['start'] - buffer_size)
                end = min(chrom_length, feat['end'] + buffer_size)
                euch_intervals.append((start, end))
        euch_intervals = merge_intervals(euch_intervals)
        # Build heterochromatin intervals as the complement of euch_intervals.
        hetero_intervals = []
        prev_end = 0
        for (start, end) in euch_intervals:
            if start > prev_end:
                hetero_intervals.append((prev_end, start))
            prev_end = end
        if prev_end < chrom_length:
            hetero_intervals.append((prev_end, chrom_length))
        # Now combine intervals with effective weight.
        combined = []
        for (start, end) in euch_intervals:
            length = end - start
            if length > 0:
                combined.append((start, end, length * bias_factor))
        for (start, end) in hetero_intervals:
            length = end - start
            if length > 0:
                combined.append((start, end, length))
        # If there are no intervals at all (e.g. no genes), use full chromosome uniformly.
        if not combined:
            combined.append((0, chrom_length, chrom_length))
        intervals[chrom] = combined
    return intervals

def choose_weighted_position(chrom_length, intervals):
    """
    Given a list of intervals (start, end, weight), choose a position randomly weighted by weight.
    """
    total = sum(weight for (_, _, weight) in intervals)
    r = random.uniform(0, total)
    cumulative = 0
    for (start, end, weight) in intervals:
        cumulative += weight
        if r <= cumulative:
            # Choose a uniform position within this interval.
            return random.randint(start, end)
    # Fallback (should not happen)
    return random.randint(0, chrom_length)

def build_gene_index(features):
    """Build sorted gene intervals with max-end augmentation for bisect lookup.

    Returns (gene_intervals, gene_starts, max_ends) where max_ends[i] is the
    maximum end coordinate among gene_intervals[0..i].  This allows
    is_in_gene to terminate early even when genes overlap.
    """
    gene_intervals = []
    for feat in features:
        feature_id = feat['name'].split(';')[0]
        if "gene" in feature_id.lower():
            gene_intervals.append((feat['start'], feat['end']))
    gene_intervals.sort()
    gene_starts = [iv[0] for iv in gene_intervals]
    # Augment with running max of end positions for early termination
    max_ends = []
    running_max = 0
    for _start, end in gene_intervals:
        running_max = max(running_max, end)
        max_ends.append(running_max)
    return gene_intervals, gene_starts, max_ends

def is_in_gene(position, gene_intervals, gene_starts, max_ends):
    """Check if position falls within any gene region using binary search.

    Uses bisect on sorted gene starts, then walks backwards with max-end
    augmentation for correct early termination even with overlapping genes.
    O(log G) typical vs O(F) for the old linear scan over all features.
    """
    if not gene_intervals:
        return False
    idx = bisect.bisect_right(gene_starts, position)
    # Walk backwards: any gene starting at or before position might contain it
    for i in range(idx - 1, -1, -1):
        if gene_intervals[i][0] <= position < gene_intervals[i][1]:
            return True
        # If the max end up to this index doesn't reach position,
        # no earlier interval can contain position either
        if max_ends[i] <= position:
            break
    return False

# --- End of new helper functions ---

def process_chromosome(args_tuple):
    """
    Process all insertion events for one chromosome.
    This function is designed to be run in parallel for each chromosome.
    It applies insertion events sequentially on the chromosome's sequence and features.
    """
    (chrom, seq_list, features, events,
     bias_intervals, global_seed, disable_genes) = args_tuple
    # TE library is shared read-only via module globals (fork COW), not passed
    # per task -- avoids pickling a private copy into every worker.
    te_by_category = _TE_BY_CATEGORY
    te_all_list = _TE_ALL_LIST
    # Compute a deterministic seed using the global seed and the chromosome name.
    chrom_seed = (global_seed if global_seed is not None else 0) + sum(ord(c) for c in chrom)
    random.seed(chrom_seed)
    nested_count = 0
    non_nested_count = 0
    chrom_length = len(seq_list)
    te_seq_cache = {}  # Cache encoded TE bytearrays to avoid re-encoding
    # Process each insertion event assigned to this chromosome.
    for event in events:
        te_class_target, te_superfamily_target = event
        if chrom_length == 0:
            continue
        # Choose insertion position: if bias is active and intervals exist, pick weighted by euchromatin bias.
        # If disable_genes is active, re-sample if the candidate falls within a gene.
        if bias_intervals:
            candidate = choose_weighted_position(chrom_length, bias_intervals)
        else:
            candidate = random.randint(0, chrom_length)
        if disable_genes:
            gene_intervals, gene_starts, max_ends = build_gene_index(features)
            attempt = 0
            max_attempts = 1000
            while is_in_gene(candidate, gene_intervals, gene_starts, max_ends) and attempt < max_attempts:
                if bias_intervals:
                    candidate = choose_weighted_position(chrom_length, bias_intervals)
                else:
                    candidate = random.randint(0, chrom_length)
                attempt += 1
            # If a gene-free position was not found after many attempts, proceed with the candidate.
        insertion_pos = candidate

        te_header, te_seq = pick_random_TE_by_category(
            te_by_category, te_all_list, te_class_target, te_superfamily_target)
        te_class, te_superfamily = extract_te_info(te_header)
        tsd_length = get_tsd_length(te_class, te_superfamily)
        strand = "+" if random.random() < 0.5 else "-"
        # Get TE sequence as bytearray, using cache to avoid re-encoding
        if strand == "-":
            cache_key = (te_header, '-')
            if cache_key not in te_seq_cache:
                te_seq_cache[cache_key] = bytearray(reverse_complement(te_seq).encode('ascii'))
            te_seq_bytes = te_seq_cache[cache_key]
        else:
            cache_key = (te_header, '+')
            if cache_key not in te_seq_cache:
                te_seq_cache[cache_key] = bytearray(te_seq.encode('ascii'))
            te_seq_bytes = te_seq_cache[cache_key]
        if tsd_length > insertion_pos:
            print(f"Warning: TSD length {tsd_length} is greater than insertion position {insertion_pos} on {chrom}. Reducing.")
            tsd_length = insertion_pos
        tsd_seq_bytes = bytearray()
        if tsd_length > 0:
            tsd_seq_bytes = seq_list[insertion_pos - tsd_length:insertion_pos]
        tsd_string = tsd_seq_bytes.decode('ascii') if tsd_length > 0 else "NA"
        # Insert TE sequence into the chromosome sequence.
        seq_list[insertion_pos:insertion_pos] = te_seq_bytes
        # If TSD exists, duplicate it after the inserted TE.
        if tsd_length > 0:
            seq_list[insertion_pos + len(te_seq_bytes):insertion_pos + len(te_seq_bytes)] = tsd_seq_bytes
        shift_amount = len(te_seq_bytes) + tsd_length

        te_start = insertion_pos
        te_end   = insertion_pos + len(te_seq_bytes)
        new_te_name = te_header
        new_te_feature = {
            'chrom': chrom,
            'start': te_start,
            'end':   te_end,
            'name':  new_te_name,
            'strand': strand,
            'tsd':   tsd_string if tsd_length > 0 else 'NA',
        }

        new_te_feature['new'] = True

        updated_features = []
        nesting_happened = False
        for feat in features:
            if feat['end'] <= insertion_pos:
                updated_features.append(feat)
            elif feat['start'] >= insertion_pos:
                feat['start'] += shift_amount
                feat['end']   += shift_amount
                updated_features.append(feat)
            else:
                # Feature spans the insertion point; split it.
                cut_feat_end = insertion_pos
                left_piece = {
                    'chrom': chrom,
                    'start': feat['start'],
                    'end':   cut_feat_end,
                    'name':  feat['name'] + f";CUT_BY:{new_te_name}",
                    'strand': feat['strand'],
                    'tsd':  feat['tsd'],
                }
                right_length = feat['end'] - cut_feat_end
                right_piece_start = insertion_pos + shift_amount
                right_piece_end   = right_piece_start + right_length
                if right_length > 0:
                    right_piece = {
                        'chrom': chrom,
                        'start': right_piece_start,
                        'end':   right_piece_end,
                        'name':  feat['name'] + f";CUT_BY:{new_te_name}",
                        'strand': feat['strand'],
                        'tsd':  feat['tsd'],
                    }
                    updated_features.append(left_piece)
                    updated_features.append(right_piece)
                else:
                    updated_features.append(left_piece)
                new_te_feature['name'] += f";NESTED_IN:{feat['name']}"
                nesting_happened = True
        features = updated_features
        features.append(new_te_feature)

        print(f"Insertion on {chrom}: Inserted TE '{te_header}' at position {insertion_pos}, strand={strand}, TSD={tsd_string}")
        if nesting_happened:
            print(f"  => Nested inside an existing feature. TE name now: {new_te_feature['name']}")
            nested_count += 1
        else:
            non_nested_count += 1
        chrom_length = len(seq_list)  # update chromosome length after insertion

    # Return the updated sequence, features, and counts for this chromosome.
    return (chrom, seq_list, features, nested_count, non_nested_count)

def main():
    args = parse_args()

    # In fix_in mode, simply ignore the euchromatin/heterochromatin bias parameters if provided.
    if args.fix_in is not None:
        if args.euch_het_buffer is not None or args.euch_het_bias is not None:
            print("Warning: --euch_het_buffer and --euch_het_bias are ignored in fix_in mode.")

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)  # Ensure numpy randomness is reproducible

    print("Reading genome FASTA ...")
    genome_raw = read_fasta(args.genome)
    genome_size = sum(len(seq) for seq in genome_raw.values())
    print(f"Genome size: {genome_size} bases")

    print("Converting genome to editable bytearrays ...")
    genome = convert_genome_to_dict_of_lists(genome_raw)
    del genome_raw  # original strings no longer needed; lengths read from `genome`

    print("Reading TE FASTA ...")
    te_raw = read_fasta(args.TE)

    print("Building TE category index ...")
    te_by_category, te_all_list = build_te_index(te_raw)
    # Expose the read-only TE library to worker processes via fork copy-on-write.
    global _TE_BY_CATEGORY, _TE_ALL_LIST
    _TE_BY_CATEGORY = te_by_category
    _TE_ALL_LIST = te_all_list

    print("Reading BED file ...")
    features_all = parse_bed(args.bed)
    features_original = features_all[:]  # keep a copy

    intact_TE_count, intact_distribution = get_intact_te_stats(features_original)

    cutpaste_set = _cc.load_cutpaste_set(args.TE_ratio_file)
    if cutpaste_set:
        _excluded = 0
        for key in list(intact_distribution):
            if key in cutpaste_set:
                intact_TE_count -= intact_distribution[key]
                del intact_distribution[key]
                _excluded += 1
        print(f"Cut-and-paste: excluded {_excluded} family(ies) present in "
              f"the genome from replicative insertion "
              f"({len(cutpaste_set)} flagged in ratios).")
        if args.fix_in is None and intact_TE_count <= 0:
            print("Warning: all intact TEs in the genome belong to "
                  "cut-and-paste families; replicative insertion rate is 0 "
                  "this generation (relocation debt, if any, still applies).")

    print(f"Total number of intact TEs from BED: {intact_TE_count}")

    print("Distribution of intact TEs by (te_class, te_superfamily):")
    for key, count in intact_distribution.items():
        print(f"  {key[0]}/{key[1]}: {count}")

    # --- Calculate total TE insertions using a Poisson distribution ---
    if args.fix_in is not None:
        lambda_value = args.fix_in * genome_size * args.generations
        total_insertions = np.random.poisson(lambda_value)
        print(f"Using fixed insertion rate: {args.fix_in} per base per generation")
        print(f"Lambda (expected insertions): {lambda_value}")
        print(f"Simulated total TE insertions to perform: {total_insertions}")
    else:
        lambda_value = args.rate * intact_TE_count * args.generations
        total_insertions = np.random.poisson(lambda_value)
        print(f"Rate (per intact TE per generation): {args.rate}")
        print(f"Generations: {args.generations}")
        print(f"Lambda (expected insertions): {lambda_value}")
        print(f"Simulated total TE insertions to perform: {total_insertions}")

    # Build insertion events.
    if args.fix_in is not None:
        if not args.TE_ratio_file:
            print("Error: --TE_ratio file must be provided if --fix_in is used.")
            sys.exit(1)
        te_ratio = {}
        try:
            with open(args.TE_ratio_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    te_class, te_superfamily, weight_str = parts[:3]
                    try:
                        weight = float(weight_str)
                    except:
                        weight = 0.0
                    te_ratio[(te_class, te_superfamily)] = weight
        except Exception as e:
            print(f"Error reading TE_ratio_file: {e}")
            sys.exit(1)
        if not te_ratio:
            print("Error: No TE ratio information loaded from TE_ratio_file.")
            sys.exit(1)
        for _k in list(te_ratio):
            if _k in cutpaste_set:
                del te_ratio[_k]
        if not te_ratio:
            print("Error: all TE ratio families are cut-and-paste; "
                  "no replicative families to insert in fix_in mode.")
            sys.exit(1)
        total_weight = sum(te_ratio.values())
        for k in te_ratio:
            te_ratio[k] /= total_weight
        insertion_events = []
        for i in range(total_insertions):
            chosen_cat = random.choices(list(te_ratio.keys()), weights=list(te_ratio.values()), k=1)[0]
            insertion_events.append(chosen_cat)
    else:
        insertion_events = []
        intact_total = sum(intact_distribution.values())
        for cat, count in intact_distribution.items():
            insertion_count = round(total_insertions * (count / intact_total))
            for i in range(insertion_count):
                insertion_events.append(cat)
        diff = total_insertions - len(insertion_events)
        if diff > 0:
            categories = list(intact_distribution.keys())
            weights = [intact_distribution[cat] for cat in categories]
            for i in range(diff):
                chosen_cat = random.choices(categories, weights=weights)[0]
                insertion_events.append(chosen_cat)
        elif diff < 0:
            for i in range(-diff):
                removal_index = random.randrange(len(insertion_events))
                insertion_events.pop(removal_index)

    # --- New TE Birth functionality ---
    if args.fix_in is None and args.birth_rate and args.birth_file and args.TE_ratio_file:
        try:
            with open(args.birth_file, 'r') as bf:
                line = bf.readline().strip()
            m = re.search(r'(\d+)\s+TEs', line)
            if m:
                number_initial_tes = int(m.group(1))
            else:
                print("Error: Could not extract number of initial TEs from birth_file.")
                sys.exit(1)
        except Exception as e:
            print(f"Error reading birth_file: {e}")
            sys.exit(1)
        number_of_born_tes = int(round(args.generations * args.birth_rate * number_initial_tes))
        print(f"Number of born TEs to insert (from birth_rate and birth_file): {number_of_born_tes}")

        te_ratio_birth = {}
        try:
            with open(args.TE_ratio_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    te_class, te_superfamily, weight_str = parts[:3]
                    try:
                        weight = float(weight_str)
                    except:
                        weight = 0.0
                    te_ratio_birth[(te_class, te_superfamily)] = weight
        except Exception as e:
            print(f"Error reading TE_ratio_file: {e}")
            sys.exit(1)
        if not te_ratio_birth:
            print("Error: No TE ratio information loaded from TE_ratio_file.")
            sys.exit(1)
        total_weight_birth = sum(te_ratio_birth.values())
        for k in te_ratio_birth:
            te_ratio_birth[k] /= total_weight_birth
        born_categories = list(te_ratio_birth.keys())
        born_weights = [te_ratio_birth[cat] for cat in born_categories]
        for i in range(number_of_born_tes):
            chosen_cat = random.choices(born_categories, weights=born_weights, k=1)[0]
            insertion_events.append(chosen_cat)
        print(f"Total insertion events after adding born TEs: {len(insertion_events)}")
    # --- End New TE Birth functionality ---

    # --- Consume cut-and-paste relocation debt from prior gen ---
    # Each excised element yields on average --cutpaste_reinsertion successful
    # re-insertions: 1.0 conserves copy number, <1.0 nets loss (failed
    # re-insertion), >1.0 nets gain (replication/repair-coupled amplification).
    # The fractional remainder is resolved by one stochastic-rounding draw per
    # family; the frac > 0 guard keeps the random stream untouched at 1.0.
    if os.path.exists(_DEBT_FILE):
        debt = _cc.read_debt(_DEBT_FILE)
        excised = added = skipped = 0
        for fam, cnt in debt.items():
            if fam in te_by_category:
                excised += cnt
                expected = cnt * args.cutpaste_reinsertion
                n = int(expected)
                frac = expected - n
                if frac > 0 and random.random() < frac:
                    n += 1
                insertion_events.extend([fam] * n)
                added += n
            else:
                skipped += cnt
        os.remove(_DEBT_FILE)
        print(f"Cut-and-paste relocation: {excised} excised element(s) -> "
              f"{added} re-inserted (reinsertion ratio "
              f"{args.cutpaste_reinsertion}); skipped {skipped} (family "
              f"extinct in library). Consumed and removed {_DEBT_FILE}.")

    random.shuffle(insertion_events)

    # --- Pre-calculate insertion position intervals if euchromatin bias is in effect (only in rate mode) ---
    use_bias = (args.fix_in is None and args.euch_het_buffer is not None and args.euch_het_bias is not None)
    bias_intervals_all = None
    if use_bias:
        print("Computing euchromatin/heterochromatin intervals based on gene features and buffer ...")
        bias_intervals_all = compute_bias_intervals(genome, features_original, args.euch_het_buffer, args.euch_het_bias)
    # --- End interval pre-calculation ---

    # Distribute insertion events among chromosomes.
    # Use the original chromosome lengths from genome_raw.
    chroms = list(genome.keys())
    chrom_lengths = [len(genome[c]) for c in chroms]
    event_assignment = {chrom: [] for chrom in chroms}
    for event in insertion_events:
        chosen_chrom = random.choices(chroms, weights=chrom_lengths, k=1)[0]
        event_assignment[chosen_chrom].append(event)

    # Partition features by chromosome.
    features_by_chrom = {chrom: [] for chrom in chroms}
    for feat in features_original:
        if feat['chrom'] in features_by_chrom:
            features_by_chrom[feat['chrom']].append(feat)
        else:
            features_by_chrom[feat['chrom']] = [feat]

    # Determine if disable_genes should be active (only if fix_in is provided).
    disable_genes_flag = args.disable_genes if args.fix_in is not None else False

    # Prepare arguments for parallel processing.
    tasks = []
    for chrom in chroms:
        seq_list = genome[chrom]
        feats = features_by_chrom.get(chrom, [])
        events = event_assignment.get(chrom, [])
        bias_intervals = None
        if use_bias and chrom in bias_intervals_all:
            bias_intervals = bias_intervals_all[chrom]
        tasks.append((chrom, seq_list, feats, events,
                       bias_intervals, args.seed, disable_genes_flag))

    # Stream output in sorted-chromosome order so peak memory never holds more
    # than one finished chromosome at a time (plus the pool's in-flight results):
    # write each chromosome's FASTA as it completes, then free it -- no whole-
    # genome `results` list and no decoded `final_genome` copy. Sorting tasks by
    # chromosome makes the ordered imap / serial map yield in the same sorted
    # order the previous implementation wrote in, so the FASTA is byte-identical.
    # Per-chromosome RNG is seeded from the chromosome name (process_chromosome),
    # so submission order never changes results.
    tasks.sort(key=lambda t: t[0])

    all_features = []
    total_nested = 0
    total_non_nested = 0
    out_fasta = f"{args.output}.fasta"
    out_bed = f"{args.output}.bed"
    tmp_fasta = out_fasta + ".tmp"
    tmp_bed = out_bed + ".tmp"
    print(f"Processing {len(chroms)} chromosomes using up to {args.max_processes} processes...")
    print(f"Writing updated FASTA to {out_fasta}")
    # Write to temp files and os.replace() only after BOTH finish, so a mid-run
    # worker failure can never leave a truncated .fasta/.bed that looks complete
    # (the streaming write would otherwise expose a partial file to a resume).
    with open(tmp_fasta, "wb") as fout:
        if args.max_processes > 1:
            # fork context: workers inherit the read-only TE-library globals via
            # copy-on-write instead of pickling a private copy in every task.
            # The COW-globals design REQUIRES fork (Linux); guard so a fork-less
            # platform fails loudly instead of workers silently seeing None.
            if "fork" not in multiprocessing.get_all_start_methods():
                sys.exit("Error: parallel mode (-m >1) requires the 'fork' start "
                         "method (Linux). Re-run with -m 1 on this platform.")
            ctx = multiprocessing.get_context("fork")
            # chunksize mirrors pool.map's auto-batching so a many-scaffold genome
            # doesn't pay per-task IPC overhead (imap defaults to chunksize=1).
            chunk = max(1, len(tasks) // (4 * args.max_processes))
            with ctx.Pool(processes=args.max_processes) as pool:
                for chrom, seq_list, feats, nested_count, non_nested_count in \
                        pool.imap(process_chromosome, tasks, chunksize=chunk):
                    _write_fasta_record(fout, chrom, seq_list)
                    all_features.extend(feats)
                    total_nested += nested_count
                    total_non_nested += non_nested_count
                    del seq_list
        else:
            # Serial: process_chromosome mutates the genome bytearray in place and
            # returns that same object, still bound by genome[chrom] and tasks[i].
            # Null both after writing so each finished chromosome is freed instead
            # of the whole (possibly ballooned) genome staying resident to the end.
            for i in range(len(tasks)):
                chrom, seq_list, feats, nested_count, non_nested_count = \
                    process_chromosome(tasks[i])
                _write_fasta_record(fout, chrom, seq_list)
                all_features.extend(feats)
                total_nested += nested_count
                total_non_nested += non_nested_count
                seq_list = None
                genome[chrom] = None
                tasks[i] = None

    all_features.sort(key=lambda x: (x['chrom'], x['start']))
    print(f"Writing updated BED to {out_bed}")
    write_bed(all_features, tmp_bed)
    os.replace(tmp_fasta, out_fasta)
    os.replace(tmp_bed, out_bed)
    total_complete = total_nested + total_non_nested
    print(f"Total TE insertions performed: {total_complete} (Nested: {total_nested}, Non-nested: {total_non_nested})")
    print("Done.")

if __name__ == "__main__":
    main()
