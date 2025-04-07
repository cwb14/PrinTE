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
      --euch_het_buffer 1000 \
      --euch_het_bias 1.1 \
      -m 4 \
      --disable_genes
"""

import argparse
import random
import re
import sys
import multiprocessing

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
                        help="Buffer (in bp) around gene features to be considered euchromatin.")
    parser.add_argument("--euch_het_bias", type=float, default=None,
                        help="Bias factor to increase probability of TE insertion in euchromatin.")
    parser.add_argument("-m", "--max_processes", type=int, default=1,
                        help="Number of chromosomes to process simultaneously. Default=1")
    parser.add_argument("--disable_genes", action="store_true",
                        help="Disable insertion into genes (only effective if --fix_in is provided).")
    return parser.parse_args()

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
    out = {}
    for chrom, seq in genome_dict.items():
        out[chrom] = list(seq)
    return out

def convert_genome_back_to_fasta(genome_dict_of_lists):
    out = {}
    for chrom, seq_list in genome_dict_of_lists.items():
        out[chrom] = "".join(seq_list)
    return out

def extract_te_info(header):
    match = re.match(r"[^#]+#([^/]+)/([^~;]+)", header)
    if match:
        te_class = match.group(1).strip()
        te_superfamily = match.group(2).strip()
        return te_class, te_superfamily
    else:
        return "unknown", "unknown"

def get_tsd_length(te_class, te_superfamily):
    import random
    tsd_mapping = {
        ('LTR', 'Copia'): lambda: 5,
        ('LTR', 'Gypsy'): lambda: 5,
        ('LTR', 'Solo'):  lambda: 5,
        ('LTR', 'Ty3'):   lambda: 5,
        ('LTR', 'unknown'): lambda: 5,

        ('DNA', 'Harbinger'): lambda: 3,
        ('DNA', 'Mariner'):   lambda: 2,

        ('DNAauto', 'Helitron'): lambda: 0,
        ('DNA', 'Helitron'):     lambda: 0,
        ('DNAnona', 'Helitron'): lambda: 0,

        ('MITE', 'Stow'):    lambda: random.randint(2, 10),
        ('MITE', 'Tourist'): lambda: random.randint(2, 10),

        ('SINE', 'tRNA'):    lambda: random.randint(5, 20),
        ('SINE', 'unknown'): lambda: random.randint(5, 20),
        ('LINE', 'L1'):      lambda: random.randint(5, 20),
        ('LINE', 'unknown'): lambda: random.randint(5, 20),

        ('DNAauto', 'CACTA'): lambda: random.randint(2, 4),
        ('DNA', 'CACTA'):     lambda: random.randint(2, 4),
        ('DNAnona', 'CACTA'): lambda: random.randint(2, 4),

        ('DNA', 'hAT'):     lambda: random.randint(5, 8),
        ('DNAnona', 'hAT'): lambda: random.randint(5, 8),
        ('DNAauto', 'hAT'): lambda: random.randint(5, 8),

        ('DNAauto', 'MuDR'): lambda: random.randint(8, 9),
        ('DNAnona', 'MULE'): lambda: random.randint(8, 9),
        ('DNA', 'MuDR'):     lambda: random.randint(8, 9),
        ('DNAauto', 'MULE'): lambda: random.randint(8, 9),
    }
    key = (te_class, te_superfamily)
    if key in tsd_mapping:
        return tsd_mapping[key]()
    else:
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
            # Write TSD in the 5th column and strand in the 6th column.
            f.write(f"{feat['chrom']}\t{feat['start']}\t{feat['end']}\t"
                    f"{feat['name']}\t{feat['tsd']}\t{feat['strand']}\n")

def pick_random_TE(te_dict):
    headers = list(te_dict.keys())
    choice = random.choice(headers)
    return choice, te_dict[choice]

def pick_random_TE_by_category(te_dict, te_class, te_superfamily):
    matching = [(h, s) for h, s in te_dict.items() if extract_te_info(h) == (te_class, te_superfamily)]
    if not matching:
        print(f"Warning: No TE found in TE library for category {te_class}/{te_superfamily}. Picking random from all TEs.")
        return pick_random_TE(te_dict)
    return random.choice(matching)

def reverse_complement(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','g':'c','c':'g','N':'N','n':'n'}
    rc = []
    for base in reversed(seq):
        rc.append(comp.get(base, base))
    return "".join(rc)

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

def count_intact_TE_count(features):
    n = len(features)
    intact_flags = [True] * n
    for i, feat in enumerate(features):
        feature_id = feat['name'].split(';')[0]
        if 'gene' in feature_id.lower():
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
    count = sum(1 for flag in intact_flags if flag)
    return count

def get_intact_TE_distribution(features):
    n = len(features)
    intact_flags = [True] * n
    for i, feat in enumerate(features):
        feature_id = feat['name'].split(';')[0]
        if 'gene' in feature_id.lower():
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
    distribution = {}
    for i, feat in enumerate(features):
        if intact_flags[i]:
            feature_id = feat['name'].split(';')[0]
            te_class, te_superfamily = extract_te_info(feature_id)
            key = (te_class, te_superfamily)
            distribution[key] = distribution.get(key, 0) + 1
    return distribution

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

def is_in_gene(position, features):
    """
    Check if the given position falls within any gene region in features.
    A feature is considered a gene if 'gene' is in its first name field.
    """
    for feat in features:
        feature_id = feat['name'].split(';')[0]
        if "gene" in feature_id.lower():
            if feat['start'] <= position < feat['end']:
                return True
    return False

# --- End of new helper functions ---

def process_chromosome(args_tuple):
    """
    Process all insertion events for one chromosome.
    This function is designed to be run in parallel for each chromosome.
    It applies insertion events sequentially on the chromosome’s sequence and features.
    """
    (chrom, seq_list, features, events, te_raw, bias_intervals, global_seed, disable_genes) = args_tuple
    # Compute a deterministic seed using the global seed and the chromosome name.
    chrom_seed = (global_seed if global_seed is not None else 0) + sum(ord(c) for c in chrom)
    random.seed(chrom_seed)
    nested_count = 0
    non_nested_count = 0
    chrom_length = len(seq_list)
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
            attempt = 0
            max_attempts = 1000
            while is_in_gene(candidate, features) and attempt < max_attempts:
                if bias_intervals:
                    candidate = choose_weighted_position(chrom_length, bias_intervals)
                else:
                    candidate = random.randint(0, chrom_length)
                attempt += 1
            # If a gene-free position was not found after many attempts, proceed with the candidate.
        insertion_pos = candidate

        te_header, te_seq = pick_random_TE_by_category(te_raw, te_class_target, te_superfamily_target)
        te_class, te_superfamily = extract_te_info(te_header)
        tsd_length = get_tsd_length(te_class, te_superfamily)
        strand = "+" if random.random() < 0.5 else "-"
        if strand == "-":
            te_seq = reverse_complement(te_seq)
        te_seq_list = list(te_seq)
        if tsd_length > insertion_pos:
            print(f"Warning: TSD length {tsd_length} is greater than insertion position {insertion_pos} on {chrom}. Reducing.")
            tsd_length = insertion_pos
        tsd_seq_list = []
        if tsd_length > 0:
            tsd_seq_list = seq_list[insertion_pos - tsd_length:insertion_pos]
        tsd_string = "".join(tsd_seq_list) if tsd_length > 0 else "NA"
        # Insert TE sequence into the chromosome sequence.
        seq_list[insertion_pos:insertion_pos] = te_seq_list
        # If TSD exists, duplicate it after the inserted TE.
        if tsd_length > 0:
            seq_list[insertion_pos + len(te_seq_list):insertion_pos + len(te_seq_list)] = tsd_seq_list
        shift_amount = len(te_seq_list) + tsd_length

        te_start = insertion_pos
        te_end   = insertion_pos + len(te_seq_list)
        new_te_name = te_header
        new_te_feature = {
            'chrom': chrom,
            'start': te_start,
            'end':   te_end,
            'name':  new_te_name,
            'strand': strand,
            'tsd':   tsd_string if tsd_length > 0 else 'NA',
        }

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
    if args.seed is not None:
        random.seed(args.seed)

    print("Reading genome FASTA ...")
    genome_raw = read_fasta(args.genome)
    genome_size = sum(len(seq) for seq in genome_raw.values())
    print(f"Genome size: {genome_size} bases")

    print("Converting genome to editable lists ...")
    genome = convert_genome_to_dict_of_lists(genome_raw)

    print("Reading TE FASTA ...")
    te_raw = read_fasta(args.TE)

    print("Reading BED file ...")
    features_all = parse_bed(args.bed)
    features_original = features_all[:]  # keep a copy

    intact_TE_count = count_intact_TE_count(features_original)
    print(f"Total number of intact TEs from BED: {intact_TE_count}")

    intact_distribution = get_intact_TE_distribution(features_original)
    print("Distribution of intact TEs by (te_class, te_superfamily):")
    for key, count in intact_distribution.items():
        print(f"  {key[0]}/{key[1]}: {count}")

    # Determine total TE insertions to perform.
    if args.fix_in is not None:
        total_insertions = int(args.fix_in * genome_size * args.generations)
        print(f"Using fixed insertion rate: {args.fix_in} per base per generation")
        print(f"Calculated total TE insertions to perform: {total_insertions}")
    else:
        total_insertions = int(args.rate * intact_TE_count * args.generations)
        print(f"Rate (per intact TE per generation): {args.rate}")
        print(f"Generations: {args.generations}")
        print(f"Total TE insertions to perform: {total_insertions}")

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

    random.shuffle(insertion_events)

    # --- Pre-calculate insertion position intervals if euchromatin bias is in effect ---
    use_bias = (args.euch_het_buffer is not None and args.euch_het_bias is not None)
    bias_intervals_all = None
    if use_bias:
        print("Computing euchromatin/heterochromatin intervals based on gene features and buffer ...")
        bias_intervals_all = compute_bias_intervals(genome, features_original, args.euch_het_buffer, args.euch_het_bias)
    # --- End interval pre-calculation ---

    # Distribute insertion events among chromosomes.
    # Use the original chromosome lengths from genome_raw.
    chroms = list(genome.keys())
    chrom_lengths = [len(genome_raw[c]) for c in chroms]
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
        tasks.append((chrom, seq_list, feats, events, te_raw, bias_intervals, args.seed, disable_genes_flag))

    print(f"Processing {len(chroms)} chromosomes using up to {args.max_processes} processes...")
    if args.max_processes > 1:
        with multiprocessing.Pool(processes=args.max_processes) as pool:
            results = pool.map(process_chromosome, tasks)
    else:
        results = list(map(process_chromosome, tasks))

    # Merge the processed chromosomes.
    final_genome = {}
    all_features = []
    total_nested = 0
    total_non_nested = 0
    for chrom, seq_list, feats, nested_count, non_nested_count in results:
        final_genome[chrom] = "".join(seq_list)
        all_features.extend(feats)
        total_nested += nested_count
        total_non_nested += non_nested_count

    all_features.sort(key=lambda x: (x['chrom'], x['start']))
    out_bed = f"{args.output}.bed"
    print(f"Writing updated BED to {out_bed}")
    write_bed(all_features, out_bed)
    out_fasta = f"{args.output}.fasta"
    print(f"Writing updated FASTA to {out_fasta}")
    with open(out_fasta, 'w') as f:
        for chrom in sorted(final_genome.keys()):
            f.write(f">{chrom}\n")
            f.write(final_genome[chrom] + "\n")
    total_complete = total_nested + total_non_nested
    print(f"Total TE insertions performed: {total_complete} (Nested: {total_nested}, Non-nested: {total_non_nested})")
    print("Done.")

if __name__ == "__main__":
    main()
