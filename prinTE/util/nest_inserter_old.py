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
      --fix_in 800 \
      -b 1e-2 \
      -bf burn_in.txt \
      --TE_ratio TE_ratio.txt
"""

import argparse
import random
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Randomly insert TE sequences into a genome, possibly nesting inside existing TEs. Output updated BED and FASTA."
    )
    parser.add_argument("--genome", required=True, help="Input genome FASTA file.")
    parser.add_argument("--TE", required=True, help="FASTA file of TEs to insert.")
    parser.add_argument("--rate", type=float, default=1e-8,
                        help="Rate of TE insertions per base per generation. Default=1e-8")
    parser.add_argument("--generations", type=int, default=1,
                        help="Number of generations to simulate. Default=1")
    parser.add_argument("--bed", required=True,
                        help="Existing BED file with TE/gene coordinates.")
    parser.add_argument("--output", required=True, help="Output prefix (for .bed and .fasta).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility.")
    # New optional parameter: if provided, use this fixed number for total insertions.
    parser.add_argument("--fix_in", type=int, default=None,
                        help="Fixed total number of TE insertions to perform (overrides rate and generations if provided).")
    # New parameters for TE birth functionality.
    parser.add_argument("-b", "--birth_rate", type=float, default=0.0,
                        help="Birth rate of new TEs. Supports scientific (e.g. 1e-2) and numeric (e.g. 10) formats. Default=0.0")
    parser.add_argument("-bf", "--birth_file",
                        help="File with burn-in genome statistics. Must contain a line like: '... 49 TEs ...'")
    parser.add_argument("--TE_ratio", dest="TE_ratio_file",
                        help="File with TE category ratios. Format: <te_class> <te_superfamily> <non-normalized ratio> per line.")
    return parser.parse_args()

def read_fasta(fasta_path):
    """
    Read a FASTA file into a dict: {header: sequence_string} if multiple sequences,
    or if it's a genome, we might keep {chrom: list_of_chars} to simplify insertions.
    This function attempts to handle multi-sequence FASTA.
    """
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous if exists
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].strip()  # drop '>'
                current_seq = []
            else:
                current_seq.append(line)
        # final
        if current_header is not None:
            sequences[current_header] = "".join(current_seq)

    return sequences

def convert_genome_to_dict_of_lists(genome_dict):
    """
    Convert {header: seq_string} -> {header: list_of_chars}
    so that we can more easily insert at arbitrary positions.
    """
    out = {}
    for chrom, seq in genome_dict.items():
        out[chrom] = list(seq)
    return out

def convert_genome_back_to_fasta(genome_dict_of_lists):
    """
    Convert {header: list_of_chars} back to {header: seq_string}.
    """
    out = {}
    for chrom, seq_list in genome_dict_of_lists.items():
        out[chrom] = "".join(seq_list)
    return out

def extract_te_info(header):
    """
    Given a TE header in the style "TEname#TEclass/TEsuperfamily"
    or "TEname#TEclass/TEsuperfamily;SomeExtraStuff" or
    with metadata appended using '~', e.g.
    "ATCOPIA31#LTR/Copia~LTRlen:271", return (te_class, te_superfamily).
    """
    # The regex looks for a '#' then captures characters until a '/' then
    # captures characters until a '~' (if present) or end-of-string.
    match = re.match(r"[^#]+#([^/]+)/([^~;]+)", header)
    if match:
        te_class = match.group(1).strip()
        te_superfamily = match.group(2).strip()
        return te_class, te_superfamily
    else:
        return "unknown", "unknown"

def get_tsd_length(te_class, te_superfamily):
    """
    Determine the TSD length based on TE class and superfamily.
    If TSD length is variable, selects a random length within the specified range.
    """
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
    """
    Parse the BED file. Return a list of dicts:
    [
      {
        'chrom': ..., 'start': ..., 'end': ...,
        'name': ..., 'strand': ..., 'tsd': ...
      },
      ...
    ]
    The current input BED format is assumed to be:
      <chromosome> <start> <end> <feature> <TSD> <strand>
    which we convert internally to:
      'chrom': <chromosome>,
      'start': <start>,
      'end': <end>,
      'name': <feature>,
      'strand': <strand>,
      'tsd': <TSD>
    All positions are assumed 0-based.
    """
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
    """
    Write out the features in BED format.
    The output format is:
      <chromosome> <start> <end> <feature> <strand> <TSD>
    """
    with open(out_bed, 'w') as f:
        for feat in features:
            f.write(f"{feat['chrom']}\t{feat['start']}\t{feat['end']}\t"
                    f"{feat['name']}\t{feat['strand']}\t{feat['tsd']}\n")

def pick_random_TE(te_dict):
    """
    Given a dict of TE_name -> TE_sequence,
    pick a random key and return (header, sequence).
    """
    headers = list(te_dict.keys())
    choice = random.choice(headers)
    return choice, te_dict[choice]

def pick_random_TE_by_category(te_dict, te_class, te_superfamily):
    """
    Given a dict of TE_name -> TE_sequence and a desired te_class and te_superfamily,
    return a random (header, sequence) whose header (parsed with extract_te_info)
    matches the desired te_class and te_superfamily.
    If no TE is found for that category, warn and return any random TE.
    """
    matching = [(h, s) for h, s in te_dict.items() if extract_te_info(h) == (te_class, te_superfamily)]
    if not matching:
        print(f"Warning: No TE found in TE library for category {te_class}/{te_superfamily}. Picking random from all TEs.")
        return pick_random_TE(te_dict)
    return random.choice(matching)

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    """
    comp = {'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','g':'c','c':'g','N':'N','n':'n'}
    rc = []
    for base in reversed(seq):
        rc.append(comp.get(base, base))
    return "".join(rc)

def partial_name_match(name1, name2):
    """
    Check if one name is a partial match of the other.
    Both names must have at least one supplemental info field (i.e. contain a semicolon).
    A partial match is defined as one name being a prefix of the other when split by ';'.
    For example:
      name1 = "ATCOPIA34#LTR/Copia~LTRlen:246;NESTED_IN:ATLINE1_12#LINE/L1;CUT_BY:Os0007_Castaway#MITE/Tourist;CUT_BY:Os3447_LTR#LTR/Copia"
      name2 = "ATCOPIA34#LTR/Copia~LTRlen:246;NESTED_IN:ATLINE1_12#LINE/L1;CUT_BY:Os0007_Castaway#MITE/Tourist"
    Here, name2 is a prefix of name1.
    """
    if ';' not in name1 or ';' not in name2:
        return False  # Both must have supplemental info.
    parts1 = name1.split(';')
    parts2 = name2.split(';')
    # Check if the smaller list is a prefix of the larger one.
    if len(parts1) <= len(parts2) and parts1 == parts2[:len(parts1)]:
        return True
    elif len(parts2) < len(parts1) and parts2 == parts1[:len(parts2)]:
        return True
    return False

def count_intact_TE_count(features):
    """
    Count the number of intact TE features from the input BED (features list) based on the following rules:
      (1) If the feature ID (i.e. the part of NAME before the first semicolon) contains 'gene' (case-insensitive),
          then it is NOT an intact TE.
      (2) If the first supplemental info (i.e. the text after the first semicolon) contains 'CUT_BY',
          then it is NOT an intact TE.
      (3) If there is at least one supplemental info field (i.e. at least one semicolon in the NAME)
          and an identical or partial (prefix) NAME occurs within 100 lines in the BED (with an identical TSD and strand),
          then it is NOT an intact TE.
      (4) If the feature ID (NAME before the first semicolon) contains '_SOLO', then it is NOT an intact TE.
    """
    n = len(features)
    intact_flags = [True] * n

    # Rule 1: Exclude features whose feature ID (before first semicolon) contains "gene"
    for i, feat in enumerate(features):
        feature_id = feat['name'].split(';')[0]
        if 'gene' in feature_id.lower():
            intact_flags[i] = False

    # Rule 4: Exclude features whose feature ID (before first semicolon) contains "_SOLO"
    for i, feat in enumerate(features):
        feature_id = feat['name'].split(';')[0]
        if "_SOLO" in feature_id:
            intact_flags[i] = False

    # Rule 2: Exclude features where the first supplemental field contains "CUT_BY"
    for i, feat in enumerate(features):
        if intact_flags[i]:
            parts = feat['name'].split(';')
            if len(parts) > 1:
                if 'CUT_BY' in parts[1]:
                    intact_flags[i] = False

    # Rule 3: For features with one or more supplemental info fields,
    # if another feature (within 100 lines before or after) has an identical TSD,
    # identical strand, and a NAME that is an exact or partial (prefix) match, mark as not intact.
    for i, feat in enumerate(features):
        if not intact_flags[i]:
            continue
        parts = feat['name'].split(';')
        if len(parts) <= 1:
            continue  # No supplemental info; skip rule 3.
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
    """
    Compute the distribution (counts) of intact TEs per (te_class, te_superfamily)
    from the BED file. Only features that are intact based on the rules (1)-(4) are used.
    Returns a dict mapping (te_class, te_superfamily) -> count.
    """
    n = len(features)
    intact_flags = [True] * n

    # Rule 1
    for i, feat in enumerate(features):
        feature_id = feat['name'].split(';')[0]
        if 'gene' in feature_id.lower():
            intact_flags[i] = False

    # Rule 4
    for i, feat in enumerate(features):
        feature_id = feat['name'].split(';')[0]
        if "_SOLO" in feature_id:
            intact_flags[i] = False

    # Rule 2
    for i, feat in enumerate(features):
        if intact_flags[i]:
            parts = feat['name'].split(';')
            if len(parts) > 1:
                if 'CUT_BY' in parts[1]:
                    intact_flags[i] = False

    # Rule 3: Check for identical or partial (prefix) NAMEs within 100 lines
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

def main():
    args = parse_args()
    # Set random seed
    if args.seed is not None:
        random.seed(args.seed)

    print("Reading genome FASTA ...")
    genome_raw = read_fasta(args.genome)  # {chrom: seq_string}
    # Convert to lists for easy insertion
    print("Converting genome to editable lists ...")
    genome = convert_genome_to_dict_of_lists(genome_raw)

    print("Reading TE FASTA ...")
    te_raw = read_fasta(args.TE)  # {header: seq_string}

    # We'll parse each TE's header to gather class/superfam info on-the-fly during insertion
    # and keep them in te_raw as is.

    print("Reading BED file ...")
    features = parse_bed(args.bed)  # list of dicts
    # Make a copy of the original BED order for counting intact TEs.
    features_original = features[:]

    # Compute the number of intact TEs based on the provided rules.
    intact_TE_count = count_intact_TE_count(features_original)
    print(f"Total number of intact TEs from BED: {intact_TE_count}")

    # Also compute the distribution of intact TEs by te_class and te_superfamily.
    intact_distribution = get_intact_TE_distribution(features_original)
    print("Distribution of intact TEs by (te_class, te_superfamily):")
    for key, count in intact_distribution.items():
        print(f"  {key[0]}/{key[1]}: {count}")

    # Compute total insertions.
    if args.fix_in is not None:
        total_insertions = args.fix_in
        print(f"Using fixed number of TE insertions: {total_insertions}")
    else:
        total_insertions = int(args.rate * intact_TE_count * args.generations)
        print(f"Rate: {args.rate}")
        print(f"Generations: {args.generations}")
        print(f"Total TE insertions to perform: {total_insertions}")

    # Determine the number of insertions per TE category based on the intact distribution.
    # We assign counts proportional to the distribution.
    insertion_events = []  # each event will be a tuple (te_class, te_superfamily)
    intact_total = sum(intact_distribution.values())
    for cat, count in intact_distribution.items():
        # Calculate insertion count for this category
        insertion_count = round(total_insertions * (count / intact_total))
        for i in range(insertion_count):
            insertion_events.append(cat)

    # Adjust if rounding made the total differ from total_insertions.
    diff = total_insertions - len(insertion_events)
    if diff > 0:
        # Add extra events randomly weighted by the intact distribution.
        categories = list(intact_distribution.keys())
        weights = [intact_distribution[cat] for cat in categories]
        for i in range(diff):
            chosen_cat = random.choices(categories, weights=weights)[0]
            insertion_events.append(chosen_cat)
    elif diff < 0:
        # Remove some events at random.
        for i in range(-diff):
            removal_index = random.randrange(len(insertion_events))
            insertion_events.pop(removal_index)

    # --- New TE Birth functionality ---
    if args.birth_rate and args.birth_file and args.TE_ratio_file:
        # Parse birth_file to extract NUMBER_OF_INITIAL_TES.
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

        # Parse TE_ratio_file to get non-normalized weights.
        te_ratio = {}  # mapping (te_class, te_superfamily) -> weight
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
        # Normalize the weights.
        total_weight = sum(te_ratio.values())
        for k in te_ratio:
            te_ratio[k] /= total_weight
        # Create lists of categories and their normalized weights.
        born_categories = list(te_ratio.keys())
        born_weights = [te_ratio[cat] for cat in born_categories]
        # Append born TE insertion events weighted by TE_ratio.
        for i in range(number_of_born_tes):
            chosen_cat = random.choices(born_categories, weights=born_weights, k=1)[0]
            insertion_events.append(chosen_cat)
        print(f"Total insertion events after adding born TEs: {len(insertion_events)}")
    # --- End New TE Birth functionality ---

    # Shuffle the insertion events to mix insertions from different categories.
    random.shuffle(insertion_events)

    # Sort existing features by chrom and start.
    features.sort(key=lambda x: (x['chrom'], x['start']))

    # Initialize counters for nested and non-nested insertions.
    nested_count = 0
    non_nested_count = 0

    # Process each insertion event.
    for insertion_index, (te_class_target, te_superfamily_target) in enumerate(insertion_events):
        # 1) Pick a random chromosome
        chroms = list(genome.keys())
        chosen_chrom = random.choice(chroms)
        chrom_length = len(genome[chosen_chrom])
        if chrom_length == 0:
            continue  # skip empty chromosomes

        # 2) Pick a random insertion position in [0, chrom_length]
        insertion_pos = random.randint(0, chrom_length)

        # 3) Pick a random TE from the TE library that matches the desired category.
        te_header, te_seq = pick_random_TE_by_category(te_raw, te_class_target, te_superfamily_target)

        # 4) (For safety) Parse the TE class and superfamily from the header.
        te_class, te_superfamily = extract_te_info(te_header)

        # 5) Determine TSD length.
        tsd_length = get_tsd_length(te_class, te_superfamily)

        # 6) Determine orientation (50% chance +, 50% -).
        strand = "+" if random.random() < 0.5 else "-"
        if strand == "-":
            te_seq = reverse_complement(te_seq)
        te_seq_list = list(te_seq)

        # 7) Handle TSD extraction; if insertion_pos < tsd_length, reduce TSD.
        if tsd_length > insertion_pos:
            print(f"Warning: TSD length {tsd_length} is greater than insertion position {insertion_pos}. Reducing.")
            tsd_length = insertion_pos

        tsd_seq_list = []
        if tsd_length > 0:
            tsd_seq_list = genome[chosen_chrom][insertion_pos - tsd_length:insertion_pos]
        tsd_string = "".join(tsd_seq_list) if tsd_length > 0 else "NA"

        # 8) Insert TE + duplicated TSD into the genome.
        # Final structure: [genome_part][TE][TSD][genome_rest]
        genome[chosen_chrom][insertion_pos:insertion_pos] = te_seq_list
        if tsd_length > 0:
            genome[chosen_chrom][insertion_pos + len(te_seq_list):insertion_pos + len(te_seq_list)] = tsd_seq_list

        shift_amount = len(te_seq_list) + tsd_length

        # 9) Create a new feature for the inserted TE.
        te_start = insertion_pos
        te_end   = insertion_pos + len(te_seq_list)  # TSD is not part of the TE proper.
        new_te_name = te_header  # will be modified below if nesting occurs.
        new_te_feature = {
            'chrom': chosen_chrom,
            'start': te_start,
            'end':   te_end,
            'name':  new_te_name,
            'strand': strand,
            'tsd':   tsd_string if tsd_length > 0 else 'NA',
        }

        # 10) Adjust existing features on the same chromosome.
        updated_features = []
        nesting_happened = False
        for feat in features:
            if feat['chrom'] != chosen_chrom:
                updated_features.append(feat)
                continue

            if feat['end'] <= insertion_pos:
                # Feature entirely before insertion site.
                updated_features.append(feat)
            elif feat['start'] >= insertion_pos:
                # Feature entirely after insertion site; shift coordinates.
                feat['start'] += shift_amount
                feat['end']   += shift_amount
                updated_features.append(feat)
            else:
                # Insertion occurs within this feature => nesting.
                cut_feat_end = insertion_pos
                left_piece = {
                    'chrom': feat['chrom'],
                    'start': feat['start'],
                    'end':   cut_feat_end,
                    'name':  feat['name'] + f";CUT_BY:{new_te_name}",
                    'strand': feat['strand'],
                    'tsd':  feat['tsd'],
                }
                right_length = feat['end'] - cut_feat_end
                right_piece_start = insertion_pos + shift_amount  # Adjust for inserted bases.
                right_piece_end   = right_piece_start + right_length

                if right_length > 0:
                    right_piece = {
                        'chrom': feat['chrom'],
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

                # Mark the new TE as nested.
                new_te_feature['name'] += f";NESTED_IN:{feat['name']}"
                nesting_happened = True

        features = updated_features

        # Finally, add the new TE feature.
        features.append(new_te_feature)

        print(f"Insertion #{insertion_index+1}: Inserted TE '{te_header}' on {chosen_chrom}:{insertion_pos}, strand={strand}, TSD={tsd_string}")
        if nesting_happened:
            print(f"  => Nested inside an existing feature. TE name now: {new_te_feature['name']}")

        # Count the type of insertion.
        if nesting_happened:
            nested_count += 1
        else:
            non_nested_count += 1

    # Sort final features.
    features.sort(key=lambda x: (x['chrom'], x['start']))

    # Convert genome back to strings.
    final_genome = convert_genome_back_to_fasta(genome)

    # Write final BED.
    out_bed = f"{args.output}.bed"
    print(f"Writing updated BED to {out_bed}")
    write_bed(features, out_bed)

    # Write final FASTA.
    out_fasta = f"{args.output}.fasta"
    print(f"Writing updated FASTA to {out_fasta}")
    with open(out_fasta, 'w') as f:
        for chrom in sorted(final_genome.keys()):
            f.write(f">{chrom}\n")
            seq = final_genome[chrom]
            f.write(seq + "\n")

    # Print summary of nested and non-nested insertions.
    total_complete = nested_count + non_nested_count
    print(f"Total TE insertions performed: {total_complete} (Nested: {nested_count}, Non-nested: {non_nested_count})")

    print("Done.")

if __name__ == "__main__":
    main()
