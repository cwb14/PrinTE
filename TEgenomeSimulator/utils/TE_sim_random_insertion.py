import os
import sys
import argparse
import random
import yaml
import numpy
from Bio import SeqIO
from Bio import Seq
from pathlib import Path

class Repeat:
    """Store the information of the TE families.
    This class has been modified by THC"""
    def __init__(self, name, subclass, superfamily, sequence, num_rep, identity, sd, indels, tsd, frag, nest, extra_attributes=None):
        self.name = name
        self.subclass = subclass
        self.superfamily = superfamily
        self.sequence = sequence
        self.num_rep = num_rep
        self.identity = identity
        self.sd = sd
        self.indels = indels
        self.tsd = tsd
        self.frag = frag
        self.nest = nest
        self.extra_attributes = extra_attributes if extra_attributes else {}

# Load params_chr from YAML file config_genome.yml in same directory (modified by THC).
def parse_random_genome_yaml(config_file):
    params_chr = yaml.load(open(config_file, 'r'), Loader=yaml.FullLoader)
    return params_chr

# Load params_chr from YAML file config_custom_genome.yml in same directory (function created by THC).
def parse_custom_genome_yaml(config_file):
    params_chr = yaml.load(open(config_file, 'r'), Loader=yaml.FullLoader)
    return params_chr

# Load existing non-repeat genome sequences (function created by THC).
def load_custom_genome(params_chr):
    chr_id = list(params_chr['chrs'].keys())
    chrs_dict = {}
    fasta = SeqIO.index(params_chr['genome_fasta'], "fasta")
    for chrid in chr_id:
        chrs_dict[chrid] = fasta[chrid].seq
    return chrs_dict

# Load collection of repeats and params for chrs simulation (modified by THC to include extra attributes).
def load_repeats_chr(params_chr):
    chr_id = list(params_chr['chrs'].keys())
    repeats_dict = {}
    
    # Parse all fasta records and extract extra attributes
    fasta_dict = {}
    with open(params_chr['rep_fasta'], 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description  # Entire header line
            # Split header on ';' to separate base name and attributes
            parts = header.split(';')
            # Extract only the base name before any space
            base_name = parts[0].lstrip('>').split()[0]
            extra_attributes = {}
            for attr in parts[1:]:
                if ':' in attr:
                    key, value = attr.split(':', 1)
                    extra_attributes[key] = value
            fasta_dict[base_name] = (str(record.seq).upper(), extra_attributes)
    
    # Read repeats from rep_list and associate with fasta sequences and extra attributes
    with open(params_chr['rep_list'], 'r') as repeats_file:
        next(repeats_file)  # Skip header line
        for line in repeats_file:
            elem = line.rstrip().split()
            name = elem[0]
            if name not in fasta_dict:
                print(f"Warning: Repeat {name} not found in fasta.")
                continue
            sequence, extra_attributes = fasta_dict[name]
            subclass = elem[1]
            superfamily = elem[2]
            num_rep = int(elem[3])
            identity = int(elem[4])
            sd = int(elem[5])
            indels = int(elem[6])
            tsd = [int(x) for x in elem[7].rstrip().split(",")]
            frag = int(elem[9])
            nest = int(elem[10])
            repeat = Repeat(name, subclass, superfamily, sequence, num_rep, identity, sd, indels, tsd, frag, nest, extra_attributes)
            repeats_dict[name] = repeat
    return repeats_dict
    
### Load other variables ###
# Calculate length of sequence of all repeats.
# sum_rep_length = sum([len(rep.sequence) * rep.num_rep for rep in repeats]).
# Calculate length of sequence that is going to be randomly generated.
# rand_seq_length = seq_length - sum_rep_length.

def generate_random_sequence(params):
    # Create DNA alphabet for random sequence.
    alphabet = ["A", "T", "G", "C"]
    gc = float(params['gc_content']) / 100
    weights = [(1 - gc) / 2, (1 - gc) / 2, gc / 2, gc / 2]
    # Generate random sequence that is going to separate the repeats.
    # abase_sequence = "".join([random.choice(alphabet) for i in range(0,params['seq_length'])]).
    # base_sequence = random.choices(alphabet, weights, params['seq_length']).
    base_sequence = "".join(numpy.random.choice(alphabet, params['seq_length'], p=weights, replace=True))
    return base_sequence

# Generate random sequence for each chromosome (function created by THC).
def generate_random_chr_sequence(params_chr):
    # Create DNA alphabet for random chr sequence.
    chrs_dict = {}
    for chr_id in list(params_chr['chrs'].keys()):
        params = params_chr['chrs'][chr_id]
        chr_seq = generate_random_sequence(params)
        chrs_dict[chr_id] = chr_seq
    return chrs_dict

# Randomly select chrs and positions to insert all the repeats (function created by THC).
def assign_chr_coord_repeats(params_chr, repeats_dict):
    # Total number of non-nested repeats.
    total_num_rep = 0
    for rep in repeats_dict.values():
        num_rep = int(rep.num_rep - rep.num_rep * (rep.nest / 100.0))
        total_num_rep += num_rep
    # Total genome bp.
    total_genome_bp = sum([params_chr['chrs'][chrid]['seq_length'] for chrid in params_chr['chrs']])
    # Create a dictionary containing chr as key and the cumulated chr length as value.
    chr_len = list(params_chr['chrs'][chrid]['seq_length'] for chrid in params_chr['chrs'])
    chr_len_cumm = list(sum(chr_len[0:x:1]) for x in range(0, len(chr_len) + 1))[1:]
    chr_id = list(params_chr['chrs'][chrid]['prefix'] for chrid in params_chr['chrs'])
    chr_cumm_dict = {}
    for i in range(0, len(chr_id)):
        ch = chr_id[i]
        end_cumm = chr_len_cumm[i]
        if i == 0:
            start_cumm = 1
        else:
            start_cumm = chr_len_cumm[i - 1] + 1
        chr_cumm_dict[ch] = (start_cumm, end_cumm)
    # Create random positions on the stitched genome and sort.
    random_coords = random.sample(range(total_genome_bp - 1), total_num_rep)
    random_coords.sort()
    # Bin the insertion positions to corresponding chr.
    random_repeats_coords = []
    for num in random_coords:
        for key, values in chr_cumm_dict.items():
            start_cumm, end_cumm = values
            if num <= end_cumm:
                random_repeats_coords.append((key, num - start_cumm + 1))
                break
    return random_repeats_coords

# Shuffle the order of TE families to be inserted.
def shuffle_repeats(repeats_dict):
    allnames = []
    allpositions = []
    for rep in repeats_dict.values():
        num_rep = int(rep.num_rep - rep.num_rep * (rep.nest / 100.0))
        names = num_rep * [rep.name]
        n_frags = int(((num_rep * rep.frag) / 100))
        positions = [0] * num_rep
        if n_frags > 0:
            sample_changes = random.sample(range(len(positions)), min(n_frags, len(positions)))
            for f in sample_changes:
                positions[f] += 1
        allnames += names
        allpositions += positions
    name_pos = list(zip(allnames, allpositions))
    random.shuffle(name_pos)
    return name_pos

# Get identity using a normal distribution (modified by THC).
def get_identity(mean, sd):
    # identity = int(numpy.random.normal(mean, sd, 1)) # int an array is deprecated in NumPy 1.25.
    # to prevent the warning message, use the following instead.
    identity = numpy.random.normal(mean, sd, 1)
    identity = int(identity[0])
    while identity > 100 or identity < 0:
        # identity = int(numpy.random.normal(mean, sd, 1)).
        identity = numpy.random.normal(mean, sd, 1)
        identity = int(identity[0])
    return identity

# Generate vector of coords for base_changes and indels.
def generate_mismatches(sequence, identity, indels, indel_size_range):
    alphabet = ["T", "G", "C", "A"]
    seq_len = len(sequence)
    seq = sequence
    # Calculate number of nucleotides that need to be changed (i.e. SNPs).
    num_changes = seq_len - int(round((seq_len * identity / 100.0)))
    # Generate a vector containing locations for SNPs.
    if num_changes > 0:
        pos_changes_vec = random.sample(range(seq_len), num_changes)
    else:
        pos_changes_vec = []
    # Calculate number of SNPs that need to be changed as indels.
    num_indels = int(round(num_changes * (indels / 100.0)))
    num_indels = min(num_indels, len(pos_changes_vec))
    # Generate a vector containing SNP locations to be changed to indels.
    if num_indels > 0:
        indel_changes_vec = random.sample(pos_changes_vec, num_indels)
    else:
        indel_changes_vec = []
    # Calculate indel sizes for each indel position.
    indel_sizes = [random.randint(indel_size_range[0], indel_size_range[1]) for _ in indel_changes_vec]
    # Replace single positions with tuples of (position, size).
    indel_changes_with_size = list(zip(indel_changes_vec, indel_sizes))
    # Generate a vector containing SNP locations (removal of indel locations).
    base_changes_vec = list(set(pos_changes_vec) - set(indel_changes_vec))
    base_changes_vec.sort()
    indel_changes_with_size.sort()
    return base_changes_vec, indel_changes_with_size

# Add SNPs.
def add_base_changes(repeat_seq, base_changes_vec):
    alphabet = ["T", "G", "C", "A"]
    repeat_seq_list = list(repeat_seq)
    for pos in base_changes_vec:
        new_base = random.choice(list(set(alphabet) - set(repeat_seq_list[pos])))
        repeat_seq_list[pos] = new_base
    new_repeat_seq = "".join(repeat_seq_list)
    return new_repeat_seq

## Add indels.
def add_indels(repeat_seq, indels_changes_with_size):
    alphabet = ["T", "G", "C", "A"]
    repeat_seq_list = list(repeat_seq)
    # Sort indels by position to avoid affecting subsequent positions.
    indels_changes_with_size_sorted = sorted(indels_changes_with_size, key=lambda x: x[0])
    offset = 0  # To keep track of the sequence length changes.
    for pos, size in indels_changes_with_size_sorted:
        pos += offset
        # Randomly choose 1 for insertion and 0 for deletion.
        if random.choice([0, 1]):
            # Insertion.
            new_bases = ''.join(random.choices(alphabet, k=size))
            repeat_seq_list.insert(pos, new_bases)
            offset += size
        else:
            # Deletion.
            del_length = min(size, len(repeat_seq_list) - pos)
            if del_length > 0:
                del repeat_seq_list[pos:pos + del_length]
                offset -= del_length
    new_repeat_seq = "".join(repeat_seq_list)
    return new_repeat_seq

# Add TSD if required (modified by THC to allow customised TSD length range).
def create_TSD(tsd_min, tsd_max, identity, indels, indel_size_range):
    alphabet = ["T", "G", "C", "A"]
    tsd_length = random.randint(tsd_min, tsd_max)
    tsd_seq_5 = "".join([random.choice(alphabet) for _ in range(tsd_length)])
    tsd_len = len(tsd_seq_5)
    base_changes_vec, indels_changes_with_size = generate_mismatches(tsd_seq_5, identity, indels, indel_size_range)
    tsd_seq_mismatches = add_base_changes(tsd_seq_5, base_changes_vec)
    tsd_seq_3 = add_indels(tsd_seq_mismatches, indels_changes_with_size)
    return tsd_seq_5, tsd_seq_3

# Fragment TE sequence.
def fragment(seq, include_tsd, fragmentation_end, tsd_seq_5, tsd_seq_3):
    frag_size = 100
    len_seq = len(seq)
    # Decide on the proportion of TE sequence to be maintained.
    if len_seq < 500:
        frag_size = random.randint(70, 99)
    else:
        frag_size = random.randint(40, 99)
    # Calculate the length of the TE sequence to be removed.
    cut_length = int(len_seq * ((100 - frag_size) / 100.0))
    
    if fragmentation_end == '5':
        # Fragment at 5' end.
        fragmented_seq = seq[cut_length:]
        if not include_tsd:
            # Remove TSD from the 5' end if present.
            fragmented_seq = fragmented_seq[len(tsd_seq_5):]
    else:
        # Fragment at 3' end.
        fragmented_seq = seq[: -cut_length]
        if not include_tsd:
            # Remove TSD from the 3' end if present.
            fragmented_seq = fragmented_seq[:-len(tsd_seq_3)]
    return fragmented_seq, frag_size, cut_length

# Generate new sequence including the repeats in the random one (modified by THC).
def generate_sequence(repeats_dict, rand_rep_pos, rand_seq, shuffled_repeats, indel_size_range):
    seq = ""
    pre_n = 0
    n = 0
    new_repeats_coord = []
    for n, m in zip(rand_rep_pos, shuffled_repeats):
        # Get sequence of repeat.
        repeat_seq = repeats_dict[m[0]].sequence

        # Get family name, subclass and superfamily.
        family = repeats_dict[m[0]].name
        subclas = repeats_dict[m[0]].subclass
        superfam = repeats_dict[m[0]].superfamily        

        # Get identity from a normal distribution.
        identity = get_identity(repeats_dict[m[0]].identity, repeats_dict[m[0]].sd)

        # Get base_changes and indels vectors and identity.
        identity_fix = identity + (100 - identity) * 0.5
        base_changes_vec, indels_changes_with_size = generate_mismatches(
            repeats_dict[m[0]].sequence, identity_fix, repeats_dict[m[0]].indels, indel_size_range
        )

        # Add mismatches to original repeat sequence.
        repeat_seq_mismatches = add_base_changes(repeat_seq, base_changes_vec)

        # Add indels to original repeat sequence.
        new_repeat_seq = add_indels(repeat_seq_mismatches, indels_changes_with_size)

        # Check if TE creates TSDs.
        if repeats_dict[m[0]].tsd != [0, 0]:
            # Generate TSD.
            TSD_min = repeats_dict[m[0]].tsd[0]
            TSD_max = repeats_dict[m[0]].tsd[1]
            tsd_seq_5, tsd_seq_3 = create_TSD(TSD_min, TSD_max, identity_fix, repeats_dict[m[0]].indels, indel_size_range)
        else:
            tsd_seq_5, tsd_seq_3 = "", ""

        # Assign sequence to a random strand.
        frag = 100  # Initiate frag size at 100% of the TE seq that has underwent identity/indel/TSD check.
        cut = 0  # Initiate cut size at 0%.
        strands = ["+", "-"]
        strand = random.choice(strands)

        # Determine fragmentation.
        if m[1] == 1:
            # Decide fragmentation direction and TSD inclusion.
            fragmentation_end = random.choice(['5', '3'])  # 50% chance.
            include_tsd = random.choice([True, False])     # 50% chance.
            # Fragment the sequence.
            new_repeat_seq_frag, frag, cut = fragment(
                new_repeat_seq, include_tsd, fragmentation_end, tsd_seq_5, tsd_seq_3
            )
            # Update TSD sequences based on fragmentation.
            if not include_tsd:
                if fragmentation_end == '5':
                    tsd_seq_5 = ""
                else:
                    tsd_seq_3 = ""
        else:
            new_repeat_seq_frag = new_repeat_seq

        # Apply strand sense.
        if strand == "-":
            new_repeat_seq_str = str(Seq.Seq(new_repeat_seq_frag).reverse_complement())
            # Reverse TSD sequences if present.
            if tsd_seq_5:
                tsd_seq_5 = str(Seq.Seq(tsd_seq_5).reverse_complement())
            if tsd_seq_3:
                tsd_seq_3 = str(Seq.Seq(tsd_seq_3).reverse_complement())
        else:
            new_repeat_seq_str = new_repeat_seq_frag

        # Concatenate TSDs and TE sequence.
        new_repeat_seq_tsd = tsd_seq_5 + new_repeat_seq_str + tsd_seq_3
        seq += rand_seq[pre_n:n] + new_repeat_seq_tsd

        # Get new repeat sequence end coordinate.
        repeat_end = len(seq) - len(tsd_seq_3)
        repeat_start = repeat_end - len(new_repeat_seq_str) + 1

        # Append to vector new data about new repeat useful for a GFF.
        new_repeats_coord.append([
            str(repeat_start),
            str(repeat_end),
            new_repeat_seq_str,
            identity,
            frag,
            strand,
            family,
            subclas,
            superfam,
            tsd_seq_5 if tsd_seq_5 else "None",
            tsd_seq_3 if tsd_seq_3 else "None"
        ])
        # Sets new end coordinate as start for next round.
        pre_n = n
    # At the last step add the remaining base sequence.
    seq += rand_seq[n:]
    # Return sequences and repeat data.
    return seq, new_repeats_coord

## Generate new sequence including the repeats in the random chr sequences (function created by THC).
def generate_genome_sequence(repeats_dict, rand_rep_pos, rand_chr_dict, shuffled_repeats, indel_size_range):
    # Create placeholders for the genome seq and TE coordinates.
    genome_dict = {}
    new_repeats_coord_dict = {}
    # Iterate through chromosomes.
    shuffled_start_index = 0
    for chromosome in rand_chr_dict.keys():
        # Capture chromosome sequence.
        rand_seq = rand_chr_dict[chromosome]
        # Capture the random position corresponding to the chromosome.
        rand_rep_pos_filter = [coord for chr_id, coord in rand_rep_pos if chr_id == chromosome]
        # Define where to slice the shuffled_repeats list.
        shuffled_end_index = shuffled_start_index + len(rand_rep_pos_filter)
        # Extract the slice of the shuffled_repeats for the corresponding chromosome.
        shuffled_repeats_sliced = shuffled_repeats[shuffled_start_index:shuffled_end_index]
        # The end of the slice becomes the start for next chromosome.
        shuffled_start_index = shuffled_end_index
        # Call generate_sequence function chromosome by chromosome.
        sequence, new_repeats_coord = generate_sequence(
            repeats_dict, rand_rep_pos_filter, rand_seq, shuffled_repeats_sliced, indel_size_range
        )
        genome_dict[chromosome] = sequence
        new_repeats_coord_dict[chromosome] = new_repeats_coord 
    # Return sequences and repeat data.
    return genome_dict, new_repeats_coord_dict

# Print final genome sequence to file (function created by THC).
def print_genome_data(genome_dict, new_repeats_coord_dict, params, out_dir, repeats_dict):
    # Setup output directory.
    file_prefix = str(params['prefix'])
    final_out = os.path.join(out_dir, f'TEgenomeSimulator_{file_prefix}_result')
    os.makedirs(final_out, exist_ok=True)
    
    # Output files.
    genome_fa = f"{file_prefix}_genome_sequence_out.fasta"
    te_fa = f"{file_prefix}_repeat_sequence_out.fasta"
    te_gff = f"{file_prefix}_repeat_annotation_out.gff"
    
    # Create genome fasta file.
    with open(os.path.join(final_out, genome_fa), "w") as fasta_out:
        for chromosome in genome_dict.keys():
            seq = str(genome_dict[chromosome])
            fasta_out.write(f">{chromosome}\n{seq}\n")
    
    # Collapse new_repeats_coord_dict to a list with chr information.
    new_repeats_coord_list = []
    for chromosome, replists in new_repeats_coord_dict.items():
        for replist in replists:
            new_repeats_coord_list.append([chromosome] + replist)   
    
    # Create repeat fasta and gff files.
    with open(os.path.join(final_out, te_fa), "w") as fasta_rep, open(os.path.join(final_out, te_gff), "w") as gff_rep:
        counts = 1
        for n in new_repeats_coord_list:
            chr_id = str(n[0])
            start = str(n[1])
            end = str(n[2])
            repeat_seq = str(n[3])
            identity = str(n[4] / 100)
            frag = str(n[5] / 100)
            strand = str(n[6])
            family = str(n[7])
            subclas = str(n[8])
            superfam = str(n[9])
            tsd_5 = str(n[10])
            tsd_3 = str(n[11])
            te_id = str(counts).zfill(7)  # prints at least 7 characters wide; i.e. at most 9,999,999 TE insertions.
            repeat_name = f">{family}_TE{te_id}#{superfam} [Location={chr_id}:{start}-{end};Identity={identity};Integrity={frag};TSD_5={tsd_5};TSD_3={tsd_3}]\n"
            repeat_sequence = f"{repeat_seq}\n"

            fasta_rep.write(repeat_name)
            fasta_rep.write(repeat_sequence)

            # Construct the attribute field with TSD information and extra attributes.
            attributes = (
                f"ID={family}_TE{te_id};"
                f"Name=TE{te_id};"
                f"Classification={superfam};"
                f"Identity={identity};"
                f"Integrity={frag};"
                f"TSD_5={tsd_5};"
                f"TSD_3={tsd_3}"
            )

            # Append extra attributes from the Repeat object, if any
            extra_attrs = repeats_dict.get(family, {}).extra_attributes
            if extra_attrs:
                extra_attrs_str = ';'.join([f"{k}:{v}" for k, v in repeats_dict[family].extra_attributes.items()])
                attributes += f";{extra_attrs_str}"

            # Write to GFF
            gff_rep.write("\t".join([
                chr_id,
                "TEgenomeSimulator",
                subclas,
                start,
                end,
                ".",
                strand,
                ".",
                attributes + "\n"
            ]))
            
            counts += 1

def main():
    # Set up argument parser.
    parser = argparse.ArgumentParser(description="Simulate random TE insertions into a genome.")
    
    # Define arguments.
    parser.add_argument('-M', '--mode', type=int, help="Mode for genome simulation (either 0 or 1).", required=True)
    parser.add_argument('-p', '--prefix', type=str, help="Prefix for output files.", required=True)
    parser.add_argument('-o', '--outdir', type=str, help="Output directory.", required=True)
    parser.add_argument('-i', '--indel_size', type=str, help="Indel size range, e.g., '1,5' for 1-5 bp. Default is '1,5'.", default="1,5")  # New argument
    
    # Parse arguments.
    args = parser.parse_args()
    mode = args.mode
    prefix = args.prefix
    out_dir = args.outdir
    indel_size_input = args.indel_size  # New argument.
    
    # Parse indel size range.
    try:
        indel_min, indel_max = map(int, indel_size_input.split(','))
        if indel_min < 1 or indel_max < indel_min:
            raise ValueError
        indel_size_range = (indel_min, indel_max)
    except:
        print("Error: --indel_size must be two integers separated by a comma, e.g., '1,5'.")
        sys.exit(1)
    
    print("\n")
    print("##############################################################")
    print("### Mutate TE sequence and perform non-overlap TE insertion###")
    print("##############################################################")
    print(f"Using mode {mode} (0 for random genome; 1 for custom genome)")
    print(f"Using indel size range: {indel_size_range[0]}-{indel_size_range[1]} bp")  # Inform about indel size range.
    
    # Config file.
    final_out = os.path.join(out_dir, f'TEgenomeSimulator_{prefix}_result')
    yml_file = f"TEgenomeSimulator_{prefix}.yml"
    print(f"Using config file {yml_file}")

    # Mode-dependent config file loading.
    if args.mode == 0:
        # Load chr parameters from yml file using parse_random_genome_yaml().        
        params_chr = parse_random_genome_yaml(os.path.join(final_out, yml_file))
    elif args.mode == 1:
        # Load chr parameters from yml file using parse_custom_genome_yaml().
        params_chr = parse_custom_genome_yaml(os.path.join(final_out, yml_file))
    else:
        print("Error: Mode must be either 0 (random genome) or 1 (custom genome).")
        sys.exit(1)

    # Set seed.
    seed = params_chr.get('seed', None)
    if seed:
        random.seed(seed)
        numpy.random.seed(seed)

    # Prepare genome for random or custom mode.
    if args.mode == 0:
        chrs_dict = generate_random_chr_sequence(params_chr)
    elif args.mode == 1:
        chrs_dict = load_custom_genome(params_chr)

    # Load repeat sequences.
    repeats_dict = load_repeats_chr(params_chr)
    # Assign TE coordinates randomly.
    repeats_coord = assign_chr_coord_repeats(params_chr, repeats_dict)
    # Shuffle the order the repeats to be inserted into the genome.
    shuffled_repeats = shuffle_repeats(repeats_dict)
    # Generate genome sequence with TE insertion.
    genome, new_repeats_coord = generate_genome_sequence(
        repeats_dict, repeats_coord, chrs_dict, shuffled_repeats, indel_size_range
    )
    # Output to fasta and gff files.
    print_genome_data(genome, new_repeats_coord, params_chr, out_dir, repeats_dict)

if __name__ == "__main__":
    main()

# END
