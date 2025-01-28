import os
import sys
import argparse
import random
import yaml
import re
import numpy
from Bio import SeqIO, Seq
from pathlib import Path
from TE_sim_random_insertion import (
    parse_random_genome_yaml,
    parse_custom_genome_yaml,
    load_repeats_chr,
    generate_mismatches,
    add_indels,
    add_base_changes,
    get_identity,
    create_TSD,
    fragment
)

# Turn GFF table into a big list where each item represents a row in GFF table.
def load_gff(gff_file):
    gff = []
    with open(gff_file) as gff_fh:
        for i in gff_fh:
            e = i.strip().split("\t")
            gff.append(e)
    return gff

# Modify the coordinates of TE loci locating downstream of the nested insertion for genome with multiple chromosomes (function created by THC).
def modify_genome_coords(offset, index, new_gff, chr_id):
    new_gff_aux = []
    for i in range(index, len(new_gff)):
        if new_gff[i][0] == chr_id:
            start = int(new_gff[i][3]) + offset
            new_gff[i][3] = str(start)
            end = int(new_gff[i][4]) + offset
            new_gff[i][4] = str(end)
    return new_gff

# Excluding Alu and SINE for nested insertion.
def filter_nonest(gff):
    c = 0
    vec_cand = []
    for i in gff:
        if "Alu" in i[8] or "SINE" in i[8]:
            pass
        else:
            vec_cand.append(c)
        c += 1
    return vec_cand

# Turn the inserted repeat FASTA file into a dictionary with TE_id as key (function created by THC).
def load_isrt_te_fa(inserted_te_fasta_file):
    isrt_te_dict = SeqIO.to_dict(SeqIO.parse(inserted_te_fasta_file, "fasta"))
    keys_copy = list(isrt_te_dict.keys())
    for key in keys_copy:
        new_key = re.sub(".*_TE", "TE", key)
        new_key = re.sub("#.*", "", new_key)
        isrt_te_dict[new_key] = isrt_te_dict.pop(key)
    return isrt_te_dict

# Generate genome with nested insertions (function modified by THC).
def generate_genome_nests(
    repeats,
    isrt_te_dict,
    gff,
    genome,
    indel_size_range  # New parameter.
):
    rep_count = []
    new_seq = ""
    nest_te_dict = {}
    
    # Helper function to format extra attributes
    def format_extra_attributes(extra_attributes):
        return ''.join([f";{key}:{value}" for key, value in extra_attributes.items()])
    
    # Create a list containing the TE family and the required copies for nested insertion.
    for i in repeats:
        rep_count += [i] * int(repeats[i].num_rep * (repeats[i].nest / 100.0))
    
    # Acquire the index of TEs from GFF file (excluding Alu and SINE).
    vec_cand = filter_nonest(gff)
    
    # Randomly pick the index of TE loci to be inserted with nested TEs.
    insert_index = random.sample(vec_cand, len(rep_count))
    
    new_gff = gff.copy()  # To avoid modifying the original GFF.
    sorted_table = sorted(zip(insert_index, rep_count))
    counter = 0
    n = 1
    
    for j, k in zip(sorted_table, rep_count):
        gff_sel = new_gff[j[0] + counter]
        chr_id = gff_sel[0]
        start = int(gff_sel[3])
        end = int(gff_sel[4])
        length = end - start + 1
        
        # Decide on the position of nested insertion in association to the length of host TE locus.
        pct_pos = random.randint(40, 60)
        ins_pos = int(round((pct_pos / 100.0) * length))
                        
        # Get family name, subclass and superfamily of nested TE.
        nest_name = repeats[k].name
        nest_subclas = repeats[k].subclass
        nest_superfam = repeats[k].superfamily 
        
        # Create SNPs, indels and TSD for the nested TE.
        nest_seq = repeats[k].sequence
        nest_identity = get_identity(repeats[k].identity, repeats[k].sd)
        nest_identity_fix = nest_identity + (100 - nest_identity) * 0.5
        nest_indels = repeats[k].indels
        base_changes_vec, indels_changes_with_size = generate_mismatches(
            nest_seq,
            nest_identity_fix,
            nest_indels,
            indel_size_range  # Pass the indel size range.
        )
        nest_seq_mismatches = add_base_changes(nest_seq, base_changes_vec)
        new_nest_seq = add_indels(nest_seq_mismatches, indels_changes_with_size)
        
        new_nest_seq_tsd = new_nest_seq
        tsd_5_len = tsd_3_len = 0
        if repeats[k].tsd != [0, 0]:
            TSD_min = repeats[k].tsd[0]
            TSD_max = repeats[k].tsd[1]
            tsd_seq_5, tsd_seq_3 = create_TSD(
                TSD_min,
                TSD_max,
                nest_identity_fix,
                nest_indels,
                indel_size_range  # Pass the indel size range.
            )
            new_nest_seq_tsd = tsd_seq_5 + new_nest_seq + tsd_seq_3
            tsd_5_len = len(tsd_seq_5)
            tsd_3_len = len(tsd_seq_3)
        
        # Fragmentation Logic (Enhanced).
        isFrag = random.choice([True, False])  # Two-thirds chance was 2 out of 3; Changed to 50% for consistency.
        new_nest_seq_tsd_frag = new_nest_seq_tsd
        if isFrag:
            # Decide fragmentation direction and TSD inclusion.
            fragmentation_end = random.choice(['5', '3'])  # 50% chance for 5' or 3'.
            include_tsd = random.choice([True, False])     # 50% chance to include TSD.
            new_nest_seq_tsd_frag, frag, cut = fragment(
                new_nest_seq_tsd,
                include_tsd,
                fragmentation_end,
                tsd_seq_5,
                tsd_seq_3
            )
            # Update TSD lengths based on fragmentation.
            if not include_tsd:
                if fragmentation_end == '5':
                    tsd_5_len = 0
                else:
                    tsd_3_len = 0
        else:
            frag = 100
        
        nest_len = len(new_nest_seq_tsd_frag)
        
        # Calculate the coordinates of the host TEs and nested TEs after insertion.
        new_end_1 = start + ins_pos
        new_start_2 = new_end_1 + nest_len 
        new_end_2 = new_start_2 + (length - ins_pos)
    
        # Apply strand sense.
        strands = ["+", "-"]
        strand = random.choice(strands)
        new_nest_seq_str = new_nest_seq_tsd_frag
        
        if strand == "-":
            new_nest_seq_str = str(Seq.Seq(new_nest_seq_tsd_frag).reverse_complement())
            # Reverse TSD sequences if present.
            if 'tsd_seq_5' in locals() and tsd_seq_5:
                tsd_seq_5 = str(Seq.Seq(tsd_seq_5).reverse_complement())
            if 'tsd_seq_3' in locals() and tsd_seq_3:
                tsd_seq_3 = str(Seq.Seq(tsd_seq_3).reverse_complement())
    
        # Prepare updated content to be put into GFF list.    
        nested_te_id = str(n).zfill(6)  # prints at least 6 characters wide; i.e. at most 999,999 nested TE insertions.
        nested_te_id = "TEn" + nested_te_id
        frag_note = ";Integrity=" + str((frag / 100))
        
        ori_seq_te_id = re.sub(".*_TE", "TE", gff_sel[8])
        ori_seq_te_id = re.sub(";Name.*", "", ori_seq_te_id)
        
        # Retrieve extra attributes from repeats_dict
        extra_attributes = repeats[k].extra_attributes if hasattr(repeats[k], 'extra_attributes') else {}
        formatted_extra_attrs = format_extra_attributes(extra_attributes)
        
        # Construct ori_line_1 without extra attributes
        ori_name_1 = [
            gff_sel[8]
            .replace(";Name", "_1;Name")
            .replace(";Clas", "_1;Clas")
            + ";Cut_at=" + str((pct_pos / 100)) + ";Cut_by=" + nest_name + "_" + nested_te_id
        ]
        ori_line_1 = gff_sel[:3] + [str(start)] + [str(new_end_1)] + gff_sel[5:8] + ori_name_1
        
        # Construct nest_name_in without extra attributes (will add extra attributes separately)
        nest_name_in = (
            f"ID={nest_name}_{nested_te_id};"
            f"Name={nested_te_id};"
            f"Classification={nest_superfam};"
            f"Identity={str((nest_identity / 100))}"
            f"{frag_note};"
            f"Nest_in={ori_seq_te_id}"
        )
        
        # Nested insertion that has undergone fragmentation does not have 5' tsd.
        if frag == 100:
            nested_line = (
                gff_sel[:2]
                + [nest_subclas]
                + [str(new_end_1 + 1 + tsd_5_len)]
                + [str(new_start_2 - tsd_3_len)]
                + [".", strand, "."]
                + [nest_name_in + formatted_extra_attrs]  # Append extra attributes only to nested insertion
            )
        else:
            nested_line = (
                gff_sel[:2]
                + [nest_subclas]
                + [str(new_end_1 + 1)]
                + [str(new_start_2 - tsd_3_len)]
                + [".", strand, "."]
                + [nest_name_in + formatted_extra_attrs]  # Append extra attributes only to nested insertion
            )
     
        # Construct ori_line_2 without extra attributes
        ori_name_2 = [
            gff_sel[8]
            .replace(";Name", "_2;Name")
            .replace(";Clas", "_2;Clas")
            + ";Cut_at=" + str((pct_pos / 100)) + ";Cut_by=" + nest_name + "_" + nested_te_id
        ]
        ori_line_2 = gff_sel[:3] + [str(new_start_2)] + [str(new_end_2 - 1)] + gff_sel[5:8] + ori_name_2
    
        n += 1
    
        # Update the entire GFF list.
        index = j[0] + counter
        new_gff.pop(index)
        new_gff_aux = new_gff[:index]
        new_gff_aux.append(ori_line_1)
        new_gff_aux.append(nested_line)
        new_gff_aux.append(ori_line_2)
        new_gff_aux += new_gff[index:]
        new_gff = new_gff_aux
        new_gff = modify_genome_coords(nest_len, index + 3, new_gff, chr_id)
        counter += 2
        
        # Update the all inserted repeat seq dictionary.      
        ori_seq_name_old = isrt_te_dict[ori_seq_te_id].description
        ori_seq_name_old_head = re.sub("-.*", "-", ori_seq_name_old)
        ori_seq_name_old_tail = re.sub(".*;I", ";I", ori_seq_name_old)
        ori_seq_name_old_tail = re.sub("]", "", ori_seq_name_old_tail)
        nest_note = ";Cut_at=" + str((pct_pos / 100)) + ";Cut_by=" + nest_name + "_" + nested_te_id + "]"
        ori_seq_name_new = ori_seq_name_old_head + str(new_end_2 - 1) + ori_seq_name_old_tail + nest_note
    
        nested_seq_name = (
            f"{nest_name}_{nested_te_id}#{nest_superfam} "
            f"[Location={chr_id}:{str(new_end_1 + 1 + tsd_5_len)}-{str(new_start_2 - tsd_3_len)};"
            f"Identity={str(nest_identity / 100)}"
            f"{frag_note};Nest_in={ori_seq_te_id}"
            f"{formatted_extra_attrs}"  # Append extra attributes only to nested insertion
        )
        
        isrt_te_dict[ori_seq_te_id].description = ori_seq_name_new
        
        # Create and update dictionary for nested TE seq.
        nest_te_dict[nested_te_id] = {
            'id': nested_seq_name,
            'seq': new_nest_seq_str[tsd_5_len:(nest_len - tsd_3_len)]
        }
                
        # Update the genome sequence.                
        seq = str(genome[chr_id].seq)
        new_seq = seq[:new_end_1] + new_nest_seq_str + seq[new_end_1:] 
        genome[chr_id].seq = new_seq
    
    return genome, isrt_te_dict, nest_te_dict, new_gff

# Print final sequence to files (function created by THC).
def print_genome_nest_data(
    genome,
    isrt_te_dict,
    nest_te_dict,
    new_gff,
    params,
    out_dir
):
    # Setup output directory.   
    file_prefix = str(params['prefix'])
    final_out = os.path.join(out_dir, f"TEgenomeSimulator_{file_prefix}_result")
    os.makedirs(final_out, exist_ok=True)  # Ensure the directory exists.
    
    # Specify output files.
    genome_fa = f"{file_prefix}_genome_sequence_out_final.fasta"
    te_fa = f"{file_prefix}_repeat_sequence_out_final.fasta"
    te_gff = f"{file_prefix}_repeat_annotation_out_final.gff"
    
    # For genome FASTA file.
    with open(os.path.join(final_out, genome_fa), "w") as fasta_out:
        for chromosome in genome:
            seq = str(genome[chromosome].seq)
            fasta_out.write(f">{chromosome}\n{seq}\n")
    
    # For all inserted TE sequences (including nested TEs).
    with open(os.path.join(final_out, te_fa), "w") as te_fa_out:
        for te in isrt_te_dict:
            header = str(isrt_te_dict[te].description)
            seq = str(isrt_te_dict[te].seq)
            te_fa_out.write(f">{header}\n{seq}\n")
        for nested_te in nest_te_dict:
            header = nest_te_dict[nested_te]['id']
            seq = nest_te_dict[nested_te]['seq']
            te_fa_out.write(f">{str(header)}\n{str(seq)}\n")
    
    # For the new GFF file.
    with open(os.path.join(final_out, te_gff), "w") as gff_out:
        for i in new_gff:
            gff_out.write("\t".join(map(str, i)) + "\n")

def main():
    # Set up argument parser.
    parser = argparse.ArgumentParser(description="Simulate nested TE insertions into a genome.")
    
    # Define arguments.
    parser.add_argument(
        '-M',
        '--mode',
        type=int,
        help="Mode for genome simulation (either 0 or 1).",
        required=True
    )
    parser.add_argument(
        '-p',
        '--prefix',
        type=str,
        help="Prefix for output files.",
        required=True
    )
    parser.add_argument(
        '-o',
        '--outdir',
        type=str,
        help="Output directory.",
        required=True
    )
    parser.add_argument(
        '-i',
        '--indel_size',
        type=str,
        help="Indel size range, e.g., '1,5' for 1-5 bp. Default is '1,5'.",
        default="1,5"  # New argument.
    )
    
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
    print("### Mutate TE sequence and perform nested TE insertion ###")
    print("##############################################################")
    print(f"Using mode {mode} (0 for random genome; 1 for custom genome)")
    print(f"Using indel size range: {indel_size_range[0]}-{indel_size_range[1]} bp")  # Inform about indel size range.
    
    # Config file.
    final_out = os.path.join(out_dir, f"TEgenomeSimulator_{prefix}_result")
    yml_file = f"TEgenomeSimulator_{prefix}.yml"
    print(f"Using config file {yml_file}")
    
    # Mode-dependent config file loading.
    if args.mode == 0:
        # Load chr parameters from YML file using parse_random_genome_yaml().        
        params_chr = parse_random_genome_yaml(os.path.join(final_out, yml_file))
    elif args.mode == 1:
        # Load chr parameters from YML file using parse_custom_genome_yaml().
        params_chr = parse_custom_genome_yaml(os.path.join(final_out, yml_file))
    else:
        print("Error: Mode must be either 0 (random genome) or 1 (custom genome).")
        sys.exit(1)
    
    # Set seed.
    seed = params_chr.get('seed', None)
    if seed:
        random.seed(seed)
        numpy.random.seed(seed)
    
    # Specify files created from previous step that generates non-overlapping insertions.
    file_prefix = str(params_chr['prefix'])
    gff_file = os.path.join(final_out, f"{file_prefix}_repeat_annotation_out.gff")
    fasta_file = os.path.join(final_out, f"{file_prefix}_genome_sequence_out.fasta")
    isrt_te_fasta = os.path.join(final_out, f"{file_prefix}_repeat_sequence_out.fasta")
    
    # Load FASTA files into dictionaries.
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    isrt_te = load_isrt_te_fa(isrt_te_fasta)
    
    # Load non-redundant TE library and the GFF file containing all TE insertions.
    repeats_dict = load_repeats_chr(params_chr)
    gff = load_gff(gff_file)
    
    # Generate nested insertion.
    genome, isrt_te_dict, nest_te_dict, new_gff = generate_genome_nests(
        repeats_dict,
        isrt_te,
        gff,
        genome,
        indel_size_range  # Pass the indel size range.
    )
    
    # Output new genome FASTA, all inserted TE FASTA, and GFF after nested insertion.
    print_genome_nest_data(genome, isrt_te_dict, nest_te_dict, new_gff, params_chr, out_dir)

if __name__ == "__main__":
    main()
# END
