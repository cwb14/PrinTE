import os
import sys
import argparse
import random
import re
from Bio import SeqIO
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""
        Prepare TE library table with simulation settings for TEgenomeSimulator.

        This script generates a table with simulated parameters for transposable element (TE) families based on a provided repeat FASTA file. 
        The simulation parameters such as copy number, sequence identity, standard deviation, indel proportion, fragmented loci proportion, 
        and nested insertion proportion can be customized via command-line arguments. If a TE superfamily is not recognized, it is 
        assigned to 'Unknown' with default settings, and a warning is issued.

        Example usage:
            python prep_sim_TE_lib.py -p Project1 -r repeats.fasta -m 10 -n 2 -o output_dir
        """
    )

    # Required arguments.
    parser.add_argument("-p", "--prefix", type=str, required=True,
                        help="Project prefix for the simulation.")
    parser.add_argument("-r", "--repeat", type=str, required=True,
                        help="Path to the repeat FASTA file.")
    parser.add_argument("-m", "--maxcp", type=int, required=True,
                        help="Maximum copy number to be simulated for a TE family.")
    parser.add_argument("-n", "--mincp", type=int, required=True,
                        help="Minimum copy number to be simulated for a TE family.")
    parser.add_argument("-o", "--outdir", type=str, required=True,
                        help="Output directory.")

    # Optional arguments for ranges.
    parser.add_argument("--idn_range", type=int, nargs=2, metavar=('MIN_IDN', 'MAX_IDN'), default=[80, 95],
                        help="Averaged sequence identity range for TE families (default: 80 95).")
    parser.add_argument("--sd_range", type=int, nargs=2, metavar=('MIN_SD', 'MAX_SD'), default=[1, 20],
                        help="Standard deviation range of averaged sequence identity (default: 1 20).")
    parser.add_argument("--indel_range", type=int, nargs=2, metavar=('MIN_INDEL', 'MAX_INDEL'), default=[5, 20],
                        help="Proportion range of INDEL to total SNP for each TE family (default: 5 20).")
    parser.add_argument("--frag_range", type=int, nargs=2, metavar=('MIN_FRAG', 'MAX_FRAG'), default=[50, 98],
                        help="Proportion range of fragmented TE loci for each TE family (default: 50 98).")
    parser.add_argument("--nest_range", type=int, nargs=2, metavar=('MIN_NEST', 'MAX_NEST'), default=[0, 30],
                        help="Proportion range of nested TE insertions for Copia or Gypsy families (default: 0 30).")

    # Random seed.
    parser.add_argument('-s', '--seed', type=int, 
                        default=1, help="Random seed (default: 1).")

    return parser.parse_args()

def main():
    args = parse_arguments()
    prefix = args.prefix
    te_fa = args.repeat
    max_cp = args.maxcp
    min_cp = args.mincp
    seed = args.seed
    out_dir = args.outdir

    idn_min, idn_max = args.idn_range
    sd_min, sd_max = args.sd_range
    indel_min, indel_max = args.indel_range
    frag_min, frag_max = args.frag_range
    nest_min, nest_max = args.nest_range

    print("\n")
    print("#########################################################")
    print("### Prepare TE library table with simulation settings ###")
    print("#########################################################")
    print(f"Using repeat FASTA file: {args.repeat}")
    print(f"Output directory set as: {args.outdir}")
    print(f"Project prefix: {args.prefix}")
    print(f"Random seed: {args.seed}")
    print("\n")

    # Set seed.
    if seed is not None:
        random.seed(seed)

    # Load TE library.
    try:
        te_lib = SeqIO.to_dict(SeqIO.parse(te_fa, "fasta"))
    except FileNotFoundError:
        print(f"Error: Repeat FASTA file '{te_fa}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading repeat FASTA file: {e}")
        sys.exit(1)

    # TE families.
    # Map clean TE family names (without metadata) to original headers
    te_family_map = {te_id.split(';')[0]: te_id for te_id in te_lib.keys()}
    te_family = list(te_family_map.keys())  # Cleaned list of TE families

    # TE superfamilies.
    te_superfamily = []
    for te_id in te_family:
        # Modified line to ignore metadata after ';'
        superfamily = re.sub(".*#", "", te_id).split(';')[0]
        te_superfamily.append(superfamily)

    # Define TE subclasses.
    subclass_dict = {
        'LTR_retrotransposon': ['LTR/Copia', 'LTR/Gypsy', 'LTR/Ty3', 'LTR/Solo', 'LTR/unknown'],
        'LINE_retrotransposon': ['LINE/unknown', 'LINE/L1'],
        'SINE_retrotransposon': ['SINE/unknown', 'SINE/tRNA'],
        'TIR_transposon': ['DNA/hAT', 'DNAnona/hAT', 'DNAauto/hAT', 
                           'DNA/CACTA', 'DNAnona/CACTA', 'DNAauto/CACTA', 
                           'DNA/Harbinger', 'DNA/MuDR', 'DNAnona/MULE', 
                           'DNAauto/MULE', 'DNA/Mariner'],
        'Helitron': ['DNA/Helitron', 'DNAnona/Helitron', 'DNAauto/Helitron'],
        'MITE': ['MITE/Stow', 'MITE/Tourist']
    }

    te_subclass = []
    unknown_sup_fams = set()

    for te_id, sup_fam in zip(te_family, te_superfamily):
        subclass_assigned = False
        for subclass, sup_fams in subclass_dict.items():
            if sup_fam in sup_fams:
                te_subclass.append(subclass)
                subclass_assigned = True
                break
        if not subclass_assigned:
            te_subclass.append('Unknown')
            unknown_sup_fams.add(te_id)
            print(f"Warning: TE superfamily '{sup_fam}' in TE family '{te_id}' is unrecognized. Assigned to 'Unknown' with TSD range 5-20.")

    # Copy number.
    print("## Randomly choosing copy number for each TE family ##")
    copy_number = [random.randint(min_cp, max_cp) for _ in te_family]
    print(f"Copy number range: {min_cp} - {max_cp}")

    # Identity.
    print("\n## Randomly choosing the averaged sequence identity for each TE family ##")
    identity = [random.randint(idn_min, idn_max) for _ in te_family]
    print(f"Averaged sequence identity range: {idn_min} - {idn_max}")

    # Standard deviation of identity.
    print("\n## Randomly choosing the standard deviation of averaged sequence identity for each TE family ##")
    sd = [random.randint(sd_min, sd_max) for _ in te_family]
    print(f"Standard deviation range: {sd_min} - {sd_max}")

    # Indel proportion.
    print("\n## Randomly choosing the proportion of INDEL to total SNP for each TE family ##")
    indel = [random.randint(indel_min, indel_max) for _ in te_family]
    print(f"INDEL proportion range: {indel_min} - {indel_max}")

    # TSD.
    print("\n## Setting the length of TSD based on TE superfamily ##")
    tsd = []
    for te_id, sup_fam in zip(te_family, te_superfamily):
        if sup_fam in subclass_dict['LTR_retrotransposon']:
            tsd_range = '5,5'
        elif sup_fam in subclass_dict['LINE_retrotransposon']:
            tsd_range = '5,20'
        elif sup_fam in subclass_dict['SINE_retrotransposon']:
            tsd_range = '5,20'
        elif sup_fam in subclass_dict['TIR_transposon']:
            if sup_fam in ['DNA/hAT', 'DNAnona/hAT', 'DNAauto/hAT']:
                tsd_range = '5,8'
            elif sup_fam in ['DNA/CACTA', 'DNAnona/CACTA', 'DNAauto/CACTA']:
                tsd_range = '2,4'
            elif sup_fam in ['DNA/Harbinger']:
                tsd_range = '3,3'
            elif sup_fam in ['DNA/MuDR', 'DNAnona/MULE', 'DNAauto/MULE']:
                tsd_range = '8,9'
            elif sup_fam in ['DNA/Mariner']:
                tsd_range = '2,2'
            else:
                tsd_range = '5,20'  # Default for any unforeseen cases within TIR_transposon.
        elif sup_fam in subclass_dict['Helitron']:
            tsd_range = '0,0'
        elif sup_fam in subclass_dict['MITE']:
            tsd_range = '2,10'
        else:
            tsd_range = '5,20'  # For Unknown.
        tsd.append(tsd_range)
    print("TSD length ranges have been set based on TE superfamilies.")

    # Length of each TE family.
    print("\n## Extracting the length of each TE family ##")
    length = [len(te_lib[te_family_map[te_id]].seq) for te_id in te_family]

    # Fragmented TE loci proportion.
    print("\n## Setting the proportion of fragmented TE loci for each TE family ##")
    fragment = [random.randint(frag_min, frag_max) for _ in te_family]
    print(f"Fragmented TE loci proportion range: {frag_min} - {frag_max}")

    # Nested TE insertions proportion.
    print("\n## Setting the proportion of nested TE insertions for each Copia or Gypsy family ##")
    nested = [random.randint(nest_min, nest_max) for _ in te_family]
    print(f"Nested TE insertions proportion range: {nest_min} - {nest_max}")

    # Create output directory.
    final_out = Path(out_dir) / f"TEgenomeSimulator_{prefix}_result"
    final_out.mkdir(parents=True, exist_ok=True)

    # Write the table.
    table_path = final_out / "TElib_sim_list.table"
    with table_path.open("w") as table_out:
        header = ["#TE_family", "superfamily", "subclass", "count", "idn", "sd", "indels", "tsd", "length", "frag", "nest"]
        table_out.write("\t".join(header) + "\n")
        for i in range(len(te_family)):
            row = [
                te_family[i],
                te_superfamily[i],
                te_subclass[i],
                str(copy_number[i]),
                str(identity[i]),
                str(sd[i]),
                str(indel[i]),
                tsd[i],
                str(length[i]),
                str(fragment[i]),
                str(nested[i])
            ]
            table_out.write("\t".join(row) + "\n")

    print(f"\nGenerated the TE library table for simulation. File saved as {table_path}\n")

if __name__ == "__main__":
    main()
    
# END.
