import os
import sys
import argparse
import re
import pandas as pd
from Bio import SeqIO
from pathlib import Path

def parse_genome_size(size_str):
    """
    Parse genome size strings like '300Mb', '500kb', '2Gb' into integer base pair counts.
    """
    size_str = size_str.lower().strip()
    match = re.match(r'^(\d+)([kmg]b)$', size_str)
    if not match:
        raise ValueError(f"Invalid genome_size format: '{size_str}'. Expected formats like '300Mb', '500kb', '2Gb'.")
    number, unit = match.groups()
    number = int(number)
    if unit == 'kb':
        return number * 1_000
    elif unit == 'mb':
        return number * 1_000_000
    elif unit == 'gb':
        return number * 1_000_000_000
    else:
        raise ValueError(f"Unknown unit in genome_size: '{unit}'.")

def main():
    # Parse arguments.
    parser = argparse.ArgumentParser(
        description="Prepare TEgenomeSimulator config file. \
                     You can provide a genome fasta file (-g), a chromosome index file (-c), or specify genome size (-gs), \
                     number of chromosomes (-cn), and GC content (-gc). If a genome fasta is provided, it takes precedence."
    )
    parser.add_argument("-p","--prefix", type=str, required=True,
                        help="Project prefix of the simulation")
    parser.add_argument("-c","--chridx", type=str,
                        help="Path to the chromosome index file (CSV with columns: chr_id, length, gc)")
    parser.add_argument("-gs", "--genome_size", type=str, default="300Mb",
                        help="Total genome size (e.g., 300Mb, 500kb, 2Gb) (default: 300Mb)")
    parser.add_argument("-cn", "--chr_number", type=int, default=5,
                        help="Number of chromosomes (default: 5)")
    parser.add_argument("-gc", "--gc_content", type=float, default=40.0,
                        help="GC content percentage for all chromosomes (default: 40)")
    parser.add_argument("-g","--genome", type=str,
                        help="Path to the genome fasta file")
    parser.add_argument("-r","--repeat", type=str, required=True,
                        help="Path to the repeat fasta file")
    parser.add_argument("-t","--table", type=str, required=True,
                        help="Path to the repeat table file")
    parser.add_argument("-s","--seed", type=str, default="1",
                        help="Seed value for simulation. Default = 1")
    parser.add_argument("-o", "--outdir", type=str, required=True,
                        help="Output directory")

    args = parser.parse_args()
    prefix = args.prefix
    chr_index = args.chridx
    genome_fa = args.genome
    te_fa = args.repeat
    te_table = args.table
    seed = args.seed
    outdir = args.outdir
    genome_size_str = args.genome_size
    chr_number = args.chr_number
    gc_content = args.gc_content

    # Ensure output directory exists.
    Path(outdir).mkdir(parents=True, exist_ok=True)

    print("\n")
    print("#############################################")
    print("### Prepare TEgenomeSimulator config file ###")
    print("#############################################")
    print(f"Using genome fasta file: {args.genome if args.genome else 'Not provided'}")
    print(f"Using repeat fasta file: {args.repeat}")
    print(f"Output directory set as: {args.outdir}")
    print("\n")

    # Initialize chromosome dictionary.
    genome_chrs = {}

    # Flag to indicate which mode is being used.
    use_genome_fa = False
    use_chridx = False
    use_alternative = False

    if genome_fa:
        use_genome_fa = True
        try:
            # Read genome fasta file and extract chromosome information.
            with open(genome_fa, 'r') as fasta_handle:
                for record in SeqIO.parse(fasta_handle, 'fasta'):
                    chr_id = record.id
                    seq_length = len(record.seq)
                    genome_chrs[chr_id] = {
                        'prefix': chr_id,
                        'seq_length': seq_length
                    }
            if not genome_chrs:
                raise ValueError("Genome fasta file is empty or improperly formatted.")
        except Exception as e:
            print(f"Error: Failed to process genome fasta file '{genome_fa}'. Reason: {e}")
            sys.exit(1)
    elif chr_index:
        use_chridx = True
        try:
            # Read chromosome index file.
            df = pd.read_csv(str(chr_index), sep=',', header=None, names=['chr_id', 'length', 'gc'])
            # Validate required columns.
            if not {'chr_id', 'length', 'gc'}.issubset(df.columns):
                raise ValueError("Chromosome index file must contain 'chr_id', 'length', and 'gc' columns.")
            # Validate data types.
            if not pd.api.types.is_numeric_dtype(df['length']):
                raise ValueError("'length' column must be numeric.")
            if not pd.api.types.is_numeric_dtype(df['gc']):
                raise ValueError("'gc' column must be numeric.")
            # Convert to dict.
            for _, row in df.iterrows():
                chr_id = row['chr_id']
                length = int(row['length'])
                gc = row['gc']
                genome_chrs[chr_id] = {
                    'prefix': chr_id,
                    'seq_length': length,
                    'gc_content': gc
                }
        except Exception as e:
            print(f"Warning: Failed to process chromosome index file '{chr_index}'. Reason: {e}")
            print("Falling back to alternative parameters: genome_size, chr_number, gc_content.\n")
            use_chridx = False  # Switch to alternative parameters.

    if not use_genome_fa and not use_chridx:
        use_alternative = True
        try:
            # Parse genome size.
            total_genome_size_bp = parse_genome_size(genome_size_str)
            # Calculate uniform chromosome length.
            uniform_chr_length = int(total_genome_size_bp / chr_number)
            # Validate GC content.
            if not (0 <= gc_content <= 100):
                raise ValueError("GC content must be between 0 and 100.")
            # Generate chromosome dictionary with uniform lengths and GC content.
            for i in range(1, chr_number + 1):
                chr_id = f'chr{i}'
                genome_chrs[chr_id] = {
                    'prefix': chr_id,
                    'seq_length': uniform_chr_length,
                    'gc_content': gc_content
                }
        except Exception as e:
            print(f"Error: Failed to use alternative parameters. Reason: {e}")
            sys.exit(1)

    # Create the config file path.
    final_out = Path(outdir) / f'TEgenomeSimulator_{prefix}_result'
    final_out.mkdir(parents=True, exist_ok=True)
    yml_name = f'TEgenomeSimulator_{prefix}.yml'
    config_path = final_out / yml_name

    with open(config_path, 'w') as config_out:
        config_out.write(f'prefix: "{prefix}"\n')
        if use_genome_fa:
            config_out.write(f'genome_fasta: "{genome_fa}"\n')
        config_out.write(f'rep_fasta: "{te_fa}"\n')
        config_out.write(f'rep_list: "{te_table}"\n')
        config_out.write(f'seed: {int(seed)}\n')
        config_out.write('chrs:\n')
        for chr_id, data in genome_chrs.items():
            config_out.write(f'  {chr_id}:\n')
            config_out.write(f'    prefix: "{data["prefix"]}"\n')
            config_out.write(f'    seq_length: {int(data["seq_length"])}\n')
            if use_chridx or use_alternative:
                # Only include gc_content if not using genome_fa
                gc_ctn = data.get('gc_content') if use_chridx or use_alternative else None
                if gc_ctn is not None:
                    config_out.write(f'    gc_content: {gc_ctn}\n')

    print(f"Generated the config file for simulation. File saved as {config_path}")
    print("\n")

if __name__ == "__main__":
    main()
