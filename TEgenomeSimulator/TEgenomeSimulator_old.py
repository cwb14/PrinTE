import os
import sys
import argparse
import subprocess
from pathlib import Path

# Function to call the prep_sim_TE_lib.py script to generate the TE library table.
def run_prep_sim_TE_lib(args, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/prep_sim_TE_lib.py')
    try:
        # Parse copy_number range.
        try:
            mincp, maxcp = map(int, args.copy_number.split(','))
        except ValueError:
            print("Error: --copy_number must be in the format 'min,max' (e.g., '1,5').")
            sys.exit(1)
        
        # Base command.
        prep_telib_command = [
            'python3', script_path,
            '-p', args.prefix,
            '-r', args.repeat,
            '-m', str(maxcp),
            '-n', str(mincp),
            '-s', str(args.seed),
            '-o', args.outdir
        ]
        
        # Add optional range parameters if provided.
        if args.idn_range:
            try:
                idn_min, idn_max = args.idn_range.split(',')
                prep_telib_command.extend(['--idn_range', idn_min, idn_max])
            except ValueError:
                print("Error: --idn_range must be in the format 'min,max' (e.g., '80,95').")
                sys.exit(1)
        
        if args.sd_range:
            try:
                sd_min, sd_max = args.sd_range.split(',')
                prep_telib_command.extend(['--sd_range', sd_min, sd_max])
            except ValueError:
                print("Error: --sd_range must be in the format 'min,max' (e.g., '1,20').")
                sys.exit(1)
        
        if args.indel_range:
            try:
                indel_min, indel_max = args.indel_range.split(',')
                prep_telib_command.extend(['--indel_range', indel_min, indel_max])
            except ValueError:
                print("Error: --indel_range must be in the format 'min,max' (e.g., '5,20').")
                sys.exit(1)
        
        if args.frag_range:
            try:
                frag_min, frag_max = args.frag_range.split(',')
                prep_telib_command.extend(['--frag_range', frag_min, frag_max])
            except ValueError:
                print("Error: --frag_range must be in the format 'min,max' (e.g., '50,98').")
                sys.exit(1)
        
        if args.nest_range:
            try:
                nest_min, nest_max = args.nest_range.split(',')
                prep_telib_command.extend(['--nest_range', nest_min, nest_max])
            except ValueError:
                print("Error: --nest_range must be in the format 'min,max' (e.g., '0,30').")
                sys.exit(1)
        
        # Run the command and capture the output.
        with open(f"{final_out}/TEgenomeSimulator.log", "w") as log_file:
            subprocess.run(prep_telib_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"TE library table generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running prep_sim_TE_lib.py: {e}")
        sys.exit(1)

# Function to call the prep_yml_config.py script for Random Genome Mode.
def run_prep_config_random(args, te_table, final_out, mode):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/prep_yml_config.py')
    try:
        # Base command.
        prep_yml_command = [
            'python3', script_path,
            '-p', args.prefix,
            '-c', args.chridx,
            '-r', args.repeat,
            '-t', te_table,
            '-s', str(args.seed),
            '-o', args.outdir
        ]

        # Run the command and append the output to the log file.
        with open(f"{final_out}/TEgenomeSimulator.log", "a") as log_file:
            subprocess.run(prep_yml_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"Config file generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
        
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running prep_yml_config.py: {e}")
        sys.exit(1)

# Function to call the prep_yml_config.py script for Custom Genome Mode.
def run_prep_config_custom(args, te_table, final_out, mode):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/prep_yml_config.py')
    try:
        # Base command.
        prep_yml_command = [
            'python3', script_path,
            '-p', args.prefix,
            '-g', args.genome,
            '-r', args.repeat,
            '-t', te_table,
            '-s', str(args.seed),
            '-o', args.outdir
        ]
        
        # Run the command and append the output to the log file.
        with open(f"{final_out}/TEgenomeSimulator.log", "a") as log_file:
            subprocess.run(prep_yml_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"Config file generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running prep_yml_config.py: {e}")
        sys.exit(1)

# Function to call TE_sim_random_insertion.py for non-overlap random TE insertion.
def run_TE_sim_random_insertion(args, final_out, mode):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/TE_sim_random_insertion.py')
    try:
        # Base command.
        prep_sim_command = [
            'python3', script_path,
            '-M', str(mode),
            '-p', args.prefix,
            '-o', args.outdir
        ]

        # Add optional indel_size parameter if provided.
        if args.indel_size:
            prep_sim_command.extend(['-i', args.indel_size])
    
        # Run the command and capture the output.
        with open(f"{final_out}/TEgenomeSimulator.log", "a") as log_file:
            subprocess.run(prep_sim_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"Genome with non-overlap random TE insertions was generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running TE_sim_random_insertion.py: {e}")
        sys.exit(1)

# Function to call TE_sim_nested_insertion.py for nested TE insertion.
def run_TE_sim_nested_insertion(args, final_out, mode):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/TE_sim_nested_insertion.py')
    try:
        # Base command.
        prep_nest_command = [
            'python3', script_path,
            '-M', str(mode),
            '-p', args.prefix,
            '-o', args.outdir
        ]

        # Add optional indel_size parameter if provided.
        if args.indel_size:
            prep_nest_command.extend(['-i', args.indel_size])
    
        # Run the command and capture the output.
        with open(f"{final_out}/TEgenomeSimulator.log", "a") as log_file:
            subprocess.run(prep_nest_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"Genome with non-overlap random and nested TE insertions was generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running TE_sim_nested_insertion.py: {e}")
        sys.exit(1)

def main():
    # Set up argument parser.
    parser = argparse.ArgumentParser(
        description="TEgenomeSimulator: Simulate TE mutation and insertion into genome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments.
    required = parser.add_argument_group('required arguments')
    required.add_argument('-p', '--prefix', type=str, help="Prefix for output files.", required=True)
    required.add_argument('-r', '--repeat', type=str, help="TE family fasta file.", required=True)
    required.add_argument('-o', '--outdir', type=str, help="Output directory.", required=True)

    # Optional arguments.
    optional = parser.add_argument_group('optional arguments')
    # Removed mode argument
    # optional.add_argument('-M', '--mode', type=mode_check, choices=[0, 1], default=0,
    #                       help="Mode for genome simulation (either 0 or 1). 0 is simulated genome. 1 is custom genome. If 1, then the --genome parameter is required.")
    optional.add_argument('-cp', '--copy_number', type=str, default='1,5',
                          help="Number of copies of each TE family, in the format 'min,max'.")
    optional.add_argument('-c', '--chridx', type=str,
                          help="Chromosome index file for custom chromosome features if using random genome mode.")
    optional.add_argument('-g', '--genome', type=str,
                          help="Genome fasta file for custom genome mode. Be sure to remove existing TEs.")
    optional.add_argument('-s', '--seed', type=int, default=1,
                          help="Random seed.")
    
    # New optional range parameters for prep_sim_TE_lib.py.
    optional.add_argument('-idn', '--idn_range', type=str, default='80,95',
                          help="Averaged sequence identity range for TE families, in the format 'min,max'.")
    optional.add_argument('-sd', '--sd_range', type=str, default='1,20',
                          help="Standard deviation range of averaged sequence identity, in the format 'min,max'.")
    optional.add_argument('-in', '--indel_range', type=str, default='5,20',
                          help="Proportion range of INDEL to total SNP for each TE family, in the format 'min,max'.")
    optional.add_argument('-fr', '--frag_range', type=str, default='50,98',
                          help="Proportion range of fragmented TE loci for each TE family, in the format 'min,max'.")
    optional.add_argument('-ne', '--nest_range', type=str, default='0,30',
                          help="Proportion range of nested TE insertions for Copia or Gypsy families, in the format 'min,max'.")
    
    # New optional parameters for prep_yml_config.py.
    optional.add_argument('-gs', '--genome_size', type=str, default='300Mb',
                          help="Total genome size (e.g., '300Mb', '500kb', '2Gb').")
    optional.add_argument('-cn', '--chr_number', type=int, default=5,
                          help="Number of chromosomes.")
    optional.add_argument('-gc', '--gc_content', type=int, default=40,
                          help="GC content percentage for all chromosomes.")
    
    # New optional parameter for TE_sim_* scripts.
    optional.add_argument('-is', '--indel_size', type=str, default='1,5',
                          help="Indel size range, in the format 'min,max'.")

    # Parse arguments.
    args = parser.parse_args()

    # Determine mode based on the presence of --genome.
    mode = 1 if args.genome else 0

    # Mode-based validation.
    if mode == 0:
        # Mode 0 can use either --chridx or --genome_size, --chr_number, --gc_content.
        if not args.chridx and not (args.genome_size and args.chr_number and args.gc_content):
            print("Error: When a custom genome is not provided, you must provide either --chridx or all of --genome_size, --chr_number, and --gc_content.")
            sys.exit(1)
    elif mode == 1:
        # Mode 1 requires --genome to be specified.
        if not args.genome:
            print("Error: When a custom genome is provided, --genome must also be specified.")
            sys.exit(1)

    # Output parsed arguments (for demonstration).
    print(f"Mode: {mode}")
    print(f"Prefix: {args.prefix}")
    print(f"Repeat: {args.repeat}")
    print(f"Copy Number Range: {args.copy_number}")
    print(f"Chromosome Index: {args.chridx}")
    print(f"Genome File: {args.genome}")
    print(f"Genome Size: {args.genome_size}")
    print(f"Chromosome Number: {args.chr_number}")
    print(f"GC Content: {args.gc_content}")
    print(f"Seed: {args.seed}")
    print(f"IDN Range: {args.idn_range}")
    print(f"SD Range: {args.sd_range}")
    print(f"INDEL Range: {args.indel_range}")
    print(f"Frag Range: {args.frag_range}")
    print(f"Nest Range: {args.nest_range}")
    print(f"Indel Size: {args.indel_size}")

    # Specify output dir for each project.
    final_out = os.path.join(args.outdir, f"TEgenomeSimulator_{args.prefix}_result")
    Path(final_out).mkdir(parents=True, exist_ok=True)

    # Call the prep_sim_TE_lib.py script to generate the TE library table.
    run_prep_sim_TE_lib(args, final_out)

    # Call the prep_yml_config.py script to generate the config file using the TE library table generated from previous step.
    te_table = os.path.join(final_out, "TElib_sim_list.table")
    if mode == 0:
        print("Mode=0, running prep_yml_config.py for Random Genome Mode.")
        if args.chridx:
            run_prep_config_random(args, te_table, final_out, mode)
        else:
            # User provided genome_size, chr_number, gc_content.
            script_path = os.path.join(os.path.dirname(__file__), 'utils/prep_yml_config.py')
            try:
                prep_yml_command = [
                    'python3', script_path,
                    '-p', args.prefix,
                    '-gs', args.genome_size,
                    '-cn', str(args.chr_number),
                    '-gc', str(args.gc_content),
                    '-r', args.repeat,
                    '-t', te_table,
                    '-s', str(args.seed),
                    '-o', args.outdir
                ]
                with open(f"{final_out}/TEgenomeSimulator.log", "a") as log_file:
                    subprocess.run(prep_yml_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
                print(f"Config file generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while running prep_yml_config.py with genome_size, chr_number, gc_content: {e}")
                sys.exit(1)
    elif mode == 1:
        print("Mode=1, running prep_yml_config.py for Custom Genome Mode.")
        run_prep_config_custom(args, te_table, final_out, mode)

    # Non-overlap random TE insertion.
    run_TE_sim_random_insertion(args, final_out, mode)

    # Nested TE insertion.
    run_TE_sim_nested_insertion(args, final_out, mode)

if __name__ == "__main__":
    main()
