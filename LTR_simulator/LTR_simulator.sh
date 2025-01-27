#!/bin/bash

# ====================================================================
# LTR_simulator.sh
# 
# A wrapper script to simulate generations of DNA mutation over genomes,
# mutate LTR files, and insert new solos and intact LTRs into the mutated genomes.
# 
# ====================================================================

# Dynamically determine the directory of the current script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# Set the path to the 'bin' directory relative to the script's location
BIN="$SCRIPT_DIR/bin"

# Default values
threads=10
mutation_rate=1.3e-8
LTR_insertion_rate=1e-11
solo_intact_ratio=0:1
unique_LTRs="$SCRIPT_DIR/maize_LTRs.fa.gz"  # Updated default unique_LTRs
# Set default CDS to the gzipped file
CDS="$SCRIPT_DIR/Amborella.AMTR1.GCF_000471905.cdhit.cds.gz"
genome_size="300Mb"
chr_number=5

# Variables to store required arguments
shared_LTRs=""
shared_LTR_count=""
generation_high=""
generation_low=""
generation_step=""

# Function to display help
show_help() {
    cat << EOF
Usage: bash LTR_simulator.sh [OPTIONS]

A wrapper script to simulate generations of DNA mutation over genomes,
mutate LTR files, and insert new solos and intact LTRs into the mutated genomes.

Required Arguments:
  -sLc, --shared_LTR_count     Number of shared LTRs to insert into the genomes.
  -h, --generation_high_range  Highest generation number to simulate.
  -l, --generation_low_range   Lowest generation number to simulate.
  -s, --generation_step        Step size for generations.

Optional Arguments:
  -sL, --shared_LTRs           Path to the shared LTRs FASTA file.
                               If not provided, shared LTRs will be generated from unique_LTRs.
  -gs, --genome_size           Total genome size (e.g., 300Mb, 500kb, 2Gb) (default: 300Mb).
  -cds, --CDS                  Path to the CDS file (default: \$CDS).
  -cn, --chr_number            Number of chromosomes (default: 5).
  -t, --threads                Number of threads to use (default: 10).
  -mr, --mutation_rate         Mutation rate (default: 1.3e-8).
  -lr, --LTR_insertion_rate    LTR insertion rate (default: 1e-11).
  -r, --solo_intact_ratio      Ratio of solo to intact LTRs (default: 0:1). 
                               Note: Solo support requires '_solo' suffix in the '-uL, --unique_LTRs' fasta headers.
  -uL, --unique_LTRs           Path to the unique LTRs FASTA file (default: \$unique_LTRs).
  --help                       Display this help message and exit.

Example:
  bash LTR_simulator.sh -sL shared_LTRs.fa -sLc 8200 -h 50000 -l 10000 -s 20000

  Alternatively:
  bash LTR_simulator.sh --shared_LTRs shared_LTRs.fa --shared_LTR_count 8200 --threads 10 \\
      --mutation_rate 1.3e-8 --LTR_insertion_rate 1e-11 \\
      --generation_high_range 1000000 --generation_low_range 10000 \\
      --generation_step 50000 --solo_intact_ratio 1:1 \\
      --unique_LTRs /path/to/custom_LTR.fa \\
      --genome_size 500Mb --CDS /path/to/CDS.fa --chr_number 5

Notes:
  - The script will simulate DNA mutations for each specified generation range.
  - LTR mutations and insertions are performed for each generation.
  - Mutating genomes is parallelized to utilize multiple threads effectively.
EOF
}

# If no arguments are provided, display help
if [[ $# -eq 0 ]]; then
    show_help
    exit 0
fi

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -gs|--genome_size)
            genome_size="$2"
            shift
            shift
            ;;
        -cds|--CDS)
            CDS="$2"
            shift
            shift
            ;;
        -cn|--chr_number)
            chr_number="$2"
            shift
            shift
            ;;
        -sL|--shared_LTRs)
            shared_LTRs="$2"
            shift
            shift
            ;;
        -sLc|--shared_LTR_count)
            shared_LTR_count="$2"
            shift
            shift
            ;;
        -t|--threads)
            threads="$2"
            shift
            shift
            ;;
        -mr|--mutation_rate)
            mutation_rate="$2"
            shift
            shift
            ;;
        -lr|--LTR_insertion_rate)
            LTR_insertion_rate="$2"
            shift
            shift
            ;;
        -h|--generation_high_range)
            generation_high="$2"
            shift
            shift
            ;;
        -l|--generation_low_range)
            generation_low="$2"
            shift
            shift
            ;;
        -s|--generation_step)
            generation_step="$2"
            shift
            shift
            ;;
        -r|--solo_intact_ratio)
            solo_intact_ratio="$2"
            shift
            shift
            ;;
        -uL|--unique_LTRs)
            unique_LTRs="$2"
            shift
            shift
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)    # unknown option
            echo "Unknown option: $1"
            echo "Use --help to see available options."
            exit 1
            ;;
    esac
done

# Check required arguments
if [[ -z "$shared_LTR_count" ]]; then
    echo "Error: shared_LTR_count (-sLc or --shared_LTR_count) is required."
    echo "Use --help to see usage."
    exit 1
fi

if [[ -z "$generation_high" ]]; then
    echo "Error: generation_high_range (-h or --generation_high_range) is required."
    echo "Use --help to see usage."
    exit 1
fi

if [[ -z "$generation_low" ]]; then
    echo "Error: generation_low_range (-l or --generation_low_range) is required."
    echo "Use --help to see usage."
    exit 1
fi

if [[ -z "$generation_step" ]]; then
    echo "Error: generation_step (-s or --generation_step) is required."
    echo "Use --help to see usage."
    exit 1
fi

# Handle default CDS: decompress if default CDS is used
default_CDS="$SCRIPT_DIR/Amborella.AMTR1.GCF_000471905.cdhit.cds.gz"
expected_default_CDS="$SCRIPT_DIR/Amborella.AMTR1.GCF_000471905.cdhit.cds.gz"

if [[ "$CDS" == "$default_CDS" ]]; then
    if [[ -f "$default_CDS" ]]; then
        echo "Decompressing default CDS file: $default_CDS"
        # Decompress the CDS file while keeping the original gzipped file
        gunzip -c "$default_CDS" > "$expected_default_CDS"
        # Update CDS to point to the decompressed file
        CDS="$expected_default_CDS"
        echo "Decompressed CDS file: $CDS"
    else
        echo "Error: Default CDS compressed file not found: $default_CDS"
        exit 1
    fi
fi

# Handle default unique_LTRs: decompress if default unique_LTRs is used
default_unique_LTRs="$SCRIPT_DIR/maize_LTRs.fa.gz"
expected_default_unique_LTRs="$SCRIPT_DIR/maize_LTRs.fa"

if [[ "$unique_LTRs" == "$default_unique_LTRs" ]]; then
    if [[ -f "$default_unique_LTRs" ]]; then
        echo "Decompressing default unique LTRs file: $default_unique_LTRs"
        # Decompress the unique_LTRs file while keeping the original gzipped file
        gunzip -c "$default_unique_LTRs" > "$expected_default_unique_LTRs"
        # Update unique_LTRs to point to the decompressed file
        unique_LTRs="$expected_default_unique_LTRs"
        echo "Decompressed unique LTRs file: $unique_LTRs"
    else
        echo "Error: Default unique LTRs compressed file not found: $default_unique_LTRs"
        exit 1
    fi
fi

# If shared_LTRs is not provided, create it using shared_ltr_mutator2.py
if [[ -z "$shared_LTRs" ]]; then
    echo "No shared_LTRs provided. Generating shared LTRs from unique_LTRs..."
    shared_LTRs="$SCRIPT_DIR/shared_maize_LTRs.fa"  # Define output path for shared_LTRs
    cmd_create_shared_LTRs="python $BIN/shared_ltr_mutator.py -fasta \"$unique_LTRs\" -generations 1 -max_perc_div 10 -out \"$shared_LTRs\" -rate 1.3e-8 -shape 0.0000001 --multiplier 10"
    echo "Creating shared LTRs: $cmd_create_shared_LTRs"
    python "$BIN/shared_ltr_mutator.py" -fasta "$unique_LTRs" -generations 1 -max_perc_div 10 -out "$shared_LTRs" -rate 1.3e-8 -shape 0.0000001 --multiplier 10
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to create shared_LTRs using shared_ltr_mutator2.py"
        exit 1
    fi
    echo "Shared LTRs created at: $shared_LTRs"
else
    # If shared_LTRs is provided, ensure it exists
    if [[ ! -f "$shared_LTRs" ]]; then
        echo "Error: Provided shared_LTRs file does not exist: $shared_LTRs"
        exit 1
    fi
fi

# Print parameters
echo "========================================"
echo "LTR Simulator Parameters:"
echo "----------------------------------------"
echo "Genome Size: $genome_size"
echo "CDS: $CDS"
echo "Chromosome Number: $chr_number"
echo "Threads: $threads"
echo "Mutation Rate: $mutation_rate"
echo "LTR Insertion Rate: $LTR_insertion_rate"
echo "Generation High Range: $generation_high"
echo "Generation Low Range: $generation_low"
echo "Generation Step: $generation_step"
echo "Solo:Intact Ratio: $solo_intact_ratio"
echo "Unique LTRs: $unique_LTRs"
echo "Shared LTRs: $shared_LTRs"
echo "Shared LTR Count: $shared_LTR_count"
echo "========================================"

# Step 1: Tag the young LTRs
unique_LTRs_tagged="${unique_LTRs%.*}_tagged.fa"
cmd_tag_unique_LTRs="python $BIN/tag_inserter.py -i \"$unique_LTRs\" -o \"$unique_LTRs_tagged\""
echo "Tagging unique LTRs: $cmd_tag_unique_LTRs"
python "$BIN/tag_inserter.py" -i "$unique_LTRs" -o "$unique_LTRs_tagged"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to tag unique LTRs."
    exit 1
fi

# Step 2: Build synthetic genome
synthetic_genome="A0.fa"
cmd_build_synthetic_genome="python $BIN/synthetic_genome.py -size \"$genome_size\" -cds \"$CDS\" -chr_number \"$chr_number\" > \"$synthetic_genome\""
echo "Building synthetic genome: $cmd_build_synthetic_genome"
python "$BIN/synthetic_genome.py" -size "$genome_size" -cds "$CDS" -chr_number "$chr_number" > "$synthetic_genome"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to build synthetic genome."
    exit 1
fi

# Step 3: Add shared LTRs
genome_with_shared_LTRs="A${shared_LTR_count}s.fa"
cmd_insert_shared_LTRs="python $BIN/shared_ltr_inserter.py -genome \"$synthetic_genome\" -LTR \"$shared_LTRs\" -n \"$shared_LTR_count\" -output \"$genome_with_shared_LTRs\""
echo "Inserting shared LTRs: $cmd_insert_shared_LTRs"
python "$BIN/shared_ltr_inserter.py" -genome "$synthetic_genome" -LTR "$shared_LTRs" -n "$shared_LTR_count" -output "$genome_with_shared_LTRs"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to insert shared LTRs."
    exit 1
fi

# Step 4: Make two genomes A and B
genome1="$genome_with_shared_LTRs"
genome2="B${shared_LTR_count}s.fa"
cmd_copy_genome="cp \"$genome_with_shared_LTRs\" \"$genome2\""
echo "Creating genome2 by copying genome1: $cmd_copy_genome"
cp "$genome_with_shared_LTRs" "$genome2"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to copy genome to create genome2."
    exit 1
fi

# Prepare genome prefixes
genome1_basename=$(basename "$genome1")
genome1_prefix="${genome1_basename%.*}"

genome2_basename=$(basename "$genome2")
genome2_prefix="${genome2_basename%.*}"

# Generate list of generations
current_gen="$generation_low"
generations=()

while true
do
    comparison=$(echo "$current_gen <= $generation_high" | bc)
    if [ "$comparison" -eq 1 ]; then
        generations+=("$current_gen")
        current_gen=$(echo "$current_gen + $generation_step" | bc)
    else
        break
    fi
done

echo "Generations to process: ${generations[@]}"
echo "========================================"

# Validate generation numbers
for gen in "${generations[@]}"
do
    if ! [[ "$gen" =~ ^[0-9]+$ ]]; then
        echo "Error: Generation number '$gen' is not a positive integer."
        exit 1
    fi
done

# Function to generate generation shorthand
generate_shorthand() {
    local gen=$1
    local shorthand=""

    if (( gen >= 1000000 )); then
        # Calculate millions with one decimal if necessary
        local millions=$(echo "scale=1; $gen / 1000000" | bc)
        # Remove trailing .0 for whole numbers
        if [[ "$millions" == *.0 ]]; then
            millions=${millions%.0}
        fi
        shorthand="${millions}m"
    elif (( gen >= 100000 )); then
        # Calculate hundreds of thousands with one decimal if necessary
        local hundreds=$(echo "scale=1; $gen / 100000" | bc)
        if [[ "$hundreds" == *.0 ]]; then
            hundreds=${hundreds%.0}
        fi
        shorthand="${hundreds}ht"
    elif (( gen >= 1000 )); then
        # Calculate thousands with one decimal if necessary
        local thousands=$(echo "scale=1; $gen / 1000" | bc)
        if [[ "$thousands" == *.0 ]]; then
            thousands=${thousands%.0}
        fi
        shorthand="${thousands}t"
    else
        shorthand="${gen}"
    fi

    echo "$shorthand"
}

# Start genome mutation steps
echo "Starting genome mutation steps..."
echo "----------------------------------------"

mutation_pids=()

for gen in "${generations[@]}"
do
    # Mutate genome1
    cmd_mutate_genome1="$BIN/ltr_mutator -fasta \"$genome1\" -rate \"$mutation_rate\" -generations \"$gen\" -mode 0 -threads \"$threads\""
    echo "Running: $cmd_mutate_genome1"
    "$BIN/ltr_mutator" -fasta "$genome1" -rate "$mutation_rate" -generations "$gen" -mode 0 -threads "$threads" > /dev/null 2>&1 &
    mutation_pids+=($!)

    # Mutate genome2
    cmd_mutate_genome2="$BIN/ltr_mutator -fasta \"$genome2\" -rate \"$mutation_rate\" -generations \"$gen\" -mode 0 -threads \"$threads\""
    echo "Running: $cmd_mutate_genome2"
    "$BIN/ltr_mutator" -fasta "$genome2" -rate "$mutation_rate" -generations "$gen" -mode 0 -threads "$threads" > /dev/null 2>&1 &
    mutation_pids+=($!)
done

# Wait for genome mutation steps to finish
echo "Waiting for genome mutation steps to finish..."
wait "${mutation_pids[@]}"
echo "Genome mutation steps completed."
echo "========================================"

# Start LTR mutation and insertion steps
echo "Starting LTR mutation and insertion steps..."
echo "----------------------------------------"

inserter_pids=()

for gen in "${generations[@]}"
do
    echo "Processing generation $gen"

    # Mutate the LTR file
    cmd_mutate_ltr="$BIN/ltr_mutator -fasta \"$unique_LTRs_tagged\" -rate \"$mutation_rate\" -generations \"$gen\" -mode 1 -threads \"$threads\""
    echo "Running: $cmd_mutate_ltr"
    "$BIN/ltr_mutator" -fasta "$unique_LTRs_tagged" -rate "$mutation_rate" -generations "$gen" -mode 1 -threads "$threads" > /dev/null 2>&1

    # Check if LTR mutation was successful
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to mutate LTRs for generation $gen."
        exit 1
    fi

    # Generate shorthand
    shorthand=$(generate_shorthand "$gen")

    # Create shorthand directory
    mkdir -p "$shorthand"

    # Prepare file names
    genome1_mutated="${genome1_prefix}_${gen}.fa"
    output_genome1="${shorthand}/${genome1_prefix}${shorthand}.fa"

    genome2_mutated="${genome2_prefix}_${gen}.fa"
    output_genome2="${shorthand}/${genome2_prefix}${shorthand}.fa"

    ltr_mutated="${unique_LTRs_tagged%.*}_${gen}.fa"

    # Insert LTRs into genome1
    cmd_insert_genome1="$BIN/ltr_inserter --genome \"$genome1_mutated\" --ltr \"$ltr_mutated\" --rate \"$LTR_insertion_rate\" --generations \"$gen\" --solo_intact_ratio \"$solo_intact_ratio\" --output \"$output_genome1\" --threads \"$threads\""
    echo "Running: $cmd_insert_genome1"
    "$BIN/ltr_inserter" --genome "$genome1_mutated" --ltr "$ltr_mutated" --rate "$LTR_insertion_rate" --generations "$gen" --solo_intact_ratio "$solo_intact_ratio" --output "$output_genome1" --threads "$threads" &
    inserter_pids+=($!)

    # Insert LTRs into genome2
    cmd_insert_genome2="$BIN/ltr_inserter --genome \"$genome2_mutated\" --ltr \"$ltr_mutated\" --rate \"$LTR_insertion_rate\" --generations \"$gen\" --solo_intact_ratio \"$solo_intact_ratio\" --output \"$output_genome2\" --threads \"$threads\""
    echo "Running: $cmd_insert_genome2"
    "$BIN/ltr_inserter" --genome "$genome2_mutated" --ltr "$ltr_mutated" --rate "$LTR_insertion_rate" --generations "$gen" --solo_intact_ratio "$solo_intact_ratio" --output "$output_genome2" --threads "$threads" &
    inserter_pids+=($!)
done

# Wait for all ltr_inserter processes to finish
echo "Waiting for LTR insertion steps to finish..."
wait "${inserter_pids[@]}"
echo "LTR mutation and insertion steps completed."
echo "========================================"

echo "All steps completed successfully."

# ====================================================================
# Cleanup Steps
# ====================================================================

echo "Starting cleanup steps..."
echo "----------------------------------------"

# Remove synthetic_genome, genome_with_shared_LTRs, genome2
echo "Removing intermediate genome files..."
rm -f "$synthetic_genome" "$genome_with_shared_LTRs" "$genome2"

# Loop through each generation to handle .fa and .txt files
for gen in "${generations[@]}"
do
    echo "Cleaning up files for generation $gen..."

    # Remove .fa files
    fa_file_A="A${shared_LTR_count}s_${gen}.fa"
    fa_file_B="B${shared_LTR_count}s_${gen}.fa"
    if [[ -f "$fa_file_A" ]]; then
        echo "Removing $fa_file_A"
        rm -f "$fa_file_A"
    fi
    if [[ -f "$fa_file_B" ]]; then
        echo "Removing $fa_file_B"
        rm -f "$fa_file_B"
    fi

    # Handle .txt files
    txt_file_A="A${shared_LTR_count}s_${gen}.txt"
    txt_file_B="B${shared_LTR_count}s_${gen}.txt"

    # Generate shorthand
    shorthand=$(generate_shorthand "$gen")

    # Define target txt files
    target_txt_A="${shorthand}/${genome1_prefix}${shorthand}.txt"
    target_txt_B="${shorthand}/${genome2_prefix}${shorthand}.txt"

    # Append contents if txt files exist
    if [[ -f "$txt_file_A" ]]; then
        echo "Appending $txt_file_A to $target_txt_A"
        cat "$txt_file_A" >> "$target_txt_A"
        echo "Removing $txt_file_A"
        rm -f "$txt_file_A"
    fi

    if [[ -f "$txt_file_B" ]]; then
        echo "Appending $txt_file_B to $target_txt_B"
        cat "$txt_file_B" >> "$target_txt_B"
        echo "Removing $txt_file_B"
        rm -f "$txt_file_B"
    fi
done

echo "Cleanup steps completed."
echo "========================================"

echo "All processes finished successfully."

exit 0

# END.
