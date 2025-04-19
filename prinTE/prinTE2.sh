#!/bin/bash
###############################################################################
# prinTE.sh
#
# A wrapper to run the TE evolution simulation pipeline in two phases,
# and then perform supplemental post-processing.
#
# Phase 1 (Burn‐in) is performed unless the user provides starting input files
# via --bed and --fasta, or if --continue is provided and previous outputs are found.
#
# New functionality and changes:
#
#  (1) New flag --continue:
#      When provided, the script searches for existing final generation outputs
#      (e.g., gen${i}_final.fasta) and resumes from the next iteration.
#
#  (2) New flag --keep_temps (or -kt):
#      When provided, temporary files are kept. Otherwise, they are removed after each loop.
#
#  (3) New user-friendly parameters --chromatin_bias_insert and --chromatin_buffer
#      to control insertion bias. (Replaces old -w/--euch-bias and -j/--euch-buffer.)
#
#  (4) New flag -cbd, --chromatin_bias_delete to specify the deletion bias for TE excision.
#
#  (5) New option --model (and shorthand -md) for per-generation post-processing
#      DNA mutation model for LTR dating (options: raw, K2P, JC69; default: K2P).
#
#  (6) New options -tk/--TE_mut_k and -tmx/--TE_mut_Mmax replace the old TE_mut_in parameter.
#
#  (7) New parallel versions of internal scripts:
#         - Use 'shared_ltr_inserter_parallel.py' instead of 'shared_ltr_inserter.py'
#         - Use 'nest_inserter_parallel2.py' instead of 'nest_inserter_parallel.py', which now adds a '-m' parameter for threads.
#           Also accepts the new optional flag --disable_genes (-dg) to disable insertion into genes.
#
#  (8) New flag -bo, --burnin_only:
#      When activated, the script will run the burn-in phase only and then exit.
#      In this case, --generation_end and --step are not required.
#
#  (9) File naming: The output files from Phase 2 now have names
#      reflecting the true generation simulated.
#
#  (10) Updated TE excision: Use "TE_exciser_parallel_dev2.py" (instead of TE_exciser2.py)
#       which supports parallel execution with -m, and accepts extra parameters:
#         --euch_het_buffer ${euch_buffer} and --euch_het_bias ${euch_bias_excise}.
#
# Directories:
#   TOOL_DIR: Directory containing this script (assumed to be TESS/prinTE)
#   BIN_DIR:  TESS/prinTE/bin
#   UTIL_DIR: TESS/prinTE/util
#
# Usage:
#   prinTE.sh [options]
#
# [A detailed options list follows...]
#
###############################################################################

# --- Determine directories ---
TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
BIN_DIR="${TOOL_DIR}/bin"
UTIL_DIR="${TOOL_DIR}/util"

# --- Temporary files array (for unzipped inputs) ---
temp_files=()
cleanup() {
  for f in "${temp_files[@]}"; do
    [ -f "$f" ] && rm -f "$f"
  done
}
trap cleanup EXIT

# --- Function: If a fasta file is gzipped, decompress it.
decompress_if_gz() {
  local file="$1"
  local base_name="$(basename "$file" .gz)"  # Remove .gz if present
  local decompressed_file="./${base_name}"

  if [[ ! -f "$file" && "$file" != *.gz ]]; then
    if [ -f "${file}.gz" ]; then
      file="${file}.gz"
    fi
  fi

  if [[ "$file" == *.gz ]]; then
    echo "Decompressing $file to ${decompressed_file}" | tee -a "$LOG" >&2
    if gunzip -c "$file" > "$decompressed_file"; then
      temp_files+=("$decompressed_file")
      echo "$decompressed_file"  # Only output the filename
    else
      echo "Error decompressing $file" | tee -a "$ERR" >&2
      exit 1
    fi
  else
    echo "$file"
  fi
}

# --- Help message function ---
print_help() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
############### REQUIRED (for full pipeline) ###############
  -ge, --generation_end      Number of generations to simulate (REQUIRED unless --burnin_only is used; must be an exact multiple of step)
  -st, --step                Generation step size (number of generations per loop; REQUIRED unless --burnin_only is used)

############### BURN-IN ###############
  -c,  --cds                 Path to CDS file (default: ${TOOL_DIR}/TAIR10.cds.fa)
  -N,  --cds_num             Number of CDS sequences to insert. (Mutually exclusive with --cds_percent)
  -P,  --cds_percent         Percent of the genome that should be CDS. (Mutually exclusive with --cds_num)
  -n,  --TE_num              Number of TE insertions in burn‐in (default: 2000) (Mutually exclusive with --TE_percent)
  -p,  --TE_percent          Percent of the genome that should be TEs in the burn‐in (Mutually exclusive with --TE_num)
  -cn, --chr_number          Number of chromosomes (default: 4)
  -sz, --size                Genome size in kb, Mb, or Gb (default: 400Mb)
  -tk, --TE_mut_k            Slope of exponential decay for TE mutation (default: 10)
  -tmx, --TE_mut_Mmax         X-limit for exponential decay function (default: 20)
  -bo, --burnin_only         Run burn-in phase only and then exit.
  -i,  --TE_lib              TE library file (default: ${TOOL_DIR}/combined_curated_TE_lib_ATOSZM_selected.fasta)
  -m,  --mutation_rate       DNA Mutation rate (default: 1.3e-8)

########## FIXED TE INDEL RATE ##########
  -F,  --fix                 Fixed insertion and deletion numbers, comma-separated. Format: insertion,deletion (e.g., 1e-9,1e-9)
  -dg, --disable_genes       Disable insertion into genes (only effective if -F/--fix is provided)

######### VARIABLE TE INDEL RATE #########  
  -ir, --insert_rate         TE insertion rate (default: 1e-8)
  -dr, --delete_rate         TE deletion rate (default: 1e-7)  
  -br, --birth_rate          TE birth rate (default: 1e-3)
  -sc, --sigma               Selection coefficient for gene insertions (default: 1.0)  
  -sf, --sel_coeff           Selection coefficient for TE excision (variable-rate only; default: 0)
                             0 = neutral; 0.1 = 2× bias; 1 = 11× bias.
  -cbi, --chromatin_bias_insert   Chromatin bias for TE insertion (default: 1.0)
  -cbd, --chromatin_bias_delete   Chromatin bias for TE deletion (default: 1.0)
  -cb,  --chromatin_buffer         Interval upstream/downstream used for chromatin bias (default: 10000)

########## GENERAL USE ##########                 
  -s,  --seed                Random seed (default: 42)
  -r,  --TE_ratio            TE ratio file (default: ${TOOL_DIR}/ratios.tsv)
  -t,  --threads             Number of threads (default: 4)
  -sr, --solo_rate           Percent chance to convert an intact TE to soloLTR (default: 95)
  -k,  --k                   TE length decay slope for TE excision (default: 10)
  -b,  --bed                 Input BED file for starting generation (skip burn-in; requires --fasta)
  -f,  --fasta               Input FASTA file for starting generation (skip burn-in; requires --bed)
  -x,  --continue            Resume simulation from the last completed generation.
  -kt, --keep_temps          Keep temporary files (default: remove them after each loop)
  -md, --model               DNA mutation model for LTR dating (choose from: raw, K2P, JC69; default: K2P)
  -ex, --ex_LTR              Exclude LTR sequences without a domain hit (passed as '-exclude_no_hits' to LTR_fasta_header_appender.py)
  -h,  --help                Display this help message and exit

Example:
  $(basename "$0") -st 1000 -ge 3000
EOF
}

#—-----------------------------------------
# if no args or any of -h, --h, -help, --help, print help and exit
if [[ $# -eq 0 \
   || "$1" == "-h"  \
   || "$1" == "--h" \
   || "$1" == "-help" \
   || "$1" == "--help" ]]; then
  print_help
  exit 0
fi
#—-----------------------------------------

# --- Log file names ---
LOG="pipeline.log"
ERR="pipeline.error"
echo "Pipeline started at $(date)" > "$LOG"
echo "Pipeline started at $(date)" > "$ERR"
# --- Modification (1): Log the command used to run the script ---
echo "Command: $0 $@" >> "$LOG"

# --- Parse command-line options ---
# Initialize flags and new parameters with defaults.
cont_flag=0
keep_temps=0
burnin_only=0
euch_bias_insert=1.0
euch_bias_excise=1.0
euch_buffer=10000
model="K2P"
ex_ltr=0
disable_genes=0
# New TE inserter mutation options defaults:
TE_mut_k=10
TE_mut_Mmax=20
sel_coeff=0

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -c|--cds)
      cds="$2"
      shift; shift;;
    -N|--cds_num)
      cds_num="$2"
      shift; shift;;
    -P|--cds_percent)
      cds_percent="$2"
      shift; shift;;
    -cn|--chr_number)
      chr_number="$2"
      shift; shift;;
    -sz|--size)
      size="$2"
      shift; shift;;
    -s|--seed)
      seed="$2"
      shift; shift;;
    -i|--TE_lib)
      TE_lib="$2"
      shift; shift;;
    -m|--mutation_rate)
      mutation_rate="$2"
      shift; shift;;
    -r|--TE_ratio)
      TE_ratio="$2"
      shift; shift;;
    -n|--TE_num)
      TE_num="$2"
      shift; shift;;
    -p|--TE_percent)
      TE_percent="$2"
      shift; shift;;
    -st|--step)
      step="$2"
      shift; shift;;
    -ge|--generation_end)
      generation_end="$2"
      shift; shift;;
    -t|--threads)
      threads="$2"
      shift; shift;;
    -ir|--insert_rate)
      insert_rate="$2"
      shift; shift;;
    -br|--birth_rate)
      birth_rate="$2"
      shift; shift;;
    -dr|--delete_rate)
      delete_rate="$2"
      shift; shift;;
    -sr|--solo_rate)
      solo_rate="$2"
      shift; shift;;
    -k|--k)
      k="$2"
      shift; shift;;
    -F|--fix)
      fix="$2"
      shift; shift;;
    -b|--bed)
      input_bed="$2"
      shift; shift;;
    -f|--fasta)
      input_fasta="$2"
      shift; shift;;
    -x|--continue)
      cont_flag=1
      shift;;
    --keep_temps|-kt)
      keep_temps=1
      shift;;
    -cbi|--chromatin_bias_insert)
      euch_bias_insert="$2"
      shift; shift;;
    -cb|--chromatin_buffer)
      euch_buffer="$2"
      shift; shift;;
    -cbd|--chromatin_bias_delete)
      euch_bias_excise="$2"
      shift; shift;;
    -sc|--sigma)
      sigma="$2"
      shift; shift;;
    -tk|--TE_mut_k)
      TE_mut_k="$2"
      shift; shift;;
    -tmx|--TE_mut_Mmax)
      TE_mut_Mmax="$2"
      shift; shift;;
    -bo|--burnin_only)
      burnin_only=1
      shift;;
    -ex|--ex_LTR)
      ex_ltr=1
      shift;;
    -dg|--disable_genes)
      disable_genes=1
      shift;;
    -sf|--sel_coeff)
      sel_coeff="$2"
      shift; shift;;
    -h|--help)
      print_help
      exit 0;;
    *)
      echo "Unknown option: $1" | tee -a "$ERR"
      print_help
      exit 1;;
  esac
done

# --- Set default values for options not provided ---
cds="${cds:-${TOOL_DIR}/TAIR10.cds.fa}"
chr_number="${chr_number:-4}"
size="${size:-400Mb}"
seed="${seed:-42}"
TE_lib="${TE_lib:-${TOOL_DIR}/combined_curated_TE_lib_ATOSZM_selected.fasta}"
mutation_rate="${mutation_rate:-1.3e-8}"
TE_ratio="${TE_ratio:-${TOOL_DIR}/ratios.tsv}"
TE_num="${TE_num:-2000}"
threads="${threads:-4}"
insert_rate="${insert_rate:-1e-8}"
birth_rate="${birth_rate:-1e-3}"
delete_rate="${delete_rate:-1e-7}"
solo_rate="${solo_rate:-95}"
k="${k:-10}"
sigma="${sigma:-1.0}"

# --- Validate mutually exclusive CDS options ---
if [[ -n "$cds_num" && -n "$cds_percent" ]]; then
  echo "Error: Provide either --cds_num (or -N) or --cds_percent (or -P), not both." | tee -a "$ERR"
  exit 1
fi

# --- Validate mutually exclusive TE insertion options ---
if [[ -z "$TE_percent" ]]; then
  TE_num="${TE_num:-2000}"
fi

# --- Validate required parameters (if not running burn-in only) ---
if [[ "$burnin_only" -eq 0 ]]; then
  if [[ -z "$step" || -z "$generation_end" ]]; then
    echo "Error: Both --step and --generation_end must be provided for the looping phase." | tee -a "$ERR"
    print_help
    exit 1
  fi

  # Ensure generation_end is an exact multiple of step.
  if (( generation_end % step != 0 )); then
    echo "Error: generation_end must be an exact multiple of step." | tee -a "$ERR"
    exit 1
  fi
fi

# --- Check for BED and FASTA input consistency ---
if [[ -n "$input_bed" || -n "$input_fasta" ]]; then
  if [[ -z "$input_bed" || -z "$input_fasta" ]]; then
    echo "Error: Both --bed and --fasta must be provided together." | tee -a "$ERR"
    print_help
    exit 1
  fi
  skip_burnin=1
  echo "User provided BED and FASTA; skipping burn-in phase." | tee -a "$LOG"
else
  skip_burnin=0
fi

# If --continue is provided, force skip burn-in.
if [[ "$cont_flag" -eq 1 ]]; then
  skip_burnin=1
fi

# --- Preprocess fasta inputs: decompress if gzipped ---
cds=$(decompress_if_gz "$cds")
TE_lib=$(decompress_if_gz "$TE_lib")
if [[ -n "$input_fasta" ]]; then
  input_fasta=$(decompress_if_gz "$input_fasta")
fi

# --- Process fixed insertion/excision numbers if provided ---
extra_fix_in=""
extra_fix_ex=""
if [[ -n "$fix" ]]; then
  IFS=',' read -r fix_in fix_ex <<< "$fix"
  extra_fix_in="--fix_in ${fix_in}"
  extra_fix_ex="--fix_ex ${fix_ex}"
fi

###############################################################################
# TE Library Processing (always executed)
###############################################################################
echo "=== TE Library Processing ===" | tee -a "$LOG"

# (A) Compute divergence info from the TE library.
cmd="python ${BIN_DIR}/seq_divergence.py -i ${TE_lib} -o lib.txt -t ${threads} --min_align 100 --max_off 20 --miu ${mutation_rate} --blast_outfmt '6 qseqid sseqid sstart send slen qstart qend qlen length nident btop'"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd >> "$LOG" 2>> "$ERR"
if [ $? -ne 0 ]; then
    echo "Error running seq_divergence.py" | tee -a "$ERR"
    exit 1
fi

# (B) Append LTR lengths to TE library.
cmd="python ${BIN_DIR}/LTR_fasta_header_appender.py -fasta ${TE_lib} -domains lib.txt -div_type none"
if [[ "$ex_ltr" -eq 1 ]]; then
    cmd+=" -exclude_no_hits"
fi
echo "Running: $cmd > lib.fa" | tee -a "$LOG"
eval $cmd > lib.fa 2>> "$ERR"
if [ $? -ne 0 ]; then
    echo "Error running LTR_fasta_header_appender.py" | tee -a "$ERR"
    exit 1
fi

# (C) Extract only intact TEs into a cleaned library
echo "=== Extracting intact TEs to lib_clean.fa ===" | tee -a "$LOG"
cmd="python ${BIN_DIR}/extract_intact_TEs.py --lib lib.fa --out_fasta lib_clean.fa"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd >> "$LOG" 2>> "$ERR"
if [ $? -ne 0 ]; then
    echo "Error running extract_intact_TEs.py" | tee -a "$ERR"
    exit 1
fi

# Now use lib_clean.fa for all downstream insertions
clean_lib="lib_clean.fa"

###############################################################################
# Phase 1: Burn-in Genome Generation (only if no external BED/FASTA provided)
###############################################################################
if [[ "$skip_burnin" -eq 0 ]]; then
    echo "=== Phase 1: Burn-in ===" | tee -a "$LOG"

    # (1a) Build synthetic genome with CDS using synthetic_genome.py.
    cmd="python ${BIN_DIR}/synthetic_genome.py -cds ${cds} -out_prefix backbone -chr_number ${chr_number} -size ${size} -seed ${seed}"
    if [[ -n "$cds_num" ]]; then
        cmd+=" -cds_num ${cds_num}"
    elif [[ -n "$cds_percent" ]]; then
        cmd+=" -cds_percent ${cds_percent}"
    fi
    echo "Running: $cmd" | tee -a "$LOG"
    eval $cmd >> "$LOG" 2>> "$ERR"
    if [ $? -ne 0 ]; then
        echo "Error running synthetic_genome.py" | tee -a "$ERR"
        exit 1
    fi

    # (1b) Insert TEs into the synthetic genome using the parallel version.
    cmd="python ${BIN_DIR}/shared_ltr_inserter_parallel.py -genome backbone.fa -TE ${clean_lib}"
    if [[ -n "$TE_percent" ]]; then
        cmd+=" -p ${TE_percent}"
    else
        cmd+=" -n ${TE_num}"
    fi
    cmd+=" -bed backbone.bed -output burnin -seed ${seed} -TE_ratio ${TE_ratio} -stat_out burnin.stat"
    cmd+=" -k ${TE_mut_k} -Mmax ${TE_mut_Mmax} -pdf_out burnin_mut_dist.pdf -m ${threads}"
    echo "Running: $cmd" | tee -a "$LOG"
    eval $cmd >> "$LOG" 2>> "$ERR"
    if [ $? -ne 0 ]; then
        echo "Error running shared_ltr_inserter_parallel.py" | tee -a "$ERR"
        exit 1
    fi
fi

###############################################################################
# Phase 2: Looping Generations
###############################################################################
echo "=== Phase 2: Looping Generations ===" | tee -a "$LOG"

# Calculate the total number of iterations.
iterations=$(( generation_end / step ))

# If --continue is set, look for the highest completed generation.
start_iter=1
if [[ "$cont_flag" -eq 1 ]]; then
  last_gen=0
  for file in gen*_final.fasta; do
    if [[ $file =~ gen([0-9]+)_final.fasta ]]; then
      num=${BASH_REMATCH[1]}
      if (( num > last_gen )); then
        last_gen=$num
      fi
    fi
  done
  if (( last_gen > 0 )); then
    start_iter=$(( last_gen / step + 1 ))
    echo "Resuming from generation $(( start_iter * step )) (iteration $start_iter) based on existing output." | tee -a "$LOG"
  else
    echo "No previous final outputs found; starting from the beginning." | tee -a "$LOG"
  fi
fi

# Generate a list of seeds for each iteration based on the provided base seed.
seed_list=($(python -c "import random; random.seed(${seed}); print(' '.join([str(random.randint(1,10000)) for _ in range(${iterations})]))"))
echo "Seed list for Phase 2 iterations: ${seed_list[@]}" | tee -a "$LOG"

# Initialize prev_lib for first generation
prev_lib="lib_clean.fa"

for (( i=start_iter; i<=iterations; i++ )); do
  # Calculate the true generation number for file naming.
  current_gen=$(( i * step ))
  current_seed=${seed_list[$((i-1))]}
  echo "----- Generation ${current_gen} using seed ${current_seed} -----" | tee -a "$LOG"
  
  # For the first iteration (if not resuming) use burn-in or user-supplied inputs.
  if [ $i -eq 1 ]; then
    if [ "$skip_burnin" -eq 1 ]; then
      if [[ -n "$input_fasta" && -n "$input_bed" ]]; then
        prev_genome="$input_fasta"
        prev_bed="$input_bed"
      else
        prev_genome="burnin.fa"
        prev_bed="burnin.bed"
      fi
    else
      prev_genome="burnin.fa"
      prev_bed="burnin.bed"
    fi
  else
    prev_genome="gen$(( (i-1) * step ))_final.fasta"
    prev_bed="gen$(( (i-1) * step ))_final.bed"
  fi

  # (2a) Mutate the genome.
  mut_prefix="gen${current_gen}_mut"
  cmd="${BIN_DIR}/ltr_mutator -fasta ${prev_genome} -rate ${mutation_rate} -generations ${step} -mode 0 -threads ${threads} -seed ${current_seed} -out_prefix ${mut_prefix}"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running ltr_mutator for generation ${current_gen}" | tee -a "$ERR"
    exit 1
  fi

  # (2b) Insert new TEs (allowing for nesting) using the parallel nest inserter.
  # Now using nest_inserter_parallel2.py and adding --disable_genes if specified.
  nest_prefix="gen${current_gen}_nest"
  cmd="python ${BIN_DIR}/nest_inserter_parallel2.py --genome ${mut_prefix}.fa --TE ${prev_lib} --generations ${step} --bed ${prev_bed} --output ${nest_prefix} --seed ${current_seed} --rate ${insert_rate} ${extra_fix_in} --TE_ratio ${TE_ratio} -bf burnin.stat --birth_rate ${birth_rate}"
  cmd+=" --euch_het_bias ${euch_bias_insert} --euch_het_buffer ${euch_buffer} -m ${threads}"
  if [[ $disable_genes -eq 1 ]]; then
    cmd+=" --disable_genes"
  fi
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running nest_inserter_parallel2.py for generation ${current_gen}" | tee -a "$ERR"
    exit 1
  fi

  # (2c) Purge some TEs and convert intact TEs to soloLTRs using TE excision.
  # Updated to use TE_exciser_parallel_dev2.py with new parameters.
  cmd="python ${BIN_DIR}/TE_exciser_parallel_dev2.py --genome ${nest_prefix}.fasta --bed ${nest_prefix}.bed --rate ${delete_rate} --generations ${step} --soloLTR_freq ${solo_rate} ${extra_fix_ex} --output gen${current_gen}_final --seed ${current_seed} --sigma ${sigma} --k ${k} --sel_coeff ${sel_coeff} -m ${threads} --euch_het_buffer ${euch_buffer} --euch_het_bias ${euch_bias_excise}"
  if [ $i -ne 1 ]; then
    cmd+=" --no_fig"
  fi
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running TE_exciser_parallel_dev2.py for generation ${current_gen}" | tee -a "$ERR"
    exit 1
  fi

  # If --keep_temps is not provided, remove intermediate files.
  if [[ "$keep_temps" -ne 1 ]]; then
    echo "Removing temporary files for generation ${current_gen}" | tee -a "$LOG"
    rm -f "${mut_prefix}.fa" "${nest_prefix}.bed" "${nest_prefix}.fasta"
  fi
  
  # (2d) Build the new per‑gen TE library
  echo "=== Extracting intact TEs for generation ${current_gen} into lib file ===" | tee -a "$LOG"
  cmd="python ${BIN_DIR}/extract_intact_TEs.py \
    --genome gen${current_gen}_final.fasta \
    --bed    gen${current_gen}_final.bed \
    --out_fasta gen${current_gen}_final.lib"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
      echo "Error running extract_intact_TEs.py for generation ${current_gen}" | tee -a "$ERR"
      exit 1
  fi

  # clean up the *previous* lib file unless the user wants to keep temps
  if [[ "$keep_temps" -ne 1 && "$prev_lib" != "lib_clean.fa" ]]; then
    rm -f "$prev_lib"
  fi

  # point to the newly built library for the next iteration
  prev_lib="gen${current_gen}_final.lib"

done

echo "Pipeline completed at $(date)" | tee -a "$LOG"

###############################################################################
# Supplemental Post-Processing (Global)
###############################################################################
echo "=== Global Post-Processing ===" | tee -a "$LOG"

# --- Modification (3): Choose the proper starting file names for post-processing ---
if [[ "$skip_burnin" -eq 1 ]]; then
    initial_bed="$input_bed"
    initial_fasta="$input_fasta"
else
    initial_bed="burnin.bed"
    initial_fasta="burnin.fa"
fi

# 1. Plot TE fraction (includes starting file and all gen*_final files)
cmd="python ${UTIL_DIR}/plot_TE_frac.py --bed \$(echo \"${initial_bed}\"; ls gen*_final.bed | sort -V) --fasta \$(echo \"${initial_fasta}\"; ls gen*_final.fasta | sort -V) --feature Intact_TE:SoloLTR:Fragmented_TE --out_prefix percent_TE"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd

# 2. Plot solo versus intact TE proportions.
cmd="python ${UTIL_DIR}/plot_solo_intact.py --bed \$(echo \"${initial_bed}\"; ls gen*_final.bed | sort -V) --out_prefix solo_intact"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd

# 3. Generate overall statistics report.
cmd="python ${UTIL_DIR}/stats_report.py --bed \$(ls gen*_final.bed | sort -V) --out_prefix stat"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd

# 4. Plot superfamily count.
cmd="python ${UTIL_DIR}/plot_superfamily_count.py"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd

# 5. Plot category bar.
cmd="python ${UTIL_DIR}/plot_category_bar.py"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd

# 6. Genome size through time.
cmd="python ${UTIL_DIR}/genome_plot.py"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd

###############################################################################
# Supplemental Post-Processing (Per-Generation Analysis)
###############################################################################
echo "=== Per-Generation Post-Processing ===" | tee -a "$LOG"

# Determine the total number of final generation files (iterations)
total_gens=$iterations

# Choose at most 4 evenly distributed generations.
selected=()
if [ "$total_gens" -le 4 ]; then
  for (( i=1; i<=total_gens; i++ )); do
    selected+=("$i")
  done
else
  # Always include the first and last iterations.
  selected+=(1)
  # Compute two intermediate indices.
  mid1=$(printf "%.0f" "$(echo "scale=2; 1 + ($total_gens - 1)/3" | bc -l)")
  mid2=$(printf "%.0f" "$(echo "scale=2; 1 + 2*($total_gens - 1)/3" | bc -l)")
  selected+=("$mid1" "$mid2" "$total_gens")
fi

# Sort the selected generations in descending order (highest to lowest)
IFS=$'\n' selected=($(sort -nr <<<"${selected[*]}"))
unset IFS

echo "Selected iterations for per-generation analysis: ${selected[@]}" | tee -a "$LOG"

for iter in "${selected[@]}"; do
  gen_num=$(( iter * step ))
  final_prefix="gen${gen_num}_final"
  echo "Processing per-generation analysis for ${final_prefix}" | tee -a "$LOG"
  
  # (a) Extract intact LTR sequences.
  cmd="python ${BIN_DIR}/extract_intact_LTR.py --bed ${final_prefix}.bed --genome ${final_prefix}.fasta --out_fasta ${final_prefix}_LTR.fasta"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd
  
  # (b) Compute sequence divergence on extracted LTR sequences.
  cmd="python ${BIN_DIR}/seq_divergence.py -i ${final_prefix}_LTR.fasta -o ${final_prefix}_LTR.tsv -t ${threads} --min_align 100 --max_off 20 --miu ${mutation_rate} --blast_outfmt '6 qseqid sseqid sstart send slen qstart qend qlen length nident btop'"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd
done

# Run ltr_dens.py once after per-generation analyses.
cmd="python ${BIN_DIR}/ltr_dens.py --model ${model} --output all_LTR_density.pdf --miu ${mutation_rate}"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd

echo "Post-processing completed at $(date)" | tee -a "$LOG"
