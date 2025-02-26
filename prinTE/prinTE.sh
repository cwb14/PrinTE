#!/bin/bash
###############################################################################
# prinTE.sh
#
# A wrapper to run the TE evolution simulation pipeline in two phases.
#
# Phase 1 (Burn‐in) is performed unless the user provides starting input files
# via --bed and --fasta, or if --continue is provided and previous outputs are found.
#
#   (1) Burn‐in (skipped if both --bed and --fasta are provided OR if --continue is provided and previous outputs exist):
#       (a) Build a synthetic genome with CDS using synthetic_genome.py.
#           Note: synthetic_genome.py accepts two additional, mutually
#           exclusive optional arguments: -cds_num and -cds_percent.
#       (b) Estimate LTR lengths with seq_divergence.py.
#       (c) Append LTR lengths to the TE library.
#       (d) Insert TEs using shared_ltr_inserter.py.
#           Note: shared_ltr_inserter.py now accepts either -n (number of
#           insertions) or -p (target percent TE) (but not both) and an extra
#           option: -stat_out burnin.stat.
#
#   (2) Looping generations:
#       For each generation “step” (until generation_end) run:
#         - ltr_mutator (simulate mutation),
#         - nest_inserter (insert new TEs, with an optional fixed insertion
#           count via --fix_in). Updated nest_inserter.py now also accepts:
#             --TE_ratio, -bf burnin.stat, --birth_rate (default 1e-3),
#             and two new parameters controlling euchromatin bias:
#             --euch_het_bias and --euch_het_buffer.
#           In this wrapper, these are provided with friendlier names:
#             --euch-bias (default: 1.2) and --euch-buffer (default: 10000).
#         - TE_exciser (purge some TEs and convert a percentage of intact TEs to
#           soloLTRs). New options: --sigma (default 1.0) and --k (default 10).
#           The first TE_exciser run will include figures; subsequent runs add
#           --no_fig.
#
# New functionality and changes:
#
#  (1) New flag --continue:
#      When provided, the script searches for existing final generation outputs
#      (e.g., gen${i}_final.fasta) and resumes from the next iteration.
#
#  (2) New flag --keep_temps (or -kt):
#      When provided, temporary files (e.g., ${mut_prefix}.fa, ${nest_prefix}.bed, and ${nest_prefix}.fasta)
#      are kept. Otherwise, these are removed at the end of each loop.
#
#  (3) New user-friendly parameters --euch-bias and --euch-buffer
#      (mapped internally to --euch_het_bias and --euch_het_buffer)
#      to control the insertion bias toward open euchromatin.
#
# The script determines the directories for the pipeline’s binaries
# (assumed to be in TOOL_DIR/bin/) and tools (TOOL_DIR) and writes its log and
# error files ("pipeline.log" and "pipeline.error") in the current directory.
#
# Usage:
#   prinTEsh [options]
#
# Options (defaults in parentheses; some options are REQUIRED):
#
#   -c,  --cds              Path to CDS file for synthetic_genome.py
#                           (default: $TOOL_DIR/TAIR10.cds.fa)
#
#   -N,  --cds_num          Number of CDS sequences to insert.
#                           (Mutually exclusive with --cds_percent)
#
#   -P,  --cds_percent      Percent of the genome that should be CDS.
#                           (Mutually exclusive with --cds_num)
#
#   -cn, --chr_number       Number of chromosomes (default: 3)
#
#   -sz, --size             Genome size in kb, Mb, or Gb (default: 100Mb)
#
#   -s,  --seed             Random seed (default: 42)
#
#   -i,  --TE_lib           TE library file (default: $TOOL_DIR/combined_curated_TE_lib_ATOSZM_selected.fasta)
#
#   -m,  --mutation_rate    Mutation rate (default: 1.3e-8)
#
#   -r,  --TE_ratio         TE ratio file (default: $TOOL_DIR/ratios.tsv)
#
#   -n,  --TE_num           Number of TE insertions in burn‐in (default: 2000)
#
#   -p,  --TE_percent       Target percent of the genome that should be TE
#                           (mutually exclusive with --TE_num)
#
#   -st, --step             Generation step (number of generations per loop; REQUIRED)
#
#   -ge, --generation_end   Final generation (REQUIRED; must be a multiple of step)
#
#   -t,  --threads          Number of threads (default: 4)
#
#   -ir, --insert_rate      TE insertion rate for nest_inserter.py (default: 1e-8)
#
#   -br, --birth_rate       Birth rate for new TEs in nest_inserter.py (default: 1e-3)
#
#   -dr, --delete_rate      TE deletion rate for TE_exciser.py (default: 1e-4)
#
#   -sr, --solo_rate        Percentage chance to convert an intact TE to soloLTR
#                           in TE_exciser.py (default: 5)
#
#   -sc, --sigma           Sigma value for TE_exciser.py (default: 1.0)
#
#   -k,  --k               TE length decay slope for TE_exciser.py (default: 10)
#
#   -F,  --fix              Fixed insertion and excision numbers, comma-separated.
#                           Format: insertion,excision (e.g., 500,200)
#
#   -b,  --bed              Input BED file for starting generation (requires --fasta)
#
#   -f,  --fasta            Input FASTA file for starting generation (requires --bed)
#
#   -x, --continue          Resume simulation from the last completed generation.
#
#   -z, --keep_temps, -kt   Keep temporary files (default: remove them after each loop).
#
#   -w, --euch-bias         Weight of bias toward euchromatin (default: 1.2)
#
#   -j, --euch-buffer       Interval upstream/downstream of genes considered euchromatin (default: 10000)
#
#   -h,  --help             Display this help message and exit.
#
# Note: Fasta inputs (CDS, TE_lib, and -fasta) may be provided as plain text or gzipped (.gz).
#
# Example:
#   $(basename "$0") -st 1000 -ge 3000
#
###############################################################################

# --- Determine directories ---
# Assume this script is stored in TOOL_DIR; BIN_DIR is assumed to be TOOL_DIR/bin
TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
BIN_DIR="$TOOL_DIR/bin"

# --- Log file names ---
LOG="pipeline.log"
ERR="pipeline.error"
echo "Pipeline started at $(date)" > "$LOG"
echo "Pipeline started at $(date)" > "$ERR"

# --- Temporary files array (for unzipped inputs) ---
temp_files=()
cleanup() {
  for f in "${temp_files[@]}"; do
    [ -f "$f" ] && rm -f "$f"
  done
}
trap cleanup EXIT

# --- Function: If a fasta file is gzipped, decompress it to a temporary file.
#     If the given file does not exist and does not end with .gz, check if a .gz version exists.
decompress_if_gz() {
  local file="$1"
  if [[ ! -f "$file" && "$file" != *.gz ]]; then
    if [ -f "${file}.gz" ]; then
      file="${file}.gz"
    fi
  fi
  if [[ "$file" == *.gz ]]; then
    local tmp
    tmp=$(mktemp)
    if gunzip -c "$file" > "$tmp"; then
      temp_files+=("$tmp")
      echo "$tmp"
    else
      echo "Error decompressing $file" | tee -a "$ERR"
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
  -c,  --cds              Path to CDS file for synthetic_genome.py
                         (default: ${TOOL_DIR}/TAIR10.cds.fa)
                         (gzipped files ending in .gz are supported)
  -N,  --cds_num          Number of CDS sequences to insert.
                         (Mutually exclusive with --cds_percent)
  -P,  --cds_percent      Percent of the genome that should be CDS.
                         (Mutually exclusive with --cds_num)
  -cn, --chr_number       Number of chromosomes (default: 3)
  -sz, --size             Genome size in kb, Mb, or Gb (default: 100Mb)
  -s,  --seed             Random seed (default: 42)
  -i,  --TE_lib           TE library file (default: ${TOOL_DIR}/combined_curated_TE_lib_ATOSZM_selected.fasta)
                         (gzipped files ending in .gz are supported)
  -m,  --mutation_rate    Mutation rate (default: 1.3e-8)
  -r,  --TE_ratio         TE ratio file (default: ${TOOL_DIR}/ratios.tsv)
  -n,  --TE_num           Number of TE insertions in burn‐in (default: 2000)
  -p,  --TE_percent       Target percent of the genome that should be TE
                         (mutually exclusive with --TE_num)
  -st, --step             Generation step (number of generations per loop; REQUIRED)
  -ge, --generation_end   Final generation (REQUIRED; must be a multiple of step)
  -t,  --threads          Number of threads (default: 4)
  -ir, --insert_rate      TE insertion rate for nest_inserter.py (default: 1e-8)
  -br, --birth_rate       Birth rate for new TEs in nest_inserter.py (default: 1e-3)
  -dr, --delete_rate      TE deletion rate for TE_exciser.py (default: 1e-4)
  -sr, --solo_rate        Percentage chance to convert an intact TE to soloLTR
                         in TE_exciser.py (default: 5)
  -sc, --sigma           Sigma value for TE_exciser.py (default: 1.0)
  -k,  --k               TE length decay slope for TE_exciser.py (default: 10)
  -F,  --fix              Fixed insertion and excision numbers, comma-separated.
                         Format: insertion,excision (e.g., 500,200)
  -b,  --bed              Input BED file for starting generation (requires --fasta)
  -f,  --fasta            Input FASTA file for starting generation (requires --bed)
                         (gzipped files ending in .gz are supported)
  -x, --continue          Resume simulation from the last completed generation.
  --keep_temps, -kt       Keep temporary files (default: remove them after each loop)
  -w, --euch-bias         Weight of bias toward euchromatin (default: 1.2)
  -j, --euch-buffer       Interval upstream/downstream of genes considered euchromatin (default: 10000)
  -h,  --help             Display this help message and exit

Example:
  $(basename "$0") -st 1000 -ge 3000
EOF
}

# --- Parse command-line options ---
# Initialize flags and new parameters with defaults.
cont_flag=0
keep_temps=0
euch_bias=1.2
euch_buffer=10000

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
    -sc|--sigma)
      sigma="$2"
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
    --keep_temps| -kt)
      keep_temps=1
      shift;;
    -w|--euch-bias)
      euch_bias="$2"
      shift; shift;;
    -j|--euch-buffer)
      euch_buffer="$2"
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
chr_number="${chr_number:-3}"
size="${size:-100Mb}"
seed="${seed:-42}"
TE_lib="${TE_lib:-${TOOL_DIR}/combined_curated_TE_lib_ATOSZM_selected.fasta}"
mutation_rate="${mutation_rate:-1.3e-8}"
TE_ratio="${TE_ratio:-${TOOL_DIR}/ratios.tsv}"
TE_num="${TE_num:-2000}"
threads="${threads:-4}"
insert_rate="${insert_rate:-1e-8}"
birth_rate="${birth_rate:-1e-3}"
delete_rate="${delete_rate:-1e-4}"
solo_rate="${solo_rate:-5}"
sigma="${sigma:-1.0}"
k="${k:-10}"

# --- Validate mutually exclusive CDS options ---
if [[ -n "$cds_num" && -n "$cds_percent" ]]; then
  echo "Error: Please provide either --cds_num (or -N) or --cds_percent (or -P), not both." | tee -a "$ERR"
  exit 1
fi

# --- Validate mutually exclusive TE insertion options ---
if [[ -z "$TE_percent" ]]; then
  TE_num="${TE_num:-2000}"
fi

# --- Validate required parameters ---
if [[ -z "$step" || -z "$generation_end" ]]; then
  echo "Error: Both --step and --generation_end must be provided." | tee -a "$ERR"
  print_help
  exit 1
fi

# Ensure generation_end is an exact multiple of step.
if (( generation_end % step != 0 )); then
  echo "Error: generation_end must be an exact multiple of step." | tee -a "$ERR"
  exit 1
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
# Phase 1: Burn-in (skip if BED/FASTA are provided or if resuming with --continue)
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

  # (1b) Get LTR lengths from the TE fasta library.
  cmd="python ${BIN_DIR}/seq_divergence.py -i ${TE_lib} -o lib.txt -t ${threads} --min_align 100 --max_off 20 --miu ${mutation_rate} --blast_outfmt '6 qseqid sseqid sstart send slen qstart qend qlen length nident btop'"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running seq_divergence.py" | tee -a "$ERR"
    exit 1
  fi

  # (1c) Append LTR lengths to the TE library.
  cmd="python ${BIN_DIR}/LTR_fasta_header_appender.py -fasta ${TE_lib} -domains lib.txt -div_type none"
  echo "Running: $cmd > lib.fa" | tee -a "$LOG"
  eval $cmd > lib.fa 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running LTR_fasta_header_appender.py" | tee -a "$ERR"
    exit 1
  fi

  # (1d) Insert TEs into the synthetic genome.
  # Use either -n (TE_num) or -p (TE_percent) and add -stat_out burnin.stat.
  cmd="python ${BIN_DIR}/shared_ltr_inserter.py -genome backbone.fa -TE lib.fa"
  if [[ -n "$TE_percent" ]]; then
    cmd+=" -p ${TE_percent}"
  else
    cmd+=" -n ${TE_num}"
  fi
  cmd+=" -bed backbone.bed -output burnin -seed ${seed} -TE_ratio ${TE_ratio} -stat_out burnin.stat"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running shared_ltr_inserter.py" | tee -a "$ERR"
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
  last_iter=0
  for file in gen*_final.fasta; do
    if [[ $file =~ gen([0-9]+)_final.fasta ]]; then
      num=${BASH_REMATCH[1]}
      if (( num > last_iter )); then
        last_iter=$num
      fi
    fi
  done
  if (( last_iter > 0 )); then
    start_iter=$(( last_iter + 1 ))
    echo "Resuming from generation $(( start_iter * step )) (iteration $start_iter) based on existing output." | tee -a "$LOG"
  else
    echo "No previous final outputs found; starting from the beginning." | tee -a "$LOG"
  fi
fi

# Generate a list of seeds for each iteration based on the provided base seed.
seed_list=($(python -c "import random; random.seed(${seed}); print(' '.join([str(random.randint(1,10000)) for _ in range(${iterations})]))"))
echo "Seed list for Phase 2 iterations: ${seed_list[@]}" | tee -a "$LOG"

for (( i=start_iter; i<=iterations; i++ )); do
  current_seed=${seed_list[$((i-1))]}
  echo "----- Generation $(( i * step )) using seed ${current_seed} -----" | tee -a "$LOG"
  
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
    prev_genome="gen$((i-1))_final.fasta"
    prev_bed="gen$((i-1))_final.bed"
  fi

  # (2a) Mutate the genome.
  mut_prefix="gen${i}_mut"
  cmd="${BIN_DIR}/ltr_mutator -fasta ${prev_genome} -rate ${mutation_rate} -generations ${step} -mode 0 -threads ${threads} -seed ${current_seed} -out_prefix ${mut_prefix}"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running ltr_mutator for generation $i" | tee -a "$ERR"
    exit 1
  fi

  # (2b) Insert new TEs (allowing for nesting) using nest_inserter.py.
  nest_prefix="gen${i}_nest"
  cmd="python ${BIN_DIR}/nest_inserter.py --genome ${mut_prefix}.fa --TE lib.fa --generations ${step} --bed ${prev_bed} --output ${nest_prefix} --seed ${current_seed} --rate ${insert_rate} ${extra_fix_in} --TE_ratio ${TE_ratio} -bf burnin.stat --birth_rate ${birth_rate}"
  # Append the new euchromatin bias parameters (mapping our friendlier names).
  cmd+=" --euch_het_bias ${euch_bias} --euch_het_buffer ${euch_buffer}"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running nest_inserter.py for generation $i" | tee -a "$ERR"
    exit 1
  fi

  # (2c) Purge some TE insertions and convert intact TEs to soloLTRs using TE_exciser.py.
  cmd="python ${BIN_DIR}/TE_exciser.py --genome ${nest_prefix}.fasta --bed ${nest_prefix}.bed --rate ${delete_rate} --generations ${step} --soloLTR_freq ${solo_rate} ${extra_fix_ex} --output gen${i}_final --seed ${current_seed} --sigma ${sigma} --k ${k}"
  # For generation 1, include figures; for later generations, add --no_fig.
  if [ $i -ne 1 ]; then
    cmd+=" --no_fig"
  fi
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running TE_exciser.py for generation $i" | tee -a "$ERR"
    exit 1
  fi

  # If --keep_temps is not provided, remove intermediate files.
  if [[ "$keep_temps" -ne 1 ]]; then
    echo "Removing temporary files for iteration $i" | tee -a "$LOG"
    rm -f "${mut_prefix}.fa" "${nest_prefix}.bed" "${nest_prefix}.fasta"
  fi

done

echo "Pipeline completed at $(date)" | tee -a "$LOG"
