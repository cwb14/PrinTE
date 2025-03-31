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
#  (3) New user-friendly parameters --euch-bias and --euch-buffer to control insertion bias.
#
#  (4) New option --model for per-generation post-processing (model options: raw, K2P, JC69; default: raw).
#
#  (5) New option -m, --TE_mut_in is replaced by two new options for shared_ltr_inserter_parallel.py:
#         --TE_mut_k and --TE_mut_Mmax, which replace the old TE_mut_in parameter.
#         Also, shared_ltr_inserter_parallel.py now takes:
#           - '-pdf_out burnin_mut_dist.pdf'
#           - '-m ${threads}' for threads.
#
#  (6) New parallel versions of internal scripts:
#         - Use 'shared_ltr_inserter_parallel.py' instead of 'shared_ltr_inserter.py'
#         - Use 'nest_inserter_parallel.py' instead of 'nest_inserter.py', which now adds a '-m' parameter for threads.
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
############### REQUIRED ###############
  -ge, --generation_end   Number of generations to simulate (REQUIRED; must be an exact multiple of step)
  -st, --step             Generation step size (number of generations per loop; REQUIRED)

############### BURN-IN ###############
  -c,  --cds              Path to CDS file (default: ${TOOL_DIR}/TAIR10.cds.fa)
  -N,  --cds_num          Number of CDS sequences to insert. (Mutually exclusive with --cds_percent)
  -P,  --cds_percent      Percent of the genome that should be CDS. (Mutually exclusive with --cds_num)
  -n,  --TE_num           Number of TE insertions in burn‐in (default: 2000) (Mutually exclusive with --TE_percent)
  -p,  --TE_percent       Percent of the genome that should be TEs in the burn‐in (Mutually exclusive with --TE_num)
  -cn, --chr_number       Number of chromosomes (default: 3)
  -sz, --size             Genome size in kb, Mb, or Gb (default: 100Mb)
  # The old option -tm, --TE_mut_in is replaced by two new options below:
  -tk, --TE_mut_k             Slope of exponential decay for TE mutation (default: 10)
  -tmx, --TE_mut_Mmax          X-limit for exponential decay function (default: 20)

########## FIXED TE INDEL RATE ##########
  -F,  --fix              Fixed insertion and deletion numbers, comma-separated. Format: insertion,deletion (e.g., 1e-9,1e-9)

######### VARIABLE TE INDEL RATE #########  
  -ir, --insert_rate      TE insertion rate (default: 1e-8)
  -dr, --delete_rate      TE deletion rate (default: 1e-4)  
  -br, --birth_rate       TE birth rate (default: 1e-3)
  -sc, --sigma            Selection coefficient for gene insertions (default: 1.0)                
  
########## GENERAL USE ##########                 
  -s,  --seed             Random seed (default: 42)
  -i,  --TE_lib           TE library file (default: ${TOOL_DIR}/combined_curated_TE_lib_ATOSZM_selected.fasta)
  -m,  --mutation_rate    DNA Mutation rate (default: 1.3e-8)
  -r,  --TE_ratio         TE ratio file (default: ${TOOL_DIR}/ratios.tsv)
  -t,  --threads          Number of threads (default: 4)
  -sr, --solo_rate        Percentage chance to convert an intact TE to soloLTR in TE_exciser.py (default: 5)
  -k,  --k                TE length decay slope for TE excision (default: 10)
  -b,  --bed              Input BED file for starting generation (skip burn-in; requires --fasta)
  -f,  --fasta            Input FASTA file for starting generation (skip burn-in; requires --bed)
  -x, --continue          Resume simulation from the last completed generation.
  --keep_temps, -kt       Keep temporary files (default: remove them after each loop)
  -w, --euch-bias         Weight of bias toward euchromatin (default: 1.0; disabled)
  -j, --euch-buffer       Interval upstream/downstream of genes considered euchromatin (default: 10000)
  --model                 DNA mutation model for LTR dating (choose from: raw, K2P, JC69; default: raw)
  -h,  --help             Display this help message and exit

Example:
  $(basename "$0") -st 1000 -ge 3000
EOF
}

# --- Parse command-line options ---
# Initialize flags and new parameters with defaults.
cont_flag=0
keep_temps=0
euch_bias=1.0
euch_buffer=10000
model="raw"
# New TE inserter mutation options defaults:
TE_mut_k=10
TE_mut_Mmax=20

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
    -w|--euch-bias)
      euch_bias="$2"
      shift; shift;;
    -j|--euch-buffer)
      euch_buffer="$2"
      shift; shift;;
    --model)
      model="$2"
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

  # (1b) Estimate LTR lengths from the TE fasta library.
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

  # (1d) Insert TEs into the synthetic genome using the parallel version.
  cmd="python ${BIN_DIR}/shared_ltr_inserter_parallel.py -genome backbone.fa -TE lib.fa"
  if [[ -n "$TE_percent" ]]; then
    cmd+=" -p ${TE_percent}"
  else
    cmd+=" -n ${TE_num}"
  fi
  cmd+=" -bed backbone.bed -output burnin -seed ${seed} -TE_ratio ${TE_ratio} -stat_out burnin.stat"
  # Add new TE mutation parameters in place of the old TE_mut_in.
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

  # (2b) Insert new TEs (allowing for nesting) using the parallel nest inserter.
  nest_prefix="gen${i}_nest"
  cmd="python ${BIN_DIR}/nest_inserter_parallel.py --genome ${mut_prefix}.fa --TE lib.fa --generations ${step} --bed ${prev_bed} --output ${nest_prefix} --seed ${current_seed} --rate ${insert_rate} ${extra_fix_in} --TE_ratio ${TE_ratio} -bf burnin.stat --birth_rate ${birth_rate}"
  cmd+=" --euch_het_bias ${euch_bias} --euch_het_buffer ${euch_buffer} -m ${threads}"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running nest_inserter_parallel.py for generation $i" | tee -a "$ERR"
    exit 1
  fi

  # (2c) Purge some TEs and convert intact TEs to soloLTRs using TE_exciser.py.
  cmd="python ${BIN_DIR}/TE_exciser.py --genome ${nest_prefix}.fasta --bed ${nest_prefix}.bed --rate ${delete_rate} --generations ${step} --soloLTR_freq ${solo_rate} ${extra_fix_ex} --output gen${i}_final --seed ${current_seed} --sigma ${sigma} --k ${k}"
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

###############################################################################
# Supplemental Post-Processing (Global)
###############################################################################
echo "=== Global Post-Processing ===" | tee -a "$LOG"

# Run the global analysis scripts from the util directory.
# 1. Plot TE fraction (includes burnin and all gen*_final files)
python ${UTIL_DIR}/plot_TE_frac.py --bed $(echo "burnin.bed"; ls gen*_final.bed | sort -V) \
  --fasta $(echo "burnin.fa"; ls gen*_final.fasta | sort -V) \
  --feature Intact_TE:SoloLTR:Fragmented_TE --out_prefix percent_TE

# 2. Plot solo versus intact TE proportions.
python ${UTIL_DIR}/plot_solo_intact.py --bed $(echo "burnin.bed"; ls gen*_final.bed | sort -V) \
  --out_prefix solo_intact

# 3. Generate overall statistics report.
python ${UTIL_DIR}/stats_report.py --bed $(ls gen*_final.bed | sort -V) \
  --out_prefix stat

# 4. Plot superfamily count.
python ${UTIL_DIR}/plot_superfamily_count.py

# 5. Plot category bar.
python ${UTIL_DIR}/plot_category_bar.py

# 6. Genome size through time.
python ${UTIL_DIR}/genome_plot.py

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

# Sort the selected generations in descending order (as per example: highest to lowest)
IFS=$'\n' selected=($(sort -nr <<<"${selected[*]}"))
unset IFS

echo "Selected generations for per-generation analysis: ${selected[@]}" | tee -a "$LOG"

for i in "${selected[@]}"; do
  final_prefix="gen${i}_final"
  echo "Processing per-generation analysis for ${final_prefix}" | tee -a "$LOG"
  
  # (a) Extract intact LTR sequences.
  python ${BIN_DIR}/extract_intact_LTR.py --bed ${final_prefix}.bed --genome ${final_prefix}.fasta --out_fasta ${final_prefix}_LTR.fasta
  
  # (b) Compute sequence divergence on extracted LTR sequences.
  python ${BIN_DIR}/seq_divergence.py -i ${final_prefix}_LTR.fasta -o ${final_prefix}_LTR.tsv -t ${threads} \
    --min_align 100 --max_off 20 --miu ${mutation_rate} --blast_outfmt '6 qseqid sseqid sstart send slen qstart qend qlen length nident btop'
  
  # Removed ltr_dens.py call from here.
done

# Run ltr_dens.py once after per-generation analyses.
echo "Running global LTR density analysis" | tee -a "$LOG"
python ${BIN_DIR}/ltr_dens.py --model ${model} --output all_LTR_density.pdf --miu ${mutation_rate}

echo "Post-processing completed at $(date)" | tee -a "$LOG"
