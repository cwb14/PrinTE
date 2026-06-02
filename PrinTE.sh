#!/bin/bash
###############################################################################
# PrinTE.sh
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
#         - Use 'shared_ltr_inserter_parallel2.py' (updated CLI: -n_intact/-p_intact plus -n_frag/-p_frag)
#         - Use 'nest_inserter_parallel.py' with '-m' for threads and optional --disable_genes (-dg).
#
#  (8) New flag -bo, --burnin_only:
#      When activated, the script will run the burn-in phase only and then exit.
#      In this case, --generation_end and --step are not required.
#
#  (9) File naming: The output files from Phase 2 now have names
#      reflecting the true generation simulated.
#
#  (10) Updated TE excision: Use "TE_exciser_parallel.py" (instead of TE_exciser2.py)
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

# --- Auto-detect OS and pick appropriate ltr_mutator binary ---
OS="$(uname)"
case "$OS" in
  Darwin)
    mutator_exec="ltr_mutator_mac"
    ;;
  Linux)
    mutator_exec="ltr_mutator"
    ;;
  *)
    echo "Error: Unsupported OS '$OS'. Only macOS (Darwin) or Linux are supported." >&2
    exit 1
    ;;
esac

# --- Verify mutator binary works; recompile from source if not ---
# Called once after LOG/ERR are initialised (see Change 2).
ensure_mutator() {
  local bin="${BIN_DIR}/${mutator_exec}"
  local src="${BIN_DIR}/ltr_mutator.cpp"

  # Try to execute the binary. Any exit that isn't a clean run
  # (loader errors, GLIBC mismatches, wrong ELF class, etc.) will
  # produce a non-zero exit before main() even starts.
  local probe_out probe_exit
  probe_out=$( "$bin" </dev/null 2>&1 )
  probe_exit=$?

  # Exit 0 or a "usage / help" style non-zero (e.g. 1) where the binary
  # actually printed something to stdout/stderr means it loaded fine.
  # We distinguish a loader failure by checking for known error phrases.
  if echo "$probe_out" | grep -qiE \
       'not found|cannot execute|exec format|no such file|illegal instruction|GLIBC|version.*required'; then
    :  # fall through to recompile
  elif [[ $probe_exit -eq 126 || $probe_exit -eq 127 ]]; then
    :  # permission denied / command not found → fall through
  else
    echo "Mutator binary OK: ${bin}" | tee -a "$LOG"
    return 0
  fi

  echo "WARNING: Pre-compiled binary '${bin}' does not run on this system." | tee -a "$LOG"
  echo "  Reason hint: $(echo "$probe_out" | head -3)" | tee -a "$LOG"
  echo "Attempting to recompile from source: ${src}" | tee -a "$LOG"

  if [[ ! -f "$src" ]]; then
    echo "Error: Source file '${src}' not found. Cannot recompile ltr_mutator." | tee -a "$ERR"
    exit 1
  fi

  if ! command -v g++ &>/dev/null; then
    echo "Error: 'g++' is not available on PATH. Cannot recompile ltr_mutator." | tee -a "$ERR"
    exit 1
  fi

  local recompile_cmd="g++ -std=c++17 -O3 -fopenmp \"${src}\" -o \"${bin}\""
  echo "Running: ${recompile_cmd}" | tee -a "$LOG"
  if eval "$recompile_cmd" >> "$LOG" 2>> "$ERR"; then
    echo "Recompile succeeded. Pipeline will use the freshly built binary." | tee -a "$LOG"
  else
    echo "Error: Recompile of ltr_mutator failed. See ${ERR} for details." | tee -a "$ERR"
    exit 1
  fi
}

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
print_usage() {
  cat <<EOF
PrinTE - a forward, per-generation simulator of transposable-element (TE) genome evolution.

Usage:
  $(basename "$0") [options]

Common workflows:
  # Build a burn-in genome only (synthetic genome + TEs), then stop:
  $(basename "$0") --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10

  # Forward-evolve for 30,000 generations, sampling every 10,000:
  $(basename "$0") -sz 135Mb -cn 5 -P 25 -itp 21 -ir 1e-6 -dr 4e-6 -ge 30000 -st 10000 -t 20

  # Start from your own genome instead of a burn-in:
  $(basename "$0") -f genome.fasta -b genome.bed -ge 20000 -st 10000

  # Resume / extend a finished run (optionally with new rates):
  $(basename "$0") -ge 50000 -st 10000 -F 7e-11,1e-11 --continue

Most-used options:
  -sz -cn -P -itp     burn-in genome size / chromosomes / % CDS / % intact-TE
  -ge -st             total generations / generations per step   (required to evolve)
  -ir -dr   |   -F    variable per-TE rates   OR   fixed per-bp rates (insertion,deletion)
  -m -TsTv -sr -k     DNA mutation rate / Ts:Tv / solo-LTR % / length-bias
  -f -b               start from an existing FASTA + BED (skip burn-in)
  --burnin_only       stop after the burn-in genome
  --continue          resume from the last completed generation
  -t                  threads

Run '$(basename "$0") -h' for the full option list.
EOF
}

print_help() {
  cat <<EOF
PrinTE - a forward, per-generation simulator of transposable-element (TE) genome evolution.
It builds (or imports) a starting genome, then repeatedly mutates the DNA, inserts new TEs,
and excises old ones, writing a genome FASTA + annotation BED + TE library per sampled generation.

Usage:
  $(basename "$0") [options]

Pipeline:  TE-library prep -> [Phase 1: burn-in genome] -> [Phase 2: generation loop] -> post-processing.
Phase 1 is skipped when you supply --fasta/--bed or --continue.

----------------------- EVOLUTION (Phase 2; required unless --burnin_only) -----------------------
  -ge,  --generation_end N     Total generations to simulate (must be a multiple of --step).
  -st,  --step N               Generations per loop. e.g. -ge 30000 -st 10000 samples 10k,20k,30k.

-------------------------------- BURN-IN GENOME (Phase 1) ----------------------------------------
  -sz,  --size SIZE            Genome size with a kb/Mb/Gb suffix (default: 400Mb).
  -cn,  --chr_number N         Number of chromosomes (default: 4). Work is parallel per chromosome.
  -c,   --cds FILE             CDS FASTA used as 'genes' (default: data/TAIR10.cds.fa, 19,621 seqs).
  -P,   --cds_percent FLOAT    Target % of genome that is CDS        (mutually exclusive with -N).
  -N,   --cds_num N            Number of CDS to insert               (mutually exclusive with -P).
  -itp, --intact_TE_percent F  Target % genome as INTACT TE bp       (default: 20; excl. -itn).
  -itn, --intact_TE_num N      Number of INTACT TE insertions        (excl. -itp).
  -ftp, --frag_TE_percent F    Target % genome as FRAGMENTED TE bp   (default: 0; excl. -ftn).
  -ftn, --frag_TE_num N        Number of FRAGMENTED TE insertions    (excl. -ftp).
  -i,   --TE_lib FILE          TE library FASTA, RepeatMasker headers (default: data/maize_rice_arab_curated_TE.lib.gz).
  -cl,  --clean_lib FILE       Pre-cleaned TE library; skips library processing (alternative to -i).
  -r,   --TE_ratio FILE        Per-superfamily insertion weights (default: ratios.tsv).
  -tk,  --TE_mut_k FLOAT       Burn-in TE age-decay slope (default: 10; larger -> TEs look younger).
  -tmx, --TE_mut_Mmax FLOAT    Max burn-in TE divergence % (default: 20; use 0 to disable TE aging).
  -mb,  --mutation_bins FILE   Explicit TE age distribution (3-col bins); overrides -tk/-tmx.
  -bo,  --burnin_only          Build the burn-in genome and exit (no Phase 2; -ge/-st not needed).

------------------------------ TE RATES (choose Variable OR Fixed) -------------------------------
  Variable -- rates scale with the evolving TE content; insertion mutations accrue over time:
  -ir,  --insert_rate RATE     Insertions per intact TE per generation (default: 1e-8).
  -dr,  --delete_rate RATE     Deletions  per intact TE per generation (default: 1e-7).
  -br,  --birth_rate RATE      Rate of reintroducing TEs from the original library (default: 1e-3).
  Fixed -- constant per-base rates:
  -F,   --fix INS,DEL          Fixed insertion,deletion per bp per generation, e.g. -F 5e-9,1e-8.
  -dg,  --disable_genes        Forbid TE insertion into genes (Fixed mode only).

------------------ SELECTION, CHROMATIN & SOLO-LTR DYNAMICS (mostly Variable mode) ---------------
  -sc,  --sigma FLOAT          Spread of per-gene selective constraint (log-normal; default: 1.0).
  -sf,  --sel_coeff FLOAT      Bias deletion toward more deleterious insertions (default: 0 = off).
  -cbi, --chromatin_bias_insert F  Euchromatin (near-gene) insertion bias (default: 1.0).
  -cbd, --chromatin_bias_delete F  Euchromatin (near-gene) deletion  bias (default: 1.0).
  -cb,  --chromatin_buffer N       bp around genes treated as euchromatin (default: 10000).
  -pb,  --promoter-boundary N      bp up/downstream of genes where a TE overlap disrupts a gene (default: 0).
  -sr,  --solo_rate FLOAT      % chance an excised intact LTR-RT leaves a solo-LTR (default: 95).
  -k,   --k FLOAT              Length-bias slope for deletion; longer TEs lost faster (default: 10; 0 = off).

----------------------------------- INPUT / RESUME ----------------------------------------------
  -f,   --fasta FILE           Start-genome FASTA (skip burn-in; requires --bed).
  -b,   --bed FILE             Start-annotation BED (skip burn-in; requires --fasta).
  -x,   --continue             Resume from the last completed generation in the working directory.

----------------------------------- DNA MUTATION ------------------------------------------------
  -m,   --mutation_rate RATE   Substitution rate per bp per generation (default: 1.3e-8).
  -TsTv,--TsTv FLOAT           Transition/transversion ratio (default: 1.0).

------------------------- GENOME-SIZE BOUNDS (optional early stop, Phase 2) ----------------------
  -mxgs,--max_size SIZE        Stop early once the projected size is confidently ABOVE SIZE (e.g. 1G).
  -mngs,--min_size SIZE        Stop early once the projected size is confidently BELOW SIZE (e.g. 100M).

------------------------------- OUTPUT / POST-PROCESSING -----------------------------------------
  -md,  --model {raw,K2P,JC69} Substitution model for LTR-RT dating (default: K2P).
  -pgs, --pergen_select N      How many generations get per-gen LTR-RT dating (default: 2 = first + last).
  -ex,  --ex_LTR               Drop library LTR-RTs that lack a detectable LTR.
  -np,  --no_postproc          Stop after the final genome; skip all plots and reports.
  -kt,  --keep_temps           Keep per-generation intermediate files.

--------------------------------------- GENERAL -------------------------------------------------
  -s,   --seed N               Random seed (default: 42).
  -t,   --threads N            Threads (default: 4).
  -h,   --help                 Show this help and exit.

Examples:
  # Burn-in only: 100 Mb, 1 chromosome, 4% CDS, 10% intact TE
  $(basename "$0") --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10

  # Variable-rate forward simulation
  $(basename "$0") -sz 135Mb -cn 5 -P 25 -itp 21 -ir 1.1e-6 -dr 4e-6 -br 1e-7 -m 7e-9 -ge 40000 -st 10000 -t 20

  # Fixed-rate run, then extend it with new rates via --continue
  $(basename "$0") -sz 113Mb -cn 20 -P 20 -itn 6000 -F 3e-11,5e-11 -ge 300000 -st 100000 -t 10
  $(basename "$0") -ge 400000 -st 100000 -F 7e-11,1e-11 --continue
EOF
}

#—-----------------------------------------
# No arguments -> concise usage.  -h/--help (any form) -> full help.
if [[ $# -eq 0 ]]; then
  print_usage
  exit 0
fi
if [[ "$1" == "-h" || "$1" == "--h" || "$1" == "-help" || "$1" == "--help" ]]; then
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

ensure_mutator   # verify / recompile ltr_mutator before anything else runs

# --- Parse command-line options ---
# Initialize flags and new parameters with defaults.
cont_flag=0
keep_temps=0
burnin_only=0
no_postproc=0
euch_bias_insert=1.0
euch_bias_excise=1.0
euch_buffer=10000
model="K2P"
ex_ltr=0
disable_genes=0
clean_lib_arg=""
# New TE inserter mutation options defaults:
TE_mut_k=10
TE_mut_Mmax=20
sel_coeff=0
TsTv=1.0
promoter_boundary=0
mutation_bins=""

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
    -cl|--clean_lib)
      clean_lib_arg="$2"
      shift; shift;;
    -m|--mutation_rate)
      mutation_rate="$2"
      shift; shift;;
    -r|--TE_ratio)
      TE_ratio="$2"
      shift; shift;;
    -itn|--intact_TE_num)
      intact_TE_num="$2"
      shift; shift;;
    -itp|--intact_TE_percent)
      intact_TE_percent="$2"
      shift; shift;;
    -ftn|--frag_TE_num)
      frag_TE_num="$2"
      shift; shift;;
    -ftp|--frag_TE_percent)
      frag_TE_percent="$2"
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
    -TsTv|--TsTv)
      TsTv="$2"
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
    -mb|--mutation_bins)
      mutation_bins="$2"
      shift; shift;;
    -bo|--burnin_only)
      burnin_only=1
      shift;;
    -ex|--ex_LTR)
      ex_ltr=1
      shift;;
    -np|--no_postproc)
      no_postproc=1
      shift;;
    -dg|--disable_genes)
      disable_genes=1
      shift;;
    -mxgs|--max_size)
      max_size="$2"
      shift; shift;;
    -mngs|--min_size)
      min_size="$2"
      shift; shift;;
    -sf|--sel_coeff)
      sel_coeff="$2"
      shift; shift;;
    -pgs|--pergen_select)
      pergen_select="$2"
      shift; shift;;
    -pb|--promoter-boundary)
      promoter_boundary="$2"
      shift; shift;;
    -md|--model)
      model="$2"
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
cds="${cds:-${TOOL_DIR}/data/TAIR10.cds.fa}"
chr_number="${chr_number:-4}"
size="${size:-400Mb}"
seed="${seed:-42}"
TE_lib="${TE_lib:-${TOOL_DIR}/data/maize_rice_arab_curated_TE.lib.gz}"
mutation_rate="${mutation_rate:-1.3e-8}"
TE_ratio="${TE_ratio:-${TOOL_DIR}/ratios.tsv}"
threads="${threads:-4}"
insert_rate="${insert_rate:-1e-8}"
birth_rate="${birth_rate:-1e-3}"
delete_rate="${delete_rate:-1e-7}"
solo_rate="${solo_rate:-95}"
k="${k:-10}"
sigma="${sigma:-1.0}"
max_size="${max_size:-}" # If not set, remains empty
min_size="${min_size:-}" # If not set, remains empty
TsTv="${TsTv:-1.0}"
pergen_select="${pergen_select:-2}"   # how many evenly spaced generations to select (incl. burnin and max)
if [[ -z "$intact_TE_num" && -z "$intact_TE_percent" ]]; then
  intact_TE_percent=20
fi
if [[ -z "$frag_TE_num" && -z "$frag_TE_percent" ]]; then
  frag_TE_percent=0
fi

# --- Validate mutually exclusive CDS options ---
if [[ -n "$cds_num" && -n "$cds_percent" ]]; then
  echo "Error: Provide either --cds_num (or -N) or --cds_percent (or -P), not both." | tee -a "$ERR"
  exit 1
fi

# --- Validate mutually exclusive TE insertion options (burn-in) ---
if [[ -n "$intact_TE_num" && -n "$intact_TE_percent" ]]; then
  echo "Error: Provide either --intact_TE_num or --intact_TE_percent, not both." | tee -a "$ERR"
  exit 1
fi
if [[ -n "$frag_TE_num" && -n "$frag_TE_percent" ]]; then
  echo "Error: Provide either --frag_TE_num or --frag_TE_percent, not both." | tee -a "$ERR"
  exit 1
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
if [[ -z "$clean_lib_arg" ]]; then
  TE_lib=$(decompress_if_gz "$TE_lib")
fi
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
# TE Library Processing
###############################################################################
if [[ -n "$clean_lib_arg" ]]; then
    clean_lib=$(decompress_if_gz "$clean_lib_arg")
    echo "=== Skipping TE Library Processing; using user-provided clean library: ${clean_lib} ===" | tee -a "$LOG"
else
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
fi

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

    # (1b) Insert TEs into the synthetic genome using the updated parallel inserter.
    cmd="python ${BIN_DIR}/shared_ltr_inserter_parallel.py -genome backbone.fa -TE ${clean_lib}"
    # INTACT controls
    if [[ -n "$intact_TE_percent" ]]; then
        cmd+=" -p_intact ${intact_TE_percent}"
    elif [[ -n "$intact_TE_num" ]]; then
        cmd+=" -n_intact ${intact_TE_num}"
    fi
    # FRAGMENTED controls
    if [[ -n "$frag_TE_percent" ]]; then
        cmd+=" -p_frag ${frag_TE_percent}"
    elif [[ -n "$frag_TE_num" ]]; then
        cmd+=" -n_frag ${frag_TE_num}"
    fi
    
    cmd+=" -TsTv ${TsTv}"
    cmd+=" -bed backbone.bed -output burnin -seed ${seed} -TE_ratio ${TE_ratio} -stat_out burnin.stat"

    if [[ -n "$mutation_bins" ]]; then
        # Use user-provided bins: script will ignore k/Mmax and disable the decay plot internally
        cmd+=" -mutation_bins ${mutation_bins} -m ${threads}"
    else
        # Default behavior: exponential decay model + PDF plot
        cmd+=" -k ${TE_mut_k} -Mmax ${TE_mut_Mmax} -pdf_out burnin_mut_dist.pdf -m ${threads}"
    fi
    echo "Running: $cmd" | tee -a "$LOG"
    eval $cmd >> "$LOG" 2>> "$ERR"
    if [ $? -ne 0 ]; then
        echo "Error running shared_ltr_inserter_parallel2.py" | tee -a "$ERR"
        exit 1
    fi
fi

if [[ "$burnin_only" -eq 1 ]]; then
  echo "Burn-in complete; exiting due to --burnin_only." | tee -a "$LOG"
  exit 0
fi

###############################################################################
# Phase 2: Looping Generations
###############################################################################


if [[ -n "$max_size" ]]; then
  raw="$max_size"
  unit="${raw: -1}"              # last character: M, G, or digit
  num="${raw%[MGmg]}"            # strip a possible M/G
  case "$unit" in
    [Mm]) max_bytes=$(( num * 1024 * 1024 )) ;;
    [Gg]) max_bytes=$(( num * 1024 * 1024 * 1024 )) ;;
    *)    # assume bytes if purely numeric
      if [[ "$raw" =~ ^[0-9]+$ ]]; then
        max_bytes=$raw
      else
        echo "Error: invalid format for --max_size: $raw" | tee -a "$ERR"
        exit 1
      fi
      ;;
  esac
  echo "Max genome size set to $raw → $max_bytes bytes" | tee -a "$LOG"
fi

if [[ -n "$min_size" ]]; then
  raw="$min_size"
  unit="${raw: -1}"              # last character: M, G, or digit
  num="${raw%[MGmg]}"            # strip a possible M/G
  case "$unit" in
    [Mm]) min_bytes=$(( num * 1024 * 1024 )) ;;
    [Gg]) min_bytes=$(( num * 1024 * 1024 * 1024 )) ;;
    *)    # assume bytes if purely numeric
      if [[ "$raw" =~ ^[0-9]+$ ]]; then
        min_bytes=$raw
      else
        echo "Error: invalid format for --min_size: $raw" | tee -a "$ERR"
        exit 1
      fi
      ;;
  esac
  echo "Min genome size set to $raw → $min_bytes bytes" | tee -a "$LOG"
fi

# If both are set, make sure the interval makes sense.
if [[ -n "$min_size" && -n "$max_size" ]]; then
  if (( min_bytes > max_bytes )); then
    echo "Error: --min_size ($min_bytes bytes) cannot be greater than --max_size ($max_bytes bytes)." | tee -a "$ERR"
    exit 1
  fi
fi

echo "=== Phase 2: Looping Generations ===" | tee -a "$LOG"

# Conserved cut-and-paste: a relocation-debt sidecar (cutpaste_debt.tsv) is
# produced by the exciser and consumed by the next generation's inserter.
# On a fresh run (NOT --continue) remove any leftover from a prior unrelated
# run so it cannot leak into generation 1. --continue intentionally keeps it.
if [[ "$cont_flag" -ne 1 ]]; then
  rm -f cutpaste_debt.tsv
fi

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
prev_lib="${clean_lib}"

last_gen_done=0
# Arrays to track genome sizes across iterations for exponential projection.
# On resume, backfill from existing gen*_final.fasta files so projection has history.
genome_size_iters=()
genome_size_bytes=()
if [[ "$cont_flag" -eq 1 && $start_iter -gt 1 ]]; then
  echo "Backfilling genome size history from existing generations for projection..." | tee -a "$LOG"
  for (( bi=1; bi<start_iter; bi++ )); do
    bf="gen$(( bi * step ))_final.fasta"
    if [[ -f "$bf" ]]; then
      if [[ "$OS" == "Darwin" ]]; then
        bsz=$(stat -f%z "$bf")
      else
        bsz=$(stat -c%s "$bf")
      fi
      genome_size_iters+=("$bi")
      genome_size_bytes+=("$bsz")
      echo "  iteration $bi (gen$(( bi * step ))): ${bsz} bytes" | tee -a "$LOG"
    fi
  done
  echo "Backfilled ${#genome_size_iters[@]} data points for projection." | tee -a "$LOG"
fi

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
        prev_genome="burnin.fasta"
        prev_bed="burnin.bed"
      fi
    else
      prev_genome="burnin.fasta"
      prev_bed="burnin.bed"
    fi
  else
    prev_genome="gen$(( (i-1) * step ))_final.fasta"
    prev_bed="gen$(( (i-1) * step ))_final.bed"
  fi

  # (2a) Mutate the genome.
  mut_prefix="gen${current_gen}_mut"

  # choose mode and bed based on whether this is the first iteration
  if [ $i -eq 1 ]; then
    mode=2
  else
    mode=3
  fi

  cmd="${BIN_DIR}/${mutator_exec} \
    -fasta ${prev_genome} \
    -bed   ${prev_bed} \
    -rate ${mutation_rate} \
    -generations ${step} \
    -mode ${mode} \
    -threads ${threads} \
    -TsTv ${TsTv} \
    -seed ${current_seed} \
    -out_prefix ${mut_prefix}"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running ltr_mutator for generation ${current_gen}" | tee -a "$ERR"
    exit 1
  fi

  # strip out column 7 from the .bed in-place
  # BUT do not do this on the user-provided BED passed via --bed
  if [[ -z "$input_bed" || "$prev_bed" != "$input_bed" ]]; then
    tmp="${prev_bed}.tmp"
    cut -f1-6 "${prev_bed}" > "$tmp" && mv "$tmp" "${prev_bed}"
  else
    echo "Skipping column 7 stripping for user-provided BED file: ${prev_bed}" | tee -a "$LOG"
  fi

  # (2b) Insert new TEs (allowing for nesting) using the parallel nest inserter.
  # Now using nest_inserter_parallel.py and adding --disable_genes if specified.
  nest_prefix="gen${current_gen}_nest"
  cmd="python ${BIN_DIR}/nest_inserter_parallel.py --genome ${mut_prefix}.fa --TE ${prev_lib} --generations ${step} --bed ${prev_bed} --output ${nest_prefix} --seed ${current_seed} --rate ${insert_rate} ${extra_fix_in} --TE_ratio ${TE_ratio} -bf burnin.stat --birth_rate ${birth_rate}"
  cmd+=" --euch_het_bias ${euch_bias_insert} --euch_het_buffer ${euch_buffer} -m ${threads}"
  if [[ $disable_genes -eq 1 ]]; then
    cmd+=" --disable_genes"
  fi
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running nest_inserter_parallel.py for generation ${current_gen}" | tee -a "$ERR"
    exit 1
  fi

  # (2c) Purge some TEs and convert intact TEs to soloLTRs using TE excision.
  # Updated to use TE_exciser_parallel.py with new parameters.
  cmd="python ${BIN_DIR}/TE_exciser_parallel.py --genome ${nest_prefix}.fasta --bed ${nest_prefix}.bed --rate ${delete_rate} --generations ${step} --soloLTR_freq ${solo_rate} ${extra_fix_ex} --output gen${current_gen}_final --seed ${current_seed} --sigma ${sigma} --k ${k} --sel_coeff ${sel_coeff} -m ${threads} --euch_het_buffer ${euch_buffer} --euch_het_bias ${euch_bias_excise} --promoter-boundary ${promoter_boundary} --TE_ratio ${TE_ratio}"

  if [ $i -ne 1 ]; then
    cmd+=" --no_fig"
  fi
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error running TE_exciser_parallel.py for generation ${current_gen}" | tee -a "$ERR"
    exit 1
  fi

  # If --keep_temps is not provided, remove intermediate files.
  if [[ "$keep_temps" -ne 1 ]]; then
    echo "Removing temporary files for generation ${current_gen}" | tee -a "$LOG"
    rm -f "${mut_prefix}.fa" "${nest_prefix}.bed" "${nest_prefix}.fasta" "backbone.fa" "backbone.cds" "backbone.bed" "lib.fa"
  fi
  
  # (2d) Build the new per‑gen TE library
  echo "=== Extracting intact TEs for generation ${current_gen} into lib file ===" | tee -a "$LOG"
  cmd="python ${BIN_DIR}/extract_intact_TEs.py \
    --genome gen${current_gen}_final.fasta \
    --bed    gen${current_gen}_final.bed \
    --weight_by ${clean_lib} --exclude_missing_ltr_len \
    --exclude_truncated \
    --out_fasta gen${current_gen}_final.lib"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
      echo "Error running extract_intact_TEs.py for generation ${current_gen}" | tee -a "$ERR"
      exit 1
  fi

  # Remove .fai index created by pysam (temp artifact)
  rm -f "gen${current_gen}_final.fasta.fai"

  # (2e) Update pipeline report
  echo "=== Step 2e: Updating insertion/deletion pipeline report ===" | tee -a "$LOG"
  cmd="python ${UTIL_DIR}/log_to_report.py -in ${LOG} -out pipeline.report"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
      echo "Error running log_to_report.py" | tee -a "$ERR"
      exit 1
  fi

  # clean up the *previous* lib file unless the user wants to keep temps
  if [[ "$keep_temps" -ne 1 && "$prev_lib" != "${clean_lib}" ]]; then
    rm -f "$prev_lib"
  fi

  # point to the newly built library for the next iteration
  prev_lib="gen${current_gen}_final.lib"
  
  # record that we successfully reached this generation
  last_gen_done=$current_gen
  
  # --- Exponential regression with empirical prediction interval ---
  # After observing >= 10% of total iterations, project the terminal genome
  # size via exponential regression (log(y) = a + bx).  Uncertainty bounds
  # use empirically calibrated prediction intervals (LOO-validated at ~93%
  # coverage) instead of the theoretical OLS PI which is miscalibrated for
  # this extrapolation task.
  #
  # The empirical error envelope is asymmetric: models almost never over-
  # predict (q97.5 ≈ 0%) but can substantially under-predict (q2.5 varies
  # from -56% at 10% observed to -25% at 50%).  This means:
  #   "too large" stops are high-confidence (tight upper bound)
  #   "too small" stops require wider margins (conservative lower bound)
  #
  # Termination requires the ENTIRE empirical interval to fall outside the
  # user-specified bounds, making false termination rare.
  if [[ -n "$max_size" || -n "$min_size" ]]; then
    fasta="gen${current_gen}_final.fasta"
    if [[ -f "$fasta" ]]; then
      if [[ "$OS" == "Darwin" ]]; then
          actual_bytes=$(stat -f%z "$fasta")  # macOS
      else
          actual_bytes=$(stat -c%s "$fasta")  # Linux
      fi

      # Track genome size for trend projection
      genome_size_iters+=("$i")
      genome_size_bytes+=("$actual_bytes")

      # Minimum data: 10% of total iterations (empirically validated threshold)
      min_n=$(( iterations / 10 ))
      [[ $min_n -lt 3 ]] && min_n=3

      if [[ ${#genome_size_iters[@]} -ge $min_n ]]; then
        projection_result=$(python -c "
import sys, math

xs = [int(x) for x in sys.argv[1].split(',')]
ys = [int(y) for y in sys.argv[2].split(',')]
target_x = int(sys.argv[3])
total_iters = int(sys.argv[4])
n = len(xs)
frac = n / total_iters

# --- Exponential regression: log(y) = a + bx ---
log_ys = [math.log(y) for y in ys]
sum_x   = sum(xs)
sum_ly  = sum(log_ys)
sum_xly = sum(x * ly for x, ly in zip(xs, log_ys))
sum_x2  = sum(x * x for x in xs)
denom   = n * sum_x2 - sum_x * sum_x

if denom == 0:
    proj = int(math.exp(sum_ly / n))
    print(f'{proj} {proj} {proj} 0.0')
    sys.exit(0)

b = (n * sum_xly - sum_x * sum_ly) / denom
a = (sum_ly - b * sum_x) / n
log_proj = a + b * target_x
proj = math.exp(log_proj)

# R-squared in log-space
ss_res = sum((ly - (a + b * x))**2 for x, ly in zip(xs, log_ys))
ss_tot = sum((ly - sum_ly / n)**2 for ly in log_ys)
r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

# --- Empirical prediction interval (calibrated from 61 trajectories) ---
# Lookup table: fraction -> (q2.5_SPE%, q97.5_SPE%) for exponential model.
# SPE = (projected - actual) / actual * 100.
# Derived via leave-one-out cross-validation (~93% coverage).
#   actual_lo = proj / (1 + q97.5/100)  [actual is at least this]
#   actual_hi = proj / (1 + q2.5/100)   [actual is at most this]
_margin_table = [
    (0.03, -84.6,  13.4),
    (0.05, -67.1,   4.2),
    (0.07, -58.3,   1.2),
    (0.10, -55.7,   0.3),
    (0.15, -50.0,   0.8),
    (0.20, -44.3,   0.4),
    (0.25, -40.3,   0.7),
    (0.30, -36.9,   0.4),
    (0.40, -30.2,   0.3),
    (0.50, -25.1,   0.1),
    (0.60, -20.5,   0.1),
    (0.70, -15.0,   0.1),
    (0.80, -10.0,   0.1),
    (0.90,  -5.0,   0.1),
]

# Interpolate margins for current fraction
q_lo, q_hi = _margin_table[-1][1], _margin_table[-1][2]
if frac < _margin_table[0][0]:
    q_lo, q_hi = _margin_table[0][1], _margin_table[0][2]
else:
    for j in range(len(_margin_table) - 1):
        f0, lo0, hi0 = _margin_table[j]
        f1, lo1, hi1 = _margin_table[j + 1]
        if f0 <= frac <= f1:
            t = (frac - f0) / (f1 - f0)
            q_lo = lo0 + t * (lo1 - lo0)
            q_hi = hi0 + t * (hi1 - hi0)
            break

# Empirical bounds on actual terminal genome size
denom_hi = 1.0 + q_hi / 100.0
denom_lo = 1.0 + q_lo / 100.0
actual_lo = int(proj / denom_hi) if denom_hi > 0 else 0
actual_hi = int(proj / denom_lo) if denom_lo > 0 else int(proj * 10)

print(f'{int(proj)} {actual_lo} {actual_hi} {r2:.4f}')
" "$(IFS=,; echo "${genome_size_iters[*]}")" "$(IFS=,; echo "${genome_size_bytes[*]}")" "$iterations" "$iterations")

        read -r proj_point proj_lo proj_hi proj_r2 <<< "$projection_result"
        frac_pct=$(( 100 * ${#genome_size_iters[@]} / iterations ))

        echo "Projection (iter $i, n=${#genome_size_iters[@]}, ${frac_pct}% observed, R²=${proj_r2}): point=${proj_point}, empirical 95% PI=[${proj_lo}..${proj_hi}]" | tee -a "$LOG"

        if [[ -n "$min_size" ]] && python3 -c "import sys; sys.exit(0 if int('$proj_hi') < int('$min_bytes') else 1)"; then
          echo "Empirical PI upper bound (${proj_hi} bytes) is below minimum (${min_bytes} bytes)." | tee -a "$LOG"
          echo "Confidently below target. Stopping at generation ${current_gen}." | tee -a "$LOG"
          break
        fi

        if [[ -n "$max_size" ]] && python3 -c "import sys; sys.exit(0 if int('$proj_lo') > int('$max_bytes') else 1)"; then
          echo "Empirical PI lower bound (${proj_lo} bytes) exceeds maximum (${max_bytes} bytes)." | tee -a "$LOG"
          echo "Confidently above target. Stopping at generation ${current_gen}." | tee -a "$LOG"
          break
        fi

        echo "Projection within bounds [${min_bytes:-0}..${max_bytes:-inf}]. Continuing." | tee -a "$LOG"
      fi
    else
      echo "Warning: expected output '$fasta' not found." | tee -a "$ERR"
    fi
  fi

done

# --- Dump genome size trajectory for model fitting ---
if [[ ${#genome_size_iters[@]} -gt 0 ]]; then
  trajectory_file="genome_size_trajectory.tsv"
  echo -e "iteration\tgenome_size_bytes" > "$trajectory_file"
  for (( ti=0; ti<${#genome_size_iters[@]}; ti++ )); do
    echo -e "${genome_size_iters[$ti]}\t${genome_size_bytes[$ti]}" >> "$trajectory_file"
  done
  echo "Genome size trajectory written to ${trajectory_file} (${#genome_size_iters[@]} data points)" | tee -a "$LOG"
fi

# The final generation's debt has no next generation to consume it (by
# design); remove the spent sidecar so the workdir is clean.
rm -f cutpaste_debt.tsv

echo "Pipeline completed at $(date)" | tee -a "$LOG"

if [[ "$no_postproc" -eq 1 ]]; then
  echo "Skipping post-processing due to --no_postproc." | tee -a "$LOG"
  exit 0
fi

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
    initial_fasta="burnin.fasta"
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
cmd="python ${UTIL_DIR}/stats_report.py --bed \$(ls gen*_final.bed burnin.bed | sort -V) --out_prefix stat"
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

# Ensure Kmer2LTR is available (clone into TOOL_DIR if missing)
if [[ ! -d "${TOOL_DIR}/Kmer2LTR" ]]; then
  echo "Cloning Kmer2LTR into ${TOOL_DIR}..." | tee -a "$LOG"
  (
    cd "${TOOL_DIR}" \
      && git clone https://github.com/cwb14/Kmer2LTR.git
  ) >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error cloning Kmer2LTR" | tee -a "$ERR"
    exit 1
  fi
fi

# Ensure Kmer2LTR is available (clone into TOOL_DIR if missing)
if [[ ! -d "${TOOL_DIR}/Kmer2LTR" ]]; then
  echo "Cloning Kmer2LTR into ${TOOL_DIR}..." | tee -a "$LOG"
  (
    cd "${TOOL_DIR}" \
      && git clone https://github.com/cwb14/Kmer2LTR.git
  ) >> "$LOG" 2>> "$ERR"
  if [ $? -ne 0 ]; then
    echo "Error cloning Kmer2LTR" | tee -a "$ERR"
    exit 1
  fi
fi

# Determine how many generations we actually ran
# (i.e. highest_gen / step)
# total_gens=$iterations
total_gens=$(( last_gen_done / step ))

# Build the list of available generations from files, treating burnin as generation 0.
# Always include 0 if burnin files exist.
gens=( )
labels=( )  # parallel array with "burnin" or "gen<GEN>_final"
if [[ -f "burnin.fasta" && -f "burnin.bed" ]]; then
  gens+=(0)
  labels+=("burnin")
fi
for f in gen*_final.fasta; do
  [[ -e "$f" ]] || continue
  if [[ $f =~ gen([0-9]+)_final\.fasta ]]; then
    gens+=("${BASH_REMATCH[1]}")
    labels+=("gen${BASH_REMATCH[1]}_final")
  fi
done

# Sort by generation and deduplicate, keeping labels in sync
if (( ${#gens[@]} == 0 )); then
  echo "No generations found for per-generation analysis; skipping." | tee -a "$LOG"
  selected_gens=()
else
  # Portable replacement for mapfile + process substitution (BSD/macOS safe)
  sorted_pairs="$( paste \
      <(printf "%s\n" "${gens[@]}") \
      <(printf "%s\n" "${labels[@]}") \
      | sort -n -u -k1,1 )"
  gens=( ); labels=( )
  while IFS=$'\t' read -r gen lab; do
    [ -n "$gen" ] || continue
    gens+=("$gen")
    labels+=("$lab")
  done <<< "$sorted_pairs"
fi

# Decide how many selections to make (at least 2 if we have ≥2 gens)
k="$pergen_select"
n="${#gens[@]}"
if (( n == 0 )); then
  selected_idx=( )
elif (( n == 1 )); then
  # Only one available—pick it
  selected_idx=(0)
else
  # Ensure k within [2, n]
  if (( k < 2 )); then k=2; fi
  if (( k > n )); then k="$n"; fi

  # Choose k evenly spaced indices in [0, n-1], always including 0 and n-1.
  # Round to nearest and deduplicate while preserving ascending order.
  # Use command substitution instead of mapfile
  idxs=( $(python - "$n" "$k" <<'PY'
import sys, math
n = int(sys.argv[1])
k = int(sys.argv[2])
# positions across [0, n-1]
step = (n - 1) / (k - 1) if k > 1 else 0
raw = [round(i * step) for i in range(k)]
# de-duplicate while preserving order
seen = set()
out = []
for v in raw:
    if v < 0: v = 0
    if v > n-1: v = n-1
    if v not in seen:
        seen.add(v)
        out.append(v)
print("\n".join(map(str, out)))
PY
) )
  selected_idx=( "${idxs[@]}" )
fi

# Build the final selection arrays
selected_gens=( )
selected_labels=( )
for idx in "${selected_idx[@]}"; do
  selected_gens+=( "${gens[$idx]}" )
  selected_labels+=( "${labels[$idx]}" )
done

# For logging: show selections in ascending by generation
echo "Selected generations (by gen number): ${selected_gens[*]}" | tee -a "$LOG"
echo "Selected labels: ${selected_labels[*]}" | tee -a "$LOG"

# Iterate in descending order (as your original loop expected)
# We’ll drive by labels to directly find the right files, and handle burnin specially.
# (No reliance on 'step' or iteration math.)
# Create an index array sorted by descending generation (BSD sort safe)
desc_idx=( $(for i in "${!selected_gens[@]}"; do echo "$i"; done | sort -nr) )

for i_idx in "${desc_idx[@]}"; do
  gen="${selected_gens[$i_idx]}"
  lab="${selected_labels[$i_idx]}"
  if [[ "$gen" -eq 0 ]]; then
    final_prefix="burnin"
  else
    final_prefix="$lab"   # e.g., gen200000_final
  fi

  echo "Processing per-generation analysis for ${final_prefix}" | tee -a "$LOG"

  # (a) Extract intact LTR sequences.
  cmd="python ${BIN_DIR}/extract_intact_LTR.py --bed ${final_prefix}.bed --genome ${final_prefix}.fasta --out_fasta ${final_prefix}_LTR.fasta"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd

  # (b1) Pull domains
  cmd="python ${BIN_DIR}/ltr_domain_puller.py ${final_prefix}_LTR.fasta.key.tsv ${final_prefix}_LTR.domain"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd

  # (b2) Kmer2LTR
  cmd="python ${TOOL_DIR}/Kmer2LTR/Kmer2LTR.py -p ${threads} -i ${final_prefix}_LTR.fasta -D ${final_prefix}_LTR.domain -o ${final_prefix}_LTR.tsv -u ${mutation_rate} --no-plot --purge-subdirs"
  echo "Running: $cmd" | tee -a "$LOG"
  eval $cmd

  # (b3) Clean up
  rm -f ${final_prefix}_LTR.domain ${final_prefix}_LTR.tsv.log ${final_prefix}_LTR.tsv.summary
done

# Run ltr_dens.py once after per-generation analyses.
cmd="python ${BIN_DIR}/ltr_dens.py --model ${model} --output all_LTR_density.pdf --miu ${mutation_rate} --gradient"
echo "Running: $cmd" | tee -a "$LOG"
eval $cmd

echo "Post-processing completed at $(date)" | tee -a "$LOG"

# END
