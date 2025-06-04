#!/bin/bash
#SBATCH --job-name=prinTE-array
#SBATCH --output=logs/output_%a.log
#SBATCH --error=logs/error_%a.log
#SBATCH --array=0-300%100   # Adjust max as needed or use sbatch --array=1-<N> externally. the starting number is the starting item in the job array (83 will start from the 83th job)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=15:00:00

# ===== Constants ====
MUTATION_RATE="1e-8"
BASE_DIR="./job_combinations_u${MUTATION_RATE}"
PRINTE_PATH="/users/bin/TESS/prinTE"
BED_FILE="/users/data/TESS/Crori1/crori.bed"
FASTA_FILE="/users/data/TESS/Crori1/crori.fa"
PARAMS_FILE="/users/data/TESS/Crori1/params_v2u${MUTATION_RATE}.list"
GE="200000"
ST="20000"
SEED="41"
THREADS="48"

# File size control
MAX_FILE_SIZE_MB=280
MAX_SIZE_DIFFERENCE_MB=20   # <--- CUSTOMIZABLE: Max allowed size diff (MB) between first and last gen_final.fasta
CHECK_INTERVAL=120
MAX_FILE_SIZE_BYTES=$((MAX_FILE_SIZE_MB * 1024 * 1024))
MAX_SIZE_DIFFERENCE_BYTES=$((MAX_SIZE_DIFFERENCE_MB * 1024 * 1024))

# Path to R script
PLOT_SCRIPT="/users/bin/TESS/prinTE/plot_divergence5.R"
REF_LIST="/users/data/TESS/Crori1/Crori1_AssemblyScaffolds.fasta.mod.pass.list"

# Create logs directory if not exists
mkdir -p logs

# ===== Read this task's parameters ====
TASK_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PARAMS_FILE})
if [[ -z "$TASK_LINE" ]]; then
    echo "No task found at line ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 0
fi

# Extract parameters from the line name
fix_first=$(echo "$TASK_LINE" | grep -oP 'fix1_\K[^_]+')
fix_second=$(echo "$TASK_LINE" | grep -oP 'fix2_\K[^_]+')
k=$(echo "$TASK_LINE" | grep -oP 'k_\K[^_]+')
solo_rate=$(echo "$TASK_LINE" | grep -oP 'solo_\K[^_]*')

JOB_DIR="$BASE_DIR/$TASK_LINE"
mkdir -p "$JOB_DIR"

cd "$JOB_DIR" || exit 1

# ===== Resume Check with oversized file check ====
has_large_file=false
for file in gen*_final.fasta; do
    if [ -f "$file" ]; then
        size=$(stat -c "%s" "$file")
        if [ "$size" -gt "$MAX_FILE_SIZE_BYTES" ]; then
            echo "?~Z| ?~O Large file detected: $file > ${MAX_FILE_SIZE_MB} MB. Will not resume."
            has_large_file=true
            break
        fi
    fi
done

# Final resume condition
if [ -f "all_LTR_density.pdf" ] && \
   [ -f "overall.solo_intact.txt" ] && \
   [ -n "$(ls -1 *_final.fasta.mod.pass.list 2>/dev/null)" ] && \
   ! $has_large_file; then
    echo "?~\~E Job already completed successfully. Skipping: $TASK_LINE"
    exit 0
else
    echo "?~_~T~D Running job: $TASK_LINE"
fi

# ===== Run prinTE command in background ====
/usr/bin/time -v bash ${PRINTE_PATH}/prinTE.sh \
  --bed ${BED_FILE} \
  --fasta ${FASTA_FILE} \
  -ge ${GE} \
  -st ${ST} \
  --fix ${fix_first},${fix_second} \
  --seed ${SEED} \
  --threads ${THREADS} \
  -k ${k} \
  --solo_rate ${solo_rate} \
  --mutation_rate ${MUTATION_RATE} \
  --continue &
PRINTE_PID=$!

ABORT=0

# ===== Monitor file sizes during execution ====
while true; do
    ps -p $PRINTE_PID > /dev/null
    if [ $? -ne 0 ]; then break; fi

    # Existing check for oversized files
    for file in gen*_final.fasta; do
        if [ -f "$file" ]; then
            size=$(stat -c "%s" "$file")
            if [ "$size" -gt "$MAX_FILE_SIZE_BYTES" ]; then
                echo "File $file exceeds limit ($MAX_FILE_SIZE_MB MB). Aborting prinTE.sh..."
                kill -9 $PRINTE_PID > /dev/null 2>&1
                wait $PRINTE_PID 2>/dev/null
                ABORT=1
                break 2
            fi
        fi
    done

    # New check: Size difference between first and last gen_final.fasta
    files=($(ls -v gen*_final.fasta 2>/dev/null))
    if [ ${#files[@]} -ge 2 ]; then
        first_file="${files[0]}"
        last_file="${files[-1]}"
        if [ -f "$first_file" ] && [ -f "$last_file" ]; then
            first_size=$(stat -c "%s" "$first_file")
            last_size=$(stat -c "%s" "$last_file")
            size_diff=$(( last_size - first_size ))
            if [ "$size_diff" -lt 0 ]; then
                size_diff=$(( -size_diff ))  # Absolute value
            fi
            if [ "$size_diff" -ge "$MAX_SIZE_DIFFERENCE_BYTES" ]; then
                echo "Size difference between first ($first_file) and last ($last_file) files exceeds ${MAX_SIZE_DIFFERENCE_MB}MB ($size_diff bytes). Aborting..."
                kill -9 $PRINTE_PID > /dev/null 2>&1
                ABORT=1
                break
            fi
        fi
    fi

    sleep $CHECK_INTERVAL
done

wait $PRINTE_PID 2>/dev/null

if [ $ABORT -eq 1 ]; then
    echo "Job aborted due to file size growth trend or oversized output."
    exit 1
fi

# ===== Post-processing only if not aborted ====
for i in *_final.fasta; do
    final_prefix=${i%.fasta}
    python ~/bin/TESS/prinTE/bin/extract_intact_LTR.py \
        --bed ${final_prefix}.bed \
        --genome ${final_prefix}.fasta \
        --out_fasta ${final_prefix}_LTR.fasta &
done
wait

for i in *_final.fasta; do
    final_prefix=${i%.fasta}
    if [ ! -f "${final_prefix}_LTR.tsv" ]; then
        nohup python ~/bin/TESS/prinTE/bin/seq_divergence.py \
            -i ${final_prefix}_LTR.fasta \
            -o ${final_prefix}_LTR.tsv \
            -t 5 \
            --min_align 100 \
            --max_off 20 \
            --miu ${MUTATION_RATE} \
            --blast_outfmt '6 qseqid sseqid sstart send slen qstart qend qlen length nident btop' &
    fi
done
wait

BIN_DIR=~/bin/TESS/prinTE/bin
python ${BIN_DIR}/ltr_dens.py --model K2P --output all_LTR_density.pdf --miu ${MUTATION_RATE}

for i in *final.fasta; do
    nohup perl ~/bin/EDTA/EDTA_raw.pl --genome $i --type LTR -t 5 &
done
wait

for i in *final.fasta.mod; do
    cd $i.EDTA.raw/LTR/
    ln -s ../../$i
    if [ ! -f "$i.out" ]; then
        nohup RepeatMasker -pa 5 -q -div 40 -cutoff 225 -lib "$i.LTRlib.fa" "$i" &
    fi
    cd ../../
done
wait

for i in */LTR/*lib.fa; do
    sh ~/bin/LTR_retriever/bin/solo_intact_ratio_wrap.sh $i ${i%.LTRlib.fa}.out
done > overall.solo_intact.txt

# Use correct path to plot script
Rscript "$PLOT_SCRIPT" "$REF_LIST" *raw/LTR/*mod.pass.list "${JOB_DIR}_all"
Rscript "$PLOT_SCRIPT" "$REF_LIST" gen160000*/LTR/*mod.pass.list gen180000*/LTR/*mod.pass.list gen200000*/LTR/*mod.pass.list "$JOB_DIR"
