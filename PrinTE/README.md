# prinTE: TE Evolution Simulation Pipeline

**prinTE** is a comprehensive wrapper script (bash) that orchestrates a two‐phase transposable element (TE) evolution simulation pipeline, followed by supplemental post‐processing and reporting. It automates:

1. **Phase 1 (Burn-in)**  
   Builds a synthetic genome, inserts an initial set of TEs, and generates “burn-in” outputs (e.g., `burnin.fa`, `burnin.bed`), unless the user provides their own starting FASTA/BED files or restarts from a previous run.

2. **Phase 2 (Looping Generations)**  
   For each generation (in user-specified steps), the pipeline:
   - Mutates the genome (LTR mutator).  
   - Inserts new TEs (nested insertion).  
   - Excises TEs (parallel TE exciser).  
   - Extracts intact TEs to form the next-generation TE library.  
   - Logs progress and checks for maximum genome size.

3. **Supplemental Post-Processing**  
   After all generations complete (or stop early due to size limits), the pipeline runs a suite of reporting and plotting utilities:
   - TE fraction over time  
   - Solo vs. intact TE proportions  
   - Superfamily counts  
   - Category bar plots  
   - Genome size trajectory  
   - Per-generation LTR extraction and divergence analyses  
   - Overall LTR density plot  

---

## Table of Contents

1. [Prerequisites](#prerequisites)  
2. [Directory Layout](#directory-layout)  
3. [Installation & Setup](#installation--setup)  
4. [Running the Pipeline](#running-the-pipeline)  
   - [Basic Full-Pipeline Invocation](#basic-full-pipeline-invocation)  
   - [Burn-In Only](#burn-in-only)  
   - [Continuing a Previous Run](#continuing-a-previous-run)  
   - [Skipping Burn-In with Custom BED/FASTA](#skipping-burn-in-with-custom-bedfasta)  
   - [Keeping Temporary Files](#keeping-temporary-files)  
   - [Adjusting Chromatin Bias & Deletion](#adjusting-chromatin-bias--deletion)  
   - [Custom TE Mutation Model](#custom-te-mutation-model)  
   - [Limiting Maximum Genome Size](#limiting-maximum-genome-size)  
5. [Options & Flags](#options--flags)  
   - [Required Arguments (Full Pipeline)](#required-arguments-full-pipeline)  
   - [Burn-In Phase Options](#burn-in-phase-options)  
   - [Fixed vs. Variable TE Indel Rates](#fixed-vs-variable-te-indel-rates)  
   - [Chromatin Bias & Buffers](#chromatin-bias--buffers)  
   - [TE Mutation Parameters (Burn-In & Phase 2)](#te-mutation-parameters-burn-in--phase-2)  
   - [General-Use Options (Phases 1 & 2)](#general-use-options-phases-1--2)  
   - [Post-Processing & Reporting](#post-processing--reporting)  
6. [Output Files](#output-files)  
7. [Example Workflows](#example-workflows)  
   - [Simulate 3000 Generations (Step 1000)](#simulate-3000-generations-step-1000)  
   - [Run Burn-In Only](#run-burn-in-only)  
   - [Resume from Generation 2000](#resume-from-generation-2000)  
   - [Start from Custom BED/FASTA](#start-from-custom-bedfasta)  
   - [Parallel TE Insertion with Custom Rates & Bias](#parallel-te-insertion-with-custom-rates--bias)  
8. [Logging & Error Handling](#logging--error-handling)  
9. [Contact & Support](#contact--support)  

---

## Prerequisites

- **Operating System**: Linux/macOS (bash-compatible shell)  
- **Python ≥ 3.6** (with modules required by internal scripts):
  - `seq_divergence.py`  
  - `LTR_fasta_header_appender.py`  
  - `extract_intact_TEs.py`  
  - `synthetic_genome.py`  
  - `shared_ltr_inserter_parallel.py`  
  - `nest_inserter_parallel.py`  
  - `TE_exciser_parallel.py`  
  - `extract_intact_LTR.py`  
  - `ltr_mutator` (compiled binary)  
  - `log_to_report.py`  
  - Plotting utilities in `util/` (e.g., `plot_TE_frac.py`, `plot_solo_intact.py`, etc.)  

- **Bash utilities**:  
  - `bash`  
  - `cut`  
  - `gunzip`  
  - `stat`  
  - `sort`  
  - `ls`  
  - `bc`  

- **GNU Coreutils**, **GNU grep**, **sort**, **ls**  
- **Optional**: `Rscript` or other plotting backends (if any Python plotting scripts invoke R)  

Before running, ensure that:

1. The entire **TESS/prinTE** repository (or folder structure) is cloned/copied.  
2. You have execute permissions on `prinTE.sh` (e.g., `chmod +x prinTE.sh`).  
3. Python dependencies for each helper script are installed (check each script for `import` statements).  

   ```bash
   pip install biopython numpy pandas matplotlib

Before running, ensure that:

1. The entire **TESS/prinTE** repository (or folder structure) is cloned/copied.
2. You have execute permissions on `prinTE.sh` (e.g., `chmod +x prinTE.sh`).
3. Python dependencies for each helper script are installed (check each script for `import` statements).

---

## Directory Layout

```text
TESS/
└── prinTE/
    ├── prinTE.sh                    # This wrapper script
    ├── bin/                         # Core simulation binaries & scripts
    │   ├── ltr_mutator              # Compiled C/C++ binary
    │   ├── seq_divergence.py
    │   ├── LTR_fasta_header_appender.py
    │   ├── extract_intact_TEs.py
    │   ├── synthetic_genome.py
    │   ├── shared_ltr_inserter_parallel.py
    │   ├── nest_inserter_parallel.py
    │   ├── TE_exciser_parallel.py
    │   ├── extract_intact_LTR.py
    │   └── ...                      # (other Python scripts)
    ├── util/                        # Supplemental post‐processing and plotting
    │   ├── log_to_report.py
    │   ├── plot_TE_frac.py
    │   ├── plot_solo_intact.py
    │   ├── stats_report.py
    │   ├── plot_superfamily_count.py
    │   ├── plot_category_bar.py
    │   ├── genome_plot.py
    │   └── ltr_dens.py
    └── README.md (→ this file)
Installation & Setup
Clone or copy the prinTE directory to your local machine:
git clone https://github.com/YourOrg/TESS.git
cd TESS/prinTE
chmod +x prinTE.sh
Install Python dependencies (if not already installed). Many helper scripts rely only on standard Python modules, but confirm:
pip install biopython numpy pandas matplotlib  # if required by scripts
Verify that bin/ltr_mutator is executable:
ls -l bin/ltr_mutator
# If not executable:
chmod +x bin/ltr_mutator
Confirm you have all necessary Python modules:
python -c "import sys; print(sys.version)"
# Try importing packages used by helper scripts, e.g.:
python - <<EOF
import argparse, subprocess, Bio
import numpy as np
EOF
Optional: Add prinTE.sh to your $PATH or create an alias:
ln -s "$(pwd)/prinTE.sh" /usr/local/bin/prinTE
Running the Pipeline
At a minimum, to run the full two‐phase pipeline, you must supply:
--step (-st): Generation step size (number of generations simulated per iteration)
--generation_end (-ge): Total number of generations to simulate (must be a multiple of step)
Everything else has sensible defaults (genome size = 400 Mb, 4 chromosomes, TE counts, etc.). Below are different usage patterns and example commands.
Basic Full-Pipeline Invocation
./prinTE.sh \
  --step 1000 \
  --generation_end 3000
This will:
Run Phase 1 (burn‐in) to build a 4-chromosome, 400 Mb synthetic genome with 2000 TE insertions (default).
Iterate Phase 2 for 3 iterations (1000, 2000, 3000), producing gen1000_final.fasta, gen2000_final.fasta, gen3000_final.fasta, and corresponding .bed & .lib files.
Generate pipeline reports and plotting files in the current directory.
Burn-In Only
To run only the burn‐in phase and exit immediately (no looping generations):
./prinTE.sh --burnin_only
Does not require --step or --generation_end.
Outputs backbone.fa, backbone.bed, burnin.fa, burnin.bed, burnin.stat, and burnin_mut_dist.pdf.
Exits after burn‐in; no Phase 2 or post‐processing is executed.
Continuing a Previous Run
If your pipeline was interrupted or you want to resume from the last completed generation, use --continue (or -x):
./prinTE.sh \
  --step 1000 \
  --generation_end 5000 \
  --continue
The script scans for gen<NUMBER>_final.fasta files in the working directory.
If e.g. gen2000_final.fasta and gen2000_final.bed exist, it will start Phase 2 at iteration 3 (generation 3000).
If no previous outputs are found, it starts from iteration 1 as usual.
Skipping Burn-In with Custom BED/FASTA
If you already have a starting genome (.fasta) and TE coordinate file (.bed) (e.g., from an external source), supply both --bed and --fasta to skip burn‐in entirely:
./prinTE.sh \
  --bed my_start.bed \
  --fasta my_start.fasta \
  --step 1000 \
  --generation_end 2000
skip_burnin is automatically enabled.
Phase 2 “first iteration” will use my_start.fasta & my_start.bed in place of burnin.fa/burnin.bed.
Note: Both --bed and --fasta must be provided together; otherwise, an error is thrown.
Keeping Temporary Files
By default, intermediate mutation and insertion FASTA/BED files (e.g., gen1000_mut.fa, gen1000_nest.fasta, etc.) are deleted after each generation. To keep them for debugging or inspection, use --keep_temps (or -kt):
./prinTE.sh \
  --step 500 \
  --generation_end 1500 \
  --keep_temps
Temporary files from each iteration will not be removed.
This is helpful if you suspect errors or want to inspect intermediate states.
Adjusting Chromatin Bias & Deletion
The pipeline now supports user‐friendly flags for insertion and deletion biases instead of the old --euch-bias & --euch-buffer:
./prinTE.sh \
  --step 1000 \
  --generation_end 3000 \
  --chromatin_bias_insert 2.0 \
  --chromatin_buffer 5000 \
  --chromatin_bias_delete 0.5
--chromatin_bias_insert (default 1.0): Multiplies TE insertion preference in euchromatic regions.
--chromatin_buffer (default 10000): Upstream/downstream buffer (in bp) around genes for insertion bias.
--chromatin_bias_delete (default 1.0): Bias for excision in euchromatin.
These parameters are passed to:
nest_inserter_parallel.py via --euch_het_bias & --euch_het_buffer.
TE_exciser_parallel.py via --euch_het_buffer & --euch_het_bias.
Custom TE Mutation Model
By default, LTR dating uses the Kimura 2‐Parameter (K2P) model. You can choose between:
raw: No correction
K2P: Kimura 2-parameter (default)
JC69: Jukes-Cantor
Use --model (or -md) to set your preference:
./prinTE.sh \
  --step 1000 \
  --generation_end 3000 \
  --model JC69
This affects the final LTR density plot (ltr_dens.py --model).
Limiting Maximum Genome Size
To stop Phase 2 early if the simulated genome exceeds a disk‐size threshold, supply --max_size:
./prinTE.sh \
  --step 1000 \
  --generation_end 5000 \
  --max_size 500M
--max_size 500M → stops if any gen<GEN>_final.fasta exceeds 500 MB (i.e., 500 × 1024 × 1024 bytes).
Acceptable units:
Suffix M or m for megabytes (× 1024²).
Suffix G or g for gigabytes (× 1024³).
No suffix → interpreted as raw bytes.
If threshold is reached, Phase 2 exits at the current generation, but post‐processing still runs on completed generations.
Options & Flags
Below is a complete list of prinTE.sh options, grouped by functionality. Each flag includes its shorthand, default value, and description.
Required Arguments (Full Pipeline)
These are only required if --burnin_only is not used.
-st, --step <INT>
Generation step size (number of generations simulated per iteration).
Example: --step 1000.
-ge, --generation_end <INT>
Total number of generations to simulate (must be an exact multiple of step).
Example: --generation_end 3000.
Burn-In Phase Options
-c, --cds <PATH>
Path to a CDS (coding sequence) FASTA file.
Default: ${TOOL_DIR}/TAIR10.cds.fa.
-N, --cds_num <INT>
Number of CDS sequences to insert (mutually exclusive with --cds_percent).
-P, --cds_percent <FLOAT>
Percent of the genome to be CDS (mutually exclusive with --cds_num).
-n, --TE_num <INT>
Number of TE insertions during burn‐in. Default: 2000.
-p, --TE_percent <FLOAT>
Percent of the genome to be TEs during burn‐in (mutually exclusive with --TE_num).
-cn, --chr_number <INT>
Number of chromosomes to simulate. Default: 4.
-sz, --size <STRING>
Genome size (e.g., 400Mb, 1Gb, 500Kb). Default: 400Mb.
-i, --TE_lib <PATH>
Input TE library FASTA for burn‐in. Default:
${TOOL_DIR}/combined_curated_TE_lib_ATOSZM_selected.fasta.
-m, --mutation_rate <FLOAT>
DNA mutation rate (per nucleotide per generation). Default: 1.3e-8.
-bo, --burnin_only
Run only the burn‐in phase and exit immediately (skipping Phase 2 and post‐processing).
-ex, --ex_LTR
Exclude LTR sequences without a domain hit in the TE library processing (passed to LTR_fasta_header_appender.py as -exclude_no_hits).
-dg, --disable_genes
Disable TE insertion into gene regions (only effective if using --fix for fixed rates).
Fixed vs. Variable TE Indel Rates
Fixed TE Indel Rate (mutually exclusive with variable options):
-F, --fix <INSERT,DELETE>
Provide fixed insertion & deletion rates as comma‐separated floats (e.g., 1e-9,1e-9).
→ Internally sets extra_fix_in="--fix_in <INSERT>" & extra_fix_ex="--fix_ex <DELETE>".
Variable TE Indel Rate:
-ir, --insert_rate <FLOAT>
TE insertion rate. Default: 1e-8.
-dr, --delete_rate <FLOAT>
TE deletion rate. Default: 1e-7.
-br, --birth_rate <FLOAT>
TE birth rate (used by nested inserter). Default: 1e-3.
-sc, --sigma <FLOAT>
Selection coefficient for gene insertions. Default: 1.0.
-sf, --sel_coeff <FLOAT>
Selection coefficient for TE excision (variable‐rate). Default: 0 (neutral).
0 = neutral (no bias)
0.1 = 2× bias
1 = 11× bias
Chromatin Bias & Buffers
-cbi, --chromatin_bias_insert <FLOAT>
Chromatin bias multiplier for TE insertion in euchromatin. Default: 1.0.
-cbd, --chromatin_bias_delete <FLOAT>
Chromatin bias multiplier for TE deletion in euchromatin. Default: 1.0.
-cb, --chromatin_buffer <INT>
Interval upstream/downstream of genes (in bp) used for chromatin bias. Default: 10000.
TE Mutation Parameters (Burn-In & Phase 2)
-tk, --TE_mut_k <INT>
Slope of exponential decay for TE mutation (used by shared_ltr_inserter_parallel.py). Default: 10.
-tmx, --TE_mut_Mmax <INT>
Maximum X on the exponential decay function for TE mutation. Default: 20.
General‐Use Options (Phases 1 & 2)
-s, --seed <INT>
Random seed for reproducibility. Default: 42.
-r, --TE_ratio <PATH>
TE ratio file (ratios of families). Default: ${TOOL_DIR}/ratios.tsv.
-t, --threads <INT>
Number of threads for parallel scripts. Default: 4.
-sr, --solo_rate <FLOAT>
Percent chance to convert an intact TE to solo LTR during excision. Default: 95.
-k, --k <INT>
TE length decay slope parameter for TE excision. Default: 10.
-b, --bed <PATH>
Input BED file for a custom starting generation (must be used with --fasta).
-f, --fasta <PATH>
Input FASTA file for a custom starting generation (must be used with --bed).
-x, --continue
Resume simulation from the last completed generation (skips burn‐in automatically).
-kt, --keep_temps
Keep intermediate temporary files instead of deleting them each loop.
-mgs, --max_size <STRING>
Maximum genome file size (e.g., 100M, 1G). If exceeded during Phase 2, the pipeline stops.
Suffix M/m → megabytes (× 2²⁰)
Suffix G/g → gigabytes (× 2³⁰)
No suffix → raw bytes
-md, --model <STRING>
Mutation model for LTR density plot: raw, K2P, or JC69. Default: K2P.
-ex, --ex_LTR
Exclude LTR sequences with no domain hits (passed to LTR_fasta_header_appender.py).
-h, --help
Display help message and exit.
Post-Processing & Reporting
These options are invoked automatically at the end of Phases 1 & 2; you generally do not need to specify them explicitly when calling prinTE.sh, but they define what is run:
Plot TE fraction over time
plot_TE_frac.py --bed <(initial_bed + sorted gen*_final.bed)> --fasta <(initial_fasta + sorted gen*_final.fasta)> --feature Intact_TE:SoloLTR:Fragmented_TE --out_prefix percent_TE
Plot solo vs. intact TE
plot_solo_intact.py --bed <(initial_bed + sorted gen*_final.bed)> --out_prefix solo_intact
Generate overall statistics report
stats_report.py --bed <(sorted gen*_final.bed)> --out_prefix stat
Plot superfamily counts
plot_superfamily_count.py
Plot category bar
plot_category_bar.py
Genome size trajectory
genome_plot.py
Per-generation LTR extraction & divergence
For a subset (≤ 4) of evenly distributed generations (first, middle, last), runs:
extract_intact_LTR.py --bed gen<GEN>_final.bed --genome gen<GEN>_final.fasta --out_fasta gen<GEN>_LTR.fasta
seq_divergence.py -i gen<GEN>_LTR.fasta -o gen<GEN>_LTR.tsv -t <threads> --min_align 100 --max_off 20 --miu <mutation_rate> --blast_outfmt '6 ...'
LTR density plot (all generations)
ltr_dens.py --model <model> --output all_LTR_density.pdf --miu <mutation_rate> --gradient
Output Files
During execution, the following files (and directories) will be generated:
Log files
pipeline.log → records stdout of each step
pipeline.error → records stderr/error messages
Phase 1 (Burn-In) Outputs
backbone.fa, backbone.bed → synthetic genome & coordinate skeleton
burnin.fa, burnin.bed → burn‐in genome and TE coordinates
burnin.stat → TE ratio & divergence statistics
burnin_mut_dist.pdf → Mutation distribution plot for burn‐in
Phase 2 (Looping Generations)
For each generation GEN (where GEN = step, 2×step, …, generation_end):
Intermediate files (unless --keep_temps is set):
gen<GEN>_mut.fa → mutated genome (before insertion)
gen<GEN>_nest.fasta, gen<GEN>_nest.bed → nested insertion input/output
Final outputs:
gen<GEN>_final.fasta → genome after excision
gen<GEN>_final.bed → TE coordinates after excision
gen<GEN>_final.lib → extracted intact TE library for next iteration
Supplemental Post-Processing
pipeline.report → summary report generated by log_to_report.py
Plots
percent_TE_*.pdf → TE fraction over time
solo_intact_*.pdf → Solo vs. intact TE proportions
stat_*.pdf / .tsv → Overall TE statistics
plot_superfamily_count_*.pdf → TE superfamily counts
plot_category_bar_*.pdf → Category bar plots
genome_plot_*.pdf → Genome size over time
gen<GEN>_LTR.fasta → Extracted LTR sequences (for selected generations)
gen<GEN>_LTR.tsv → LTR divergence statistics (for selected generations)
all_LTR_density.pdf → Combined LTR density plot
Temporary Files (removed by default; kept if --keep_temps)
*.tmp files created during cut -f1-6 on BED files
*.gz decompressed FASTA/BED, if original inputs were gzipped
Example Workflows
Below are several common use‐cases. Each example is a ready‐to‐copy command block.
1. Simulate 3000 Generations (Step 1000)
cd /path/to/TESS/prinTE

# Basic full pipeline, 3 total iterations: 1000, 2000, 3000
./prinTE.sh \
  --step 1000 \
  --generation_end 3000
What happens:
Build synthetic genome (4 chr, 400Mb) with 2000 TEs.
Run nested insertion & excision for generational steps: 1000, 2000, 3000.
Output:
gen1000_final.fasta, gen2000_final.fasta, gen3000_final.fasta
Per-phase TE libraries: gen1000_final.lib, etc.
Final plots & reports in current directory.
2. Run Burn-In Only
# Only run Phase 1 (burn-in), then exit
./prinTE.sh --burnin_only
What happens:
Creates: backbone.fa, backbone.bed, burnin.fa, burnin.bed, burnin.stat, burnin_mut_dist.pdf.
Then exits; no Phase 2 or post‐processing.
3. Resume from Generation 2000
Assume you previously ran 3000 generations but crashed halfway through generation 2000. You have:
gen1000_final.fasta
gen1000_final.bed
gen2000_final.fasta
gen2000_final.bed
To resume:
./prinTE.sh \
  --step 1000 \
  --generation_end 5000 \
  --continue
What happens:
Scans for existing gen*_final.fasta: finds gen2000_final.fasta.
Calculates start_iter = (2000/1000) + 1 = 3.
Continues with iteration 3 → generation 3000, then 4000, 5000.
All Phase 2 steps, plus post‐processing at the end.
4. Start from Custom BED/FASTA
If you have a pre‐constructed reference:
# Provide your own starting genome & TE coordinates
./prinTE.sh \
  --bed ./custom_start.bed \
  --fasta ./custom_start.fasta \
  --step 500 \
  --generation_end 1500
What happens:
skip_burnin=1 (no synthetic genome is built).
Phase 2 iteration 1 uses custom_start.fasta & custom_start.bed.
Iterations at generations 500, 1000, 1500 produce gen500_final.fasta, etc.
5. Parallel TE Insertion with Custom Rates & Bias
./prinTE.sh \
  --step 2000 \
  --generation_end 8000 \
  --threads 8 \
  --insert_rate 5e-8 \
  --delete_rate 5e-7 \
  --birth_rate 2e-3 \
  --chromatin_bias_insert 2.5 \
  --chromatin_buffer 20000 \
  --chromatin_bias_delete 0.8 \
  --seed 1234 \
  --keep_temps \
  --model raw \
  --max_size 1G
What happens:
Burn‐in with:
TE insertion rate = 5×10⁻⁸
TE deletion rate = 5×10⁻⁷
TE birth rate = 2×10⁻³
Chromatin bias multiplier = 2.5, buffer = 20 kb, deletion bias = 0.8
Phase 2 (4 iterations: 2000, 4000, 6000, 8000) with 8 threads.
Keeps intermediate gen<GEN>_mut.fa, gen<GEN>_nest.fasta, etc. for debugging.
Uses “raw” LTR mutation model in final density plot (no correction).
Stops early if gen<GEN>_final.fasta exceeds 1 GB.
Logging & Error Handling
pipeline.log: Standard output from all commands, with timestamps prefixed.
pipeline.error: Standard error from all commands.
On any non‐zero exit code from a sub‐step (e.g., failure in seq_divergence.py), prinTE.sh will log the error, print a message to pipeline.error, and exit 1 immediately.
A trap cleanup EXIT ensures that any decompressed temporary files (from gzipped inputs) are removed, unless you requested --keep_temps.
If an unknown option is passed, the script prints an “Unknown option” error to both terminal and pipeline.error, displays the help message, and exits.
