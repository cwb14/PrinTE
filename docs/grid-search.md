## Grid search: finding parameters without bias

Finding the insertion and deletion rates that make a simulation match a
target genome, without hand-tuning by eye.

[Back to the README](../README.md)

---

Both tracks need the same thing: the TE insertion/deletion rates (and solo-LTR and length-bias
settings) that make a simulation match a target. Rather than guess, PrinTE's `grid/` tools
**search the parameter space reproducibly** and score every simulation against a reference.

The search space is four dimensions: **insertion rate × deletion rate × solo ratio (`sr`) ×
length bias (`k`)**. Two ways to cover it:

- **`printe-grid`** (`printe.grid.guided_search`) - *active learning*: round 0 explores with
  Latin-hypercube sampling; later rounds train a random-forest surrogate on the results so far
  and propose the most informative next batch. This converges on the target with far fewer
  simulations, and it is what the published analysis used.
- **`nextflow run . --mode sweep`** - a fixed grid, scattered across whatever executor you
  point it at. Use this when you want to cover a defined space rather than let a surrogate
  choose, or when you want the portability. See [the Nextflow pipeline](nextflow.md).

**What the grid writes to disk:** every sampled parameter combination becomes its own directory
named `insertion_rates_<x>_deletion_rates_<y>_solo_ratio_<s>_length_bias_<k>/`, each holding that
run's `gen<N>_final.fasta/.bed/.lib`. The re-annotation and scoring steps below iterate over
these directories automatically.

**You will need a cleaned TE library** (`--clean-lib`). PrinTE writes `lib_clean.fa` into the
working directory whenever it processes a `--TE_lib`; reuse that file. To build one directly:

```bash
# Clean a TE library to pure ACGT and pre-process it once; reuse as --clean-lib / -cl.
zcat data/maize_rice_arab_curated_TE.lib.gz \
  | python -m printe.util.fix_non_ATGC - > lib_clean.fa
```

### Step 1 - seed and launch the search

```bash
python -m printe.grid.guided_search init \
    --samples 100 --seed 7 \
    --ins-start 1e-9 --ins-end 1e-17 --del-start 1e-9 --del-end 1e-17 \
    --sr-start 5 --sr-end 95 --sr-step 5 --k-start 0 --k-end 10 --k-step 1 \
    --printe-script ./PrinTE.sh \
    --ge 5400000 --st 54000 --mut 2e-9 --tstv 2.0 \
    --mxgs 211730673 --mngs 211530673 \
    --bed burnin.bed --fasta burnin.fasta \
    --clean-lib lib_clean.fa --ratios src/printe/data/ratios_ltr_only.tsv \
    --threads 250
```

`--mxgs`/`--mngs` bracket the target genome size in bytes (≈ the real genome's bp); the search
prioritizes parameter combinations whose terminal genome lands inside that band. Set `--ge`/`--st`
so the run spans your branch (a 5.4-My branch at 100-yr generations is `--ge 5400000`).

`init` writes the launcher into `slurm_grid/`. **Choose local or SLURM:**

- **Local (default).** `init` wrote `slurm_grid/submit_array.sh`. Launch this round, then hand the
  whole loop to the orchestrator (run both under `nohup` so they survive logout):

  ```bash
  nohup bash slurm_grid/submit_array.sh > submit.log 2>&1 &
  nohup python -m printe.grid.run_loop --target-hits 200 --max-rounds 500 \
        --samples 100 --explore-frac 0.2 > orchestrator.log 2>&1 &
  ```

  `run_loop.py` harvests results, proposes the next batch, and re-launches each round itself -
  run it from the directory containing `search_state.json` (written by `init`). It only drives
  the **local** launcher.

- **SLURM.** Add `--scheduler slurm --slurm-partition <p> --slurm-cpus 48 --slurm-mem 120G
  --slurm-time 24:00:00` to the `init` command; it writes `slurm_grid/submit_array.sbatch`, which
  you submit with `sbatch slurm_grid/submit_array.sbatch`. On SLURM, drive the rounds manually
  (`guided_search.py harvest` → `next` → `sbatch`) rather than with `run_loop.py`.

### Step 2 - re-annotate the simulated genomes

The grid runs skip post-processing for speed, so the simulated genomes are **not** dated
automatically. Annotate each one's LTR-RTs with the **same tool you used to benchmark**
([synLTR](#annotating-ltr-rts-and-building-a-species-specific-library)), producing in every
simulation directory a `<genome>_r1_ltr.tsv` table (**column 1 = `chrom:start-end`, column 7 =
p-distance** - exactly what the scorer reads). `scripts/discover_ltr.sh` is a template that loops
over every simulation directory and runs the annotation; edit its target FASTA name and paths to
match your run. Build the **reference** table the same way on your real genome.

### Step 3 - score and visualize

```bash
python -m printe.grid.build_composite_matrix \
    --ref-tsv real_genome_LTR.tsv --ref-fasta real_genome.fa \
    --gen-prefix gen5400000_final.fasta \
    --exp-tsv-name gen5400000_final_r1_ltr.tsv \
    --dist-metric rms --alphas 0 0 10 1 \
    --output composite_matrix.tsv --threads 64

# Best-fitting parameters = the row with the lowest Composite (column 20).
(head -1 composite_matrix.tsv; tail -n +2 composite_matrix.tsv | sort -t$'\t' -k20,20n) | head

# Visualize how each parameter affects the fit.
python -m printe.grid.contour_plot --input composite_matrix.tsv
```

- **`real_genome.fa` / `real_genome_LTR.tsv`** are your real target assembly and *its* LTR-RT
  table, produced by the **same** re-annotation command as the simulations (so the columns match).
- **`--gen-prefix`** must equal the FASTA filename inside each simulation directory;
  **`--exp-tsv-name`** is the basename of the per-simulation LTR TSV inside each directory. If your
  re-annotation wrote a different name (e.g. `<prefix>_kmer2ltr_dedup`), either rename it in each
  directory or pass that name to `--exp-tsv-name`.
- **Reading the winner:** the first four columns of the top row are your best-fit
  `insertion_rate`, `deletion_rate`, `solo_ratio`, `length_bias`. `--alphas` weights the four
  distance terms (count, length, genome size, age distribution); `0 0 10 1` weights genome size
  most heavily and age distribution next.

`contour_plot.py` writes pairwise contour maps of the Composite score and prints a sensitivity
report (parameter ranges, correlations, and confidence intervals for the optimum).

---
