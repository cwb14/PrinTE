## The Nextflow pipeline

Runs PrinTE on a laptop, a Slurm cluster, or AWS Batch from the same command, and
scatters a parameter sweep across whatever the executor gives it.

[Back to the README](../README.md)

---

### Try it

You need Nextflow itself first - `mamba install -y -c bioconda nextflow`. Then, from any
directory:

```bash
nextflow run cwb14/PrinTE -profile test,conda
```

That pulls the pipeline from GitHub, so you do not need to clone anything. It runs the
**simulate** path end to end on the bundled test fixtures in a couple of minutes: burn-in,
then two generations. CI runs the same thing on every push as `-profile test,local`. The
**sweep** path is not covered by CI, so run it once yourself before you trust it with a
large grid.

The examples below use `nextflow run .`, which works when you are inside a clone. Replace
the `.` with `cwb14/PrinTE` to run from anywhere without cloning.

### The two modes

**`--mode simulate`** (the default) runs one parameter set: burn-in, then the generation
loop, then post-processing if you ask for it.

```bash
nextflow run . -profile conda \
    --te_lib my_TE.lib --size 135Mb --chr_number 5 \
    --intact_te_percent 21 --generation_end 30000 --step 10000 \
    --insert_rate 1e-6 --delete_rate 4e-6 --outdir results
```

**`--mode sweep`** builds the parameter grid, runs one simulation per combination, scores
each against a reference genome, and gathers the scores into one table.

```
BUILD START ─┐
             ├─▶ SIMULATE (× N combos) ─▶ SCORE (× N) ─▶ AGGREGATE ─▶ CONTOUR
EMIT_GRID  ──┘
```

```bash
nextflow run . -profile singularity,slurm --mode sweep \
    --fasta burnin.fasta --bed burnin.bed --clean_lib lib_clean.fa \
    --ins-start 1e-9 --ins-end 1e-11 --ins-count 5 \
    --del-start 1e-9 --del-end 1e-11 --del-count 5 \
    --ref_tsv real_genome_LTR.tsv --ref_fasta real_genome.fa \
    --outdir sweep_results
```

Each simulation lands in `results/simulations/insertion_rates_<i>_deletion_rates_<d>_solo_ratio_<s>_length_bias_<k>/`,
the same directory naming the Python front-ends use, so the two are interchangeable
downstream.

**Scoring needs an LTR-RT table per simulation.** The sweep skips post-processing for
speed, so nothing dates the simulated LTR-RTs for you. Annotate each simulation with the
same tool you benchmarked with (see [Annotating LTR-RTs](ltr-annotation.md)) and leave
the table beside the final FASTA before running `SCORE`.

### Profiles

Combine one environment profile with one executor profile.

| profile | what it does |
|---|---|
| `conda` | resolves `environment.yml` per process |
| `docker` | uses `ghcr.io/cwb14/printe`, running as your uid so bind-mounted output is not root-owned |
| `singularity` | same image through Apptainer or Singularity, with automounts |
| `local` | runs on the current machine |
| `slurm` | submits jobs, rate-limited to 10/min |
| `awsbatch` | see [AWS Batch](aws-batch.md) |
| `test` | tiny fixtures, capped at 2 cpus |

Slurm needs a partition, and possibly an account, which are site-specific and so are not
baked in:

```bash
nextflow run . -profile singularity,slurm \
    -process.queue mypartition -process.clusterOptions '-A myaccount'
```

### Resources and retries

Processes carry `process_low`, `process_medium`, or `process_high` labels; the ceilings
live in `nf/conf/base.config`, and memory doubles on a retry. See
[Troubleshooting](troubleshooting.md#how-much-memory-do-i-need) for what to budget.

`-resume` works as usual. Because `SIMULATE` is keyed on the parameter combination,
adding combinations to a sweep re-runs only the new ones.

Each run writes `pipeline_info/` with a timeline, an execution report, a trace, and the
DAG.

### What this does not replace

`guided_search.py` stays the tool for active-learning search. Its loop is iterative -
harvest the finished runs, fit a random-forest surrogate, propose the next batch, relaunch
- and that does not map onto a DAG without leaning on preview features for no benefit to
anyone. Use the Nextflow sweep when you want to cover a fixed grid and the portability
that comes with it; use `guided_search.py` when you want the surrogate to spend your
compute for you. See [Grid search](grid-search.md).
