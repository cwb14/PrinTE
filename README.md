# PrinTE: a forward simulator of transposable-element genome evolution

[![CI](https://github.com/cwb14/PrinTE/actions/workflows/ci.yml/badge.svg)](https://github.com/cwb14/PrinTE/actions/workflows/ci.yml)
[![Container](https://img.shields.io/badge/container-ghcr.io-blue)](https://github.com/cwb14/PrinTE/pkgs/container/printe)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.10-brightgreen)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

**PrinTE** evolves a genome one generation at a time. Starting from a real or synthetic
genome, it repeatedly **mutates** the DNA, **inserts** new transposable elements (TEs), and
**excises** old ones, writing a genome FASTA, an annotation BED, and an evolved TE library
for every generation you sample. It is built for plant-genome-scale data and runs on a laptop
or a cluster.

```
                 Present
Ancestor        ╭─▶▬▶─▶▬▶─▶▬◀─▶─▶▬▶─▶▬▶─▶▬▶  Expansion
─▶▬▶──▶▬◀──▶▬▶──┤ ⟶ forward time
                ╰─▶▬▶──▶▬▶──                 Contraction

▶▬▶ LTR-RT · ▶▬◀ DNA transposon (TIR) · ▶ solo-LTR · ── DNA
```

PrinTE is four tools in one:

1. **A forward TE simulator** that supports both **fixed** (constant) and **variable**
   (composition-dependent) rates of TE insertion and deletion, with realistic
   nesting, solo-LTR formation, cut-and-paste vs copy-and-paste transposition,
   selection against gene disruption, and chromatin bias.
2. **A grid search** that finds the evolutionary parameters which best reproduce a real
   genome - in an unbiased, reproducible way - rather than hand-tuning by eye.
3. **A benchmarking tool** because PrinTE knows exactly where every TE is, you can run your
   TE-annotation pipeline on a simulated genome and measure how well it recovers the truth.
4. **Utilities for two directions of time** reconstructing a clade's *past* and evolving it
   to the present, or projecting a present-day genome into the *future*.

> **PrinTE works from your TE annotations, so its conclusions are only as good as those
> annotations.** Before drawing biological inferences, **benchmark your annotation pipeline**
> once (see [Benchmark your TE annotations](docs/benchmarking.md)). This applies
> to **both** tracks below.

> **Scope note - PrinTE is LTR-RT-centric.** PrinTE simulates *all* TE classes, but its
> genome-size dynamics (solo-LTR formation, length-biased loss, LTR-RT dating) and its
> grid-search **scoring** are built around **LTR retrotransposons** since theyre the dominant driver of
> plant genome-size evolution. Non-LTR TEs are simulated and tracked, but the size signal you
> reconstruct and fit is an LTR-RT signal. The annotation, library-building, and ancestral-
> reconstruction steps below therefore focus on LTR-RTs (see
> [Annotating LTR-RTs](docs/ltr-annotation.md)).

---

## Install

PrinTE relies on **conda** and **mamba** - we recommend installing them via
[Miniforge](https://github.com/conda-forge/miniforge/releases) (pick the installer for your
OS and CPU). None of the routes below need root.

| | |
|---|---|
| **container** | `apptainer pull printe.sif docker://ghcr.io/cwb14/printe:latest`<br>**Usually the least trouble on a shared cluster** - Apptainer needs no daemon and no privileges. |
| **source** | `git clone https://github.com/cwb14/PrinTE && cd PrinTE`<br>`mamba env create -f environment.yml && conda activate PrinTE` |
| **pip** | `pip install printe` (code only - `make fetch-data` pulls the reference libraries) |

A **bioconda** recipe is in `conda/`, but it has not been submitted yet, so
`mamba install -c bioconda printe` will not work until it merges.

From a clone, `bash PrinTE.sh` and the installed `printe` command are the **same program**,
so every example below works either way.

Track 1 additionally needs R (`mamba env create -f environment-r.yml`) and an LTR-RT
discovery tool. Track 2 additionally needs
[TEgenomeSimulator](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator).

PrinTE clones [Kmer2LTR](https://github.com/cwb14/Kmer2LTR) into `~/.cache/printe/` the
first time it post-processes an evolve run, so that step needs network access once.
Pre-clone it there to run offline.

## Quick start

```bash
# 1) Build a burn-in genome and stop: 100 Mb, 1 chromosome, 4% CDS, 10% intact TE.
#    Writes burnin.fasta, burnin.bed, burnin.stat to the current directory.
printe --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10

# 2) Evolve for 30,000 generations, sampling every 10,000 (variable rates).
printe -sz 135Mb -cn 5 -P 25 -itp 21 -ir 1e-6 -dr 4e-6 -ge 30000 -st 10000 -t 20

# 3) Start from your OWN genome instead of a burn-in.
#    Supply -i so new insertions match your clade.
printe -f genome.fasta -b genome.bed -i my_TE.lib -ge 20000 -st 10000
```

Run `printe` with no arguments for a short menu, or `printe -h` for every option.

To run many parameter sets at once, or on a cluster or AWS, use the Nextflow pipeline:

```bash
nextflow run . -profile test,conda            # tiny end-to-end run on the bundled fixtures
nextflow run . -profile singularity,slurm --mode sweep \
    --fasta burnin.fasta --bed burnin.bed \
    --ref_tsv real_LTR.tsv --ref_fasta real_genome.fa
```

## Documentation

**Start with [Benchmark your annotations](docs/benchmarking.md).** Everything downstream
assumes you know how good your TE calls are, and that page is how you find out.

Then, depending on what you are doing:

| | |
|---|---|
| [Annotating LTR-RTs](docs/ltr-annotation.md) | Build a **species-specific** library so insertions resemble your own clade's elements |
| [Grid search](docs/grid-search.md) | Find the rates that reproduce a target genome - **this is what the paper used** |
| [Track 1 - past to present](docs/track1-past-to-present.md) | Reconstruct an ancestor and evolve it to the present |
| [Track 2 - present to future](docs/track2-present-to-future.md) | Project a genome forward and scan for stability |
| [Nextflow](docs/nextflow.md) and [AWS Batch](docs/aws-batch.md) | Running sweeps on a cluster or on AWS |
| [Custom TE libraries](docs/custom-te-library.md) | Assembling a library by hand |
| [Troubleshooting](docs/troubleshooting.md) | Read `pipeline.error` first. Also memory, offline runs, and the compiler |

## Two ways to use PrinTE

PrinTE answers questions in two directions of time. Pick the track that matches your data and
hypothesis.

| | **Past → present** | **Present → future** |
|---|---|---|
| **You start from** | a time-scaled tree + TE annotations for a clade | one present-day genome |
| **Build the burn-in with** | **PrinTE** (from a reconstructed ancestor) | **TEgenomeSimulator** (a digital replica) |
| **You evolve toward** | the **real** present-day genome | an unknown future |
| **You can score the result?** | **Yes** - against the real genome | **No** - the future hasn't happened |
| **Typical question** | Did TE activity **speed up or slow down** over time? | Under which parameters does the genome stay **stable**? |
| **Typical horizon** | the full branch length (millions of years) | short (~1 My) |

- **[Track 1 - Past → present](docs/track1-past-to-present.md)** reconstructs the genomic and TE state
  of a common ancestor, builds that ancestor as a burn-in, and evolves it forward until the
  simulation looks like the real genome. Because you have the real genome as ground truth,
  you can **score** each simulation and **test increasing- vs decreasing-rate hypotheses**.
- **[Track 2 - Present → future](docs/track2-present-to-future.md)** takes a digital replica of a real
  genome (built with [TEgenomeSimulator](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator)),
  feeds it to PrinTE, and scans for the parameter regimes that **maintain** the genome over a
  short horizon. There is no future genome to score against, so this is best framed as a
  **genome-stability** question.

Both tracks share the same engine ([the core simulator](#the-core-simulator)) and the same
optimizer ([grid search](docs/grid-search.md)).

---

## The core simulator

A PrinTE run has two phases:

- **Phase 1 - burn-in.** Build a starting genome: a synthetic genome seeded with genes and TEs.
  Skipped if you supply your own genome (`--fasta`/`--bed`) or resume a run (`--continue`).
- **Phase 2 - generation loop.** For each step, *mutate → insert → excise → rebuild the TE
  library*, then write `gen<N>_final.{fasta,bed,lib}`.

PrinTE models both **copy-and-paste** (Type 1) and **cut-and-paste** (Type 2) transposition.
Families flagged `cutpaste` in the `transposition` column of `ratios.tsv` are excluded from
replicative amplification, but real cut-and-paste transposons often do not precisely conserve copy 
number in a 1:1 way. PrinTE controls the balance between excision and successful re-insertion 
is set by `--cutpaste_reinsertion` (default 1.0: each excision is repaid by exactly one re-insertion, 
conserving copy number).  Set the ratio <1.0 for a net loss (failed re-insertion) or >1.0 for net 
amplification (insertions without excision). 

### Phase 1 - the burn-in genome

PrinTE's burn-in is flexible enough to mimic the composition of almost any real genome. It does
**not** create solo-LTRs in the burn-in - those arise during Phase 2.

**Inputs:**

1. **CDS FASTA (`--cds`)** - the **genes** inserted into the synthetic genome. The default is
   *Arabidopsis thaliana* TAIR10 CDS (19,621 sequences). To request more genes than that
   (`--cds_percent` / `--cds_num`), supply your own CDS FASTA. Each CDS is inserted at most once.

2. **TE library FASTA (`--TE_lib`)** - the **TEs** to insert. Headers must use **RepeatMasker
   format**:

   ```
   >[name]#[class]/[superfamily]
   >Os2670#MITE/Tourist
   ```

   The supported `#[class]/[superfamily]` suffixes are listed in `ratios.tsv`. TEs are sampled
   with replacement, so a single element can be inserted many times - **larger libraries give
   more diverse TE landscapes**. Control abundance with `--intact_TE_percent`/`--intact_TE_num`
   (intact) and `--frag_TE_percent`/`--frag_TE_num` (fragmented). In the burn-in, genes and TEs
   are laid out in tandem (no nesting) at random distances ≥ 20 bp apart.

3. **TE ratios file (`--TE_ratio`, default `ratios.tsv`)** - the relative frequency of each
   superfamily. Columns:

   1. `class`
   2. `superfamily`
   3. `weight` of **intact** TE
   4. `weight` of **fragmented** TE
   5. `transposition` (optional) - `cutpaste` for conserved Type 2 DNA transposons; blank or
      `copypaste` for Type 1 (the default).

   ```text
   # class   superfamily   intact   fragmented   transposition
   DNA       Helitron      0.07     0.07         copypaste
   MITE      DTH           0.07     0.07         cutpaste
   LTR       Gypsy         0.19     0.19         copypaste
   ```

   With the row above and `--intact_TE_percent 20 --frag_TE_percent 10`:
   intact `LTR/Gypsy` → `20 × 0.19 = 3.8%` of the genome; fragmented `LTR/Gypsy` →
   `10 × 0.19 = 1.9%`. A ready-made `ratios_ltr_only.tsv` is included for LTR-only studies.

**Key burn-in parameters:**

- **`--size` / `--chr_number`** - total genome size (must carry a `kb`/`Mb`/`Gb` suffix, e.g.
  `135Mb`) and chromosome count. Work is parallelized per chromosome, so more chromosomes
  generally run faster (but avoid oversplitting).
- **`--TE_mut_k` / `--TE_mut_Mmax`** - the age (substitution-divergence) distribution of TEs in
  the burn-in. TE mutations are clock-like, but old elements are preferentially lost, so the
  *observed* distribution skews young. `--TE_mut_k` sets the decay slope (larger → more young
  TEs); `--TE_mut_Mmax` sets the maximum divergence (use `--TE_mut_Mmax 0` to disable aging).
  PrinTE writes the chosen distribution to `burnin_mut_dist.pdf`.
- **`--mutation_bins`** - supply an explicit TE age distribution (a 3-column `start end frequency`
  file) instead of the decay model. **This is how the past→present track injects a reconstructed
  ancestral age distribution** (overrides `--TE_mut_k`/`--TE_mut_Mmax`).
- **`--burnin_only`** - build the burn-in and stop (no `-ge`/`-st` needed).

### Phase 2 - evolving the genome

PrinTE offers two rate models. Choose one.

#### Fixed rates

Constant TE insertion and deletion **per base of genome per generation**.

```bash
--fix 5e-9,1e-8     # insertion=5e-9, deletion=1e-8 per bp per generation
```

- **`--disable_genes`** - forbid TE insertion into genes (fixed mode only).

#### Variable rates

Insertion and deletion scale with the genome's current TE content. Rates are specified **per
intact TE per generation** and converted internally to per-base rates; accumulated TE mutations
are inherited by new copies.

- **`--insert_rate` / `--delete_rate`** - insertions/deletions per intact TE per generation.
- **`--birth_rate`** - occasionally resamples from the **original** library (`--TE_lib`),
  modeling horizontal transfer or revival of an extinct lineage. Without it, only existing
  intact TEs propagate.

**Selective constraint** (variable mode). PrinTE treats genes as selectively constrained and
TEs + intergenic DNA as neutral.

- **`--sigma`** - how unevenly constraint is spread across genes (log-normal). Low → most genes
  similarly constrained; high → a few highly constrained, the rest mild. Visualized in
  `lognormal_distribution.pdf`.
- **`--sel_coeff`** - by default PrinTE does not purge TEs by fitness; set this to make more
  deleterious (gene-disrupting) insertions more likely to be deleted.
- **`--promoter-boundary`** - bp up/downstream of a gene where an overlapping TE still counts as
  disrupting that gene.

**Chromatin bias.** TEs may insert/delete more readily in euchromatin. PrinTE approximates
euchromatin as the neighborhood of genes.

- **`--chromatin_buffer`** - bp around genes treated as euchromatin.
- **`--chromatin_bias_insert` / `--chromatin_bias_delete`** - euchromatin bias for insertion / deletion.

**LTR / length dynamics** (both modes).

- **`--solo_rate`** - when an intact LTR-RT is removed, the % chance it leaves a **solo-LTR**
  rather than disappearing entirely.
- **`--k`** - length bias of deletion: longer elements (more prone to unequal recombination) are
  removed faster. `--k 0` disables it. See `weighted_candidate_selection.pdf`.

In short: use **fixed** for simple constant-rate experiments, and **variable** for dynamics
shaped by genome composition, selection, and chromatin.

### General options

- **`--generation_end` / `--step`** - total generations and generations per step.
  `--generation_end 3000 --step 1000` samples generations 1000, 2000, 3000.
- **`--mutation_rate` / `--TsTv`** - DNA substitution rate per bp per generation, and the
  transition/transversion ratio.
- **`--max_size` / `--min_size`** - stop Phase 2 early once a regression projection of the
  terminal genome size is *confidently* above (`--max_size`) or below (`--min_size`) a bound
  (e.g. `1G`, `100M`). Handy for parameter scans that would otherwise blow up or collapse.
- **`--fasta` / `--bed`** - start from your own genome + PrinTE-format BED (skips the burn-in;
  see [Bringing your own genome](#bringing-your-own-genome)).
- **`--continue`** - resume from the last completed generation in the working directory. Combine
  with a larger `--generation_end` (and new rates) to **extend** a finished run.
- **`--model`** - substitution model for LTR-RT dating (`raw`, `K2P`, `JC69`; default `K2P`).
- **`--pergen_select`** - how many generations get the (slower) per-generation LTR-RT dating
  analysis (default 2 = first + last).
- **`--ex_LTR`** - drop library LTR-RTs that lack a detectable LTR.
- **`--no_postproc` / `--keep_temps`** - skip all plots/reports / keep intermediate files.
- **`--seed` / `--threads`** - reproducibility and parallelism.

### Outputs

For each sampled generation, the primary outputs are the **genome**, its **annotation**, and
the **evolved TE library**:

```bash
ls gen40000_final.*
gen40000_final.fasta   # genome
gen40000_final.bed     # gene + TE coordinates (PrinTE BED, see below)
gen40000_final.lib     # FASTA of intact TEs eligible to transpose next generation
```

The `.lib` is a FASTA whose headers are the bare feature IDs of the surviving intact TEs (e.g.
`>tuteh_AC183372_584#LTR/unknown~LTRlen:126`). You can inspect it (`grep '>' gen40000_final.lib`)
or feed it to a later run via `--clean_lib`.

The **BED** has six tab-separated columns - `chromosome`, `start`, `end`, `feature_ID`, `TSD`,
`strand` (coordinates are 0-based, half-open):

```text
chr1  1857719  1859314  tuteh_AC183372_584#LTR/unknown~LTRlen:126;CUT_BY:Os2721#DNAnona/hAT   CATTC  +
chr1  1859314  1859734  Os2721#DNAnona/hAT;NESTED_IN:tuteh_AC183372_584#LTR/unknown~LTRlen:126 CTTCCG +
chr1  1862094  1863435  TE_00016444_FRAG#MITE/DTA                                              NA     -
chr1  1869245  1869502  anysaf_AC211487_11211#LTR/Ty3~LTRlen:257_SOLO                          GCGCG  -
```

- `Os2721#DNAnona/hAT` is **nested inside** `tuteh_AC183372_584#LTR/unknown` (which is therefore
  tagged `CUT_BY`).
- `_FRAG` marks a fragmented element, `_SOLO` a solo-LTR. Neither is intact, so neither can
  transpose, but both remain available for deletion.

Other useful files:

- **`gen<N>_mut.txt`** - the per-generation mutation report (mutations introduced this step,
  and the accumulated substitution distance from the burn-in genome).
- **`burnin.stat`** - composition of the burn-in (bp, % genes, % intact/fragmented TE).
- **`pipeline.report`** - insertions (nested/non-nested) and deletions per generation.
- **`pipeline.log` / `pipeline.error`** - full run log and error log. Both are always created;
  on success `pipeline.error` holds only a start banner, so **if a run dies, read
  `pipeline.error` for the traceback.**
- **Figures** - `percent_TE.pdf` (TE fraction over time), `solo_intact.pdf` (solo:intact ratio),
  `stat_intact_plot.pdf`/`stat_frag_plot.pdf` (superfamily counts), `genome_size_plot.pdf`,
  `all_LTR_density.pdf` (LTR-RT age distribution across generations).

If `--birth_rate` was used, count horizontally acquired TEs with:

```bash
grep 'Number of born TEs to insert' pipeline.log
```

### Bringing your own genome

To skip the burn-in (`--fasta genome.fasta --bed genome.bed`), your BED must be in the same
**PrinTE BED** format PrinTE writes (shown above): 6 tab-separated columns
(`chromosome  start  end  feature_ID  TSD  strand`), 0-based half-open coordinates, where:

- `feature_ID` ends in `#class/superfamily` matching a row in `ratios.tsv` (genes are named
  `gene1`, `gene2`, … with no `#`);
- **LTR-RTs must carry `~LTRlen:<N>`** (the LTR length in bp) or they will not be treated as
  intact LTR elements;
- `TSD` is the target-site-duplication sequence, or `NA` if none;
- the easiest non-nested starting line looks like:
  `chr1  20000  25000  myCopia_1#LTR/Copia~LTRlen:430  ATGCA  +`.

The two supported ways to obtain a valid BED without hand-editing are: **(a)** let PrinTE build a
burn-in (Track 1), or **(b)** import a TEgenomeSimulator genome with `util/gff_to_bed.py`
(Track 2). If you omit `-i/--TE_lib` here, new insertions are drawn from the bundled default
library - supply `-i your_TE.lib` (RepeatMasker headers) so they match your organism.

---
## Requirements & supported systems

Supported on **Linux** and **macOS**. Tested on macOS 14.6.1 (Apple M2 Pro, 16 GB RAM) and
Ubuntu 22.04.5 LTS (AMD EPYC 7763). Requires ~2 GB RAM minimum; 16 GB+, 4+ cores, and a network
connection (to clone Kmer2LTR) are recommended. The bundled `ltr_mutator` binary is rebuilt
automatically from source (`g++ -std=c++17 -O3 -fopenmp`) if the precompiled one won't run.
