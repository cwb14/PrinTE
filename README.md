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

Pick **one** of the three routes below and paste the whole block. Each one ends with the
same check, so you know it worked before you go any further. None of them needs root, and
after any of them the `printe` command works the same way.

### 1. Container - least trouble on a shared cluster

Nothing to compile and nothing to resolve. Most HPC systems already have Apptainer
(check with `apptainer --version`); if yours does not, the first line installs it into
your own space.

```bash
mamba install -y -c conda-forge apptainer     # skip if `apptainer --version` already works

apptainer pull printe.sif docker://ghcr.io/cwb14/printe:latest

./printe.sif --version
```

You should see `PrinTE 1.0.0`. **The image is the program** - run it as `./printe.sif`
followed by the usual options. Nothing was installed anywhere else, `ls` shows you the
whole of PrinTE, and `rm printe.sif` removes it. Apptainer makes your current directory
visible inside the image, so results land where you ran the command.

In the examples below, wherever you see `printe`, type `./printe.sif` instead (or the full
path to it, if you are working in a different directory).

<details>
<summary><b>Optional:</b> type <code>printe</code> from anywhere instead of <code>./printe.sif</code></summary>

This is the only step in these instructions that changes anything outside the current
folder. It writes one small launcher and adds one line to `~/.bashrc`. Re-running it is
harmless, and the last line tells you how to undo it.

```bash
mkdir -p ~/.local/bin
printf '#!/bin/bash\nexec "%s/printe.sif" "$@"\n' "$PWD" > ~/.local/bin/printe
chmod +x ~/.local/bin/printe
grep -qxF 'export PATH="$HOME/.local/bin:$PATH"' ~/.bashrc \
  || echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
export PATH="$HOME/.local/bin:$PATH"

printe --version        # to undo later:  rm ~/.local/bin/printe
```

The launcher points at that exact `printe.sif`, so do not move the image afterwards. If
you later install PrinTE from source as well, delete this launcher first - it sits ahead
of your conda environment on `PATH` and would shadow the other copy.

</details>

### 2. From source - if you want to read or change the code

```bash
git clone https://github.com/cwb14/PrinTE.git
cd PrinTE
mamba env create -f environment.yml
conda activate PrinTE
pip install -e .

printe --version
```

You should see `PrinTE 1.0.0`. The `pip install -e .` line is the one that puts `printe`
on your PATH - without it the quick start below will not run. You can also call
`bash PrinTE.sh` from inside the clone instead; it is the same program.

Remember to `conda activate PrinTE` in every new terminal.

If `printe --version` fails with `apptainer: not found`, you have a leftover wrapper from
route 1 shadowing this install. Run `command -v printe` to confirm, then
`rm ~/.local/bin/printe`.

### 3. Nextflow - if you want to run many parameter sets, or use a cluster

```bash
mamba install -y -c bioconda nextflow

nextflow run cwb14/PrinTE -profile test,conda
```

That downloads the pipeline from GitHub and runs a small test simulation, so you do not
need to clone anything. See [the Nextflow pipeline](docs/nextflow.md) for real runs.

### Not available yet

`pip install printe` and `mamba install -c bioconda printe` do **not** work yet - the
package is not on PyPI and the bioconda recipe in `conda/` has not been submitted. Use one
of the three routes above.

### Extras, only if you need them

Track 1 (past to present) also needs R: `mamba env create -f environment-r.yml`, plus an
LTR-RT discovery tool. Track 2 also needs
[TEgenomeSimulator](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator).

The first time PrinTE post-processes an evolve run it clones
[Kmer2LTR](https://github.com/cwb14/Kmer2LTR) into `~/.cache/printe/` to date the LTR-RTs,
so that one step needs internet. A burn-in never does, and `--no_postproc` skips it.

## Quick start

These work whichever route you installed with. Run them in a directory you can write to.

```bash
# 1) Build a starting genome and stop: 100 Mb, 1 chromosome, 4% genes, 10% intact TE.
#    Takes a few minutes. Writes burnin.fasta, burnin.bed and burnin.stat here.
printe --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10

# 2) Evolve a genome for 30,000 generations, saving it every 10,000.
printe -sz 135Mb -cn 5 -P 25 -itp 21 -ir 1e-6 -dr 4e-6 -ge 30000 -st 10000 -t 20

# 3) Start from YOUR genome instead of a simulated one. The BED must be in PrinTE's
#    format (see docs/file-formats.md), and -i should be your clade's TE library.
printe -f genome.fasta -b genome.bed -i my_TE.lib -ge 20000 -st 10000
```

Run `printe` on its own for a short menu, or `printe -h` for every option.

**Your first run will sit for several minutes on `=== TE Library Processing ===`.** That is
normal. PrinTE is comparing every TE in the library against itself to work out how old each
one is, and it only has to do that once. It writes the result as `lib_clean.fa` in your
working directory, so pass that back on every later run and the step is skipped entirely:

```bash
printe --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10 -cl lib_clean.fa
```

**If something fails,** look at `pipeline.error` in the directory you ran from. On a
successful run it contains only a start banner, so anything else in it is the real problem.
[Troubleshooting](docs/troubleshooting.md) covers the common ones.

## What comes with PrinTE

You do not have to supply anything to try it - these ship with the tool and are used by
default. Look at them before you build your own, so you can see the formats expected.

| file | what it is |
|---|---|
| `ratios.tsv` | how often each TE superfamily is inserted, and which ones cut-and-paste. The default `-r`. |
| `ratios_ltr_only.tsv` | the same, restricted to LTR retrotransposons |
| `maize_rice_arab_curated_TE.lib.gz` | the TE library new insertions are drawn from. The default `-i`. |
| `TAIR10.cds.fa.gz` | *Arabidopsis* CDS, used as the "genes" in a simulated genome. The default `-c`. |
| `TAIR10.pep.fa.gz` | proteins, for classifying LTR-RTs during annotation |
| `ltr-db.fa.gz` | a larger LTR-RT exemplar database. Not bundled - `make fetch-data` downloads it. |

To see where they are on your system:

```bash
python -m printe.paths                       # source or pip install
apptainer exec printe.sif python -m printe.paths     # if you installed the container
```

That prints one line per file. Add a filename to get just that path, which you can use
straight away:

```bash
less "$(python -m printe.paths ratios.tsv)"
zcat "$(python -m printe.paths TAIR10.cds.fa.gz)" | head
```

**To edit one**, copy it out first, change your copy, and pass it back:

```bash
cp "$(python -m printe.paths ratios.tsv)" my_ratios.tsv
# edit my_ratios.tsv, then
printe --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10 -r my_ratios.tsv
```

The same applies to the TE library and CDS: `-i my_TE.lib` and `-c my_cds.fa`. For a
species-specific library, see [Annotating LTR-RTs](docs/ltr-annotation.md).

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
