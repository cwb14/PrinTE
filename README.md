# PrinTE: a forward simulator of transposable-element genome evolution

**PrinTE** evolves a genome one generation at a time. Starting from a real or synthetic
genome, it repeatedly **mutates** the DNA, **inserts** new transposable elements (TEs), and
**excises** old ones, writing a genome FASTA, an annotation BED, and an evolved TE library
for every generation you sample. It is built for plant-genome-scale data and runs on a laptop
or a cluster.

PrinTE is four tools in one:

1. **A forward TE simulator** that supports both **fixed** (constant) and **variable**
   (composition-dependent) rates of TE insertion and deletion, with realistic
   nesting, solo-LTR formation, cut-and-paste vs copy-and-paste transposition,
   selection against gene disruption, and chromatin bias.
2. **A grid search** that finds the evolutionary parameters which best reproduce a real
   genome — in an unbiased, reproducible way — rather than hand-tuning by eye.
3. **A benchmarking tool**: because PrinTE knows exactly where every TE is, you can run your
   TE-annotation pipeline on a simulated genome and measure how well it recovers the truth.
4. **Utilities for two directions of time**: reconstructing a clade's *past* and evolving it
   to the present, or projecting a present-day genome into the *future*.

> **PrinTE works from your TE annotations, so its conclusions are only as good as those
> annotations.** Before drawing biological inferences, **benchmark your annotation pipeline**
> once (see [Benchmark your TE annotations](#benchmark-your-te-annotations-first)). This applies
> to **both** tracks below.

> **Scope note — PrinTE is LTR-RT-centric.** PrinTE simulates *all* TE classes, but its
> genome-size dynamics (solo-LTR formation, length-biased loss, LTR-RT dating) and its
> grid-search **scoring** are built around **LTR retrotransposons** — the dominant driver of
> plant genome-size evolution. Non-LTR TEs are simulated and tracked, but the size signal you
> reconstruct and fit is an LTR-RT signal. The annotation, library-building, and ancestral-
> reconstruction steps below therefore focus on LTR-RTs (see
> [Annotating LTR-RTs](#annotating-ltr-rts-and-building-a-species-specific-library)).

---

## Two ways to use PrinTE

PrinTE answers questions in two directions of time. Pick the track that matches your data and
hypothesis.

| | **Past → present** | **Present → future** |
|---|---|---|
| **You start from** | a time-scaled tree + TE annotations for a clade | one present-day genome |
| **Build the burn-in with** | **PrinTE** (from a reconstructed ancestor) | **TEgenomeSimulator** (a digital replica) |
| **You evolve toward** | the **real** present-day genome | an unknown future |
| **You can score the result?** | **Yes** — against the real genome | **No** — the future hasn't happened |
| **Typical question** | Did TE activity **speed up or slow down** over time? | Under which parameters does the genome stay **stable**? |
| **Typical horizon** | the full branch length (millions of years) | short (~1 My) |

- **[Track 1 — Past → present](#track-1--past--present)** reconstructs the genomic and TE state
  of a common ancestor, builds that ancestor as a burn-in, and evolves it forward until the
  simulation looks like the real genome. Because you have the real genome as ground truth,
  you can **score** each simulation and **test increasing- vs decreasing-rate hypotheses**.
- **[Track 2 — Present → future](#track-2--present--future)** takes a digital replica of a real
  genome (built with [TEgenomeSimulator](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator)),
  feeds it to PrinTE, and scans for the parameter regimes that **maintain** the genome over a
  short horizon. There is no future genome to score against, so this is best framed as a
  **genome-stability** question.

Both tracks share the same engine ([the core simulator](#the-core-simulator)) and the same
optimizer ([grid search](#grid-search-finding-parameters-without-bias)).

---

## Installation

PrinTE relies on **conda** and **mamba**. We recommend installing them via
[Miniforge](https://github.com/conda-forge/miniforge/releases) (pick the installer for your OS
and CPU). Then:

```bash
git clone https://github.com/cwb14/PrinTE.git
mamba env create -f PrinTE/env.yml
conda activate PrinTE
```

This installs everything the **core simulator** and **grid search** need. PrinTE additionally
clones [Kmer2LTR](https://github.com/cwb14/Kmer2LTR.git) (used to date LTR-RTs) into the PrinTE
directory the first time it runs the post-processing of an *evolve* run — not during a burn-in.
If you work offline, pre-clone it: `git clone https://github.com/cwb14/Kmer2LTR.git` inside the
PrinTE directory.

**Track 1 (past → present)** additionally needs **R** with the packages `ape`, `phytools`, and
`optparse` for ancestral reconstruction (install these yourself, e.g. `mamba install -c
conda-forge r-ape r-phytools r-optparse`), and an **LTR-RT discovery tool** of your choice
(e.g. LTRharvest/LTR_retriever) to annotate genomes. **Track 2** additionally needs
[TEgenomeSimulator](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator).

---

## Quick start

```bash
# 1) Build a burn-in genome only: 100 Mb, 1 chromosome, 4% CDS, 10% intact TE.
#    Writes burnin.fasta, burnin.bed, burnin.stat, burnin_mut_dist.pdf to the current directory.
bash PrinTE/PrinTE.sh --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10

# 2) Evolve that genome for 30,000 generations, sampling every 10,000 (variable rates).
bash PrinTE/PrinTE.sh -sz 135Mb -cn 5 -P 25 -itp 21 -ir 1e-6 -dr 4e-6 -ge 30000 -st 10000 -t 20

# 3) Start from your OWN genome instead of a burn-in (see "Bringing your own genome" below).
#    Supply -i so new insertions match your clade; otherwise the bundled default library is used.
bash PrinTE/PrinTE.sh -f genome.fasta -b genome.bed -i my_TE.lib -ge 20000 -st 10000
```

Run `bash PrinTE/PrinTE.sh` for a concise menu, or `bash PrinTE/PrinTE.sh -h` for the full
option list.

---

## The core simulator

A PrinTE run has two phases:

- **Phase 1 — burn-in.** Build a starting genome: a synthetic genome seeded with genes and TEs.
  Skipped if you supply your own genome (`--fasta`/`--bed`) or resume a run (`--continue`).
- **Phase 2 — generation loop.** For each step, *mutate → insert → excise → rebuild the TE
  library*, then write `gen<N>_final.{fasta,bed,lib}`.

PrinTE models both **copy-and-paste** (Type 1) and **cut-and-paste** (Type 2) transposition.
Families flagged `cutpaste` in the `transposition` column of `ratios.tsv` are conserved: they
are excluded from replicative amplification, and when an *intact* element is excised it
relocates (an element of the same superfamily is re-inserted the next generation, keeping copy
number constant); a *fragmented* copy (`_FRAG`, `_SOLO`, or one cut by a nested insertion) is
simply lost. Every other family is copy-and-paste and can amplify. Older 4-column `ratios.tsv`
files (no `transposition` column) behave as all copy-and-paste, exactly as before.

### Phase 1 — the burn-in genome

PrinTE's burn-in is flexible enough to mimic the composition of almost any real genome. It does
**not** create solo-LTRs in the burn-in — those arise during Phase 2.

**Inputs:**

1. **CDS FASTA (`--cds`)** — the **genes** inserted into the synthetic genome. The default is
   *Arabidopsis thaliana* TAIR10 CDS (19,621 sequences). To request more genes than that
   (`--cds_percent` / `--cds_num`), supply your own CDS FASTA. Each CDS is inserted at most once.

2. **TE library FASTA (`--TE_lib`)** — the **TEs** to insert. Headers must use **RepeatMasker
   format**:

   ```
   >[name]#[class]/[superfamily]
   >Os2670#MITE/Tourist
   ```

   The supported `#[class]/[superfamily]` suffixes are listed in `ratios.tsv`. TEs are sampled
   with replacement, so a single element can be inserted many times — **larger libraries give
   more diverse TE landscapes**. Control abundance with `--intact_TE_percent`/`--intact_TE_num`
   (intact) and `--frag_TE_percent`/`--frag_TE_num` (fragmented). In the burn-in, genes and TEs
   are laid out in tandem (no nesting) at random distances ≥ 20 bp apart.

3. **TE ratios file (`--TE_ratio`, default `ratios.tsv`)** — the relative frequency of each
   superfamily. Columns:

   1. `class`
   2. `superfamily`
   3. `weight` of **intact** TE
   4. `weight` of **fragmented** TE
   5. `transposition` (optional) — `cutpaste` for conserved Type 2 DNA transposons; blank or
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

- **`--size` / `--chr_number`** — total genome size (must carry a `kb`/`Mb`/`Gb` suffix, e.g.
  `135Mb`) and chromosome count. Work is parallelized per chromosome, so more chromosomes
  generally run faster (but avoid oversplitting).
- **`--TE_mut_k` / `--TE_mut_Mmax`** — the age (substitution-divergence) distribution of TEs in
  the burn-in. TE mutations are clock-like, but old elements are preferentially lost, so the
  *observed* distribution skews young. `--TE_mut_k` sets the decay slope (larger → more young
  TEs); `--TE_mut_Mmax` sets the maximum divergence (use `--TE_mut_Mmax 0` to disable aging).
  PrinTE writes the chosen distribution to `burnin_mut_dist.pdf`.
- **`--mutation_bins`** — supply an explicit TE age distribution (a 3-column `start end frequency`
  file) instead of the decay model. **This is how the past→present track injects a reconstructed
  ancestral age distribution** (overrides `--TE_mut_k`/`--TE_mut_Mmax`).
- **`--burnin_only`** — build the burn-in and stop (no `-ge`/`-st` needed).

### Phase 2 — evolving the genome

PrinTE offers two rate models. Choose one.

#### Fixed rates

Constant TE insertion and deletion **per base of genome per generation**.

```bash
--fix 5e-9,1e-8     # insertion=5e-9, deletion=1e-8 per bp per generation
```

- **`--disable_genes`** — forbid TE insertion into genes (fixed mode only).

#### Variable rates

Insertion and deletion scale with the genome's current TE content. Rates are specified **per
intact TE per generation** and converted internally to per-base rates; accumulated TE mutations
are inherited by new copies.

- **`--insert_rate` / `--delete_rate`** — insertions/deletions per intact TE per generation.
- **`--birth_rate`** — occasionally resamples from the **original** library (`--TE_lib`),
  modeling horizontal transfer or revival of an extinct lineage. Without it, only existing
  intact TEs propagate.

**Selective constraint** (variable mode). PrinTE treats genes as selectively constrained and
TEs + intergenic DNA as neutral.

- **`--sigma`** — how unevenly constraint is spread across genes (log-normal). Low → most genes
  similarly constrained; high → a few highly constrained, the rest mild. Visualized in
  `lognormal_distribution.pdf`.
- **`--sel_coeff`** — by default PrinTE does not purge TEs by fitness; set this to make more
  deleterious (gene-disrupting) insertions more likely to be deleted.
- **`--promoter-boundary`** — bp up/downstream of a gene where an overlapping TE still counts as
  disrupting that gene.

**Chromatin bias.** TEs may insert/delete more readily in euchromatin. PrinTE approximates
euchromatin as the neighborhood of genes.

- **`--chromatin_buffer`** — bp around genes treated as euchromatin.
- **`--chromatin_bias_insert` / `--chromatin_bias_delete`** — euchromatin bias for insertion / deletion.

**LTR / length dynamics** (both modes).

- **`--solo_rate`** — when an intact LTR-RT is removed, the % chance it leaves a **solo-LTR**
  rather than disappearing entirely.
- **`--k`** — length bias of deletion: longer elements (more prone to unequal recombination) are
  removed faster. `--k 0` disables it. See `weighted_candidate_selection.pdf`.

In short: use **fixed** for simple constant-rate experiments, and **variable** for dynamics
shaped by genome composition, selection, and chromatin.

### General options

- **`--generation_end` / `--step`** — total generations and generations per step.
  `--generation_end 3000 --step 1000` samples generations 1000, 2000, 3000.
- **`--mutation_rate` / `--TsTv`** — DNA substitution rate per bp per generation, and the
  transition/transversion ratio.
- **`--max_size` / `--min_size`** — stop Phase 2 early once a regression projection of the
  terminal genome size is *confidently* above (`--max_size`) or below (`--min_size`) a bound
  (e.g. `1G`, `100M`). Handy for parameter scans that would otherwise blow up or collapse.
- **`--fasta` / `--bed`** — start from your own genome + PrinTE-format BED (skips the burn-in;
  see [Bringing your own genome](#bringing-your-own-genome)).
- **`--continue`** — resume from the last completed generation in the working directory. Combine
  with a larger `--generation_end` (and new rates) to **extend** a finished run.
- **`--model`** — substitution model for LTR-RT dating (`raw`, `K2P`, `JC69`; default `K2P`).
- **`--pergen_select`** — how many generations get the (slower) per-generation LTR-RT dating
  analysis (default 2 = first + last).
- **`--ex_LTR`** — drop library LTR-RTs that lack a detectable LTR.
- **`--no_postproc` / `--keep_temps`** — skip all plots/reports / keep intermediate files.
- **`--seed` / `--threads`** — reproducibility and parallelism.

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

The **BED** has six tab-separated columns — `chromosome`, `start`, `end`, `feature_ID`, `TSD`,
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

- **`gen<N>_mut.txt`** — the per-generation mutation report (mutations introduced this step,
  and the accumulated substitution distance from the burn-in genome).
- **`burnin.stat`** — composition of the burn-in (bp, % genes, % intact/fragmented TE).
- **`pipeline.report`** — insertions (nested/non-nested) and deletions per generation.
- **`pipeline.log` / `pipeline.error`** — full run log and error log. Both are always created;
  on success `pipeline.error` holds only a start banner, so **if a run dies, read
  `pipeline.error` for the traceback.**
- **Figures** — `percent_TE.pdf` (TE fraction over time), `solo_intact.pdf` (solo:intact ratio),
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
library — supply `-i your_TE.lib` (RepeatMasker headers) so they match your organism.

---

## Benchmark your TE annotations first

PrinTE simulates from TE annotations, so **the forward framework assumes high-quality
annotations**. Conveniently, PrinTE is also the tool to *test* that assumption: it knows the
exact location of every TE it places, so a simulated genome is a ground-truth benchmark for any
annotation pipeline. **Do this once before either track.**

1. **Build a genome** whose TE landscape mimics your real one. `--burnin_only` writes
   `burnin.fasta` (genome) and `burnin.bed` (the ground truth) to the current directory:

   ```bash
   bash PrinTE/PrinTE.sh --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10
   ```

2. **Annotate the LTR-RTs** in `burnin.fasta`. We recommend
   [synLTR](#annotating-ltr-rts-and-building-a-species-specific-library) (below): with
   `--out_prefix burnin` it writes `burnin_r1_ltr.tsv`, whose column 1 is the `chrom:start-end`
   locus — exactly what `-pass_scn` reads. Any LTR finder works as long as you hand `bedtools.py`
   an SCN or PASS file (an SCN line has ≥12 whitespace columns with start/end in columns 1–2 and
   the sequence name last; a PASS line starts `chr:start..end`). Use the **same pipeline on your
   real genomes** in Track 1/2.

3. **Compare predictions to the truth** with `util/bedtools.py`:

   ```bash
   python PrinTE/util/bedtools.py \
       -pass_scn burnin_r1_ltr.tsv \
       -bed burnin.bed --printe-intact \
       -r 0.9
   ```

   `--printe-intact` subsets the truth BED to intact LTR-RTs so false negatives stay meaningful;
   `-r 0.9` requires reciprocal 90% overlap (each call and the true element cover ≥90% of the
   other — lower it, e.g. `-r 0.5`, for looser positional agreement). Expected output:

   ```text
   Subset BED via --printe-intact: 412 intact LTR-RT entries
   Overlapping entries: 389 (388 unique)
   Entries unique to SCN/PASS file: 14
   Entries unique to BED file: 24

   Precision = 0.9653
   Recall    = 0.9417        # = sensitivity = detection rate = TP / (TP+FN)
   F1        = 0.9534
   FDR       = 0.0347
   ```

As a rule of thumb, **recall and F1 ≥ 0.90** at `-r 0.9` indicate a pipeline you can trust;
below ~0.8, fix the annotation (or simplify the simulation target) before interpreting any
evolutionary result.

---

## Annotating LTR-RTs and building a species-specific library

Two LTR-RT operations recur in PrinTE: **annotating** LTR-RTs in a genome (to benchmark your
calls and to score simulations), and **building a species-specific LTR-RT library** so that
simulated insertions resemble the real elements of your genome. The companion tools
[synLTR](https://github.com/cwb14/synLTR) (annotation) and
[Kmer2LTR](https://github.com/cwb14/Kmer2LTR) (library building + LTR dating) do this and benchmark well for PrinTE's purposes.

```bash
git clone https://github.com/cwb14/synLTR.git
git clone https://github.com/cwb14/Kmer2LTR.git    # also auto-cloned by PrinTE/synLTR when needed
```

**Annotate LTR-RTs in a genome** (real or simulated). `--max-rounds 1` does a single, non-nested
pass — the right setting for PrinTE; `--proteins` is a protein FASTA used to classify elements
(the bundled `PrinTE/data/TAIR10.pep.fa.gz` works for most plant genomes):

```bash
bash synLTR/module2/ltrharvest_wrapper2.sh \
    --genome mygenome.fa --proteins PrinTE/data/TAIR10.pep.fa.gz \
    --max-rounds 1 --threads 64 --out_prefix mygenome
```

Two of its outputs matter here (with `--out_prefix mygenome`):

- **`mygenome_r1_ltr.tsv`** — the LTR-RT table. **Column 1 is the `chrom:start-end` locus and
  column 7 is the LTR–LTR p-distance** (the age proxy). Pass this to `bedtools.py -pass_scn` for
  [benchmarking](#benchmark-your-te-annotations-first) and to the grid scorer (as `--ref-tsv` and
  the per-simulation `--exp-tsv-name`).
- **`mygenome_depth0_ltr.fa`** — a FASTA of the (unnested) LTR-RT sequences, with
  RepeatMasker-style headers `>chrom:start-end#LTR/superfamily/...`.

(synLTR needs `mafft` and `trimal` on PATH and clones Kmer2LTR on first run, so that run needs
network access.)

**Build a species-specific library** from that FASTA. `Kmer2LTR --make-perfect-ltr-rt 5p`
rebuilds each element with two *identical* LTRs (as at the moment of insertion) and appends the
`~LTRlen:<N>` tag PrinTE needs. `5p` (use the 5′ LTR) is the natural choice:

```bash
python Kmer2LTR/Kmer2LTR.py -i mygenome_depth0_ltr.fa --make-perfect-ltr-rt 5p
# -> mygenome_depth0_ltr.LTRs.alns.perfect_5p.fa
```

Feed that file to PrinTE as a ready-made library: its headers already use the
`name#class/superfamily~LTRlen:<N>` convention, so PrinTE uses them as-is (skips library processing):

```bash
bash PrinTE/PrinTE.sh ... --clean_lib mygenome_depth0_ltr.LTRs.alns.perfect_5p.fa
```

**Why this matters biologically:** the simulation now inserts the *actual* LTR-RTs of your species
 (their real sequences, length distribution, and superfamily mix) so the simulated TE landscape,
and the genome-size trajectory it drives, reflect your genome's own elements rather than a generic
library. (`--make-perfect-ltr-rt` **requires** a mode; a bare `--make-perfect-ltr-rt` errors.
Kmer2LTR only appends `~LTRlen:N`; the `#class/superfamily` comes from synLTR's headers, so the
chain works because synLTR emits RepeatMasker-style names.)

---

## Grid search: finding parameters without bias

Both tracks need the same thing: the TE insertion/deletion rates (and solo-LTR and length-bias
settings) that make a simulation match a target. Rather than guess, PrinTE's `grid/` tools
**search the parameter space reproducibly** and score every simulation against a reference.

The search space is four dimensions: **insertion rate × deletion rate × solo ratio (`sr`) ×
length bias (`k`)**. Two front-ends are provided:

- **`grid/gridsearch.py`** — random/grid sampling of the space (simple, exhaustive-ish).
- **`grid/guided_search.py`** — *active learning*: round 0 explores with Latin-hypercube
  sampling; later rounds train a random-forest surrogate on the results so far and propose the
  most informative next batch. This converges on the target with far fewer simulations.

**What the grid writes to disk:** every sampled parameter combination becomes its own directory
named `insertion_rates_<x>_deletion_rates_<y>_solo_ratio_<s>_length_bias_<k>/`, each holding that
run's `gen<N>_final.fasta/.bed/.lib`. The re-annotation and scoring steps below iterate over
these directories automatically.

**You will need a cleaned TE library** (`--clean-lib`). PrinTE writes `lib_clean.fa` into the
working directory whenever it processes a `--TE_lib`; reuse that file. To build one directly:

```bash
# Clean a TE library to pure ACGT and pre-process it once; reuse as --clean-lib / -cl.
zcat PrinTE/data/maize_rice_arab_curated_TE.lib.gz \
  | python PrinTE/data/fix_non_ATGC.py - > lib_clean.fa
```

### Step 1 — seed and launch the search

```bash
python PrinTE/grid/guided_search.py init \
    --samples 100 --seed 7 \
    --ins-start 1e-9 --ins-end 1e-17 --del-start 1e-9 --del-end 1e-17 \
    --sr-start 5 --sr-end 95 --sr-step 5 --k-start 0 --k-end 10 --k-step 1 \
    --printe-script PrinTE/PrinTE.sh \
    --ge 5400000 --st 54000 --mut 2e-9 --tstv 2.0 \
    --mxgs 211730673 --mngs 211530673 \
    --bed burnin.bed --fasta burnin.fasta \
    --clean-lib lib_clean.fa --ratios PrinTE/ratios_ltr_only.tsv \
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
  nohup python PrinTE/grid/run_loop.py --target-hits 200 --max-rounds 500 \
        --samples 100 --explore-frac 0.2 > orchestrator.log 2>&1 &
  ```

  `run_loop.py` harvests results, proposes the next batch, and re-launches each round itself —
  run it from the directory containing `search_state.json` (written by `init`). It only drives
  the **local** launcher.

- **SLURM.** Add `--scheduler slurm --slurm-partition <p> --slurm-cpus 48 --slurm-mem 120G
  --slurm-time 24:00:00` to the `init` command; it writes `slurm_grid/submit_array.sbatch`, which
  you submit with `sbatch slurm_grid/submit_array.sbatch`. On SLURM, drive the rounds manually
  (`guided_search.py harvest` → `next` → `sbatch`) rather than with `run_loop.py`.

### Step 2 — re-annotate the simulated genomes

The grid runs skip post-processing for speed, so the simulated genomes are **not** dated
automatically. Annotate each one's LTR-RTs with the **same tool you used to benchmark**
([synLTR](#annotating-ltr-rts-and-building-a-species-specific-library)), producing in every
simulation directory a `<genome>_r1_ltr.tsv` table (**column 1 = `chrom:start-end`, column 7 =
p-distance** — exactly what the scorer reads). `grid/discover_ltr.sh` is a template that loops
over every simulation directory and runs the annotation; edit its target FASTA name and paths to
match your run. Build the **reference** table the same way on your real genome.

### Step 3 — score and visualize

```bash
python PrinTE/grid/build_composite_matrix.py \
    --compare-script PrinTE/grid/compare_genomes.py \
    --ref-tsv real_genome_LTR.tsv --ref-fasta real_genome.fa \
    --gen-prefix gen5400000_final.fasta \
    --exp-tsv-name gen5400000_final_r1_ltr.tsv \
    --dist-metric rms --alphas 0 0 10 1 \
    --output composite_matrix.tsv --threads 64

# Best-fitting parameters = the row with the lowest Composite (column 20).
(head -1 composite_matrix.tsv; tail -n +2 composite_matrix.tsv | sort -t$'\t' -k20,20n) | head

# Visualize how each parameter affects the fit.
python PrinTE/grid/contour_plot.py --input composite_matrix.tsv
```

- **`real_genome.fa` / `real_genome_LTR.tsv`** are your real target assembly and *its* LTR-RT
  table, produced by the **same** re-annotation command as the simulations (so the columns match).
- **`--gen-prefix`** must equal the FASTA filename inside each simulation directory;
  **`--exp-tsv-name`** is the basename of the per-simulation LTR TSV inside each directory. If your
  re-annotation wrote a different name (e.g. `<prefix>_kmer2ltr_dedup`), either rename it in each
  directory or pass that name to `--exp-tsv-name`.
- Pass `--compare-script PrinTE/grid/compare_genomes.py` explicitly — the script's built-in
  default points at a differently named file.
- **Reading the winner:** the first four columns of the top row are your best-fit
  `insertion_rate`, `deletion_rate`, `solo_ratio`, `length_bias`. `--alphas` weights the four
  distance terms (count, length, genome size, age distribution); `0 0 10 1` weights genome size
  most heavily and age distribution next.

`contour_plot.py` writes pairwise contour maps of the Composite score and prints a sensitivity
report (parameter ranges, correlations, and confidence intervals for the optimum).

---

## Track 1 — Past → present

**Goal:** reconstruct the genome and TE state of a clade's common ancestor, build it as a
burn-in, and evolve it forward until the simulation resembles the **real** present-day genome.
Because the real genome is ground truth, you can score the fit and ask whether TE activity
**increased or decreased** along the branch.

You should already have a **time-scaled Newick tree** (`tree.nwk`) and **each tip's genome
FASTA**. Run [synLTR](#annotating-ltr-rts-and-building-a-species-specific-library) on each tip's
genome with `--out_prefix <tip>` to produce, per tip, an LTR-RT table `<tip>_r1_ltr.tsv` (column
7 = p-distance/age) and an LTR-RT FASTA `<tip>_depth0_ltr.fa`. So each tip has three files (
genome, table, LTR-RT FASTA) which drive the three reconstructions below. (Distance trees can be
time-scaled with a tool such as [PATHd8](https://www2.math.su.se/PATHd8/).)

### Step 1 — reconstruct the ancestor

Two R scripts run `phytools::fastAnc` over the tree. **File-naming/placement matters:**

- `ancestral_reconstruction_gs.R` looks for each tip's file as `<tip_label><suffix>` **in the
  same directory as the Newick** (`<tip>` is the exact tip string in the tree).
- `ancestral_reconstruction_ltr_age.R` looks for `<tip_label><suffix>` **in the current working
  directory** — so run it from where those files live.

```bash
# (a) Ancestral GENOME size  -> needs one genome FASTA per tip, named <tip>.fa, next to tree.nwk.
Rscript PrinTE/grid/ancestral_reconstruction_gs.R \
    --newick tree.nwk --suffix .fa --out ancestral_genome_size --label_nodes

# (b) Ancestral LTR-RT size  -> needs each tip's LTR-RT FASTA (synLTR's <tip>_depth0_ltr.fa)
#     next to tree.nwk; the script sums their bp.
Rscript PrinTE/grid/ancestral_reconstruction_gs.R \
    --newick tree.nwk --suffix _depth0_ltr.fa --out ancestral_ltrrt_size --label_nodes

# (c) Ancestral LTR-RT AGE distribution  -> run from the dir holding each tip's LTR table
#     (synLTR's <tip>_r1_ltr.tsv; TAB-separated, col 7 = p-distance/age in 0-1).
Rscript PrinTE/grid/ancestral_reconstruction_ltr_age.R \
    --newick tree.nwk --suffix _r1_ltr.tsv
```

(`gs.R` needs `samtools` on PATH and builds `<file>.fai` automatically. `ltr_age.R` calls
`grid/summary_stats.py` via `python` — keep it beside the R script, or pass `--summary_path` —
and bins ages with defaults `--bins 50 --bin_max 0.15`.)

**Find your node.** `gs.R` writes `ancestral_genome_size.contMap.pdf` (node numbers show only
with `--label_nodes`) and `ancestral_genome_size.ancestral_genome_sizes.tsv`. `ltr_age.R` writes
`tree_with_node_ids.pdf` and `node_map.tsv` (each internal node → its full set of descendant
tips). Open a PDF (or `grep` your tip names in `node_map.tsv`), find the internal node whose
descendants are exactly your clade, and **use that node number consistently below** (here, `20`).

### Step 2 — turn the reconstruction into PrinTE inputs

```bash
# Ancestral genome size and intact-LTR-RT size are column 6 (bp_est) of each gs.R table.
gs=$(awk -F'\t' '$1==20{print $6}' ancestral_genome_size.ancestral_genome_sizes.tsv)
ltr=$(awk -F'\t' '$1==20{print $6}' ancestral_ltrrt_size.ancestral_genome_sizes.tsv)

# -sz must carry an Mb suffix; -itp is the intact-LTR-RT percent of the genome.
sz=$(awk -v g="$gs" 'BEGIN{printf "%.2fMb", g/1e6}')          # e.g. 180.27Mb
itp=$(awk -v g="$gs" -v l="$ltr" 'BEGIN{printf "%.3f", 100*l/g}')  # e.g. 16.629
echo "-sz $sz   -itp $itp"

# Build the age-distribution file for -mb from ltr_age.R's ancestral_summary_bins.tsv.
# Columns there are: stat_type, node, clade_size, clade_signature, descendant_tips,
#                    bin_index, bin_start, bin_end, est, CI_lower, CI_upper.
# Pull your node's bins (start,end,count = cols 7,8,9) and round the counts:
awk -F'\t' '$2==20' ancestral_summary_bins.tsv | cut -f7-9 \
  | awk -F'\t' 'BEGIN{OFS="\t"}{ $3=sprintf("%.0f",$3); print }' > ancestral.LTR.bins

# Normalize counts to frequencies -> 3-column "bin_start bin_end frequency" file that -mb expects:
awk -v OFS='\t' '{a[NR]=$0; v[NR]=$3; s+=$3} END{for(i=1;i<=NR;i++) print a[i], v[i]/s}' \
  ancestral.LTR.bins | cut -f1,2,4 > ancestral.LTR.bins.freq
```

### Step 3 — build the ancestor as a burn-in

```bash
bash PrinTE/PrinTE.sh --burnin_only \
    -sz "$sz" -itp "$itp" -ftp 50 \
    -mb ancestral.LTR.bins.freq \
    -m 2e-9 -TsTv 2.0 --ex_LTR \
    -cl lib_clean.fa -r PrinTE/ratios_ltr_only.tsv -t 200
```

`-sz`, `-itp`, and `-mb` come from Step 2. `-ftp` is the **fragmented**-TE percent — it is *not*
reconstructed above (the gs.R/ltr_age.R steps cover intact LTR-RTs only), so estimate it from
your clade's present-day fragment + solo load, or use it as filler/deletion substrate. `-cl
lib_clean.fa` is the cleaned TE library; for a more realistic ancestor, build a
**species-specific** LTR-RT library from a close present-day relative
(see [Annotating LTR-RTs](#annotating-ltr-rts-and-building-a-species-specific-library)) and pass it here.

### Step 4 — benchmark, evolve, and fit

Optionally [benchmark](#benchmark-your-te-annotations-first) the burn-in (run your LTR pipeline on
`burnin.fasta`; a high recall on a genome with a known TE landscape means later mismatches reflect
**biology**, not annotation error). Then use the [grid search](#grid-search-finding-parameters-without-bias)
to find the insertion/deletion rates (and `sr`, `k`) whose terminal genome best matches the real
one.

**Interpret.** The best-fit rates are your estimate of TE turnover along the branch. Compare them
to the ancestral and present-day TE loads to decide whether activity rose or fell — and use the
ancestral reconstruction to motivate which hypothesis (increasing vs decreasing rate) to test.

---

## Track 2 — Present → future

**Goal:** take a real present-day genome and project it forward. There is no future genome to
compare against, so instead of fitting, you **scan for stability**: under which PrinTE parameters
does a faithful replica of the genome stay roughly unchanged (in size, TE content, and age
distribution) over a short horizon (~1 My)?

**1. Build a digital replica** of your genome by running
[TEgenomeSimulator](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator) yourself. It
emits a genome FASTA (`TEgenomeSimulator.fa`) and a TE GFF3 (`TEgenomeSimulator.gff`). Unlike
Track 1, TEgenomeSimulator, not PrinTE, creates the burn-in here.

**2. Convert it to PrinTE format** and (optionally) confirm it loads:

```bash
python PrinTE/util/gff_to_bed.py TEgenomeSimulator.gff TEgenomeSimulator.bed
# Sanity check: one short evolve step should run cleanly.
bash PrinTE/PrinTE.sh -f TEgenomeSimulator.fa -b TEgenomeSimulator.bed -i my_TE.lib -ge 100000 -st 100000 -t 20
```

**3. Scan for the stability regime.** Use the same [grid search](#grid-search-finding-parameters-without-bias),
but with the **replica genome as the reference**, a **tight size band** around the current size,
and a **short horizon**. For example, with a current genome of ~180,000,000 bp:

```bash
python PrinTE/grid/guided_search.py init --samples 100 --seed 7 \
    --ins-start 1e-10 --ins-end 1e-13 --del-start 1e-10 --del-end 1e-13 \
    --sr-start 5 --sr-end 95 --sr-step 10 --k-start 0 --k-end 10 --k-step 2 \
    --printe-script PrinTE/PrinTE.sh --ge 1000000 --st 100000 --mut 1.3e-8 --tstv 2.0 \
    --mxgs 181800000 --mngs 178200000 \
    --bed TEgenomeSimulator.bed --fasta TEgenomeSimulator.fa \
    --clean-lib lib_clean.fa --ratios PrinTE/ratios.tsv --threads 64
```

A run **"held steady"** if its genome size and TE content stayed inside your tolerance band over
the horizon. The lightest-weight read-out is each run's `genome_size_plot.pdf` and `percent_TE.pdf`
(flat = stable); for a single number, score each run against the replica with
`build_composite_matrix.py` (low Composite = stayed close to the start). Because there is no ground
truth, treat Track 2 as **hypothesis generation** about genome stability — report *which rate
regimes maintain the genome* and which let it expand or erode, not a single historical rate.

---

## Preparing a custom TE library

PrinTE ships with a curated maize/rice/Arabidopsis library
(`data/maize_rice_arab_curated_TE.lib.gz`) and an LTR-RT exemplar database (`data/ltr-db.fa.gz`).
For LTR-RT studies, the best library is usually a **species-specific** one built from your own
genome — see [Annotating LTR-RTs](#annotating-ltr-rts-and-building-a-species-specific-library).
To assemble a library by hand instead, a few helpers are provided:

- **`util/TE_lib_stitcher.py`** — reassemble split `LTR` + `Internal` library entries into
  intact LTR-RTs.
- **`data/fasta_to_RepeatMasker.py`** — normalize heterogeneous headers into
  `>name#class/superfamily` form.
- **`data/fix_non_ATGC.py`** — replace non-ACGT characters with `T` so downstream k-mer/alignment
  steps don't choke (streams large files; writes to stdout):

  ```bash
  zcat custom.lib.gz | python PrinTE/data/fix_non_ATGC.py - > my_TE.lib
  ```

Provide the result via `--TE_lib my_TE.lib` (PrinTE processes it and writes a cleaned
`lib_clean.fa` you can reuse) or, if already cleaned/pre-processed, via `--clean_lib lib_clean.fa`
to skip library processing.

---

## Requirements & supported systems

Supported on **Linux** and **macOS**. Tested on macOS 14.6.1 (Apple M2 Pro, 16 GB RAM) and
Ubuntu 22.04.5 LTS (AMD EPYC 7763). Requires ~2 GB RAM minimum; 16 GB+, 4+ cores, and a network
connection (to clone Kmer2LTR) are recommended. The bundled `ltr_mutator` binary is rebuilt
automatically from source (`g++ -std=c++17 -O3 -fopenmp`) if the precompiled one won't run.
