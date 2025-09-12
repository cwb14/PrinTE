# PrinTE: TE Evolution Simulation Pipeline

**PrinTE** is a comprehensive pipeline to orchestrates a two‐phase transposable element (TE) evolution simulation pipeline, followed by supplemental post‐processing and reporting. It automates:

1. **Phase 1 (Burn-in)**  
   Builds a synthetic genome, inserts an initial set of TEs, and generates “burn-in” outputs (e.g., `burnin.fa`, `burnin.bed`), unless the user provides their own starting FASTA/BED files or restarts from a previous run.

2. **Phase 2 (Looping Generations)**  
   For each generation (in user-specified steps), the pipeline:
   - Mutates the genome.  
   - Inserts new TEs.  
   - Excises TEs.  
   - Extracts intact TEs to form the next-generation TE library.  
   - Logs progress and checks for maximum genome size.

3. **Supplemental Post-Processing**  
   After all generations complete, the pipeline runs a suite of reporting and plotting utilities:
   - TE fraction over time  
   - Fragmented vs. intact TE proportions  
   - Superfamily counts  
   - Category bar plots  
   - Genome size trajectory  
   - Per-generation LTR extraction and divergence analyses  
   - Overall LTR density plot  

PrinTE treats all TEs as Type 1 (copy-and-paste). It supports simulating Type 2 (cut-and-paste) TEs, but does not currently model the cut-and-paste aspect of their biology. 

---

## Installation

PrinTE relies on **conda** and **mamba**. We recommend installing them via [Miniforge](https://github.com/conda-forge/miniforge/releases) (choose the installer for your OS and CPU architecture).

Then run:

```bash
git clone https://github.com/cwb14/PrinTE.git
mamba env create -f PrinTE/env.yml
conda activate PrinTE
```

PrinTE installs and uses [Kmer2LTR](https://github.com/cwb14/Kmer2LTR.git) for LTR-RT dating. 

---

## Phase 1 (Burn-in) - Inputs and Parameters  
   - PrinTE's initial (**burn-in**) utility is quite flexible, allowing simulations that mirror the composition of a wide range any real genome.  
   - A caveat is that PrinTE does not simulate solo-LTRs in burn-in genomes. For solo-LTRs, be sure to proceed to Phase 2 (looping generattions).

### Inputs

1. **CDS FASTA (`--cds`)**  
   Specifies the **genes** to insert into the synthetic genome.  
   - By default, PrinTE uses *Arabidopsis thaliana* TAIR10 CDS, which contains 19,621 sequences.  
   - If you request more CDS sequences than this (via `--cds_percent` or `--cds_num`), you must provide an alternative CDS FASTA file.

2. **TE Library FASTA (`--TE_lib`)**  
   Specifies the **transposable elements** to insert into the synthetic genome.  
   - Default: `maize_rice_arab_curated_TE.lib.gz`  
   - FASTA headers must follow **RepeatMasker format**, e.g.:  
     ```
     >[name]#[class]/[superfamily]
     >Os2670#MITE/Tourist
     ```  
   - Supported `#[class]/[superfamily]` suffixes are listed in `ratios.tsv`.
   - Users can control the abundance of **intact TEs** in the burn-in genome using `--intact_TE_percent` or `--intact_TE_num` and **fragmented TEs** using `--frag_TE_num` or `--frag_TE_percent`.
   - Unlike genes, where each CDS is inserted a maximum of 1 time, TEs are selected randomly from the library and may be inserted any number of times in ratios matching `ratios.tsv`. So, larger TE libraries yield more diverse TE landscapes.
   - Note that in the burn-in, there are no nested insertions. All genes and TEs are laid out in tandem at random distances appart (minimum 20bp).

3. **TE Ratios File (`--TE_ratio ratios.tsv`)**  
   Defines the relative frequency of each TE superfamily in the genome.  
   - Columns:  
     1. `class`  
     2. `superfamily`  
     3. `weight` of **intact TE** (probability of inserting intact TE)
     4. `weight` of **fragmented TE** (probability of inserting fragmented TE)
   - PrinTE ships with `ratios.tsv` for this purpose.  
   - Users can adjust the `weight` values (column 3 & 4) with to tune the TE landscape to their need.

    ```bash
    # TE_class   TE_superfamily   intact_frequency	fragmented_frequency
    DNA         Helitron          0.03	0.07
    SINE        unknown           0.05	0.06
    LTR         Gypsy             0.19	0.21
    ```
    For example, with `ratios.tsv` above, and the options  `--intact_TE_percent 20 --frag_TE_percent 10`:  
        - **Intact `DNA/Helitron`** -> `20 * 0.03 = 0.6%` of the genome.  
        - **Fragmented `DNA/Helitron`** -> `10 * 0.07 = 0.7%` of the genome.  

### Parameters

1. **`--chr_number`**  
   Number of chromosomes in the burn-in genome.  
   - Parallelization occurs on a per-chromosome basis, so more chromosomes generally improve runtimes.  
   - Avoid oversplitting the genome into too many small chromosomes.  

2. **`--size`**  
   Total size (bp) of the burn-in genome.  

3. **`--TE_mut_k`** and **`--TE_mut_Mmax`**  
   Control the distribution of TE substitution mutations in the burn-in genome.  
   - TE Mutations are clock-like, but older TEs are more likely to have been deleted or fragmented than younger. This means that the observed mutation landscape appears asymetric, with most intact TEs appear young, while older TEs are often lost.  
     - **`--TE_mut_k`** sets the *decay rate*:  
       - Large values -> steep decay (most TEs with very low mutation rates).  
       - Small values -> flatter distribution (higher mutation rates more likely).  
     - **`--TE_mut_Mmax`** sets the *ceiling* for mutation percentage.  
   - **Use `--TE_mut_Mmax 0` to disables TE mutations in the burn-in.**  
   - PrinTE generates a PDF of the mutation decay function: `burnin_mut_dist.pdf`.  

4. **`--burnin_only`**  
   Runs only Phase 1 (burn-in).  
   - Generates the synthetic genome and stops before evolution steps (Phase 2).


  ## Phase 2 (Looping Generations) - Parameters

With the burn-in genome created, Phase 2 specifies how to evolve it across generations.  
PrinTE implements two approaches for TE evolution:

1. **Fixed rates** – constant TE insertion and deletion rates **per base of the genome per generation**.  
2. **Variable rates** – The **per base per generation** insertion and deletion rates depend on the TE landscape of the previous generation.

In both approaches, PrinTE dynamically updates the TE library:  
- Only **intact TEs** are eligible for insertion into the next generation. Accumulated mutations are inherited to the transposed copy.   
- The exception is when `--birth_rate` is specified (Variable approach), which allows reintroduction of TEs from the original library.

---

### Fixed Approach

The **Fixed** method is simpler: you specify constant TE insertion and deletion rates.  

- **`--fix`**  
  TE insertion and deletion rates **per base per generation**.  
  - Format: `--fix <insertion_rate>,<deletion_rate>`  
  - Example:  
    ```bash
    --fix 5e-9,1e-8
    ```
    -> insertion rate = `5e-9`, deletion rate = `1e-8`.  

- **`--disable_genes`**  
  Prevents TE insertions into genes.

---

### Variable Approach

The **Variable** method modulates TE dynamics based on the composition of the previous generation. TE mutations accrue.
Insertion and deletion rates are **per intact TE per generation**, and PrinTE internally converts them into **per-base** rates.

#### Core Parameters
- **`--insert_rate`**  
  TE insertion rate **per intact TE per generation**.  

- **`--delete_rate`**  
  TE deletion rate **per intact TE per generation**.  

- **`--birth_rate`**  
  Rate of horizontally acquired TEs.  
  - Normally, only intact TEs are propagated.  
  - With `--birth_rate`, PrinTE occasionally samples from the **original TE library** (`--TE_lib`), simulating horizontal transfer or reintroduction of extinct TE lineages.  

#### Selective Constraints
The genome is modeled as:  
1. **Selectively constrained sequence** (genes).  
2. **Neutral sequence** (TEs + intergenic space).  

- **`--sigma`**  
  Controls how unevenly selective constraints are distributed across genes (via a log-normal distribution).  
  - Low values -> constraints clustered (most genes similarly constrained).  
  - High values -> uneven constraints (a few highly constrained, most moderately constrained).  
  - PrinTE outputs `lognormal_distribution.pdf` to visualize the distribution.  

- **`--sel_coeff`**  
  By default, PrinTE does not purge TEs based on fitness.  
  - With this parameter, insertions with higher deleterious effects are **more likely to be deleted**.  

#### Chromatin Effects
TE insertions may occur more readily in euchromatin than heterochromatin.  
PrinTE approximates this by treating **genes as euchromatin** and **non-gene regions as heterochromatin**.

- **`--chromatin_buffer`**  
  Expands the euchromatin region around genes.  

- **`--chromatin_bias_insert`**  
  Increases the probability of TE **insertion** into euchromatin.  

- **`--chromatin_bias_delete`**  
  Increases the probability of TE **deletion** in euchromatin.  

---

In summary:  
- Use **Fixed** for simple, constant-rate simulations.  
- Use **Variable** for more realistic dynamics shaped by genome composition, selective constraints, and chromatin context.

   
## General Parameters

These options apply across both Phase 1 and Phase 2 of the pipeline.

- **`--generation_end`**
  The total number of generations to simulate.

- **`--step`**
  The number of generations to simulate in a single step.
  - E.g., `--generation_end 3000 --step 1000` simulates generation 1000, 2000, and 3000.

- **`--mutation_rate`**  
  Rate of DNA substitutions applied to the genome.  

- **`--max_size`**  
  Terminates Phase 2 early if the genome size exceeds this threshold, even if `--generation_end` has not been reached.  

- **`--seed`**  
  Random seed for reproducibility.  

- **`--threads`**  
  Number of threads to use for parallel execution.  

- **`--solo_rate`**  
  When an LTR-RT is selected for deletion, it may be removed entirely or partially, leaving behind a solo-LTR.  
  - This parameter controls the frequency of solo-LTR formation.  

- **`--k`**  
  Longer TEs are more likely to undergo unequal or illegitimate recombination, leading to their deletion.  
  - **`--k`** controls the slope of this length–bias curve.  
  - Use `--k 0` to disable length-biased deletion.
  - See `weighted_candidate_selection.pdf` to visualize the deletion weight by TE lenght dustribution. 

- **`--fasta`, `--bed`**
  Provide genome fasta and PrinTE-formatted bed to skip the burn-in phase, using these files instead.

- **`--continue`**
  Tells PrinTE to look for existing outputs in the working directory and pickup where it left off. 

- **`--keep_temps`**  
  Retain intermediate temporary files instead of cleaning them up automatically.  

- **`--model`**  
  Substitution model of DNA evolution used for dating LTR-RTs.  

- **`--TsTv`**  
  Transition/transversion ratio applied during substitution modeling.  

- **`--ex_LTR`**  
  From the TE library, exclude LTR-RTs that lack detectable LTRs.


  ---

## Usage
See help menu.
```bash
bash PrinTE/PrinTE.sh
```

Create a blank-slate genome (burnin.fa). 
```bash
bash PrinTE/PrinTE.sh --burnin_only --cds_percent 0 --TE_percent 0 --chr_number 1 --size 100Mb
```


Variable-rate method.
```bash
bash PrinTE/PrinTE.sh -cn 5 -sz 135Mb -tmx 5 -m 7e-9 -P 25 -p 21 -br 1e-7 -ir 1.1e-6 -dr 4e-6 -cbi 1.1 -cbd 1.0 -cb 500 -t 20 -k 2 -ge 40000 -st 10000
```

Fixed-rate method.
```bash
bash PrinTE/PrinTE.sh -mgs 1500M -P 20 -n 6000 -cn 20 -sz 113Mb -ge 300000 -st 100000 -t 10 -k 0 -kt -F 3.0e-11,5e-11 -m 1.3e-8 -sr 95 

# We can increase '-ge' and add '--continue' to pickup where we left off, adding more generations with new TE insertion/deletion rates. 
bash PrinTE/PrinTE.sh -mgs 1500M -P 20 -n 6000 -cn 20 -sz 113Mb -ge 400000 -st 100000 -t 10 -k 0 -kt -F 7.0e-11,1e-11 -m 1.3e-8 -sr 95 --continue
```
---

## Outputs

The primary outputs are:  
(1) The **genome fasta** (`gen[generation_number]_final.fasta`).  
(2) The **cooresponding bed** (`gen[generation_number]_final.bed`) showing gene and TE coordinates in the genome.  
(3) The **evolved TE library** (`gen[generation_number]_final.lib`).  
```bash
ls gen40000_final.*
gen40000_final.fasta
gen40000_final.bed
gen40000_final.lib
```

The `gen40000_final.bed` file looks like this:
```bash
chr1    1852998 1854855 gene1316        NA      +
chr1    1857719 1859314 tuteh_AC183372_584#LTR/unknown~LTRlen:126;CUT_BY:Os2721#DNAnona/hAT     CATTC   +
chr1    1859314 1859734 Os2721#DNAnona/hAT;NESTED_IN:tuteh_AC183372_584#LTR/unknown~LTRlen:126  CTTCCG  +
chr1    1859740 1860095 tuteh_AC183372_584#LTR/unknown~LTRlen:126;CUT_BY:Os2721#DNAnona/hAT     CATTC   +
chr1    1860701 1862054 gene1321        NA      +
chr1 	  1862094	1863435	TE_00016444_FRAG#MITE/DTA	NA	-
chr1	  1869245	1869502	anysaf_AC211487_11211#LTR/Ty3~LTRlen:257_SOLO	GCGCG	-
```
- Columns are `chromosome`, `start`, `end`, `feature_ID`, `target site duplication sequence (TSD)`, and `strand`.
- `Os2721#DNAnona/hAT` is nested inside `tuteh_AC183372_584#LTR/unknown`.
- `TE_00016444` has `_FRAG` in its feature_ID, indicating it was added into the burn-in using `--frag_TE_num` or `--frag_TE_percent`.
   - Neither `Os2721#DNAnona/hAT` or `TE_00016444_FRAG` are intact, so theyre not eligable for transpoition but avaible for deletion.
- `anysaf_AC211487_11211#LTR/Ty3` has the `_SOLO` tag on its `feature_ID`, so its a solo-LTR and also not viable for transposition.

---

`gen[generation_number]_mut.txt` is useful for looking at empirical mutation data:
```bash
cat gen40000_mut.txt                
Genome size:  136266136
Total mutations: 9487
Recurrent mutations: 1
Mutations / site:    6.962111261e-05
Mutations / site * 2:    0.0001392422252
Non-recurrent Mutations / site:    6.961377403e-05
Non-recurrent Mutations / site * 2:    0.0001392275481
Accumulated Mutations / site:    0.000265089172
```
- Lets assume I ran PrinTE with `-ge 400000 -st 100000`.  
- This means that, between `gen40000_mut.txt` and the preceeding simulation (`gen30000_mut.txt`), PrinTE introduced `9487` random mutations and only `1` was reccurent. 
- `Accumulated Mutations / site` shows us the frequency of mutations up to this point in the simulation (ie, distance from the original `burn-in` genome).

---

`burnin.stat` tells the composition of the original burnin genome:
```bash
The burn-in genome is 135000000bp in length with 19621 genes (21.21%) and 68300 TEs (51.16%).
  INTACT:  insertions=18300, TE bp inserted=27009913, (20.01% of genome)
  FRAG   : insertions=50000, TE bp inserted=42062617, (31.16% of genome)
```
---

`pipeline.report` tells us the number of insertions and deletions per generation simulated:
```bash
Generation	TE_inserts(nest/nonnest)	Calculated_TE_deletions	Actual_TE_deletions
10000	154(69/85)	2	2
20000	189(75/114)	13	13
30000	201(85/116)	23	23
40000	168(70/98)	16	16
```
---
If the variable-rate method was used with `--birth_rate`, I can check how many TEs were inserted horizontally:
```bash
cat pipeline.log | grep 'Number of born TEs to insert' 
Number of born TEs to insert (from birth_rate and birth_file): 14
Number of born TEs to insert (from birth_rate and birth_file): 14
Number of born TEs to insert (from birth_rate and birth_file): 14
Number of born TEs to insert (from birth_rate and birth_file): 14
```
---
- `all_LTR_density.pdf` shows a density distribution of LTR-RT ages through the simulation.  
- `genome_size_plot.pdf` plots the change in genome size.  
- `percent_TE.pdf` plots the change in TE abundance through generation simulated.   
- `stat_intact_plot.pdf` and `stat_frag_plot.pdf` show superfamily counts for each simulation.  
- `solo_intact.pdf` plots the changing solo-LTR ratio through generations.  
