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
   - Solo vs. intact TE proportions  
   - Superfamily counts  
   - Category bar plots  
   - Genome size trajectory  
   - Per-generation LTR extraction and divergence analyses  
   - Overall LTR density plot  

---

## Installation

PrinTE relies on **conda** and **mamba**. We recommend installing them via [Miniforge](https://github.com/conda-forge/miniforge/releases) (choose the installer for your OS and CPU architecture).

Then run:

```bash
git clone https://github.com/cwb14/PrinTE.git
mamba env create -f PrinTE/env.yml
conda activate PrinTE
```
---

## Phase 1 - Inputs and Parameters 
Allowing customization of the sequence composition for the initial (**burn-in**) genome.

### Inputs

1. **CDS FASTA (`--cds`)**  
   Specifies the **genes** to insert into the synthetic genome.  
   - By default, PrinTE uses *Arabidopsis thaliana* TAIR10 CDS, which contains 19,621 sequences.  
   - If you request more CDS sequences than this (via `--cds_percent` or `--cds_num`), you must provide an alternative CDS FASTA file.

2. **TE Library FASTA (`--TE_lib`)**  
   Specifies the **transposable elements** to insert into the synthetic genome.  
   - Default: `combined_curated_TE_lib_ATOSZM_selected.fasta`  
   - FASTA headers must follow **RepeatMasker format**, e.g.:  
     ```
     >[name]#[class]/[superfamily]
     >Os2670#MITE/Tourist
     ```  
   - Supported `#[class]/[superfamily]` suffixes are listed in `ratios.tsv`.
   - Users can control the abundance of TEs in the burn-in genome using `--TE_percent` or `--TE_num`.

3. **TE Ratios File (`--TE_ratio ratios.tsv`)**  
   Defines the relative frequency of each TE superfamily in the genome.  
   - Columns:  
     1. `class`  
     2. `superfamily`  
     3. `weight` (probability of insertion)
   - PrinTE ships with `ratios.tsv` for this purpose.  
   - Users can adjust the `weight` values (column 3) to tune the TE landscape to their needs.

### Parameters

1. **`--chr_number`**  
   Number of chromosomes in the burn-in genome.  
   - Parallelization occurs on a per-chromosome basis, so more chromosomes generally improve runtimes.  
   - Avoid oversplitting the genome into too many small chromosomes.  

2. **`--size`**  
   Total size (bp) of the burn-in genome.  

3. **`--TE_mut_k` and `--TE_mut_Mmax`**  
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
   
