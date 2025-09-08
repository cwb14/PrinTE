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
