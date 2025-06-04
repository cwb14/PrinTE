# TESS: TE Evolution Simulation Suite

**Version**: 1.00  

---

## Components

- **TESS.sh**  
  - Entry-point script that dispatches subcommands (`RandSeqInsert`, `PrinTE`, `TEgenomeSimulator`).  
  - Manages `--verbose` logging and `--force-update` flags.  
  - Clones or updates the external TEgenomeSimulator repository at runtime.  
  - Decompresses input FASTA files automatically (if needed).  

- **RandSeqInsert/**  
  - `RandSeqInsert.py`: Python script that simulates TIR insertions.  
  - Detailed usage: [RandSeqInsert README](https://github.com/cwb14/TESS/blob/main/RandSeqInsert/README.md)

- **TEgenomeSimulator/**  
  - Populated automatically by cloning [`https://github.com/Plant-Food-Research-Open/TEgenomeSimulator`](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator) when `TESS.sh TEgenomeSimulator` is run for the first time.  
  - Contains scripts, test data, and dependencies needed to run `TEgenomeSimulator.py`.  
  - Detailed usage: [TEgenomeSimulator GitHub](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator/tree/main)

- **PrinTE/**  
  - prinTE.sh: Bash wrapper that performs forward TE simulations.  
  - Usage instructions are embedded at the top of [`prinTE/prinTE.sh`](https://github.com/cwb14/TESS/blob/main/prinTE/prinTE.sh).  

---

## Prerequisites

Before installing or running **TESS**, ensure you have:

1. **mamba** (or **conda**) for environment management.  
2. **Git** (for cloning TEgenomeSimulator).  
3. **C++ compiler** supporting C++17 (e.g., `g++`).  
4. **Python 3.9** (to match the environment instructions below).  
5. Standard Unix utilities (`bash`, `gunzip`, etc.).  

---

## Installation

### 1. Create Conda Environment

Use **mamba** (recommended) or **conda** to create and activate a new environment named `TESS`:

```bash
mamba create -n TESS -c conda-forge biopython matplotlib seaborn numpy pandas python=3.9 pyyaml scipy setuptools rmblast graphviz r-viridis

# Activate the environment:
conda activate TESS
```

### 2. Clone TESS Repository

```bash
git clone https://github.com/cwb14/TESS.git
chmod +x ./TESS/TESS.sh
```

---

## Usage

```bash
./TESS/TESS.sh <command> [options]
```

### Global Flags

- `--verbose`  
  Enable detailed logging (prints INFO messages to stderr).  
- `--force-update`  
  Force removal and re-clone of the `TEgenomeSimulator` repository (if it already exists).  

### Available Commands

| Command             | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `RandSeqInsert`     | Run `RandSeqInsert.py` to simulate TIR insertions.                          |
| `TEgenomeSimulator` | Clone or update TEgenomeSimulator, then run `TEgenomeSimulator.py`.          |
| `PrinTE`            | Run `prinTE/prinTE.sh` to simulate LTR mutation/insertion (forward evolution).  |

---

### 1. RandSeqInsert

Invokes the Python script under `RandSeqInsert/`. Detailed usage can be found at:

[RandSeqInsert README](https://github.com/cwb14/TESS/blob/main/RandSeqInsert/README.md)

Basic examples:

```bash
# Help menu.
./TESS/TESS.sh RandSeqInsert -h

# Create a blank-slate fake genome (burnin.fa) for use with RandSeqInsert. 
./TESS/TESS.sh PrinTE --burnin_only --cds_percent 0 --TE_percent 0 --chr_number 1 --size 100Mb
rm pipeline.* TAIR10.cds.fa combined_curated_TE_lib_ATOSZM_selected.fasta lib* backbone.* burnin_mut_dist.pdf burnin.bed

# Run RandSeqInsert.
./TESS/TESS.sh RandSeqInsert -i burnin.fa -is 50 -it 3 -d ./TESS/RandSeqInsert/lib/TIR/maize/ -w 0.7 --tsd 9 --track --visual
```

---

### 2. TEgenomeSimulator

Clones (if necessary) the external TEgenomeSimulator repository and executes `TEgenomeSimulator.py`. Detailed usage can be found at:

[TEgenomeSimulator GitHub](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator/tree/main)

Basic examples:

```bash
# Help menu.
./TESS/TESS.sh TEgenomeSimulator -h

# Random Synthesized Genome mode.
./TESS/TESS.sh TEgenomeSimulator -M 0 -p random_mode -c TESS/TEgenomeSimulator/test/input/random_genome_chr_index.csv -r TESS/TEgenomeSimulator/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta -m 5 -n 1 -o random_mode_test
```

---
### 3. PrinTE

The `PrinTE` command launches `prinTE/prinTE.sh`. Please refer to the top of the script for all supported options:

[prinTE/prinTE.sh](https://github.com/cwb14/TESS/blob/main/prinTE/prinTE.sh)

Basic examples:

```bash
# Help menu.
./TESS/TESS.sh PrinTE -h

# Variable-rate method.
./TESS/TESS.sh PrinTE -cn 5 -sz 135Mb -tmx 5 -m 7e-9 -P 25 -p 21 -br 1e-7 -ir 1.1e-6 -dr 4e-6 -cbi 1.1 -cbd 1.0 -cb 500 -t 20 -k 2 -ge 40000 -st 10000

# What was the gene/TE landscape of the starting genome?
cat burnin.stat 

# How many TEs were inserted and deleted?
cat pipeline.report

# How many TEs were inserted due to vertical acquisition?
cat pipeline.log | grep 'Number of born TEs to insert'

# Fixed-rate method.
./TESS/TESS.sh PrinTE -mgs 1500M -P 20 -n 6000 -cn 20 -sz 113Mb -ge 300000 -st 100000 -t 10 -k 0 -kt -F 3.0e-11,5e-11 -m 1.3e-8 -sr 95 

# How'd it do?
cat burnin.stat 
cat pipeline.report

# We can increase '-ge' and add '--continue' to pickup where we resume with more generations. 
./TESS/TESS.sh PrinTE -mgs 1500M -P 20 -n 6000 -cn 20 -sz 113Mb -ge 400000 -st 100000 -t 10 -k 0 -kt -F 7.0e-11,1e-11 -m 1.3e-8 -sr 95 --continue
```

---

## Troubleshooting

### Permission Denied

- Confirm `TESS.sh` is executable:
  ```bash
  chmod +x ./TESS/TESS.sh
  ```
- Ensure you have write permissions in this directory if re-cloning `TEgenomeSimulator`.

### C++ Compiler Errors (`ltr_mutator`)

- Ensure you have a modern C++ compiler supporting C++17 (e.g., `g++` ≥ 7.0).  
- The compile command is:
  ```bash
  g++ -std=c++17 -fopenmp -O3 \
    -o ./TESS/prinTE/bin/ltr_mutator \
    ./TESS/prinTE/bin/ltr_mutator.cpp
  ```
- If errors persist, inspect the source code and flags, or install a newer `g++`.

### Python Not Found

- TESS expects either `python3` or `python` on your `$PATH`.  
- Install Python 3.9 if missing, or create a symlink:
  ```bash
  ln -s /usr/bin/python3 /usr/local/bin/python
  ```

---

## Configuration & Customization

### Force Update

To force-reclone TEgenomeSimulator (e.g., after pulling new changes upstream), add `--force-update`:

```bash
./TESS/TESS.sh TEgenomeSimulator --force-update
```

---

## Citation

If you use **TESS** in your research, please cite:

> Benson CW, Chen T-H, Lu T, Angelin-Bonnet O, Thomson SJ, Deng CH, Ou S. (2025).  
> **TESS: A Forward Simulation Framework for Studying the Role of Transposable Elements in Genome Expansion and Contraction.**

---

## License

This project is released under the **GNU GPL v3** license. See the full text in `LICENSE.md`.
