# TESS: A Forward Simulation Framework for Studying the Role of Transposable Elements in Genome Expansion and Contraction

**Version**: 1.00  
**Repository**: [https://github.com/cwb14/TESS/](https://github.com/cwb14/TESS/)

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

- **PrinTE/**  
  - prinTE.sh: Bash wrapper that performs forward TE simulations.  
  - Usage instructions are embedded at the top of [`prinTE/prinTE.sh`](https://github.com/cwb14/TESS/blob/main/prinTE/prinTE.sh).  

- **TEgenomeSimulator/**  
  - Populated automatically by cloning [`https://github.com/Plant-Food-Research-Open/TEgenomeSimulator`](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator) when `TESS.sh TEgenomeSimulator` is run for the first time.  
  - Contains scripts, test data, and dependencies needed to run `TEgenomeSimulator.py`.  
  - Detailed usage: [TEgenomeSimulator GitHub](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator/tree/main)

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
cd TESS
chmod +x TESS.sh
```

### 3. Compile `ltr_mutator` (PrinTE)

Inside the `prinTE/bin/` directory, compile the C++ source:

```bash
cd prinTE/bin/
g++ -std=c++17 -fopenmp -O3 -o ltr_mutator ltr_mutator.cpp
```

Confirm that `ltr_mutator` exists and is executable:

```bash
ls -l ltr_mutator
# -rwxr-xr-x 1 user group <size> <date> ltr_mutator
```

---

## Usage

```bash
./TESS.sh <command> [options]
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
| `PrinTE`            | Run `prinTE/prinTE.sh` to simulate LTR mutation/insertion (forward evolution).  |
| `TEgenomeSimulator` | Clone or update TEgenomeSimulator, then run `TEgenomeSimulator.py`.          |

---

### 1. RandSeqInsert

Invokes the Python script under `RandSeqInsert/`. Refer to the full documentation here:

[RandSeqInsert README](https://github.com/cwb14/TESS/blob/main/RandSeqInsert/README.md)

Basic examples:

```bash
# Minimal usage
./TESS.sh RandSeqInsert

# Pass additional flags to RandSeqInsert.py, e.g., specifying length or output:
./TESS.sh RandSeqInsert --length 500 --output out.fasta
```

---

### 2. PrinTE

The `PrinTE` command launches `prinTE/prinTE.sh`, which runs the compiled `ltr_mutator` binary. Refer to the top of the script for all supported options:

[prinTE/prinTE.sh](https://github.com/cwb14/TESS/blob/main/prinTE/prinTE.sh)

Basic examples:

```bash
# Run forward TE simulation with input genome and output directory:
./TESS.sh PrinTE -i input_genome.fasta -o ltr_output/

# View full list of PrinTE options:
./TESS.sh PrinTE -h
```

---

### 3. TEgenomeSimulator

Clones (if necessary) the external TEgenomeSimulator repository and executes `TEgenomeSimulator.py`. Detailed usage can be found at:

[TEgenomeSimulator GitHub](https://github.com/Plant-Food-Research-Open/TEgenomeSimulator/tree/main)

Basic examples:

```bash
# First time: clone and run using default test inputs:
./TESS.sh TEgenomeSimulator

# Enable verbose logging:
./TESS.sh TEgenomeSimulator --verbose

# Force re-cloning of the TEgenomeSimulator repo:
./TESS.sh TEgenomeSimulator --force-update

# Pass additional TEgenomeSimulator.py arguments:
./TESS.sh TEgenomeSimulator \
  -M 0 -p test_random_mode \
  -c TEgenomeSimulator/test/input/random_genome_chr_index.csv \
  -r TEgenomeSimulator/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta \
  -m 5 -n 1 -o ./test/output
  --config TEgenomeSimulator/configs/tess_config.yaml \
  --genome TEgenomeSimulator/reference.fasta \
  --output sims/
```

Internally, TESS will:

1. Check for the `TEgenomeSimulator/` directory.  
2. If missing or `--force-update` is set:
   ```bash
   git clone --depth=1 https://github.com/Plant-Food-Research-Open/TEgenomeSimulator TEgenomeSimulator
   ```
3. Look for a compressed test FASTA at:
   ```bash
   TEgenomeSimulator/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz
   ```
   - If found and not decompressed, run `gunzip`.  
4. Resolve any relative paths in your arguments to absolute paths.  
5. Execute:
   ```bash
   python3 TEgenomeSimulator/TEgenomeSimulator.py <cleaned-args>
   ```

---

## Examples

```bash
# 1. Run RandSeqInsert with minimal parameters:
./TESS.sh RandSeqInsert

# 2. Run PrinTE, specifying an input and output directory:
./TESS.sh PrinTE -i my_genome.fasta -o ltr_results/

# 3. Run TEgenomeSimulator with verbose logging:
./TESS.sh TEgenomeSimulator --verbose --genome genomes/human_chr1.fasta --output sims/
```

---

## Troubleshooting

### Permission Denied

- Confirm `TESS.sh` is executable:
  ```bash
  chmod +x TESS.sh
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
./TESS.sh TEgenomeSimulator --force-update
```

---

## Citation

If you use **TESS** in your research, please cite:

> Benson CW, Chen T-H, Lu T, Angelin-Bonnet O, Thomson SJ, Deng CH, Ou S. (2025).  
> **TESS: A Forward Simulation Framework for Studying the Role of Transposable Elements in Genome Expansion and Contraction.**

---

## License

This project is released under the **GNU GPL v3** license. See the full text in `LICENSE.md`.
