# RandSeqGen

[![license](https://img.shields.io/github/license/lutianyu2001/RandSeqGen.svg)](https://github.com/lutianyu2001/RandSeqGen/blob/master/LICENSE)

RandSeqGen is a high-performance Python tool for generating random DNA sequences. It supports both purely random sequence generation and reference-based sequence generation with customizable random base insertion.

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Usage](#usage)
- [Examples](#examples)
- [Output Format](#output-format)
- [Modes](#modes)
- [Performance](#performance)
- [License](#license)

## Features

- Generate purely random DNA sequences or combine reference sequences with random bases
- Multi-processing support for efficient sequence generation
- Support for various sequence length formats (number, kb, mb)

## Prerequisites
- Python >=3.8
- BioPython

## Usage

```bash
python RandSeqGen.py [-h] [-v] -l LENGTH [-n NUMBER] [-b BATCH] [-p PROCESSOR] [-o OUTPUT] [-r REFERENCE] [--ratio RATIO] [--verbose]
```

### Program Information

- `-h, --help`
  - Show help message and exit
- `-v, --version`
  - Show version information and exit

### Required Arguments

- `-l, --length`
  - Length of each sequence (e.g., 100, 1kb, 1mb)

### Optional Arguments

- `-n, --number` (default: 1)
  - Number of sequences in each file

- `-b, --batch` (default: 1)
  - Number of generated FASTA files

- `-p, --processor` (default: available cores - 2)
  - Number of processor cores to use

- `-o, --output` (default: "RandSeqGen-Result")
  - Output directory path

- `-r, --reference`
  - Path to reference library directory
  - If not specified, the program will generate pure random sequences
  - Built-in reference libraries options:
    - `TIR`: Use built-in TIR reference library

- `--ratio` (default: 0.2)
  - Ratio of random bases in reference mode

- `--verbose`
  - Enable detailed progress output

## Examples

### Generate Random Sequences

```bash
# Generate 100 sequences of 1kb each in 2 batches
python RandSeqGen.py -l 1kb -n 100 -b 2
```

### Generate Reference-Based Sequences

```bash
# Generate sequences using built-in TIR reference library with 20% random bases
python RandSeqGen.py -l 10kb -n 100 -b 2 -r TIR
```

## Output Format

The program generates FASTA files in the specified output directory:
- Format: `sequences_batch_X.fa` where X is the batch number
- Each sequence has a unique identifier: `seq_Y_batch_X_len_Z`
  - Y: sequence number
  - X: batch number
  - Z: sequence length

## Modes

### Random Mode

Generates sequences using random choices from predefined 4 DNA nucleotides (A, T, G, C)

### Reference Mode
1. Loads reference sequences from specified directory
2. Generates random segment distribution
3. Constructs sequences with the following steps:
    1. Inserting reference segments at calculated positions
    2. Filling gaps with randomly generated DNA bases
    3. Ensuring total sequence length matches specified requirement

## Performance

- Parallel processing for every sequence generation and file I/O
- Binary search optimization for reference sequence selection
- Minimal disk I/O with in-memory processing

## License

This project is licensed under the LGPL-2.1 License - see the [LICENSE](LICENSE) file for details.
