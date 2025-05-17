# RandSeqInsert

[![license](https://img.shields.io/github/license/lutianyu2001/RandSeqInsert.svg)](https://github.com/lutianyu2001/RandSeqInsert/blob/master/LICENSE)

RandSeqInsert is a high-performance Python tool for inserting random DNA segments from reference libraries into existing sequences. It allows for customizable insertion of references with precise control over insertion parameters.

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

- Insert reference sequences from libraries into existing DNA sequences
- Control insertion frequency and distribution
- Multi-processing support for efficient sequence processing
- Support for various sequence length formats (number, kb, mb)
- Tracking of inserted references (optional)

## Prerequisites
- Python >=3.8
- BioPython

## Usage

```bash
python RandSeqInsert.py [-h] [-v] -i INPUT -ins INSERTION [-it ITERATION] [-b BATCH] [-p PROCESSOR] [-o OUTPUT] [-r REFERENCE] [-w WEIGHT] [--track] [--filter_n] [--verbose]
```

### Program Information

- `-h, --help`
    - Show help message and exit
- `-v, --version`
    - Show version information and exit

### Arguments

- `-i, --input`
    - Path to input sequence file
- `-is, --insertion`
    - Number of insertions to perform in each sequence

- `-it, --iteration` (default: 1)
    - Number of iterations to process each sequence

- `-b, --batch` (default: 1)
    - Number of batches to process

- `-p, --processor` (default: available cores - 2)
    - Number of processor cores to use

- `-o, --output` (default: "RandSeqInsert-Result")
    - Output directory path

- `-r, --reference`
    - Path to reference library directory
    - Can be specified multiple times for different libraries
    - Built-in reference libraries options:
        - `TIR`: Use built-in TIR reference library

- `-w, --weight` (default: equal weights)
    - Weight for each reference library
    - Must match number of reference libraries

- `--track`
    - Enable tracking of inserted references

- `--filter_n`
    - Avoid using reference sequences with N

## Examples

### Basic Usage

```sh
# Insert 10 references from the built-in TIR/maize library into sequences in input.fa
python RandSeqInsert.py -i input.fa -is 10 -r lib/TIR/maize
```

#### Track References

```sh
# Enable tracking of inserted references
python RandSeqInsert.py -i input.fa -is 10 -r lib/TIR/maize --track --tsd 3
```

```sh
# Enable tracking of inserted references
python RandSeqInsert.py -i test_all_dna_nt_5kb.fa -is 10 -r lib/TIR/maize --track --tsd 3 --visual
```

### Advanced Examples

#### Multiple Reference Libraries with Weights

```sh
# Insert 10 references per sequence using two libraries with different weights
python RandSeqInsert.py -i input.fa -is 10 -it 10 -r lib/TIR/maize -w 0.8 -r lib/TIR/rice -w 0.2 --track
```

#### Multiple Iterations

```sh
# Process each sequence 3 times with 10 insertions each time
python RandSeqInsert.py -i input.fa -is 10 -it 3 -r lib/TIR/maize
```

#### Filter References with N

```sh
# Filter out reference sequences containing N
python RandSeqInsert.py -i input.fa -is 10 -r lib/TIR/maize --filter_n
```

#### Multiple Batches with Multi-processing

```sh
# Process in 5 batches using 8 processor cores
python RandSeqInsert.py -i input.fa -is 10 -b 5 -p 8 -r lib/TIR/maize
```

#### Comprehensive Example

```sh
# Complex example with multiple libraries, tracking, and verbose output
python RandSeqInsert.py -i input.fa -is 20 -it 2 -b 3 -p 12 -o custom_output \
  -r lib/TIR/maize -w 0.6 -r lib/TIR/rice -w 0.4 \
  --track --filter_n
```

## Output Format

The program generates FASTA files in the specified output directory:
- Format: `sequences_batch_X.fa` where X is the batch number
- If tracking is enabled: `references_batch_X.fa` with details of inserted references
- Each sequence retains its original ID with additional metadata

### Tracking Output Format

When `--track` is enabled, a separate file is generated with details about each inserted reference:
- Position in the original sequence
- Source reference library
- Reference sequence ID
- Length of the inserted reference

## Modes

### Standard Mode

Inserts reference sequences at random positions in the input sequences

### Tracking Mode

Same as Standard Mode but with additional tracking of inserted references

## Performance

- Parallel processing for efficient sequence processing
- Binary search optimization for reference sequence selection
- Minimal disk I/O with in-memory processing

## License

[TODO]
