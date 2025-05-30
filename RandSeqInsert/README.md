# RandSeqInsert

[![license](https://img.shields.io/github/license/lutianyu2001/RandSeqInsert.svg)](https://github.com/lutianyu2001/RandSeqInsert/blob/master/LICENSE)

RandSeqInsert is a high-performance Python tool for simulating transposable element insertions in genomic sequences. Built around an AVL tree-based algorithm with event sourcing architecture, it enables precise modeling of complex nested insertions and provides comprehensive sequence reconstruction capabilities.

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Core Features](#core-features)
- [Examples](#examples)
- [Output Formats](#output-formats)
- [Advanced Features](#advanced-features)
- [Performance](#performance)
- [License](#license)

## Features

- **AVL Tree-Based Architecture**: Efficient O(log n) insertion operations with automatic tree balancing
- **Event Sourcing System**: Selective tracking of nested insertion events for complex scenario reconstruction
- **Target Site Duplication (TSD) Modeling**: Biologically accurate TSD generation with configurable mutations
- **Nested Insertion Support**: Simulation of donor-to-donor insertions with fragment tracking
- **Sequence Reconstruction**: Three reconstruction modes (full, clean, event history)
- **Dual Visualization**: Tree structure and event graph visualizations in Graphviz DOT format
- **Multiple Output Formats**: FASTA, BED, and specialized reconstruction files
- **High-Performance Processing**: Multi-core support with memory-efficient operations
- **Flexible Donor Libraries**: Support for custom and built-in transposon libraries

## Prerequisites

- Python ≥3.8
- BioPython
- NumPy
- (Optional) Graphviz for visualization rendering

## Installation

```bash
git clone https://github.com/lutianyu2001/RandSeqInsert.git
cd RandSeqInsert
pip install -r requirements.txt
```

## Usage

```bash
python RandSeqInsert.py [-h] [-v] -i INPUT -is INSERTION [-it ITERATION] [-b BATCH] 
                        [-p PROCESSORS] [-o OUTPUT] [-d DONOR [DONOR ...]] 
                        [-w WEIGHT [WEIGHT ...]] [-l LIMIT] [--seed SEED]
                        [--tsd TSD_LENGTH] [--track] [--visual] [--recursive] 
                        [--filter_n] [--debug]
```

### Core Arguments

- `-i, --input` **[Required]**
    - Input sequence file in FASTA format
- `-is, --insert` **[Required]**
    - Number of insertions per sequence (supports 1k, 1m notation)
- `-it, --iteration` (default: 1)
    - Number of insertion iterations per sequence
- `-b, --batch` (default: 1)
    - Number of independent result files to generate
- `-p, --processors` (default: CPU cores - 2)
    - Number of processors for parallel processing
- `-o, --output` (default: "RandSeqInsert-Result")
    - Output directory path

### Donor Library Arguments

- `-d, --donor` **[Required]**
    - Donor sequence library file(s) or directory paths
    - Multiple libraries supported
    - Built-in libraries: `TIR/rice`, `TIR/maize`
- `-w, --weight`
    - Weights for donor libraries (must match number of libraries)
- `-l, --limit`
    - Maximum donor sequence length to load

### Feature Flags

- `--tsd TSD_LENGTH`
    - Enable Target Site Duplication with specified length
- `--track`
    - Track and save used donor sequences with reconstruction
- `--visual`
    - Generate Graphviz DOT files for tree and event visualization
- `--recursive`
    - Use recursive insertion method (default: iterative)
- `--filter_n`
    - Filter out donor sequences containing N bases
- `--debug`
    - Enable debug mode with detailed information
- `--seed SEED`
    - Random seed for reproducible results

## Core Features

### AVL Tree Architecture

RandSeqInsert uses a balanced binary tree structure for efficient sequence manipulation:
- **O(log n) insertion complexity** maintaining performance for large sequences
- **Automatic balancing** through tree rotations
- **Memory efficient** with node-based sequence representation

### Event Sourcing for Nested Insertions

Advanced tracking system for complex insertion scenarios:
- **Selective recording** of only nested (donor-to-donor) insertions
- **Complete reconstruction** of fragmented donor sequences
- **Event history** preservation for temporal analysis

### Target Site Duplication (TSD) Modeling

Biologically accurate simulation of insertion signatures:
- **Configurable TSD length** based on transposon type
- **Independent 5' and 3' mutations** with SNP and InDel support
- **Realistic mutation rates** for authentic simulation

## Examples

### Basic Insertion Simulation

```bash
# Insert 100 TIR elements into genome sequences
python RandSeqInsert.py -i genome.fa -is 100 -d TIR/maize
```

### Complex Nested Insertion with TSD

```bash
# Simulate nested insertions with TSD and tracking
python RandSeqInsert.py -i genome.fa -is 50 -it 3 \
  -d TIR/maize -d TIR/rice -w 0.7 -w 0.3 \
  --tsd 9 --track --visual
```

### Multiple Iterations and Batches

```bash
# Generate 5 independent datasets with multiple iterations
python RandSeqInsert.py -i genome.fa -is 20 -it 5 -b 5 \
  -d custom_library.fa --track --seed 12345
```

### Benchmarking Setup

```bash
# Create ground truth dataset for annotation tool benchmarking
python RandSeqInsert.py -i reference.fa -is 1000 \
  -d comprehensive_TE_lib.fa --tsd 5 --track --visual \
  --filter_n --debug -o benchmark_dataset
```

### Large-Scale Simulation

```bash
# High-throughput simulation with multi-processing
python RandSeqInsert.py -i large_genome.fa -is 5k -it 2 -b 10 \
  -p 16 -d TIR/maize -d TIR/rice -w 0.6 -w 0.4 \
  --tsd 7 --track --recursive
```

## Output Formats

### Primary Outputs

- **`sequences_batch_X.fa`**: Modified sequences with insertions
- **`used_donors_batch_X.fa`**: Active donor sequences (if `--track` enabled)
- **`reconstructed_donors_batch_X.fa`**: Full reconstructed sequences
- **`clean_reconstructed_donors_batch_X.fa`**: Clean reconstructed sequences
- **`donors_batch_X.bed`**: BED format annotation with insertion coordinates

### Visualization Outputs (with `--visual`)

- **`visualization/seqid_tree_visual.dot`**: Sequence tree structure
- **`visualization/seqid_event_visual.dot`**: Insertion event relationships

### BED File Format

```
chr1    1000    1500    donor_123;TIR_element    ATCG    +
chr1    2000    2300    donor_456;LTR_element    GCTA    +
```

Columns: chromosome, start, end, name, TSD_sequence, strand

## Advanced Features

### Sequence Reconstruction Modes

1. **Full Reconstruction**: Complete sequences including all nested content
2. **Clean Reconstruction**: Original donor sequences with nested elements removed  
3. **Event History**: Step-by-step sequence states through insertion events

### Event Sourcing Architecture

- **Selective tracking** reduces memory overhead
- **Complete reconstruction** of complex nested scenarios
- **Temporal analysis** capabilities for evolutionary studies

### Dual Visualization System

- **Tree Visualization**: Hierarchical structure with position annotations
- **Event Graph**: Temporal relationships and nesting patterns
- **Interactive exploration** of complex insertion scenarios

### Performance Optimizations

- **Multi-core processing** for large-scale simulations
- **Memory-efficient** tree-based operations
- **Incremental balancing** maintains performance
- **Chunked processing** for very large sequences

## Performance

### Complexity
- **Insertion**: O(log n) per operation
- **Balancing**: Automatic with O(log n) overhead
- **Memory**: Linear with sequence length plus tree structure

### Benchmarks

### Scalability

## Use Cases

### Method Development
- **Benchmarking** TE annotation tools
- **Ground truth generation** with known insertion histories
- **Algorithm validation** for structural variant detection

### Evolutionary Studies
- **Genome size evolution** modeling
- **TE accumulation** simulation over time
- **Nested insertion** impact analysis

### Population Genomics
- **Pan-genome** TE variation modeling
- **Insertion polymorphism** simulation
- **Population-specific** TE landscapes
