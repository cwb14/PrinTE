## Track 2 - Present → future

Take a replica of a present-day genome and scan for the parameter regimes
that keep it stable over a short horizon.

[Back to the README](../README.md)

---

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
python -m printe.util.gff_to_bed TEgenomeSimulator.gff TEgenomeSimulator.bed
# Sanity check: one short evolve step should run cleanly.
printe -f TEgenomeSimulator.fa -b TEgenomeSimulator.bed -i my_TE.lib -ge 100000 -st 100000 -t 20
```

**3. Scan for the stability regime.** Use the same [grid search](#grid-search-finding-parameters-without-bias),
but with the **replica genome as the reference**, a **tight size band** around the current size,
and a **short horizon**. For example, with a current genome of ~180,000,000 bp:

```bash
python -m printe.grid.guided_search init --samples 100 --seed 7 \
    --ins-start 1e-10 --ins-end 1e-13 --del-start 1e-10 --del-end 1e-13 \
    --sr-start 5 --sr-end 95 --sr-step 10 --k-start 0 --k-end 10 --k-step 2 \
    --printe-script ./PrinTE.sh --ge 1000000 --st 100000 --mut 1.3e-8 --tstv 2.0 \
    --mxgs 181800000 --mngs 178200000 \
    --bed TEgenomeSimulator.bed --fasta TEgenomeSimulator.fa \
    --clean-lib lib_clean.fa --ratios src/printe/data/ratios.tsv --threads 64
```

A run **"held steady"** if its genome size and TE content stayed inside your tolerance band over
the horizon. The lightest-weight read-out is each run's `genome_size_plot.pdf` and `percent_TE.pdf`
(flat = stable); for a single number, score each run against the replica with
`build_composite_matrix.py` (low Composite = stayed close to the start). Because there is no ground
truth, treat Track 2 as **hypothesis generation** about genome stability - report *which rate
regimes maintain the genome* and which let it expand or erode, not a single historical rate.

---
