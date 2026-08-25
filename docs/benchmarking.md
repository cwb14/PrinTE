## Benchmark your TE annotations first

How to measure your TE-annotation pipeline against a genome whose
TE positions are known exactly. Do this once before either track.

[Back to the README](../README.md)

---

PrinTE simulates from TE annotations, so **the forward framework assumes high-quality
annotations**. Conveniently, PrinTE is also the tool to *test* that assumption because it knows the
exact location of every TE it places, so a simulated genome is a ground-truth benchmark for any
annotation pipeline. **Do this once before either track.**

1. **Build a genome** whose TE landscape mimics your real one. `--burnin_only` writes
   `burnin.fasta` (genome) and `burnin.bed` (the ground truth) to the current directory:

   ```bash
   printe --burnin_only -sz 100Mb -cn 1 -P 4 -itp 10
   ```

2. **Annotate the LTR-RTs** in `burnin.fasta`. We recommend
   [synLTR](#annotating-ltr-rts-and-building-a-species-specific-library) (below): with
   `--out_prefix burnin` it writes `burnin_r1_ltr.tsv`, whose column 1 is the `chrom:start-end`
   locus - exactly what `-pass_scn` reads. Any LTR finder works as long as you hand `bedtools.py`
   an SCN or PASS file (an SCN line has ≥12 whitespace columns with start/end in columns 1–2 and
   the sequence name last; a PASS line starts `chr:start..end`). Use the **same pipeline on your
   real genomes** in Track 1/2.

3. **Compare predictions to the truth** with `util/bedtools.py`:

   ```bash
   python -m printe.util.bedtools \
       -pass_scn burnin_r1_ltr.tsv \
       -bed burnin.bed --printe-intact \
       -r 0.9
   ```

   `--printe-intact` subsets the truth BED to intact LTR-RTs so false negatives stay meaningful;
   `-r 0.9` requires reciprocal 90% overlap (each call and the true element cover ≥90% of the
   other - lower it, e.g. `-r 0.5`, for looser positional agreement). Expected output:

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
