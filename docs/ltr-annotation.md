## Annotating LTR-RTs and building a species-specific library

Annotating LTR-RTs in a genome, and building a species-specific LTR-RT
library so simulated insertions resemble your own elements.

[Back to the README](../README.md)

---

Two LTR-RT operations recur in PrinTE: **annotating** LTR-RTs in a genome (to benchmark your
calls and to score simulations), and **building a species-specific LTR-RT library** so that
simulated insertions resemble the real elements of your genome. The companion tools
[synLTR](https://github.com/cwb14/synLTR) (annotation) and
[Kmer2LTR](https://github.com/cwb14/Kmer2LTR) (library building + LTR dating) do this and benchmark well for PrinTE's purposes.

```bash
git clone https://github.com/cwb14/synLTR.git
git clone https://github.com/cwb14/Kmer2LTR.git    # also auto-cloned by PrinTE/synLTR when needed
```

**Annotate LTR-RTs in a genome** (real or simulated). `--max-rounds 1` does a single, non-nested
pass - the right setting for PrinTE; `--proteins` is a protein FASTA used to classify elements
(the bundled `data/TAIR10.pep.fa.gz` works for most plant genomes):

```bash
bash synLTR/module2/ltrharvest_wrapper2.sh \
    --genome mygenome.fa --proteins data/TAIR10.pep.fa.gz \
    --max-rounds 1 --threads 64 --out_prefix mygenome
```

Two of its outputs matter here (with `--out_prefix mygenome`):

- **`mygenome_r1_ltr.tsv`** - the LTR-RT table. **Column 1 is the `chrom:start-end` locus and
  column 7 is the LTR–LTR p-distance** (the age proxy). Pass this to `bedtools.py -pass_scn` for
  [benchmarking](#benchmark-your-te-annotations-first) and to the grid scorer (as `--ref-tsv` and
  the per-simulation `--exp-tsv-name`).
- **`mygenome_depth0_ltr.fa`** - a FASTA of the (unnested) LTR-RT sequences, with
  RepeatMasker-style headers `>chrom:start-end#LTR/superfamily/...`.

(synLTR needs `mafft` and `trimal` on PATH and clones Kmer2LTR on first run, so that run needs
network access.)

**Build a species-specific library** from that FASTA. `Kmer2LTR --make-perfect-ltr-rt 5p`
rebuilds each element with two *identical* LTRs (as at the moment of insertion) and appends the
`~LTRlen:<N>` tag PrinTE needs. `5p` (use the 5′ LTR) is the natural choice:

```bash
python Kmer2LTR/Kmer2LTR.py -i mygenome_depth0_ltr.fa --make-perfect-ltr-rt 5p
# -> mygenome_depth0_ltr.LTRs.alns.perfect_5p.fa
```

Feed that file to PrinTE as a ready-made library: its headers already use the
`name#class/superfamily~LTRlen:<N>` convention, so PrinTE uses them as-is (skips library processing):

```bash
printe ... --clean_lib mygenome_depth0_ltr.LTRs.alns.perfect_5p.fa
```

**Why this matters biologically:** the simulation now inserts the *actual* LTR-RTs of your species
 (their real sequences, length distribution, and superfamily mix) so the simulated TE landscape,
and the genome-size trajectory it drives, reflect your genome's own elements rather than a generic
library. (`--make-perfect-ltr-rt` **requires** a mode; a bare `--make-perfect-ltr-rt` errors.
Kmer2LTR only appends `~LTRlen:N`; the `#class/superfamily` comes from synLTR's headers, so the
chain works because synLTR emits RepeatMasker-style names.)

---
