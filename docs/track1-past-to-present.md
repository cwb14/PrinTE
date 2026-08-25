## Track 1 - Past → present

Reconstruct a clade's common ancestor, build it as a burn-in, and evolve it
forward until the simulation looks like the real present-day genome.

[Back to the README](../README.md)

---

**Goal:** reconstruct the genome and TE state of a clade's common ancestor, build it as a
burn-in, and evolve it forward until the simulation resembles the **real** present-day genome.
Because the real genome is ground truth, you can score the fit and ask whether TE activity
**increased or decreased** along the branch.

You should already have a **time-scaled Newick tree** (`tree.nwk`) and **each tip's genome
FASTA**. Run [synLTR](#annotating-ltr-rts-and-building-a-species-specific-library) on each tip's
genome with `--out_prefix <tip>` to produce, per tip, an LTR-RT table `<tip>_r1_ltr.tsv` (column
7 = p-distance/age) and an LTR-RT FASTA `<tip>_depth0_ltr.fa`. So each tip has three files (
genome, table, LTR-RT FASTA) which drive the three reconstructions below. (Distance trees can be
time-scaled with a tool such as [PATHd8](https://www2.math.su.se/PATHd8/).)

### Step 1 - reconstruct the ancestor

Two R scripts run `phytools::fastAnc` over the tree. **File-naming/placement matters:**

- `ancestral_reconstruction_gs.R` looks for each tip's file as `<tip_label><suffix>` **in the
  same directory as the Newick** (`<tip>` is the exact tip string in the tree).
- `ancestral_reconstruction_ltr_age.R` looks for `<tip_label><suffix>` **in the current working
  directory** - so run it from where those files live.

```bash
# (a) Ancestral GENOME size  -> needs one genome FASTA per tip, named <tip>.fa, next to tree.nwk.
Rscript R/ancestral_reconstruction_gs.R \
    --newick tree.nwk --suffix .fa --out ancestral_genome_size --label_nodes

# (b) Ancestral LTR-RT size  -> needs each tip's LTR-RT FASTA (synLTR's <tip>_depth0_ltr.fa)
#     next to tree.nwk; the script sums their bp.
Rscript R/ancestral_reconstruction_gs.R \
    --newick tree.nwk --suffix _depth0_ltr.fa --out ancestral_ltrrt_size --label_nodes

# (c) Ancestral LTR-RT AGE distribution  -> run from the dir holding each tip's LTR table
#     (synLTR's <tip>_r1_ltr.tsv; TAB-separated, col 7 = p-distance/age in 0-1).
Rscript R/ancestral_reconstruction_ltr_age.R \
    --newick tree.nwk --suffix _r1_ltr.tsv
```

(`gs.R` needs `samtools` on PATH and builds `<file>.fai` automatically. `ltr_age.R` calls
`grid/summary_stats.py` via `python` - keep it beside the R script, or pass `--summary_path` -
and bins ages with defaults `--bins 50 --bin_max 0.15`.)

**Find your node.** `gs.R` writes `ancestral_genome_size.contMap.pdf` (node numbers show only
with `--label_nodes`) and `ancestral_genome_size.ancestral_genome_sizes.tsv`. `ltr_age.R` writes
`tree_with_node_ids.pdf` and `node_map.tsv` (each internal node → its full set of descendant
tips). Open a PDF (or `grep` your tip names in `node_map.tsv`), find the internal node whose
descendants are exactly your clade, and **use that node number consistently below** (here, `20`).

### Step 2 - turn the reconstruction into PrinTE inputs

```bash
# Ancestral genome size and intact-LTR-RT size are column 6 (bp_est) of each gs.R table.
gs=$(awk -F'\t' '$1==20{print $6}' ancestral_genome_size.ancestral_genome_sizes.tsv)
ltr=$(awk -F'\t' '$1==20{print $6}' ancestral_ltrrt_size.ancestral_genome_sizes.tsv)

# -sz must carry an Mb suffix; -itp is the intact-LTR-RT percent of the genome.
sz=$(awk -v g="$gs" 'BEGIN{printf "%.2fMb", g/1e6}')          # e.g. 180.27Mb
itp=$(awk -v g="$gs" -v l="$ltr" 'BEGIN{printf "%.3f", 100*l/g}')  # e.g. 16.629
echo "-sz $sz   -itp $itp"

# Build the age-distribution file for -mb from ltr_age.R's ancestral_summary_bins.tsv.
# Columns there are: stat_type, node, clade_size, clade_signature, descendant_tips,
#                    bin_index, bin_start, bin_end, est, CI_lower, CI_upper.
# Pull your node's bins (start,end,count = cols 7,8,9) and round the counts:
awk -F'\t' '$2==20' ancestral_summary_bins.tsv | cut -f7-9 \
  | awk -F'\t' 'BEGIN{OFS="\t"}{ $3=sprintf("%.0f",$3); print }' > ancestral.LTR.bins

# Normalize counts to frequencies -> 3-column "bin_start bin_end frequency" file that -mb expects:
awk -v OFS='\t' '{a[NR]=$0; v[NR]=$3; s+=$3} END{for(i=1;i<=NR;i++) print a[i], v[i]/s}' \
  ancestral.LTR.bins | cut -f1,2,4 > ancestral.LTR.bins.freq
```

### Step 3 - build the ancestor as a burn-in

```bash
# the bundled ratios file, wherever PrinTE was installed
RATIOS=$(python -m printe.paths ratios_ltr_only.tsv)

printe --burnin_only \
    -sz "$sz" -itp "$itp" -ftp 50 \
    -mb ancestral.LTR.bins.freq \
    -m 2e-9 -TsTv 2.0 --ex_LTR \
    -cl lib_clean.fa -r "$RATIOS" -t 200
```

`-sz`, `-itp`, and `-mb` come from Step 2. `-ftp` is the **fragmented**-TE percent - it is *not*
reconstructed above (the gs.R/ltr_age.R steps cover intact LTR-RTs only), so estimate it from
your clade's present-day fragment + solo load, or use it as filler/deletion substrate. `-cl
lib_clean.fa` is the cleaned TE library; for a more realistic ancestor, build a
**species-specific** LTR-RT library from a close present-day relative
(see [Annotating LTR-RTs](#annotating-ltr-rts-and-building-a-species-specific-library)) and pass it here.

### Step 4 - benchmark, evolve, and fit

Optionally [benchmark](#benchmark-your-te-annotations-first) the burn-in (run your LTR pipeline on
`burnin.fasta`; a high recall on a genome with a known TE landscape means later mismatches reflect
**biology**, not annotation error). Then use the [grid search](#grid-search-finding-parameters-without-bias)
to find the insertion/deletion rates (and `sr`, `k`) whose terminal genome best matches the real
one.

**Interpret.** The best-fit rates are your estimate of TE turnover along the branch. Compare them
to the ancestral and present-day TE loads to decide whether activity rose or fell - and use the
ancestral reconstruction to motivate which hypothesis (increasing vs decreasing rate) to test.

---
