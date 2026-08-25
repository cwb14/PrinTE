# Bundled data

## What PrinTE uses by default

| file | used as | size |
|---|---|---|
| `maize_rice_arab_curated_TE.lib.gz` | default `--TE_lib`, the elements new insertions are drawn from | 4.4 MB |
| `TAIR10.cds.fa.gz` | default `--cds`, the genes placed in a burn-in genome. 19,621 sequences | 5.2 MB |
| `TAIR10.pep.fa.gz` | protein set for classifying LTR-RTs during annotation. Works for most plant genomes | 5.1 MB |

## Source libraries

Intact LTR-RT exemplars, from Shujun Ou, 2025-07-23:

| file | source |
|---|---|
| `athrep.updated.fasta_062024.ltr.full.gz` | [oushujun/TAIR12-TE](https://github.com/oushujun/TAIR12-TE/tree/main/data/TE_library) |
| `rice7.0.0.liban.ltr.full.gz` | [oushujun/riceTElib](https://github.com/oushujun/riceTElib) |
| `maizeTE02052020.ltr.full.gz` | [oushujun/MTEC](https://github.com/oushujun/MTEC) |

The rice and Arabidopsis libraries ship partitioned into `Internal` and `LTR` entries.
They were reconstituted here into intact LTR-RTs as `[LTR][INTERNAL][LTR]`; see
`printe.util.merge_ltr_parts` for the tool that does this. The maize library is already
intact on download and needed no reconstitution.

These three are kept for provenance. Nothing reads them automatically; they are the inputs
that `maize_rice_arab_curated_TE.lib.gz` was built from.

## Not in this directory

`ltr-db.fa.gz`, the LTR-RT exemplar database, is 38 MB and no code path reads it, so it
ships as a release asset rather than a repository file:

```bash
make fetch-data          # downloads it into ~/.cache/printe/data
```

`printe.paths` finds it there, or wherever `$PRINTE_DATA` points.
