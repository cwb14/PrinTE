## Preparing a custom TE library

Assembling a TE library by hand, and the helpers that normalise headers and
sequence characters.

[Back to the README](../README.md)

---

PrinTE ships with a curated maize/rice/Arabidopsis library
(`maize_rice_arab_curated_TE.lib.gz`, the default `--TE_lib`) and, as a release asset,
an LTR-RT exemplar database (`ltr-db.fa.gz` - fetch it with `make fetch-data`).
For LTR-RT studies, the best library is usually a **species-specific** one built from your own
genome - see [Annotating LTR-RTs](#annotating-ltr-rts-and-building-a-species-specific-library).
To assemble a library by hand instead, a few helpers are provided:

- **`util/TE_lib_stitcher.py`** - reassemble split `LTR` + `Internal` library entries into
  intact LTR-RTs.
- **`data/fasta_to_RepeatMasker.py`** - normalize heterogeneous headers into
  `>name#class/superfamily` form.
- **`data/fix_non_ATGC.py`** - replace non-ACGT characters with `T` so downstream k-mer/alignment
  steps don't choke (streams large files; writes to stdout):

  ```bash
  zcat custom.lib.gz | python -m printe.util.fix_non_ATGC - > my_TE.lib
  ```

Provide the result via `--TE_lib my_TE.lib` (PrinTE processes it and writes a cleaned
`lib_clean.fa` you can reuse) or, if already cleaned/pre-processed, via `--clean_lib lib_clean.fa`
to skip library processing.

---
