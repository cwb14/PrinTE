# Changelog

Notable changes to PrinTE. Format follows [Keep a Changelog](https://keepachangelog.com/).

## [1.0.0]

First packaged release. The simulator itself is unchanged: the burn-in and generation loop
produce byte-identical output to the pre-packaging code for the same seed.

### Added

- Installable package: `pip install -e .` from a clone puts `printe`, `printe-grid`,
  `printe-score` and `printe-benchmark` on PATH. Not on PyPI or bioconda yet; the recipe
  is in `conda/`.
- Container images at `ghcr.io/cwb14/printe`, built and pushed on tag. Apptainer
  definition in `containers/` for clusters without a Docker daemon.
- Nextflow pipeline (`main.nf`) with `--mode simulate` and `--mode sweep`, and profiles for
  conda, Docker, Singularity, local, Slurm, and AWS Batch.
- Test suite: unit tests plus end-to-end runs on small fixtures in `tests/data/`.
- GitHub Actions running lint, tests on Linux and macOS, a container build, and the
  Nextflow test profile.
- `Makefile` carrying the `ltr_mutator` build recipes for both platforms.
- `--version` on `PrinTE.sh`.
- `--title` / bare `--title` on the plotting scripts.
- `docs/`, including a troubleshooting page and AWS Batch deployment notes.

### Changed

- `env.yml` is now `environment.yml`, and it gained `scikit-learn` and `minimap2`, which
  the grid search and `plot_indel.py` import but which were never declared. The optional R
  stack moved to `environment-r.yml`.
- `bin/`, `util/`, and `grid/` moved into an importable `printe` package under `src/`.
  `bash PrinTE.sh` works exactly as before from a clone.
- `ratios.tsv` and `ratios_ltr_only.tsv` moved to `src/printe/data/`.
- The R scripts moved to `R/`.
- `ltr_mutator` is compiled on first use instead of being shipped as a binary. The
  prebuilt Linux binary did not load on RHEL 8.
- Kmer2LTR is cloned into `~/.cache/printe/` rather than into the installation directory,
  so PrinTE works from a read-only install.
- Figures no longer carry titles by default, since journals strip them. Pass `--title` to
  restore the old text.
- `plot_category_bar.py`, `plot_superfamily_count.py`, and `genome_plot.py` accept input
  and output paths instead of hardcoding them. Defaults match the previous behaviour.

### Fixed

- `data/fasta_to_RepeatMasker.py` did not parse at all. An example alias file had been
  pasted into the header without comment markers, so the file raised `SyntaxError` on
  import while the README told people to run it.
- Starting from `--fasta`/`--bed` failed on the first generation, because the insertion
  step was always passed `-bf burnin.stat` even when no burn-in had run. This broke the
  documented Track 2 workflow. TE birth is now skipped in that case, with a note in the
  log.
- `build_composite_matrix.py` defaulted `--compare-script` to `compare_genomes2.py`, which
  has never existed in the repository.
- `TE_lib_stitcher.py` contained two complete programs concatenated into one file, so it
  ran both and always exited non-zero. Split into `TE_lib_stitcher.py` and
  `merge_ltr_parts.py`. Note that the former pairs on the first `LTR` or `I` anywhere in a
  header, so it silently skips names containing an I; prefer `merge_ltr_parts.py`.
- Thirteen plotting scripts imported `matplotlib.pyplot` without selecting a non-interactive
  backend, which fails on a headless node.
- `genome_plot.py` executed its whole body on import.
- A malformed weight in `ratios.tsv` silently disabled that TE family for the entire run.
  It now warns on stderr; the fallback value is unchanged.
- `plot_pipeline_report.R` reported saving a different filename than it wrote.

### Removed

- Five superseded modules that nothing referenced: the serial ancestors of the three
  simulator steps, `intact_LTR_extractor.py`, and `ltr_mutator_random_gen.cpp`.
- `grid/run_array.sh` and `grid/generate_params.sh`, which were hardcoded to a filesystem
  that no longer exists and were superseded by `guided_search.py`.
- `util/plot_TE.py`, `util/plot_benchmark_metrics.py`, `util/plot_divergence5.R`, and
  `util/pipeline_report_rate.py`.
- `grid/gridsearch.py`. Fixed-grid sampling is now the Nextflow sweep; `guided_search.py`
  remains the active-learning front end and is what the published analysis used.
- `data/ltr-db.fa.gz` moved to a release asset; `make fetch-data` retrieves it.
