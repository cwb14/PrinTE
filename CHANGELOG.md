# Changelog

Notable changes to PrinTE. Format follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

## [1.0.3] - 2026-09-29

1.0.1 and 1.0.2 were tagged without changelog entries or version bumps; their changes are
listed here.

### Changed

- Post-processing is off by default. Pass `--postproc` (`-pp`) to date LTR-RTs and make the
  plots and reports. `--no_postproc` still works and is now the default.
- `--birth_rate` defaults to 1e-8 instead of 1e-3.
- One environment: the R stack moved from `environment-r.yml`, now removed, into
  `environment.yml`. Most users need the R scripts, and `ancestral_reconstruction_ltr_age.R`
  calls PrinTE's Python anyway. Update an existing env with
  `mamba env update -n PrinTE -f environment.yml --prune`.
- Nothing is written to `~/.cache` any more. `make fetch-data` downloads `ltr-db.fa.gz` into
  `data/`, Kmer2LTR is cloned into the PrinTE directory, and `ltr_mutator` is built into its
  `bin/`. `PRINTE_CACHE`, when set, holds `bin/` and `Kmer2LTR/` instead.
- The container carries Kmer2LTR and the R scripts (`/opt/printe/R`), so post-processing
  needs no network and clones nothing at runtime.
- Kmer2LTR is fetched when a run starts rather than when post-processing begins, so a missing
  network or an unwritable directory fails at once instead of after the simulation. Runs
  started together, such as a job array or a sweep, share one clone.
- A Kmer2LTR that fails, or a clone at another commit than the pinned one, now prints a
  warning instead of passing silently.

### Fixed

- On a fresh setup, post-processing dated no LTR-RTs yet exited 0: it cloned Kmer2LTR's
  `main` branch, a rewrite without the `Kmer2LTR.py` that PrinTE runs. Kmer2LTR is now
  pinned to a commit on its `legacy` branch.
- `ltr_dens.py` misread Kmer2LTR's 16-column output: pandas made the four surplus columns
  the index and shifted every name four places, so `all_LTR_density.pdf` plotted LTR
  coordinates in place of divergences. It now reads the first 12 columns by position.
- `ltr_mutator` built through `make` always landed in PrinTE's own `bin/`, even when
  `PRINTE_MUTATOR_DIR` said to put it elsewhere, so the run could not find it.
- An `ltr_mutator` that no longer runs, such as one built on a machine with another C library,
  is rebuilt; `make` used to call it up to date and the run failed at the first generation.
- `--version` reports the release; 1.0.1 and 1.0.2 still printed `PrinTE 1.0.0`.

### Added

- `manuscript_figures/`: the paper's seven figures and the scripts that drew them.

### Upgrading from 1.0.0

PrinTE no longer looks in `~/.cache/printe`, so what 1.0.0 put there is ignored:

- `ltr-db.fa.gz`: move it into `data/`, or point `PRINTE_DATA` at its directory.
- A Kmer2LTR clone: to keep using it, `export PRINTE_CACHE=~/.cache/printe`. If it is a clone
  of `main`, PrinTE says so and prints the `git checkout` that fixes it.
- A read-only install that relied on the cache fallback now stops at startup until
  `PRINTE_CACHE` names a writable directory.

Otherwise `rm -r ~/.cache/printe`.

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
- The container could not be run as a command. `printe` was the Docker `CMD`, and
  arguments replace the CMD, so `./printe.sif --version` or `docker run printe --version`
  dropped `printe` and tried to exec `--version`. It is part of the `ENTRYPOINT` now, so
  the image works as `./printe.sif <options>`.
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
