# Contributing to PrinTE

Thanks for your interest. Bug reports and pull requests are both welcome.

## Reporting a bug

Open an issue using the bug template. The single most useful thing you can include is the
tail of `pipeline.error` from the failed run: on success that file holds only a start
banner, so anything else in it is the real failure.

## Setting up a development environment

```bash
git clone https://github.com/cwb14/PrinTE.git && cd PrinTE
mamba env create -f environment.yml && conda activate PrinTE
pip install -e '.[dev]'
pre-commit install
make ltr-mutator
```

Then:

```bash
make lint        # ruff + shellcheck, both blocking in CI
make test-fast   # unit tests, seconds
make test        # adds the end-to-end runs, a few minutes
```

## Tests

Fixtures live in `tests/data/` and are tiny synthetic files - a 21-sequence TE library, 100
CDS, a small GFF. Keep them that way: no real genomes, and nothing that makes the suite
take longer than a coffee. `tests/data/README.md` records how each was generated so they
can be rebuilt.

Tests that invoke the pipeline are marked `slow`. `make test-fast` skips them, CI runs both
on Linux and only the fast ones on macOS.

## Changing simulation behaviour

**A green test suite is not sufficient.** The tests check that outputs are well formed, not
that they are scientifically correct. If your change touches the simulator, include a
before-and-after comparison on a fixed seed in the pull request:

```bash
# on main, then on your branch, in separate directories
printe -sz 300kb -cn 2 -P 5 -itp 12 -ftp 5 -cl tests/data/tiny_TE.lib \
       -c tests/data/tiny_cds.fa -r tests/data/tiny_ratios.tsv \
       -ir 1e-5 -dr 1e-5 -ge 200 -st 100 -t 4 -s 7 --no_postproc
md5sum gen200_final.fasta gen200_final.bed
```

If the hashes change, say so and explain why. If they should not have changed, that is a
bug in the change.

## Design notes

Two decisions that look like oversights but are not. Please do not "fix" them without
raising an issue first.

**`PrinTE.sh` does not use `set -e`.** It handles errors explicitly with `if [ $? -ne 0 ]`
after each step, which lets it print a useful message naming the step that failed. Under
`set -e` the shell would exit before those handlers ran, and the failure would become
silent. Two shellcheck rules are disabled at the top of the file for the same reason:
SC2181 flags that idiom, and SC2086 flags the deliberate word splitting inside the `eval`
command strings.

**`ltr_mutator` is compiled on first use, not shipped.** A prebuilt binary in the
repository breaks the moment someone's libc differs, which is exactly what happened before
1.0.0. `ensure_mutator` in `PrinTE.sh` builds it if it is missing or will not load, and the
`Makefile` carries the Linux and macOS recipes.

## Style

`ruff` with `E,F,I,B,UP` is the whole style guide; run `make lint`. Beyond that, match the
file you are editing rather than any external convention. The codebase does not use
`logging`, does not annotate types everywhere, and keeps domain acronyms capitalised in
module names (`TE_exciser.py`). That is deliberate and consistent.
