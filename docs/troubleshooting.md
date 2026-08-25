## Troubleshooting

Common problems, in roughly the order people hit them.

[Back to the README](../README.md)

---

### A run failed. Where is the error?

`pipeline.error`, in the directory you ran from. On a successful run it holds only a
start banner, so anything else in it is the actual failure. `pipeline.log` has the full
transcript including every command PrinTE ran, which is what you want when a step
produced the wrong answer rather than no answer.

```bash
tail -40 pipeline.error
grep 'Running:' pipeline.log | tail -5     # the last command that started
```

### How much memory do I need?

Every genome-touching step loads the whole genome into memory, and the excision step
hands a copy to each worker process. Budget roughly **four times your genome size**:

| genome | comfortable |
|---|---|
| 150 Mb (Arabidopsis-scale) | 4 GB |
| 400 Mb | 8 GB |
| 1 Gb | 8-16 GB |
| 3 Gb (maize-scale) | 16-32 GB |

More threads (`-t`) means more concurrent chromosome copies, so on a memory-tight node
lower `-t` before you lower anything else. Under Nextflow the `process_high` label
retries with double the memory, which covers most first-attempt OOMs.

### `printe` runs the wrong PrinTE, or fails with `apptainer`/`run-singularity` not found

You installed PrinTE more than one way. The container route's optional launcher lives in
`~/.local/bin/printe` and that directory sits at the front of your `PATH`, so it wins over
the `printe` that `pip install -e .` puts in your conda environment. Depending on whether
Apptainer is on your `PATH` you will either get an error like

```
printe: line 2: exec: apptainer: not found
/usr/bin/env: 'run-singularity': No such file or directory
```

or, more confusingly, no error at all - just the container's PrinTE rather than the source
copy you meant to run.

Find out which one you are getting:

```bash
command -v printe
```

If that prints `~/.local/bin/printe` and you meant to use the source or conda install,
delete the launcher:

```bash
rm ~/.local/bin/printe
```

If you did mean to use the container, the message means Apptainer is not available in this
shell - load the module that provides it, or `mamba install -y -c conda-forge apptainer`.
You can always skip the launcher entirely and run the image directly as
`/full/path/to/printe.sif --burnin_only ...`.

### ltr_mutator will not build

PrinTE compiles the point-mutation binary on first use rather than shipping one, so a
clone works on any libc. It needs a C++17 compiler with OpenMP.

```bash
mamba install cxx-compiler          # Linux, no root needed
brew install libomp                 # macOS, alongside Xcode command line tools
make ltr-mutator                    # build it explicitly
```

If you already have a binary somewhere, point at it and skip the build:

```bash
export PRINTE_MUTATOR=/path/to/ltr_mutator
```

### Can I run without network access?

The simulation itself, yes. Post-processing clones
[Kmer2LTR](https://github.com/cwb14/Kmer2LTR) to date the LTR-RTs, which needs network
the first time. Either pre-clone it:

```bash
mkdir -p ~/.cache/printe && git clone https://github.com/cwb14/Kmer2LTR.git ~/.cache/printe/Kmer2LTR
```

or skip that phase entirely with `--no_postproc`. A burn-in never needs the network.

Set `PRINTE_CACHE` to move that directory somewhere else, which is what the container
does so it can write into a read-only image.

### Where are the bundled TE libraries?

Ask PrinTE:

```bash
python -m printe.paths                                # source install
apptainer exec printe.sif python -m printe.paths      # container
```

That lists every bundled file and where it resolves on your system. Anything shown as
`not found` is just not present in the way you installed: the reference libraries come
with a clone and with the container, but a bare `pip install` of the package gets the
code and the small `.tsv` tables only.

If you have the libraries elsewhere, point PrinTE at that directory:

```bash
export PRINTE_DATA=/path/to/your/data
```

`ltr-db.fa.gz` is different again - it is 38 MB and nothing reads it automatically, so it
ships as a release asset rather than in the repository. `make fetch-data` downloads it
into `~/.cache/printe/data`, where PrinTE looks by default.

### `mamba env create -f env.yml` fails

The file is `environment.yml` now. `env.yml` was removed in 1.0.0.

### Cloning is slow

The repository carries a long history. For a working copy you only need the tip:

```bash
git clone --depth 1 https://github.com/cwb14/PrinTE.git
```

### Do two runs with the same --seed give identical output?

Yes, on the same machine with the same thread count. `--seed` fixes both the parameter
draws and the per-worker streams. Changing `-t` changes how work is split across
processes and so changes the output, which is expected: rerun a comparison with the
same `-t` you used the first time.

### The grid search cannot find sbatch

`guided_search.py` and `run_loop.py` submit jobs from wherever you launch them, so they
need the Slurm client commands on PATH, i.e. a login node rather than inside a compute
allocation. `run_loop.py` checks for this and says so; the others will fail at submit
time. `--scheduler local` runs everything on the current machine instead.

### A TE family never gets inserted

Check `ratios.tsv` for that family. A row whose weight fails to parse is treated as
weight 0, which silently disables the family; since 1.0.0 that also prints a warning to
stderr, so search `pipeline.error` for `unparseable weight`. A family present in your
library but absent from `ratios.tsv` gets no weight either.
