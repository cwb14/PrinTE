"""Locate the reference libraries and config tables that ship with PrinTE.

Data files are looked for in three places, in order: $PRINTE_DATA, the repo's data/
directory next to the package, then the user cache that 'make fetch-data' writes into.
A pip install carries the code but not the multi-MB libraries, so the cache is what
makes an installed PrinTE usable without a clone.
"""

import os
from pathlib import Path

_PKG = Path(__file__).resolve().parent


def _cache_root():
    xdg = os.environ.get("XDG_CACHE_HOME")
    return Path(xdg) if xdg else Path.home() / ".cache"


def _candidates(name, subdir):
    env = os.environ.get("PRINTE_DATA")
    if env:
        yield Path(env) / name
    # src/printe/ -> src/ -> repo root
    yield _PKG.parent.parent / subdir / name
    yield _cache_root() / "printe" / subdir / name


def _resolve(name, subdir):
    tried = []
    for c in _candidates(name, subdir):
        if c.is_file():
            return c
        tried.append(str(c))
    raise FileNotFoundError(
        "could not find {!r}. Looked in:\n  {}\n"
        "Run 'make fetch-data', or set PRINTE_DATA to the directory holding it.".format(
            name, "\n  ".join(tried)
        )
    )


def data_file(name):
    """Path to a bundled reference library, e.g. maize_rice_arab_curated_TE.lib.gz."""
    return _resolve(name, "data")


def config_file(name):
    """Path to a bundled config table, e.g. ratios.tsv.

    These are small enough to ship inside the wheel, so the packaged copy wins.
    """
    local = _PKG / "data" / name
    if local.is_file():
        return local
    return _resolve(name, "data")


# Everything that ships with PrinTE, so `python -m printe.paths` can list it.
CONFIG_FILES = ("ratios.tsv", "ratios_ltr_only.tsv")
DATA_FILES = (
    "maize_rice_arab_curated_TE.lib.gz",
    "TAIR10.cds.fa.gz",
    "TAIR10.pep.fa.gz",
    "athrep.updated.fasta_062024.ltr.full.gz",
    "rice7.0.0.liban.ltr.full.gz",
    "maizeTE02052020.ltr.full.gz",
    "ltr-db.fa.gz",
)


def find(name):
    """Resolve a bundled file whether it is a config table or a reference library."""
    return config_file(name) if name in CONFIG_FILES else data_file(name)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        print(find(sys.argv[1]))
    else:
        for name in CONFIG_FILES + DATA_FILES:
            try:
                print(f"{name:38s} {find(name)}")
            except FileNotFoundError:
                # fetch-data only pulls ltr-db; the rest ship with a clone or the image
                hint = "make fetch-data" if name == "ltr-db.fa.gz" else "use a clone or the container"
                print(f"{name:38s} not found  ({hint})")
