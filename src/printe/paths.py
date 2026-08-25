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


if __name__ == "__main__":
    import sys

    print(data_file(sys.argv[1]) if len(sys.argv) > 1 else _PKG)
