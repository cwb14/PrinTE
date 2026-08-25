"""Console entry points.

'printe' is a thin exec of PrinTE.sh rather than a reimplementation, so the installed
command and the 'bash PrinTE.sh' invocation in the paper cannot drift apart.
"""

import os
import site
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()


def _candidates():
    # Working from a checkout: src/printe/_cli.py, with the script at the repo root.
    yield _HERE.parents[2] / "PrinTE.sh"
    # Installed: pyproject puts it in <prefix>/share/printe. sys.prefix covers venvs
    # and conda envs, base_prefix the case where a venv inherits from the base env,
    # and USER_BASE a `pip install --user`.
    for base in (sys.prefix, sys.base_prefix, site.USER_BASE):
        if base:
            yield Path(base) / "share" / "printe" / "PrinTE.sh"
    # Container image layout.
    yield Path("/opt/printe/PrinTE.sh")


def _pipeline():
    env = os.environ.get("PRINTE_SCRIPT")
    if env:
        if Path(env).is_file():
            return Path(env)
        sys.exit(f"printe: PRINTE_SCRIPT is set to {env}, which is not a file")

    tried = []
    for cand in _candidates():
        if cand.is_file():
            return cand
        tried.append(str(cand))
    sys.exit(
        "printe: cannot find PrinTE.sh. Looked in:\n  "
        + "\n  ".join(tried)
        + "\nSet PRINTE_SCRIPT to its path."
    )


def main():
    # Answer --version here rather than letting PrinTE.sh shell out to `python`,
    # which in an installed layout need not be the interpreter holding the package.
    if len(sys.argv) == 2 and sys.argv[1] in ("-v", "--version"):
        from printe import __version__

        print(f"PrinTE {__version__}")
        return
    os.execv("/bin/bash", ["bash", str(_pipeline()), *sys.argv[1:]])
