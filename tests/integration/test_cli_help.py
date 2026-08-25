"""Every shipped CLI must answer --help with exit 0.

Cheap, and it catches the whole class of breakage where a module is renamed, an
import goes stale, or a parser is misconfigured.
"""

import pkgutil
import subprocess
import sys

import pytest

import printe

# Library modules with no command line of their own.
NO_CLI = {
    "printe._cli",
    "printe.cutpaste_common",
    "printe.grid.grid_utils",
    "printe.grid.slurm_introspect",
    "printe.grid.slurm_sbatch",
    "printe.util._figure",
    "printe.paths",
}


def _modules():
    for m in pkgutil.walk_packages(printe.__path__, "printe."):
        if not m.ispkg and m.name not in NO_CLI:
            yield m.name


@pytest.mark.parametrize("mod", sorted(_modules()))
def test_help_exits_zero(mod):
    p = subprocess.run(
        [sys.executable, "-m", mod, "--help"], capture_output=True, text=True, timeout=600
    )
    assert p.returncode == 0, f"{mod} --help exited {p.returncode}\n{p.stderr[-600:]}"
    assert p.stdout.strip(), f"{mod} --help printed nothing"


def test_printe_sh_help(repo):
    p = subprocess.run(["bash", str(repo / "PrinTE.sh"), "--help"],
                       capture_output=True, text=True)
    assert p.returncode == 0
    assert "--generation_end" in p.stdout


def test_printe_sh_with_no_arguments_prints_usage(repo):
    p = subprocess.run(["bash", str(repo / "PrinTE.sh")], capture_output=True, text=True)
    assert p.returncode == 0
    assert "Common workflows" in p.stdout
