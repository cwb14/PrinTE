import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"
FIXTURES = ("tiny_TE.lib", "tiny_cds.fa", "tiny_ratios.tsv")


@pytest.fixture(scope="session")
def repo():
    return REPO


@pytest.fixture
def tiny_lib():
    return DATA / "tiny_TE.lib"


@pytest.fixture
def tiny_cds():
    return DATA / "tiny_cds.fa"


@pytest.fixture
def tiny_ratios():
    return DATA / "tiny_ratios.tsv"


def run_printe(cwd, *args, timeout=1800):
    """Invoke PrinTE.sh in cwd. On failure, surface pipeline.error - that is where
    the actual traceback lands, and a bare exit code says nothing useful."""
    proc = subprocess.run(
        ["bash", str(REPO / "PrinTE.sh"), *(str(a) for a in args)],
        cwd=str(cwd), capture_output=True, text=True, timeout=timeout,
    )
    if proc.returncode != 0:
        err = cwd / "pipeline.error"
        raise AssertionError(
            f"PrinTE.sh exited {proc.returncode}\n"
            f"--- stderr ---\n{proc.stderr[-2000:]}\n"
            f"--- pipeline.error ---\n"
            f"{err.read_text()[-2000:] if err.is_file() else '(not written)'}"
        )
    return proc


def _seed_dir(d):
    for name in FIXTURES:
        shutil.copy(DATA / name, d / name)
    return d


@pytest.fixture
def rundir(tmp_path):
    return _seed_dir(tmp_path)


@pytest.fixture(scope="session")
def burnin_run(tmp_path_factory):
    """One burn-in shared across the read-only tests that inspect its output."""
    d = _seed_dir(tmp_path_factory.mktemp("burnin"))
    run_printe(
        d, "--burnin_only", "-sz", "60kb", "-cn", "1", "-P", "4", "-itp", "10",
        "-cl", "tiny_TE.lib", "-c", "tiny_cds.fa", "-r", "tiny_ratios.tsv",
        "-t", "2", "-s", "42",
    )
    return d
