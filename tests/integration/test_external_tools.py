"""Where PrinTE puts the tools it builds or fetches: its own directory, or $PRINTE_CACHE."""

import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from printe import ltr_dens
from tests.conftest import LOCATION_VARS, _seed_dir, run_printe

EVOLVE = (
    "-sz", "80kb", "-cn", "2", "-P", "5", "-itp", "25", "-ftp", "5",
    "-cl", "tiny_TE.lib", "-c", "tiny_cds.fa", "-r", "tiny_ratios.tsv",
    "-ir", "1e-5", "-dr", "1e-5", "-ge", "200", "-st", "100", "-t", "2", "-s", "7",
)
BURNIN = (
    "--burnin_only", "-sz", "60kb", "-cn", "1", "-P", "4", "-itp", "10",
    "-cl", "tiny_TE.lib", "-c", "tiny_cds.fa", "-r", "tiny_ratios.tsv", "-t", "2", "-s", "42",
)
# For runs that stop at the Kmer2LTR check, before any mutation: anything that passes
# PrinTE's probe of the mutator will do.
TRUE = shutil.which("true")

FAKE_GIT = """#!/bin/sh
# Stands in for git: logs each call, and "clones" a minimal legacy Kmer2LTR the way real
# git does - the destination appears at once and is filled in after a pause.
echo "$*" >> "{log}"
if [ "$1" = clone ]; then
  for dest; do :; done
  mkdir -p "$dest"
  sleep 3
  : > "$dest/Kmer2LTR.py"
fi
"""


def clean_env(strip_mutator=False, **extra):
    """This process's environment minus anything that hands PrinTE a mutator or a place
    to put things, so each test decides those itself."""
    env = {k: v for k, v in os.environ.items() if k not in (*LOCATION_VARS, "PRINTE_MUTATOR")}
    if strip_mutator:
        keep = [p for p in env["PATH"].split(os.pathsep) if not (Path(p) / "ltr_mutator").exists()]
        python = shutil.which("python", path=env["PATH"])
        if python and str(Path(python).parent) not in keep:
            pytest.skip("ltr_mutator sits next to python on PATH, so it cannot be hidden")
        env["PATH"] = os.pathsep.join(keep)
    env.update(extra)
    return env


@pytest.fixture
def built_mutator(repo):
    exe = repo / "bin" / "ltr_mutator"
    if not os.access(exe, os.X_OK):
        pytest.skip("needs a built ltr_mutator (make ltr-mutator)")
    return str(exe)


@pytest.fixture
def fake_git(tmp_path):
    """(directory to put first on PATH, log file that exists once git has been called)"""
    bindir, log = tmp_path / "fakebin", tmp_path / "git_calls.log"
    bindir.mkdir()
    (bindir / "git").write_text(FAKE_GIT.format(log=log))
    (bindir / "git").chmod(0o755)
    return bindir, log


def with_fake_git(env, fake):
    env["PATH"] = f"{fake[0]}{os.pathsep}{env['PATH']}"
    return env


def github_reachable():
    try:
        proc = subprocess.run(
            ["git", "ls-remote", "--exit-code", "https://github.com/cwb14/Kmer2LTR.git", "HEAD"],
            capture_output=True, timeout=30, env={**os.environ, "GIT_TERMINAL_PROMPT": "0"},
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0


@pytest.mark.slow
def test_mutator_is_built_under_printe_cache(rundir, tmp_path):
    cache = tmp_path / "cache"
    run_printe(rundir, *EVOLVE, "--no_postproc",
               env=clean_env(strip_mutator=True, PRINTE_CACHE=str(cache)))
    assert (cache / "bin" / "ltr_mutator").is_file()
    assert (rundir / "gen200_final.fasta").is_file()


@pytest.mark.slow
def test_mutator_that_no_longer_runs_is_rebuilt(rundir, tmp_path):
    # A checkout shared by machines with different C libraries: the binary on disk is
    # newer than its source, so make alone would call it up to date.
    stale = tmp_path / "cache" / "bin" / "ltr_mutator"
    stale.parent.mkdir(parents=True)
    stale.write_text("#!/bin/sh\necho \"ltr_mutator: version 'GLIBC_2.99' not found\" >&2\nexit 1\n")
    stale.chmod(0o755)
    run_printe(rundir, *EVOLVE, "--no_postproc",
               env=clean_env(strip_mutator=True, PRINTE_CACHE=str(tmp_path / "cache")))
    assert (rundir / "gen200_final.fasta").is_file()


@pytest.mark.skipif(os.geteuid() == 0, reason="root writes through directory permissions")
def test_kmer2ltr_that_cannot_be_fetched_fails_before_simulating(rundir, tmp_path):
    cache = tmp_path / "read_only"
    cache.mkdir(mode=0o500)
    try:
        with pytest.raises(AssertionError, match="Kmer2LTR"):
            run_printe(rundir, *EVOLVE, "--postproc",
                       env=clean_env(PRINTE_CACHE=str(cache), PRINTE_MUTATOR=TRUE))
    finally:
        cache.chmod(0o700)
    assert not (rundir / "burnin.fasta").exists()


def test_incompatible_kmer2ltr_fails_before_simulating(rundir, tmp_path):
    # What a clone of Kmer2LTR's main branch looks like to PrinTE: no Kmer2LTR.py.
    (tmp_path / "cache" / "Kmer2LTR" / "src").mkdir(parents=True)
    with pytest.raises(AssertionError, match="Kmer2LTR"):
        run_printe(rundir, *EVOLVE, "--postproc",
                   env=clean_env(PRINTE_CACHE=str(tmp_path / "cache"), PRINTE_MUTATOR=TRUE))
    assert not (rundir / "burnin.fasta").exists()


def test_argument_errors_are_reported_before_fetching_kmer2ltr(rundir, tmp_path, fake_git):
    env = with_fake_git(clean_env(PRINTE_CACHE=str(tmp_path / "cache"), PRINTE_MUTATOR=TRUE),
                        fake_git)
    with pytest.raises(AssertionError):
        run_printe(rundir, *EVOLVE, "--postproc", "-b", "genome.bed", env=env)  # --bed without --fasta
    assert not fake_git[1].exists(), "git ran before the arguments were checked"


@pytest.mark.slow
@pytest.mark.parametrize("args", [BURNIN, EVOLVE, (*EVOLVE, "--no_postproc")],
                         ids=["burnin_only", "default", "no_postproc"])
def test_runs_without_postprocessing_never_fetch_kmer2ltr(rundir, tmp_path, fake_git,
                                                         built_mutator, args):
    env = with_fake_git(clean_env(PRINTE_CACHE=str(tmp_path / "cache"),
                                  PRINTE_MUTATOR=built_mutator), fake_git)
    run_printe(rundir, *args, env=env)
    assert not fake_git[1].exists(), fake_git[1].read_text()


@pytest.mark.slow
def test_runs_started_together_share_one_kmer2ltr_clone(tmp_path, fake_git, built_mutator):
    cache = tmp_path / "cache"
    env = with_fake_git(clean_env(PRINTE_CACHE=str(cache), PRINTE_MUTATOR=built_mutator),
                        fake_git)
    runs = []
    for name in ("a", "b"):
        (tmp_path / name).mkdir()
        runs.append(_seed_dir(tmp_path / name))
    with ThreadPoolExecutor(2) as pool:
        first = pool.submit(run_printe, runs[0], *EVOLVE, "--postproc", env=env)
        time.sleep(1)  # the second run reaches the check while the first is still cloning
        second = pool.submit(run_printe, runs[1], *EVOLVE, "--postproc", env=env)
        first.result()
        second.result()
    assert (cache / "Kmer2LTR" / "Kmer2LTR.py").is_file()
    assert [p.name for p in cache.iterdir()] == ["Kmer2LTR"], "a partial clone was left behind"


@pytest.mark.slow
def test_postprocessing_dates_the_ltr_rts(rundir, tmp_path, built_mutator):
    if not github_reachable():
        pytest.skip("clones Kmer2LTR from GitHub, which is not reachable")
    cache = tmp_path / "cache"
    run_printe(rundir, *EVOLVE, "--postproc",
               env=clean_env(PRINTE_CACHE=str(cache), PRINTE_MUTATOR=built_mutator))
    assert (cache / "Kmer2LTR" / "Kmer2LTR.py").is_file()
    for label in ("burnin", "gen200_final"):
        # What the density plot reads: K2P distances, not columns shifted out of place.
        k2p = ltr_dens.read_ltr_tsv(str(rundir / f"{label}_LTR.tsv"))["K2P_d"]
        assert len(k2p), f"{label}_LTR.tsv is empty"
        assert k2p.between(0, 1).all(), f"{label}: K2P_d {k2p.tolist()}"
    assert (rundir / "all_LTR_density.pdf").is_file()
