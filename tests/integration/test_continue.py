import pytest

from tests.conftest import run_printe

pytestmark = pytest.mark.slow

BASE = [
    "-sz", "80kb", "-cn", "1", "-P", "5", "-itp", "12",
    "-cl", "tiny_TE.lib", "-c", "tiny_cds.fa", "-r", "tiny_ratios.tsv",
    "-ir", "1e-5", "-dr", "1e-5", "-st", "100", "-t", "2", "-s", "7",
    "--no_postproc",
]


def test_continue_extends_a_finished_run(rundir):
    run_printe(rundir, *BASE, "-ge", "100")
    assert (rundir / "gen100_final.fasta").is_file()
    assert not (rundir / "gen200_final.fasta").exists()

    run_printe(rundir, *BASE, "-ge", "200", "--continue")
    assert (rundir / "gen200_final.fasta").is_file()
    # the first step's output must survive the resume
    assert (rundir / "gen100_final.fasta").is_file()
