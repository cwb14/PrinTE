"""project_terminal signals a size breach to PrinTE.sh through its exit code:
10 = confidently above --mxgs, 11 = confidently below --mngs, 0 = keep going."""

import subprocess
import sys


def run(*args):
    return subprocess.run(
        [sys.executable, "-m", "printe.project_terminal", *(str(a) for a in args)],
        capture_output=True, text=True,
    )


BOUNDS = ("--clamp-lo", 1, "--clamp-hi", 10_000_000_000)


def test_exit_10_when_growing_past_the_maximum():
    p = run("--iters", "0,1,2,3", "--sizes", "100,140,190,260", "--target", 10,
            "--mxgs", 250, *BOUNDS)
    assert p.returncode == 10, p.stdout + p.stderr


def test_exit_11_when_shrinking_past_the_minimum():
    p = run("--iters", "0,1,2,3", "--sizes", "260,190,140,100", "--target", 10,
            "--mngs", 110, *BOUNDS)
    assert p.returncode == 11, p.stdout + p.stderr


def test_exit_0_while_inside_the_band():
    p = run("--iters", "0,1,2,3", "--sizes", "200,201,199,200", "--target", 10,
            "--mxgs", 400, "--mngs", 100, *BOUNDS)
    assert p.returncode == 0, p.stdout + p.stderr


def test_exit_0_when_there_is_too_little_data_to_fit():
    p = run("--iters", "0", "--sizes", "100", "--target", 10, "--mxgs", 50, *BOUNDS)
    assert p.returncode == 0


def test_no_breach_when_the_last_step_moved_back_toward_the_band():
    # over the max, but shrinking again - not a confident breach
    p = run("--iters", "0,1,2,3", "--sizes", "100,400,380,300", "--target", 10,
            "--mxgs", 250, *BOUNDS)
    assert p.returncode == 0, p.stdout + p.stderr
