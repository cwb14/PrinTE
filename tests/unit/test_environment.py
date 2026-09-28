"""environment.yml must carry everything the R/ scripts load, so one env runs all of PrinTE."""

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

# Ship with every R installation, so they never appear in environment.yml.
BASE_R = {
    "base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods",
    "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils",
}


def r_packages():
    """Packages R/*.R loads, by library()/require() or a pkg:: reference."""
    pkgs = set()
    for script in (REPO / "R").glob("*.R"):
        for line in script.read_text().splitlines():
            code = line.split("#", 1)[0]
            pkgs.update(re.findall(r"\b(?:library|require)\(\s*([A-Za-z][\w.]*)\s*\)", code))
            pkgs.update(re.findall(r"\b([A-Za-z][\w.]*)::", code))
    return pkgs - BASE_R


def declared_packages():
    text = (REPO / "environment.yml").read_text()
    return set(re.findall(r"^\s*-\s*([A-Za-z0-9][\w.\-]*)", text, re.MULTILINE))


def test_environment_declares_every_r_package():
    declared = declared_packages()
    missing = sorted(
        p for p in r_packages()
        if f"r-{p.lower()}" not in declared and f"bioconductor-{p.lower()}" not in declared
    )
    assert not missing, f"R/ loads packages that environment.yml does not install: {missing}"


@pytest.mark.needs_r
@pytest.mark.skipif(shutil.which("Rscript") is None, reason="needs Rscript on PATH")
def test_every_r_package_loads():
    loads = "; ".join(f"library({p})" for p in sorted(r_packages()))
    proc = subprocess.run(
        ["Rscript", "-e", f"suppressPackageStartupMessages({{{loads}}})"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
