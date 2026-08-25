import shutil
from pathlib import Path

import pytest

from tests.conftest import FIXTURES, run_printe

pytestmark = pytest.mark.slow

DATA = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture(scope="module")
def evolved(tmp_path_factory):
    d = tmp_path_factory.mktemp("evolve")
    for name in FIXTURES:
        shutil.copy(DATA / name, d / name)
    run_printe(
        d, "-sz", "80kb", "-cn", "2", "-P", "5", "-itp", "12", "-ftp", "5",
        "-cl", "tiny_TE.lib", "-c", "tiny_cds.fa", "-r", "tiny_ratios.tsv",
        "-ir", "1e-5", "-dr", "1e-5", "-ge", "200", "-st", "100",
        "-t", "2", "-s", "7", "--no_postproc",
    )
    return d


def chrom_lengths(fasta):
    lengths, name = {}, None
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        else:
            lengths[name] += len(line.strip())
    return lengths


def test_one_genome_and_annotation_per_sampled_generation(evolved):
    for gen in (100, 200):
        for ext in ("fasta", "bed"):
            f = evolved / f"gen{gen}_final.{ext}"
            assert f.is_file(), f"missing gen{gen}_final.{ext}"
            assert f.stat().st_size > 0


def test_only_the_final_generation_keeps_its_library(evolved):
    # Each step's .lib feeds the next one and is then removed, unless --keep_temps.
    assert (evolved / "gen200_final.lib").is_file()
    assert not (evolved / "gen100_final.lib").exists()


def test_unsampled_generations_are_not_written(evolved):
    assert not (evolved / "gen150_final.fasta").exists()
    assert not (evolved / "gen300_final.fasta").exists()


def test_mutation_report_written_per_step(evolved):
    for gen in (100, 200):
        assert (evolved / f"gen{gen}_mut.txt").is_file()


def test_pipeline_report_is_written(evolved):
    assert (evolved / "pipeline.report").is_file()


def test_library_headers_are_bare_feature_ids(evolved):
    heads = [
        ln for ln in (evolved / "gen200_final.lib").read_text().splitlines()
        if ln.startswith(">")
    ]
    assert heads
    for h in heads:
        assert "#" in h, f"library header without a class: {h}"
        assert " " not in h.strip(), f"whitespace in header: {h}"


def test_bed_stays_inside_the_genome(evolved):
    lengths = chrom_lengths(evolved / "gen200_final.fasta")
    for line in (evolved / "gen200_final.bed").read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        assert f[0] in lengths, line
        assert 0 <= int(f[1]) < int(f[2]) <= lengths[f[0]], line


def test_pipeline_error_holds_only_the_start_banner_on_success(evolved):
    # The README tells users that anything beyond the banner means a real failure.
    text = (evolved / "pipeline.error").read_text()
    assert "Traceback" not in text
    assert "Error running" not in text
