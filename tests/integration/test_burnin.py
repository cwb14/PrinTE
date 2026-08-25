import pytest

pytestmark = pytest.mark.slow


def test_writes_the_documented_outputs(burnin_run):
    for name in ("burnin.fasta", "burnin.bed", "burnin.stat"):
        f = burnin_run / name
        assert f.is_file(), f"missing {name}"
        assert f.stat().st_size > 0, f"{name} is empty"


def test_genome_is_acgt_only(burnin_run):
    seq = "".join(
        ln.strip()
        for ln in (burnin_run / "burnin.fasta").read_text().splitlines()
        if not ln.startswith(">")
    )
    assert set(seq) <= set("ACGT"), sorted(set(seq) - set("ACGT"))


def test_genome_is_about_the_requested_size(burnin_run):
    # -sz 60kb was requested; TE insertion grows it, so allow generous headroom.
    seq = "".join(
        ln.strip()
        for ln in (burnin_run / "burnin.fasta").read_text().splitlines()
        if not ln.startswith(">")
    )
    assert 40_000 < len(seq) < 200_000, len(seq)


def test_stat_reports_a_plausible_te_fraction(burnin_run):
    stat = (burnin_run / "burnin.stat").read_text()
    pct = float(stat.split("TE bp inserted=")[1].split("(")[1].split("%")[0])
    assert 0.0 < pct < 40.0, stat


def test_chromosome_count_matches_the_request(burnin_run):
    heads = [
        ln for ln in (burnin_run / "burnin.fasta").read_text().splitlines()
        if ln.startswith(">")
    ]
    assert len(heads) == 1, heads
