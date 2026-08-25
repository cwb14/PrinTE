import pytest

from printe.paths import config_file


def rows(name):
    return [
        ln.split("\t")
        for ln in config_file(name).read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]


@pytest.mark.parametrize("name", ["ratios.tsv", "ratios_ltr_only.tsv"])
def test_shipped_ratios_parse(name):
    got = rows(name)
    assert got
    for r in got:
        assert len(r) >= 4, r
        float(r[2])
        float(r[3])


@pytest.mark.parametrize("name", ["ratios.tsv", "ratios_ltr_only.tsv"])
def test_weights_are_non_negative(name):
    for r in rows(name):
        assert float(r[2]) >= 0 and float(r[3]) >= 0, r


def test_transposition_column_only_uses_known_values():
    for r in rows("ratios.tsv"):
        if len(r) > 4 and r[4].strip():
            assert r[4].strip().lower() in ("cutpaste", "copypaste"), r


def test_default_ratios_weights_sum_to_about_one():
    total = sum(float(r[2]) for r in rows("ratios.tsv"))
    assert total == pytest.approx(1.0, abs=0.02)


def test_ltr_only_file_puts_all_its_weight_on_ltr():
    got = rows("ratios_ltr_only.tsv")
    total = sum(float(r[2]) for r in got)
    ltr = sum(float(r[2]) for r in got if r[0] == "LTR")
    assert total == pytest.approx(1.0, abs=0.02)
    assert ltr == pytest.approx(total, abs=1e-9)


def test_fixture_ratios_cover_every_family_in_the_fixture_library(tiny_lib, tiny_ratios):
    import re

    fams = set()
    for line in tiny_lib.read_text().splitlines():
        if line.startswith(">"):
            m = re.match(r">[^#]+#([^/]+)/([^~\s]+)", line)
            assert m, line
            fams.add((m.group(1), m.group(2)))
    declared = {
        (r[0], r[1])
        for r in [
            ln.split("\t")
            for ln in tiny_ratios.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        ]
    }
    assert fams <= declared, fams - declared
