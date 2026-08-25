import pytest

from printe import paths


def test_config_file_finds_bundled_ratios():
    p = paths.config_file("ratios.tsv")
    assert p.is_file()
    assert p.read_text().lstrip().startswith("#")


def test_ltr_only_ratios_also_bundled():
    assert paths.config_file("ratios_ltr_only.tsv").is_file()


def test_env_override_wins(tmp_path, monkeypatch):
    fake = tmp_path / "maize_rice_arab_curated_TE.lib.gz"
    fake.write_bytes(b"")
    monkeypatch.setenv("PRINTE_DATA", str(tmp_path))
    assert paths.data_file("maize_rice_arab_curated_TE.lib.gz") == fake


def test_missing_file_names_every_location_searched(tmp_path, monkeypatch):
    monkeypatch.setenv("PRINTE_DATA", str(tmp_path))
    with pytest.raises(FileNotFoundError) as e:
        paths.data_file("nope.fa.gz")
    msg = str(e.value)
    assert str(tmp_path) in msg
    assert "make fetch-data" in msg
