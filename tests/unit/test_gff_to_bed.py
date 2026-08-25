import subprocess
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"


def convert(tmp_path):
    out = tmp_path / "out.bed"
    proc = subprocess.run(
        [sys.executable, "-m", "printe.util.gff_to_bed", str(DATA / "sample.gff"), str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    return [ln for ln in out.read_text().splitlines() if ln.strip()]


def test_output_is_six_column_printe_bed(tmp_path):
    for line in convert(tmp_path):
        f = line.split("\t")
        assert len(f) == 6, f
        assert int(f[1]) < int(f[2]), f
        assert f[5] in ("+", "-"), f


def test_every_gff_record_is_converted(tmp_path):
    assert len(convert(tmp_path)) == 5


def test_feature_ids_carry_class_and_superfamily(tmp_path):
    for line in convert(tmp_path):
        assert "#" in line.split("\t")[3]


def test_partial_elements_are_marked_frag(tmp_path):
    frag = [ln for ln in convert(tmp_path) if "_FRAG" in ln.split("\t")[3]]
    # Integrity 0.62 and 0.45 are the two non-intact records in the fixture
    assert len(frag) == 2


def test_tsd_is_a_sequence_or_na(tmp_path):
    for line in convert(tmp_path):
        tsd = line.split("\t")[4]
        assert tsd == "NA" or set(tsd.upper()) <= set("ACGTN"), tsd
