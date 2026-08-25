"""log_to_report turns pipeline.log into the per-generation pipeline.report table.
The phrasings below are what PrinTE.sh's steps actually print."""

import subprocess
import sys


def build(tmp_path, log_text):
    log = tmp_path / "pipeline.log"
    log.write_text(log_text)
    out = tmp_path / "pipeline.report"
    proc = subprocess.run(
        [sys.executable, "-m", "printe.util.log_to_report", "-in", str(log), "-out", str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    return [ln for ln in out.read_text().splitlines() if ln.strip()]


LOG = (
    "Total TE insertions performed: 4 (Nested: 1, Non-nested: 3)\n"
    "Calculated number of TE excisions: 2\n"
    "Selected 2 removal events\n"
    "Updated FASTA written to gen100_final.fasta\n"
    "Total TE insertions performed: 6 (Nested: 2, Non-nested: 4)\n"
    "Calculated number of TE excisions: 1\n"
    "Selected 1 removal events\n"
    "Updated FASTA written to gen200_final.fasta\n"
)


def test_one_row_per_sampled_generation(tmp_path):
    rows = build(tmp_path, LOG)
    assert len(rows) >= 3  # header plus two generations


def test_generations_appear_in_order(tmp_path):
    body = "\n".join(build(tmp_path, LOG))
    assert "100" in body and "200" in body
    assert body.index("100") < body.index("200")


def test_insertion_counts_are_carried_through(tmp_path):
    body = "\n".join(build(tmp_path, LOG))
    assert "4" in body and "6" in body


def test_empty_log_does_not_crash(tmp_path):
    build(tmp_path, "")
