import subprocess
import sys


def run(tmp_path, **over):
    args = {
        "--ins-start": "1e-9", "--ins-end": "1e-11", "--ins-count": "3",
        "--del-start": "1e-9", "--del-end": "1e-11", "--del-count": "2",
        "--sr-start": "90", "--sr-end": "95", "--sr-step": "5",
        "--k-start": "0", "--k-end": "1", "--k-step": "1",
    }
    args.update(over)
    out = tmp_path / "combos.csv"
    flat = [x for kv in args.items() for x in kv]
    p = subprocess.run(
        [sys.executable, "-m", "printe.grid.emit_combos", *flat, "--out", str(out)],
        capture_output=True, text=True,
    )
    assert p.returncode == 0, p.stderr
    return [ln for ln in out.read_text().splitlines() if ln.strip()]


def test_header_matches_what_the_workflow_reads(tmp_path):
    assert run(tmp_path)[0] == "ins,del,sr,k,dirname"


def test_emits_the_full_cartesian_product(tmp_path):
    assert len(run(tmp_path)) - 1 == 3 * 2 * 2 * 2


def test_counts_are_honoured(tmp_path):
    rows = run(tmp_path, **{"--ins-count": "5", "--del-count": "4"})
    assert len(rows) - 1 == 5 * 4 * 2 * 2


def test_dirname_matches_the_grid_search_convention(tmp_path):
    from printe.grid.grid_utils import dir_name_from_combo

    for row in run(tmp_path)[1:]:
        ins, dele, sr, k, dirname = row.split(",")
        assert dirname == dir_name_from_combo(ins, dele, sr, k)


def test_rates_span_the_requested_range(tmp_path):
    rates = {r.split(",")[0] for r in run(tmp_path)[1:]}
    vals = sorted(float(x) for x in rates)
    assert vals[0] == 1e-11
    assert vals[-1] == 1e-9
