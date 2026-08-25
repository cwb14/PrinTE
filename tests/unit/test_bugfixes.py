"""Things that shipped broken, plus a few greps so they do not come back."""

import ast
import subprocess
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "src" / "printe"


def test_fasta_to_repeatmasker_parses():
    # An SO.txt example was pasted into the header without comment markers, so the
    # file raised SyntaxError on import while the README told people to run it.
    ast.parse((SRC / "util" / "fasta_to_RepeatMasker.py").read_text())


def test_compare_script_default_points_at_a_real_file():
    from printe.grid.build_composite_matrix import build_parser

    default = build_parser().get_default("compare_script")
    assert "compare_genomes2" not in default
    assert Path(default).is_file()


def test_headless_plotters_select_agg_before_importing_pyplot():
    for path in sorted(SRC.rglob("*.py")):
        src = path.read_text()
        if "import matplotlib.pyplot" not in src:
            continue
        assert 'matplotlib.use("Agg")' in src, f"{path.name}: no Agg backend"
        assert src.index('matplotlib.use("Agg")') < src.index("import matplotlib.pyplot"), (
            f"{path.name}: Agg selected after pyplot import"
        )


def test_no_absolute_home_paths_in_shipped_source():
    hits = subprocess.run(
        ["grep", "-rn", "-e", "/home/chris", "-e", "/users/bin", "-e", "/users/data", str(SRC)],
        capture_output=True, text=True,
    )
    assert hits.stdout == "", hits.stdout


def test_no_bare_except_clauses():
    for path in sorted(SRC.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ExceptHandler):
                assert node.type is not None, f"{path.name}:{node.lineno} bare except"


def test_modules_do_not_execute_on_import():
    # genome_plot used to run its whole body at import time and sys.exit().
    proc = subprocess.run(
        [sys.executable, "-c", "import printe.util.genome_plot"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr


def test_stitcher_split_into_two_working_programs():
    for name in ("TE_lib_stitcher.py", "merge_ltr_parts.py"):
        src = (SRC / "util" / name).read_text()
        assert src.count("#!/usr/bin/env python3") == 1, f"{name}: more than one program"
        assert src.count('if __name__ ==') == 1, f"{name}: more than one entry point"
