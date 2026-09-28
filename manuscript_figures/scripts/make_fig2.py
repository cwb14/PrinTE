#!/usr/bin/env python3
"""Build final Fig. 2 (benchmarking PrinTE) for Nature Communications.

Wraps the approved plotting script merged_plot28.py (imported read-only) and
  * switches all text to Arial (TrueType, editable) with ~7 pt text at final size,
  * uses mean +/- SD (n = 3 replicate genomes) for panel D (was 1.96 x s.e.m.),
  * scales the vector PDF to 180 mm width,
  * writes the Source Data blocks (build/source_data/fig2).
"""
import argparse
import glob
import locale
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ncomms_style  # noqa: E402
from srcdata import SourceData  # noqa: E402

V10 = "/home/chris/data/TEGE/benchmarking/v2/v3/v4/v5/v6/v7/v8/v9/v10"
SRC = f"{V10}/v11/v12"
MM = 72 / 25.4
# pre-scale sizes; x ~0.745 at 180 mm -> ~7 pt text, 8 pt panel letters
RC = {"font.size": 9.5, "axes.labelsize": 10, "axes.titlesize": 10.5,
      "legend.fontsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9}
LETTER_PT = 10.75


def shell_glob(pattern):
    """Same order as bash's glob under en_US.UTF-8 (the order used for the published panel B)."""
    files = glob.glob(pattern)
    try:
        locale.setlocale(locale.LC_COLLATE, "en_US.UTF-8")
        return sorted(files, key=locale.strxfrm)
    except locale.Error:
        return sorted(files)


def per_replicate_A():
    """Per-replicate panel A values (PrinTE Kmer2LTR summaries; SLiM TiTv=1 minimap2 stats)."""
    rows = []
    for ds in ["TiTv1", "TiTv10"]:
        for rep in ["rep1", "rep2", "rep3"]:
            d = f"{V10}/{ds}/{rep}"
            for fn in os.listdir(d):
                m = re.match(r"gen(\d+)_Kmer2LTR\.summary$", fn)
                if not m:
                    continue
                g = int(m.group(1))
                vals = dict(ln.strip().split("\t") for ln in open(f"{d}/{fn}") if "\t" in ln)
                nltr = sum(1 for ln in open(f"{d}/gen{g}_Kmer2LTR") if ln.strip() and not ln.startswith("#"))
                rows.append(dict(dataset=f"PrinTE {ds}", replicate=rep, generation=g, n_LTR_RTs=nltr,
                                 p_dist=float(vals["raw_d"]), JC69=float(vals["JC69_d"]), K2P=float(vals["K2P_d"])))
    for rep in ["rep1", "rep2", "rep3"]:
        d = f"{V10}/v11/TiTv1/{rep}"
        for fn in os.listdir(d):
            m = re.match(r"gen(\d+)_1_vs_gen\d+_2\.stat$", fn)
            if not m:
                continue
            t = open(f"{d}/{fn}").read()
            mp = int(re.search(r"mapped bases: (\d+)", t).group(1))
            sb = int(re.search(r"substitutions: (\d+)", t).group(1))
            rows.append(dict(dataset="SLiM/minimap2 TiTv1", replicate=rep, generation=int(m.group(1)),
                             mapped_bp=mp, substitutions=sb, p_dist=(sb / mp if mp else 0.0)))
    return pd.DataFrame(rows).sort_values(["dataset", "generation", "replicate"]).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src-dir", default=SRC, help="dir with merged_plot28.py and its inputs")
    ap.add_argument("--metrics-error", choices=["sd", "se", "ci95"], default="sd",
                    help="panel D error bars (default sd; published draft used ci95)")
    ap.add_argument("--width-mm", type=float, default=180.0, help="final width in mm")
    ap.add_argument("-o", "--out", default="fig2/fig2.pdf", help="output PDF")
    ap.add_argument("-v", "--verbose", action="store_true", help="verbose")
    a = ap.parse_args()

    import fitz
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    src = Path(a.src_dir)
    sys.dont_write_bytecode = True  # never write __pycache__ into the read-only source dir
    sys.path.insert(0, str(src))
    import merged_plot28 as mp  # sets its own rcParams on import; we override below

    ncomms_style.arial()
    plt.rcParams.update(RC)
    mp.panel_label.__defaults__ = (3, LETTER_PT)

    dens = shell_glob(str(src / "*.LTRs.alns.results"))
    efg = [str(src / f) for f in ("s8.tsv", "s9_subset2.tsv", "s10.tsv")]
    native = HERE.parent / "build" / "fig2" / "fig2_native.pdf"
    native.parent.mkdir(parents=True, exist_ok=True)
    sys.argv = ["merged_plot28.py", "--div-csv", str(src / "divergence_aggregated_k2p_TiTv10_new.csv"),
                "--metrics-csv", str(src / "repeatmasker.per_replicate_metrics.csv"),
                "--threads-csv", str(src / "threads_perf.csv"), "--metrics-error", a.metrics_error,
                "--panelA-alpha", "0.6", "--efg-tsv", *efg, "--density-in", *dens, "--out", str(native)]
    mp.main()

    # scale the vector page to the final width
    doc = fitz.open(native)
    w = a.width_mm * MM
    h = doc[0].rect.height * w / doc[0].rect.width
    out = fitz.open()
    out.new_page(width=w, height=h).show_pdf_page(fitz.Rect(0, 0, w, h), doc, 0)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.save(a.out, garbage=4, deflate=True)
    print(f"wrote {a.out}: {a.width_mm:.0f} x {h / MM:.1f} mm (scale {w / doc[0].rect.width:.3f})", file=sys.stderr)

    # ------------------------------------------------------------ Source Data
    sd = SourceData("fig2")
    div = pd.read_csv(src / "divergence_aggregated_k2p_TiTv10_new.csv")
    div.insert(2, "n_replicates", np.where(div["generation"] == 0, 0, 3))
    sd.add("Fig. 2A (summary)", "Mean and SD (sample SD, n = 3 replicate simulated genomes per point; generation 0 is "
           "the fixed origin, not data) of LTR-RT divergence by substitution model. SLiM/minimap2 rows: 0 = no "
           "whole-genome alignment obtained (>= 4.5 Mgen; 8.0-8.5 Mgen not run). expected_divergence_mean = red "
           "dashed line.", div)
    per = per_replicate_A()
    sd.add("Fig. 2A (replicates)", "Per-replicate values behind each point: PrinTE divergence pooled over all intact "
           "LTR-RTs of one replicate genome (n_LTR_RTs); SLiM/minimap2 TiTv=1 p-distance = substitutions / mapped "
           "bases of the whole-genome alignment (0 = unaligned).", per)
    rowsB = []
    for f in dens:
        name = os.path.basename(f).split(".", 1)[0]
        k2p = pd.read_csv(f, sep="\t", header=None).iloc[:, 10].astype(float)
        summ = dict(ln.split() for ln in open(f + ".summary") if ln.strip())
        rowsB += [{"genome": name, "LTR_RT": i + 1, "K2P_divergence": v, "genome_pooled_K2P": float(summ["K2P_d"])}
                  for i, v in enumerate(k2p)]
    sd.add("Fig. 2B", "K2P divergence of each intact LTR-RT (Kmer2LTR) per simulated genome; curves are Gaussian "
           "KDEs (Silverman bandwidth). genome_pooled_K2P = grey dashed line. Axis shows 0-0.20 "
           "(athrep1-9 visible).", pd.DataFrame(rowsB))
    tsd = pd.DataFrame({
        "simulator": ["SimulaTE"] * 3 + ["ReplicaTE"] * 3 + ["PrinTE no TSD"] * 3 + ["PrinTE"] * 3,
        "replicate": [1, 2, 3] * 4,
        "LTR_RTs_inserted": [682, 675, 616, 25, 24, 28, 1000, 1000, 1000, 1000, 1000, 1000],
        "detection_rate": np.array([7.77, 8.00, 8.28, 12.00, 20.83, 10.71, 8.10, 7.70, 6.10, 78.20, 74.00, 78.80]) / 100})
    tsd["LTR_RTs_detected"] = (tsd["detection_rate"] * tsd["LTR_RTs_inserted"]).round().astype(int)
    sd.add("Fig. 2C", "LTR_retriever detection rate = detected / inserted intact LTR-RTs per replicate simulated "
           "genome (n = 3 per simulator). Box plots: centre line = median; box = 25th-75th percentiles; whiskers "
           "= min-max (1.5 x IQR rule); dots = replicates. LTR_RTs_detected = rate x inserted.", tsd)
    met = pd.read_csv(src / "repeatmasker.per_replicate_metrics.csv")
    sd.add("Fig. 2D", "RepeatMasker per-base performance per replicate simulated genome (500 Mb; n = 3 per "
           "condition; TP+FP+FN+TN = bp). Bars = mean; error bars = SD; dots = replicates.", met)
    thr = pd.read_csv(src / "threads_perf.csv")
    sd.add("Fig. 2E", "PrinTE wall-clock runtime and peak memory (GNU time) vs CPU threads; one run per thread "
           "count (n = 1, no error bars).", thr)
    for panel, f, var in [("Fig. 2F", efg[1], "deletion rate"), ("Fig. 2G", efg[0], "deletion length bias k"),
                          ("Fig. 2H", efg[2], "solo LTR formation rate")]:
        p = mp.load_param_tsv(Path(f)).rename(columns={"metadata": "simulation_parameters",
                                                         "gen_million": "generation_millions",
                                                         "size_mb": "genome_size_Mb"})
        sd.add(panel, f"Genome size across generations under varying {var}; one simulation per setting (n = 1). "
               "Size parsed from the file-size listing of each generation's FASTA.", p)
    sd.save()
    print(f"source data blocks: {len(sd.blocks)}", file=sys.stderr)


if __name__ == "__main__":
    main()
