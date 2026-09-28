#!/usr/bin/env python3
"""Build the final Fig. 7 (Nature Communications) — one page, 180 mm wide.

Derived from /data2/chris/fungi/best/plots/panel_plot.py (page 1, run
2026-08-19 as `panel_plot.py --k2p-counts -o panel_plot_counts.pdf -v`).
Data, computed values, the starred best fits and the best-5% set definition
are unchanged; this script only re-lays the page out at print size and adds:

  * Arial everywhere (ncomms_style), TrueType-embedded; panel A is re-drawn by
    fig7_contmap.R through cairo_pdf(family="Arial") at 1:1 print size and
    embedded as a vector XObject (PyMuPDF show_pdf_page, no scaling).
  * data points on the bar charts M-P: the best-5% simulations behind each
    whisker (M, N, P) and the per-generation values of the best run (O).
  * Source Data blocks (srcdata.SourceData -> build/source_data/fig7/).

Panels: A tree | B-D ancestral K2P | E-L sweep landscapes + K2P (Crori1,
PuccoNC29_1, myrtle_rust up, myrtle_rust down) | M-P best-fit summary.
Nothing is written under /data2/chris/fungi or /home/chris/data.
"""

import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, "/data2/chris/manuscript/v2/scripts")
import ncomms_style  # noqa: E402

ncomms_style.use()

import matplotlib._mathtext as _mathtext  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FormatStrFormatter, MaxNLocator  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402

# Super/subscripts (bp^-1, log_10, 10^-11 tick labels) are drawn at
# SHRINK_FACTOR x the parent size; matplotlib's 0.7 puts a 6-pt tick
# exponent at 4.2 pt (< 5 pt minimum). 0.84 keeps every script >= 5.0 pt
# while tick labels stay 6 pt and axis labels 7 pt.
MATH_SCRIPT_SCALE = 0.84
if not hasattr(_mathtext, "SHRINK_FACTOR"):
    sys.exit("[error] matplotlib._mathtext.SHRINK_FACTOR not found "
             "(matplotlib API changed); cannot enforce >= 5 pt scripts")
_mathtext.SHRINK_FACTOR = MATH_SCRIPT_SCALE

# -----------------------------------------------------------------------------
# Print-size typography (points) and geometry (inches)
# -----------------------------------------------------------------------------
FS_LETTER = 8.0      # bold panel letters
FS_LABEL = 7.0       # axis labels
FS_TICK = 6.0        # tick labels
FS_LEGEND = 6.0      # legends
FS_CBAR = 6.0        # colourbar label + ticks
FS_SPECIES = 6.5     # species titles above the E-L pairs
PHYLO_POINTSIZE = 7.0  # R base pointsize; tree text = 0.85 * 7 = 5.95 pt

MM = 1 / 25.4
FIG_W = 180 * MM                 # 7.0866 in, exact
X_LEFT = 0.50                    # left edge of the first contour axes
X_RIGHT = FIG_W - 0.05           # right edge of the last K2P / P axes
CB_PAD, CB_W = 0.04, 0.06        # contour -> colourbar gap, colourbar width
CB_TO_K2P = 0.87                 # colourbar right edge -> K2P axes left
PAIR_GAP = 0.56                  # K2P axes right -> next contour axes left
ROW_H = 1.45                     # height of E-L axes
Q_H = 1.30                       # height of M-P axes
Q_LEFT, Q_GAP = 0.62, 0.62       # M-P left margin and inter-axes gap
TOP_PAD = 0.06                   # page top -> panel A top
BANNER_TO_ROW1 = 0.50            # D axes bottom -> row-1 axes top
ROW_GAP = 0.53                   # row-1 axes bottom -> row-2 axes top
ROW2_TO_Q = 0.47                 # row-2 axes bottom -> M-P axes top
BOTTOM_PAD = 0.58                # M-P axes bottom -> page bottom
PHYLO_X0 = 0.10                  # tree content left edge
PHYLO_PDF_W, PHYLO_PDF_H = 4.75, 2.10   # R cairo_pdf page (inches, 1:1)
LETTER_DX = 0.47                 # letter left of E/G/I/K/M-P axes (in)
LETTER_DX_K2P = 0.40             # letter left of F/H/J/L and B-D axes (in)
LETTER_DY = 0.05                 # letter baseline above axes top (in)

DOT_S = 2.5                      # scatter marker area (pt^2)
DOT_FC, DOT_EC, DOT_LW = "white", "#303030", 0.3
SWARM_HALF = 0.30                # max horizontal spread (bar width = 0.7)

K2P_DIST_IDX = 10
REF_COLOR = "#888888"
HIGHLIGHT_COLOR = "#FFD700"
PANEL_LABELS = "ABCDEFGHIJKLMNOP"

PHYLO_R_SUFFIX = ".contMap.pdf"
SCRIPTS = "/data2/chris/manuscript/v2/scripts"
MSCOREFONTS = "/home/chris/mamba/pkgs/mscorefonts-0.0.1-3/fonts"
READ_ONLY_ROOTS = ("/data2/chris/fungi", "/home/chris/data")


def log(msg, verbose=True):
    if verbose:
        print(f"[fig7] {msg}", file=sys.stderr)


# -----------------------------------------------------------------------------
# Inputs (identical to panel_plot.py build_rows)
# -----------------------------------------------------------------------------
def build_rows(root):
    """Per-scenario inputs, in panel_plot.py's canonical order
    [pucco, crori, mrust_up, mrust_down] (B-D colours are by this index)."""
    up = (f"{root}/mrust/long_up/insertion_rates_2.3143996843e-12_"
          f"deletion_rates_3.97591921428e-15_solo_ratio_85_length_bias_0")
    down_run = (f"{up}/short_down2/insertion_rates_1.2107737124e-12_"
                f"deletion_rates_1.23629969705e-11_solo_ratio_80_length_bias_8")
    pucco_run = (f"{root}/pucco/insertion_rates_1.32510498906e-11_"
                 f"deletion_rates_3.3507256332e-14_solo_ratio_40_length_bias_5")
    crori_run = (f"{root}/crori/insertion_rates_6.07522942927e-12_"
                 f"deletion_rates_1.14325788651e-12_solo_ratio_35_length_bias_3")
    return [
        {"key": "pucco", "row_short": "PuccoNC29_1", "color": "#4477AA",
         "species": r"$\mathit{Puccinia\ coronata}$ f. sp. $\mathit{avenae}$",
         "contour": f"{root}/pucco/composite_matrix_rms_0_0_10_1.tsv",
         "run_dir": pucco_run,
         "ref_ltr": f"{root}/pucco/PuccoNC29_1_ltr_kmer2ltr_dedup",
         "sim_ltr": f"{pucco_run}/gen260000_final.fasta_r1_ltr.tsv",
         "best": {"insertion_rate": 1.32510498906e-11,
                  "deletion_rate": 3.3507256332e-14,
                  "solo_ratio": 40, "length_bias": 5}},
        {"key": "crori", "row_short": "Crori1", "color": "#228833",
         "species": r"$\mathit{Cronartium\ ribicola}$",
         "contour": f"{root}/crori/composite_matrix_rms_0_0_10_1.tsv",
         "run_dir": crori_run,
         "ref_ltr": f"{root}/crori/Crori1_ltr_kmer2ltr_dedup",
         "sim_ltr": f"{crori_run}/gen5400000_final.fasta_r1_ltr.tsv",
         "best": {"insertion_rate": 6.07522942927e-12,
                  "deletion_rate": 1.14325788651e-12,
                  "solo_ratio": 35, "length_bias": 3}},
        {"key": "mrust_up", "row_short": "myrtle_rust ↑", "color": "#AA3377",
         "species": r"$\mathit{Austropuccinia\ psidii}$ (increase)",
         "contour": f"{root}/mrust/long_up/composite_matrix_rms_0_0_10_1.tsv",
         "run_dir": up,
         "ref_ltr": f"{root}/mrust/myrtle_rust_ltr_kmer2ltr_dedup3",
         "sim_ltr": f"{up}/gen113000000_final_r1_ltr2.tsv",
         "best": {"insertion_rate": 2.3143996843e-12,
                  "deletion_rate": 3.97591921428e-15,
                  "solo_ratio": 85, "length_bias": 0}},
        {"key": "mrust_down", "row_short": "myrtle_rust ↓", "color": "#EE6677",
         "species": r"$\mathit{Austropuccinia\ psidii}$ (decrease)",
         "contour": f"{up}/short_down2/composite_matrix_rms_0_0_10_1.tsv",
         "run_dir": down_run,
         "ref_ltr": f"{root}/mrust/myrtle_rust_ltr_kmer2ltr_dedup",
         "sim_ltr": f"{down_run}/blast_depth0_ltr.tsv",
         "best": {"insertion_rate": 1.2107737124e-12,
                  "deletion_rate": 1.23629969705e-11,
                  "solo_ratio": 80, "length_bias": 8}},
    ]


# -----------------------------------------------------------------------------
# Readers (unchanged logic)
# -----------------------------------------------------------------------------
def read_k2p_dists(path):
    """Finite, non-negative K2P_d values (column index 10; '#' lines skipped)."""
    vals = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) <= K2P_DIST_IDX:
                continue
            try:
                v = float(parts[K2P_DIST_IDX])
            except ValueError:
                continue
            if np.isfinite(v) and v >= 0:
                vals.append(v)
    return np.asarray(vals, dtype=float)


def fd_bin_count(data, xmin, xmax):
    """Tightened Freedman-Diaconis bin count, clamped to [20, 200]."""
    n = data.size
    if n < 2:
        return 20
    q75, q25 = np.percentile(data, [75, 25])
    iqr = q75 - q25
    span = xmax - xmin
    if span <= 0:
        return 20
    if iqr > 0:
        nb = max(1, int(np.ceil(span / (1.0 * iqr * n ** (-1.0 / 3.0)))))
    else:
        nb = max(1, int(np.ceil(1.0 + np.log2(n))))
    return max(20, min(nb, 200))


def read_k2p_bins(path):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["start", "end", "count"], comment="#")
    return (df["start"].to_numpy(float), df["end"].to_numpy(float),
            df["count"].to_numpy(float))


def load_solo_intact(run_dir):
    """Read-only version of compute_solo_intact_by_gen.load_or_build: reads the
    sidecar TSV and never builds/writes it (the run dirs are read-only)."""
    path = os.path.join(run_dir, "solo_intact_by_gen.tsv")
    if not os.path.exists(path):
        sys.exit(f"[error] missing {path}; build it with "
                 "plots/compute_solo_intact_by_gen.py (not done here: read-only)")
    t = pd.read_csv(path, sep="\t").sort_values("generation")
    return [(int(g), int(s), int(i))
            for g, s, i in t[["generation", "solo_count", "intact_count"]].values]


# -----------------------------------------------------------------------------
# Panel helpers
# -----------------------------------------------------------------------------
def add_axes_in(fig, x0, y0, w, h):
    """Axes from inch coordinates (origin bottom-left)."""
    fw, fh = fig.get_size_inches()
    return fig.add_axes([x0 / fw, y0 / fh, w / fw, h / fh])


def panel_letter(fig, ax, letter, dx_in):
    """8-pt bold letter, dx_in left of the axes, baseline LETTER_DY above it."""
    fw, fh = fig.get_size_inches()
    b = ax.get_position()
    fig.text(b.x0 - dx_in / fw, b.y1 + LETTER_DY / fh, letter,
             fontsize=FS_LETTER, fontweight="bold", ha="left", va="bottom")


def plot_contour(ax, tsv_path, best, zoom_decades=2.0, ceiling_pct=50.0,
                 levels=18, grid_res=180, log_eps=1e-300):
    """Composite landscape over log10(ins) x log10(del) — panel_plot.py logic."""
    df = pd.read_csv(tsv_path, sep="\t")
    for c in ("insertion_rate", "deletion_rate", "Composite"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.loc[(df.insertion_rate > 0) & (df.deletion_rate > 0)
                & df.Composite.notna()]
    if len(df) < 3:
        raise ValueError(f"{tsv_path}: <3 valid points after filtering")
    x_full = np.log10(np.maximum(df.insertion_rate.values, log_eps))
    y_full = np.log10(np.maximum(df.deletion_rate.values, log_eps))
    z_full = df.Composite.values.astype(float)
    bx = np.log10(max(best["insertion_rate"], log_eps))
    by = np.log10(max(best["deletion_rate"], log_eps))
    if zoom_decades and zoom_decades > 0:
        x_lo, x_hi = max(x_full.min(), bx - zoom_decades), min(x_full.max(), bx + zoom_decades)
        y_lo, y_hi = max(y_full.min(), by - zoom_decades), min(y_full.max(), by + zoom_decades)
    else:
        x_lo, x_hi, y_lo, y_hi = x_full.min(), x_full.max(), y_full.min(), y_full.max()
    xi, yi = np.linspace(x_lo, x_hi, grid_res), np.linspace(y_lo, y_hi, grid_res)
    Xi, Yi = np.meshgrid(xi, yi)
    pts = np.column_stack((x_full, y_full))
    Zi = griddata(pts, z_full, (Xi, Yi), method="linear")
    nan = np.isnan(Zi)
    if nan.any():
        Zi[nan] = griddata(pts, z_full, (Xi[nan], Yi[nan]), method="nearest")
    z_min, z_max = float(np.nanmin(Zi)), float(np.nanmax(Zi))
    if z_min == z_max:
        z_max = z_min + 1e-6
    if ceiling_pct and 0 < ceiling_pct < 100:
        z_top = max(min(float(np.percentile(z_full, ceiling_pct)), z_max), z_min + 1e-6)
    else:
        z_top = z_max
    Zi = np.minimum(Zi, z_top)
    cf = ax.contourf(Xi, Yi, Zi, levels=np.linspace(z_min, z_top, levels),
                     cmap="viridis_r")
    if x_lo <= bx <= x_hi and y_lo <= by <= y_hi:
        ax.plot(bx, by, marker="*", markersize=8, markerfacecolor=HIGHLIGHT_COLOR,
                markeredgecolor="black", markeredgewidth=0.5, zorder=10)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.set_xlabel(r"$\log_{10}$(insertion rate, bp$^{-1}$ gen$^{-1}$)", fontsize=FS_LABEL)
    ax.set_ylabel(r"$\log_{10}$(deletion rate, bp$^{-1}$ gen$^{-1}$)", fontsize=FS_LABEL,
                  labelpad=2)
    ax.tick_params(labelsize=FS_TICK)
    for s in ax.spines.values():          # filled field: keep the full frame
        s.set_visible(True)
    return cf


def plot_k2p_bins(ax, bins_path, color, label, xmax, show_x):
    """Pre-binned ancestral K2P histogram (B-D). Returns the plotted table."""
    starts, ends, counts = read_k2p_bins(bins_path)
    edges = np.concatenate([starts, [ends[-1]]])
    x_lo = float(starts.min())
    x_hi = float(xmax) if xmax is not None else float(ends.max())
    pad = 0.02 * (x_hi - x_lo)
    ax.hist(0.5 * (starts + ends), bins=edges, weights=counts, alpha=0.55,
            color=color, edgecolor="white", linewidth=0.2, label=label)
    ax.set_xlim(min(x_lo - pad, -pad), x_hi + pad)
    if show_x:
        ax.set_xlabel("K2P divergence", fontsize=FS_LABEL)
    else:
        ax.tick_params(axis="x", labelbottom=False)
    ax.set_ylabel("Count", fontsize=FS_LABEL, labelpad=2)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.tick_params(labelsize=FS_TICK)
    ax.legend(fontsize=FS_LEGEND, loc="upper right", handlelength=1.0,
              borderpad=0.1, borderaxespad=0.1)
    return pd.DataFrame({"bin_left": starts, "bin_right": ends,
                         "count": counts.astype(int)})


def plot_k2p_counts(ax, ref_path, sim_path, sim_color, xmax):
    """Reference vs simulated K2P histogram, raw counts (--k2p-counts mode).
    Returns (table as plotted, n_ref_total, n_sim_total)."""
    ref, sim = read_k2p_dists(ref_path), read_k2p_dists(sim_path)
    if ref.size == 0 or sim.size == 0:
        raise ValueError(f"empty K2P data: ref={ref.size} sim={sim.size}")
    all_vals = np.concatenate([ref, sim])
    x_lo = float(all_vals.min())
    x_hi = float(xmax) if xmax is not None else float(all_vals.max())
    if x_hi <= x_lo:
        x_hi = x_lo + 1e-6
    pad = 0.02 * (x_hi - x_lo)
    lo, hi = min(x_lo - pad, -pad), x_hi + pad
    nb = fd_bin_count(all_vals, lo, hi)
    n_ref, edges, _ = ax.hist(ref, bins=nb, range=(lo, hi), alpha=0.55,
                              color=REF_COLOR, linewidth=0, label="Reference")
    n_sim, _, _ = ax.hist(sim, bins=nb, range=(lo, hi), alpha=0.55,
                          color=sim_color, linewidth=0, label="Simulated")
    ax.set_xlim(lo, hi)
    ax.set_xlabel("K2P divergence", fontsize=FS_LABEL)
    ax.set_ylabel("Count", fontsize=FS_LABEL, labelpad=2)
    ax.tick_params(labelsize=FS_TICK)
    ax.legend(fontsize=FS_LEGEND, loc="upper right", handlelength=1.0,
              borderpad=0.1, borderaxespad=0.1)
    tab = pd.DataFrame({"bin_left": edges[:-1], "bin_right": edges[1:],
                        "reference_count": n_ref.astype(int),
                        "simulated_count": n_sim.astype(int)})
    return tab, ref.size, sim.size


# -----------------------------------------------------------------------------
# M-P: best-fit summary (values unchanged from panel_plot.py)
# -----------------------------------------------------------------------------
Q_PARAMS = ("insertion_rate", "deletion_rate", "solo_ratio", "length_bias")
Q_PARAM_LABELS = (r"Insertion rate (bp$^{-1}$ gen$^{-1}$)",
                  r"Deletion rate (bp$^{-1}$ gen$^{-1}$)",
                  "Solo:intact ratio (%)", "Length bias")
Q_PARAM_LOG = (True, True, False, False)
Q_TOP_PCT = 5.0
Q_SHARED_LOG = frozenset(("insertion_rate", "deletion_rate"))


def q_load_summary(tsv_path, best_params):
    """panel_plot.py _q_load_summary, unchanged, plus the matrix and the
    best-5% rows themselves: threshold = Composite.quantile(0.05) over all
    rows; keep Composite <= threshold; drop imputed failures
    (exp_genome_size == 0). Bar = the starred run (best_params)."""
    df = pd.read_csv(tsv_path, sep="\t")
    for c in (*Q_PARAMS, "Composite"):
        if c not in df.columns:
            raise ValueError(f"{tsv_path}: missing column {c!r}")
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["Composite", *Q_PARAMS])
    threshold = float(df["Composite"].quantile(Q_TOP_PCT / 100.0))
    top = df[df["Composite"] <= threshold]
    if "exp_genome_size" in top.columns:
        scored = top[pd.to_numeric(top["exp_genome_size"], errors="coerce").fillna(0) != 0]
        if len(scored) == 0:
            raise ValueError(f"{tsv_path}: top-5% set is all imputed failures")
        top = scored
    best = {p: float(best_params[p]) for p in Q_PARAMS}
    lo = {p: float(top[p].min()) for p in Q_PARAMS}
    hi = {p: float(top[p].max()) for p in Q_PARAMS}
    return best, lo, hi, top, df


def starred_mask(df, best):
    """Rows of a composite matrix equal to the starred parameter set."""
    return (np.isclose(df.insertion_rate, best["insertion_rate"], rtol=1e-9, atol=0)
            & np.isclose(df.deletion_rate, best["deletion_rate"], rtol=1e-9, atol=0)
            & (df.solo_ratio == best["solo_ratio"])
            & (df.length_bias == best["length_bias"]))


def swarm_offsets(ax, x0, values, seed):
    """Horizontal offsets (data units) for a beeswarm strip at x = x0.

    Greedy swarm in display space (marker diameter from DOT_S): each point,
    in value order, takes the offset closest to x0 that does not overlap an
    already placed point, within +/-SWARM_HALF. Points that cannot fit (dense
    ties) get a seeded uniform jitter inside the same band."""
    fig = ax.figure
    pt = 72.0 / fig.dpi
    vals = np.asarray(values, float)
    disp = ax.transData.transform(np.column_stack([np.full(vals.size, x0), vals]))
    y = disp[:, 1] * pt
    unit = (ax.transData.transform([[x0 + 1, vals[0]]])[0, 0] - disp[0, 0]) * pt
    d = np.sqrt(DOT_S) * 1.05                 # marker diameter (pt) + margin
    lim = SWARM_HALF * unit
    cands = [0.0]
    for k in range(1, int(lim / (d / 2)) + 1):
        cands += [k * d / 2, -k * d / 2]
    rng = np.random.default_rng(seed)
    placed_x, placed_y = [], []
    offs = np.zeros(vals.size)
    for i in np.argsort(y, kind="stable"):
        px, py = np.asarray(placed_x), np.asarray(placed_y)
        near = np.abs(py - y[i]) < d if py.size else np.zeros(0, bool)
        chosen = None
        for c in cands:
            if not near.any() or np.all((px[near] - c) ** 2 + (py[near] - y[i]) ** 2 >= d * d):
                chosen = c
                break
        if chosen is None:
            chosen = rng.uniform(-lim, lim)
        offs[i] = chosen
        placed_x.append(chosen)
        placed_y.append(y[i])
    return offs / unit


def plot_panel_q(axes, rows, verbose=False):
    """M-P bars (starred best fit) + min-max whiskers + overlaid points.
    Returns (points_long, summary, o_per_generation) DataFrames."""
    n = len(rows)
    bests = {p: np.zeros(n) for p in Q_PARAMS}
    lows = {p: np.zeros(n) for p in Q_PARAMS}
    highs = {p: np.zeros(n) for p in Q_PARAMS}
    points = {p: [] for p in Q_PARAMS}       # per scenario: (values, starred flags)
    n_pts = {p: [0] * n for p in Q_PARAMS}
    o_rows = []
    colors, labels = [], []
    for i, row in enumerate(rows):
        best, lo, hi, top, _df = q_load_summary(row["contour"], row["best"])
        star = starred_mask(top, row["best"]).to_numpy()
        for p in Q_PARAMS:
            bests[p][i], lows[p][i], highs[p][i] = best[p], lo[p], hi[p]
            points[p].append((top[p].to_numpy(float), star))
            n_pts[p][i] = len(top)
        # Panel O: realized solo:intact (%) of the best run; bar = endpoint
        # generation, whiskers = per-generation min/max, points = generations.
        table = load_solo_intact(row["run_dir"])
        ratios = [100.0 * s / it for _, s, it in table if it > 0]
        _g, s_end, i_end = table[-1]
        bar = 100.0 * s_end / i_end if i_end else 0.0
        bests["solo_ratio"][i] = bar
        lows["solo_ratio"][i] = min(ratios) if ratios else bar
        highs["solo_ratio"][i] = max(ratios) if ratios else bar
        o_vals = np.array([100.0 * s / it for _, s, it in table if it > 0])
        o_star = np.zeros(o_vals.size, bool)
        o_star[-1] = True                      # endpoint generation = bar
        points["solo_ratio"][i] = (o_vals, o_star)
        n_pts["solo_ratio"][i] = o_vals.size
        for (g, s, it) in table:
            o_rows.append({"scenario": row["row_short"], "generation": g,
                           "solo_count": s, "intact_count": it,
                           "solo_intact_pct": 100.0 * s / it if it else np.nan})
        log(f"M-P {row['row_short']}: best-5% n={len(top)} (starred in set: "
            f"{bool(star.any())}); O generations n={o_vals.size}", verbose)
        colors.append(row["color"])
        labels.append(row["row_short"])

    shared_lo = np.concatenate([lows[p] for p in Q_PARAMS if p in Q_SHARED_LOG])
    shared_hi = np.concatenate([highs[p] for p in Q_PARAMS if p in Q_SHARED_LOG])
    shared_h = np.concatenate([bests[p] for p in Q_PARAMS if p in Q_SHARED_LOG])
    shared_pos = shared_lo[shared_lo > 0]
    shared_bottom = (shared_pos.min() if shared_pos.size
                     else float(shared_h[shared_h > 0].min())) / 10.0
    shared_top = float(shared_hi.max()) * 3.0

    x = np.arange(n)
    err_kw = dict(ecolor="black", lw=0.6, capthick=0.6, zorder=4)
    long_rows, summ_rows = [], []
    for ax, p, lbl, logy in zip(axes, Q_PARAMS, Q_PARAM_LABELS, Q_PARAM_LOG):
        h, lo, hi = bests[p], lows[p], highs[p]
        yerr = np.clip(np.vstack([h - lo, hi - h]), 0.0, None)
        if logy:
            bottom, top_ = shared_bottom, shared_top
            ax.bar(x, h - bottom, bottom=bottom, color=colors, edgecolor="black",
                   linewidth=0.4, width=0.7, yerr=yerr, capsize=2,
                   error_kw=err_kw, zorder=2)
            ax.set_yscale("log")
            ax.set_ylim(bottom, top_)
        else:
            ax.bar(x, h, color=colors, edgecolor="black", linewidth=0.4,
                   width=0.7, yerr=yerr, capsize=2, error_kw=err_kw, zorder=2)
            y_top = max(float(hi.max()), float(h.max()))
            ax.set_ylim(0.0, y_top * 1.15 if y_top > 0 else 1.0)
        ax.set_xlim(-0.6, n - 0.4)
        for i in range(n):
            vals, star = points[p][i]
            offs = swarm_offsets(ax, x[i], vals, seed=1000 + 10 * i + Q_PARAMS.index(p))
            ax.scatter(x[i] + offs, vals, s=DOT_S, facecolor=DOT_FC,
                       edgecolor=DOT_EC, linewidth=DOT_LW, zorder=5,
                       clip_on=False)
            metric = "solo_intact_pct" if p == "solo_ratio" else p
            for v, st in zip(vals, star):
                long_rows.append({"panel": "MNOP"[Q_PARAMS.index(p)],
                                  "scenario": labels[i], "metric": metric,
                                  "value": v, "is_starred_best": bool(st)})
            summ_rows.append({"panel": "MNOP"[Q_PARAMS.index(p)],
                              "scenario": labels[i], "metric": metric,
                              "bar_value": h[i], "whisker_min": lo[i],
                              "whisker_max": hi[i], "n_points": n_pts[p][i]})
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=FS_TICK, rotation=45, ha="right",
                           rotation_mode="anchor")
        ax.set_ylabel(lbl, fontsize=FS_LABEL, labelpad=2)
        ax.tick_params(labelsize=FS_TICK)
        ax.grid(True, axis="y", linestyle=":", linewidth=0.4, alpha=0.5, zorder=0)
        ax.set_axisbelow(True)
    return pd.DataFrame(long_rows), pd.DataFrame(summ_rows), pd.DataFrame(o_rows)


# -----------------------------------------------------------------------------
# Panel A (R contMap -> cairo_pdf Arial) + vector embed
# -----------------------------------------------------------------------------
def newick_tips(newick):
    txt = open(newick).read()
    return re.findall(r"[(,]\s*([^():,;\s]+)\s*:", txt), txt.strip()


def build_phylo(args, verbose):
    """Run fig7_contmap.R in the build dir with a private fontconfig (Arial)."""
    tips, _ = newick_tips(args.phylo_newick)
    nwk_dir = os.path.dirname(os.path.abspath(args.phylo_newick))
    missing = [t for t in tips
               if not os.path.exists(os.path.join(nwk_dir, t + args.phylo_suffix + ".fai"))]
    if missing:   # the R script would write .fai files into the data dir
        sys.exit(f"[error] missing .fai for {missing} in {nwk_dir}; refusing to "
                 "let samtools write there (read-only)")
    os.makedirs(args.build_dir, exist_ok=True)
    conf = os.path.join(args.build_dir, "fonts.conf")
    with open(conf, "w") as fh:
        fh.write(
            '<?xml version="1.0"?>\n<!DOCTYPE fontconfig SYSTEM "fonts.dtd">\n'
            "<!-- Private fontconfig for the Fig. 7A R/cairo run only (installs\n"
            "     nothing): system + conda configs, plus the mscorefonts cache. -->\n"
            "<fontconfig>\n"
            '  <include ignore_missing="yes">/etc/fonts/fonts.conf</include>\n'
            f'  <include ignore_missing="yes">{os.path.dirname(os.path.dirname(args.rscript))}'
            "/etc/fonts/fonts.conf</include>\n"
            f"  <dir>{MSCOREFONTS}</dir>\n"
            f"  <cachedir>{os.path.join(args.build_dir, 'fontconfig-cache')}</cachedir>\n"
            "</fontconfig>\n")
    stem = os.path.basename(args.phylo_pdf)[: -len(PHYLO_R_SUFFIX)]
    cmd = [args.rscript, args.phylo_script, "--suffix", args.phylo_suffix,
           "--newick", os.path.abspath(args.phylo_newick), "--out", stem,
           "--pdf_width", f"{PHYLO_PDF_W}", "--pdf_height", f"{PHYLO_PDF_H}",
           "--pointsize", f"{PHYLO_POINTSIZE}", "--tree_lwd", "2.4",
           "--bar_lwd", "2.0", "--circle_cex", "1.3"]
    env = dict(os.environ, FONTCONFIG_FILE=conf,
               PATH=os.path.dirname(args.rscript) + os.pathsep + os.environ.get("PATH", ""))
    log("R: " + " ".join(cmd), verbose)
    res = subprocess.run(cmd, cwd=args.build_dir, env=env, capture_output=True, text=True)
    with open(os.path.join(args.build_dir, stem + ".R.log"), "w") as fh:
        fh.write(res.stdout + res.stderr)
    if res.returncode != 0 or not os.path.exists(args.phylo_pdf):
        sys.exit(f"[error] fig7_contmap.R failed (exit {res.returncode}):\n{res.stderr[-2000:]}")


def check_phylo_fonts(pdf):
    import fitz
    doc = fitz.open(pdf)
    fonts = doc[0].get_fonts(full=True)
    doc.close()
    bad = [f[3] for f in fonts if "Arial" not in f[3] or f[1] in ("n/a", "")]
    if not fonts or bad:
        sys.exit(f"[error] panel A fonts not embedded Arial: {bad or 'no text'} "
                 "(fontconfig did not resolve Arial?)")


def phylo_content_bbox(phylo_pdf, pad_pt=1.0):
    """Tight bbox of ink (drawings + text) on page 0 (fitz coords)."""
    import fitz
    src = fitz.open(phylo_pdf)
    try:
        page = src[0]
        bb = None
        for dr in page.get_drawings():
            r = dr.get("rect")
            fill = dr.get("fill")
            # cairo_pdf paints a white page background: not ink
            if (r is None or (dr.get("color") is None and fill is not None
                              and all(c > 0.999 for c in fill))):
                continue
            bb = r if bb is None else bb | r
        for tb in page.get_text("blocks"):
            r = fitz.Rect(tb[:4])
            bb = r if bb is None else bb | r
        prect = fitz.Rect(page.rect)
    finally:
        src.close()
    return fitz.Rect(bb.x0 - pad_pt, bb.y0 - pad_pt, bb.x1 + pad_pt, bb.y1 + pad_pt) & prect


def embed_phylo(out_pdf, phylo_pdf, clip, x0_in, top_in):
    """Place the clipped tree at 1:1 scale with its top-left at (x0, top)."""
    import fitz
    rect = fitz.Rect(x0_in * 72, top_in * 72, x0_in * 72 + clip.width,
                     top_in * 72 + clip.height)
    doc, src = fitz.open(out_pdf), fitz.open(phylo_pdf)
    try:
        doc[0].show_pdf_page(rect, src, 0, clip=clip)
        tmp = out_pdf + ".tmp"
        doc.save(tmp, garbage=3, deflate=True)
    finally:
        doc.close()
        src.close()
    os.replace(tmp, out_pdf)
    return rect


# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
def make_figure(rows, ancestral_specs, args, clip, verbose):
    """Lay out page 1 at print size. Returns (fig, source-data dict, A rect)."""
    wc = wk = ((X_RIGHT - X_LEFT - PAIR_GAP) / 2 - (CB_PAD + CB_W + CB_TO_K2P)) / 2
    phylo_h = clip.height / 72
    phylo_w = clip.width / 72
    fig_h = (TOP_PAD + phylo_h + BANNER_TO_ROW1 + ROW_H + ROW_GAP + ROW_H
             + ROW2_TO_Q + Q_H + BOTTOM_PAD)
    fig = plt.figure(figsize=(FIG_W, fig_h))
    sd = {}

    # ---- banner: A (tree, embedded after save) + B-D -----------------------
    banner_top = fig_h - TOP_PAD
    banner_bot = banner_top - phylo_h
    fig.text(0.03 / FIG_W, banner_top / fig_h, "A", fontsize=FS_LETTER,
             fontweight="bold", ha="left", va="top")
    bd_x0 = PHYLO_X0 + phylo_w + 0.52
    bd_w = X_RIGHT - bd_x0
    bd_top = banner_top - 0.13          # room for letter B
    bd_gap = 0.17
    bd_h = (bd_top - banner_bot - 2 * bd_gap) / 3
    sd["BD"] = []
    for j, spec in enumerate(ancestral_specs):
        y0 = bd_top - (j + 1) * bd_h - j * bd_gap
        ax = add_axes_in(fig, bd_x0, y0, bd_w, bd_h)
        tab = plot_k2p_bins(ax, spec["path"], spec["color"], spec["label_prefix"],
                            args.k2p_xmax, show_x=(j == 2))
        tab.insert(0, "genome", spec["label_prefix"])
        sd["BD"].append(tab)
        panel_letter(fig, ax, PANEL_LABELS[1 + j], LETTER_DX_K2P)

    # ---- rows 1-2: contour + K2P pairs --------------------------------------
    row1_top = banner_bot - BANNER_TO_ROW1
    row_tops = (row1_top, row1_top - ROW_H - ROW_GAP)
    pairs = [(rows[0], rows[1], "EFGH"), (rows[2], rows[3], "IJKL")]
    sd["contour"], sd["k2p"] = {}, {}
    for (ra, rb, letters), top in zip(pairs, row_tops):
        y0 = top - ROW_H
        x = X_LEFT
        for r, (lc, lk) in ((ra, letters[:2]), (rb, letters[2:])):
            axc = add_axes_in(fig, x, y0, wc, ROW_H)
            cf = plot_contour(axc, r["contour"], r["best"], args.zoom_decades,
                              args.ceiling_pct)
            cax = add_axes_in(fig, x + wc + CB_PAD, y0, CB_W, ROW_H)
            cb = fig.colorbar(cf, cax=cax)
            cb.set_label("Composite score", fontsize=FS_CBAR, labelpad=2)
            cb.ax.tick_params(labelsize=FS_CBAR, width=0.5, length=2)
            cb.locator = MaxNLocator(nbins=5)
            cb.formatter = FormatStrFormatter("%.2f")
            cb.update_ticks()
            cb.outline.set_linewidth(0.5)
            panel_letter(fig, axc, lc, LETTER_DX)
            fw, fh = fig.get_size_inches()
            fig.text(x / fw, (top + LETTER_DY) / fh, r["species"],
                     fontsize=FS_SPECIES, ha="left", va="bottom")
            xk = x + wc + CB_PAD + CB_W + CB_TO_K2P
            axk = add_axes_in(fig, xk, y0, wk, ROW_H)
            tab, n_ref, n_sim = plot_k2p_counts(axk, r["ref_ltr"], r["sim_ltr"],
                                                r["color"], args.k2p_xmax)
            panel_letter(fig, axk, lk, LETTER_DX_K2P)
            sd["contour"][lc] = r
            sd["k2p"][lk] = (r, tab, n_ref, n_sim)
            log(f"{lc}/{lk} {r['row_short']}: K2P ref n={n_ref} sim n={n_sim} "
                f"bins={len(tab)}", verbose)
            x = xk + wk + PAIR_GAP

    # ---- row 3: M-P -----------------------------------------------------------
    q_w = (X_RIGHT - Q_LEFT - 3 * Q_GAP) / 4
    q_y0 = row_tops[1] - ROW_H - ROW2_TO_Q - Q_H
    q_axes = [add_axes_in(fig, Q_LEFT + j * (q_w + Q_GAP), q_y0, q_w, Q_H)
              for j in range(4)]
    sd["q_points"], sd["q_summary"], sd["o_gen"] = plot_panel_q(q_axes, rows, verbose)
    for ax, letter in zip(q_axes, "MNOP"):
        panel_letter(fig, ax, letter, LETTER_DX + 0.08)
    return fig, sd, (PHYLO_X0, TOP_PAD)


# -----------------------------------------------------------------------------
# Source Data
# -----------------------------------------------------------------------------
def write_source_data(sd, rows, args, verbose):
    from srcdata import SourceData
    out = SourceData("fig7")
    stem = os.path.join(args.build_dir,
                        os.path.basename(args.phylo_pdf)[: -len(PHYLO_R_SUFFIX)])

    tips = pd.read_csv(stem + ".tip_genome_sizes.tsv", sep="\t")
    tips["genome_size_Mb"] = tips["genome_bp"] / 1e6
    out.add("Fig. 7A (tip genome sizes)", f"Genome size of the {len(tips)} tips (Mb; sum of FASTA "
            "sequence lengths), the trait mapped onto the tree (n = 18 genomes).",
            tips[["tip", "genome_size_Mb", "genome_bp"]])
    _, nwk = newick_tips(args.phylo_newick)
    out.add("Fig. 7A (newick)", "Time-scaled species tree (Newick; branch lengths in "
            "millions of generations, root-to-tip 266).", pd.DataFrame({"newick": [nwk]}))
    nodes = pd.read_csv(stem + ".ancestral_genome_sizes.tsv", sep="\t")
    nodes = pd.DataFrame({
        "node": nodes["node"], "representative_tips": nodes["representative_tips"],
        "ancestral_size_Mb": nodes["bp_est"] / 1e6,
        "CI95_lower_Mb": nodes["bp_CI_lower"] / 1e6,
        "CI95_upper_Mb": nodes["bp_CI_upper"] / 1e6})
    out.add("Fig. 7A (ancestral nodes)", f"Reconstructed ancestral genome size at the {len(nodes)} "
            "internal nodes (phytools fastAnc on log10 size, back-transformed; Mb, "
            "95% CI).", nodes)

    bd = pd.concat(sd["BD"], ignore_index=True)
    tot = bd.groupby("genome", sort=False)["count"].sum().to_dict()
    out.add("Fig. 7B–D", "Ancestral LTR-RT K2P-divergence distribution of each "
            "focal genome: per-bin counts reconstructed at the ancestral node (phytools "
            "fastAnc, rounded); 50 bins = [0, 1e-7] plus 49 bins of width 0.00306 to "
            "0.15. Totals: " + ", ".join(f"{k} n = {v:,}" for k, v in tot.items()) + ".",
            bd)

    for letter in "EFGHIJKL":
        if letter in sd["contour"]:
            r = sd["contour"][letter]
            best, _lo, _hi, top, df = q_load_summary(r["contour"], r["best"])
            blk = pd.DataFrame({
                "insertion_rate": df.insertion_rate, "deletion_rate": df.deletion_rate,
                "solo_ratio": df.solo_ratio, "length_bias": df.length_bias,
                "Composite": df.Composite,
                "completed": pd.to_numeric(df.exp_genome_size, errors="coerce").fillna(0) != 0,
                "in_best5pct": df.index.isin(top.index),
                "is_starred": starred_mask(df, r["best"]).to_numpy()})
            n_done = int(blk.completed.sum())
            if int(blk.is_starred.sum()) != 1:
                sys.exit(f"[error] starred run not unique in {r['contour']}")
            out.add(f"Fig. 7{letter}", f"{r['row_short']} parameter sweep, every "
                    f"simulation (n = {len(blk):,}; rates in bp^-1 gen^-1; Composite "
                    f"score, lower = better). {n_done} completed; the other "
                    f"{len(blk) - n_done:,} did not complete and carry the penalty "
                    f"score (max + 0.1). The landscape is interpolated from all "
                    f"{len(blk):,} points; in_best5pct = the n = "
                    f"{int(blk.in_best5pct.sum())} simulations behind panels M, N, P.",
                    blk)
        else:
            r, tab, n_ref, n_sim = sd["k2p"][letter]
            out.add(f"Fig. 7{letter}", f"{r['row_short']} K2P divergence of LTR-RTs, "
                    f"counts per bin as plotted ({len(tab)} bins of width "
                    f"{tab.bin_right.iloc[0] - tab.bin_left.iloc[0]:.5f}); reference "
                    f"n = {n_ref:,} LTR-RTs ({int(tab.reference_count.sum()):,} within "
                    f"the plotted range <= {tab.bin_right.iloc[-1]:.3f}), simulated n = "
                    f"{n_sim:,} ({int(tab.simulated_count.sum()):,} within range).", tab)

    pts, summ = sd["q_points"], sd["q_summary"]
    npt = summ.set_index(["panel", "scenario"])["n_points"]
    desc_n = "; ".join(
        f"{p}: " + ", ".join(f"{s} {npt[(p, s)]}" for s in summ.scenario.unique())
        for p in "MNOP")
    out.add("Fig. 7M–P (points)", "Every overlaid point. M, N (bp^-1 gen^-1) and P: best-5% "
            "simulations (Composite <= 5th percentile of all simulations in the sweep, "
            "non-completed simulations excluded); O: realized solo:intact ratio (%) of "
            "the best-fit simulation at each saved generation (is_starred_best marks the "
            f"final generation = bar). n per panel and scenario: {desc_n}.", pts)
    out.add("Fig. 7M–P (summary)", "Bar = starred best-fit simulation (O: its final generation); "
            "whiskers = min to max of the overlaid points (centre = bar value, not a "
            "mean); n_points = number of points.", summ)
    out.add("Fig. 7O", "Per-generation solo and intact LTR-RT counts of the best "
            "run and realized solo:intact ratio (%) (n generations: "
            + ", ".join(f"{s} {g}" for s, g in
                        sd["o_gen"].groupby("scenario", sort=False).size().items())
            + ").", sd["o_gen"])
    out.save()
    log(f"source data -> {out.dir} ({len(out.blocks)} blocks)", verbose)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Build final Fig. 7 (single page, 180 mm, Arial TrueType) "
                    "+ its Source Data blocks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--root", default="/data2/chris/fungi/best",
                   help="read-only data root (pucco/, crori/, mrust/, *.bins)")
    p.add_argument("-o", "--output", default="/data2/chris/manuscript/v2/fig7/fig7.pdf",
                   help="output PDF")
    p.add_argument("--build-dir", default="/data2/chris/manuscript/v2/build/fig7",
                   help="scratch dir for the R tree PDF, fonts.conf, logs")
    p.add_argument("--phylo-regen", action="store_true",
                   help="re-run fig7_contmap.R even if the tree PDF exists")
    p.add_argument("--phylo-script", default=f"{SCRIPTS}/fig7_contmap.R",
                   help="R contMap script (Arial/cairo copy)")
    p.add_argument("--phylo-newick", default="/home/chris/data/fungi/v2/subset_time.nwk",
                   help="time-scaled Newick (FASTA .fai files alongside)")
    p.add_argument("--phylo-suffix", default=".fa", help="FASTA suffix per tip")
    p.add_argument("--rscript",
                   default=os.path.join(os.path.dirname(sys.executable), "Rscript"),
                   help="Rscript binary (needs phytools, optparse, cairo)")
    p.add_argument("--k2p-xmax", type=float, default=0.25,
                   help="common K2P axis upper bound")
    p.add_argument("--zoom-decades", type=float, default=2.0,
                   help="contour view half-width (log10 units) around the star")
    p.add_argument("--ceiling-pct", type=float, default=50.0,
                   help="contour colour ceiling percentile of Composite")
    p.add_argument("--no-source-data", action="store_true",
                   help="skip writing build/source_data/fig7/")
    p.add_argument("-v", "--verbose", action="store_true", help="per-step progress")
    a = p.parse_args()
    a.phylo_pdf = os.path.join(a.build_dir, "fig7_tree" + PHYLO_R_SUFFIX)
    return a


def main():
    args = parse_args()
    v = args.verbose
    for d in (args.output, args.build_dir):
        if os.path.abspath(d).startswith(READ_ONLY_ROOTS):
            sys.exit(f"[error] refusing to write under a read-only data root: {d}")
    if not os.path.isdir(args.root):
        sys.exit(f"[error] --root not found: {args.root}")
    rows = build_rows(args.root)
    ancestral_specs = [
        {"label_prefix": "PuccoNC29_1", "color": rows[0]["color"],
         "path": os.path.join(args.root, "pucco.ancestral.LTR.bins")},
        {"label_prefix": "Crori1", "color": rows[1]["color"],
         "path": os.path.join(args.root, "Crori1.ancestral.LTR.bins")},
        {"label_prefix": "myrtle_rust", "color": rows[2]["color"],
         "path": os.path.join(args.root, "mrust.ancestral.LTR.bins")},
    ]
    missing = [f"{r['row_short']} {k}: {r[k]}" for r in rows
               for k in ("contour", "ref_ltr", "sim_ltr", "run_dir") if not os.path.exists(r[k])]
    missing += [s["path"] for s in ancestral_specs if not os.path.exists(s["path"])]
    if missing:
        sys.exit("[error] missing inputs:\n  " + "\n  ".join(missing))
    try:
        import fitz  # noqa: F401
    except ImportError:
        sys.exit("[error] PyMuPDF (fitz) is required for the vector tree embed")

    log("start", True)
    if args.phylo_regen or not os.path.exists(args.phylo_pdf):
        build_phylo(args, v)
    check_phylo_fonts(args.phylo_pdf)
    clip = phylo_content_bbox(args.phylo_pdf)
    log(f"tree PDF {args.phylo_pdf}: content {clip.width / 72:.2f} x "
        f"{clip.height / 72:.2f} in (embedded 1:1)", v)

    # Display order: Crori1, PuccoNC29_1, myrtle ↑, myrtle ↓ (as panel_plot.py)
    display_rows = [rows[1], rows[0], rows[2], rows[3]]
    fig, sd, (px, ptop) = make_figure(display_rows, ancestral_specs, args, clip, v)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    fig.savefig(args.output)
    w_in, h_in = fig.get_size_inches()
    plt.close(fig)
    rect = embed_phylo(args.output, args.phylo_pdf, clip, px, ptop)
    log(f"panel A embedded at {tuple(round(c, 1) for c in rect)} pt", v)
    if not args.no_source_data:
        write_source_data(sd, display_rows, args, v)
    print(f"[ok] wrote {args.output}: {w_in * 25.4:.1f} x {h_in * 25.4:.1f} mm")


if __name__ == "__main__":
    main()
