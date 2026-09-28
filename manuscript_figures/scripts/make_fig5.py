#!/usr/bin/env python3
"""Build final Fig. 5 (A. psidii TE composition approximation) for Nature Communications.

Rebuilds slide 14 of fig5.pptx (V21.1.1) at 180 mm width, fully vector (default --mode vector):

  A  TE loci per superfamily        original: RepeatMasker copy number (1,181,685 loci)
                                    simulated: features in the replica GFF (1,251,338 loci)
  B  TE bases per superfamily       original: summed hit length (981,930,851 bp)
                                    simulated: summed feature length (793,915,337 bp)
  C  divergence (p-distance) = 1 - identity to the library consensus, both genomes
  D  integrity = element length / library consensus length (values > 1 not plotted)
     original: length_ratio of every RepeatMasker hit (1,045,728 of 1,181,685 loci)
     simulated: resampled per locus from the per-family integrity lists with the simulator's own
     rule (TE_sim_random_insertion.fragment_m2, random.seed(--seed)); the bar at integrity 1.0 is
     set to 70,504 loci read off the published panel (the only number not recomputed here).

Inputs: TEgenomeSimulator mode-2 scan of the original genome (--scan; RepeatMasker 4.1.8/4.1.5,
NOT 4.1.9) and the replica annotation GFF (--gff). Stacked-bar bins follow ggplot's geom_histogram
(binwidth 0.02, right-closed (c - 0.01, c + 0.01]); fills are the published palette at alpha 0.8,
pre-blended over white. Also writes the Source Data blocks and a preview PNG.

--mode bitmap rebuilds the earlier interim composite instead: the published 1200-dpi ggplot panels
cropped to the data marks, with vector axes/labels laid over them (kept for reference).
Content change vs the slide: simulated loci total typo 1,251,335 -> 1,251,338.
"""
import argparse
import ast
import io
import json
import math
import random
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ncomms_style  # noqa: E402
from srcdata import SourceData  # noqa: E402

ROOT = HERE.parent
PPTX = Path("/data2/chris/manuscript/fig5/fig5.pptx")
SCAN = ROOT / "build/fig5/tegs_scan418/TEgenomeSimulator_psidii_scan418_result"
GFF = Path("/data2/chris/fungi/best/forward_march/psidii_run2_repeat_annotation_out_final.gff")
Image.MAX_IMAGE_PIXELS = None

W, H = 180.0, 110.0  # figure size, mm
PPI = 1200           # bitmap mode: resolution at final size (never above source)
FS_LETTER, FS_TITLE, FS_AXIS, FS_TICK, FS_LEG, FS_PCT = 8, 6.5, 6.5, 6, 6, 5.5
FRAME, TICKC, STRIP_FILL = "#A6A6A6", "#333333", "#DADADA"

SF = ["DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron", "LINE/Deceiver", "LINE/Tad1",
      "LTR/Copia", "LTR/Gypsy", "LTR/unknown", "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", "MITE/DTT", "unknown"]
PAL = ["#CC9999", "#CC6666", "#993344", "#996699", "#996633", "#666699", "#336699", "#0099AA", "#99CCCC",
       "#006655", "#669966", "#669922", "#99CC99", "#999966", "#CCCC99", "#CCCC77", "#CCCC11"]
ALPHA = 0.8
RGB = np.array([[int(h[i:i + 2], 16) for i in (1, 3, 5)] for h in PAL], float)
BLEND = ALPHA * RGB + (1 - ALPHA) * 255                      # bitmap-mode colour matching
FILL = dict(zip(SF, [tuple(c / 255) for c in BLEND]))        # published look, no PDF transparency

BINW = 0.02
FUZZ = 1e-8 * BINW
PCT_MIN = 5.0        # published rule: label slices > 5 %
PCT_R = 0.66         # radius (units of R) of the pie % labels
ONE_BIN_SIM = 70504  # D simulated, bin centred on 1.0: value read off the published panel

# slide-14 picture -> panel (bitmap mode; verified against ppt/slides/_rels/slide14.xml.rels)
IMG = {"A_orig": "image26", "A_sim": "image23", "B_orig": "image24", "B_sim": "image25",
       "C_orig": "image22", "C_sim": "image21", "D_orig": "image28", "D_sim": "image27"}
C_XLIM, C_YLIM = (-0.03, 0.60), (0, 120000)
D_XLIM, D_YLIM = (-0.03, 1.03), (0, 90000)
TOTALS = {"A_orig": 1181685, "A_sim": 1251338, "B_orig": 981930851, "B_sim": 793915337}
# layout (mm from the figure's top-left)
DIA, PIE_Y = 28.0, 12.0
PIE_X = {"A_orig": 10.5, "A_sim": 40.5, "B_orig": 82.5, "B_sim": 112.5}
CY0, CH, CW = 50.5, 48.5, 27.0
CX = {"C_orig": 11.0, "C_sim": 40.0}
DX0, DW, DH = 84.0, 60.0, 22.6
LEG_X, LEG_Y, LEG_ROW = 148.0, 33.0, 3.05


def log(msg, verbose=True):
    if verbose:
        print(msg, file=sys.stderr, flush=True)


def bink(x):
    """ggplot geom_histogram bin index k (bin centre = k * BINW, bin = (centre-0.01, centre+0.01])."""
    return np.ceil((np.asarray(x, float) - FUZZ - BINW / 2) / BINW).astype(int)


# ------------------------------------------------------------------ exact data
def exact_tables(scan, gff, seed, verbose=False):
    """Per-superfamily loci/bases and per-bin, per-superfamily counts for all four panels."""
    scan = Path(scan)
    fam = pd.read_csv(scan / "summarise_repeatmasker_out_family.txt", sep="\t",
                      usecols=["original_family", "superfamily", "copynumber", "total_length", "integrity_lst"])
    po = pd.read_csv(scan / "processed_output.temp", sep="\t", usecols=["superfamily", "identity", "length_ratio"])
    log(f"  original: {len(fam)} families, {len(po)} RepeatMasker hits", verbose)
    t = {"A_orig": fam.groupby("superfamily").copynumber.sum(),
         "B_orig": fam.groupby("superfamily").total_length.sum()}
    t["C_orig"] = po.assign(k=bink(1 - po.identity)).groupby(["k", "superfamily"]).size()
    keep = po[po.length_ratio <= 1]
    t["D_orig"] = keep.assign(k=bink(keep.length_ratio)).groupby(["k", "superfamily"]).size()
    log(f"  original integrity: {len(keep)} loci plotted, {len(po) - len(keep)} with integrity > 1 dropped", verbose)

    rows, sf_of, loci, bases, cdiv = Counter(), {}, Counter(), Counter(), Counter()
    with open(gff) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            at = f[8]
            family = at.split(";", 1)[0][3:].rsplit("_TE", 1)[0]
            a = dict(kv.split("=", 1) for kv in at.split(";") if "=" in kv)
            sf = a["Classification"]
            rows[family] += 1
            sf_of[family] = sf
            loci[sf] += 1
            bases[sf] += int(f[4]) - int(f[3]) + 1
            cdiv[(int(math.ceil((1 - float(a["Identity"]) - FUZZ - BINW / 2) / BINW)), sf)] += 1
    t["A_sim"] = pd.Series(loci)
    t["B_sim"] = pd.Series(bases)
    t["C_sim"] = pd.Series(cdiv)
    log(f"  simulated: {sum(rows.values())} loci in {len(rows)} families, {sum(bases.values())} bp", verbose)

    # D simulated: one integrity per locus, drawn from its family's integrity list (fragment_m2)
    lists = {r.original_family: ast.literal_eval(r.integrity_lst) for r in fam.itertuples()}
    random.seed(seed)
    draw, dropped = Counter(), 0
    for family in fam.original_family:              # draw order = summary-file order
        n = rows.get(family, 0)
        if not n:
            continue
        lst, sf = lists[family], sf_of[family]
        for _ in range(n):
            val = random.choice(lst)
            if val > 1:
                dropped += 1
                continue
            draw[(int(math.ceil((val - FUZZ - BINW / 2) / BINW)), sf)] += 1
    log(f"  simulated integrity (seed {seed}): {sum(draw.values())} loci plotted, {dropped} > 1 dropped", verbose)
    d_sim = pd.Series(draw)
    # single manual value: the bar at 1.0 is set to the published total, split in the published proportions
    pub = pd.read_csv(ROOT / "build/fig5/inputs/psidii_fig5_PNG_digitized.tsv", sep="\t")
    pub = pub[(pub.panel == "D_integrity") & (pub.dataset == "simulated") &
              (pd.to_numeric(pub.bin_center, errors="coerce") == 1.0)]
    if pub.empty:
        sys.exit("ERROR: no published bin-1.0 proportions in psidii_fig5_PNG_digitized.tsv")
    share = pub.set_index("superfamily").value_est / pub.value_est.sum()
    exact = share * ONE_BIN_SIM
    n = np.floor(exact).astype(int)
    for sf in (exact - n).sort_values(ascending=False).index[:ONE_BIN_SIM - int(n.sum())]:
        n[sf] += 1
    k1 = bink(1.0)
    one = pd.Series({(k1, sf): int(c) for sf, c in n.items() if c})
    log(f"  simulated integrity: bin 1.0 {int(d_sim[k1].sum())} -> {int(one.sum())} (published panel), "
        f"split over {len(one)} superfamilies", verbose)
    t["D_sim"] = pd.concat([d_sim.drop(index=k1, level=0), one])
    for key, want in TOTALS.items():
        got = int(t[key].sum())
        if got != want:
            sys.exit(f"ERROR: {key} total {got:,} != expected {want:,}")
    return t


def pct_check(t, cal, verbose=False):
    """The four pies must reproduce the published percentages."""
    for key, name in (("A_orig", "image26"), ("A_sim", "image23"), ("B_orig", "image24"), ("B_sim", "image25")):
        s = t[key] / t[key].sum() * 100
        for sf, txt in cal["pie_labels_in_png"][name].items():
            if f"{s[sf]:.2f}%" != txt:
                sys.exit(f"ERROR: {key} {sf}: {s[sf]:.2f}% != published {txt}")
        log(f"  {key}: published % labels reproduced ({', '.join(cal['pie_labels_in_png'][name].values())})", verbose)


# ------------------------------------------------------------------ bitmap mode
def load_pngs(pptx):
    """Read the 8 slide-14 PNGs from the pptx (read-only) as uint8 RGB arrays (lazy)."""
    z = zipfile.ZipFile(pptx)
    rels = z.read("ppt/slides/_rels/slide14.xml.rels").decode()
    missing = set(IMG.values()) - set(re.findall(r"media/(image\d+)\.png", rels))
    if missing:
        sys.exit(f"ERROR: {missing} not referenced by slide 14 of {pptx}")
    return lambda name: np.asarray(Image.open(io.BytesIO(z.read(f"ppt/media/{name}.png"))).convert("RGB"))


def resample(arr, w_mm, h_mm, ppi=PPI):
    """Downsample (Lanczos) so the bitmap is ~ppi at its final printed size (never upsample)."""
    tw = min(arr.shape[1], int(round(w_mm / 25.4 * ppi)))
    th = min(arr.shape[0], int(round(h_mm / 25.4 * ppi)))
    if (tw, th) == (arr.shape[1], arr.shape[0]):
        return arr
    return np.asarray(Image.fromarray(arr).resize((tw, th), Image.LANCZOS))


def hist_crop(a, cal, xlim, ylim):
    """Crop the plot-area interior to xlim x ylim (data units); pad white outside the ggplot panel.
    Returns (crop, extent) with extent = exact data coordinates of the crop's outer pixel edges."""
    x0, x1, y0, yt, tv = cal["x0"], cal["x1"], cal["y0"], cal["yt"], cal["tv"]
    ci0, ci1, ri0, ri1 = cal["interior"]
    sx, sy = x1 - x0, (y0 - yt) / tv          # px per x unit, px per count
    c_lo = int(round(x0 + 0.5 + xlim[0] * sx)); c_hi = int(round(x0 + 0.5 + xlim[1] * sx))
    r_lo = int(round(y0 + 0.5 - ylim[1] * sy)); r_hi = int(round(y0 + 0.5 - ylim[0] * sy))
    out = np.full((r_hi - r_lo, c_hi - c_lo, 3), 255, np.uint8)
    a0, a1 = max(c_lo, ci0), min(c_hi, ci1 + 1)
    b0, b1 = max(r_lo, ri0), min(r_hi, ri1 + 1)
    out[b0 - r_lo:b1 - r_lo, a0 - c_lo:a1 - c_lo] = a[b0:b1, a0:a1]
    ext = ((c_lo - x0 - 0.5) / sx, (c_hi - x0 - 0.5) / sx, (y0 + 0.5 - r_hi) / sy, (y0 + 0.5 - r_lo) / sy)
    return out, ext


def _ring_separators(a, cx, cy, R, rf, n=720000):
    """Angles (rad, clockwise from 12 o'clock) of the midpoints of white runs on a ring of radius rf*R."""
    t = (np.arange(n) + 0.5) / n * 2 * np.pi
    r = rf * R
    px = a[np.round(cy - r * np.cos(t) - 0.5).astype(int), np.round(cx + r * np.sin(t) - 0.5).astype(int)]
    wh = (px > 235).all(1)
    k0 = int(np.flatnonzero(~wh)[0])
    d = np.diff(np.r_[0, np.roll(wh, -k0).astype(int), 0])
    st, en = np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1
    return np.mod(((st + en) / 2 + k0 + 0.5) / n * 2 * np.pi, 2 * np.pi)


def pie_clean(a, cal, verbose=False):
    """Crop the pie disc and paint out the raster % labels (text pixels -> flat slice colour).
    Returns (crop, extent, labels); labels = {colour index: (x, y) of the label ink centre in units
    of R, origin at the disc centre, y up}. Separator lines are never touched."""
    cx, cy, R = cal["cx"], cal["cy"], cal["R"]
    pad = int(0.02 * R)
    r0, r1 = int(cy - R) - pad, int(cy + R) + pad + 1
    c0, c1 = int(cx - R) - pad, int(cx + R) + pad + 1
    crop = a[r0:r1, c0:c1].copy()
    rays = np.sort(np.concatenate([_ring_separators(a, cx, cy, R, rf) for rf in (0.45, 0.97)]))
    groups = []
    for r in rays:  # same separator seen on both rings -> one ray
        if groups and r - groups[-1][-1] < np.radians(0.3):
            groups[-1].append(r)
        else:
            groups.append([r])
    rays = np.array([np.mean(g) for g in groups])
    nr = len(rays)

    def geom(ys, xs):
        dx = c0 + xs + 0.5 - cx
        dy = cy - (r0 + ys + 0.5)
        rr = np.hypot(dx, dy)
        th = np.mod(np.arctan2(dx, dy), 2 * np.pi)
        dth = np.abs((th[:, None] - rays[None, :] + np.pi) % (2 * np.pi) - np.pi)
        dist = np.where(dth < np.pi / 2, rr[:, None] * np.sin(dth), np.inf).min(1)
        return dx, dy, rr, th, dist

    ys, xs = np.nonzero(crop.min(2) > 200)
    dx, dy, rr, th, dist = geom(ys, xs)
    txt = (rr > 0.55 * R) & (rr < 0.95 * R) & (dist > 20)
    ys, xs, dx, dy, th = ys[txt], xs[txt], dx[txt], dy[txt], th[txt]
    sl = np.searchsorted(rays, th) % nr
    labels = {}
    for s in np.unique(sl):
        m = sl == s
        if m.sum() < 500:
            continue
        lo, hi = rays[s - 1], rays[s]
        mid = np.mod(lo + np.mod(hi - lo, 2 * np.pi) / 2, 2 * np.pi)
        smp = a[int(round(cy - 0.45 * R * np.cos(mid) - 0.5)), int(round(cx + 0.45 * R * np.sin(mid) - 0.5))]
        ci = int(np.argmin(((BLEND - smp) ** 2).sum(1)))
        by, bx = np.mgrid[ys[m].min() - 25:ys[m].max() + 26, xs[m].min() - 25:xs[m].max() + 26]
        by, bx = by.ravel(), bx.ravel()
        _, _, brr, bth, bdist = geom(by, bx)
        ok = (np.searchsorted(rays, bth) % nr == s) & (bdist > 20) & (brr < 0.97 * R)
        crop[by[ok], bx[ok]] = np.round(BLEND[ci]).astype(np.uint8)
        labels[ci] = ((dx[m].min() + dx[m].max()) / 2 / R, (dy[m].min() + dy[m].max()) / 2 / R)
        log(f"    % label in {SF[ci]} slice at ({labels[ci][0]:+.3f}, {labels[ci][1]:+.3f}) R; "
            f"{m.sum()} text px painted out ({ok.sum()} px box)", verbose)
    ext = ((c0 - cx) / R, (c1 - cx) / R, (cy - r1) / R, (cy - r0) / R)
    return crop, ext, labels


# ------------------------------------------------------------------ figure furniture
def mm_axes(fig, x, y, w, h):
    """Axes at (x, y) = upper-left corner in mm from the figure's top-left, size w x h mm."""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def fig_text(fig, x, y, s, **kw):
    fig.text(x / W, 1 - y / H, s, **kw)


def kfmt(v):
    return "0" if v == 0 else f"{int(v / 1000)}k"


def hist_axes(fig, x, y, w, h, xlim, ylim, xticks, xlabels, yticks, ylabels, facet, right=False):
    ax = mm_axes(fig, x, y, w, h)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    for s in ax.spines.values():
        s.set_visible(True); s.set_color(FRAME); s.set_linewidth(0.6)
    ax.tick_params(colors=TICKC, labelcolor="black", labelsize=FS_TICK, width=0.5, length=2.2, pad=1.5)
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels)
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels)
    ax.text(0.985 if right else 0.04, 0.965, facet, transform=ax.transAxes,
            ha="right" if right else "left", va="top", fontsize=FS_TITLE)
    return ax


def draw_bars(ax, counts):
    """Stacked bars from a Series indexed by (bin index, superfamily); 01_ on top, as published."""
    tab = counts.unstack(fill_value=0) if isinstance(counts.index, pd.MultiIndex) else counts
    tab = tab.reindex(columns=[s for s in SF if s in tab.columns], fill_value=0)
    centres = np.array(tab.index, float) * BINW
    bottom = np.zeros(len(tab))
    for sf in reversed(list(tab.columns)):           # bottom-up: 17_unknown first, 01_DNA/DTA last
        vals = tab[sf].to_numpy(float)
        ax.bar(centres, vals, bottom=bottom, width=BINW, color=FILL[sf], edgecolor="white",
               linewidth=0.25, align="center", zorder=2)
        bottom += vals


def draw_pie(ax, counts, verbose=False):
    """Pie of one panel: wedges clockwise from 12 o'clock in superfamily order, white borders,
    white % labels on slices > 5 %. Returns the label positions for reporting."""
    from matplotlib.patches import Wedge
    s = counts.reindex(SF).fillna(0)
    frac = (s / s.sum()).to_numpy()
    start = 90.0
    for sf, f in zip(SF, frac):
        if f <= 0:
            continue
        ax.add_patch(Wedge((0, 0), 1.0, start - f * 360, start, facecolor=FILL[sf], edgecolor="white",
                           linewidth=0.3, zorder=2))
        if 100 * f > PCT_MIN:
            mid = np.radians(start - f * 180)
            ax.text(PCT_R * np.cos(mid), PCT_R * np.sin(mid), f"{100 * f:.2f}%", ha="center", va="center",
                    color="white", fontsize=FS_PCT, fontweight="bold", zorder=3)
        start -= f * 360
    ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02); ax.set_aspect("equal"); ax.axis("off")


def chrome(fig):
    """Everything that is the same in both modes: pie titles, strips, axis titles, legend, letters."""
    from matplotlib.patches import Rectangle
    titles = {"A_orig": ("Original", "loci"), "A_sim": ("Simulated", "loci"),
              "B_orig": ("Original", "bases"), "B_sim": ("Simulated", "bases")}
    for key, xmm in PIE_X.items():
        lab, unit = titles[key]
        fig_text(fig, xmm + DIA / 2, PIE_Y - 1.2, f"{lab}\n({TOTALS[key]:,} {unit})", ha="center", va="bottom",
                 fontsize=FS_TITLE, linespacing=1.15)
    for x, text in ((5.0, "Total TE loci"), (77.0, "Total TE bases")):
        ax = mm_axes(fig, x, PIE_Y, 3.6, DIA)
        ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=STRIP_FILL, edgecolor=FRAME, lw=0.6))
        ax.text(0.5, 0.5, text, rotation=90, ha="center", va="center", fontsize=FS_TITLE, transform=ax.transAxes)
        ax.axis("off")
    fig_text(fig, CX["C_orig"] + CW + 1.0, CY0 + CH + 5.2, "Divergence (p-distance)", ha="center", va="top",
             fontsize=FS_AXIS)
    fig_text(fig, DX0 + DW / 2, CY0 + CH + 5.2, "Integrity", ha="center", va="top", fontsize=FS_AXIS)
    for i, (sf, col) in enumerate(zip(SF, PAL)):
        y = LEG_Y + i * LEG_ROW
        ax = mm_axes(fig, LEG_X, y, 2.6, 2.1)
        ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=FILL[sf], lw=0))
        ax.axis("off")
        fig_text(fig, LEG_X + 3.4, y + 1.05, f"{i + 1:02d}_{sf}", ha="left", va="center", fontsize=FS_LEG)
    for letter, x, y in (("A", 0.5, 2.0), ("B", 72.5, 2.0), ("C", 0.5, 45.5), ("D", 72.5, 45.5)):
        fig_text(fig, x, y, letter, ha="left", va="top", fontsize=FS_LETTER, fontweight="bold")


CT = np.arange(0, 0.51, 0.1)
DT = np.arange(0, 1.01, 0.1)
CY_T, DY_T = np.arange(0, 100001, 25000), np.arange(0, 80001, 20000)


# ------------------------------------------------------------------ the two builds
def build_vector(fig, t, verbose=False):
    for key, xmm in PIE_X.items():
        ax = mm_axes(fig, xmm, PIE_Y, DIA, DIA)
        draw_pie(ax, t[key], verbose)
    for key, x in CX.items():
        ax = hist_axes(fig, x, CY0, CW, CH, C_XLIM, C_YLIM, CT, [f"{v:.1f}" for v in CT],
                       CY_T, [kfmt(v) for v in CY_T] if key == "C_orig" else [],
                       "Original" if key == "C_orig" else "Simulated")
        draw_bars(ax, t[key])
        if key == "C_orig":
            ax.set_ylabel("TE loci", fontsize=FS_AXIS, labelpad=2)
    for key, y in (("D_orig", CY0), ("D_sim", CY0 + CH - DH)):
        ax = hist_axes(fig, DX0, y, DW, DH, D_XLIM, D_YLIM, DT,
                       [f"{v:.1f}" for v in DT] if key == "D_sim" else [], DY_T, [kfmt(v) for v in DY_T],
                       "Original" if key == "D_orig" else "Simulated", right=True)
        draw_bars(ax, t[key])
        ax.set_ylabel("TE loci", fontsize=FS_AXIS, labelpad=2)


def build_bitmap(fig, a, cal, verbose=False):
    get = load_pngs(a.pptx)
    for key, xmm in PIE_X.items():
        name = IMG[key]
        log(f"  pie {key} ({name})", verbose)
        arr, ext, labs = pie_clean(get(name), cal["pie"][name], verbose)
        pad = DIA / 2 * (ext[1] - 1)
        ax = mm_axes(fig, xmm - pad, PIE_Y - pad, DIA + 2 * pad, DIA + 2 * pad)
        ax.imshow(resample(arr, DIA + 2 * pad, DIA + 2 * pad), extent=ext, interpolation="none", zorder=0)
        ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal"); ax.axis("off")
        want = cal["pie_labels_in_png"][name]
        got = {SF[ci]: xy for ci, xy in labs.items()}
        if set(got) != set(want):
            sys.exit(f"ERROR: {name}: found % labels in {sorted(got)}, expected {sorted(want)}")
        for sf, (x, y) in got.items():  # same angle as the PNG label; PCT_R keeps 5.5-pt text on the disc
            f = PCT_R / np.hypot(x, y)
            ax.text(x * f, y * f, want[sf], ha="center", va="center", color="white", fontsize=FS_PCT,
                    fontweight="bold")
    for key, x in CX.items():
        arr, ext = hist_crop(get(IMG[key]), cal["hist"][IMG[key]], C_XLIM, C_YLIM)
        ax = hist_axes(fig, x, CY0, CW, CH, C_XLIM, C_YLIM, CT, [f"{v:.1f}" for v in CT],
                       CY_T, [kfmt(v) for v in CY_T] if key == "C_orig" else [],
                       "Original" if key == "C_orig" else "Simulated")
        ax.imshow(resample(arr, CW, CH), extent=ext, aspect="auto", interpolation="none", zorder=0)
        if key == "C_orig":
            ax.set_ylabel("TE loci", fontsize=FS_AXIS, labelpad=2)
        log(f"  C {key} crop {arr.shape[1]}x{arr.shape[0]} px", verbose)
    for key, y in (("D_orig", CY0), ("D_sim", CY0 + CH - DH)):
        arr, ext = hist_crop(get(IMG[key]), cal["hist"][IMG[key]], D_XLIM, D_YLIM)
        ax = hist_axes(fig, DX0, y, DW, DH, D_XLIM, D_YLIM, DT,
                       [f"{v:.1f}" for v in DT] if key == "D_sim" else [], DY_T, [kfmt(v) for v in DY_T],
                       "Original" if key == "D_orig" else "Simulated", right=True)
        ax.imshow(resample(arr, DW, DH), extent=ext, aspect="auto", interpolation="none", zorder=0)
        ax.set_ylabel("TE loci", fontsize=FS_AXIS, labelpad=2)
        log(f"  D {key} crop {arr.shape[1]}x{arr.shape[0]} px", verbose)


# ------------------------------------------------------------------ source data
def long_hist(t, key, dataset):
    s = t[key]
    df = s.rename("value").reset_index()
    df.columns = ["k", "superfamily", "value"]
    df.insert(0, "dataset", dataset)
    df["bin_left"] = (df.k * BINW - BINW / 2).round(2)
    df["bin_right"] = (df.k * BINW + BINW / 2).round(2)
    return df[["dataset", "superfamily", "bin_left", "bin_right", "value"]]


def source_data(t, seed):
    sd = SourceData("fig5")
    for panel, keys, unit, desc in (
        ("Fig. 5A", ("A_orig", "A_sim"), "TE loci",
         "Number of annotated TE loci per TE superfamily. Original: RepeatMasker copy number per family, summed by "
         "superfamily (total 1,181,685 loci). Simulated: features in the TEgenomeSimulator replica annotation "
         "(total 1,251,338 loci)."),
        ("Fig. 5B", ("B_orig", "B_sim"), "bases (bp)",
         "Bases (bp) covered by annotated TE loci per TE superfamily, summed over loci without merging overlaps. "
         "Original total 981,930,851 bp; simulated total 793,915,337 bp.")):
        df = pd.concat([t[k].rename("value").rename_axis("superfamily").reset_index().assign(dataset=ds)
                        for k, ds in zip(keys, ("original", "simulated"))])
        sd.add(panel, desc, df[["dataset", "superfamily", "value"]].reset_index(drop=True))
    dfc = pd.concat([long_hist(t, "C_orig", "original"), long_hist(t, "C_sim", "simulated")], ignore_index=True)
    sd.add("Fig. 5C", "Number of TE loci per divergence bin and TE superfamily. Divergence (p-distance) = 1 - identity "
           "to the TE library consensus (RepeatMasker substitutions + indels for the original genome; the per-copy "
           "identity assigned by the simulator for the replica). Bins are (bin_left, bin_right], width 0.02. "
           "Original n = 1,181,685 loci; simulated n = 1,251,338 loci.", dfc)
    dfd = pd.concat([long_hist(t, "D_orig", "original"), long_hist(t, "D_sim", "simulated")], ignore_index=True)
    tot = dfd.groupby("dataset").value.sum()
    sd.add("Fig. 5D", "Number of TE loci per integrity bin and TE superfamily. Integrity = element length / TE library "
           "consensus length; values > 1 are not plotted. Bins are (bin_left, bin_right], width 0.02. Original: "
           f"length ratio of every RepeatMasker hit, n = {tot['original']:,} of 1,181,685 loci. Simulated: integrity "
           f"drawn per locus from the original per-family integrity lists (random seed {seed}), n = {tot['simulated']:,} "
           f"of 1,251,338 loci, with the bin centred on 1.0 set to {ONE_BIN_SIM:,} loci as in the published panel.", dfd)
    sd.save()


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default=str(ROOT / "fig5" / "fig5.pdf"), help="output PDF (default v2/fig5/fig5.pdf)")
    ap.add_argument("--mode", choices=("vector", "bitmap"), default="vector",
                    help="vector: redraw every mark from the data (default); bitmap: published panels as images")
    ap.add_argument("--scan", default=str(SCAN), help="TEgenomeSimulator mode-2 scan dir of the original genome")
    ap.add_argument("--gff", default=str(GFF), help="simulated replica annotation GFF")
    ap.add_argument("--seed", type=int, default=1, help="random seed for the simulated integrity draw (default 1)")
    ap.add_argument("--pptx", default=str(PPTX), help="--mode bitmap: source deck (read-only; slide 14 PNGs)")
    ap.add_argument("--inputs", default=str(ROOT / "build" / "fig5" / "inputs"), help="calibrations + published tables")
    ap.add_argument("--png-dpi", type=int, default=200, help="preview PNG resolution (0 = none; default 200)")
    ap.add_argument("--no-source-data", action="store_true", help="skip writing Source Data blocks")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-panel progress and sanity checks")
    a = ap.parse_args()
    v = a.verbose
    cal = json.loads((Path(a.inputs) / "calibrations.json").read_text())

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ncomms_style.use()
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    t = None
    if a.mode == "vector":
        log(f"make_fig5: exact data from {a.scan} and {a.gff}")
        t = exact_tables(a.scan, a.gff, a.seed, v)
        pct_check(t, cal, v)
        build_vector(fig, t, v)
    else:
        log(f"make_fig5: published panels from {a.pptx}")
        build_bitmap(fig, a, cal, v)
    chrome(fig)

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=600)
    log(f"make_fig5: wrote {out} ({a.mode} mode)")
    if a.png_dpi:
        import fitz
        png = ROOT / "build" / "fig5" / f"fig5_{a.png_dpi}dpi.png"
        fitz.open(out)[0].get_pixmap(dpi=a.png_dpi).save(str(png))
        log(f"make_fig5: preview {png}")
    if not a.no_source_data:
        if t is None:
            t = exact_tables(a.scan, a.gff, a.seed, v)
        source_data(t, a.seed)
        log("make_fig5: Source Data -> build/source_data/fig5")
    log("make_fig5: done")


if __name__ == "__main__":
    main()
