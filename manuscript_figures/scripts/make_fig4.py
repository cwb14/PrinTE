#!/usr/bin/env python3
"""Build final Fig. 4 (LTR-RT landscape across Pucciniomycotina) for Nature Communications.

Rebuilds the Inkscape composite Figure_4_v2 (original plotting scripts lost) from recovered data,
re-laid-out at 180 mm width with Arial text >= 5 pt, and writes its Source Data blocks.

  A  IQ-TREE 2.4.0 ML tree (69 BUSCO loci, UFBoot2 x1000 shown); tip dots = Table S2 assembly size
  B  PCA (centred, unit variance) of 4 LTR metrics; n = 22 genomes with LTR content > 1 %
  C  Pearson r (two-sided) + OLS fit and 95 % CI of the mean vs assembly size; same n = 22
  D  divergence (1 - LTR identity) of EDTA structurally intact LTR-RTs; bins of 0.001; top-5 families

Inputs (copied read-only into --inputs; see build/fig4/memo_fig4.md for provenance):
  supp_table2.xlsx, Pucciniomycotina_dedup_busco_partitions.txt.treefile, tip_order.txt,
  3 x *.TEanno.gff3.intactLTR (Name, classification, ltr_identity, motif, TSD).
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ncomms_style  # noqa: E402
from srcdata import SourceData  # noqa: E402

ROOT = HERE.parent
W = 180.0  # figure width, mm
H = 170.0  # figure height, mm
PT = 25.4 / 72  # mm per point

# ---- font sizes (pt, final size)
FS_LETTER, FS_AXIS, FS_TICK, FS_TIP, FS_SUPPORT = 8, 6.5, 5.5, 5.5, 5
FS_STRIP, FS_LEG, FS_LEGT, FS_CLASS, FS_SPECIES = 6, 5.5, 6, 6.5, 6.5
GREY_TICK, GREY_AXIS, GREY_SUPPORT = "#4D4D4D", "#333333", "#666666"

METRICS = [  # plotted name (as in original), Table S2 column
    ("LTR_percent", "LTR percentage (%)"),
    ("Intact.count", "Intact count"),
    ("Solo.count", "Solo count"),
    ("Solo.Intact.ratio", "Solo:Intact ratio"),
]
METRIC_COL = dict(zip([m for m, _ in METRICS], ["#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"]))  # ggplot hue
EXPECTED_C = {"LTR_percent": "LTR_percent (r = 0.81, p = 5.0e-06)",       # titles printed on Figure_4_v2
              "Intact.count": "Intact.count (r = 0.76, p = 3.9e-05)",
              "Solo.count": "Solo.count (r = 0.98, p = 2.5e-16)",
              "Solo.Intact.ratio": "Solo.Intact.ratio (r = 0.65, p = 1.2e-03)"}
EXPECTED_PCA = (72.99, 19.17)

FAM_GROUPS = ["Top 1", "Top 2", "Top 3", "Top 4", "Top 5", "Other"]
FAM_COL = ["#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#7F7F7F"]  # as in Figure_4_v2
D_GENOMES = [  # (tree tip, intactLTR file, title (mathtext), expected n plotted, expected n total)
    ("Crori1", "Crori1_AssemblyScaffolds.fasta.mod.EDTA.TEanno.gff3.intactLTR",
     r"$\it{Cronartium\ ribicola}$", 1498, 1499),
    ("PuccoNC29_1", "PuccoNC29_1_AssemblyScaffolds.dedup.fasta.mod.EDTA.TEanno.gff3.intactLTR",
     r"$\it{Puccinia\ coronata}$ f. sp. $\it{avenae}$", 1118, 1118),
    ("myrtle_rust", "GCA_902702905.1_Austropuccinia_psidii_genomic.fna.mod.EDTA.TEanno.gff3.intactLTR",
     r"$\it{Austropuccinia\ psidii}$", 15093, 15228),
]
D_SPECIES = {"Crori1": "Cronartium ribicola Cypress_4 v1.0 (Crori1)",
             "PuccoNC29_1": "Puccinia coronata f. sp. avenae 12NC29 (PuccoNC29_1)",
             "myrtle_rust": "Austropuccinia psidii GCA_902702905.1 (myrtle_rust)"}
B_LABELS = {"Crori1": "C. ribicola", "PuccoNC29_1": "P. coronata", "myrtle_rust": "A. psidii"}
CLASS_STRIPS = [("Microbotryomycetes", "Rhoto_IFO1236_1", "Hyabl1", "#D9F2D0"),  # tip ranges/colours
                ("Pucciniomycetes", "Pgt_Ug99_C1", "Sepsp1", "#F2CFEE")]       # as in Figure_4_v2
BIN = 0.001


def log(msg, verbose, always=False):
    if verbose or always:
        print(msg, file=sys.stderr)


# ------------------------------------------------------------------ data
def read_s2(path):
    s2 = pd.read_excel(path, header=1)
    s2 = s2[s2["ID"].notna()].copy()
    need = ["Genome", "Class", "ID", "Assembly size"] + [c for _, c in METRICS]
    miss = [c for c in need if c not in s2.columns]
    if miss:
        sys.exit(f"ERROR: {path} lacks columns {miss}")
    s2["size_Mb"] = s2["Assembly size"].astype(float) / 1e6
    return s2.set_index("ID", drop=False)


def parse_newick(s):
    """Minimal Newick parser -> nested dicts {name, label, bl, children}."""
    s = s.strip().rstrip(";")
    pos = 0
    tok_re, bl_re = re.compile(r"[^:,();]*"), re.compile(r":([-0-9.eE+]+)")

    def node():
        nonlocal pos
        n = {"children": [], "name": None, "label": None, "bl": 0.0}
        if s[pos] == "(":
            pos += 1
            while True:
                n["children"].append(node())
                if s[pos] == ",":
                    pos += 1
                elif s[pos] == ")":
                    pos += 1
                    break
                else:
                    raise ValueError(f"bad newick at {pos}")
        m = tok_re.match(s, pos)
        pos = m.end()
        if n["children"]:
            n["label"] = m.group(0) or None
        else:
            n["name"] = m.group(0)
        m = bl_re.match(s, pos)
        if m:
            n["bl"], pos = float(m.group(1)), m.end()
        return n

    root = node()
    if pos != len(s):
        raise ValueError("trailing characters in newick")
    return root


def layout_tree(root, order):
    """x = root-to-node distance; tip y = row in `order`; internal y = mean of children (ggtree)."""
    row = {t: i for i, t in enumerate(order)}
    nodes = []

    def walk(n, x):
        n["x"] = x + (n["bl"] if n is not root else 0.0)
        if not n["children"]:
            if n["name"] not in row:
                sys.exit(f"ERROR: tip {n['name']} missing from tip_order.txt")
            n["y"], n["rows"] = row[n["name"]], [row[n["name"]]]
        else:
            for c in n["children"]:
                walk(c, n["x"])
            n["y"] = float(np.mean([c["y"] for c in n["children"]]))
            n["rows"] = sorted(r for c in n["children"] for r in c["rows"])
            if n["rows"] != list(range(n["rows"][0], n["rows"][-1] + 1)):
                sys.exit("ERROR: tip_order.txt is not a valid drawing order for this topology")
        nodes.append(n)

    walk(root, 0.0)
    ntip = sum(1 for n in nodes if not n["children"])
    if ntip != len(order):
        sys.exit(f"ERROR: tree has {ntip} tips, tip_order.txt {len(order)}")
    return nodes


def pearson_ols(x, y):
    from scipy import stats
    import statsmodels.api as sm
    res = stats.pearsonr(x, y)
    ci = res.confidence_interval(0.95)
    n = len(x)
    r = res.statistic
    fit = sm.OLS(y, sm.add_constant(x)).fit()
    sci = fit.conf_int(0.05)[1]
    return dict(n=n, r=r, r_CI95_low=ci.low, r_CI95_high=ci.high, t=r * np.sqrt(n - 2) / np.sqrt(1 - r * r),
                df=n - 2, P_two_sided=res.pvalue, P_Bonferroni_4=min(1.0, 4 * res.pvalue), R2=fit.rsquared,
                slope_per_Mb=fit.params[1], slope_CI95_low=sci[0], slope_CI95_high=sci[1],
                intercept=fit.params[0]), fit


def pca_scaled(X):
    """prcomp(X, center=TRUE, scale.=TRUE); sign: largest-|loading| of each PC positive (matches R here)."""
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    _, s, vt = np.linalg.svd(Z, full_matrices=False)
    load = vt.T
    sign = np.sign(load[np.abs(load).argmax(0), np.arange(load.shape[1])])
    load = load * sign
    var = s ** 2 / (len(X) - 1)
    return Z @ load, load, 100 * var / var.sum()


def nrd0(x):
    """R bw.nrd0 (Silverman's rule of thumb)."""
    sd = np.std(x, ddof=1)
    lo = min(sd, np.subtract(*np.percentile(x, [75, 25])) / 1.34) or sd
    return 0.9 * lo * len(x) ** -0.2


def kde_reflect0(x, grid, bw):
    """Gaussian KDE reflected at 0 (divergence >= 0)."""
    k = lambda g: np.exp(-0.5 * ((g[:, None] - x[None, :]) / bw) ** 2).sum(1)
    return (k(grid) + k(-grid)) / (len(x) * bw * np.sqrt(2 * np.pi))


def lab_gradient(t):
    """ggplot2 scale_colour_gradient(low='blue', high='red'): interpolation in CIE Lab (D65)."""
    M = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    wp = M @ np.ones(3)
    d = 6 / 29
    lin = lambda c: np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    gam = lambda c: np.where(c <= 0.0031308, 12.92 * c, 1.055 * np.clip(c, 0, None) ** (1 / 2.4) - 0.055)
    f = lambda v: np.where(v > d ** 3, np.cbrt(v), v / (3 * d * d) + 4 / 29)
    finv = lambda v: np.where(v > d, v ** 3, 3 * d * d * (v - 4 / 29))

    def to_lab(rgb):
        fx, fy, fz = f(M @ lin(np.asarray(rgb, float)) / wp)
        return np.array([116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)])

    lo, hi = to_lab([0, 0, 1]), to_lab([1, 0, 0])
    out = []
    for tt in np.atleast_1d(t):
        L, a, b = lo + (hi - lo) * tt
        fy = (L + 16) / 116
        xyz = np.array([finv(fy + a / 500), finv(fy), finv(fy - b / 200)]) * wp
        out.append(np.clip(gam(np.linalg.solve(M, xyz)), 0, 1))
    return np.array(out)


def text_w_mm(s, size, style="normal", weight="normal"):
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath
    p = TextPath((0, 0), s, size=size, prop=FontProperties(family="Arial", style=style, weight=weight))
    return p.get_extents().width * PT


# ------------------------------------------------------------------ figure helpers
def ax_mm(fig, x, y, w, h):
    """Axes from mm box (x, y measured from the top-left corner)."""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def F(x, y):
    """mm (from top-left) -> figure fraction."""
    return x / W, 1 - y / H


def data_to_mm(fig, ax, x, y):
    fx, fy = fig.transFigure.inverted().transform(ax.transData.transform((x, y)))
    return fx * W, (1 - fy) * H


def ftext(fig, x, y, s, **kw):
    fx, fy = F(x, y)
    return fig.text(fx, fy, s, **kw)


def fline(fig, xs, ys, **kw):
    from matplotlib.lines import Line2D
    fx, fy = zip(*[F(x, y) for x, y in zip(xs, ys)])
    ln = Line2D(fx, fy, transform=fig.transFigure, **kw)
    fig.add_artist(ln)
    return ln


def letter(fig, x, y, s):
    ftext(fig, x, y, s, fontsize=FS_LETTER, fontweight="bold", ha="left", va="top")


def style_ticks(ax, box=True):
    ax.tick_params(labelsize=FS_TICK, labelcolor=GREY_TICK, color=GREY_AXIS, length=1.8, width=0.5, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(box)
        sp.set_color(GREY_AXIS)
        sp.set_linewidth(0.5)


# ------------------------------------------------------------------ panels
def panel_A(fig, nodes, order, s2, v):
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Rectangle
    sp = 2.3                      # tip spacing, mm (6.5 pt)
    top, x0, w = 4.5, 3.0, 42.0   # tree axes box, mm
    ax = ax_mm(fig, x0, top, w, sp * len(order))
    tips = {n["name"]: n for n in nodes if not n["children"]}
    xmax = max(n["x"] for n in tips.values())
    ax.set_xlim(0, xmax)
    ax.set_ylim(len(order) - 0.5, -0.5)
    segs = []
    for n in nodes:
        for c in n["children"]:
            segs.append([(n["x"], c["y"]), (c["x"], c["y"])])                  # horizontal
        if n["children"]:
            ys = [c["y"] for c in n["children"]]
            segs.append([(n["x"], min(ys)), (n["x"], max(ys))])                # vertical
    ax.add_collection(LineCollection(segs, colors="black", linewidths=0.6, capstyle="projecting", clip_on=False))
    # UFBoot support (second field of SH-aLRT/UFBoot) left of and above each non-root internal node;
    # thin white halo (outlined copy underneath) keeps digits legible where a vertical branch passes behind
    from matplotlib import patheffects
    halo = patheffects.withStroke(linewidth=1.0, foreground="white")
    from matplotlib.transforms import Bbox
    renderer = fig.canvas.get_renderer()

    def ink(t):   # digit ink box (cap height ~0.61 of the layout box, which includes descent + leading)
        b = t.get_window_extent(renderer)
        return Bbox.from_extents(b.x0, b.y0 + 0.18 * b.height, b.x1, b.y0 + 0.79 * b.height)

    placed, flipped = [], []
    for n in sorted([n for n in nodes if n["children"] and n["label"]], key=lambda n: -n["x"]):
        uf = n["label"].split("/")[-1]
        kw = dict(fontsize=FS_SUPPORT, ha="right", clip_on=False)
        t = ax.text(n["x"] - 0.03, n["y"] - 0.16, uf, va="bottom", **kw)         # above-left (as original)
        bb = ink(t)
        if any(bb.overlaps(o) for o in placed):                                    # ink clash -> try below-left
            t2 = ax.text(n["x"] - 0.03, n["y"] + 0.16, uf, va="top", **kw)
            if any(ink(t2).overlaps(o) for o in placed):
                t2.remove()
            else:
                t.remove()
                t, bb = t2, ink(t2)
                flipped.append(uf)
        t.set(color=GREY_SUPPORT, zorder=5)                                        # real (editable) text
        ax.text(*t.get_position(), uf, va=t.get_va(), color="white", zorder=4, path_effects=[halo], **kw)
        placed.append(bb)
    nsup = len(placed)
    # tip dots coloured by assembly size (ggplot blue->red Lab gradient over the 45 tips)
    size = s2.loc[order, "size_Mb"].values
    vmin, vmax = size.min(), size.max()
    cols = lab_gradient((size - vmin) / (vmax - vmin))
    label_x = x0 + w + 1.6   # mm
    for i, t in enumerate(order):
        tx = tips[t]["x"]
        ax.plot([tx], [i], "o", ms=2.3, mfc=cols[i], mec=cols[i], mew=0.3, clip_on=False, zorder=3)
        dx, dy = data_to_mm(fig, ax, tx, i)
        fline(fig, [dx + 0.7, label_x - 0.35], [dy, dy], color="black", lw=0.4, ls=(0, (0.6, 1.8)),
              dash_capstyle="round")                                   # aligned-label leader (ggtree style)
        ftext(fig, label_x, dy, t, fontsize=FS_TIP, ha="left", va="center_baseline")
    # tree x axis (substitutions per site), as in the original: ticks only
    for s_ in ("left", "right", "top"):
        ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_bounds(0, xmax)
    ax.set_yticks([])
    ax.set_xticks([0, 0.5, 1.0, 1.5])
    ax.set_xticklabels(["0.0", "0.5", "1.0", "1.5"])
    style_ticks(ax, box=False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_position(("outward", 1.5))
    ax.patch.set_visible(False)
    # colour bar "Assembly Size (Mb)" (vector: 256 stacked rectangles)
    cb_x, cb_y, cb_w, cb_h = 3.0, 9.0, 2.6, 15.0
    ftext(fig, cb_x, cb_y - 1.2, "Assembly Size (Mb)", fontsize=FS_LEGT, ha="left", va="bottom")
    cax = ax_mm(fig, cb_x, cb_y, cb_w, cb_h)
    cax.set_xlim(0, 1)
    cax.set_ylim(vmin, vmax)
    edges = np.linspace(vmin, vmax, 257)
    gcols = lab_gradient((edges[:-1] + np.diff(edges) / 2 - vmin) / (vmax - vmin))
    for e0, c in zip(edges[:-1], gcols):   # each slice runs to the top: later slices cover seams
        cax.add_patch(Rectangle((0, e0), 1, vmax - e0, fc=c, ec="none", lw=0))
    brk = [250, 500, 750, 1000]
    for b in brk:
        cax.plot([0, 0.2], [b, b], color="white", lw=0.4)
        cax.plot([0.8, 1], [b, b], color="white", lw=0.4)
    cax.set_xticks([])
    cax.yaxis.tick_right()
    cax.set_yticks(brk)
    cax.set_yticklabels([str(b) for b in brk])
    cax.tick_params(axis="y", length=0, labelsize=FS_LEG, pad=1.2)
    for sp_ in cax.spines.values():
        sp_.set_visible(False)
    # geometry for class strips / connectors
    widths = {t: text_w_mm(t, FS_TIP) for t in order}
    first_conn = min(order.index(t) for t, *_ in D_GENOMES)
    below = max(widths[t] for t in order[first_conn:])
    chan0 = label_x + below + 1.6
    strip_x = max(label_x + max(widths.values()) + 1.2, chan0 + 3.0 + 1.4)
    for name, t0, t1, col in CLASS_STRIPS:
        y0 = data_to_mm(fig, ax, 0, order.index(t0) - 0.45)[1]
        y1 = data_to_mm(fig, ax, 0, order.index(t1) + 0.45)[1]
        fig.add_artist(Rectangle(F(strip_x, y1), 3.4 / W, (y1 - y0) / H, transform=fig.transFigure,
                                 fc=col, ec="none", lw=0))
        ftext(fig, strip_x + 1.7, (y0 + y1) / 2, name, fontsize=FS_CLASS, style="italic", rotation=90,
              ha="center", va="center")
    geo = dict(label_x=label_x, chan0=chan0, strip_right=strip_x + 3.4, widths=widths,
               tip_mm={t: data_to_mm(fig, ax, tips[t]["x"], i)[1] for i, t in enumerate(order)},
               tree_bottom=top + sp * len(order))
    log(f"A: {len(order)} tips, {nsup} UFBoot labels ({len(flipped)} moved below their branch to avoid a clash), size range {vmin:.2f}-{vmax:.2f} Mb, "
        f"labels at {label_x:.1f} mm, channel {chan0:.1f} mm, strips {strip_x:.1f} mm", v)
    return geo


def panel_C(fig, sub, x_left, v):
    import statsmodels.api as sm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FixedLocator, FuncFormatter
    stats_rows, fits = [], {}
    gap, right = 8.0, 2.0
    w = (W - right - x_left - gap) / 2
    hs, hp, top = 3.4, 38.5, 4.5
    x = sub["size_Mb"].values
    xg = np.linspace(x.min(), x.max(), 80)                   # ggplot stat_smooth grid
    rx = np.ptp(x)
    xl = (x.min() - 0.05 * rx, x.max() + 0.05 * rx)
    fmt = {"LTR_percent": (lambda t, _: f"{t:.0f}%", [0, 25, 50, 75, 100]),
           "Intact.count": (lambda t, _: f"{t:.0f}", [0, 2500, 5000, 7500]),
           "Solo.count": (lambda t, _: f"{t / 1000:.0f}k", [0, 50000, 100000, 150000, 200000]),
           "Solo.Intact.ratio": (lambda t, _: f"{t:.0f}", [0, 10, 20, 30])}
    fixed = {"LTR_percent": (0, 100), "Intact.count": (0, 9000)}  # limits used in Figure_4_v2
    for k, (m, col) in enumerate(METRICS):
        r_, c_ = divmod(k, 2)
        ax_x = x_left + c_ * (w + gap)
        ax_y = top + r_ * (hs + hp + 4.2) + hs
        ax = ax_mm(fig, ax_x, ax_y, w, hp)
        y = sub[col].astype(float).values
        st, fit = pearson_ols(x, y)
        pr = fit.get_prediction(sm.add_constant(xg)).summary_frame(alpha=0.05)
        lo, hi, mu = pr["mean_ci_lower"].values, pr["mean_ci_upper"].values, pr["mean"].values
        ax.fill_between(xg, lo, hi, color="#ADD8E6", alpha=0.4, lw=0, zorder=1)
        ax.plot(xg, mu, color="black", lw=1.0, ls=(0, (1, 3)), dash_capstyle="round", zorder=2)
        ax.scatter(x, y, s=3.3 ** 2, c=METRIC_COL[m], alpha=0.7, lw=0, zorder=3)
        if m in fixed:
            a, b = fixed[m]
        else:
            a, b = min(y.min(), lo.min()), max(y.max(), hi.max())
        ax.set_ylim(a - 0.05 * (b - a), b + 0.05 * (b - a))
        ax.set_xlim(*xl)
        f_, ticks = fmt[m]
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        ax.yaxis.set_major_formatter(FuncFormatter(f_))
        ax.xaxis.set_major_locator(FixedLocator([0, 250, 500, 750, 1000]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda t, _: f"{t:,.0f}"))
        if r_ == 0:
            ax.tick_params(labelbottom=False)
        style_ticks(ax, box=True)
        title = f"{m} (r = {st['r']:.2f}, p = {st['P_two_sided']:.1e})"
        if title != EXPECTED_C[m]:
            sys.exit(f"ERROR: panel C title drift: {title!r} != {EXPECTED_C[m]!r}")
        fig.add_artist(Rectangle(F(ax_x, ax_y), w / W, hs / H, transform=fig.transFigure, fc="#D9D9D9",
                                 ec=GREY_AXIS, lw=0.5))
        ftext(fig, ax_x + w / 2, ax_y - hs / 2, title, fontsize=FS_STRIP, ha="center", va="center_baseline")
        stats_rows.append(dict(metric=m, table_S2_column=col, **st, panel_title=title))
        fits[m] = (xg, lo, hi, mu)
    bottom = top + 2 * (hs + hp) + 4.2
    ftext(fig, x_left + w + gap / 2, bottom + 4.2, "Genome size (Mb)", fontsize=FS_AXIS, ha="center", va="top")
    # legend row "Metric"
    ly = bottom + 9.3
    xx = x_left + 3
    ftext(fig, xx, ly, "Metric", fontsize=FS_LEGT, ha="left", va="center_baseline")
    xx += text_w_mm("Metric", FS_LEGT) + 2.5
    for m, _ in METRICS:
        fx, fy = F(xx, ly - 0.55)
        fig.add_artist(Line2D([fx], [fy], transform=fig.transFigure, marker="o", ms=3.3, mfc=METRIC_COL[m],
                              mec="none", alpha=0.7))
        ftext(fig, xx + 1.6, ly, m, fontsize=FS_LEG, ha="left", va="center_baseline")
        xx += 1.6 + text_w_mm(m, FS_LEG) + 3.2
    log(f"C: facets {w:.1f} x {hp:.1f} mm, n = {len(x)}", v)
    return pd.DataFrame(stats_rows), ly


def panel_B(fig, sub, top, v):
    from matplotlib.patches import FancyArrowPatch
    X = sub[[c for _, c in METRICS]].astype(float).values
    scores, load, pve = pca_scaled(X)
    if (round(pve[0], 2), round(pve[1], 2)) != EXPECTED_PCA:
        sys.exit(f"ERROR: PCA variance drift {pve[:2]} != {EXPECTED_PCA}")
    x0, w, h = 9.5, 47.5, 38.0
    ax = ax_mm(fig, x0, top + 4.5, w, h)
    arrows = 3.0 * load[:, :2]                 # loadings x 3, as drawn in Figure_4_v2
    ids = list(sub.index)
    lab = {i: (scores[ids.index(i), 0], scores[ids.index(i), 1]) for i in B_LABELS}
    allx = np.r_[scores[:, 0], arrows[:, 0], 0]
    ally = np.r_[scores[:, 1], arrows[:, 1], 0, [lab[i][1] + 0.3 for i in lab]]
    ax.set_xlim(allx.min() - 0.05 * np.ptp(allx), allx.max() + 0.05 * np.ptp(allx))
    ax.set_ylim(ally.min() - 0.05 * np.ptp(ally), ally.max() + 0.05 * np.ptp(ally))
    ax.scatter(scores[:, 0], scores[:, 1], s=2.7 ** 2, c="black", lw=0, zorder=3)
    for (m, _), (ax_, ay_) in zip(METRICS, arrows):
        ax.add_patch(FancyArrowPatch((0, 0), (ax_, ay_), arrowstyle="->,head_length=2.4,head_width=1.3",
                                     mutation_scale=1, color="blue", lw=0.75, zorder=2, shrinkA=0, shrinkB=0))
    # loading labels next to the arrow tips; LTR_percent needs a leader (its tip sits by the P. pachyrhizi points)
    pos = {"Solo.Intact.ratio": ((arrows[3, 0] + 0.12, arrows[3, 1]), "left"),
           "Solo.count": ((arrows[2, 0] + 0.14, arrows[2, 1] + 0.12), "left"),
           "Intact.count": ((arrows[1, 0] + 0.14, arrows[1, 1]), "left"),
           "LTR_percent": ((2.95, -0.52), "left")}
    for m, ((tx, ty), ha) in pos.items():
        ax.text(tx, ty, m, color="blue", fontsize=FS_LEG, ha=ha, va="center", zorder=4)
    lx, ly = arrows[0]
    ax.plot([lx + 0.08, 2.9], [ly - 0.02, -0.52], color="blue", lw=0.35, zorder=2)
    for i, name in B_LABELS.items():
        px, py = lab[i]
        ax.text(px - 0.14, py + 0.15, name, color="red", fontsize=FS_SPECIES, style="italic", ha="right",
                va="center", zorder=4)
    ax.set_xticks([-2, 0, 2, 4])
    ax.set_yticks([-1, 0, 1, 2, 3])
    style_ticks(ax, box=True)
    ax.set_xlabel(f"PC1 ({pve[0]:.2f}%)", fontsize=FS_AXIS, labelpad=1.5)
    ax.set_ylabel(f"PC2 ({pve[1]:.2f}%)", fontsize=FS_AXIS, labelpad=1.5)
    sc = pd.DataFrame({"ID": ids, "genome": sub["Genome"].values, "PC1": scores[:, 0], "PC2": scores[:, 1],
                       "labelled_as": [B_LABELS.get(i, "") for i in ids]})
    ld = pd.DataFrame({"row": [m for m, _ in METRICS], "PC1": load[:, 0], "PC2": load[:, 1],
                       "arrow_end_PC1 (3 x loading)": arrows[:, 0], "arrow_end_PC2 (3 x loading)": arrows[:, 1]})
    ld = pd.concat([ld, pd.DataFrame([{"row": "variance explained (%)", "PC1": pve[0], "PC2": pve[1]}])],
                   ignore_index=True)
    log(f"B: PC1 {pve[0]:.2f}% PC2 {pve[1]:.2f}%", v)
    return sc, ld, pve


def panel_D(fig, inputs, top, x_left, verbose):
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
    gap, ah = 9.0, 29.5
    aw = (W - 1.5 - x_left - 2 * gap) / 3
    rows, geo = [], {}
    grid = np.linspace(0, 0.101, 600)
    for k, (tip, fn, title, n_exp, n_tot) in enumerate(D_GENOMES):
        d = pd.read_csv(inputs / fn, sep="\t", header=None,
                        names=["family", "classification", "ltr_identity", "motif", "tsd"])
        if len(d) != n_tot:
            sys.exit(f"ERROR: {fn}: {len(d)} elements, expected {n_tot}")
        d["divergence"] = 1 - d["ltr_identity"].astype(float)
        d["bin"] = np.floor(d["divergence"] / BIN).astype(int)   # float64 floor, as in Figure_4_v2
        d["plotted"] = d["bin"] <= 100                           # bins [0, 0.101): 0.00-0.10 axis
        top5 = d["family"].value_counts().head(5).index.tolist()
        d["top5_group"] = d["family"].map({f: f"Top {i + 1}" for i, f in enumerate(top5)}).fillna("Other")
        p = d[d["plotted"]]
        if len(p) != n_exp:
            sys.exit(f"ERROR: {tip}: {len(p)} plotted elements, expected {n_exp}")
        ax_x = x_left + k * (aw + gap)
        ax = ax_mm(fig, ax_x, top + 6.5, aw, ah)
        edges = np.arange(0, 102) * BIN
        base = np.zeros(101)
        for g, col in zip(FAM_GROUPS, FAM_COL):
            cnt = np.bincount(p.loc[p["top5_group"] == g, "bin"], minlength=101)[:101]
            ax.fill_between(edges, np.r_[base, base[-1]], np.r_[base + cnt, (base + cnt)[-1]], step="post",
                            color=col, lw=0, zorder=1)
            base = base + cnt
        bw = nrd0(p["divergence"].values)
        dens = kde_reflect0(p["divergence"].values, grid, bw) * len(p) * BIN
        ax.plot(grid, dens, color="black", lw=1.0, zorder=3, solid_capstyle="round")
        ax.set_xlim(-0.005, 0.106)
        ymax = max(base.max(), dens.max())
        ax.set_ylim(-0.04 * ymax, 1.05 * ymax)
        ax.xaxis.set_major_locator(FixedLocator([0, 0.02, 0.04, 0.06, 0.08, 0.10]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda t, _: f"{t:.2f}"))
        ax.yaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 5, 10], integer=True))
        style_ticks(ax, box=False)
        ax.tick_params(length=0, pad=1.2)
        ax.set_xlabel("LTR divergence", fontsize=FS_AXIS, labelpad=2)
        ax.set_ylabel("Count", fontsize=FS_AXIS, labelpad=1.5)
        ax.patch.set_visible(False)
        ftext(fig, ax_x, top + 4.3, title, fontsize=FS_SPECIES, ha="left", va="baseline")
        geo[tip] = dict(title_x=ax_x, title_top=top + 4.3 - FS_SPECIES * PT * 0.75)
        out = d[["family", "classification", "ltr_identity", "divergence", "top5_group", "plotted"]].copy()
        out["divergence"] = out["divergence"].round(4)         # identity has 4 decimals; drop float noise
        out.insert(5, "bin_start", np.where(d["plotted"], (d["bin"] * BIN).round(3), np.nan))
        out.insert(0, "genome", D_SPECIES[tip])
        rows.append(out)
        log(f"D: {tip}: n = {len(d)} intact, {len(p)} plotted, top5 {top5}, KDE bw {bw:.5f}, "
            f"peak bar {base.max():.0f}, peak curve {dens.max():.1f}", verbose)
    # legend "LTR Family": one row under the three histograms
    ly = top + 6.5 + ah + 10.0
    xx = x_left
    ftext(fig, xx, ly, "LTR Family", fontsize=FS_LEGT, ha="left", va="center_baseline")
    xx += text_w_mm("LTR Family", FS_LEGT) + 2.5
    for g, col in zip(FAM_GROUPS, FAM_COL):
        fig.add_artist(Rectangle(F(xx, ly + 0.55), 2.2 / W, 2.2 / H, transform=fig.transFigure, fc=col, ec="none"))
        ftext(fig, xx + 3.0, ly, g, fontsize=FS_LEG, ha="left", va="center_baseline")
        xx += 3.0 + text_w_mm(g, FS_LEG) + 3.0
    log(f"D: axes {aw:.1f} x {ah:.1f} mm, legend row at {ly:.1f} mm", verbose)
    return pd.concat(rows, ignore_index=True), geo


def connectors(fig, geoA, geoD, top_row2):
    """Dashed leaders from tree tips to the panel D titles (routing as in Figure_4_v2)."""
    kw = dict(color="black", lw=0.55, ls=(0, (4, 3)))
    xs = {"Crori1": geoA["chan0"], "PuccoNC29_1": geoA["chan0"] + 1.5, "myrtle_rust": geoA["chan0"] + 3.0}
    levels = {"myrtle_rust": top_row2 + 0.3, "PuccoNC29_1": top_row2 + 1.9}
    for tip, xv in xs.items():
        y_tip = geoA["tip_mm"][tip]
        x_start = geoA["label_x"] + geoA["widths"][tip] + 0.8
        g = geoD[tip]
        if tip == "Crori1":   # into the title from the left
            yb = g["title_top"] + 1.0
            fline(fig, [x_start, xv, xv, g["title_x"] - 0.8], [y_tip, y_tip, yb, yb], **kw)
        else:                 # over, then down onto the start of the title
            yl, xt = levels[tip], g["title_x"] + 1.0
            fline(fig, [x_start, xv, xv, xt, xt], [y_tip, y_tip, yl, yl, g["title_top"] - 0.6], **kw)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", default=str(ROOT / "build/fig4/inputs"), help="dir with the copied inputs")
    ap.add_argument("-o", "--out", default=str(ROOT / "fig4/fig4.pdf"), help="output PDF (180 mm wide)")
    ap.add_argument("--png", default=str(ROOT / "build/fig4/fig4_200dpi.png"),
                    help="200-dpi preview render ('' to skip)")
    ap.add_argument("--no-source-data", action="store_true", help="skip writing build/source_data/fig4")
    ap.add_argument("--extract-tip-order", metavar="PDF",
                    help="write tip_order.txt from an original figure PDF (tip label y positions) and exit")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-panel progress and sanity checks")
    a = ap.parse_args()
    inputs = Path(a.inputs)

    if a.extract_tip_order:
        import fitz
        s2 = read_s2(inputs / "supp_table2.xlsx")
        pg = fitz.open(a.extract_tip_order)[0]
        tips = []
        for b in pg.get_text("dict")["blocks"]:
            for ln in b.get("lines", []):
                for s in ln["spans"]:
                    t = s["text"].replace("ﬁ", "fi").strip()
                    if t in s2.index and s["bbox"][0] < 260:
                        tips.append(((s["bbox"][1] + s["bbox"][3]) / 2, t))
        order = [t for _, t in sorted(tips)]
        (inputs / "tip_order.txt").write_text("\n".join(order) + "\n")
        print(f"wrote {inputs / 'tip_order.txt'} ({len(order)} tips)", file=sys.stderr)
        return

    for f in ["supp_table2.xlsx", "Pucciniomycotina_dedup_busco_partitions.txt.treefile", "tip_order.txt"] + \
             [g[1] for g in D_GENOMES]:
        if not (inputs / f).exists():
            sys.exit(f"ERROR: missing input {inputs / f}")
    log("fig4: start", a.verbose, always=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ncomms_style.use()
    plt.rcParams.update({"axes.unicode_minus": True, "pdf.compression": 9})

    s2 = read_s2(inputs / "supp_table2.xlsx")
    order = [ln.strip() for ln in (inputs / "tip_order.txt").read_text().splitlines() if ln.strip()]
    newick = (inputs / "Pucciniomycotina_dedup_busco_partitions.txt.treefile").read_text().strip()
    nodes = layout_tree(parse_newick(newick), order)
    miss = [t for t in order if t not in s2.index]
    if miss:
        sys.exit(f"ERROR: tips without Table S2 row: {miss}")
    sub = s2[s2["LTR percentage (%)"].astype(float) > 1].copy()
    if len(sub) != 22:
        sys.exit(f"ERROR: expected 22 genomes with LTR content > 1 %, got {len(sub)}")

    def build(verbose):
        fig = plt.figure(figsize=(W / 25.4, H / 25.4))
        geoA = panel_A(fig, nodes, order, s2, verbose)
        letter(fig, 0.0, 0.6, "A")
        c_left = geoA["strip_right"] + 9.0
        statsC, c_bottom = panel_C(fig, sub, c_left, verbose)
        letter(fig, c_left - 8.3, 0.6, "C")
        row2 = max(geoA["tree_bottom"] + 6.5, c_bottom + 2.5)
        letter(fig, 0.0, row2 + 0.3, "B")
        scB, ldB, pveB = panel_B(fig, sub, row2, verbose)
        d_left = geoA["chan0"] + 11.5
        letter(fig, geoA["chan0"] - 4.6, row2 + 0.3, "D")
        elemD, geoD = panel_D(fig, inputs, row2, d_left, verbose)
        connectors(fig, geoA, geoD, row2)
        return fig, statsC, scB, ldB, pveB, elemD

    # pass 1 on a tall canvas (everything is placed in mm from the top), then crop the height to the content
    global H
    H = 240.0
    fig, *_ = build(False)
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    plt.close(fig)
    H = round(H - bb.y0 * 25.4 + 1.5, 1)
    if H > 240:
        sys.exit(f"ERROR: content needs {H} mm > 240 mm")
    fig, statsC, scB, ldB, pveB, elemD = build(a.verbose)

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="pdf", metadata={"Title": "Figure 4", "Creator": "make_fig4.py"})
    plt.close(fig)
    log(f"wrote {out} ({W:.0f} x {H:.0f} mm)", a.verbose, always=True)
    if a.png:
        import fitz
        Path(a.png).parent.mkdir(parents=True, exist_ok=True)
        fitz.open(out)[0].get_pixmap(dpi=200).save(a.png)
        log(f"wrote {a.png}", a.verbose)

    if a.no_source_data:
        return
    sd = SourceData("fig4")
    tipsA = s2.loc[order, ["ID", "Genome", "Class", "size_Mb"]].rename(
        columns={"Genome": "genome", "Class": "class (Table S2)", "size_Mb": "assembly_size_Mb"})
    tipsA.insert(0, "tip_order_top_to_bottom", range(1, len(order) + 1))
    chunks = [newick[i:i + 32000] for i in range(0, len(newick), 32000)]
    tipsA["newick (chunks in consecutive rows)"] = chunks + [np.nan] * (len(tipsA) - len(chunks))
    sd.add("Fig. 4A", f"Tips of the maximum-likelihood tree (n = {len(order)} genomes; IQ-TREE 2.4.0, 69 BUSCO "
           "proteins, 63,603 aa sites) with assembly size (Mb, Table S2; dot colour). Newick internal labels = "
           "SH-aLRT/UFBoot support (%); the figure shows UFBoot (1,000 replicates); branch lengths in "
           "substitutions per site.", tipsA.reset_index(drop=True))
    sd.add("Fig. 4B (scores)", "PC scores from PCA (centred, scaled to unit variance) of LTR percentage, intact "
           "count, solo count and solo:intact ratio for the n = 22 genomes with LTR content > 1%.", scB)
    sd.add("Fig. 4B (loadings)", "PCA loadings (unitless) of the 4 metrics on PC1/PC2; arrows drawn at 3 x "
           f"loading. Last row: % variance explained (PC3 {pveB[2]:.2f}%, PC4 {pveB[3]:.2f}%; not plotted).", ldB)
    cC = sub[["ID", "Genome", "size_Mb"] + [c for _, c in METRICS]].copy()
    cC.columns = ["ID", "genome", "genome_size_Mb"] + [m for m, _ in METRICS]
    sd.add("Fig. 4C", "Assembly size (Mb) and LTR metrics (LTR_percent = % of assembly; counts = number of "
           "elements; ratio = solo:intact) per genome; n = 22 genomes with LTR content > 1%.",
           cC.reset_index(drop=True))
    sd.add("Fig. 4C (statistics)", "Two-sided Pearson correlation (t test, df = n - 2; 95% CI by Fisher z) and OLS "
           "fit of each metric on genome size (Mb); P_Bonferroni_4 = P x 4 (not applied in the figure). Band = "
           "95% CI of the fitted mean.", statsC)
    sd.add("Fig. 4D", "Structurally intact LTR-RTs (EDTA) per genome; divergence = 1 - LTR identity (unitless); "
           "family = EDTA TE family ID (unassigned elements keep their own coordinates); top5_group = five most "
           "abundant families per genome; bin_start = left edge of the 0.001-wide histogram bin. Plotted "
           "(divergence < 0.101): n = 1,498 of 1,499 (C. ribicola), 1,118 of 1,118 (P. coronata 12NC29), "
           "15,093 of 15,228 (A. psidii).", elemD)
    sd.save()
    log(f"source data: {len(sd.blocks)} blocks -> {sd.dir}", a.verbose, always=True)


if __name__ == "__main__":
    main()
