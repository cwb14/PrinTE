"""Shared matplotlib style for Nature Communications final figures.

- Arial (registered from the mscorefonts package cache; nothing installed)
- TrueType (Type 42) embedding so text stays editable in Illustrator
- 7 pt text, 8 pt bold upper-case panel letters
- Okabe-Ito colour-blind-safe palette
"""
from pathlib import Path

import matplotlib as mpl
from matplotlib import font_manager

MM = 1 / 25.4  # inches per mm
FULL_W = 180 * MM  # double-column width
HALF_W = 88 * MM   # single-column width
ARIAL_DIR = Path("/home/chris/mamba/pkgs/mscorefonts-0.0.1-3/fonts")
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def _fix_mathtext_xheight():
    """matplotlib reads x-height from Arial's PCLT table with a wrong unit scale (~2x too big),
    which floats super/subscripts (bp^-1, R^2, log_10) far too high. Use the measured 'x'
    glyph instead - matplotlib's own fallback for fonts without a PCLT table."""
    import matplotlib._mathtext as mt

    def get_xheight(self, fontname, fontsize, dpi):
        return self.get_metrics(fontname, mpl.rcParams["mathtext.default"], "x", fontsize, dpi).iceberg

    mt.TruetypeFonts.get_xheight = get_xheight


def arial():
    """Register Arial and set only font-related rcParams (for scripts with their own tuned style)."""
    _fix_mathtext_xheight()
    for f in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
        p = ARIAL_DIR / f
        if not p.exists():
            raise FileNotFoundError(f"Arial font missing: {p}")
        font_manager.fontManager.addfont(str(p))
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.sf": "Arial",
        "mathtext.cal": "Arial",
        "mathtext.default": "regular",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def use(base_size=7):
    """Register Arial and apply the full house style. Call once before plotting."""
    arial()
    mpl.rcParams.update({
        "font.size": base_size,
        "axes.titlesize": base_size,
        "axes.labelsize": base_size,
        "xtick.labelsize": base_size - 1,
        "ytick.labelsize": base_size - 1,
        "legend.fontsize": base_size - 1,
        "legend.title_fontsize": base_size - 1,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "patch.linewidth": 0.5,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 600,
        "figure.dpi": 150,
    })


def panel_label(ax, letter, x=-0.18, y=1.06, fig=None, size=8):
    """Bold upper-case panel letter at the upper-left of an axes (axes coords)."""
    target = fig if fig is not None else ax
    kw = dict(transform=target.transFigure if fig is not None else ax.transAxes)
    target.text(x, y, letter, fontsize=size, fontweight="bold", va="bottom", ha="left", **kw)
