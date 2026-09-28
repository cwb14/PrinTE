#!/usr/bin/env python3
"""Main-text figure - PrinTE parameter recovery test.

Promoted from Fig. S2 (v4) to a main-text figure for v5.

Six panels:
  A  searched parameter space and convergence onto the generating rates (hero panel)
  B  recovered vs generating insertion / deletion rate
  C  LTR-RT divergence distribution, true genome vs best-fit simulation
  D  held-out (zero-weight) genome summaries
  E  identifiability: best-5% interval as a fraction of the searched range
  F  sensitivity: variance in the composite score explained by each parameter

Changes relative to make_figS15.py (the supplementary version):
  * panel A enlarged - it carries the main result; B-F shrunk to make room
  * B, D, E, F transposed so the categorical variable is on x and the numeric
    variable is on y (A and C already had that orientation)
  * row-to-row spacing tightened
  * panel F axis label "Solo CV R^2" -> "CV R^2" ('solo' is reserved for solo LTRs)
  * axis labels harmonised with the other main-text figures (Fig. 6 in particular):
    rates carry "bp^-1 gen^-1", the divergence histogram y-axis is "Count"

All explanatory prose lives in the figure legend text, not in the panels.

Nature Communications final version (manuscript/v2/scripts/make_fig6.py), changes vs
make_fig_recovery.py:
  * Arial (TrueType, editable) instead of DejaVu Sans; panel letters 8 pt
  * D: the n = 33 best-5% simulations overlaid as dots on each bar
  * F: the 5 per-fold CV R^2 values overlaid as dots; error bars = sample SD (ddof=1)
  * B: the n = 33 best-5% simulations overlaid as dots on each whisker
  * writes Source Data blocks (build/source_data/fig6) and a single PDF
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold, cross_val_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ncomms_style  # noqa: E402
from srcdata import SourceData  # noqa: E402

ap = argparse.ArgumentParser(description="Build final Fig. 6 (parameter recovery test).")
ap.add_argument("--root", default="/data2/chris/PrinTE_R2", help="PrinTE_R2 data root")
ap.add_argument("-o", "--out", default="fig6/fig6.pdf", help="output PDF")
ap.add_argument("--ddof", type=int, default=1, help="ddof for panel F SD across folds (1 = sample SD)")
ap.add_argument("-v", "--verbose", action="store_true", help="print per-panel statistics")
args = ap.parse_args()

ROOT = args.root
MATRIX = f"{ROOT}/composite_matrix_rms_0_0_10_1.tsv"
REF_TSV = f"{ROOT}/burnin/gen3000000_final_LTRs_depth0_clean_ltr.tsv"
BEST_DIR = ("insertion_rates_2.83307534071e-11_deletion_rates_1.03245721314e-11"
            "_solo_ratio_15_length_bias_7")
BEST_TSV = f"{ROOT}/{BEST_DIR}/gen3000000_final_LTRs_depth0_clean_ltr.tsv"
for _p in (MATRIX, REF_TSV, BEST_TSV):
    if not os.path.exists(_p):
        sys.exit(f"ERROR: input not found: {_p}")
C_PTS = "0.15"       # overlaid data points
rng = np.random.default_rng(1)

# generating ("true") parameters: -F 3e-11,1e-11 with PrinTE defaults sr=95, k=10
TRUE = dict(insertion_rate=3e-11, deletion_rate=1e-11, solo_ratio=95, length_bias=10)
BOUNDS = dict(insertion_rate=(1e-14, 1e-10), deletion_rate=(1e-14, 1e-10),
              solo_ratio=(5, 95), length_bias=(0, 10))

C_INS = "#0072B2"    # blue         - insertion rate / best-fit simulation
C_DEL = "#D55E00"    # vermillion   - deletion rate
C_HELD = "#009E73"   # bluish green - held-out features
C_GREY = "#595959"   # dark grey    - unconstrained / weighted feature
C_LGREY = "#BFBFBF"  # light grey   - true genome, early-terminated runs
C_STAR = "#FFC20A"   # gold         - best fit (matches Fig. 6)

RATE_UNIT = "bp$^{-1}$ gen$^{-1}$"   # as in Fig. 6E,G,I,K,M,N
ROT = 30                             # x tick label rotation, as in Fig. 6M-P
VMAX = 0.70                          # composite-score colour cap (~80th pct of
                                     # the scored runs; the tail reaches 2.03, so
                                     # the top tick is labelled ">=" rather than
                                     # drawing an extend arrow on the colour bar)

ncomms_style.arial()
plt.rcParams.update({
    "font.size": 7,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 6.8,
    "legend.fontsize": 6.2,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ---------------------------------------------------------------- load data
df = pd.read_csv(MATRIX, sep="\t")
imp = df["exp_genome_size"] == 0
real = df[~imp]
best = df.loc[df["Composite"].idxmin()]
thr = df["Composite"].quantile(0.05)
top = df[df["Composite"] <= thr]

ref_div = pd.read_csv(REF_TSV, sep="\t")["K2P_d"].astype(float).values
sim_div = pd.read_csv(BEST_TSV, sep="\t")["K2P_d"].astype(float).values

print(f"combinations={len(df)}  scored={len(real)}  best5%={len(top)} (thr={thr:.4g})")
print(f"ref LTR-RTs={len(ref_div)}  best-fit LTR-RTs={len(sim_div)}")
N_TOP = len(top)


def strip(ax, x, vals, half=0.1, s=5, color=C_PTS, alpha=0.55, z=5):
    """Overlay individual data points as a jittered strip centred on x."""
    vals = np.asarray(vals, float)
    ax.scatter(x + rng.uniform(-half, half, len(vals)), vals, s=s, color=color, alpha=alpha,
               lw=0, zorder=z, clip_on=False)

# ---------------------------------------------------------------- layout
# 2 x 3 grid.  Column 0 is ~1.8x the others so panel A (the main result) is the
# largest panel; the row gap is kept short so A sits close to D.
FIGW, FIGH = 7.1, 5.6
L, R, T, B_ = 0.095, 0.988, 0.955, 0.135
fig = plt.figure(figsize=(FIGW, FIGH))
gs = fig.add_gridspec(2, 3,
                      width_ratios=[2.28, 1.34, 1.38],
                      height_ratios=[1.40, 1.00],
                      wspace=0.55, hspace=0.272,
                      left=L, right=R, top=T, bottom=B_)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[0, 2])
axD = fig.add_subplot(gs[1, 0])
axE = fig.add_subplot(gs[1, 1])
axF = fig.add_subplot(gs[1, 2])

letters = []

# =========================================================== A: search space
axA.scatter(np.log10(df.loc[imp, "insertion_rate"]),
            np.log10(df.loc[imp, "deletion_rate"]),
            s=6, c=C_LGREY, lw=0, alpha=0.9, zorder=1)
vmax = VMAX
sc = axA.scatter(np.log10(real["insertion_rate"]), np.log10(real["deletion_rate"]),
                 c=real["Composite"], cmap="viridis_r", s=11, lw=0,
                 vmin=real["Composite"].min(), vmax=vmax, zorder=2)
tx, ty = np.log10(TRUE["insertion_rate"]), np.log10(TRUE["deletion_rate"])
axA.axvline(tx, color="k", lw=0.6, ls=(0, (4, 2)), zorder=3)
axA.axhline(ty, color="k", lw=0.6, ls=(0, (4, 2)), zorder=3)
axA.set_xlim(-14.3, -9.7)
axA.set_ylim(-14.3, -9.7)
axA.set_xticks([-14, -13, -12, -11, -10])
axA.set_yticks([-14, -13, -12, -11, -10])
axA.set_xlabel(f"log$_{{10}}$ (insertion rate, {RATE_UNIT})")
axA.set_ylabel(f"log$_{{10}}$ (deletion rate, {RATE_UNIT})")
cb = fig.colorbar(sc, ax=axA, fraction=0.042, pad=0.026)   # blunt at both ends
cb.set_label("Composite score", fontsize=6.6, labelpad=1.2)
cb_ticks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
cb.set_ticks(cb_ticks)
cb.set_ticklabels([("\u2265" if t == VMAX else "") + f"{t:.1f}" for t in cb_ticks])
cb.ax.tick_params(labelsize=6.0, width=0.6, length=2, pad=1.1)
cb.outline.set_linewidth(0.6)
n_capped = int((real["Composite"] > VMAX).sum())
print(f"  A colour scale capped at {VMAX}: {n_capped}/{len(real)} scored runs "
      f"share the darkest colour (max score {real['Composite'].max():.2f})")

axIn = inset_axes(axA, width="45%", height="45%", loc="lower left",
                  bbox_to_anchor=(0.045, 0.045, 1, 1), bbox_transform=axA.transAxes)
axIn.scatter(np.log10(real["insertion_rate"]), np.log10(real["deletion_rate"]),
             c=real["Composite"], cmap="viridis_r", s=20, lw=0,
             vmin=real["Composite"].min(), vmax=vmax)
axIn.axvline(tx, color="k", lw=0.6, ls=(0, (4, 2)))
axIn.axhline(ty, color="k", lw=0.6, ls=(0, (4, 2)))
axIn.plot(tx, ty, marker="+", ms=10, mew=1.5, color="k", zorder=6)
axIn.plot(np.log10(best["insertion_rate"]), np.log10(best["deletion_rate"]),
          marker="*", ms=12, color=C_STAR, mec="k", mew=0.6, zorder=7)
axIn.set_xlim(-10.78, -10.32)
axIn.set_ylim(-11.42, -10.58)
axIn.set_xticks([-10.7, -10.4])
axIn.set_yticks([-11.2, -10.8])
axIn.tick_params(labelsize=6.0, width=0.5, length=1.8, pad=1.2)
for s in axIn.spines.values():
    s.set_linewidth(0.6)
    s.set_color("0.35")
axIn.set_facecolor("white")
axIn.patch.set_alpha(1.0)
mark_inset(axA, axIn, loc1=2, loc2=4, fc="none", ec="0.55", lw=0.5, zorder=4)

leg = axA.legend(handles=[
    Line2D([], [], marker="+", ls="", color="k", ms=6, mew=1.2, label="True genome"),
    Line2D([], [], marker="*", ls="", color=C_STAR, mec="k", mew=0.5, ms=7.5, label="Best fit"),
    Line2D([], [], marker="o", ls="", color=C_LGREY, ms=3.4, label="Terminated"),
], loc="upper left", frameon=True, framealpha=0.9, edgecolor="none", fontsize=6.5,
    handletextpad=0.3, borderpad=0.18, labelspacing=0.24, handlelength=1.0,
    borderaxespad=0.28)
leg.get_frame().set_facecolor("white")
letters.append((axA, "A"))

# ====================================================== B: rate recovery
rows = [("Insertion", "insertion_rate", C_INS),
        ("Deletion", "deletion_rate", C_DEL)]
xpos = [0.0, 1.0]
# plotted as log10 on a linear axis so the axis matches panel A and Fig. 6E,G,I,K
for x, (name, key, col) in zip(xpos, rows):
    lo_s, hi_s = BOUNDS[key]
    axB.bar(x, np.log10(hi_s) - np.log10(lo_s), bottom=np.log10(lo_s), width=0.135,
            color=C_LGREY, lw=0, zorder=1)
    lo, hi = np.log10(top[key].min()), np.log10(top[key].max())
    bf, tv = np.log10(best[key]), np.log10(TRUE[key])
    axB.errorbar([x], [bf], yerr=[[bf - lo], [hi - bf]], fmt="none",
                 ecolor=col, elinewidth=1.8, capsize=2.8, capthick=1.4, zorder=4)
    strip(axB, x + 0.2, np.log10(top[key]), half=0.06, s=4, color=col, alpha=0.6)
    axB.plot([x], [tv], marker="_", ms=13, mew=1.8, color="k", zorder=5)
    axB.plot([x], [bf], marker="o", ms=4.6, color=col, mec="white", mew=0.7, zorder=6)
    err = 100 * (best[key] - TRUE[key]) / TRUE[key]
    axB.annotate(f"{err:+.1f}%", xy=(x, np.log10(hi_s)), xytext=(0, 3.5),
                 textcoords="offset points", ha="center", va="bottom", fontsize=6.8,
                 fontweight="bold", color=col)
axB.set_ylim(-14.85, -9.35)
axB.set_xlim(-0.68, 1.68)
axB.set_xticks(xpos)
axB.set_xticklabels([r[0] for r in rows])
axB.set_ylabel(f"log$_{{10}}$ (rate, {RATE_UNIT})")
axB.set_yticks([-14, -13, -12, -11, -10])
axB.grid(axis="y", color="0.88", lw=0.5, zorder=0)
axB.set_axisbelow(True)
lgB = axB.legend(handles=[
    Patch(fc=C_LGREY, ec="none", label="Searched range"),
    Line2D([], [], marker="_", ls="", color="k", ms=7, mew=1.5, label="True genome"),
    Line2D([], [], color="0.35", lw=1.8, marker="o", ms=4, mec="white", mew=0.6,
           label="Best fit ± best 5%"),
], loc="lower center", frameon=True, framealpha=0.92, edgecolor="none", ncol=1,
    handletextpad=0.4, borderpad=0.16, labelspacing=0.24, handlelength=1.3,
    borderaxespad=0.3)
lgB.get_frame().set_facecolor("white")
letters.append((axB, "B"))

# ============================================ C: divergence distributions
hi = np.nanpercentile(np.concatenate([ref_div, sim_div]), 99.7)
bins = np.linspace(0, hi, 51)
axC.hist(ref_div, bins=bins, color=C_LGREY, label="True genome", zorder=1)
axC.hist(sim_div, bins=bins, histtype="step", color=C_INS, lw=1.0,
         label="Best fit", zorder=2)
axC.set_xlabel("K2P divergence")
axC.set_ylabel("Count")
axC.set_xlim(0, hi)
axC.set_xticks([0, 0.05, 0.10, 0.15, 0.20])
axC.set_ylim(0, max(np.histogram(ref_div, bins=bins)[0]) * 1.16)
axC.legend(loc="upper right", frameon=False, handletextpad=0.5,
           borderpad=0.0, labelspacing=0.3, borderaxespad=0.3)
letters.append((axC, "C"))

# ================================================= D: held-out features
feats = [
    ("LTR-RT count", "ref_ltr_rt_count", "exp_ltr_rt_count", C_HELD),
    ("LTR-RT length", "ref_cumulative_length", "exp_cumulative_length", C_HELD),
    ("Genome size", "ref_genome_size", "exp_genome_size", C_GREY),
]
xp = np.arange(len(feats)).astype(float)
axD.axhspan(-5, 5, color="0.925", zorder=0, lw=0)
for x, (name, refc, expc, col) in zip(xp, feats):
    dev = 100 * (best[expc] - best[refc]) / best[refc]
    rel = 100 * (top[expc] - top[refc]) / top[refc]
    lo, hi_ = rel.min(), rel.max()
    axD.bar(x, dev, width=0.34, color=col, alpha=0.75, zorder=2)
    axD.errorbar([x], [dev], yerr=[[dev - lo], [hi_ - dev]], fmt="none",
                 ecolor="0.25", elinewidth=0.9, capsize=2.2, capthick=0.9, zorder=4)
    strip(axD, x, rel, half=0.11)
    axD.text(x, 12.4, f"{dev:+.1f}%", ha="center", va="top",
             fontsize=6.8, fontweight="bold", color=col, zorder=5)
    print(f"  D {name}: best={dev:+.2f}%  top5% [{lo:+.2f}, {hi_:+.2f}]")
axD.axhline(0, color="k", lw=0.7, zorder=3)
axD.set_xticks(xp)
axD.set_xticklabels([f[0] for f in feats], rotation=ROT, ha="right",
                    rotation_mode="anchor")
axD.set_ylim(-37, 14)
axD.set_yticks([-30, -20, -10, 0, 10])
axD.set_xlim(-0.62, len(feats) - 0.38)
axD.set_ylabel("Deviation from\ntrue genome (%)")
axD.legend(handles=[
    Patch(fc=C_HELD, ec="none", label="Held out"),
    Patch(fc=C_GREY, ec="none", label="Scored"),
], loc="lower right", frameon=False, ncol=1, handletextpad=0.4, borderpad=0.0,
    labelspacing=0.28, handlelength=1.1, borderaxespad=0.35)
axD.grid(axis="y", color="0.85", lw=0.5, zorder=1)
axD.set_axisbelow(True)
letters.append((axD, "D"))

# ================================================= E: identifiability
def frac_of_range(key, logscale):
    lo_s, hi_s = BOUNDS[key]
    lo, hi_ = top[key].min(), top[key].max()
    if logscale:
        return 100 * np.log10(hi_ / lo) / np.log10(hi_s / lo_s)
    return 100 * (hi_ - lo) / (hi_s - lo_s)

items = [("Insertion rate", frac_of_range("insertion_rate", True), C_INS),
         ("Deletion rate", frac_of_range("deletion_rate", True), C_DEL),
         ("Solo LTR rate", frac_of_range("solo_ratio", False), C_GREY),
         ("Length bias", frac_of_range("length_bias", False), C_GREY)]
xp = np.arange(len(items)).astype(float)
for x, (name, val, col) in zip(xp, items):
    axE.bar(x, val, width=0.62, color=col, zorder=2)
    axE.text(x, val + 3.5, f"{val:.0f}%", ha="center", va="bottom", fontsize=6.8,
             fontweight="bold", color=col, zorder=3)
axE.set_xticks(xp)
axE.set_xticklabels([i[0] for i in items], rotation=ROT, ha="right",
                    rotation_mode="anchor")
axE.set_ylim(0, 128)
axE.set_yticks([0, 25, 50, 75, 100])
axE.set_xlim(-0.66, len(items) - 0.34)
axE.set_ylabel("Best-5% interval\n(% of range)")
axE.grid(axis="y", color="0.88", lw=0.5, zorder=1)
axE.set_axisbelow(True)
letters.append((axE, "E"))
for name, val, _ in items:
    print(f"  E {name}: {val:.1f}% of searched range")

# ================================================= F: sensitivity
params = ["insertion_rate", "deletion_rate", "solo_ratio", "length_bias"]
disp = {"insertion_rate": "Insertion rate", "deletion_rate": "Deletion rate",
        "solo_ratio": "Solo LTR rate", "length_bias": "Length bias"}
cols = {"insertion_rate": C_INS, "deletion_rate": C_DEL,
        "solo_ratio": C_GREY, "length_bias": C_GREY}
Xc = []
for p in params:
    v = real[p].astype(float).values
    Xc.append(np.log10(v) if "rate" in p else v)
X = np.column_stack(Xc)
y = real["Composite"].astype(float).values
ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
X, y = X[ok], y[ok]
cv = KFold(n_splits=5, shuffle=True, random_state=0)
rf = RandomForestRegressor(n_estimators=400, oob_score=True, random_state=0, n_jobs=-1)
joint = cross_val_score(rf, X, y, cv=cv, scoring="r2")
rf.fit(X, y)
perm = permutation_importance(rf, X, y, n_repeats=20, random_state=0,
                              n_jobs=-1, scoring="r2")
single, single_sd, folds = [], [], []
for i in range(len(params)):
    s = cross_val_score(RandomForestRegressor(n_estimators=300, random_state=0,
                                              n_jobs=-1),
                        X[:, [i]], y, cv=cv, scoring="r2")
    single.append(s.mean())
    single_sd.append(s.std(ddof=args.ddof))
    folds.append(s)
print(f"  joint CV R2 = {joint.mean():.4f} +/- {joint.std():.4f}  OOB={rf.oob_score_:.4f}")
for p, s, sd, pm, ps in zip(params, single, single_sd, perm.importances_mean,
                            perm.importances_std):
    print(f"  {p}: single-parameter CV R2={s:.4f}+/-{sd:.4f}  perm={pm:.4f}+/-{ps:.4f}")

xp = np.arange(len(params)).astype(float)
for x, p, s, sd, fv in zip(xp, params, single, single_sd, folds):
    axF.bar(x, s, width=0.62, color=cols[p], alpha=0.75, zorder=2)
    axF.errorbar([x], [s], yerr=[[sd], [sd]], fmt="none", ecolor="0.25",
                 elinewidth=0.9, capsize=2.2, capthick=0.9, zorder=4)
    strip(axF, x, fv, half=0.16, s=7)
    if s >= 0:
        axF.text(x, max(s + sd, fv.max()) + 0.030, f"{s:.2f}", ha="center", va="bottom",
                 fontsize=6.8, fontweight="bold", color=cols[p], zorder=5)
    else:
        axF.text(x, min(s - sd, fv.min()) - 0.030, f"{s:.2f}", ha="center", va="top",
                 fontsize=6.8, fontweight="bold", color=cols[p], zorder=5)
axF.axhline(0, color="k", lw=0.7, zorder=3)
axF.set_xticks(xp)
axF.set_xticklabels([disp[p] for p in params], rotation=ROT, ha="right",
                    rotation_mode="anchor")
axF.set_ylim(-0.32, 0.95)
axF.set_yticks([-0.25, 0, 0.25, 0.5, 0.75])
axF.set_xlim(-0.66, len(params) - 0.34)
axF.set_ylabel("CV $R^2$")
axF.grid(axis="y", color="0.88", lw=0.5, zorder=1)
axF.set_axisbelow(True)
letters.append((axF, "F"))

# --------------------------------------------------------- panel letters
# left-aligned with each panel's y-axis furniture (label + tick labels), so the
# letters line up down the figure without hand-tuned offsets.
fig.canvas.draw()
rend = fig.canvas.get_renderer()
inv = fig.transFigure.inverted()

def furniture_x0(ax):
    xs = [ax.get_window_extent(renderer=rend).x0]
    lb = ax.yaxis.get_label()
    if lb.get_text():
        xs.append(lb.get_window_extent(renderer=rend).x0)
    for t in ax.get_yticklabels():
        if t.get_text():
            xs.append(t.get_window_extent(renderer=rend).x0)
    return min(xs)

col_x0 = {}
for ax, letter in letters:
    key = round(ax.get_position().x0, 3)
    x0 = furniture_x0(ax)
    col_x0[key] = min(col_x0.get(key, x0), x0)
for ax, letter in letters:
    key = round(ax.get_position().x0, 3)
    xfig = inv.transform((col_x0[key], 0))[0]
    fig.text(max(0.004, xfig - 0.012), ax.get_position().y1 + 0.008, letter,
             fontsize=8, fontweight="bold", va="bottom", ha="left")

# ------------------------------------------------- overflow check (verification)
fig.canvas.draw()
rend = fig.canvas.get_renderer()
problems = []
panels = [("A", axA), ("B", axB), ("C", axC), ("D", axD), ("E", axE), ("F", axF)]
for name, ax in panels:
    ab = ax.get_window_extent(renderer=rend)
    items_ = [(t, "text:%r" % t.get_text()[:22]) for t in ax.texts]
    lg = ax.get_legend()
    if lg is not None:
        items_.append((lg, "legend"))
    for art, what in items_:
        bb = art.get_window_extent(renderer=rend)
        if bb.x0 < ab.x0 - 0.5 or bb.x1 > ab.x1 + 0.5 or bb.y0 < ab.y0 - 0.5 or bb.y1 > ab.y1 + 0.5:
            over = max(ab.x0 - bb.x0, bb.x1 - ab.x1, ab.y0 - bb.y0, bb.y1 - ab.y1)
            problems.append(f"  panel {name}: {what} overflows axes by {over:.1f} px")
# axis labels must fit along the side they label
for name, ax in panels:
    ab = ax.get_window_extent(renderer=rend)
    yl = ax.yaxis.get_label()
    if yl.get_text():
        h = yl.get_window_extent(renderer=rend).height
        if h > ab.height:
            problems.append(f"  panel {name}: ylabel is {h - ab.height:.0f} px taller than the axes")
    xl = ax.xaxis.get_label()
    if xl.get_text():
        w = xl.get_window_extent(renderer=rend).width
        if w > ab.width:
            problems.append(f"  panel {name}: xlabel is {w - ab.width:.0f} px wider than the axes")
# neighbouring panels (incl. tick labels) must not collide
tb = {n: ax.get_tightbbox(rend) for n, ax in panels}
for a, b in [("A", "B"), ("B", "C"), ("D", "E"), ("E", "F")]:
    if tb[a].x1 > tb[b].x0 + 0.5:
        problems.append(f"  panels {a}/{b} overlap horizontally by {tb[a].x1 - tb[b].x0:.1f} px")
for a, b in [("A", "D"), ("B", "E"), ("C", "F")]:
    if tb[b].y1 > tb[a].y0 + 0.5:
        problems.append(f"  panels {a}/{b} overlap vertically by {tb[b].y1 - tb[a].y0:.1f} px")
cb_bb = cb.ax.get_tightbbox(rend)
if cb_bb.x1 > tb["B"].x0 + 0.5:
    problems.append(f"  colorbar/panel B overlap by {cb_bb.x1 - tb['B'].x0:.1f} px")
print("OVERFLOW CHECK:", "clean" if not problems else "")
for pr in problems:
    print(pr)

os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
fig.savefig(args.out, dpi=600)
print("wrote", args.out)

# ---------------------------------------------------------------- Source Data
sd = SourceData("fig6")
status = np.where(imp, "terminated", "scored")
blkA = df[["insertion_rate", "deletion_rate", "solo_ratio", "length_bias", "Composite"]].copy()
blkA.insert(0, "simulation", np.arange(1, len(df) + 1))
blkA["status"] = status
blkA["in_best_5pct"] = df["Composite"] <= thr
blkA["is_best_fit"] = df.index == best.name
sd.add("Fig. 6A", f"All {len(df)} parameter combinations evaluated by the Grid Search Method "
       f"({int((~imp).sum())} scored, {int(imp.sum())} terminated before finishing; terminated runs carry a "
       f"fixed penalty score and are shown grey). Rates in bp^-1 gen^-1. Best 5% = Composite <= "
       f"{thr:.5f} (5th percentile of all {len(df)}), n = {N_TOP}.", blkA)
rowsB = []
for name, key in [("Insertion rate", "insertion_rate"), ("Deletion rate", "deletion_rate")]:
    rowsB.append({"parameter": name, "true_value": TRUE[key], "best_fit": best[key],
                  "best5pct_min": top[key].min(), "best5pct_max": top[key].max(),
                  "n_best5pct": N_TOP, "searched_min": BOUNDS[key][0], "searched_max": BOUNDS[key][1],
                  "best_fit_deviation_pct": 100 * (best[key] - TRUE[key]) / TRUE[key]})
sd.add("Fig. 6B", f"Recovered vs true rates (bp^-1 gen^-1; plotted as log10). Point = best fit; whiskers = "
       f"min-max over the best 5% of simulations (n = {N_TOP}; individual values in Fig. 6A block, "
       f"in_best_5pct = TRUE).", pd.DataFrame(rowsB))
blkC = pd.DataFrame({"genome": ["true genome"] * len(ref_div) + ["best-fit simulation"] * len(sim_div),
                     "K2P_divergence": np.concatenate([ref_div, sim_div])})
sd.add("Fig. 6C", f"K2P divergence of every intact LTR-RT (true genome n = {len(ref_div)}; best-fit "
       f"simulation n = {len(sim_div)}). Histogram: 50 equal bins from 0 to {hi:.5f} (99.7th percentile "
       f"of pooled values); values above are not drawn.", blkC)
blkD = top[["insertion_rate", "deletion_rate", "solo_ratio", "length_bias", "Composite"]].copy()
for name, refc, expc, _ in feats:
    blkD[f"{name} true"] = top[refc]
    blkD[f"{name} simulated"] = top[expc]
    blkD[f"{name} deviation (%)"] = 100 * (top[expc] - top[refc]) / top[refc]
blkD["is_best_fit"] = top.index == best.name
blkD = blkD.sort_values("Composite")
sd.add("Fig. 6D", f"Deviation (%) = 100 x (simulated - true) / true for each of the n = {N_TOP} best-5% "
       f"simulations (dots). Bar = best-fit simulation (is_best_fit = TRUE); whiskers = min-max of the "
       f"{N_TOP}. LTR-RT count and length were held out of the composite score; genome size was scored.", blkD)
rowsE = []
for (name, val, _), key, lg in zip(items, ["insertion_rate", "deletion_rate", "solo_ratio", "length_bias"],
                                   [True, True, False, False]):
    rowsE.append({"parameter": name, "best5pct_min": top[key].min(), "best5pct_max": top[key].max(),
                  "searched_min": BOUNDS[key][0], "searched_max": BOUNDS[key][1],
                  "scale": "log10" if lg else "linear", "interval_width_pct_of_range": val})
sd.add("Fig. 6E", f"Width of the best-5% interval (n = {N_TOP} simulations) as % of the searched range "
       f"(log10 scale for rates). Single values; no error bars.", pd.DataFrame(rowsE))
blkF = pd.DataFrame({"parameter": [disp[p] for p in params]})
for k in range(5):
    blkF[f"fold{k + 1}_R2"] = [f[k] for f in folds]
blkF["mean_R2"] = single
blkF[f"SD_R2 (ddof={args.ddof})"] = single_sd
sd.add("Fig. 6F", f"Single-parameter cross-validated R^2 of the composite score (random forest, 300 trees; "
       f"5-fold CV, shuffled, random_state 0) on the n = {len(y)} scored simulations. Bar = mean of the 5 "
       f"out-of-fold R^2 values (dots); error bars = SD across folds.", blkF)
sd.save()
print("source data blocks:", len(sd.blocks))
