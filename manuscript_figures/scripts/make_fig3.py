#!/usr/bin/env python3
"""Build final Fig. 3 (TE-driven genome contraction/expansion) for Nature Communications.

Rebuilds slide 7 of fig3.pptx (V3.3.1) from the original source data as ONE vector PDF,
180 mm wide, all text Arial (TrueType, editable):
  A  chr ribbon plot + per-generation pies      (riparian3.py / plot_TE_frac2.py logic, imported read-only)
  B,C LTR-RT K2P density (seaborn KDE)          (ltr_dens.py load_data + identical kdeplot calls)
  D,E Tekay RT trees                            (make_fig3_trees.R: ggtree circular, as LTR_phylo7.R)
  F,J TE insertion/deletion counts and rates    (Rscript2.R --collapse, redrawn in matplotlib)
  G  genome size, H TE counts per superfamily   (genome_plot.py / plot_superfamily_count.py)
  I  LTR-RT K2P density, A. thaliana replica    (ltr_dens.py)
Source scripts are exec'd from their original locations without writing anything there.
Also writes the Source Data blocks (build/source_data/fig3) and a 200-dpi PNG preview.
"""
import argparse
import glob
import os
import re
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True  # never drop __pycache__ next to (read-only) source scripts

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator, NullFormatter  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import ncomms_style  # noqa: E402
from srcdata import SourceData  # noqa: E402

V29 = ("/home/chris/data/TEGE/v2/v3/v4/v5/v6/v7/v8/v9/v10/v11/v12/v13/v14/v15/v16/v17/v18/v19/v20/"
       "v21/v22/v23/v24/v25/v26/v27/v28/v29")
CE = f"{V29}/figures/contraction_expansion"   # 10 M-generation contraction -> expansion run
RIP = f"{CE}/riparian"                         # synteny anchors between generations
VM = f"{V29}/variable_mode4/4e-6"               # A. thaliana replica, variable-rate mode
TK = f"{V29}/LTR_phylogeny/Tekay2"             # Tekay phylogenies
PRINTE = f"{V29}/TESS/prinTE"
TREE_D = f"{TK}/contraction/both_4000000_LTR.fa.rexdb-plant.cls.pep.RT.aln.nwk"
TREE_E = f"{TK}/expansion/both_4000000_LTR_subset.fa.rexdb-plant.cls.pep.RT.aln.nwk"
RSCRIPT = "/home/chris/bin/mambaforge/envs/synLTR/bin/Rscript"

FW, FH = 180.0, 152.0          # figure size (mm)
TXT, LAB, LEG = 6.0, 7.0, 6.0  # tick / axis-label / legend font sizes (pt)
RATE_LAB = 7.2  # rate-axis title: its 0.7x mathtext superscripts must stay >= 5 pt
INS, DEL = "#7570b3", "#d95f02"  # Rscript2.R insertion / deletion colours
FEAT_KEYS = ["Intact TE", "SoloLTR", "Fragmented TE"]
TREE_MM_PER_UNIT = 6.2         # both trees share one radial scale (published: 6.2-6.4 mm/unit at 180 mm)
TREE_RADIUS = 4.5              # ggtree radial limit (> max root-to-tip of both trees)
TREE_LW_PT = 0.4


def log(msg, verbose=True):
    if verbose:
        print(msg, flush=True)


def load_script(path, name):
    """Exec an original plotting script as a module (no bytecode cache, nothing written)."""
    if not os.path.isfile(path):
        sys.exit(f"ERROR: missing source script {path}")
    mod = types.ModuleType(name)
    mod.__file__ = path
    exec(compile(Path(path).read_text(), path, "exec"), mod.__dict__)
    return mod


def need(*paths):
    for p in paths:
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing input {p}")


# ----------------------------------------------------------------------------- layout helpers (mm, top-left)
def ax_mm(fig, x, y, w, h, **kw):
    return fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH], **kw)


def text_mm(fig, x, y, s, **kw):
    return fig.text(x / FW, 1 - y / FH, s, **kw)


def letter(fig, x, y, s):
    text_mm(fig, x, y, s, fontsize=8, fontweight="bold", va="top", ha="left")


def box(ax, lw=0.6):
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(lw)


def sig4(v):
    return float(f"{v:.4g}")


def gfmt(v, _pos=None):
    return f"{v:g}"


# ----------------------------------------------------------------------------- panel A
def pies_table(ptf, fai_len, verbose):
    """Per-generation bp of each feature (plot_TE_frac2.process_bed_file) and genome length (all.fai)."""
    beds = [("burnin.bed", 0)] + [(f"gen{i}000000_final.bed", i * 1_000_000) for i in range(1, 11)]
    rows = []
    for bed, gen in beds:
        path = f"{CE}/expansion/{bed}"
        need(path)
        acc = "gen0000000" if gen == 0 else f"gen{gen}"
        acc = acc if acc in fai_len else f"gen{gen:07d}"
        g = fai_len[acc]
        fl = ptf.process_bed_file(path)
        row = dict(generation=gen, bed_file=bed, genome_bp=g)
        for k in FEAT_KEYS:
            row[f"{k.replace(' ', '_')}_bp"] = fl[k]
        for k in FEAT_KEYS:
            row[f"{k.replace(' ', '_')}_pct"] = fl[k] / g * 100
        row["Non-TE_pct"] = 100 - sum(row[f"{k.replace(' ', '_')}_pct"] for k in FEAT_KEYS)
        rows.append(row)
    df = pd.DataFrame(rows)
    df[[c for c in df.columns if c.endswith("_pct")]] = df[[c for c in df.columns if c.endswith("_pct")]].round(4)
    rep = f"{CE}/expansion/pipeline_all.report.revised.cleaned"
    if os.path.exists(rep):  # cross-check fai lengths vs FASTA bp counted by calculate_rate.py
        r = pd.read_csv(rep, sep="\t").drop_duplicates("Generation").set_index("Generation")["genome_size"]
        bad = [(g, s, r[g]) for g, s in zip(df.generation, df.genome_bp) if g in r.index and r[g] != s]
        if bad:
            print(f"WARNING: genome length mismatch fai vs report: {bad}", file=sys.stderr)
        else:
            log(f"  A pies: fai genome lengths match FASTA bp counts for {sum(df.generation.isin(r.index))} "
                "generations", verbose)
    return df


def panel_A(fig, rip, ptf, chrom, chrom_label, verbose):
    need(f"{RIP}/all.fai", f"{RIP}/all.recip.anchors.coords", f"{RIP}/names.txt")
    chrom_info = rip.parse_chrom_size(f"{RIP}/all.fai")
    order = rip.read_accession_order(f"{RIP}/names.txt", chrom_info)
    fai_len = {a: sum(v.values()) for a, v in chrom_info.items()}
    blocks = rip.merge_blocks(rip.parse_syntenic_file(f"{RIP}/all.recip.anchors.coords"))  # default merge
    one = {a: {chrom: chrom_info[a][chrom]} for a in order}
    pos, _, _, ypos = rip.compute_chromosome_positions(one, order, alignment="center", acc_gap=50)
    offset = 0.07 * 50
    idx = {a: i for i, a in enumerate(order)}
    polys = []
    for b in blocks:  # same selection/geometry as riparian3.draw_ribbons (--disable_paralogs)
        if abs(idx[b["acc1"]] - idx[b["acc2"]]) != 1 or b["seq1"] != chrom or b["seq2"] != chrom:
            continue
        x1s, x1e, y1, L1 = pos[b["acc1"]][chrom]
        x2s, x2e, y2, L2 = pos[b["acc2"]][chrom]
        r1s = b["start1"] / L1 * (x1e - x1s) + x1s
        r1e = b["end1"] / L1 * (x1e - x1s) + x1s
        r2s = b["start2"] / L2 * (x2e - x2s) + x2s
        r2e = b["end2"] / L2 * (x2e - x2s) + x2s
        y1e, y2e = (y1 - offset, y2 + offset) if y1 > y2 else (y1 + offset, y2 - offset)
        if b["strand"] == "+":
            polys.append([(r1s, y1e), (r1e, y1e), (r2e, y2e), (r2s, y2e)])
        else:
            polys.append([(r1s, y1e), (r1e, y1e), (r2s, y2e), (r2e, y2e)])
    log(f"  A ribbon: {chrom}, {len(polys)} syntenic ribbons between adjacent generations", verbose)

    x0s = [pos[a][chrom][0] for a in order]
    x1s_ = [pos[a][chrom][1] for a in order]
    xmin, xmax = min(x0s), max(x1s_)
    ax = ax_mm(fig, 9.5, 3.5, 34.5, 76.0)
    ax.set_xlim(xmin - 0.02 * (xmax - xmin), xmax + 0.02 * (xmax - xmin))
    ax.set_ylim(-30, 506)
    ax.axis("off")
    col = matplotlib.colors.TABLEAU_COLORS["tab:blue"]  # riparian3: chr1 and chr2 both tab:blue
    for a in order:
        xs, xe, y, _ = pos[a][chrom]
        ax.plot([xs, xe], [y, y], color=col, lw=2.2, solid_capstyle="round", zorder=2)
    ax.add_collection(PolyCollection(polys, facecolor="grey", edgecolor="none", alpha=0.4, zorder=3))
    for a in order:  # generation labels (millions) at the left end of each bar
        xs, _, y, _ = pos[a][chrom]
        g = int(a.replace("gen", ""))
        ax.annotate(f"{g // 1_000_000}", (xs, y), xytext=(-3.5, 0), textcoords="offset points",
                    ha="right", va="center", fontsize=TXT)
    text_mm(fig, 1.2, 3.5 + 76.0 * (506 - 250) / 536, "Millions of Generations", rotation=90,
            ha="left", va="center", fontsize=LAB)
    # 50 Mb scale bar and chromosome label (riparian3 geometry)
    y10 = min(ypos.values())
    xs10 = pos[order[-1]][chrom][0]
    ysb = y10 - offset * 2.5
    ax.plot([xs10, xs10 + 50e6], [ysb, ysb], color="black", lw=1.0, solid_capstyle="butt")
    ax.text(xs10 + 25e6, ysb - offset * 1.2, "50 Mb", ha="center", va="top", fontsize=TXT)
    xm = (pos[order[-1]][chrom][0] + pos[order[-1]][chrom][1]) / 2
    ax.text(xm, ysb - offset * 1.2, chrom_label, ha="center", va="top", fontsize=LAB)

    # pies (plot_TE_frac2.plot_pie_charts colours/orientation), right of each bar
    pies = pies_table(ptf, fai_len, verbose)
    cmap = plt.get_cmap("tab10")
    cmap_f = {k: cmap(i) for i, k in enumerate(FEAT_KEYS)}
    cmap_f["Non-TE"] = "gray"
    d = 5.6  # pie diameter (mm)
    to_fig = ax.transData + fig.transFigure.inverted()
    for a, (_, r) in zip(order, pies.iterrows()):
        _, xe, y, _ = pos[a][chrom]
        fx, fy = to_fig.transform((xe, y))
        pax = fig.add_axes([fx + 1.3 / FW, fy - d / 2 / FH, d / FW, d / FH])
        labels, sizes = [], []
        for k in FEAT_KEYS:
            p = r[f"{k.replace(' ', '_')}_pct"]
            if p > 0:
                labels.append(k)
                sizes.append(p)
        if 100 - sum(sizes) > 0:
            labels.append("Non-TE")
            sizes.append(100 - sum(sizes))
        pax.pie(sizes, startangle=90, counterclock=False, colors=[cmap_f[k] for k in labels])
        pax.set_aspect("equal")
        pax.axis("off")
    handles = [Patch(facecolor=cmap_f[k], edgecolor="none", label=k) for k in FEAT_KEYS + ["Non-TE"]]
    fig.legend(handles=handles, title="Features", loc="upper left", fontsize=LEG, title_fontsize=LEG,
               bbox_to_anchor=(41.5 / FW, 1 - 26.0 / FH), handlelength=1.0, handleheight=1.0,
               handletextpad=0.4, labelspacing=0.3, borderaxespad=0, frameon=False, alignment="left")
    lens = pd.DataFrame([dict(generation=int(a.replace("gen", "")), chromosome=s, length_bp=L,
                              plotted_in_ribbon=(s == chrom))
                         for a in order for s, L in sorted(chrom_info[a].items())])
    return pies, lens, len(polys)


# ----------------------------------------------------------------------------- panels B, C, I
def density_panel(fig, ld, rect, files, highlight, xmax_manual, show_xmax, yticks, lw, lw_hi, verbose, tag):
    """ltr_dens.create_density_plot semantics: per-curve Gaussian KDE (Scott bw), clip=(0, x_max),
    plasma gradient 0.15-0.85, highlighted generation = red dashed."""
    need(*files)
    data, hl = ld.load_data(files, "K2P", highlight)
    if xmax_manual is not None:
        x_max = xmax_manual
    else:
        m = data["distance"].max()
        x_max = m * 1.1 if m > 0 else 1
    ax = ax_mm(fig, *rect)
    ax.set_xlim(0, x_max)
    sources = sorted(data["source"].unique(),
                     key=lambda s: int(s.replace("Generation", "").strip()) if "Generation" in s else float("inf"))
    colors = {s: plt.cm.plasma(v) for s, v in zip(sources, np.linspace(0.15, 0.85, len(sources)))}
    import seaborn as sns
    for s in sources:
        sub = data[data["source"] == s]["distance"]
        if s == hl:
            sns.kdeplot(sub, ax=ax, color="red", linestyle="--", linewidth=lw_hi, fill=False,
                        clip=(0, x_max), warn_singular=False)
        else:
            sns.kdeplot(sub, ax=ax, color=colors[s], linewidth=lw, fill=False, clip=(0, x_max),
                        warn_singular=False)
        if verbose:
            x = sub.dropna()
            print(f"  {tag} {s}: n={len(x)} scott_bw={x.std(ddof=1) * len(x) ** -0.2:.4f} "
                  f"clip=(0,{x_max:.4f})")
    ax.set_xlim(0, show_xmax if show_xmax is not None else x_max)
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.xaxis.set_major_formatter(FuncFormatter(gfmt))
    ax.set_xlabel("LTR divergence (K2P)", fontsize=LAB, labelpad=1.5)
    ax.set_ylabel("Density", fontsize=LAB, labelpad=2)
    box(ax)
    handles = [Patch(facecolor="red" if s == hl else colors[s], edgecolor="black", linewidth=0.4,
                     label=s.replace("Generation", "").strip()) for s in sources]
    ax.legend(handles=handles, title="Generation", loc="upper right", fontsize=LEG, title_fontsize=LEG,
              handlelength=1.4, handleheight=0.8, handletextpad=0.5, labelspacing=0.25,
              borderaxespad=0.4, frameon=False)
    long = pd.DataFrame({"generation": data["source"].str.replace("Generation", "").str.strip().astype(int),
                         "K2P_divergence": data["distance"]})
    return long.sort_values("generation", kind="stable").reset_index(drop=True), x_max


# ----------------------------------------------------------------------------- panels F, J
def indel_table(report):
    need(report)
    t = pd.read_csv(report, sep="\t")
    t["TE_insertions"] = t["TE_inserts(nest/nonnest)"].str.extract(r"^(\d+)")[0].astype(int)
    return t


def indel_panel(fig, rect, t, step, xlim, xticks, yticks, rticks, rlabels, ylabel, lw):
    """Rscript2.R --collapse: counts (thousands, solid) and per-bp rates (dashed) on a secondary axis.
    Same axis mapping as ggplot (primary 0..max count, secondary = rate * max count / max rate);
    the rate axis is shown per generation (report rate per step / step), as relabelled on the slide."""
    g = t["Generation"] / 1e6
    ins_k, del_k = t["TE_insertions"] / 1e3, t["Actual_TE_deletions"] / 1e3
    pmax = max(ins_k.max(), del_k.max())
    rmax = max(t["insertion_rate"].max(), t["deletion_rate"].max())
    sec = pmax / rmax
    ax = ax_mm(fig, *rect)
    ax.plot(g, del_k, color=DEL, lw=lw, solid_joinstyle="round", zorder=3)
    ax.plot(g, ins_k, color=INS, lw=lw, solid_joinstyle="round", zorder=3)
    dash = (0, (3.0, 2.2))
    ax.plot(g, t["deletion_rate"] * sec, color=DEL, lw=lw, ls=dash, zorder=3)
    ax.plot(g, t["insertion_rate"] * sec, color=INS, lw=lw, ls=dash, zorder=3)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, pmax)
    ax.xaxis.set_major_locator(FixedLocator(xticks))
    ax.xaxis.set_minor_locator(FixedLocator([(a + b) / 2 for a, b in zip(xticks[:-1], xticks[1:])]))
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.yaxis.set_minor_locator(FixedLocator([(a + b) / 2 for a, b in zip(yticks[:-1], yticks[1:])]))
    ax.xaxis.set_major_formatter(FuncFormatter(gfmt))
    ax.yaxis.set_major_formatter(FuncFormatter(gfmt))
    ax.grid(True, which="major", color="#EBEBEB", lw=0.45)
    ax.grid(True, which="minor", color="#EBEBEB", lw=0.25)
    ax.set_axisbelow(True)
    ax.tick_params(which="both", length=0, pad=1.5)
    box(ax, 0.6)
    sec_ax = ax.secondary_yaxis("right", functions=(lambda y: y / sec / step, lambda r: r * sec * step))
    sec_ax.yaxis.set_major_locator(FixedLocator(rticks))
    sec_ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: rlabels.get(round(v, 20), "")))
    sec_ax.yaxis.set_minor_formatter(NullFormatter())
    sec_ax.tick_params(which="both", length=0, pad=1.5, labelsize=TXT)
    sec_ax.spines["right"].set_visible(False)
    sec_ax.set_ylabel("Rate (bp$^{-1}$\ngeneration$^{-1}$)", fontsize=RATE_LAB, rotation=270, va="bottom",
                      labelpad=2.5, linespacing=1.1)
    ax.set_xlabel("Generations (millions)", fontsize=LAB, labelpad=1.5)
    ax.set_ylabel(ylabel, fontsize=LAB, labelpad=2, linespacing=1.1)
    h = [Line2D([], [], color=INS, lw=lw, label="Insertion count"),
         Line2D([], [], color=DEL, lw=lw, label="Deletion count"),
         Line2D([], [], color=INS, lw=lw, ls=dash, label="Insertion rate"),
         Line2D([], [], color=DEL, lw=lw, ls=dash, label="Deletion rate")]
    ax.legend(handles=h, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=LEG,
              handlelength=2.0, handletextpad=0.4, columnspacing=1.0, labelspacing=0.2,
              borderaxespad=0.3, borderpad=0, frameon=False)
    out = pd.DataFrame({
        "generation": t["Generation"],
        "TE_insertions_count": t["TE_insertions"],
        "TE_deletions_count": t["Actual_TE_deletions"],
        "insertion_rate_per_bp_per_generation": (t["insertion_rate"] / step).map(sig4),
        "deletion_rate_per_bp_per_generation": (t["deletion_rate"] / step).map(sig4),
        "insertion_rate_per_bp_per_step_as_in_report": t["insertion_rate"],
        "deletion_rate_per_bp_per_step_as_in_report": t["deletion_rate"],
        "genome_size_bp": t["genome_size"],
    })
    return out


# ----------------------------------------------------------------------------- panels G, H
def panel_G(fig, rect):
    t = indel_table(f"{VM}/pipeline.report.revised.tsv")
    df = pd.DataFrame({"generation": t["Generation"].astype(int), "genome_size_bp": t["genome_size"].astype(int)})
    df = df.sort_values("generation").reset_index(drop=True)
    ax = ax_mm(fig, *rect)
    ax.set_xticks(df["generation"], minor=True)  # genome_plot.py: a (dashed) grid line at every generation
    ax.grid(True, which="both", axis="x", linestyle="--", alpha=0.6, lw=0.25)
    ax.grid(True, which="major", axis="y", linestyle="--", alpha=0.6, lw=0.25)
    ax.plot(df["generation"], df["genome_size_bp"] / 1e6, marker="o", ls="-", color="black", ms=2.2, lw=0.5,
            mew=0, clip_on=False)
    ax.set_xlim(0, 1.0e6)
    ax.xaxis.set_major_locator(FixedLocator([0, 1_000_000]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{int(v):,}"))
    ax.tick_params(axis="x", which="minor", length=0)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.set_ylim(122.5, 141.5)
    ax.set_xlabel("Generations", fontsize=LAB, labelpad=1.5)
    ax.set_ylabel("Genome Size (Mb)", fontsize=LAB, labelpad=2)
    box(ax)
    df["genome_size_Mb"] = df["genome_size_bp"] / 1e6
    return df


def panel_H(fig, rect_i, rect_f, leg_xy):
    def read(fn):
        need(fn)
        d = pd.read_csv(fn, sep="\t")
        long = d.melt(id_vars=["TE_class/TE_superfamily"], var_name="Generation", value_name="Count")
        long["Generation"] = long["Generation"].str.extract(r"gen(\d+)_final_Count")[0].astype(int)
        return long
    intact, frag = read(f"{VM}/stat_intact.tsv"), read(f"{VM}/stat_frag.tsv")
    classes = sorted(set(intact["TE_class/TE_superfamily"]) | set(frag["TE_class/TE_superfamily"]))
    import seaborn as sns
    pal = sns.color_palette("tab20", len(classes))  # plot_superfamily_count.py colours
    cmap = {c: pal[i] for i, c in enumerate(classes)}
    out, axes = [], []
    for long, rect, title, yt in [(intact, rect_i, "Intact TEs", [0, 2, 4, 6, 8, 10]),
                                  (frag, rect_f, "Fragmented TEs", [0, 1, 2])]:
        pv = long.pivot(index="Generation", columns="TE_class/TE_superfamily", values="Count").fillna(0).sort_index()
        ax = ax_mm(fig, *rect)
        axes.append(ax)
        bottom = np.zeros(len(pv))
        for c in pv.columns:  # stacked bars, width 0.8 x 10,000-generation step (pandas bar width=0.8)
            ax.bar(pv.index, pv[c] / 1e3, bottom=bottom, width=8000, color=cmap[c], linewidth=0, label=c)
            bottom += pv[c].values / 1e3
        ax.set_xlim(0, 1.006e6)
        ax.set_ylim(0, bottom.max() * 1.05)
        ax.xaxis.set_major_locator(FixedLocator([0, 1_000_000]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{int(v):,}"))
        ax.yaxis.set_major_locator(FixedLocator(yt))
        ax.set_title(title, fontsize=LAB, loc="left", pad=1.5)
        ax.set_xlabel("Generations", fontsize=LAB, labelpad=1.5)
        ax.set_ylabel("TE count (thousands)", fontsize=LAB, labelpad=2)
        box(ax)
        long2 = long.rename(columns={"TE_class/TE_superfamily": "superfamily", "Generation": "generation",
                                     "Count": "count"})
        long2.insert(0, "TE_state", title.split()[0].lower())
        out.append(long2)
    fig.align_ylabels(axes)
    handles = [Patch(facecolor=cmap[c], edgecolor="none", label=c) for c in classes]
    fig.legend(handles=handles, title="TE class/superfamily", loc="upper left", fontsize=LEG,
               title_fontsize=LEG, bbox_to_anchor=(leg_xy[0] / FW, 1 - leg_xy[1] / FH), handlelength=1.0,
               handleheight=1.0, handletextpad=0.4, labelspacing=0.22, borderaxespad=0, frameon=False,
               alignment="left")
    return pd.concat(out, ignore_index=True).sort_values(["TE_state", "generation", "superfamily"],
                                                         ascending=[False, True, True]).reset_index(drop=True)


# ----------------------------------------------------------------------------- panels D, E
def newick_tips(nwk):
    return re.findall(r"[(,]([^(),:;]+):", nwk)


def render_trees(build, rscript, verbose):
    need(TREE_D, TREE_E, rscript, HERE / "make_fig3_trees.R")
    out = [str(build / "fig3D_tree.pdf"), str(build / "fig3E_tree.pdf")]
    cmd = [rscript, str(HERE / "make_fig3_trees.R"), TREE_D, out[0], TREE_E, out[1],
           str(TREE_MM_PER_UNIT), str(TREE_RADIUS), str(TREE_LW_PT)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: tree rendering failed:\n{r.stderr}")
    if verbose:
        for ln in r.stdout.strip().splitlines():
            print("  " + ln)
    return out


def tree_tables():
    tips, chunks = [], []
    for panel, name, fn in [("D", "contraction", TREE_D), ("E", "expansion", TREE_E)]:
        nwk = Path(fn).read_text().strip()
        for lab in newick_tips(nwk):
            tag = lab.split("#")[0]
            cat = "unique" if tag.startswith("uniq") else "shared" if tag.startswith("shared") else "ungrouped"
            tips.append(dict(panel=panel, tree=name, tip_label=lab, category=cat))
        for i in range(0, len(nwk), 32000):
            chunks.append(dict(panel=panel, tree=name, chunk_index=i // 32000 + 1, newick_chunk=nwk[i:i + 32000]))
    return pd.DataFrame(tips), pd.DataFrame(chunks)


def place_trees(pdf_in, pdf_out, tree_pdfs, anchors):
    """Insert the text-free ggtree PDFs as vector Form XObjects at 1:1 scale; anchors = (x, y) mm of bbox TL."""
    import fitz
    doc = fitz.open(pdf_in)
    page = doc[0]
    k = 72 / 25.4
    boxes = []
    for fn, (x, y) in zip(tree_pdfs, anchors):
        src = fitz.open(fn)
        bb = fitz.Rect()
        for d in src[0].get_drawings():
            bb |= d["rect"]
        bb = (bb + (-1, -1, 1, 1)) & src[0].rect  # stroke-width margin
        tgt = fitz.Rect(x * k, y * k, x * k + bb.width, y * k + bb.height)
        page.show_pdf_page(tgt, src, 0, clip=bb, keep_proportion=True, overlay=True)
        boxes.append((tgt, bb))
    doc.save(pdf_out, garbage=4, deflate=True)
    return boxes


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default=str(ROOT / "fig3" / "fig3.pdf"), help="output PDF (default: %(default)s)")
    ap.add_argument("--build-dir", default=str(ROOT / "build" / "fig3"),
                    help="intermediates: tree PDFs, pre-tree PDF, PNG preview (default: %(default)s)")
    ap.add_argument("--rscript", default=RSCRIPT, help="Rscript with ggtree/ape (default: synLTR env)")
    ap.add_argument("--ribbon-chrom", default="chr2",
                    help="chromosome drawn in panel A (published crop = chr2; default: %(default)s)")
    ap.add_argument("--ribbon-label", default="Chr1",
                    help="label under the panel A ribbon (authors label the plotted chromosome Chr1; "
                         "Source Data renames chr1<->chr2 to match; default: %(default)s)")
    ap.add_argument("--published-labels", action="store_true",
                    help="use the slide's exact wording where it is inaccurate ('Chr1' for the chr2 ribbon; "
                         "'TE InDel (Counts)' for counts in thousands)")
    ap.add_argument("--no-source-data", action="store_true", help="skip writing Source Data blocks")
    ap.add_argument("--png-dpi", type=int, default=200, help="preview PNG resolution (0 = none; default 200)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-panel statistics and sanity checks")
    a = ap.parse_args()
    v = a.verbose
    build = Path(a.build_dir)
    build.mkdir(parents=True, exist_ok=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    print(f"make_fig3: building {a.out}")

    rip = load_script(f"{RIP}/synLTR/bin/riparian3.py", "riparian3")
    ptf = load_script(f"{CE}/plot_TE_frac2.py", "plot_TE_frac2")
    ld = load_script(f"{PRINTE}/bin/ltr_dens.py", "ltr_dens")

    ncomms_style.use()
    plt.rcParams.update({"axes.spines.top": True, "axes.spines.right": True, "xtick.labelsize": TXT,
                         "ytick.labelsize": TXT, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
                         "xtick.major.size": 2.0, "ytick.major.size": 2.0, "axes.labelpad": 2.0})
    fig = plt.figure(figsize=(FW / 25.4, FH / 25.4))

    # ---- top row
    chrom_label = "Chr1" if a.published_labels else a.ribbon_label
    count_label = "TE InDel\n(Counts)" if a.published_labels else "TE InDels\n(thousands)"
    letter(fig, 0, 0, "A")
    pies, lens, n_rib = panel_A(fig, rip, ptf, a.ribbon_chrom, chrom_label, v)
    print("  panel A done")

    letter(fig, 60.5, 0, "B")
    fB = sorted(glob.glob(f"{CE}/temp/gen[0-9]*_LTR.tsv")) + [f"{CE}/temp/burnin_LTR.tsv"]
    dB, _ = density_panel(fig, ld, (69.5, 3.0, 47.5, 30.5), fB, f"{CE}/temp/gen5000000_final_LTR.tsv", 0.2, None,
                          [0, 10, 20, 30], 0.8, 1.2, v, "B")
    letter(fig, 60.5, 42.0, "C")
    fC = sorted(glob.glob(f"{CE}/expansion/temp/gen[0-9]*_LTR.tsv"))
    dC, _ = density_panel(fig, ld, (69.5, 45.0, 47.5, 30.5), fC, f"{CE}/expansion/temp/gen5000000_final_LTR.tsv",
                          0.2, None, [0, 20, 40, 60, 80], 0.8, 1.2, v, "C")
    print("  panels B, C done")

    letter(fig, 121.5, 0, "D")
    letter(fig, 121.5, 20.0, "E")
    k = TREE_MM_PER_UNIT
    sb = ax_mm(fig, 157.0, 30.0, 0.5 * k, 3.0)  # tree scale bar (0.5 substitutions per site)
    sb.set_xlim(0, 0.5 * k)
    sb.set_ylim(0, 1)
    sb.plot([0, 0.5 * k], [0.25, 0.25], color="black", lw=0.8, solid_capstyle="butt", clip_on=False)
    sb.text(0.25 * k, 0.45, "0.5", ha="center", va="bottom", fontsize=TXT)
    sb.axis("off")
    fig.legend(handles=[Patch(facecolor="#1f78b4", label="Shared"), Patch(facecolor="#e31a1c", label="Unique")],
               loc="upper left", bbox_to_anchor=(158.0 / FW, 1 - 3.0 / FH), fontsize=LEG, handlelength=1.0,
               handleheight=1.0, handletextpad=0.4, labelspacing=0.3, borderaxespad=0, frameon=False)

    letter(fig, 121.5, 45.5, "F")
    tF = indel_table(f"{CE}/expansion/pipeline_all.report.revised.cleaned")
    sF = indel_panel(fig, (130.5, 55.0, 31.5, 20.5), tF, 1e6, (0, 10), [0, 2.5, 5, 7.5, 10], [0, 5, 10, 15],
                     [1e-11, 3e-11, 5e-11], {1e-11: "1e-11", 3e-11: "3e-11", 5e-11: "5e-11"}, count_label, 0.8)
    print("  panel F done")

    # ---- bottom row
    letter(fig, 0, 86.0, "G")
    sG = panel_G(fig, (11.0, 90.0, 44.0, 51.0))
    letter(fig, 60.5, 86.0, "H")
    sH = panel_H(fig, (69.5, 91.0, 24.5, 19.0), (69.5, 123.0, 24.5, 18.0), (99.0, 92.5))
    print("  panels G, H done")
    letter(fig, 121.5, 86.0, "I")
    fI = sorted(glob.glob(f"{VM}/gen[0-9]*_LTR.tsv"))
    dI, xmaxI = density_panel(fig, ld, (130.5, 89.5, 31.5, 19.5), fI, None, None, 0.06,
                              [0, 10, 20, 30, 40, 50, 60], 0.8, 1.2, v, "I")
    letter(fig, 121.5, 117.5, "J")
    tJ = indel_table(f"{VM}/pipeline.report.revised.tsv")
    sJ = indel_panel(fig, (130.5, 126.5, 31.5, 17.5), tJ, 1e4, (0, 1), [0, 0.25, 0.5, 0.75, 1], [0, 0.1, 0.2],
                     [5e-11, 1e-10, 1.5e-10], {5e-11: "5e-11", 1e-10: "1e-10", 1.5e-10: "1.5e-10"}, count_label,
                     0.8)
    print("  panels I, J done")

    pre = build / "fig3_no_trees.pdf"
    fig.savefig(pre, dpi=600)
    plt.close(fig)

    tree_pdfs = render_trees(build, a.rscript, v)
    boxes = place_trees(str(pre), a.out, tree_pdfs, [(129.0, 1.0), (127.5, 21.5)])
    for (tgt, bb), p in zip(boxes, "DE"):
        log(f"  {p} tree placed at {tgt.x0 / 72 * 25.4:.1f},{tgt.y0 / 72 * 25.4:.1f} mm, "
            f"{bb.width / 72 * 25.4:.1f} x {bb.height / 72 * 25.4:.1f} mm (1:1, {TREE_MM_PER_UNIT} mm per unit)", v)
    print("  panels D, E placed")

    if a.png_dpi:
        import fitz
        png = build / f"fig3_{a.png_dpi}dpi.png"
        fitz.open(a.out)[0].get_pixmap(dpi=a.png_dpi).save(str(png))
        log(f"  preview: {png}", v)

    if not a.no_source_data:
        write_source_data(pies, lens, n_rib, dB, dC, dI, xmaxI, sF, sG, sH, sJ, a.ribbon_chrom, chrom_label)
    print(f"make_fig3: done -> {a.out}")


def ncount(d):
    return "; ".join(f"{g:,}: n = {n:,}" for g, n in d.groupby("generation").size().items())


def write_source_data(pies, lens, n_rib, dB, dC, dI, xmaxI, sF, sG, sH, sJ, chrom, chrom_label):
    sd = SourceData("fig3")
    sd.add("Fig. 3A (pies)",
           "Genome composition per generation of the 10 M-generation contraction (0-5 M) then expansion "
           "(5-10 M) simulation: bp and % of genome covered by intact TEs, solo LTRs and fragmented TEs "
           "(from each generation's BED annotation); Non-TE = 100% minus these; genome_bp = total sequence "
           "length (bp) of that generation's genome FASTA. One simulated genome per generation "
           "(n = 11 generations).", pies)
    sd.add("Fig. 3A (chromosome size)",
           f"Chromosome lengths (bp) of the simulated genome (4 chromosomes) at each generation (0-10 M, step "
           f"1 M). The ribbon plot shows {chrom} (labelled '{chrom_label}'); its {n_rib:,} grey ribbons are "
           "syntenic gene-anchor blocks between adjacent generations (all.recip.anchors.coords), "
           "merged within 10 kb (the 24.9 MB anchor file is not reproduced here).", lens)
    sd.add("Fig. 3B",
           "Per-element LTR divergence (K2P distance between the two LTRs of each intact LTR-RT) during genome "
           "contraction, generations 0-5 M; curves are Gaussian KDEs (Scott bandwidth, clipped to 0-0.2, each "
           "curve normalised separately). " + ncount(dB) + ".", dB)
    sd.add("Fig. 3C",
           "Per-element LTR divergence (K2P) during genome expansion, generations 5-10 M (5 M identical to "
           "panel B); Gaussian KDEs (Scott bandwidth, clipped to 0-0.2, per-curve normalisation). "
           + ncount(dC) + ".", dC)
    tips, chunks = tree_tables()
    nd = tips[tips.panel == "D"].category.value_counts().to_dict()
    ne = tips[tips.panel == "E"].category.value_counts().to_dict()
    sd.add("Fig. 3D-E (tips)",
           "Tip labels of the Tekay RT-domain phylogenies (VeryFastTree -lg -gamma on TEsorter RT protein "
           "alignments) and their category: shared = copies from the burn-in (shared) library, unique = insertions "
           "made after the burn-in. Trees were built from two replicate runs of 4,000,000 generations each "
           "(prinTE --generation_end 4000000; decrease: -F 2e-11,1e-11; increase: -F 2e-10,1e-10; increase "
           f"subsampled to 20% with seqtk). D (decrease): n = {sum(nd.values())} tips ({nd.get('shared', 0)} "
           f"shared, {nd.get('unique', 0)} unique); E (increase): n = {sum(ne.values())} tips "
           f"({ne.get('shared', 0)} shared, {ne.get('unique', 0)} unique).", tips)
    sd.add("Fig. 3D-E (newick)",
           "Newick strings of the two trees (branch lengths = substitutions per site; internal labels = SH-like "
           "support), split into consecutive chunks of <= 32,000 characters (concatenate by chunk_index).", chunks)
    sd.add("Fig. 3F",
           "TE insertions and deletions per 1 M-generation step (counts; plotted in thousands) and rates per bp "
           "per generation (per-step rate / 1e6) across the 10 M generations. The 5,000,000 row appears twice: "
           "the second copy repeats the 6,000,000 values so the plotted lines step at 5 M, as in the original "
           "plotting table.", sF.assign(note=np.where(sF.duplicated("generation", keep="first"),
                                                     "repeat of 6,000,000 values (plotted step at 5 M)", "")))
    sd.add("Fig. 3G",
           "Genome size (bp and Mb) of the A. thaliana replica forward-simulated for 1 M generations with the "
           "variable-rate method, every 10,000 generations (n = 100 time points, one simulated genome).", sG)
    sd.add("Fig. 3H",
           "Number of intact and fragmented TEs per class/superfamily (counts; plotted in thousands) every "
           "10,000 generations of the A. thaliana replica simulation (13 superfamilies x 100 time points).", sH)
    sd.add("Fig. 3I",
           "Per-element LTR divergence (K2P) of intact LTR-RTs in the A. thaliana replica simulation; Gaussian "
           f"KDEs (Scott bandwidth, clipped to 0-{xmaxI:.4f}, per-curve normalisation; x axis shown to 0.06). "
           + ncount(dI) + ".", dI)
    sd.add("Fig. 3J",
           "TE insertions and deletions per 10,000-generation step (counts; plotted in thousands) and rates per "
           "bp per generation (per-step rate / 1e4) in the A. thaliana replica simulation (n = 100 steps).", sJ)
    sd.save()
    print(f"  source data: {len(sd.blocks)} blocks -> {sd.dir}")


if __name__ == "__main__":
    main()
