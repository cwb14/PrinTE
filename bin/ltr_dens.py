#!/usr/bin/env python3
import glob
import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# Expected columns for the new headerless format:
# <LTR-RT_name>  <LTR_LEN>  <ALN_LEN>  <substitions> <transitions>  <transversions>
# <raw_d> <raw_T> <JC69_d> <JC69_T> <K2P_d>  <K2P_T>
NEW_FORMAT_COLS = [
    "LTR_RT_name", "LTR_LEN", "ALN_LEN", "substitutions", "transitions", "transversions",
    "raw_d", "raw_T", "JC69_d", "JC69_T", "K2P_d", "K2P_T"
]

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot density distributions from gen[int]_LTR.tsv files.")
    parser.add_argument("--model", required=True, choices=["raw", "K2P", "JC69"],
                        help="Model to use (choose from raw, K2P, JC69).")
    parser.add_argument("--miu", required=True, type=float,
                        help="Mutation rate (used in the plot title and for time conversion).")
    parser.add_argument("--output", default="density_plot.pdf",
                        help="Output PDF file name.")
    parser.add_argument("--gradient", action="store_true",
                        help="If provided, color lines using a distinctive gradient.")
    parser.add_argument("--xmax", type=float, default=None,
                        help="Manually set the upper x-axis limit.")
    return parser.parse_args()

def find_files():
    pattern = "gen[0-9]*_LTR.tsv"
    files = glob.glob(pattern)
    burnin_file = "burnin_LTR.tsv"
    if os.path.exists(burnin_file):
        files.append(burnin_file)
    if not files:
        raise FileNotFoundError("No files found matching pattern 'gen[int]_LTR.tsv' or 'burnin_LTR.tsv'")
    # Sort numerically by generation if possible, keeping non-matching (e.g., burnin) stable
    def gen_key(f):
        b = os.path.basename(f)
        m = re.search(r'gen(\d+)', b)
        return (0, int(m.group(1))) if m else (1, b)
    return sorted(files, key=gen_key)

def extract_source_label(filename):
    base = os.path.basename(filename)
    if base == "burnin_LTR.tsv":
        return "Generation 0"
    m = re.search(r'gen(\d+)', base)
    return f"Generation {m.group(1)}" if m else base

def read_ltr_tsv(fpath):
    """
    Try reading with a header; if the required columns aren't there,
    re-read as headerless using NEW_FORMAT_COLS.
    """
    # First attempt: assume there might be a header
    try:
        df_try = pd.read_csv(fpath, sep="\t")
        if set(["raw_d","raw_T","JC69_d","JC69_T","K2P_d","K2P_T"]).issubset(df_try.columns):
            return df_try
    except Exception:
        pass

    # Second attempt: headerless new format
    df = pd.read_csv(fpath, sep="\t", header=None, names=NEW_FORMAT_COLS, engine="python")
    # Ensure numeric columns are numeric
    num_cols = [c for c in NEW_FORMAT_COLS if c.endswith("_d") or c.endswith("_T") or c in ["LTR_LEN","ALN_LEN","substitutions","transitions","transversions"]]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def load_data(files, model):
    model = model.strip()
    d_col = f"{model}_d"
    t_col = f"{model}_T"

    data_list = []
    for f in files:
        try:
            df = read_ltr_tsv(f)
            # Guard: ensure required cols exist
            if d_col not in df.columns or t_col not in df.columns:
                print(f"Skipping {f}: missing required columns {d_col} or {t_col}.")
                continue

            # Some rows may be empty/zero; keep them but coerce to float
            df = df[[d_col, t_col]].copy()
            df.rename(columns={d_col: "distance", t_col: "time"}, inplace=True)

            # Distance can be NaN if parsing failed; drop NaN-only distances
            before = len(df)
            df = df[pd.to_numeric(df["distance"], errors="coerce").notna()]
            after = len(df)
            if after == 0:
                print(f"Skipping {f}: no valid numeric distances after parsing.")
                continue
            if after < before:
                print(f"Note: {f} had {before-after} rows with non-numeric distance and were dropped.")

            df["source"] = extract_source_label(f)
            data_list.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")

    if not data_list:
        raise ValueError("No valid data loaded from files.")
    return pd.concat(data_list, ignore_index=True)

def create_density_plot(data, miu, model, output, use_gradient=False, manual_xmax=None):
    plt.figure(figsize=(10, 6))
    ax = plt.gca()

    # X limit
    if manual_xmax is not None and manual_xmax > 0:
        x_max = manual_xmax
    else:
        max_distance = data['distance'].max()
        x_max = max(max_distance * 1.1, 1e-6)  # avoid zero-width axis
    ax.set_xlim(0, x_max)

    sources = list(data['source'].unique())

    # Sort sources numerically by generation; non-numeric go to the end
    def gen_val(s):
        try:
            return int(s.replace("Generation", "").strip())
        except Exception:
            return float('inf')
    sorted_sources = sorted(sources, key=gen_val)

    # Colors
    if use_gradient:
        cmap = plt.cm.plasma
        colors = {s: cmap(val) for s, val in zip(sorted_sources, np.linspace(0.15, 0.85, len(sorted_sources)))}
    else:
        palette = sns.color_palette("husl", len(sorted_sources))
        colors = {s: palette[i] for i, s in enumerate(sorted_sources)}

    # Plot KDE per source
    for source in sorted_sources:
        subset = data.loc[data['source'] == source, "distance"]
        # Skip if all zeros (singular KDE); plot will error or be flat—guard it.
        if subset.nunique(dropna=True) <= 1:
            # Draw a tiny vertical line to indicate a degenerate distribution
            x = float(subset.iloc[0]) if len(subset) else 0.0
            if np.isfinite(x):
                ax.vlines(x, 0, 1, colors=colors[source], linestyles="dotted", label=None)
            continue
        sns.kdeplot(subset, label=source, color=colors[source],
                    fill=False, clip=(0, x_max), warn_singular=False)

    # Axis labels
    plt.xlabel(f"{model} genetic distance")
    plt.ylabel("Density")

    # Top axis: time in MYA using miu
    def forward(x):
        return (x / miu) / 1e6
    def inverse(x):
        return x * miu * 1e6
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel("Time (MYA)")

    plt.title(f"Density Plot for {model} distances\nMutation Rate (miu) = {miu}")

    # Legend (sorted by generation, label with numeric gen if possible)
    legend_items = []
    for s in sorted_sources:
        try:
            gen_label = str(int(s.replace("Generation", "").strip()))
        except ValueError:
            gen_label = s
        legend_items.append(Patch(facecolor=colors[s], edgecolor='black', label=gen_label))
    ax.legend(handles=legend_items, title="Generation", frameon=True)

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print(f"Plot saved to {output}")

def main():
    args = parse_arguments()
    files = find_files()
    data = load_data(files, args.model)
    create_density_plot(data, args.miu, args.model, args.output,
                        use_gradient=args.gradient, manual_xmax=args.xmax)

if __name__ == "__main__":
    main()
