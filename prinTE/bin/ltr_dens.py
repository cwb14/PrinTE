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

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot density distributions from gen[int]_LTR.tsv files.")
    parser.add_argument("--model", required=True, choices=["raw", "K2P", "JC69"],
                        help="Model to use (choose from raw, K2P, JC69).")
    parser.add_argument("--miu", required=True, type=float,
                        help="Mutation rate (used in the plot title and for time conversion).")
    parser.add_argument("--output", default="density_plot.pdf",
                        help="Output PDF file name.")
    parser.add_argument("--gradiant", action="store_true",
                        help="If provided, color lines using a gradient (lower generations darker, higher lighter).")
    parser.add_argument("--xmax", type=float, default=None,
                        help="Manually set the upper x-axis limit.")
    return parser.parse_args()

def find_files():
    # Find all files matching the pattern gen[int]_LTR.tsv.
    pattern = "gen[0-9]*_LTR.tsv"
    files = glob.glob(pattern)
    
    # Check for burnin_LTR.tsv (representing generation 0)
    burnin_file = "burnin_LTR.tsv"
    if os.path.exists(burnin_file):
        files.append(burnin_file)
    
    if not files:
        raise FileNotFoundError("No files found matching pattern 'gen[int]_LTR.tsv' or 'burnin_LTR.tsv'")
    return files

def extract_source_label(filename):
    base = os.path.basename(filename)
    if base == "burnin_LTR.tsv":
        return "Generation 0"
    m = re.search(r'gen(\d+)', base)
    return f"Generation {m.group(1)}" if m else base

def load_data(files, model):
    data_list = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t')
            d_col = f"{model}_d"
            t_col = f"{model}_T"
            if d_col not in df.columns or t_col not in df.columns:
                print(f"Skipping {f}: missing required columns {d_col} or {t_col}.")
                continue
            source_label = extract_source_label(f)
            df = df[[d_col, t_col]].rename(columns={d_col: "distance", t_col: "time"})
            df['source'] = source_label
            data_list.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")
    if not data_list:
        raise ValueError("No valid data loaded from files.")
    return pd.concat(data_list, ignore_index=True)

def create_density_plot(data, miu, model, output, use_gradient=False, manual_xmax=None):
    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    
    # Determine x-axis limit: use manual limit if provided, else 10% above maximum.
    if manual_xmax is not None:
        x_max = manual_xmax
    else:
        max_distance = data['distance'].max()
        x_max = max_distance * 1.1 if max_distance > 0 else 1
    ax.set_xlim(0, x_max)
    
    # Get unique sources.
    sources = data['source'].unique()
    
    # Create a mapping from source to generation number for sorting.
    source_to_gen = {}
    for s in sources:
        try:
            # Remove "Generation" and extra spaces, then convert to int.
            gen = int(s.replace("Generation", "").strip())
        except ValueError:
            gen = float('inf')
        source_to_gen[s] = gen
    
    if use_gradient:
        # Sort sources by generation number.
        sorted_sources = sorted(sources, key=lambda s: source_to_gen[s])
        # Create a gradient from dark to light using a colormap.
        # Here we use np.linspace to get values between 0.3 (dark) and 0.7 (light).
        cmap = plt.cm.viridis
        colors = {s: cmap(val) for s, val in zip(sorted_sources, np.linspace(0.3, 0.7, len(sorted_sources)))}
    else:
        # Use the standard husl palette.
        palette = sns.color_palette("husl", len(sources))
        colors = {s: palette[i] for i, s in enumerate(sources)}
    
    # Plot density curves.
    for source in sources:
        subset = data[data['source'] == source]
        sns.kdeplot(subset['distance'], label=source, color=colors[source],
                    fill=False, clip=(0, x_max), warn_singular=False)
    
    plt.xlabel(f"{model} genetic distance")
    plt.ylabel("Density")
    
    # Create secondary x-axis converting genetic distance to time (MYA).
    def forward(x):
        return (x / miu) / 1e6
    def inverse(x):
        return x * miu * 1e6
    
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel("Time (MYA)")
    
    plt.title(f"Density Plot for {model} distances\nMutation Rate (miu) = {miu}")
    
    # Sort legend entries numerically by generation.
    # Create sorted legend entries with colored squares
    legend_items = []
    for s in sorted(sources, key=lambda s: source_to_gen[s]):
        try:
            gen_label = str(int(s.replace("Generation", "").strip()))
        except ValueError:
            gen_label = s
        legend_items.append(Patch(facecolor=colors[s], edgecolor='black', label=gen_label))

    ax.legend(handles=legend_items, title="Generation")    
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print(f"Plot saved to {output}")

def main():
    args = parse_arguments()
    files = find_files()
    data = load_data(files, args.model)
    create_density_plot(data, args.miu, args.model, args.output,
                        use_gradient=args.gradiant, manual_xmax=args.xmax)

if __name__ == "__main__":
    main()
