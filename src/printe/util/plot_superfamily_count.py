#!/usr/bin/env python3
"""
plot_superfamily_count.py

Usage:
    python -m printe.util.plot_superfamily_count

Description:
    - Reads 'stat_intact.tsv' and 'stat_frag.tsv'.
    - Extracts TE counts per generation.
    - Generates stacked bar plots, ensuring the same colors are used across both plots.
    - Saves two PDFs: 'stat_intact_plot.pdf' and 'stat_frag_plot.pdf'.
"""

import argparse

import matplotlib
import pandas as pd

matplotlib.use("Agg")  # for non-interactive / cluster use
import matplotlib.pyplot as plt
import seaborn as sns

from ._figure import add_title_args, resolve_title


def read_te_data(file_path):
    """Reads the TE data from a TSV file and formats it for plotting."""
    df = pd.read_csv(file_path, sep='\t')

    # Convert to long format for easier plotting
    df_long = df.melt(id_vars=['TE_class/TE_superfamily'], 
                      var_name='Generation', 
                      value_name='Count')

    # Extract the integer generation number from column names
    df_long['Generation'] = df_long['Generation'].str.extract(r'gen(\d+)_final_Count').astype(int)
    
    return df_long

def plot_stacked_bar(df_long, output_file, color_map, title=None):
    """Generates and saves a stacked bar plot from the TE data."""
    # Pivot to get TE categories stacked per generation
    pivot_df = df_long.pivot(index='Generation', columns='TE_class/TE_superfamily', values='Count').fillna(0)
    
    # Sort generations numerically
    pivot_df = pivot_df.sort_index()

    # Set up the plot
    fig, ax = plt.subplots(figsize=(10, 7))
    pivot_df.plot(kind='bar', stacked=True, width=0.8, color=[color_map[cat] for cat in pivot_df.columns], ax=ax)

    # Formatting
    ax.set_xlabel("Generation", fontsize=14)
    ax.set_ylabel("TE Count", fontsize=14)
    if title:
        ax.set_title(title, fontsize=16, pad=20)
    
    # Improve legend readability
    legend = ax.legend(title="TE class/superfamily", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10, title_fontsize=12, frameon=False)
    plt.setp(legend.get_texts(), fontsize=9)

    # Save as PDF
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches="tight")
    print(f"Saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Stacked superfamily counts from stats_report.py output."
    )
    parser.add_argument("--intact", default="stat_intact.tsv",
                        help="Intact TE table (default: %(default)s)")
    parser.add_argument("--frag", default="stat_frag.tsv",
                        help="Fragmented TE table (default: %(default)s)")
    parser.add_argument("--out_prefix", default="stat",
                        help="Writes <prefix>_intact_plot.pdf and <prefix>_frag_plot.pdf "
                             "(default: %(default)s)")
    add_title_args(parser)
    args = parser.parse_args()

    intact_data = read_te_data(args.intact)
    frag_data = read_te_data(args.frag)

    # Get unique TE categories
    unique_classes = sorted(set(intact_data["TE_class/TE_superfamily"].unique()) | set(frag_data["TE_class/TE_superfamily"].unique()))

    # Assign consistent colors using seaborn palette
    color_palette = sns.color_palette("tab20", len(unique_classes))
    color_map = {cat: color_palette[i] for i, cat in enumerate(unique_classes)}

    # Generate plots
    _t = resolve_title(args, "Stacked TE Count Across Generations")
    plot_stacked_bar(intact_data, f"{args.out_prefix}_intact_plot.pdf", color_map, _t)
    plot_stacked_bar(frag_data, f"{args.out_prefix}_frag_plot.pdf", color_map, _t)

if __name__ == "__main__":
    main()
