import matplotlib
import pandas as pd

matplotlib.use("Agg")  # for non-interactive / cluster use
import argparse
import os
import sys

import matplotlib.pyplot as plt

from ._figure import add_title_args, resolve_title


def plot_grouped_line(stat_file, output_pdf, title=None):
    # Read the data
    df = pd.read_csv(stat_file, sep='\t')
    
    # Melt dataframe for plotting
    df_melted = df.melt(id_vars=['Category'], var_name='Group', value_name='Count')
    df_melted['Group'] = df_melted['Group'].str.extract(r'gen(\d+)_')[0]
    
    # Convert group to categorical with sorted order
    df_melted['Group'] = pd.Categorical(df_melted['Group'], categories=sorted(df_melted['Group'].unique()), ordered=True)
    
    # Set plot aesthetics
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot line for each category
    for category in df['Category']:
        subset = df_melted[df_melted['Category'] == category]
        ax.plot(subset['Group'], subset['Count'], marker='o', linestyle='-', label=category)
    
    ax.set_xlabel("Groups")
    ax.set_ylabel("Counts")
    if title:
        ax.set_title(title)
    ax.legend(title="Category")
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Save to PDF
    plt.tight_layout()
    plt.savefig(output_pdf, format='pdf', dpi=300)
    print(f"Plot saved to {output_pdf}")

def main():
    parser = argparse.ArgumentParser(
        description="Plot TE category counts across generations from stats_report.py output."
    )
    parser.add_argument("--input", default="stat_overall.tsv",
                        help="Category table written by stats_report.py (default: %(default)s)")
    parser.add_argument("--output", default="grouped_line_plot.pdf",
                        help="Output PDF (default: %(default)s)")
    add_title_args(parser)
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Error: {args.input} not found.")
    plot_grouped_line(args.input, args.output, resolve_title(args, "Grouped Line Plot"))


if __name__ == "__main__":
    main()
