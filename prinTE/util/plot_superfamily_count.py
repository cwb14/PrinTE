#!/usr/bin/env python3
"""
plot_te_stacked_bar.py

Usage:
    python plot_te_stacked_bar.py

Description:
    - Reads 'stat_intact.tsv' and 'stat_frag.tsv'.
    - Extracts TE counts per generation.
    - Generates stacked bar plots, ensuring the same colors are used across both plots.
    - Saves two PDFs: 'stat_intact_plot.pdf' and 'stat_frag_plot.pdf'.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

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

def plot_stacked_bar(df_long, output_file, color_map):
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
    ax.set_title(f"Stacked TE Count Across Generations", fontsize=16, pad=20)
    
    # Improve legend readability
    legend = ax.legend(title="TE class/superfamily", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10, title_fontsize=12, frameon=False)
    plt.setp(legend.get_texts(), fontsize=9)

    # Save as PDF
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches="tight")
    print(f"Saved: {output_file}")

def main():
    # Read both files
    intact_data = read_te_data("stat_intact.tsv")
    frag_data = read_te_data("stat_frag.tsv")

    # Get unique TE categories
    unique_classes = sorted(set(intact_data["TE_class/TE_superfamily"].unique()) | set(frag_data["TE_class/TE_superfamily"].unique()))

    # Assign consistent colors using seaborn palette
    color_palette = sns.color_palette("tab20", len(unique_classes))
    color_map = {cat: color_palette[i] for i, cat in enumerate(unique_classes)}

    # Generate plots
    plot_stacked_bar(intact_data, "stat_intact_plot.pdf", color_map)
    plot_stacked_bar(frag_data, "stat_frag_plot.pdf", color_map)

if __name__ == "__main__":
    main()
