#!/usr/bin/env python3
import glob
import argparse
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot density distributions from gen[int]_LTR.tsv files.")
    parser.add_argument("--model", required=True, choices=["raw", "K2P", "JC69"],
                        help="Model to use (choose from raw, K2P, JC69).")
    parser.add_argument("--miu", required=True, type=float,
                        help="Mutation rate (used in the plot title and for time conversion).")
    parser.add_argument("--output", default="density_plot.pdf",
                        help="Output PDF file name.")
    return parser.parse_args()

def find_files():
    # Find all files that match the pattern gen[int]_LTR.tsv.
    pattern = "gen[0-9]*_LTR.tsv"
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError("No files found matching pattern 'gen[int]_LTR.tsv'")
    return files

def extract_source_label(filename):
    # Extract the numeric part from the filename and format as "Generation <number>".
    base = os.path.basename(filename)
    m = re.search(r'gen(\d+)', base)
    return f"Generation {m.group(1)}" if m else base

def load_data(files, model):
    # For each file, load the data and extract the required columns.
    data_list = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t')
            d_col = f"{model}_d"
            t_col = f"{model}_T"
            if d_col not in df.columns or t_col not in df.columns:
                print(f"Skipping {f}: missing required columns {d_col} or {t_col}.")
                continue
            # Extract a simplified source label.
            source_label = extract_source_label(f)
            df = df[[d_col, t_col]].rename(columns={d_col: "distance", t_col: "time"})
            df['source'] = source_label
            data_list.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")
    if not data_list:
        raise ValueError("No valid data loaded from files.")
    return pd.concat(data_list, ignore_index=True)

def create_density_plot(data, miu, model, output):
    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    
    # Calculate an appropriate x-axis maximum based on the data.
    max_distance = data['distance'].max()
    # Set x_max to 10% above the maximum value in the data.
    x_max = max_distance * 1.1 if max_distance > 0 else 1
    ax.set_xlim(0, x_max)
    
    # Plot a density curve for each source file.
    sources = data['source'].unique()
    palette = sns.color_palette("husl", len(sources))
    
    for i, source in enumerate(sources):
        subset = data[data['source'] == source]
        sns.kdeplot(subset['distance'], label=source, color=palette[i], fill=False, clip=(0, x_max), warn_singular=False)
    
    plt.xlabel(f"{model} genetic distance")
    plt.ylabel("Density")
    
    # Create a secondary x-axis converting genetic distance to time (MYA)
    # Conversion: time (years) = genetic_distance / miu, then convert years to MYA.
    def forward(x):
        return (x / miu) / 1e6
    def inverse(x):
        return x * miu * 1e6
    
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel("Time (MYA)")
    
    plt.title(f"Density Plot for {model} distances\nMutation Rate (miu) = {miu}")
    plt.legend(title="Generation")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print(f"Plot saved to {output}")

def main():
    args = parse_arguments()
    files = find_files()
    data = load_data(files, args.model)
    create_density_plot(data, args.miu, args.model, args.output)

if __name__ == "__main__":
    main()
