import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def plot_grouped_line(stat_file, output_pdf):
    # Read the data
    df = pd.read_csv(stat_file, sep='\t')
    
    # Extract groups from column names
    groups = [col.split('_')[0][3:] for col in df.columns[1:]]
    
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
    ax.set_title("Grouped Line Plot")
    ax.legend(title="Category")
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Save to PDF
    plt.tight_layout()
    plt.savefig(output_pdf, format='pdf', dpi=300)
    print(f"Plot saved to {output_pdf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate grouped line plot from stat_overall.tsv.")
    args = parser.parse_args()
    
    stat_file = "stat_overall.tsv"
    output_pdf = "grouped_line_plot.pdf"
    
    if not os.path.exists(stat_file):
        print(f"Error: {stat_file} not found.")
    else:
        plot_grouped_line(stat_file, output_pdf)
