#!/usr/bin/env python3
"""
plot_te_stats.py

Usage:
    python plot_te_stats.py --infile stats2.txt [--outfile output.pdf]

Description:
    Reads a TE stats file with sample data blocks like:

        SAMPLE: 1
        Total number of intact TEs from BED: 1000
        Distribution of intact TEs by (te_class, te_superfamily):
          LTR/Copia: 168
          MITE/Stow: 28
          ...

    It then plots a stacked bar graph for the breakdown of categories per sample,
    and overlays a line graph of the total TEs per sample.
    
    The output is saved as a publication quality PDF.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import re

def parse_file(infile):
    """
    Parses the input file into a list of sample dictionaries.
    
    Each dictionary has:
       - 'sample': sample number (int)
       - 'total': total count (int)
       - 'categories': a dict mapping category names to counts
    """
    samples = []
    current_sample = None

    with open(infile, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            # Look for a new sample block
            if line.startswith("SAMPLE:"):
                # Save the previous sample (if any)
                if current_sample is not None:
                    samples.append(current_sample)
                # Initialize new sample
                m = re.search(r"SAMPLE:\s*(\d+)", line)
                if m:
                    sample_id = int(m.group(1))
                    current_sample = {'sample': sample_id, 'total': None, 'categories': {}}
            # Look for the total count line
            elif line.startswith("Total number of intact TEs from BED:"):
                m = re.search(r"Total number of intact TEs from BED:\s*(\d+)", line)
                if m:
                    total = int(m.group(1))
                    current_sample['total'] = total
            # Skip the header for the distribution section
            elif "Distribution of intact TEs by" in line:
                continue
            # Blank lines are skipped
            elif line.strip() == "":
                continue
            else:
                # Expecting lines like "  LTR/Copia: 168"
                m = re.match(r"^\s*(.+):\s*(\d+)", line)
                if m:
                    category = m.group(1).strip()
                    count = int(m.group(2))
                    current_sample['categories'][category] = count

        # Append the last sample if it exists
        if current_sample is not None:
            samples.append(current_sample)
    return samples

def main():
    parser = argparse.ArgumentParser(
        description="Plot TE statistics as stacked bars with an overlaid total line plot."
    )
    parser.add_argument('--infile', required=True, help="Input file with TE statistics")
    parser.add_argument('--outfile', default="output.pdf", help="Output PDF file (default: output.pdf)")
    args = parser.parse_args()

    # Parse the data file
    samples = parse_file(args.infile)

    # Sort samples by sample number
    samples.sort(key=lambda s: s['sample'])
    sample_ids = [s['sample'] for s in samples]
    total_counts = [s['total'] for s in samples]

    # Get the union of all categories across all samples
    all_categories = set()
    for s in samples:
        all_categories.update(s['categories'].keys())
    # Sorting the categories (alphabetically; modify as needed)
    all_categories = sorted(all_categories)

    # Prepare data for the stacked bar plot:
    # For each category, build a list of counts per sample (fill missing with 0)
    category_data = {}
    for cat in all_categories:
        category_data[cat] = [s['categories'].get(cat, 0) for s in samples]

    # Begin plotting
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(12, 8))

    import itertools
    import matplotlib.colors as mcolors

    # Use multiple colormaps to generate a large, vibrant palette
    cmap1 = plt.get_cmap("tab20")    # Standard vibrant
    cmap2 = plt.get_cmap("Paired")   # More contrast, avoids dark colors
    cmap3 = plt.get_cmap("Set3")     # Light pastel colors

    # Get colors from all three colormaps
    colors_combined = list(cmap1.colors) + list(cmap2.colors) + list(cmap3.colors)

    # Manually filter out dark colors (optional)
    filtered_colors = [c for c in colors_combined if mcolors.rgb_to_hsv(c)[2] > 0.5]

    # Ensure we have enough colors and cycle through them if needed
    color_pickers = itertools.cycle(filtered_colors)
    colors = [next(color_pickers) for _ in range(len(all_categories))]

    # Plot stacked bars
    bottoms = np.zeros(len(samples))
    for i, cat in enumerate(all_categories):
        counts = category_data[cat]
        ax.bar(sample_ids, counts, bottom=bottoms, color=colors[i], label=cat)
        bottoms += np.array(counts)

    # Overlay the line graph for total counts
    ax.plot(sample_ids, total_counts, color='black', marker='o',
            linestyle='-', linewidth=2, label='Total TEs')

    # Labeling the axes and title
    ax.set_xlabel("Sample", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)
    ax.set_title("bash TEvo/wrapper3.sh -P 10 -cn 1 -sz 1Mb --TE_num 1000 -st 1000 -ge 100000 -sr 2 -F 100,100", fontsize=16)

    # Place the legend outside the plot if needed
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=10)

    plt.tight_layout()
    plt.savefig(args.outfile, format='pdf')
    print(f"Plot saved to {args.outfile}")

if __name__ == '__main__':
    main()
