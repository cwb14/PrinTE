import argparse
import os
import re
import sys

import matplotlib
from Bio import SeqIO

matplotlib.use("Agg")  # for non-interactive / cluster use
import matplotlib.pyplot as plt
import pandas as pd

from ._figure import add_title_args, resolve_title


def calculate_genome_size(fasta_file):
    total_size = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))
    return total_size


def auto_scale_bp(values):
    """Return (divisor, unit_label) based on the magnitude of the values."""
    max_val = max(values)
    if max_val >= 1e9:
        return 1e9, "Gb"
    elif max_val >= 1e6:
        return 1e6, "Mb"
    elif max_val >= 1e3:
        return 1e3, "kb"
    return 1, "bp"


def auto_scale_generation(values):
    """Return (divisor, unit_label) based on the magnitude of generation counts."""
    max_val = max(values)
    if max_val >= 1e9:
        return 1e9, "Generation count (billions)"
    elif max_val >= 1e6:
        return 1e6, "Generation count (millions)"
    elif max_val >= 1e3:
        return 1e3, "Generation count (thousands)"
    elif max_val >= 1e2:
        return 1e2, "Generation count (hundreds)"
    return 1, "Generation count"


def main():
    parser = argparse.ArgumentParser(description="Plot genome size across generations.")
    parser.add_argument("reference", nargs="?", default=None,
                        help="Optional reference genome FASTA for a horizontal line.")
    parser.add_argument("--zero-y", action="store_true",
                        help="Start the Y-axis at zero.")
    parser.add_argument("--flip-x", action="store_true",
                        help="Flip the X-axis tick labels. If a reference genome is "
                             "also provided, plot it as a red point at generation 0 "
                             "instead of a dashed line.")
    parser.add_argument("--output", default="genome_size_plot.pdf",
                        help="Output PDF (default: %(default)s)")
    parser.add_argument("--no-legend", action="store_true",
                        help="Hide the legend.")
    add_title_args(parser)
    args = parser.parse_args()

    reference_fasta = args.reference
    if reference_fasta and not os.path.isfile(reference_fasta):
        print(f"Error: reference file '{reference_fasta}' not found.")
        sys.exit(1)

    # Find all files matching 'gen[int]_final.fasta'
    fasta_files = [f for f in os.listdir() if re.match(r"gen\d+_final\.fasta$", f)]

    if not fasta_files:
        print("No files matching 'gen[int]_final.fasta' found in the current directory.")
        sys.exit(1)

    # Extract genome size for each file
    genome_data = []
    for fasta_file in fasta_files:
        match = re.search(r"gen(\d+)_final\.fasta", fasta_file)
        if match:
            genome_id = int(match.group(1))
            genome_size = calculate_genome_size(fasta_file)
            genome_data.append((genome_id, genome_size))

    genome_data.sort()
    df = pd.DataFrame(genome_data, columns=["Generation", "Genome Size"])

    # Auto-scale axes
    bp_divisor, bp_unit = auto_scale_bp(df["Genome Size"].tolist())
    gen_divisor, gen_label = auto_scale_generation(df["Generation"].tolist())

    df["Genome Size Scaled"] = df["Genome Size"] / bp_divisor
    df["Generation Scaled"] = df["Generation"] / gen_divisor

    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(
        df["Generation Scaled"],
        df["Genome Size Scaled"],
        marker="o",
        linestyle="-",
        color="black",
        markersize=8,
    )

    # Reference genome handling
    ref_size_scaled = None
    if reference_fasta:
        ref_size = calculate_genome_size(reference_fasta)
        ref_size_scaled = ref_size / bp_divisor

        if args.flip_x:
            # Plot reference as a red point at generation 0 (rightmost after flip)
            tick_positions = list(df["Generation Scaled"])
            if len(tick_positions) >= 2:
                step = tick_positions[1] - tick_positions[0]
            else:
                step = tick_positions[0]
            ref_x = tick_positions[-1] + step
            # Dashed red line connecting last data point to reference point
            last_x = tick_positions[-1]
            last_y = df["Genome Size Scaled"].iloc[-1]
            plt.plot([last_x, ref_x], [last_y, ref_size_scaled],
                     linestyle="--", color="red", linewidth=1.5, zorder=4)
            plt.plot(ref_x, ref_size_scaled, marker="o", color="red", markersize=8,
                     zorder=5, label=f"Reference ({ref_size_scaled:.2f} {bp_unit})")
            if not args.no_legend:
                plt.legend(fontsize=12)
        else:
            # Plot reference as a horizontal dashed line
            plt.axhline(y=ref_size_scaled, color="red", linestyle="--", linewidth=1.5,
                         label=f"Reference ({ref_size_scaled:.2f} {bp_unit})")
            if not args.no_legend:
                plt.legend(fontsize=12)

    # ── Y-axis limits (single authoritative block) ──────────────────────
    # 1. Collect the extreme plotted values.
    y_data_max = df["Genome Size Scaled"].max()
    y_data_min = df["Genome Size Scaled"].min()
    if ref_size_scaled is not None:
        y_data_max = max(y_data_max, ref_size_scaled)
        y_data_min = min(y_data_min, ref_size_scaled)

    # 2. Determine the bottom boundary first (--zero-y may override).
    y_bottom = 0 if args.zero_y else y_data_min

    # 3. Compute padding as 5% of the *visible* range (bottom → data max),
    #    so headroom stays visually consistent regardless of --zero-y.
    visible_range = y_data_max - y_bottom
    padding = visible_range * 0.05 if visible_range > 0 else y_data_max * 0.05

    # 4. Apply: pull bottom down (unless pinned at 0), push top up.
    if not args.zero_y:
        y_bottom -= padding
    y_top = y_data_max + padding

    plt.ylim(bottom=y_bottom, top=y_top)

    if args.flip_x:
        gen_label = gen_label.replace("Generation count", "Generations in past count")
    plt.xlabel(gen_label, fontsize=14)
    plt.ylabel(f"Genome size ({bp_unit})", fontsize=14)
    _t = resolve_title(args, "Genome Size Distribution")
    if _t:
        plt.title(_t, fontsize=16)

    # X-axis ticks: flip labels if requested
    tick_positions = list(df["Generation Scaled"])
    if args.flip_x:
        flipped_labels = list(reversed(tick_positions))
        if reference_fasta:
            tick_positions.append(ref_x)
            flipped_labels.append(0)
        plt.xticks(tick_positions, [f"{v:g}" for v in flipped_labels])
    else:
        plt.xticks(tick_positions)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.savefig(args.output, bbox_inches="tight")
    print(f"Genome size plot saved as '{args.output}'.")



if __name__ == "__main__":
    main()
