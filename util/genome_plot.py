import os
import re
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd


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


import argparse

parser = argparse.ArgumentParser(description="Plot genome size across generations.")
parser.add_argument("reference", nargs="?", default=None,
                    help="Optional reference genome FASTA for a horizontal line.")
parser.add_argument("--zero-y", action="store_true",
                    help="Start the Y-axis at zero.")
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

# Optional reference genome horizontal line
if reference_fasta:
    ref_size = calculate_genome_size(reference_fasta)
    ref_size_scaled = ref_size / bp_divisor
    plt.axhline(y=ref_size_scaled, color="red", linestyle="--", linewidth=1.5,
                label=f"Reference ({ref_size_scaled:.2f} {bp_unit})")
    plt.legend(fontsize=12)

# Ensure some headroom above the highest visible element
y_max = df["Genome Size Scaled"].max()
if reference_fasta:
    y_max = max(y_max, ref_size_scaled)
plt.ylim(top=y_max * 1.05)

plt.xlabel(gen_label, fontsize=14)
plt.ylabel(f"Genome size ({bp_unit})", fontsize=14)
plt.title("Genome Size Distribution", fontsize=16)
plt.xticks(df["Generation Scaled"])
plt.grid(True, linestyle="--", alpha=0.6)
if args.zero_y:
    plt.ylim(bottom=0)
plt.savefig("genome_size_plot.pdf", bbox_inches="tight")
print("Genome size plot saved as 'genome_size_plot.pdf'.")
