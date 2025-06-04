import os
import re
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd

# Function to calculate genome size
def calculate_genome_size(fasta_file):
    total_size = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))
    return total_size

# Find all files matching 'gen[int]_final.fasta'
fasta_files = [f for f in os.listdir() if re.match(r'gen\d+_final\.fasta$', f)]

# Extract genome size for each file
genome_data = []
for fasta_file in fasta_files:
    match = re.search(r'gen(\d+)_final\.fasta', fasta_file)
    if match:
        genome_id = int(match.group(1))  # Extract numeric part
        genome_size = calculate_genome_size(fasta_file)
        genome_data.append((genome_id, genome_size))

# Sort by genome id
genome_data.sort()

# Convert to DataFrame
df = pd.DataFrame(genome_data, columns=["Genome ID", "Genome Size"])

# Plot
plt.figure(figsize=(8, 6))
plt.plot(df["Genome ID"], df["Genome Size"], marker="o", linestyle="-", color="black", markersize=8)
plt.xlabel("Genome ID", fontsize=14)
plt.ylabel("Genome Size (bp)", fontsize=14)
plt.title("Genome Size Distribution", fontsize=16)
plt.xticks(df["Genome ID"])
plt.grid(True, linestyle="--", alpha=0.6)
plt.savefig("genome_size_plot.pdf", bbox_inches="tight")

print("Genome size plot saved as 'genome_size_plot.pdf'.")
