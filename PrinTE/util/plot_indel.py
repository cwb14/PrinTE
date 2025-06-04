#!/usr/bin/env python3
import os
import sys
import re
import subprocess
import argparse
import matplotlib.pyplot as plt

#####################
# Utility Functions #
#####################

def run_command(cmd):
    """Run a shell command and check for errors."""
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def index_genome(genome):
    """Run samtools faidx on a genome if the .fai file doesn't exist."""
    fai_file = genome + ".fai"
    if os.path.exists(fai_file):
        print(f"Index file {fai_file} exists; using it.")
    else:
        run_command(["samtools", "faidx", genome])
    return fai_file

def run_minimap2(ref, query):
    """Run minimap2 alignment if output PAF doesn't exist."""
    prefix = f"{os.path.splitext(ref)[0]}_{os.path.splitext(query)[0]}"
    paf_file = prefix + ".paf"
    if os.path.exists(paf_file):
        print(f"PAF file {paf_file} exists; using it.")
    else:
        cmd = ["minimap2", "-cx", "asm5", ref, query]
        print(f"Aligning {ref} -> {query}")
        with open(paf_file, "w") as outf:
            subprocess.run(cmd, stdout=outf, check=True)
    return paf_file

def parse_fai(fai_file):
    """
    Parse a .fai file and return a list of (chromosome, length) tuples.
    The order of lines is preserved.
    """
    chroms = []
    with open(fai_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                chroms.append((parts[0], int(parts[1])))
    return chroms

def parse_cigar(cigar_str):
    """
    Parse a CIGAR string (from the cg:Z: tag).
    Returns a list of (length, op) tuples.
    """
    return [(int(num), op) for num, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)]

def compute_indels(cigar_tuples):
    """
    Walk through the CIGAR tuples to compute indel events.
    Returns a list of dicts with:
      - "type": "I" or "D"
      - "pos": the position on the reference (relative to the alignment start)
      - "size": the indel size
    """
    ref_pos = 0
    events = []
    for length, op in cigar_tuples:
        if op in ("M", "=", "X"):
            ref_pos += length
        elif op == "I":
            events.append({"type": "I", "pos": ref_pos, "size": length})
        elif op == "D":
            events.append({"type": "D", "pos": ref_pos, "size": length})
            ref_pos += length
    return events

def parse_paf(paf_file):
    """
    Parse a PAF file and return a list of alignment events.
    Each event dict will include:
      - "ref": reference chromosome name
      - "pos": position on reference (absolute, i.e. alignment start offset + event pos)
      - "size": indel size
      - "type": "I" or "D"
    """
    events = []
    with open(paf_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 13:
                continue
            ref_name = fields[0]
            ref_start = int(fields[2])
            # Look for the cg:Z: tag containing the CIGAR string
            cigar = None
            for tag in fields[12:]:
                if tag.startswith("cg:Z:"):
                    cigar = tag.split(":", 2)[2]
                    break
            if cigar:
                tuples = parse_cigar(cigar)
                aln_events = compute_indels(tuples)
                # Offset the event position by the alignment's start on the reference.
                for ev in aln_events:
                    ev["pos"] += ref_start
                    ev["ref"] = ref_name
                events.extend(aln_events)
    return events

###########################
# Plotting Implementation #
###########################

def compute_chrom_offsets(chrom_list, gap=1e6):
    """
    Given a list of (chrom, length) tuples, compute the cumulative offset for each chromosome.
    Returns a dict mapping chromosome name to (offset, length).
    """
    offsets = {}
    current = 0
    for chrom, length in chrom_list:
        offsets[chrom] = (current, length)
        current += length + gap
    total_length = current - gap  # remove last gap
    return offsets, total_length

def plot_genomes(genome_fai_dict, genome_order, align_events_by_ref, output_pdf):
    """
    Plot each genome on its own horizontal row.
    For each genome, draw its chromosomes using their cumulative offsets.
    Then, for each alignment event assigned to that genome (by reference),
    plot an arrow (up for insertions, down for deletions) at the proper x position.
    """
    n = len(genome_order)
    row_height = 20
    gap_y = 15

    # Create a mapping: for each genome, compute chromosome offsets.
    genome_offsets = {}
    total_lengths = {}
    for genome in genome_order:
        chrom_list = genome_fai_dict[genome]
        offsets, tot = compute_chrom_offsets(chrom_list)
        genome_offsets[genome] = offsets
        total_lengths[genome] = tot

    # Determine overall x range as the maximum total length.
    overall_x = max(total_lengths.values())

    fig, ax = plt.subplots(figsize=(12, 2*n))
    y_ticks = []
    y_labels = []

    # Plot each genome as a horizontal bar.
    for idx, genome in enumerate(genome_order):
        y = idx * (row_height + gap_y)
        y_ticks.append(y + row_height/2)
        y_labels.append(os.path.basename(genome))
        offsets = genome_offsets[genome]
        # Draw each chromosome segment as a broken_barh.
        for chrom, (start, length) in offsets.items():
            ax.broken_barh([(start, length)], (y, row_height), facecolors="gray", edgecolors="black")
            # Label the chromosome (centered in its segment)
            ax.text(start + length/2, y + row_height/2, chrom, ha="center", va="center", fontsize=8, color="white")

        # Now, if this genome has alignment events (as the reference in a pair),
        # plot them.
        if genome in align_events_by_ref:
            for ev in align_events_by_ref[genome]:
                offsets_for_genome = genome_offsets[genome]
                # Only plot if the event's chromosome is found.
                if ev["ref"] in offsets_for_genome:
                    start_offset, _ = offsets_for_genome[ev["ref"]]
                    x_event = start_offset + ev["pos"]
                    if ev["type"] == "I":
                        # Upward arrow for insertion
                        ax.annotate(f"+{ev['size']}", xy=(x_event, y + row_height/2),
                                    xytext=(x_event, y + row_height/2 + 10),
                                    arrowprops=dict(facecolor="blue", width=2),
                                    color="blue", fontsize=6, ha="center")
                    elif ev["type"] == "D":
                        # Downward arrow for deletion
                        ax.annotate(f"-{ev['size']}", xy=(x_event, y + row_height/2),
                                    xytext=(x_event, y + row_height/2 - 10),
                                    arrowprops=dict(facecolor="red", width=2),
                                    color="red", fontsize=6, ha="center")

    ax.set_xlim(0, overall_x * 1.05)
    ax.set_ylim(-gap_y, n*(row_height + gap_y))
    ax.set_xlabel("Genomic coordinate")
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_title("Genome Alignment Presence/Absence Plot")
    plt.tight_layout()
    plt.savefig(output_pdf, format="pdf")
    print(f"Plot saved to {output_pdf}")

##################
# Main Execution #
##################

def main():
    parser = argparse.ArgumentParser(description="Genome alignment and indel plotting pipeline")
    parser.add_argument("--list", required=True, help="File containing list of genome fasta files (one per line)")
    parser.add_argument("--output", default="genome_alignment.pdf", help="Output PDF file for the plot")
    args = parser.parse_args()

    # Read genomes list
    with open(args.list) as f:
        genomes = [line.strip() for line in f if line.strip()]
    if len(genomes) < 2:
        sys.exit("Need at least two genomes for alignment.")

    # Index genomes and parse FAI files (preserving chromosome order)
    genome_fai_dict = {}
    for genome in genomes:
        fai_file = index_genome(genome)
        genome_fai_dict[genome] = parse_fai(fai_file)

    # For each adjacent genome pair, run minimap2 (if needed) and parse events.
    # We assign the events to the reference genome in the pair.
    align_events_by_ref = {g: [] for g in genomes}
    for i in range(len(genomes) - 1):
        ref = genomes[i]
        query = genomes[i+1]
        paf_file = run_minimap2(ref, query)
        events = parse_paf(paf_file)
        # Append these events under the reference genome key.
        align_events_by_ref[ref].extend(events)

    # Now plot genomes with their corresponding alignment events.
    plot_genomes(genome_fai_dict, genomes, align_events_by_ref, args.output)

if __name__ == "__main__":
    main()
