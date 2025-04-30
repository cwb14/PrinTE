#!/usr/bin/env python3
import argparse
import os
import re
import matplotlib.pyplot as plt

def parse_line(line):
    """Parse a BED line into its columns."""
    parts = line.strip().split("\t")
    if len(parts) < 6:
        return None
    chrom, start, end, name, tsd, strand = parts[:6]
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'name': name,
        'tsd': tsd,
        'strand': strand
    }

def parse_attributes(name):
    """Split the NAME into feature_ID and additional attributes."""
    if ";" in name:
        parts = name.split(";")
        feature_id = parts[0]
        additional = parts[1:]
    else:
        feature_id = name
        additional = []
    return feature_id, additional

def extract_TE_info(feature_id):
    """Extract TE information from feature_ID."""
    try:
        te_name, rest = feature_id.split("#", 1)
        te_class, te_superfamily = rest.split("/", 1)
        te_superfamily = te_superfamily.split("~")[0]  # Remove junk
        return te_name, te_class, te_superfamily
    except Exception:
        return None, None, None

def process_bed_file(bed_file):
    """Process a BED file to count SoloLTRs and IntactLTRs."""
    overall_counts = {"SoloLTR": 0, "Intact TE": 0}
    intact_TE_counts = {}

    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            record = parse_line(line)
            if not record:
                continue
            feature_id, additional = parse_attributes(record['name'])

            # Exclude fragments
            if "_FRAG" in feature_id:
                continue

            _, te_class, te_superfamily = extract_TE_info(feature_id)

            # Count SoloLTRs
            if "_SOLO" in feature_id:
                overall_counts["SoloLTR"] += 1

            # Count IntactLTRs (where TE_class is "LTR")
            elif te_class == "LTR":
                overall_counts["Intact TE"] += 1
                key = f"{te_class}/{te_superfamily}"
                intact_TE_counts[key] = intact_TE_counts.get(key, 0) + 1

    return overall_counts, intact_TE_counts

def extract_generation(filename):
    """Extract the generation number from a filename."""
    base = os.path.basename(filename)
    if base.lower() == "burnin.bed":
        return 0
    match = re.search(r"gen(\d+)_final", base)
    return int(match.group(1)) if match else None

def main():
    parser = argparse.ArgumentParser(
        description="Plot SoloLTR/IntactLTR ratio across generations from BED files."
    )
    parser.add_argument("--bed", nargs="+", required=True, help="Input BED file(s)")
    parser.add_argument("--out_prefix", required=True, help="Output prefix for PDF file")
    args = parser.parse_args()

    data = []

    for bed_file in args.bed:
        gen = extract_generation(bed_file)
        if gen is None:
            print(f"Warning: Could not extract generation from {bed_file}. Skipping.")
            continue

        overall_counts, intact_TE_counts = process_bed_file(bed_file)
        solo_count = overall_counts["SoloLTR"]
        intact_ltr_count = overall_counts["Intact TE"]

        # Compute ratio, avoid division by zero
        ratio = solo_count / intact_ltr_count if intact_ltr_count > 0 else 0
        data.append((gen, ratio))

    # Sort data by generation
    data.sort()
    generations, ratios = zip(*data) if data else ([], [])

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(generations, ratios, marker="o", linestyle="-", color="black", linewidth=2)
    ax.set_xlabel("Generation", fontsize=14)
    ax.set_ylabel("SoloLTR/IntactLTR Ratio", fontsize=14)
    ax.set_title("SoloLTR/IntactLTR Ratio Across Generations", fontsize=16)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.set_xticks(generations)  # Ensure integer ticks

    plt.tight_layout()
    out_pdf = args.out_prefix + ".pdf"
    plt.savefig(out_pdf, format="pdf")
    print(f"Plot saved to {out_pdf}")

if __name__ == "__main__":
    main()
