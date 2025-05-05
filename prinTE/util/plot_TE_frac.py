#!/usr/bin/env python3
import argparse
import os
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_line(line):
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
    if ";" in name:
        parts = name.split(";")
        feature_id = parts[0]
        additional = parts[1:]
    else:
        feature_id = name
        additional = []
    return feature_id, additional

def process_bed_file(bed_file):
    feature_lengths = {
        "Intact gene": 0,
        "Fragmented gene": 0,
        "Intact TE": 0,
        "SoloLTR": 0,
        "Fragmented TE": 0
    }

    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            record = parse_line(line)
            if not record:
                continue
            feature_id, additional = parse_attributes(record['name'])

            if feature_id.startswith("gene"):
                category = "Intact gene" if not additional else "Fragmented gene"
            elif "_SOLO" in feature_id:
                category = "SoloLTR"
            elif "_FRAG" in feature_id:
                category = "Fragmented TE"
            else:
                category = "Fragmented TE" if additional and "CUT_BY" in additional[0] else "Intact TE"

            try:
                length = int(record['end']) - int(record['start'])
                feature_lengths[category] += length
            except ValueError:
                continue

    return feature_lengths

def get_genome_length(fasta_file):
    total_length = 0
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                continue
            total_length += len(line.strip())
    return total_length

def extract_gen_number(label):
    if "burnin" in label.lower():
        return 0
    match = re.search(r'gen(\d+)_', label)
    return int(match.group(1)) if match else None

def main():
    parser = argparse.ArgumentParser(
        description=("Classify BED entries from multiple files into gene and TE categories, "
                     "compute feature percentage over genome generations using a corresponding FASTA, "
                     "and produce summary stats and a publication-quality PDF plot."),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--bed", nargs="+", required=True, help="Input BED file(s).")
    parser.add_argument("--fasta", nargs="+", required=True, help="Input FASTA file(s).")
    parser.add_argument("--feature", required=True,
                        help="Colon-separated feature categories (e.g., 'Intact_TE:SoloLTR').")
    parser.add_argument("--out_prefix", required=True, help="Output prefix for stats files and plot (PDF).")
    parser.add_argument("--ymax", type=int, default=100, help="Maximum y-axis limit for the plot. Default is 100.")
    parser.add_argument("--vert_grid", type=int, default=None,
                        help="Interval at which to draw vertical grid lines on the plot (e.g., every 10 generations).")

    args = parser.parse_args()

    if len(args.bed) != len(args.fasta):
        raise ValueError("The number of BED files must equal the number of FASTA files.")

    features_to_plot = args.feature.split(":")
    feature_keys = [feat.replace("_", " ") for feat in features_to_plot]

    generation_data = []

    for bed_file, fasta_file in zip(args.bed, args.fasta):
        file_base = os.path.splitext(os.path.basename(bed_file))[0]
        feature_lengths = process_bed_file(bed_file)
        genome_length = get_genome_length(fasta_file)

        per_feature_percent = {}
        total_bp = 0
        for key in feature_keys:
            bp = feature_lengths.get(key, 0)
            per_feature_percent[key] = (bp / genome_length * 100) if genome_length > 0 else 0
            total_bp += bp

        overall_percent = (total_bp / genome_length * 100) if genome_length > 0 else 0
        # Attempt to extract the generation number from the file name.
        gen_num = extract_gen_number(file_base)
        if gen_num is None:
            # If no generation number can be extracted, assume generation zero.
            gen_num = 0

        generation_data.append((gen_num, file_base, overall_percent, per_feature_percent))

    generation_data.sort(key=lambda x: x[0])
    x_vals = [x[0] for x in generation_data]
    x_tick_labels = [str(x[0]) for x in generation_data]

    stacked_data = {key: [] for key in feature_keys}
    overall_percents = []
    for gen in generation_data:
        overall_percents.append(gen[2])
        per_feature = gen[3]
        for key in feature_keys:
            stacked_data[key].append(per_feature.get(key, 0))

    plt.figure(figsize=(8, 6))

    # Add vertical grid lines *before* plotting to keep them in the background
    if args.vert_grid:
        for idx, x in enumerate(x_vals):
            if idx % args.vert_grid == 0 and idx != 0:
                plt.axvline(x=x, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # Plot the data on top of grid
    if len(feature_keys) > 1:
        stack_arrays = [stacked_data[key] for key in feature_keys]
        for stack in stack_arrays:
            assert len(stack) == len(x_vals)
        plt.stackplot(x_vals, *stack_arrays, labels=feature_keys, zorder=2)
        plt.plot(x_vals, overall_percents, color='black', linewidth=1.5, zorder=3)
        plt.legend(loc='best')
    else:
        plt.plot(x_vals, overall_percents, marker='o', linestyle='-', zorder=2)

    plt.xlabel("Generation")
    plt.ylabel("Feature Percentage")
    plt.title(f"Feature Percentage Over Generations: {', '.join(feature_keys)}")
    plt.xticks(x_vals, x_tick_labels, rotation=45)
    plt.xlim(min(x_vals), max(x_vals))
    plt.ylim(0, args.ymax)

    plt.tight_layout()
    plt.savefig(args.out_prefix + ".pdf")
    plt.close()

    # Generate pie charts per generation
    pie_output = args.out_prefix + "_pies.pdf"
    plot_pie_charts(generation_data, feature_keys, pie_output)

def plot_pie_charts(generation_data, feature_keys, out_pdf):
    non_te_label = "Non-TE"
    
    # ─── build a consistent color map ───
    cmap = plt.get_cmap('tab10')
    color_map = {label: cmap(i) for i, label in enumerate(feature_keys)}
    color_map[non_te_label] = "gray"  # fixed color for Non-TE

    with PdfPages(out_pdf) as pdf:
        for gen_num, file_base, overall_percent, per_feature in generation_data:
            fig, ax = plt.subplots(figsize=(6, 6))

            labels = []
            sizes = []
            for key in feature_keys:
                pct = per_feature.get(key, 0)
                if pct > 0:
                    labels.append(key)
                    sizes.append(pct)

            total_feat = sum(sizes)
            remaining = 100.0 - total_feat
            if remaining > 0:
                labels.append(non_te_label)
                sizes.append(remaining)

            if not sizes:
                continue

            # consistent slice coloring
            slice_colors = [color_map[label] for label in labels]

            wedges, _ = ax.pie(
                sizes,
                startangle=90,
                counterclock=False,
                colors=slice_colors
            )

            ax.set_title(f"{file_base} (Generation {gen_num})")
            ax.legend(wedges, labels, title="Features", loc="center left", bbox_to_anchor=(1, 0.5))

            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)


if __name__ == "__main__":
    main()
