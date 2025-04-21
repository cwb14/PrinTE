#!/usr/bin/env python3

# With prinTE bed input, it plots seq length density for intact TE, intact gene, SoloLTR, Fragmented TE, and fragmented gene. 
# With prinTE fasta input, it plots seq length density for each unique TE superfamily. 

import sys
import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

# ------------------------------------------------------------------------------
# Helper functions for BED classification (from your provided script)
# ------------------------------------------------------------------------------
def parse_line(line):
    parts = line.strip().split("\t")
    if len(parts) < 6:
        return None
    chrom, start, end, name, tsd, strand = parts[:6]
    return {
        'chrom': chrom,
        'start': int(start),
        'end': int(end),
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

def extract_TE_info(feature_id):
    try:
        te_name, rest = feature_id.split("#", 1)
        te_class, te_superfamily = rest.split("/", 1)
        te_superfamily = te_superfamily.split("~")[0]
        return te_name, te_class, te_superfamily
    except Exception:
        return None, None, None

def classify_bed_records(lines):
    # First pass: gene vs TE candidates
    for rec in lines:
        fid = rec['feature_id']
        add = rec['additional']
        if fid.startswith("gene"):
            rec['category'] = "Intact gene" if not add else "Fragmented gene"
        else:
            if "_SOLO" in fid:
                rec['category'] = "SoloLTR"
            elif add and "CUT_BY" in add[0]:
                rec['category'] = "Fragmented TE"
            else:
                rec['category'] = "Potential intact TE"

    # Rule (3): demote some Potential intact TE to Fragmented TE
    for i, rec in enumerate(lines):
        if rec['category'] != "Potential intact TE" or not rec['additional']:
            continue
        for j in range(max(0, i-100), min(len(lines), i+101)):
            if i == j:
                continue
            other = lines[j]
            if other['feature_id'].startswith("gene"):
                continue
            if (other['tsd'] == rec['tsd'] and
                other['strand'] == rec['strand'] and
                (rec['name'].startswith(other['name']) or other['name'].startswith(rec['name']))):
                rec['category'] = "Fragmented TE"
                break

    # Final pass: promote remaining Potential intact TE
    for rec in lines:
        if rec['category'] == "Potential intact TE":
            rec['category'] = "Intact TE"

    return lines

def parse_bed_with_category(bed_file):
    lines = []
    with open(bed_file) as f:
        for L in f:
            if L.startswith("#") or not L.strip():
                continue
            rec = parse_line(L)
            if not rec:
                continue
            fid, add = parse_attributes(rec['name'])
            rec['feature_id'] = fid
            rec['additional'] = add
            lines.append(rec)

    # classify each record in place
    classify_bed_records(lines)

    # build DataFrame of length + category
    df = pd.DataFrame(lines)
    df['length'] = df['end'] - df['start']
    return df[['length', 'category']]

# ------------------------------------------------------------------------------
# FASTA parser that extracts TE_superfamily as category
# ------------------------------------------------------------------------------
def parse_fasta_with_superfamily(fasta_file):
    recs = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        feature_id, _ = parse_attributes(record.id)
        _, te_class, te_superfamily = extract_TE_info(feature_id)
        cat = te_superfamily if te_superfamily is not None else "Unknown"
        recs.append({'length': len(record.seq), 'category': cat})
    df = pd.DataFrame(recs)
    return df

# ------------------------------------------------------------------------------
# Plotting function
# ------------------------------------------------------------------------------
def plot_density(df, output_pdf, title, hue_col):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.kdeplot(
        data=df,
        x="length",
        hue=hue_col,
        common_norm=False,   # so each category is normalized separately
        fill=True,
        linewidth=2
    )
    plt.title(title, fontsize=16)
    plt.xlabel("Length (bp)", fontsize=14)
    plt.ylabel("Density", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_pdf, format="pdf")
    print(f"Saved density plot to {output_pdf}")

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main(input_file, output_pdf):
    ext = os.path.splitext(input_file)[1].lower()
    if ext == ".bed":
        df = parse_bed_with_category(input_file)
        title = "Feature Length Density (BED)"
        hue_col = "category"
    elif ext in (".fa", ".fasta"):
        df = parse_fasta_with_superfamily(input_file)
        title = "Sequence Length Density (FASTA)"
        hue_col = "category"
    else:
        sys.stderr.write("Error: Unsupported file format. Use .bed, .fa, or .fasta\n")
        sys.exit(1)

    plot_density(df, output_pdf, title, hue_col)

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Plot length‐density shaded by category (BED) or TE_superfamily (FASTA)."
    )
    p.add_argument("input", help="Input .bed | .fa | .fasta file")
    p.add_argument("output", help="Output PDF filename")
    args = p.parse_args()
    main(args.input, args.output)
