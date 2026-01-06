#!/usr/bin/env python3
import argparse
import math
import re
from typing import Tuple, List

import matplotlib
matplotlib.use("Agg")  # for non-interactive / cluster use
import matplotlib.pyplot as plt

import numpy as np
from scipy.stats import wasserstein_distance, ks_2samp


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare reference vs experimental genomic datasets using multiple metrics."
    )
    parser.add_argument(
        "--ref-tsv", required=True,
        help="Reference TSV file."
    )
    parser.add_argument(
        "--exp-tsv", required=True,
        help="Experimental TSV file."
    )
    parser.add_argument(
        "--ref-fasta", required=True,
        help="Reference genome FASTA (reference.fa)."
    )
    parser.add_argument(
        "--exp-fasta", required=True,
        help="Experimental genome FASTA (experimental.fa)."
    )
    parser.add_argument(
        "--alphas", "-a", nargs=4, type=float, required=True,
        metavar=("ALPHA_COUNT", "ALPHA_LENGTH", "ALPHA_DIST", "ALPHA_GENOME"),
        help="Four weights (α1..α4) for sample-count, cumulative-length, "
             "distribution distance, and genome-size distances."
    )
    parser.add_argument(
        "--pdist-col", type=int, default=7,
        help="1-based column index of p-distance (default: 7)."
    )
    parser.add_argument(
        "--plot-pdf", default=None,
        help="Output PDF filename for the p-distance distribution plot (optional)."
    )
    parser.add_argument(
        "--ignore-header", action="store_true",
        help="If set, skip the first line of each TSV as header."
    )
    parser.add_argument(
        "--dist-metric",
        choices=["wasserstein", "ks", "rms"],
        default="wasserstein",
        help=(
            "Distance metric for p-distance distributions: "
            "'wasserstein' (default), 'ks' (Kolmogorov–Smirnov), or "
            "'rms' (root-mean-square difference between binned densities)."
        ),
    )
    return parser.parse_args()


#interval_regex = re.compile(r"^(.+):(\d+)-(\d+)$")
interval_regex = re.compile(r"^(.+?):(\d+)-(\d+)")


def parse_interval(interval_str: str) -> int:
    """
    Parse 'chrom:start-end' and return interval length.
    Assumes inclusive coordinates -> length = end - start + 1.
    """
    m = interval_regex.match(interval_str.strip())
    if not m:
        raise ValueError(f"Could not parse interval from: {interval_str}")
    start = int(m.group(2))
    end = int(m.group(3))
    if end < start:
        raise ValueError(f"End < start in interval: {interval_str}")
    # inclusive genomic interval
    return end - start + 1


def read_tsv_metrics(tsv_path: str, pdist_col_1based: int, ignore_header: bool
                     ) -> Tuple[int, int, List[float]]:
    """
    Return:
      sample_count: number of (non-empty, non-comment) lines
      cumulative_length: sum of interval lengths from column1
      p_dists: list of p-distance values from given column (1-based index)
    """
    sample_count = 0
    cumulative_length = 0
    p_dists: List[float] = []

    col_idx = pdist_col_1based - 1
    with open(tsv_path) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ignore_header and i == 0:
                # treat first non-comment line as header
                continue

            parts = line.split("\t")
            if len(parts) <= col_idx:
                # Not enough columns; skip
                continue

            # Sample count
            sample_count += 1

            # Cumulative length (from column 1)
            try:
                length = parse_interval(parts[0])
                cumulative_length += length
            except Exception as e:
                raise ValueError(
                    f"Error parsing interval in file {tsv_path} line {i+1}: {e}"
                )

            # p-distance
            pdist_str = parts[col_idx]
            if pdist_str not in ("NA", "", "NaN", "nan", "NAN"):
                try:
                    p_dists.append(float(pdist_str))
                except ValueError:
                    # If header or weird text survived, skip
                    pass

    return sample_count, cumulative_length, p_dists


def fasta_genome_size(fasta_path: str) -> int:
    """
    Sum the lengths of all sequences in a FASTA.
    Ignore header lines (starting with '>') and whitespace.
    """
    total = 0
    with open(fasta_path) as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total


def relative_difference(ref_value: float, exp_value: float) -> float:
    """
    Return |ref - exp| / ref, or 0 if ref == 0.
    This is dimensionless and can exceed 1 if exp is more than 2x ref.
    """
    if ref_value == 0:
        return 0.0
    return abs(ref_value - exp_value) / float(ref_value)


def compute_distribution_distance(ref_pdists: List[float],
                                  exp_pdists: List[float],
                                  metric: str) -> float:
    """
    Compute a scalar distance between two p-distance distributions.

    All returned distances are **dimensionless** and intended to be
    comparable in scale to |ref - exp| / ref, i.e.:
      - 0   → no difference
      - ~0.1→ ~10% difference
      - values > 1 allowed for very large differences

    Normalization details:
      - 'wasserstein':
          Compute 1D Wasserstein distance W (same units as p-distance),
          then normalize by the combined p-distance range R:
              R = max(all pdists) - min(all pdists)
              d = W / R
          If R == 0 (all values identical), return 0.

      - 'ks':
          Kolmogorov–Smirnov statistic is already dimensionless in [0,1],
          measuring the maximum CDF difference. We return it unchanged.

      - 'rms':
          Bin both distributions over the same range with 'bins' bins,
          using numpy.histogram(..., density=True). Let RMS_raw be the
          root-mean-square of the bin-wise density differences.

          For B bins spanning range R with uniform width and density=True,
          the **maximum possible** RMS (RMS_max) occurs when each
          distribution puts all its mass in a different single bin:
              RMS_max = sqrt(2 * B) / R

          We return:
              d = RMS_raw / RMS_max

          This yields a dimensionless quantity that is on a comparable
          scale to the other distances. Due to numerical and sampling
          effects, the returned value may slightly exceed 1 in extreme
          cases, which is acceptable.
    """
    if not ref_pdists or not exp_pdists:
        return float("nan")

    ref_arr = np.array(ref_pdists)
    exp_arr = np.array(exp_pdists)

    # Combined range of p-distances, used for normalization of
    # Wasserstein and RMS-based metrics.
    data_min = float(min(ref_arr.min(), exp_arr.min()))
    data_max = float(max(ref_arr.max(), exp_arr.max()))
    data_range = data_max - data_min

    if metric == "wasserstein":
        if data_range == 0.0:
            # All values identical -> no effective distributional difference
            return 0.0
        w_raw = float(wasserstein_distance(ref_arr, exp_arr))
        # Normalize by p-distance range to obtain a dimensionless quantity.
        return w_raw / data_range

    elif metric == "ks":
        stat, _ = ks_2samp(ref_arr, exp_arr)
        # KS statistic is already dimensionless in [0,1].
        return float(stat)

    elif metric == "rms":
        if data_range == 0.0:
            # All values identical -> no distributional difference
            return 0.0

        # Bin both distributions over the same range and compute the
        # RMS of the bin-wise density differences, then normalize by
        # the theoretical maximum RMS for the chosen binning and range.
        bins = 50
        hist_ref, bin_edges = np.histogram(
            ref_arr, bins=bins, range=(data_min, data_max), density=True
        )
        hist_exp, _ = np.histogram(
            exp_arr, bins=bin_edges, density=True
        )
        diff = hist_ref - hist_exp
        rms_raw = float(np.sqrt(np.mean(diff ** 2)))

        # Theoretical maximum RMS for B bins across range R, assuming
        # each distribution has all its mass in a different single bin
        # and using density=True:
        #   RMS_max = sqrt(2 * B) / R
        B = float(bins)
        rms_max = math.sqrt(2.0 * B) / data_range

        if rms_max == 0.0:
            return 0.0

        return rms_raw / rms_max

    else:
        raise ValueError(f"Unknown distribution metric: {metric}")


def plot_pdistance_distributions(ref_pdists: List[float],
                                 exp_pdists: List[float],
                                 dist_value: float,
                                 metric: str,
                                 pdf_path: str):
    """
    Save a PDF showing the density/histogram of p-distance distributions
    for reference and experimental datasets, annotated with the chosen
    (normalized, dimensionless) distribution distance.
    """
    if not ref_pdists or not exp_pdists:
        print("Warning: one of the p-distance lists is empty; skipping plot.")
        return

    ref_arr = np.array(ref_pdists)
    exp_arr = np.array(exp_pdists)

    plt.figure(figsize=(6, 4))

    # Histograms as density estimates
    bins = 50
    plt.hist(ref_arr, bins=bins, density=True, alpha=0.4, label="Reference")
    plt.hist(exp_arr, bins=bins, density=True, alpha=0.4, label="Experimental")

    metric_label = {
        "wasserstein": "Wasserstein (normalized)",
        "ks": "Kolmogorov–Smirnov",
        "rms": "RMS (normalized binned density)"
    }.get(metric, metric)

    plt.xlabel("p-distance")
    plt.ylabel("Density")
    plt.title(f"p-distance distributions\n{metric_label} distance = {dist_value:.6f}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(pdf_path)
    plt.close()


def main():
    args = parse_args()

    # --- TSV metrics ---
    ref_count, ref_cumlen, ref_pdists = read_tsv_metrics(
        args.ref_tsv, args.pdist_col, args.ignore_header
    )
    exp_count, exp_cumlen, exp_pdists = read_tsv_metrics(
        args.exp_tsv, args.pdist_col, args.ignore_header
    )

    # --- Genome sizes ---
    ref_gsize = fasta_genome_size(args.ref_fasta)
    exp_gsize = fasta_genome_size(args.exp_fasta)

    # --- Distance metrics ---
    d_count = relative_difference(ref_count, exp_count)
    d_cumlen = relative_difference(ref_cumlen, exp_cumlen)

    d_dist = compute_distribution_distance(ref_pdists, exp_pdists, args.dist_metric)

    d_gsize = relative_difference(ref_gsize, exp_gsize)

    alpha_count, alpha_len, alpha_dist, alpha_gsize = args.alphas

    # Composite metric: sum(alpha_i * distance_i)
    component_values = np.array([d_count, d_cumlen, d_dist, d_gsize], dtype=float)
    alpha_values = np.array([alpha_count, alpha_len, alpha_dist, alpha_gsize],
                            dtype=float)

    composite = float(np.nansum(alpha_values * component_values))

    metric_label = {
        "wasserstein": "Wasserstein (normalized)",
        "ks": "Kolmogorov–Smirnov",
        "rms": "RMS (normalized binned density)"
    }.get(args.dist_metric, args.dist_metric)

    # --- Print results ---
    print("=== Reference vs Experimental comparison ===")
    print(f"Reference TSV:   {args.ref_tsv}")
    print(f"Experimental TSV:{args.exp_tsv}")
    print(f"Reference FASTA: {args.ref_fasta}")
    print(f"Experimental FASTA:{args.exp_fasta}")
    print()

    print("--- Raw values ---")
    print(f"Sample count:    ref={ref_count}, exp={exp_count}")
    print(f"Cumulative length (bp):   ref={ref_cumlen}, exp={exp_cumlen}")
    print(f"Genome size (bp):         ref={ref_gsize},  exp={exp_gsize}")
    if not math.isnan(d_dist):
        # Note: d_dist is now a normalized, dimensionless distance.
        print(f"{metric_label} distance (p-distance distributions): {d_dist:.6g}")
    else:
        print(f"{metric_label} distance: NaN (insufficient p-distance data)")
    print()

    print("--- Per-metric distances (dimensionless) ---")
    print(f"d_sample_count (|Δ| / ref):      {d_count:.6g}")
    print(f"d_cumulative_length (|Δ| / ref): {d_cumlen:.6g}")
    print(f"d_genome_size (|Δ| / ref):       {d_gsize:.6g}")
    print(f"d_distribution ({metric_label}): {d_dist:.6g}")
    print()

    print("--- Weights (α) ---")
    print(f"α_sample_count = {alpha_count}")
    print(f"α_cumulative_length = {alpha_len}")
    print(f"α_distribution = {alpha_dist}")
    print(f"α_genome_size = {alpha_gsize}")
    print()

    print("--- Composite metric ---")
    print(f"Composite = Σ α_i * d_i = {composite:.6g}")
    print()

    # --- Plot, if requested ---
    if args.plot_pdf:
        if math.isnan(d_dist):
            print("Skipping PDF plot because distribution distance is NaN.")
        else:
            plot_pdistance_distributions(ref_pdists, exp_pdists, d_dist,
                                         args.dist_metric, args.plot_pdf)
            print(f"Distribution plot written to: {args.plot_pdf}")


if __name__ == "__main__":
    main()
