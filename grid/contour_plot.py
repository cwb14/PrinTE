#!/usr/bin/env python3
"""
Make contour plots of Composite metric vs. pairs of parameters
from a random grid search sample.

Requires:
  - pandas
  - matplotlib

Example usage:
  python plot_sensitivity_contour.py \
      --input params.tsv \
      --metric Composite \
      --output sensitivity_contours.pdf \
      --min_insertion_rate 0.2e-12 --max_insertion_rate 0.1e-11 \
      --min_solo_ratio 30 --max_solo_ratio 70 \
      --min_deletion_rate 0.06e-12 --max_deletion_rate 0.06e-11
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.backends.backend_pdf import PdfPages


def parse_args():
    parser = argparse.ArgumentParser(
        description="Contour plots of quality metric vs. parameter pairs"
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Input table (tab-separated) with parameters and metric"
    )
    parser.add_argument(
        "--metric", "-m", default="Composite",
        help="Name of quality metric column (default: Composite)"
    )
    parser.add_argument(
        "--output", "-o", default="contour_plot.pdf",
        help="Output multi-page PDF filename (default: contour_plot.pdf)"
    )
    parser.add_argument(
        "--levels", "-L", type=int, default=20,
        help="Number of contour levels (default: 20)"
    )

    # --- Parameter range filters ---
    # insertion_rate
    parser.add_argument(
        "--min_insertion_rate", type=float, default=None,
        help="Minimum insertion_rate to include (optional)"
    )
    parser.add_argument(
        "--max_insertion_rate", type=float, default=None,
        help="Maximum insertion_rate to include (optional)"
    )

    # deletion_rate
    parser.add_argument(
        "--min_deletion_rate", type=float, default=None,
        help="Minimum deletion_rate to include (optional)"
    )
    parser.add_argument(
        "--max_deletion_rate", type=float, default=None,
        help="Maximum deletion_rate to include (optional)"
    )

    # solo_ratio
    parser.add_argument(
        "--min_solo_ratio", type=float, default=None,
        help="Minimum solo_ratio to include (optional)"
    )
    parser.add_argument(
        "--max_solo_ratio", type=float, default=None,
        help="Maximum solo_ratio to include (optional)"
    )

    # length_bias
    parser.add_argument(
        "--min_length_bias", type=float, default=None,
        help="Minimum length_bias to include (optional)"
    )
    parser.add_argument(
        "--max_length_bias", type=float, default=None,
        help="Maximum length_bias to include (optional)"
    )

    return parser.parse_args()


def pearson_ci(r, n, alpha=0.05):
    """
    Compute confidence interval for Pearson r using Fisher z-transform.
    Returns (lower, upper). If n < 4, returns (nan, nan).
    """
    if n < 4 or np.isnan(r):
        return np.nan, np.nan

    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = 1.96  # 95% CI
    lo_z, hi_z = z - z_crit * se, z + z_crit * se
    return np.tanh(lo_z), np.tanh(hi_z)


def main():
    args = parse_args()

    # Columns we will use for plotting
    param_pairs = [
        ("insertion_rate", "deletion_rate"),
        ("insertion_rate", "solo_ratio"),
        ("insertion_rate", "length_bias"),
        ("deletion_rate", "solo_ratio"),
        ("deletion_rate", "length_bias"),
        ("solo_ratio", "length_bias"),
    ]
    all_param_cols = sorted({p for pair in param_pairs for p in pair})

    # Read table (tab-separated by default)
    df_raw = pd.read_csv(args.input, sep="\t")

    # Basic checks
    missing_cols = [c for c in (all_param_cols + [args.metric]) if c not in df_raw.columns]
    if missing_cols:
        raise ValueError(f"Column(s) {missing_cols} not found in input file.")

    # Convert relevant columns to numeric (handles strings like '0.1e-12')
    df = df_raw.copy()
    numeric_cols_to_force = set(all_param_cols + [args.metric])
    for col in numeric_cols_to_force:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # --- Apply global parameter filters before plotting and correlations ---
    mask = np.ones(len(df), dtype=bool)

    # Helper: apply one parameter's min/max filter if provided
    def apply_filter(col_name, min_val, max_val, current_mask):
        col = df[col_name]
        if min_val is not None:
            current_mask &= (col >= min_val)
        if max_val is not None:
            current_mask &= (col <= max_val)
        return current_mask

    mask = apply_filter("insertion_rate",
                        args.min_insertion_rate, args.max_insertion_rate, mask)
    mask = apply_filter("deletion_rate",
                        args.min_deletion_rate, args.max_deletion_rate, mask)
    mask = apply_filter("solo_ratio",
                        args.min_solo_ratio, args.max_solo_ratio, mask)
    mask = apply_filter("length_bias",
                        args.min_length_bias, args.max_length_bias, mask)

    # Subset after filters
    df = df[mask].copy()

    if df.shape[0] < 3:
        raise ValueError(
            "Not enough valid points after filtering to make contour plots (need >= 3)."
        )

    # Metric (already numeric by construction)
    metric_col = args.metric
    df_metric = df[metric_col]
    # === Optimal parameter CI (profile-based, empirical 95%) ===

    # Find best (minimum) composite score
    best_idx = df_metric.idxmin()
    best_row = df.loc[best_idx]
    best_score = df_metric.loc[best_idx]

    # Define 95% threshold as top 5% best scores
    threshold = df_metric.quantile(0.05)

    df_top = df[df_metric <= threshold]

    print("\n=== Optimal parameter estimates (profile-based 95% CI) ===")
    print(f"Best Composite score: {best_score:.6g}")
    print(f"95% threshold (5th percentile): {threshold:.6g}\n")

    for param in ["insertion_rate", "deletion_rate", "solo_ratio", "length_bias"]:
        best_val = best_row[param]
        lo = df_top[param].min()
        hi = df_top[param].max()

        print(
            f"{param:20s}: "
            f"{best_val:.6g} "
            f"[95% CI: {lo:.6g}, {hi:.6g}] "
            f"(n={df_top.shape[0]})"
        )

    print("==========================================================")


    # --- Simple sensitivity summary (correlation) ---
    print("=== Sensitivity summary (Pearson correlation w/ metric; after filtering) ===")
    print("    (R² = % variance in metric explained by a single parameter; linear, univariate)")
    numeric_cols = df.select_dtypes(include=[np.number]).columns

    params_of_interest = [
        "insertion_rate",
        "deletion_rate",
        "solo_ratio",
        "length_bias",
    ]

    for col in params_of_interest:
        valid = ~(df_metric.isna() | df[col].isna())
        n = int(valid.sum())

        if n < 4:
            r = np.nan
            lo, hi = np.nan, np.nan
            r2 = np.nan
            r2_lo, r2_hi = np.nan, np.nan
        else:
            r = df_metric[valid].corr(df[col][valid])
            lo, hi = pearson_ci(r, n)

            # Univariate variance explained
            r2 = r ** 2

            # CI for R² by transforming r CI endpoints
            r2_endpoints = np.array([lo ** 2, hi ** 2], dtype=float)
            r2_lo = np.nanmin(r2_endpoints)
            r2_hi = np.nanmax(r2_endpoints)

        print(
            f"{col:20s}: r = {r: .4f}  "
            f"[95% CI: {lo: .4f}, {hi: .4f}]  "
            f"R² = {r2: .4f} ({r2*100:5.1f}%)  "
            f"[95% CI: {r2_lo: .4f}, {r2_hi: .4f}]  "
            f"(n={n})"
        )

    print("==================================================================")

    # --- Multi-page PDF with contour plots for each parameter pair ---
    figures_made = 0

    with PdfPages(args.output) as pdf:
        for x_param, y_param in param_pairs:
            x = df[x_param]
            y = df[y_param]
            z = df[metric_col]

            # Drop rows with NaNs for this pair + metric
            valid = ~(x.isna() | y.isna() | z.isna())
            x = x[valid].values
            y = y[valid].values
            z = z[valid].values

            if len(x) < 3:
                print(
                    f"Skipping plot for ({x_param}, {y_param}) "
                    f"– not enough valid points (n={len(x)})"
                )
                continue

            triang = Triangulation(x, y)

            fig, ax = plt.subplots(figsize=(6, 5))

            # Filled contour
            contour = ax.tricontourf(
                triang, z,
                levels=args.levels
            )

            # Overlay points: size shrunk to 1/5 (from 10 to 2)
            ax.scatter(x, y, s=2, alpha=0.7, edgecolor="k", linewidth=0.3)

            # Axis labels & title
            ax.set_xlabel(x_param)
            ax.set_ylabel(y_param)
            ax.set_title(
                f"{metric_col} vs {x_param} & {y_param}\n(lower is better)"
            )

            # Colorbar
            cbar = fig.colorbar(contour, ax=ax)
            cbar.set_label(metric_col)

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            figures_made += 1

    if figures_made == 0:
        raise ValueError(
            "No plots were generated (all parameter pairs had < 3 valid points)."
        )

    print(f"Saved {figures_made} contour plots to multi-page PDF: {args.output}")


if __name__ == "__main__":
    main()

