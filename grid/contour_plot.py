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
      --min_deletion_rate 0.06e-12 --max_deletion_rate 0.06e-11 \
      --log_params \
      --highlight_top_pct 5
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import griddata


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
    parser.add_argument(
        "--grid_res", type=int, default=200,
        help="Resolution of interpolation grid (default: 200). Higher = smoother but slower."
    )

    # --- Parameter range filters ---
    parser.add_argument(
        "--min_insertion_rate", type=float, default=None,
        help="Minimum insertion_rate to include (optional)"
    )
    parser.add_argument(
        "--max_insertion_rate", type=float, default=None,
        help="Maximum insertion_rate to include (optional)"
    )

    parser.add_argument(
        "--min_deletion_rate", type=float, default=None,
        help="Minimum deletion_rate to include (optional)"
    )
    parser.add_argument(
        "--max_deletion_rate", type=float, default=None,
        help="Maximum deletion_rate to include (optional)"
    )

    parser.add_argument(
        "--min_solo_ratio", type=float, default=None,
        help="Minimum solo_ratio to include (optional)"
    )
    parser.add_argument(
        "--max_solo_ratio", type=float, default=None,
        help="Maximum solo_ratio to include (optional)"
    )

    parser.add_argument(
        "--min_length_bias", type=float, default=None,
        help="Minimum length_bias to include (optional)"
    )
    parser.add_argument(
        "--max_length_bias", type=float, default=None,
        help="Maximum length_bias to include (optional)"
    )

    parser.add_argument(
        "--no_points", action="store_true",
        help="Do not overlay individual sample points on contour plots."
    )

    # --- Log scaling options ---
    parser.add_argument(
        "--log_x", action="store_true",
        help="Plot x-axes on log10 scale (for rate-like params)."
    )
    parser.add_argument(
        "--log_y", action="store_true",
        help="Plot y-axes on log10 scale (for rate-like params)."
    )
    parser.add_argument(
        "--log_params", action="store_true",
        help="Automatically use log10 axes for any parameter containing 'rate'."
    )
    parser.add_argument(
        "--log_epsilon", type=float, default=1e-300,
        help="Small value to guard log10(0); values <= 0 are dropped anyway. Default: 1e-300"
    )

    # --- Highlight best region ---
    parser.add_argument(
        "--highlight_top_pct", type=float, default=None,
        help=(
            "Highlight region corresponding to best (lowest-metric) PCT%% of runs "
            "(e.g., 5 for best 5%%). If omitted, no highlighting."
        )
    )
    parser.add_argument(
        "--highlight_hatch", type=str, default="////",
        help="Hatch pattern for highlighted region overlay (default: '////')."
    )
    parser.add_argument(
        "--highlight_alpha", type=float, default=0.0,
        help=(
            "Fill alpha for highlight overlay (default: 0.0). "
            "Keep 0.0 if you only want hatch without fill."
        )
    )
    parser.add_argument(
        "--highlight_linewidth", type=float, default=1.2,
        help="Line width for highlight boundary (default: 1.2)."
    )

    # --- Report scope control ---
    parser.add_argument(
        "--report_on_filtered", action="store_true",
        help=(
            "Compute summary statistics (correlations, optimal CIs) on the "
            "filtered subset instead of the full dataset. By default, stats "
            "are computed on the full (unfiltered) data so they remain stable "
            "regardless of plot zoom level."
        )
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


def apply_filter(df, col_name, min_val, max_val, current_mask):
    col = df[col_name]
    if min_val is not None:
        current_mask &= (col >= min_val)
    if max_val is not None:
        current_mask &= (col <= max_val)
    return current_mask


def print_report(df, metric_col, label, params_of_interest):
    """Print parameter ranges, optimal estimates, and correlations."""
    df_metric = df[metric_col]

    print(f"\n=== Parameter ranges ({label}) ===")
    for param in params_of_interest:
        vals = df[param].dropna()
        if len(vals) == 0:
            print(f"{param:20s}: no valid values")
        else:
            print(
                f"{param:20s}: "
                f"min = {vals.min():.6g}, "
                f"max = {vals.max():.6g}, "
                f"(n={len(vals)})"
            )
    print("========================================")

    # === Optimal parameter CI (profile-based, empirical 95%) ===
    best_idx = df_metric.idxmin()
    best_row = df.loc[best_idx]
    best_score = df_metric.loc[best_idx]

    threshold = df_metric.quantile(0.05)
    df_top = df[df_metric <= threshold]

    print(f"\n=== Optimal parameter estimates – profile-based 95% CI ({label}) ===")
    print(f"Best {metric_col} score: {best_score:.6g}")
    print(f"95% threshold (5th percentile): {threshold:.6g}\n")

    for param in params_of_interest:
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

    # --- Sensitivity summary (correlation) ---
    print(f"=== Sensitivity summary – Pearson r ({label}) ===")
    print("    (R² = % variance in metric explained by a single parameter; linear, univariate)")

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

            r2 = r ** 2
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
    params_of_interest = [
        "insertion_rate",
        "deletion_rate",
        "solo_ratio",
        "length_bias",
    ]

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

    # =========================================================================
    # Report on FULL dataset (before range filters) – stable across zoom levels
    # =========================================================================
    metric_col = args.metric

    if not args.report_on_filtered:
        print_report(df, metric_col, "full dataset – unfiltered", params_of_interest)

    # --- Apply global parameter filters (for plotting) ---
    mask = np.ones(len(df), dtype=bool)
    mask = apply_filter(df, "insertion_rate",
                        args.min_insertion_rate, args.max_insertion_rate, mask)
    mask = apply_filter(df, "deletion_rate",
                        args.min_deletion_rate, args.max_deletion_rate, mask)
    mask = apply_filter(df, "solo_ratio",
                        args.min_solo_ratio, args.max_solo_ratio, mask)
    mask = apply_filter(df, "length_bias",
                        args.min_length_bias, args.max_length_bias, mask)

    # Subset after filters
    df_filtered = df[mask].copy()

    if df_filtered.shape[0] < 3:
        raise ValueError(
            "Not enough valid points after filtering to make contour plots (need >= 3)."
        )

    # If user explicitly wants stats on the filtered subset, print those too
    if args.report_on_filtered:
        print_report(df_filtered, metric_col, "filtered subset", params_of_interest)
    else:
        # Still print the filtered ranges so the user knows the plot extent
        print(f"\n--- Plot range after filtering: {df_filtered.shape[0]} points ---")
        for param in params_of_interest:
            vals = df_filtered[param].dropna()
            if len(vals):
                print(f"  {param:20s}: [{vals.min():.6g}, {vals.max():.6g}]")
        print("---------------------------------------------------")

    # --- Determine highlight threshold from FULL dataset (stable) ---
    highlight_enabled = args.highlight_top_pct is not None
    if highlight_enabled:
        if not (0 < args.highlight_top_pct < 100):
            raise ValueError("--highlight_top_pct must be between 0 and 100 (exclusive).")
        # Threshold always computed on full data so it means the same thing
        # regardless of zoom
        highlight_threshold = df[metric_col].quantile(args.highlight_top_pct / 100.0)
        n_best = int((df[metric_col] <= highlight_threshold).sum())
        print(
            f"\n=== Highlighting best {args.highlight_top_pct:.3g}% region on plots ===\n"
            f"Highlight threshold ({args.highlight_top_pct:.3g}th percentile of {metric_col}): "
            f"{highlight_threshold:.6g}  (n={n_best} points in full dataset)\n"
            "==============================================================="
        )
    else:
        highlight_threshold = None

    # --- Multi-page PDF with contour plots for each parameter pair ---
    figures_made = 0

    # Collect filter bounds per parameter (used to define the view window)
    param_filter_bounds = {
        "insertion_rate": (args.min_insertion_rate, args.max_insertion_rate),
        "deletion_rate":  (args.min_deletion_rate,  args.max_deletion_rate),
        "solo_ratio":     (args.min_solo_ratio,     args.max_solo_ratio),
        "length_bias":    (args.min_length_bias,    args.max_length_bias),
    }

    with PdfPages(args.output) as pdf:
        for x_param, y_param in param_pairs:
            # ---- FULL dataset for interpolation ----
            x_full_s = df[x_param]
            y_full_s = df[y_param]
            z_full_s = df[metric_col]

            use_log_x = args.log_x or (args.log_params and ("rate" in x_param))
            use_log_y = args.log_y or (args.log_params and ("rate" in y_param))

            valid_full = ~(x_full_s.isna() | y_full_s.isna() | z_full_s.isna())
            if use_log_x:
                valid_full &= (x_full_s > 0)
            if use_log_y:
                valid_full &= (y_full_s > 0)

            x_full = x_full_s[valid_full].values.astype(float)
            y_full = y_full_s[valid_full].values.astype(float)
            z_full = z_full_s[valid_full].values.astype(float)

            if use_log_x:
                x_full = np.log10(np.maximum(x_full, args.log_epsilon))
            if use_log_y:
                y_full = np.log10(np.maximum(y_full, args.log_epsilon))

            if len(x_full) < 3:
                print(
                    f"Skipping plot for ({x_param}, {y_param}) "
                    f"– not enough valid points in full data (n={len(x_full)})"
                )
                continue

            # ---- View window from filter bounds ----
            x_lo_bound, x_hi_bound = param_filter_bounds[x_param]
            y_lo_bound, y_hi_bound = param_filter_bounds[y_param]

            def to_plot_coord(val, use_log, fallback):
                if val is None:
                    return fallback
                return np.log10(max(val, args.log_epsilon)) if use_log else val

            x_view_lo = to_plot_coord(x_lo_bound, use_log_x, x_full.min())
            x_view_hi = to_plot_coord(x_hi_bound, use_log_x, x_full.max())
            y_view_lo = to_plot_coord(y_lo_bound, use_log_y, y_full.min())
            y_view_hi = to_plot_coord(y_hi_bound, use_log_y, y_full.max())

            # ---- Interpolate FULL data onto grid spanning the VIEW window ----
            grid_res = args.grid_res
            xi = np.linspace(x_view_lo, x_view_hi, grid_res)
            yi = np.linspace(y_view_lo, y_view_hi, grid_res)
            Xi, Yi = np.meshgrid(xi, yi)

            pts_full = np.column_stack((x_full, y_full))

            Zi = griddata(pts_full, z_full, (Xi, Yi), method="linear")

            nan_mask = np.isnan(Zi)
            if nan_mask.any():
                Zi[nan_mask] = griddata(
                    pts_full, z_full,
                    (Xi[nan_mask], Yi[nan_mask]),
                    method="nearest"
                )

            # ---- Filtered points for scatter overlay ----
            x_filt_s = df_filtered[x_param]
            y_filt_s = df_filtered[y_param]
            z_filt_s = df_filtered[metric_col]

            valid_filt = ~(x_filt_s.isna() | y_filt_s.isna() | z_filt_s.isna())
            if use_log_x:
                valid_filt &= (x_filt_s > 0)
            if use_log_y:
                valid_filt &= (y_filt_s > 0)

            x_filt = x_filt_s[valid_filt].values.astype(float)
            y_filt = y_filt_s[valid_filt].values.astype(float)
            z_filt = z_filt_s[valid_filt].values.astype(float)

            if use_log_x:
                x_filt = np.log10(np.maximum(x_filt, args.log_epsilon))
            if use_log_y:
                y_filt = np.log10(np.maximum(y_filt, args.log_epsilon))

            if highlight_enabled:
                best_mask = z_filt <= highlight_threshold
                n_best_pair = int(best_mask.sum())
            else:
                best_mask = None
                n_best_pair = 0

            # ---- Plot ----
            fig, ax = plt.subplots(figsize=(6, 5))

            # Color scale spans the view window's interpolated range
            z_min, z_max = float(Zi.min()), float(Zi.max())
            if z_min == z_max:
                z_max = z_min + 1e-6
            contour_levels = np.linspace(z_min, z_max, args.levels)

            contour = ax.contourf(
                Xi, Yi, Zi,
                levels=contour_levels,
                extend="both",
            )

            # Highlight best region
            if highlight_enabled and highlight_threshold is not None:
                # Only draw if the threshold falls within the view's z range;
                # otherwise the entire view is already "best" (or none of it is).
                if n_best_pair >= 1 and z_min < highlight_threshold < z_max:
                    ax.contourf(
                        Xi, Yi, Zi,
                        levels=[z_min, highlight_threshold],
                        alpha=args.highlight_alpha,
                        hatches=[args.highlight_hatch]
                    )
                    ax.contour(
                        Xi, Yi, Zi,
                        levels=[highlight_threshold],
                        linewidths=args.highlight_linewidth
                    )

            # Scatter filtered points only
            if not args.no_points:
                ax.scatter(
                    x_filt, y_filt,
                    s=2, alpha=0.7,
                    edgecolor="k", linewidth=0.3
                )

            if highlight_enabled and (best_mask is not None) and (n_best_pair >= 1):
                ax.scatter(
                    x_filt[best_mask], y_filt[best_mask],
                    s=8, alpha=0.9,
                    edgecolor="k", linewidth=0.4
                )

            ax.set_title(f"{metric_col} vs {x_param} & {y_param}\n(lower is better)")
            ax.set_xlabel(f"log10({x_param})" if use_log_x else x_param)
            ax.set_ylabel(f"log10({y_param})" if use_log_y else y_param)

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

    print(f"\nSaved {figures_made} contour plots to multi-page PDF: {args.output}")


if __name__ == "__main__":
    main()
