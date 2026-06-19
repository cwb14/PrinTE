#!/usr/bin/env python3
"""
Make contour plots of Composite metric vs. pairs of parameters
from a random grid search sample.

Requires:
  - pandas
  - matplotlib
  - scipy
  - scikit-learn  (for multiple R²)

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
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression


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
        "--highlight_color", type=str, default="#FFD700",
        help=(
            "Color for the highlight boundary and fill "
            "(default: '#FFD700' gold). Use any matplotlib color string."
        )
    )
    parser.add_argument(
        "--highlight_linewidth", type=float, default=2.0,
        help="Line width for highlight boundary (default: 2.0)."
    )
    parser.add_argument(
        "--highlight_linestyle", type=str, default="--",
        help="Line style for highlight boundary (default: '--' dashed)."
    )
    parser.add_argument(
        "--highlight_fill_alpha", type=float, default=0.12,
        help=(
            "Fill alpha for highlight region overlay (default: 0.12). "
            "Set to 0.0 to disable the fill entirely."
        )
    )
    parser.add_argument(
        "--highlight_marker_color", type=str, default="#FFD700",
        help=(
            "Color for highlighted scatter points "
            "(default: '#FFD700' gold)."
        )
    )

    # --- Color scale ceiling ---
    parser.add_argument(
        "--color_ceiling", type=float, default=None,
        help=(
            "Cap the color scale at this composite score value. "
            "The full color gradient is compressed into [min, ceiling]; "
            "any data point above the ceiling is rendered at the top color (yellow). "
            "Useful when a small fraction of high-scoring runs would otherwise "
            "dominate the scale and wash out variation near the optimum."
        )
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


def r2_ci_from_r_ci(r_lo, r_hi):
    """
    Compute R² CI from the endpoints of a Pearson r CI.
    If the r CI spans zero, the minimum possible R² is 0.
    """
    if np.isnan(r_lo) or np.isnan(r_hi):
        return np.nan, np.nan

    # If the CI for r contains zero, R² can be as low as 0
    if r_lo <= 0 <= r_hi:
        r2_lo = 0.0
    else:
        r2_lo = min(r_lo ** 2, r_hi ** 2)

    r2_hi = max(r_lo ** 2, r_hi ** 2)
    return r2_lo, r2_hi


def apply_filter(df, col_name, min_val, max_val, current_mask):
    col = df[col_name]
    if min_val is not None:
        current_mask &= (col >= min_val)
    if max_val is not None:
        current_mask &= (col <= max_val)
    return current_mask


# Display labels for printed reports (column names stay unchanged internally)
DISPLAY_NAMES = {"solo_ratio": "solo_rate"}


def disp(name):
    """Return the human-facing label for a column name."""
    return DISPLAY_NAMES.get(name, name)


def print_report(df, metric_col, label, params_of_interest):
    """Print parameter ranges, optimal estimates, correlations, and multiple R²."""
    df_metric = df[metric_col]

    print(f"\n=== Parameter ranges ({label}) ===")
    for param in params_of_interest:
        vals = df[param].dropna()
        if len(vals) == 0:
            print(f"{disp(param):20s}: no valid values")
        else:
            print(
                f"{disp(param):20s}: "
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
            f"{disp(param):20s}: "
            f"{best_val:.6g} "
            f"[95% CI: {lo:.6g}, {hi:.6g}] "
            f"(n={df_top.shape[0]})"
        )

    print("==========================================================")

    # --- Sensitivity summary (Pearson + Spearman) ---
    print(f"\n=== Sensitivity summary – Pearson & Spearman ({label}) ===")
    print(
        "    Pearson r  = linear correlation "
        "(R² = fraction of variance explained, linear & univariate)"
    )
    print(
        "    Spearman ρ = rank correlation "
        "(detects monotonic nonlinear relationships)"
    )
    print()

    for col in params_of_interest:
        valid = ~(df_metric.isna() | df[col].isna())
        n = int(valid.sum())

        if n < 4:
            r = np.nan
            r_lo, r_hi = np.nan, np.nan
            r2 = np.nan
            r2_lo, r2_hi = np.nan, np.nan
            rho = np.nan
            sp_pval = np.nan
        else:
            metric_valid = df_metric[valid].values
            col_valid = df[col][valid].values

            # Pearson
            r = float(pd.Series(metric_valid).corr(pd.Series(col_valid)))
            r_lo, r_hi = pearson_ci(r, n)
            r2 = r ** 2
            r2_lo, r2_hi = r2_ci_from_r_ci(r_lo, r_hi)

            # Spearman
            rho, sp_pval = spearmanr(col_valid, metric_valid)

        print(
            f"  {disp(col):20s}:  "
            f"Pearson r = {r: .4f}  "
            f"[95% CI: {r_lo: .4f}, {r_hi: .4f}]  "
            f"R² = {r2: .4f} ({r2*100:5.1f}%)  "
            f"[95% CI: {r2_lo: .4f}, {r2_hi: .4f}]"
        )
        print(
            f"  {'':20s}   "
            f"Spearman ρ = {rho: .4f}  "
            f"(p = {sp_pval:.2e})  "
            f"(n={n})"
        )

    # --- Multiple R² (joint linear model) ---
    valid_all = np.ones(len(df), dtype=bool)
    for col in params_of_interest:
        valid_all &= ~df[col].isna()
    valid_all &= ~df_metric.isna()
    n_joint = int(valid_all.sum())

    if n_joint >= len(params_of_interest) + 1:
        X = df.loc[valid_all, params_of_interest].values
        y = df_metric[valid_all].values

        model = LinearRegression().fit(X, y)
        r2_multi = model.score(X, y)

        # Adjusted R²
        p = X.shape[1]
        r2_adj = 1.0 - (1.0 - r2_multi) * (n_joint - 1) / (n_joint - p - 1)

        print(f"\n=== Multiple linear regression ({label}) ===")
        print(f"  Predictors: {', '.join(disp(p) for p in params_of_interest)}")
        print(f"  R²          = {r2_multi:.4f} ({r2_multi*100:.1f}%)")
        print(f"  Adjusted R² = {r2_adj:.4f} ({r2_adj*100:.1f}%)")
        print(f"  (n={n_joint}, p={p})")
        print()
        print("  Coefficients:")
        for name, coef in zip(params_of_interest, model.coef_):
            print(f"    {disp(name):20s}: {coef: .6g}")
        print(f"    {'(intercept)':20s}: {model.intercept_: .6g}")
    else:
        print(
            f"\n=== Multiple linear regression ({label}) ===\n"
            f"  Skipped: not enough complete cases (n={n_joint})\n"
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
                print(f"  {disp(param):20s}: [{vals.min():.6g}, {vals.max():.6g}]")
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

            # Color scale spans the view window's interpolated range,
            # optionally capped at --color_ceiling.
            z_min, z_max = float(Zi.min()), float(Zi.max())
            if z_min == z_max:
                z_max = z_min + 1e-6

            if args.color_ceiling is not None and args.color_ceiling > z_min:
                z_max_plot = args.color_ceiling
            else:
                z_max_plot = z_max

            contour_levels = np.linspace(z_min, z_max_plot, args.levels)

            contour = ax.contourf(
                Xi, Yi, Zi,
                levels=contour_levels,
                extend="both",
            )

            # Highlight best region — publication-quality contrasting style
            if highlight_enabled and highlight_threshold is not None:
                if n_best_pair >= 1 and z_min < highlight_threshold < z_max:
                    # Subtle semi-transparent fill to tint the best region
                    if args.highlight_fill_alpha > 0:
                        ax.contourf(
                            Xi, Yi, Zi,
                            levels=[z_min, highlight_threshold],
                            colors=[args.highlight_color],
                            alpha=args.highlight_fill_alpha,
                        )

                    # Bold contrasting boundary line
                    cs_line = ax.contour(
                        Xi, Yi, Zi,
                        levels=[highlight_threshold],
                        colors=[args.highlight_color],
                        linewidths=args.highlight_linewidth,
                        linestyles=args.highlight_linestyle,
                    )

            # Scatter filtered points only
            if not args.no_points:
                ax.scatter(
                    x_filt, y_filt,
                    s=2, alpha=0.7, color="0.35",
                    edgecolor="k", linewidth=0.3,
                    zorder=3,
                )

            if highlight_enabled and (best_mask is not None) and (n_best_pair >= 1):
                ax.scatter(
                    x_filt[best_mask], y_filt[best_mask],
                    s=12, alpha=0.95,
                    facecolor=args.highlight_marker_color,
                    edgecolor="k", linewidth=0.4,
                    zorder=4,
                    label=f"Top {args.highlight_top_pct:.3g}%",
                )
                ax.legend(
                    loc="best", fontsize=8, framealpha=0.85,
                    edgecolor="0.6", fancybox=False,
                )

            ax.set_title(f"{metric_col} vs {x_param} & {y_param}\n(lower is better)")
            ax.set_xlabel(f"log10({x_param})" if use_log_x else x_param)
            ax.set_ylabel(f"log10({y_param})" if use_log_y else y_param)

            cbar = fig.colorbar(contour, ax=ax)
            cbar_label = metric_col
            if args.color_ceiling is not None and args.color_ceiling > z_min:
                cbar_label += f" (ceiling={args.color_ceiling:g})"
            cbar.set_label(cbar_label)

            # Draw the highlight threshold on the colorbar for reference
            if highlight_enabled and highlight_threshold is not None:
                if z_min < highlight_threshold < z_max_plot:
                    cbar.ax.axhline(
                        y=highlight_threshold, color=args.highlight_color,
                        linewidth=1.5, linestyle=args.highlight_linestyle,
                    )

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
