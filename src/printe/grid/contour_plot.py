#!/usr/bin/env python3
"""
Make contour plots of a Composite metric vs. pairs of parameters from an
adaptive surrogate-guided search (guided_search.py / run_loop.py). Samples are
NOT a uniform random grid: they cluster near the genome-size target and the
insertion/deletion rates are drawn in log10 space, so the sensitivity report
treats the data accordingly.

Each panel interpolates the metric over the two plotted parameters (the other
two vary freely across the points), so the smooth surface is an interpolation
of a 4-D cloud onto 2-D rather than a true 2-parameter response surface — read
it qualitatively. Imputed/failed runs (early-terminated sims, penalized to the
worst score) ARE included in the surface so explored-but-failed regions render
in the worst color; they are excluded from the stats and the highlighted points.

Requires:
  - pandas
  - matplotlib
  - scipy
  - scikit-learn  (RandomForest sensitivity)

Example usage:
  python contour_plot.py \
      --input composite_matrix.tsv \
      --metric Composite \
      --output sensitivity_contours.pdf \
      --min_insertion_rate 0.2e-12 --max_insertion_rate 0.1e-11 \
      --min_solo_ratio 30 --max_solo_ratio 70 \
      --min_deletion_rate 0.06e-12 --max_deletion_rate 0.06e-11 \
      --log_params \
      --highlight_top_pct 5
"""

import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")  # for non-interactive / cluster use
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import griddata
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold, cross_val_score


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
        help="Resolution of interpolation grid (default: 200). Higher = finer raster."
    )
    parser.add_argument(
        "--cmap", type=str, default="viridis_r",
        help=(
            "Matplotlib colormap for the metric surface (default: 'viridis_r', so "
            "the worst/most-penalized scores are deep purple and the best are "
            "yellow). Use 'viridis' to flip (worst = yellow)."
        )
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
        "--highlight_linestyle", type=str, default="--",
        help="Line style for the colorbar threshold tick (default: '--' dashed)."
    )
    parser.add_argument(
        "--highlight_marker_color", type=str, default="#D7191C",
        help=(
            "Color for the top-PCT%% scatter points and the colorbar threshold "
            "tick (default: '#D7191C' red, which contrasts with viridis at both "
            "ends)."
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
            "Compute summary statistics (RandomForest sensitivity, best-5% "
            "ranges) on the filtered subset instead of the full dataset. By "
            "default, stats are computed on the full (unfiltered) data so they "
            "remain stable regardless of plot zoom level. Note: imputed/failed "
            "rows are always excluded from the sensitivity model."
        )
    )

    return parser.parse_args()


def imputed_mask(df):
    """
    Boolean mask of synthetic/imputed rows.

    build_composite_matrix.py penalizes simulations that never produced a
    terminal genome (missing FASTA/TSV) with metric = max(real) + 0.1 and
    writes their raw size columns as 0. Those rows are not real genome
    comparisons, so the sensitivity analysis excludes them. Detection uses
    exp_genome_size == 0 when that column is present; otherwise nothing is
    flagged (the input may be a generic metric table).
    """
    if "exp_genome_size" in df.columns:
        return pd.to_numeric(df["exp_genome_size"], errors="coerce").fillna(0) == 0
    return pd.Series(False, index=df.index)


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

    print(f"\n=== Near-optimal parameter ranges – min/max over best 5% of runs ({label}) ===")
    print("    NOTE: a descriptive range of the best-scoring runs, NOT a confidence")
    print("    interval — unlike a CI it tends to WIDEN with more sampling. 'best' is")
    print("    a single stochastic simulation (one draw), so treat it as indicative.")
    print(f"Best {metric_col} score (single run): {best_score:.6g}")
    print(f"Best-5% threshold (5th percentile of {metric_col}): {threshold:.6g}\n")

    for param in params_of_interest:
        best_val = best_row[param]
        lo = df_top[param].min()
        hi = df_top[param].max()

        print(
            f"{disp(param):20s}: "
            f"{best_val:.6g} "
            f"[best-5% range: {lo:.6g}, {hi:.6g}] "
            f"(n={df_top.shape[0]})"
        )

    print("==========================================================")

    # --- Nonlinear sensitivity (RandomForest, cross-validated) ---
    # Composite is a distance-to-reference, hence U-shaped / non-monotone in the
    # rate parameters. Linear (Pearson) and rank (Spearman) statistics are blind
    # to that and badly understate sensitivity (a U-shape has ~zero linear AND
    # ~zero rank correlation). We instead fit the same model family the guided
    # search uses — a RandomForest on log10-rates — and report a cross-validated
    # R² plus permutation importances.
    print(f"\n=== Sensitivity summary – RandomForest, cross-validated ({label}) ===")

    imp = imputed_mask(df)
    n_imp = int(imp.sum())
    df_real = df[~imp]

    feat_names = []
    x_cols = []
    for col in params_of_interest:
        v = pd.to_numeric(df_real[col], errors="coerce").values.astype(float)
        if "rate" in col:  # log10 the rate params (matches plot axes + surrogate)
            v = np.log10(np.where(v > 0, v, np.nan))
            feat_names.append(f"log10({disp(col)})")
        else:
            feat_names.append(disp(col))
        x_cols.append(v)
    X = np.column_stack(x_cols)
    y = pd.to_numeric(df_real[metric_col], errors="coerce").values.astype(float)

    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    X, y = X[ok], y[ok]
    n_rf = len(y)

    if n_imp:
        print(f"  Excluded {n_imp} imputed/failed sim(s); using {n_rf} real simulations.")
    else:
        print(f"  n = {n_rf} simulations.")

    if n_rf < 30:
        print("  Skipped: need >= 30 real simulations for a stable RandomForest.")
    else:
        print("    Importances are permutation-based (mean drop in R² when a feature")
        print("    is shuffled) and describe the *explored* region, which the adaptive")
        print("    search concentrated near the optimum — not the global space.")
        print()

        rf = RandomForestRegressor(
            n_estimators=400, oob_score=True, random_state=0, n_jobs=-1
        )
        cv = KFold(n_splits=min(5, n_rf), shuffle=True, random_state=0)
        cv_r2 = cross_val_score(rf, X, y, cv=cv, scoring="r2")
        rf.fit(X, y)

        print(f"  Joint R² (5-fold CV): {cv_r2.mean(): .4f}  (+/- {cv_r2.std():.4f})")
        print(f"  Joint R² (OOB):       {rf.oob_score_: .4f}")
        print(f"  (n={n_rf}, p={len(feat_names)})")
        print()

        perm = permutation_importance(
            rf, X, y, n_repeats=20, random_state=0, n_jobs=-1, scoring="r2"
        )

        print("  Per-parameter sensitivity (ranked by permutation importance):")
        print(f"    {'parameter':22s} {'perm. importance':>20s}   {'solo CV R²':>11s}")
        for i in np.argsort(perm.importances_mean)[::-1]:
            solo = cross_val_score(
                RandomForestRegressor(n_estimators=300, random_state=0, n_jobs=-1),
                X[:, [i]], y, cv=cv, scoring="r2",
            ).mean()
            print(
                f"    {feat_names[i]:22s} "
                f"{perm.importances_mean[i]: .4f} +/- {perm.importances_std[i]:.4f}   "
                f"{solo: 11.4f}"
            )
        print()
        print("    perm. importance = joint contribution (accounts for redundancy);")
        print("    solo CV R² = variance each parameter explains alone (can be < 0).")

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

    # The plotted SURFACE uses every sampled run, including imputed/failed sims:
    # those were explored and penalized to the worst score (metric = max(real)+0.1),
    # so painting them in the worst color shows which regions were tried and failed.
    # Real rows still drive the STATS, the scatter points, and the top-PCT% markers
    # (imputed values are placeholders, not true genome comparisons).
    df_real = df[~imputed_mask(df)]
    df_filt_real = df_filtered[~imputed_mask(df_filtered)]
    n_imp_plot = len(df) - len(df_real)
    if n_imp_plot:
        print(f"\n(Surface shows {n_imp_plot} imputed/failed rows at the worst "
              f"color; stats and highlighted points use {len(df_real)} real sims.)")

    # --- Determine highlight threshold from the REAL dataset (stable) ---
    highlight_enabled = args.highlight_top_pct is not None
    if highlight_enabled:
        if not (0 < args.highlight_top_pct < 100):
            raise ValueError("--highlight_top_pct must be between 0 and 100 (exclusive).")
        # Threshold computed on the real data so "best X%" is independent of zoom
        # and of how many sims failed.
        highlight_threshold = df_real[metric_col].quantile(args.highlight_top_pct / 100.0)
        n_best = int((df_real[metric_col] <= highlight_threshold).sum())
        print(
            f"\n=== Highlighting best {args.highlight_top_pct:.3g}% region on plots ===\n"
            f"Highlight threshold ({args.highlight_top_pct:.3g}th percentile of {metric_col}): "
            f"{highlight_threshold:.6g}  (n={n_best} real simulations)\n"
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

    def to_plot_coord(val, use_log, fallback):
        if val is None:
            return fallback
        return np.log10(max(val, args.log_epsilon)) if use_log else val

    def get_coords(frame, x_param, y_param, use_log_x, use_log_y):
        """Numeric, validity-filtered, optionally-log10 (x, y, z) arrays."""
        xs = pd.to_numeric(frame[x_param], errors="coerce")
        ys = pd.to_numeric(frame[y_param], errors="coerce")
        zs = pd.to_numeric(frame[metric_col], errors="coerce")
        good = ~(xs.isna() | ys.isna() | zs.isna())
        if use_log_x:
            good &= xs > 0
        if use_log_y:
            good &= ys > 0
        xv = xs[good].values.astype(float)
        yv = ys[good].values.astype(float)
        zv = zs[good].values.astype(float)
        if use_log_x:
            xv = np.log10(xv)
        if use_log_y:
            yv = np.log10(yv)
        return xv, yv, zv

    with PdfPages(args.output) as pdf:
        for x_param, y_param in param_pairs:
            use_log_x = args.log_x or (args.log_params and ("rate" in x_param))
            use_log_y = args.log_y or (args.log_params and ("rate" in y_param))

            # All sampled rows (incl. imputed/failed at their worst-color penalty)
            # build the interpolated surface, so explored-but-failed regions read as
            # worst. Only real, range-filtered rows are scattered as points on top.
            x_surf, y_surf, z_surf = get_coords(df, x_param, y_param,
                                                use_log_x, use_log_y)
            x_filt, y_filt, z_filt = get_coords(df_filt_real, x_param, y_param,
                                                use_log_x, use_log_y)

            if len(x_surf) < 3:
                print(f"Skipping ({x_param}, {y_param}) – only {len(x_surf)} points.")
                continue

            # ---- View window (filter bounds, else real-data range) ----
            x_lo_bound, x_hi_bound = param_filter_bounds[x_param]
            y_lo_bound, y_hi_bound = param_filter_bounds[y_param]

            x_view_lo = to_plot_coord(x_lo_bound, use_log_x, float(x_surf.min()))
            x_view_hi = to_plot_coord(x_hi_bound, use_log_x, float(x_surf.max()))
            y_view_lo = to_plot_coord(y_lo_bound, use_log_y, float(y_surf.min()))
            y_view_hi = to_plot_coord(y_hi_bound, use_log_y, float(y_surf.max()))

            # ---- Interpolate onto a grid spanning the view window ----
            grid_res = args.grid_res
            xi = np.linspace(x_view_lo, x_view_hi, grid_res)
            yi = np.linspace(y_view_lo, y_view_hi, grid_res)
            Xi, Yi = np.meshgrid(xi, yi)

            pts = np.column_stack((x_surf, y_surf))
            Zi = griddata(pts, z_surf, (Xi, Yi), method="linear")
            nan_mask = np.isnan(Zi)
            if nan_mask.any():
                Zi[nan_mask] = griddata(
                    pts, z_surf, (Xi[nan_mask], Yi[nan_mask]), method="nearest"
                )

            if highlight_enabled:
                best_mask = z_filt <= highlight_threshold
                n_best_pair = int(best_mask.sum())
            else:
                best_mask = None
                n_best_pair = 0

            # ---- Plot ----
            fig, ax = plt.subplots(figsize=(6, 5))

            z_min, z_max = float(Zi.min()), float(Zi.max())
            if z_min == z_max:
                z_max = z_min + 1e-6
            if args.color_ceiling is not None and args.color_ceiling > z_min:
                z_max_plot = args.color_ceiling
            else:
                z_max_plot = z_max

            contour_levels = np.linspace(z_min, z_max_plot, args.levels)
            contour = ax.contourf(Xi, Yi, Zi, levels=contour_levels,
                                  cmap=args.cmap, extend="both")

            # Scatter real, range-filtered points
            if not args.no_points:
                ax.scatter(
                    x_filt, y_filt,
                    s=2, alpha=0.7, color="0.35",
                    edgecolor="k", linewidth=0.3, zorder=3,
                )

            if highlight_enabled and (best_mask is not None) and (n_best_pair >= 1):
                ax.scatter(
                    x_filt[best_mask], y_filt[best_mask],
                    s=12, alpha=0.95,
                    facecolor=args.highlight_marker_color,
                    edgecolor="k", linewidth=0.4, zorder=4,
                    label=f"Top {args.highlight_top_pct:.3g}%",
                )
                ax.legend(loc="best", fontsize=8, framealpha=0.85,
                          edgecolor="0.6", fancybox=False)

            # No title — caption belongs in the manuscript (lower metric = better).
            ax.set_xlabel(f"log10({x_param})" if use_log_x else x_param)
            ax.set_ylabel(f"log10({y_param})" if use_log_y else y_param)

            cbar = fig.colorbar(contour, ax=ax)
            cbar_label = metric_col
            if args.color_ceiling is not None and args.color_ceiling > z_min:
                cbar_label += f" (ceiling={args.color_ceiling:g})"
            cbar.set_label(cbar_label)

            # Mark the top-PCT% threshold on the colorbar (same color as the points)
            if highlight_enabled and highlight_threshold is not None:
                if z_min < highlight_threshold < z_max_plot:
                    cbar.ax.axhline(
                        y=highlight_threshold, color=args.highlight_marker_color,
                        linewidth=1.5, linestyle=args.highlight_linestyle,
                    )

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            figures_made += 1

    if figures_made == 0:
        raise ValueError(
            "No plots were generated (no parameter pair had >= 3 real points)."
        )

    print(f"\nSaved {figures_made} contour plots to multi-page PDF: {args.output}")


if __name__ == "__main__":
    main()
