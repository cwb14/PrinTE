#!/usr/bin/env python3
"""Exponential trajectory projection + hard-breach decision for PrinTE.

Used by PrinTE.sh's per-iteration loop. Given the observed genome-size
trajectory (including the iteration-0 anchor), it (1) projects the terminal
size via log-linear regression, clamped to a sane range, and (2) decides
whether the run has hard-breached a size bound while moving further out of it.

On a breach it prints a standard `Projection (...)` line plus the matching
bound line (the exact strings grid/guided_search.py harvests) and exits 10
(too_large) or 11 (too_small). Otherwise it prints nothing and exits 0.
"""
import argparse
import math
import sys


def fit_project(xs: list[float], ys: list[float], target: float,
                clamp_lo: float, clamp_hi: float) -> tuple[int, int, int, float]:
    """Fit log(y)=a+b*x over >=2 points, project to `target`, clamp the point.

    Returns (point, lo, hi, r2). lo/hi are a nominal +-0 interval around the
    clamped point (the calibrated PI is not meaningful at <min_n points; the
    breach decision does not use it). r2 is the log-space fit quality.
    """
    n = len(xs)
    log_ys = [math.log(y) for y in ys]
    sum_x = sum(xs)
    sum_ly = sum(log_ys)
    sum_xly = sum(x * ly for x, ly in zip(xs, log_ys))
    sum_x2 = sum(x * x for x in xs)
    denom = n * sum_x2 - sum_x * sum_x

    # Fit coefficients once (degenerate denom == 0 => flat fit at the mean).
    if denom == 0:
        b = 0.0
        a = sum_ly / n
    else:
        b = (n * sum_xly - sum_x * sum_ly) / denom
        a = (sum_ly - b * sum_x) / n
    proj = math.exp(a + b * target)
    point = int(max(clamp_lo, min(clamp_hi, proj)))

    # r2 in log space (1.0 when there is no variance to explain, e.g. 2 points)
    mean_ly = sum_ly / n
    ss_res = sum((ly - (a + b * x)) ** 2 for x, ly in zip(xs, log_ys))
    ss_tot = sum((ly - mean_ly) ** 2 for ly in log_ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return point, point, point, r2


def _emit_projection_line(i, n, target, point, lo, hi, r2):
    frac_pct = int(100 * n / target) if target else 100
    print(f"Projection (iter {i}, n={n}, {frac_pct}% observed, R²={r2:.4f}): "
          f"point={point}, empirical 95% PI=[{lo}..{hi}]")


def main():
    p = argparse.ArgumentParser(description="PrinTE projection + hard-breach check")
    p.add_argument("--iters", required=True, help="comma-separated iteration indices (incl. 0 anchor)")
    p.add_argument("--sizes", required=True, help="comma-separated byte sizes, aligned with --iters")
    p.add_argument("--target", type=int, required=True, help="terminal iteration count")
    p.add_argument("--mxgs", type=int, default=None)
    p.add_argument("--mngs", type=int, default=None)
    p.add_argument("--clamp-lo", type=int, required=True)
    p.add_argument("--clamp-hi", type=int, required=True)
    args = p.parse_args()

    xs = [float(x) for x in args.iters.split(",")]
    ys = [float(y) for y in args.sizes.split(",")]
    if len(xs) < 2 or len(xs) != len(ys):
        sys.exit(0)  # not enough data for a direction; let the caller continue

    last = ys[-1]
    prev = ys[-2]
    i = int(xs[-1])
    n = len(xs)
    point, lo, hi, r2 = fit_project(xs, ys, args.target, args.clamp_lo, args.clamp_hi)

    if args.mxgs is not None and last > args.mxgs and last > prev:
        _emit_projection_line(i, n, args.target, point, lo, hi, r2)
        print(f"Empirical PI lower bound ({point} bytes) exceeds maximum ({args.mxgs} bytes).")
        print(f"Confidently above target (hard breach). Stopping at iteration {i}.")
        sys.exit(10)

    if args.mngs is not None and last < args.mngs and last < prev:
        _emit_projection_line(i, n, args.target, point, lo, hi, r2)
        print(f"Empirical PI upper bound ({point} bytes) is below minimum ({args.mngs} bytes).")
        print(f"Confidently below target (hard breach). Stopping at iteration {i}.")
        sys.exit(11)

    sys.exit(0)


if __name__ == "__main__":
    main()
