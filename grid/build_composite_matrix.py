#!/usr/bin/env python3
import os
import re
import glob
import csv
import sys
import subprocess
import argparse
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

DIR_REGEX = re.compile(
    r"^insertion_rates_(?P<insertion_rate>[^_]+)"
    r"_deletion_rates_(?P<deletion_rate>[^_]+)"
    r"_solo_ratio_(?P<solo_ratio>[^_]+)"
    r"_length_bias_(?P<length_bias>[^_]+)$"
)

# Allow floats and scientific notation
FLOAT_RE = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
INT_RE = r"(\d+)"

# --- Patterns for raw values ---
RAW_PATTERNS = {
    "ref_sample_count":   re.compile(r"Sample count:\s+ref=" + INT_RE),
    "exp_sample_count":   re.compile(r"Sample count:\s+ref=\d+,\s+exp=" + INT_RE),
    "ref_cumulative_length": re.compile(r"Cumulative length \(bp\):\s+ref=" + INT_RE),
    "exp_cumulative_length": re.compile(r"Cumulative length \(bp\):\s+ref=\d+,\s+exp=" + INT_RE),
    "ref_genome_size":    re.compile(r"Genome size \(bp\):\s+ref=" + INT_RE),
    "exp_genome_size":    re.compile(r"Genome size \(bp\):\s+ref=\d+,\s+exp=" + INT_RE),
    "raw_distribution":   re.compile(
        r"(?:RMS|Kolmogorov|Wasserstein)[^\n]*?distance[^\n]*?:\s*" + FLOAT_RE
    ),
}

# --- Patterns for per-metric distances ---
METRIC_PATTERNS = {
    "d_sample_count":      re.compile(r"d_sample_count.*?"      + FLOAT_RE),
    "d_cumulative_length": re.compile(r"d_cumulative_length.*?" + FLOAT_RE),
    "d_genome_size":       re.compile(r"d_genome_size.*?"       + FLOAT_RE),
    "d_distribution":      re.compile(r"d_distribution.*?"      + FLOAT_RE),
    "Composite":           re.compile(r"Composite\s*=\s*Σ[^=]*=\s*" + FLOAT_RE),
}

# --- Patterns for weights ---
WEIGHT_PATTERNS = {
    "alpha_sample_count":      re.compile(r"α_sample_count\s*=\s*"      + FLOAT_RE),
    "alpha_cumulative_length": re.compile(r"α_cumulative_length\s*=\s*" + FLOAT_RE),
    "alpha_genome_size":       re.compile(r"α_genome_size\s*=\s*"       + FLOAT_RE),
    "alpha_distribution":      re.compile(r"α_distribution\s*=\s*"      + FLOAT_RE),
}


def parse_compare_output(stdout_text):
    """Extract numeric metrics from compare_genomes2.py output."""
    metrics = {}

    # Raw values
    for key, pattern in RAW_PATTERNS.items():
        m = pattern.search(stdout_text)
        if not m:
            return None
        try:
            metrics[key] = float(m.group(1))
        except ValueError:
            return None

    # Per-metric distances
    for key, pattern in METRIC_PATTERNS.items():
        m = pattern.search(stdout_text)
        if not m:
            return None
        try:
            metrics[key] = float(m.group(1))
        except ValueError:
            return None

    # Weights
    for key, pattern in WEIGHT_PATTERNS.items():
        m = pattern.search(stdout_text)
        if not m:
            return None
        try:
            metrics[key] = float(m.group(1))
        except ValueError:
            return None

    return metrics


# ---------------------------------------------------------------------------
# Progress tracker
# ---------------------------------------------------------------------------
class ProgressTracker:
    """Thread-safe progress counter with ETA."""

    def __init__(self, total):
        self.total = total
        self.completed = 0
        self.failed = 0
        self.lock = threading.Lock()
        self.start_time = time.monotonic()

    def update(self, success=True):
        with self.lock:
            self.completed += 1
            if not success:
                self.failed += 1
            self._print()

    def _print(self):
        elapsed = time.monotonic() - self.start_time
        pct = self.completed / self.total * 100
        if self.completed > 0:
            eta = elapsed / self.completed * (self.total - self.completed)
            m, s = divmod(int(eta), 60)
            eta_str = f"{m}:{s:02d}"
        else:
            eta_str = "?"
        fail_str = f" | failed {self.failed}" if self.failed else ""
        print(
            f"\rRunning: {self.completed}/{self.total} ({pct:.1f}%)"
            f" | ETA {eta_str}{fail_str}   ",
            end="",
            file=sys.stderr,
            flush=True,
        )

    def finish(self):
        elapsed = time.monotonic() - self.start_time
        m, s = divmod(int(elapsed), 60)
        print(
            f"\rCompleted {self.total} jobs in {m}:{s:02d}"
            f" ({self.total - self.failed} ok, {self.failed} failed).   ",
            file=sys.stderr,
        )


def run_compare_job(job, progress=None):
    """
    Run compare_genomes2.py for one simulation directory.

    job is a tuple:
      (d, insertion_rate, deletion_rate, solo_ratio, length_bias,
       exp_fasta, exp_tsv,
       compare_script, ref_tsv, ref_fasta, dist_metric, alphas)
    """
    (d,
     insertion_rate, deletion_rate, solo_ratio, length_bias,
     exp_fasta, exp_tsv,
     compare_script, ref_tsv, ref_fasta, dist_metric, alphas) = job

    cmd = [
        sys.executable,  # same python interpreter
        compare_script,
        "--ref-tsv", ref_tsv,
        "--exp-tsv", exp_tsv,
        "--ref-fasta", ref_fasta,
        "--exp-fasta", exp_fasta,
        "--dist-metric", dist_metric,
        "--alphas",
        *(str(a) for a in alphas),
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(
            f"\n{os.path.basename(compare_script)} failed for {d} "
            f"(exit {result.returncode}); stderr:\n{result.stderr}",
            file=sys.stderr,
        )
        if progress:
            progress.update(success=False)
        return None

    metrics = parse_compare_output(result.stdout)
    if metrics is None:
        print(
            f"\nCould not parse metrics from {os.path.basename(compare_script)} output for {d}",
            file=sys.stderr,
        )
        if progress:
            progress.update(success=False)
        return None

    if progress:
        progress.update(success=True)

    return [
        insertion_rate,
        deletion_rate,
        solo_ratio,
        length_bias,
        # Raw values (ref / exp)
        int(metrics["ref_sample_count"]),
        int(metrics["exp_sample_count"]),
        int(metrics["ref_cumulative_length"]),
        int(metrics["exp_cumulative_length"]),
        int(metrics["ref_genome_size"]),
        int(metrics["exp_genome_size"]),
        metrics["raw_distribution"],
        # Per-metric distances
        metrics["d_sample_count"],
        metrics["d_cumulative_length"],
        metrics["d_genome_size"],
        metrics["d_distribution"],
        # Weights
        metrics["alpha_sample_count"],
        metrics["alpha_cumulative_length"],
        metrics["alpha_genome_size"],
        metrics["alpha_distribution"],
        # Composite
        metrics["Composite"],
    ]


# Column indices for the numeric fields that need imputation
# (offsets within a row, starting after the 4 param columns)
_METRIC_OFFSET = 4   # first metric column index
_NUM_METRICS = 17     # columns from ref_ltr_rt_count .. Composite


def main():
    parser = argparse.ArgumentParser(
        description="Run compare_genomes2.py over all simulation dirs and build a composite matrix."
    )
    parser.add_argument(
        "--compare-script",
        default="compare_genomes2.py",
        help="Path to compare_genomes2.py (default: compare_genomes2.py in $PWD)",
    )
    parser.add_argument(
        "--ref-tsv",
        required=True,
        help="Reference TSV path (e.g. ../Athal.fa_Kmer2LTR_TSD_class_purge)",
    )
    parser.add_argument(
        "--ref-fasta",
        required=True,
        help="Reference FASTA path (e.g. ../Athal.fa)",
    )
    parser.add_argument(
        "--output",
        default="composite_matrix.tsv",
        help="Output TSV file (default: composite_matrix.tsv)",
    )
    parser.add_argument(
        "--pattern",
        default="insertion_rates_*_deletion_rates_*_solo_ratio_*_length_bias_*",
        help="Glob pattern for simulation dirs (default: insertion_rates_*_deletion_rates_*_solo_ratio_*_length_bias_*)",
    )
    parser.add_argument(
        "--gen-prefix",
        default="gen5100000_final.fasta",
        help="Base name of generated FASTA/TSV prefix (default: gen5100000_final.fasta)",
    )
    parser.add_argument(
        "--dist-metric",
        choices=["wasserstein", "ks", "rms"],
        default="wasserstein",
        help="Distribution distance metric to use in compare_genomes2.py "
             "(default: wasserstein).",
    )
    parser.add_argument(
        "--alphas",
        nargs=4,
        type=float,
        metavar=("ALPHA_SAMPLE", "ALPHA_LENGTH", "ALPHA_GENOME", "ALPHA_DISTRIBUTION"),
        default=[0.0, 0.0, 10.0, 1.0],
        help="Four alpha weights passed to compare_genomes2.py via --alphas "
             "(default: 0 0 10 1).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of simulation dirs to process in parallel (default: 1)",
    )

    args = parser.parse_args()

    # Find all candidate directories
    candidate_dirs = sorted(
        d for d in glob.glob(args.pattern) if os.path.isdir(d)
    )

    if not candidate_dirs:
        print("No matching directories found.", file=sys.stderr)
        sys.exit(1)

    rows = []  # successfully parsed rows
    skipped_param_rows = []  # dirs missing prereq files but with valid params
    jobs = []  # jobs for parallel compare runs
    skipped_count = 0

    # First pass: gather parameters, check required files, build jobs list
    for d in candidate_dirs:
        base = os.path.basename(d)
        m = DIR_REGEX.match(base)
        if not m:
            skipped_count += 1
            continue

        insertion_rate = m.group("insertion_rate")
        deletion_rate = m.group("deletion_rate")
        solo_ratio = m.group("solo_ratio")
        length_bias = m.group("length_bias")

        exp_fasta = os.path.join(d, args.gen_prefix)
        base = os.path.basename(exp_fasta)

        # Remove extension (.fasta/.fa) if present
        base_noext = re.sub(r"\.(?:fasta|fa)$", "", base)

        # Remove trailing "_final" if present
        base_noext = re.sub(r"_final$", "", base_noext)

        exp_tsv = os.path.join(d, base_noext + "_ltr_r1_kmer2ltr_dedup")

        missing = False
        if not os.path.exists(exp_fasta):
            missing = True
        if not os.path.exists(exp_tsv):
            missing = True

        if missing:
            skipped_count += 1
            skipped_param_rows.append([
                insertion_rate,
                deletion_rate,
                solo_ratio,
                length_bias,
            ])
            continue

        # Build a job tuple for this directory
        jobs.append((
            d,
            insertion_rate,
            deletion_rate,
            solo_ratio,
            length_bias,
            exp_fasta,
            exp_tsv,
            args.compare_script,
            args.ref_tsv,
            args.ref_fasta,
            args.dist_metric,
            args.alphas,
        ))

    print(
        f"Found {len(candidate_dirs)} dirs: "
        f"{len(jobs)} to process, {skipped_count} skipped (missing files/bad name).",
        file=sys.stderr,
    )

    # Second pass: run compare_genomes2.py over all jobs, possibly in parallel
    if jobs:
        max_workers = max(1, args.threads if args.threads and args.threads >= 1 else 1)
        progress = ProgressTracker(len(jobs))

        if max_workers == 1:
            for job in jobs:
                row = run_compare_job(job, progress=progress)
                if row is not None:
                    rows.append(row)
        else:
            print(f"Using {max_workers} threads.", file=sys.stderr)
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_job = {
                    executor.submit(run_compare_job, job, progress): job
                    for job in jobs
                }
                for future in as_completed(future_to_job):
                    row = future.result()
                    if row is not None:
                        rows.append(row)

        progress.finish()

    # If we have at least one valid row, compute max metrics and
    # impute values for skipped_param_rows as max + 0.1
    if rows and skipped_param_rows:
        # For imputation, only use the distance / composite columns
        # (indices relative to row list)
        #   idx  4-5:  ref/exp LTR-RT count       (int, skip for imputation)
        #   idx  6-7:  ref/exp cumulative length   (int, skip for imputation)
        #   idx  8-9:  ref/exp genome size         (int, skip for imputation)
        #   idx 10:    raw_distribution            (skip)
        #   idx 11-14: d_sample_count .. d_distribution
        #   idx 15-18: alpha weights               (skip)
        #   idx 19:    Composite

        # For imputed rows: raw counts get 0, distances get max+0.1,
        # weights mirror the alphas used.
        offset = 0.1

        imputed_d_sample      = max(r[11] for r in rows) + offset
        imputed_d_cumlength   = max(r[12] for r in rows) + offset
        imputed_d_genome      = max(r[13] for r in rows) + offset
        imputed_d_distribution = max(r[14] for r in rows) + offset
        imputed_composite     = max(r[19] for r in rows) + offset

        for params in skipped_param_rows:
            rows.append(params + [
                0, 0,           # ref/exp LTR-RT count (unknown)
                0, 0,           # ref/exp cumulative length (unknown)
                0, 0,           # ref/exp genome size (unknown)
                0.0,            # raw distribution distance (unknown)
                imputed_d_sample,
                imputed_d_cumlength,
                imputed_d_genome,
                imputed_d_distribution,
                args.alphas[0], # alpha_sample_count
                args.alphas[1], # alpha_cumulative_length
                args.alphas[2], # alpha_genome_size
                args.alphas[3], # alpha_distribution
                imputed_composite,
            ])

    elif not rows and skipped_param_rows:
        print(
            "Warning: no successful compare_genomes2.py runs; "
            "imputing all metrics as 0.1 for skipped dirs.",
            file=sys.stderr,
        )
        for params in skipped_param_rows:
            rows.append(params + [
                0, 0, 0, 0, 0, 0, 0.0,
                0.1, 0.1, 0.1, 0.1,
                args.alphas[0], args.alphas[1], args.alphas[2], args.alphas[3],
                0.1,
            ])

    # Write output matrix
    header = [
        "insertion_rate",
        "deletion_rate",
        "solo_ratio",
        "length_bias",
        # Raw counts (ref / exp)
        "ref_ltr_rt_count",
        "exp_ltr_rt_count",
        "ref_cumulative_length",
        "exp_cumulative_length",
        "ref_genome_size",
        "exp_genome_size",
        # Raw distribution distance
        "raw_distribution_distance",
        # Per-metric distances
        "d_sample_count",
        "d_cumulative_length",
        "d_genome_size",
        "d_distribution",
        # Weights
        "alpha_sample_count",
        "alpha_cumulative_length",
        "alpha_genome_size",
        "alpha_distribution",
        # Composite
        "Composite",
    ]

    with open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)

    print(f"Wrote {len(rows)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
