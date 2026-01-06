#!/usr/bin/env python3
import os
import re
import glob
import csv
import sys
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

DIR_REGEX = re.compile(
    r"^insertion_rates_(?P<insertion_rate>[^_]+)"
    r"_deletion_rates_(?P<deletion_rate>[^_]+)"
    r"_solo_ratio_(?P<solo_ratio>[^_]+)"
    r"_length_bias_(?P<length_bias>[^_]+)$"
)

# Allow floats and scientific notation
FLOAT_RE = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"

# Updated to use the new density metric line: d_distribution (...)
METRIC_PATTERNS = {
    "d_sample_count":      re.compile(r"d_sample_count.*?"      + FLOAT_RE),
    "d_cumulative_length": re.compile(r"d_cumulative_length.*?" + FLOAT_RE),
    "d_genome_size":       re.compile(r"d_genome_size.*?"       + FLOAT_RE),
    # This line looks like:
    #   d_distribution (Kolmogorov–Smirnov): 0.134263
    #   d_distribution (Wasserstein): 0.00976438
    #   d_distribution (RMS (binned density)): 2.66063
    "d_distribution":      re.compile(r"d_distribution.*?"      + FLOAT_RE),
    # Assuming Composite is still printed in compare_genomes2.py
    "Composite":           re.compile(r"Composite\s*=\s*Σ[^=]*=\s*" + FLOAT_RE),
}


def parse_compare_output(stdout_text):
    """Extract numeric metrics from compare_genomes2.py output."""
    metrics = {}
    for key, pattern in METRIC_PATTERNS.items():
        m = pattern.search(stdout_text)
        if not m:
            return None  # missing something, treat as failure
        try:
            metrics[key] = float(m.group(1))
        except ValueError:
            return None
    return metrics


def run_compare_job(job):
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

    # Note: stderr prints may interleave when using multiple threads, which is fine for a log.
    print(f"Running ({os.path.basename(d)}): {' '.join(cmd)}", file=sys.stderr)

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(
            f"{os.path.basename(compare_script)} failed for {d} "
            f"(exit {result.returncode}); stderr:\n{result.stderr}",
            file=sys.stderr,
        )
        return None

    metrics = parse_compare_output(result.stdout)
    if metrics is None:
        print(
            f"Could not parse metrics from {os.path.basename(compare_script)} output for {d}",
            file=sys.stderr,
        )
        return None

    return [
        insertion_rate,
        deletion_rate,
        solo_ratio,
        length_bias,
        metrics["d_sample_count"],
        metrics["d_cumulative_length"],
        metrics["d_genome_size"],
        metrics["d_distribution"],
        metrics["Composite"],
    ]


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

    # First pass: gather parameters, check required files, build jobs list
    for d in candidate_dirs:
        base = os.path.basename(d)
        m = DIR_REGEX.match(base)
        if not m:
            # Not exactly in the expected format; skip quietly
            print(f"Skipping {d}: name doesn't match expected pattern", file=sys.stderr)
            continue

        insertion_rate = m.group("insertion_rate")
        deletion_rate = m.group("deletion_rate")
        solo_ratio = m.group("solo_ratio")
        length_bias = m.group("length_bias")

        exp_fasta = os.path.join(d, args.gen_prefix)
        # Derive exp_tsv from exp_fasta, e.g.
        #   gen5100000_final.fasta -> gen5100000_ltrharvest_kmer2ltr
        #   gen5100000.fasta       -> gen5100000_ltrharvest_kmer2ltr
        base = os.path.basename(exp_fasta)

        # Remove extension (.fasta/.fa) if present
        base_noext = re.sub(r"\.(?:fasta|fa)$", "", base)

        # Remove trailing "_final" if present
        base_noext = re.sub(r"_final$", "", base_noext)

        exp_tsv = os.path.join(d, base_noext + "_ltrharvest_kmer2ltr_dedup")

        missing = False
        if not os.path.exists(exp_fasta):
            print(f"Skipping {d}: missing {exp_fasta}", file=sys.stderr)
            missing = True
        if not os.path.exists(exp_tsv):
            print(f"Skipping {d}: missing {exp_tsv}", file=sys.stderr)
            missing = True

        if missing:
            # store just the parameter part; metrics will be filled later
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

    # Second pass: run compare_genomes2.py over all jobs, possibly in parallel
    if jobs:
        if args.threads is None or args.threads < 1:
            max_workers = 1
        else:
            max_workers = args.threads

        if max_workers == 1:
            # Serial execution (original behavior)
            for job in jobs:
                row = run_compare_job(job)
                if row is not None:
                    rows.append(row)
        else:
            print(f"Processing {len(jobs)} jobs with {max_workers} threads.", file=sys.stderr)
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_job = {executor.submit(run_compare_job, job): job for job in jobs}
                for future in as_completed(future_to_job):
                    row = future.result()
                    if row is not None:
                        rows.append(row)

    # If we have at least one valid row, compute max metrics and
    # impute values for skipped_param_rows as max + 0.1
    if rows and skipped_param_rows:
        # metric indices in rows: 4..8
        max_d_sample_count = max(r[4] for r in rows)
        max_d_cumulative_length = max(r[5] for r in rows)
        max_d_genome_size = max(r[6] for r in rows)
        max_d_distribution = max(r[7] for r in rows)
        max_composite = max(r[8] for r in rows)

        offset = 0.1

        imputed_sample = max_d_sample_count + offset
        imputed_cumulative_length = max_d_cumulative_length + offset
        imputed_genome_size = max_d_genome_size + offset
        imputed_distribution = max_d_distribution + offset
        imputed_composite = max_composite + offset

        for params in skipped_param_rows:
            rows.append(params + [
                imputed_sample,
                imputed_cumulative_length,
                imputed_genome_size,
                imputed_distribution,
                imputed_composite,
            ])

    elif not rows and skipped_param_rows:
        # Edge case: everything was skipped due to missing files.
        # We can still emit rows with a small constant (0.1) so the
        # matrix isn't empty. Adjust if you prefer a different behavior.
        print(
            "Warning: no successful compare_genomes2.py runs; "
            "imputing all metrics as 0.1 for skipped dirs.",
            file=sys.stderr,
        )
        for params in skipped_param_rows:
            rows.append(params + [0.1, 0.1, 0.1, 0.1, 0.1])

    # Write output matrix
    header = [
        "insertion_rate",
        "deletion_rate",
        "solo_ratio",
        "length_bias",
        "d_sample_count",
        "d_cumulative_length",
        "d_genome_size",
        "d_distribution",   # generic for ks/wasserstein/rms
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
