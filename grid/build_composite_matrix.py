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

METRIC_PATTERNS = {
    "d_sample_count":      re.compile(r"d_sample_count.*?"      + FLOAT_RE),
    "d_cumulative_length": re.compile(r"d_cumulative_length.*?" + FLOAT_RE),
    "d_genome_size":       re.compile(r"d_genome_size.*?"       + FLOAT_RE),
    "d_distribution":      re.compile(r"d_distribution.*?"      + FLOAT_RE),
    "Composite":           re.compile(r"Composite\s*=\s*Σ[^=]*=\s*" + FLOAT_RE),
}


def parse_compare_output(stdout_text):
    metrics = {}
    for key, pattern in METRIC_PATTERNS.items():
        m = pattern.search(stdout_text)
        if not m:
            return None
        try:
            metrics[key] = float(m.group(1))
        except ValueError:
            return None
    return metrics


def run_compare_job(job):
    (d,
     insertion_rate, deletion_rate, solo_ratio, length_bias,
     exp_fasta, exp_tsv,
     compare_script, ref_tsv, ref_fasta, dist_metric, alphas) = job

    cmd = [
        sys.executable,
        compare_script,
        "--ref-tsv", ref_tsv,
        "--exp-tsv", exp_tsv,
        "--ref-fasta", ref_fasta,
        "--exp-fasta", exp_fasta,
        "--dist-metric", dist_metric,
        "--alphas",
        *(str(a) for a in alphas),
    ]

    print(f"Running ({os.path.basename(d)}): {' '.join(cmd)}", file=sys.stderr)

    result = subprocess.run(cmd, capture_output=True, text=True)

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
        insertion_rate, deletion_rate, solo_ratio, length_bias,
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
    parser.add_argument("--compare-script", default="compare_genomes2.py")
    parser.add_argument("--ref-tsv", required=True)
    parser.add_argument("--ref-fasta", required=True)
    parser.add_argument("--output", default="composite_matrix.tsv")
    parser.add_argument(
        "--pattern",
        default="insertion_rates_*_deletion_rates_*_solo_ratio_*_length_bias_*",
    )
    parser.add_argument("--gen-prefix", default="gen5100000_final.fasta")
    parser.add_argument(
        "--tsv-suffix",
        default="_ltrharvest_kmer2ltr_dedup",
        help="Suffix appended to the FASTA base name (minus extension and '_final') "
             "to form the experimental TSV filename. "
             "Ignored when --exp-tsv-name is provided.",
    )
    parser.add_argument(
        "--exp-tsv-name",
        default=None,
        help="Exact base name of the experimental TSV file inside each simulation dir "
             "(e.g. ltrharvest_r1_kmer2ltr_dedup). "
             "When provided, overrides --tsv-suffix and any name derived from --gen-prefix.",
    )
    parser.add_argument(
        "--dist-metric", choices=["wasserstein", "ks", "rms"], default="wasserstein",
    )
    parser.add_argument(
        "--alphas", nargs=4, type=float,
        metavar=("ALPHA_SAMPLE", "ALPHA_LENGTH", "ALPHA_GENOME", "ALPHA_DISTRIBUTION"),
        default=[0.0, 0.0, 10.0, 1.0],
    )
    parser.add_argument("--threads", type=int, default=1)

    args = parser.parse_args()

    candidate_dirs = sorted(
        d for d in glob.glob(args.pattern) if os.path.isdir(d)
    )

    if not candidate_dirs:
        print("No matching directories found.", file=sys.stderr)
        sys.exit(1)

    rows = []
    skipped_param_rows = []
    jobs = []

    for d in candidate_dirs:
        base = os.path.basename(d)
        m = DIR_REGEX.match(base)
        if not m:
            print(f"Skipping {d}: name doesn't match expected pattern", file=sys.stderr)
            continue

        insertion_rate = m.group("insertion_rate")
        deletion_rate  = m.group("deletion_rate")
        solo_ratio     = m.group("solo_ratio")
        length_bias    = m.group("length_bias")

        exp_fasta = os.path.join(d, args.gen_prefix)

        # Determine exp_tsv:
        #   --exp-tsv-name  ->  use exactly that filename (no derivation)
        #   otherwise       ->  strip ext + "_final" from gen-prefix, append --tsv-suffix
        if args.exp_tsv_name:
            exp_tsv = os.path.join(d, args.exp_tsv_name)
        else:
            base_name  = os.path.basename(exp_fasta)
            base_noext = re.sub(r"\.(?:fasta|fa)$", "", base_name)
            base_noext = re.sub(r"_final$", "", base_noext)
            exp_tsv    = os.path.join(d, base_noext + args.tsv_suffix)

        missing = False
        if not os.path.exists(exp_fasta):
            print(f"Skipping {d}: missing {exp_fasta}", file=sys.stderr)
            missing = True
        if not os.path.exists(exp_tsv):
            print(f"Skipping {d}: missing {exp_tsv}", file=sys.stderr)
            missing = True

        if missing:
            skipped_param_rows.append([insertion_rate, deletion_rate, solo_ratio, length_bias])
            continue

        jobs.append((
            d,
            insertion_rate, deletion_rate, solo_ratio, length_bias,
            exp_fasta, exp_tsv,
            args.compare_script, args.ref_tsv, args.ref_fasta,
            args.dist_metric, args.alphas,
        ))

    if jobs:
        max_workers = max(1, args.threads or 1)
        if max_workers == 1:
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

    if rows and skipped_param_rows:
        offset = 0.1
        imputed = [
            max(r[4] for r in rows) + offset,
            max(r[5] for r in rows) + offset,
            max(r[6] for r in rows) + offset,
            max(r[7] for r in rows) + offset,
            max(r[8] for r in rows) + offset,
        ]
        for params in skipped_param_rows:
            rows.append(params + imputed)
    elif not rows and skipped_param_rows:
        print(
            "Warning: no successful compare_genomes2.py runs; "
            "imputing all metrics as 0.1 for skipped dirs.",
            file=sys.stderr,
        )
        for params in skipped_param_rows:
            rows.append(params + [0.1, 0.1, 0.1, 0.1, 0.1])

    header = [
        "insertion_rate", "deletion_rate", "solo_ratio", "length_bias",
        "d_sample_count", "d_cumulative_length", "d_genome_size",
        "d_distribution", "Composite",
    ]

    with open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)

    print(f"Wrote {len(rows)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
