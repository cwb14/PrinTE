#!/usr/bin/env python3
"""
Parse a TE log file and generate a TSV report with columns:
Generation, TE_inserts(nest/nonnest), Calculated_TE_deletions, Actual_TE_deletions
Usage:
    python log_to_report.py -in logfile.txt -out report.tsv
"""
import re
import argparse

def parse_log(infile):
    # Patterns
    p_total = re.compile(r"Total TE insertions performed:\s*(\d+)\s*\(Nested:\s*(\d+),\s*Non-nested:\s*(\d+)\)")
    p_calc = re.compile(r"Calculated number of TE excisions:\s*(\d+)")
    p_sel = re.compile(r"Selected\s*(\d+)\s*removal events")
    p_gen = re.compile(r"Updated FASTA written to gen(\d+)_final\.fasta")

    records = []
    current = None

    for line in infile:
        # Total inserts
        m = p_total.search(line)
        if m:
            # Start new record
            current = {
                'inserts': int(m.group(1)),
                'nested': int(m.group(2)),
                'nonnest': int(m.group(3)),
                'deletions_calc': None,
                'deletions_actual': 0,
                'generation': None,
            }
            continue

        if current is None:
            continue

        # Calculated deletions
        m = p_calc.search(line)
        if m:
            current['deletions_calc'] = int(m.group(1))
            continue

        # Selected removal events
        m = p_sel.search(line)
        if m:
            current['deletions_actual'] += int(m.group(1))
            continue

        # Generation identifier
        m = p_gen.search(line)
        if m:
            current['generation'] = int(m.group(1))
            # Only add if we have all fields
            if (current['deletions_calc'] is not None
                and current['generation'] is not None):
                records.append(current)
            current = None
            continue

    return records


def write_report(records, outfile):
    # Write TSV header
    outfile.write("Generation\tTE_inserts(nest/nonnest)\tCalculated_TE_deletions\tActual_TE_deletions\n")
    for rec in records:
        inserts_str = f"{rec['inserts']}({rec['nested']}/{rec['nonnest']})"
        outfile.write(f"{rec['generation']}\t{inserts_str}\t{rec['deletions_calc']}\t{rec['deletions_actual']}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert TE log to TSV report")
    parser.add_argument('-in', dest='infile', required=True,
                        help='Path to the log file')
    parser.add_argument('-out', dest='outfile', required=True,
                        help='Path to the output TSV file')
    args = parser.parse_args()

    try:
        with open(args.infile, 'r') as inf:
            records = parse_log(inf)
    except FileNotFoundError:
        print(f"Error: input file '{args.infile}' not found.")
        return

    if not records:
        print("Warning: no records found in log file.")

    with open(args.outfile, 'w') as outf:
        write_report(records, outf)

    print(f"Report written to {args.outfile}")

if __name__ == '__main__':
    main()
