#!/usr/bin/env python3
import sys, os, csv, re

# Rate reflect itteration rate. I need to add a parameter to accept step size so it can output generation rate instead. 

def format_scientific(x):
    """Format x in minimal scientific notation, e.g. 4.5e-2, 1.3e-42."""
    s = "{:.1e}".format(x)           # one significant decimal
    s = s.replace("E", "e").replace("+", "")
    # strip leading zeros in exponent: e-02 -> e-2, e+02 -> e2
    s = re.sub(r"e([+-])0*([0-9]+)", r"e\1\2", s)
    return s

def count_bp(fasta_path):
    """Sum lengths of all sequence lines in a FASTA."""
    total = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total

def main(report_path):
    # read input
    with open(report_path) as inf:
        reader = csv.DictReader(inf, delimiter="\t")
        # extend header
        out_fields = reader.fieldnames + ["genome_size", "insertion_rate", "deletion_rate"]
        writer = csv.DictWriter(sys.stdout, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for row in reader:
            gen = row["Generation"]
            fasta = f"gen{gen}_final.fasta"
            if not os.path.isfile(fasta):
                sys.stderr.write(f"Warning: {fasta} not found; skipping rates for generation {gen}\n")
                size = ""
                irate = ""
                drate = ""
            else:
                size = count_bp(fasta)
                # TE_inserts is the number before the "("
                te_insert = int(row["TE_inserts(nest/nonnest)"].split("(")[0])
                del_actual = int(row["Actual_TE_deletions"])
                irate = format_scientific(te_insert / size)
                drate = format_scientific(del_actual / size)

            row["genome_size"]      = size
            row["insertion_rate"]   = irate
            row["deletion_rate"]    = drate
            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python calculate_rates.py pipeline.report")
    main(sys.argv[1])
