#!/usr/bin/env python3
import sys, os, csv, re

def format_scientific(x):
    """Format x in minimal scientific notation, e.g. 4.5e-2, 1.3e-42."""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return ""
    if x == 0:
        return "0"
    s = "{:.1e}".format(x)  # one decimal in mantissa
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

def parse_generation_value(row, gen_field="Generation"):
    """
    Parse the generation value from the 'Generation' column.
    Expected to be an integer like 10000, 20000, ...
    """
    if gen_field not in row:
        return None
    s = str(row[gen_field]).strip()
    if s == "":
        return None
    try:
        return int(s)
    except ValueError:
        return None

def main(report_path):
    with open(report_path) as inf:
        reader = csv.DictReader(inf, delimiter="\t")
        if not reader.fieldnames:
            sys.exit("Error: input report appears to have no header.")

        out_fields = reader.fieldnames + [
            "genome_size",
            "insertion_rate", "deletion_rate",                 # per-iteration (counts / bp)
            "generations_in_iteration",                        # derived from Generation deltas
            "TE_inserts_per_generation", "TE_deletions_per_generation",  # per-generation counts
            "insertion_rate_per_generation", "deletion_rate_per_generation"  # per-generation (counts/bp/gen)
        ]

        writer = csv.DictWriter(sys.stdout, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        prev_gen = 0  # first line subtracts 0, per your spec

        for row in reader:
            gen = parse_generation_value(row, "Generation")
            if gen is None:
                sys.stderr.write("Warning: could not parse Generation for a row; leaving per-generation fields blank.\n")
                gen_span = ""
            else:
                gen_span = gen - prev_gen
                if gen_span <= 0:
                    sys.stderr.write(
                        f"Warning: non-increasing Generation encountered (prev={prev_gen}, current={gen}); "
                        f"setting generations_in_iteration to blank for this row.\n"
                    )
                    gen_span = ""
                prev_gen = gen

            row["generations_in_iteration"] = gen_span

            # FASTA assumed to be tied to the Generation value
            fasta = f"gen{row.get('Generation','')}_final.fasta"

            if not os.path.isfile(fasta):
                sys.stderr.write(f"Warning: {fasta} not found; skipping genome_size/rates for generation {row.get('Generation','')}\n")
                row["genome_size"] = ""
                row["insertion_rate"] = ""
                row["deletion_rate"] = ""
                row["TE_inserts_per_generation"] = ""
                row["TE_deletions_per_generation"] = ""
                row["insertion_rate_per_generation"] = ""
                row["deletion_rate_per_generation"] = ""
                writer.writerow(row)
                continue

            size = count_bp(fasta)
            row["genome_size"] = size

            # TE_inserts is the number before the "("
            te_insert = ""
            try:
                te_insert = int(str(row["TE_inserts(nest/nonnest)"]).split("(")[0].strip())
            except Exception:
                te_insert = ""

            del_actual = ""
            try:
                del_actual = int(str(row["Actual_TE_deletions"]).strip())
            except Exception:
                del_actual = ""

            # Per-iteration rates
            row["insertion_rate"] = format_scientific(te_insert / size) if te_insert != "" else ""
            row["deletion_rate"]  = format_scientific(del_actual / size) if del_actual != "" else ""

            # Per-generation counts and rates (divide by gen_span)
            if gen_span == "" or gen_span == 0:
                row["TE_inserts_per_generation"] = ""
                row["TE_deletions_per_generation"] = ""
                row["insertion_rate_per_generation"] = ""
                row["deletion_rate_per_generation"] = ""
            else:
                if te_insert != "":
                    ins_pg = te_insert / gen_span
                    row["TE_inserts_per_generation"] = ins_pg
                    row["insertion_rate_per_generation"] = format_scientific(ins_pg / size)
                else:
                    row["TE_inserts_per_generation"] = ""
                    row["insertion_rate_per_generation"] = ""

                if del_actual != "":
                    del_pg = del_actual / gen_span
                    row["TE_deletions_per_generation"] = del_pg
                    row["deletion_rate_per_generation"] = format_scientific(del_pg / size)
                else:
                    row["TE_deletions_per_generation"] = ""
                    row["deletion_rate_per_generation"] = ""

            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python calculate_rates.py pipeline.report")
    main(sys.argv[1])
