process AGGREGATE {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scores

    output:
    path 'composite_matrix.tsv'

    script:
    """
    python - <<'PY'
    import glob, sys
    rows, header = [], None
    for path in sorted(glob.glob("*.score.tsv")):
        lines = [l.rstrip("\\n") for l in open(path) if l.strip()]
        if not lines:
            continue
        if header is None:
            header = lines[0]
        rows.extend(lines[1:] if lines[0] == header else lines)
    if header is None:
        sys.exit("no score files to aggregate")
    with open("composite_matrix.tsv", "w") as out:
        out.write(header + "\\n")
        out.write("\\n".join(rows) + "\\n")
    print(f"aggregated {len(rows)} rows")
    PY
    """
}
