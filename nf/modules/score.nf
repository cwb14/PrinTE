process SCORE {
    label 'process_medium'
    tag "${dirname}"
    publishDir "${params.outdir}/scores", mode: 'copy'

    input:
    tuple val(dirname), val(ins), val(del), val(sr), val(k),
          path(fasta), path(bed), path(lib), path(exp_tsv)
    path ref_tsv
    path ref_fasta

    output:
    path "${dirname}.score.tsv"

    script:
    """
    python -m printe.grid.compare_genomes \\
        --ref-tsv ${ref_tsv} --ref-fasta ${ref_fasta} \\
        --exp-tsv ${exp_tsv} --exp-fasta ${fasta} \\
        --dist-metric ${params.dist_metric} \\
        --alphas ${params.alphas} > ${dirname}.score.tsv
    """
}
