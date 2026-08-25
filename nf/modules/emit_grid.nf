process EMIT_GRID {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    output:
    path 'combos.csv'

    script:
    """
    python -m printe.grid.emit_combos \\
        --ins-start ${params.ins_start} --ins-end ${params.ins_end} --ins-count ${params.ins_count} \\
        --del-start ${params.del_start} --del-end ${params.del_end} --del-count ${params.del_count} \\
        --sr-start ${params.sr_start} --sr-end ${params.sr_end} --sr-step ${params.sr_step} \\
        --k-start ${params.k_start} --k-end ${params.k_end} --k-step ${params.k_step} \\
        --out combos.csv
    """
}
