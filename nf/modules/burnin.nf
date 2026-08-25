process BURNIN {
    label 'process_medium'
    publishDir "${params.outdir}/burnin", mode: 'copy'

    input:
    path lib
    path cds
    path ratios

    output:
    tuple path('burnin.fasta'), path('burnin.bed'), emit: start
    path 'burnin.stat'
    path 'pipeline.log'

    script:
    def libflag = params.clean_lib ? "-cl ${lib}" : "-i ${lib}"
    def cdsflag = cds.name != 'NO_CDS' ? "-c ${cds}" : ''
    def ratflag = ratios.name != 'NO_RATIOS' ? "-r ${ratios}" : ''
    """
    printe --burnin_only \\
        -sz ${params.size} -cn ${params.chr_number} \\
        -P ${params.cds_percent} \\
        -itp ${params.intact_te_percent} -ftp ${params.frag_te_percent} \\
        ${libflag} ${cdsflag} ${ratflag} \\
        -m ${params.mutation_rate} -TsTv ${params.tstv} \\
        -t ${task.cpus} -s ${params.seed}
    """
}
