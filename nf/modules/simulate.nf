process SIMULATE {
    label 'process_high'
    tag "${dirname}"
    publishDir path: { "${params.outdir}/simulations/${dirname}" }, mode: 'copy'

    input:
    tuple val(ins), val(del), val(sr), val(k), val(dirname), path(fasta), path(bed)
    path lib
    path ratios

    output:
    tuple val(dirname), val(ins), val(del), val(sr), val(k),
          path("gen${params.generation_end}_final.fasta"),
          path("gen${params.generation_end}_final.bed"),
          path("gen${params.generation_end}_final.lib"), emit: final_gen
    path 'pipeline.log'
    path 'pipeline.report', optional: true

    script:
    def libflag = params.clean_lib ? "-cl ${lib}" : "-i ${lib}"
    def ratflag = ratios.name != 'NO_RATIOS' ? "-r ${ratios}" : ''
    def post    = params.postproc ? '' : '--no_postproc'
    """
    printe \\
        -f ${fasta} -b ${bed} \\
        ${libflag} ${ratflag} \\
        -ir ${ins} -dr ${del} -sr ${sr} -k ${k} \\
        -ge ${params.generation_end} -st ${params.step} \\
        -m ${params.mutation_rate} -TsTv ${params.tstv} \\
        -t ${task.cpus} -s ${params.seed} ${post}
    """
}
