process CONTOUR {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path matrix

    output:
    path '*.pdf', optional: true
    path '*.txt', optional: true

    script:
    """
    python -m printe.grid.contour_plot --input ${matrix}
    """
}
