#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BURNIN    } from './nf/modules/burnin'
include { EMIT_GRID } from './nf/modules/emit_grid'
include { SIMULATE  } from './nf/modules/simulate'
include { SCORE     } from './nf/modules/score'
include { AGGREGATE } from './nf/modules/aggregate'
include { CONTOUR   } from './nf/modules/contour'

def helpMessage() {
    log.info """
    PrinTE - forward simulation of transposable-element genome evolution

    Usage:
      nextflow run . -profile <conda|docker|singularity>,<local|slurm|awsbatch> [options]

    Modes:
      --mode simulate   one parameter set, burn-in through evolve   (default)
      --mode sweep      scatter over a parameter grid, then score every run
                        against a reference genome and aggregate

    Starting genome (pick one):
      --fasta, --bed        your own genome and PrinTE-format BED
      --size, --chr_number  build a burn-in instead

    TE library (pick one):
      --te_lib     RepeatMasker-headed FASTA, preprocessed by PrinTE
      --clean_lib  already preprocessed, skips that step

    Sweep also needs:
      --ref_tsv, --ref_fasta   the real genome to score against

    Try it:
      nextflow run . -profile test,conda
    """.stripIndent()
}

workflow SIMULATE_ONE {
    take:
    start
    lib
    ratios

    main:
    combo = Channel.of(
        [params.insert_rate, params.delete_rate, params.solo_rate, params.length_bias, 'single']
    )
    SIMULATE(combo.combine(start), lib, ratios)

    emit:
    SIMULATE.out.final_gen
}

workflow SWEEP {
    take:
    start
    lib
    ratios

    main:
    combos = EMIT_GRID()
        .splitCsv(header: true)
        .map { r -> [r.ins, r.'del', r.sr, r.k, r.dirname] }

    SIMULATE(combos.combine(start), lib, ratios)

    // Scoring needs each simulation's LTR-RT table, and dating those is an external
    // annotation step (see docs/ltr-annotation.md). The sweep therefore scores only the
    // simulations whose table is already present, and skips the rest. Annotate, then
    // re-run the same command with -resume to pick up scoring without re-simulating.
    scored = SIMULATE.out.final_gen
        .map { d, i, dl, s, k, fa, bed, lb ->
            [d, i, dl, s, k, fa, bed, lb,
             file("${params.outdir}/simulations/${d}/${params.exp_tsv_name}")]
        }
        .filter { it[8].exists() }

    SCORE(scored, file(params.ref_tsv), file(params.ref_fasta))
    AGGREGATE(SCORE.out.collect())
    CONTOUR(AGGREGATE.out)

    emit:
    AGGREGATE.out
}

workflow {
    if (params.help) {
        helpMessage()
        return
    }

    if (!params.clean_lib && !params.te_lib) {
        error("Provide --clean_lib or --te_lib")
    }
    lib = file(params.clean_lib ?: params.te_lib)

    // Placeholder files stand in for optional inputs so the process signatures stay
    // fixed; each process checks the name and drops the flag.
    cds    = params.cds    ? file(params.cds)    : file("${projectDir}/assets/NO_CDS")
    ratios = params.ratios ? file(params.ratios) : file("${projectDir}/assets/NO_RATIOS")

    if (params.fasta && params.bed) {
        start = Channel.of([file(params.fasta), file(params.bed)])
    } else {
        BURNIN(lib, cds, ratios)
        start = BURNIN.out.start
    }

    if (params.mode == 'sweep') {
        if (!params.ref_tsv || !params.ref_fasta) {
            error("--mode sweep needs --ref_tsv and --ref_fasta to score against. See docs/nextflow.md for the annotation step that comes first.")
        }
        SWEEP(start, lib, ratios)
    } else {
        SIMULATE_ONE(start, lib, ratios)
    }
}
