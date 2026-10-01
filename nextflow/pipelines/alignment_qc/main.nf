#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

include { ALIGNMENT_QC } from './workflows/alignment_qc.nf'

workflow HPRC_ALIGNMENT_QC {
    
    if (!params.input && !params.graph) {
        error "Either --input text file or --graph must be provided"
    }
    
    if (params.input) {
        ch_input = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitText()
            .map { row ->
                [row]
            }
    } else {
        ch_input = Channel.of([params.graph])
    }
    
    config = [
        panacus_template:         params.panacus_template,
    ]
    
    ALIGNMENT_QC (
        ch_input,
        config
    )
}

workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'failed'}"
    log.info "Results published to: ${params.outdir}"
}

workflow {
    HPRC_ALIGNMENT_QC()
}
