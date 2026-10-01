#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PANACUS_STATS } from '../subworkflows/panacus/main.nf'

workflow ALIGNMENT_QC {
    take:
    ch_databases // channel: [ core_dbname, assembly_accession ]
    config // map: database connection configuration

    main:
    ch_reports = Channel.empty()

    if (params.run_panacus_stats) {
        PANACUS_STATS(
            ch_databases,
            config,
        )
        ch_reports = ch_reports.mix(PANACUS_STATS.out.stats.map { it -> ['mapping_stats', it] })
    }

    emit:
    reports = ch_reports
}
