#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    HPRC Ensembl vs CAT Annotation Comparison Pipeline
========================================================================================
    Compares Ensembl and CAT gene annotations for HPRC assemblies

    Usage:
        nextflow run main.nf --input assemblies.csv --outdir results

    Required Parameters:
        --input             CSV file with columns: assembly_accession,sample_name
        --outdir            Output directory for results

    Optional Parameters:
        --ensembl_cache_dir         Directory to cache Ensembl GFFs (default: ./cache/ensembl)
        --cat_cache_dir             Directory to cache CAT GFFs (default: ./cache/cat)
        --assembly_reports_dir      Directory with assembly reports (optional)
        --comparison_script         Path to comparison Python script
        --max_assemblies            Limit number of assemblies to process (for testing)
========================================================================================
*/

// Default parameters - MUST be defined before includes to avoid warnings
params.input = null
params.outdir = './results'
params.ensembl_cache_dir = './cache/ensembl'
params.cat_cache_dir = './cache/cat'
params.assembly_reports_dir = './cache/assembly_reports'
params.comparison_script = "${projectDir}/../../../packages/hprc-qc-annotation/src/hprc_qc_annotation/hprc_ensembl_cat_overlap.py"
params.max_assemblies = null
params.help = false

include { ENSEMBL_CAT_COMPARISON } from './workflows/ensembl_cat_comparison'

// Help message
def helpMessage() {
    log.info"""
    ============================================
    HPRC Ensembl vs CAT Comparison Pipeline
    ============================================

    Usage:
        nextflow run main.nf --input assemblies.csv --outdir results

    Required arguments:
        --input               CSV file with 'assembly_accession' and 'sample_name' columns
        --outdir              Output directory (default: ${params.outdir})

    Optional arguments:
        --ensembl_cache_dir   Cache directory for Ensembl GFFs (default: ${params.ensembl_cache_dir})
        --cat_cache_dir       Cache directory for CAT GFFs (default: ${params.cat_cache_dir})
        --assembly_reports_dir  Directory containing assembly reports (optional)
        --comparison_script   Path to comparison script (default: auto-detected)
        --max_assemblies      Limit number of assemblies (for testing, e.g., --max_assemblies 3)

    Example CSV format:
        assembly_accession,sample_name
        GCA_018466835.2,HG00408_pat
        GCA_041900255.1,HG00408_mat

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (!params.input) {
    log.error "ERROR: --input parameter is required"
    helpMessage()
    exit 1
}

if (!file(params.comparison_script).exists()) {
    log.error "ERROR: Comparison script not found at: ${params.comparison_script}"
    log.error "Please specify correct path with --comparison_script"
    exit 1
}

// Create cache directories if they don't exist
file(params.ensembl_cache_dir).mkdirs()
file(params.cat_cache_dir).mkdirs()

workflow {
    // Read input CSV
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.assembly_accession || !row.sample_name) {
                error "CSV must contain 'assembly_accession' and 'sample_name' columns"
            }
            tuple(row.assembly_accession, row.sample_name)
        }
        .take(params.max_assemblies ?: -1)  // Limit if specified
        .set { assemblies_ch }

    // Count assemblies
    assemblies_ch.count().view { count ->
        """
        ============================================
        Processing ${count} assemblies
        ============================================
        """
    }

    // Run comparison workflow
    ENSEMBL_CAT_COMPARISON(assemblies_ch)

    // Collect results
    ENSEMBL_CAT_COMPARISON.out.rbh.collect().view { results ->
        """
        ============================================
        Completed! RBH files created: ${results.size()}
        Results in: ${params.outdir}/results/
        ============================================
        """
    }
}

workflow.onComplete {
    log.info ""
    log.info "Pipeline execution summary"
    log.info "=========================="
    log.info "Completed at : ${workflow.complete}"
    log.info "Duration     : ${workflow.duration}"
    log.info "Success      : ${workflow.success}"
    log.info "Results dir  : ${params.outdir}"
    log.info "Work dir     : ${workflow.workDir}"
    log.info ""

    if (!workflow.success) {
        log.error "Pipeline failed. Check logs in ${params.outdir}/logs/"
    }
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorReport}"
}
