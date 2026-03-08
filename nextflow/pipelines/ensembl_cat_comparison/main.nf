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
        --ensg_lookup       Path to transcript ID to ENSG ID lookup table

    Optional Parameters:
        --ensembl_cache_dir         Directory to cache Ensembl GFFs (default: ./cache/ensembl)
        --cat_cache_dir             Directory to cache CAT GFFs (default: ./cache/cat)
        --assembly_reports_dir      Directory with assembly reports (optional)
        --comparison_script         Path to comparison Python script
        --max_assemblies            Limit number of assemblies to process (for testing)
        --gencode_version           GENCODE version for GRCh38 reference (default: 47)
        --gencode_gtf               Path to pre-downloaded GENCODE GTF (skips download)
        --mapping_stats_genes       Path to CSV with Ensembl gene projection rates
        --mapping_stats_transcripts Path to CSV with Ensembl transcript projection rates
========================================================================================
*/

// Default parameters - MUST be defined before includes to avoid warnings
params.input = null
params.ensg_lookup = null
params.outdir = './results'
params.ensembl_cache_dir = './cache/ensembl'
params.cat_cache_dir = './cache/cat'
params.assembly_reports_dir = './cache/assembly_reports'
params.comparison_script = "${projectDir}/../../../packages/hprc-qc-annotation/src/hprc_qc_annotation/hprc_ensembl_cat_overlap.py"
params.max_assemblies = null
params.help = false

// GENCODE reference for GRCh38 divergence analysis
params.gencode_version = '47'
params.gencode_gtf = null

// Ensembl projection rate data (pre-computed)
params.mapping_stats_genes = null
params.mapping_stats_transcripts = null

include { FETCH_FILES } from './subworkflows/fetch_files/main'

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
        --ensg_lookup         Path to transcript ID to ENSG ID lookup table

    Optional arguments:
        --ensembl_cache_dir   Cache directory for Ensembl GFFs (default: ${params.ensembl_cache_dir})
        --cat_cache_dir       Cache directory for CAT GFFs (default: ${params.cat_cache_dir})
        --assembly_reports_dir  Directory containing assembly reports (optional)
        --comparison_script   Path to comparison script (default: auto-detected)
        --max_assemblies      Limit number of assemblies (for testing, e.g., --max_assemblies 3)

    GENCODE/GRCh38 reference:
        --gencode_version     GENCODE version number (default: 47)
        --gencode_gtf         Path to pre-downloaded GENCODE GTF (skips download)

    Projection rates (SUPP-A):
        --mapping_stats_genes       CSV with per-assembly gene projection rates
        --mapping_stats_transcripts CSV with per-assembly transcript projection rates

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
        Using GENCODE v${params.gencode_version} as GRCh38 reference
        ============================================
        """
    }


    // Run comparison workflow
    FETCH_FILES(assemblies_ch)

    // Cache-only summary
    FETCH_FILES.out.paired.count().view { n -> "Cache population complete for ${n} pairs" }

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
