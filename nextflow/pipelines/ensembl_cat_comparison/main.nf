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
params.comparison_script = "${projectDir}/bin/hprc_ensembl_cat_overlap.py"
params.max_assemblies = null
params.help = false

// GENCODE reference for GRCh38 divergence analysis
params.gencode_version = '47'
params.gencode_gtf = null

// Ensembl projection rate data (pre-computed)
params.mapping_stats_genes = null
params.mapping_stats_transcripts = null

// Optional CAT GFF filtering branch, used for no-Kinnex/no-ULC sensitivity runs
params.filter_cat_gff = false
params.filter_cat_filter_name = 'filtered_cat'
params.filter_cat_exclude_biotypes = 'unknown_likely_coding'
params.filter_cat_exclude_sources = ''
params.filter_cat_exclude_transcript_modes = ''
params.filter_cat_exclude_gene_name_regex = ''
params.filter_cat_exclude_gene_name_file = ''
params.filter_cat_exclude_hyphenated_gene_names = false
params.filter_cat_exclude_attr_nonempty = ''
params.filter_cat_exclude_attr_equals = ''
params.filter_cat_exclude_attr_in = ''
params.filter_cat_exclude_attr_in_file = ''

// Optional QC endpoints. Keep defaults equivalent to the historical full run.
params.run_transcript_concordance = true
params.run_coding_integrity = true
params.run_gene_presence = true
params.run_multi_mapping = true
params.run_grch38_divergence = true
params.run_gene_transcript_counts = true
params.run_cat_gene_transcript_counts = true
params.run_gff_feature_metrics = true

// Optional aggregation endpoints. For Panel G no-ULC runs, only intron-chain
// aggregation is needed in addition to transcript concordance and count inputs.
params.run_aggregate = true
params.aggregate_gene_presence = true
params.aggregate_sankey = true
params.aggregate_coding_integrity = true
params.aggregate_transcript_counts = true
params.aggregate_divergence = true
params.aggregate_concordance_vs_ref = true
params.aggregate_intron_chain_by_biotype = true

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

    Sensitivity / endpoint controls:
        --filter_cat_gff true
        --filter_cat_filter_name no_unknown_likely_coding
        --filter_cat_exclude_biotypes unknown_likely_coding
        --filter_cat_exclude_hyphenated_gene_names true
        --filter_cat_exclude_attr_nonempty collapsed_gene_ids,collapsed_gene_names
        --filter_cat_exclude_attr_equals extra_paralog=True
        --filter_cat_exclude_attr_in transcript_class=paralog|possible_paralog
        --filter_cat_exclude_gene_name_file readthrough_gene_names.txt
        --filter_cat_exclude_attr_in_file source_gene_id=patch_gene_ids.txt
        --filter_cat_exclude_transcript_modes augPB,exRef,Liftoff
        --filter_cat_exclude_sources Liftoff
        --run_<endpoint> false       Disable selected QC endpoints for focused reruns
        --aggregate_<endpoint> false Disable selected aggregation endpoints

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
if (!params.ensg_lookup) {
    log.error "ERROR: --ensg_lookup parameter is required"
    helpMessage()
    exit 1
} else if (!file(params.ensg_lookup).exists()) {
    log.error "ERROR: ENSG lookup file not found at: ${params.ensg_lookup}"
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
        Using GENCODE v${params.gencode_version} as GRCh38 reference
        ============================================
        """
    }

    // Create channel for lookup file
    ensg_lookup_ch = Channel.fromPath(params.ensg_lookup, checkIfExists: true)

    // Copy mapping stats to output if provided
    if (params.mapping_stats_genes) {
        Channel.fromPath(params.mapping_stats_genes, checkIfExists: true)
            .set { mapping_genes_ch }
    }
    if (params.mapping_stats_transcripts) {
        Channel.fromPath(params.mapping_stats_transcripts, checkIfExists: true)
            .set { mapping_transcripts_ch }
    }

    // Run comparison workflow
    ENSEMBL_CAT_COMPARISON(assemblies_ch, ensg_lookup_ch)

    // Collect results
    ENSEMBL_CAT_COMPARISON.out.rbh.collect().view { results ->
        """
        ============================================
        Completed! RBH files created: ${results.size()}
        Results in: ${params.outdir}/results/
        Intermediate spreadsheets in: ${params.outdir}/intermediate_spreadsheets/
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
