include { FETCH_FILES }      from '../subworkflows/fetch_files/main'
include { RUN_COMPARISON }   from '../modules/run_comparison/main'
include { QC_METRICS }       from '../subworkflows/qc_metrics/main'

workflow ENSEMBL_CAT_COMPARISON {
    take:
    assemblies_ch  // Channel of [assembly_accession, sample_name]

    main:
    // Fetch GFF files and assembly reports
    FETCH_FILES(assemblies_ch)

    // FETCH_FILES returns [accession, sample, ensembl_gff, cat_gff, assembly_report]
    // Pass directly to RUN_COMPARISON which expects the same tuple structure
    RUN_COMPARISON(FETCH_FILES.out.paired)

    // Prepare data for QC metrics
    // Combine GFF files with comparison results
    // FETCH_FILES.out.paired: [accession, sample, ensembl_gff, cat_gff, assembly_report]
    // RUN_COMPARISON.out.rbh: [accession, rbh.tsv]
    // RUN_COMPARISON.out.all_pairs: [accession, all_pairs.tsv]

    qc_input = FETCH_FILES.out.paired
        .map { accession, sample, ensembl_gff, cat_gff, assembly_report ->
            tuple(accession, sample, ensembl_gff, cat_gff)
        }
        .join(
            RUN_COMPARISON.out.rbh.map { accession, rbh -> tuple(accession, rbh) },
            by: 0
        )
        .join(
            RUN_COMPARISON.out.all_pairs.map { accession, all_pairs -> tuple(accession, all_pairs) },
            by: 0
        )
        .map { accession, sample, ensembl_gff, cat_gff, rbh_tsv, all_pairs_tsv ->
            tuple(accession, sample, ensembl_gff, cat_gff, rbh_tsv, all_pairs_tsv)
        }

    // Run QC metrics analyses
    QC_METRICS(qc_input)

    emit:
    rbh = RUN_COMPARISON.out.rbh
    all_pairs = RUN_COMPARISON.out.all_pairs
    logs = RUN_COMPARISON.out.log
    transcript_concordance = QC_METRICS.out.transcript_concordance
    coding_integrity = QC_METRICS.out.coding_integrity
    gene_presence = QC_METRICS.out.gene_presence
    multi_mapping = QC_METRICS.out.multi_mapping
}
