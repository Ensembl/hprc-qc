include { FETCH_FILES }          from '../subworkflows/fetch_files/main'
include { FETCH_GENCODE_GTF }   from '../modules/fetch_gencode_gtf/main'
include { RENAME_ENSEMBL_GFF }  from '../modules/rename_ensembl_gff/main'
include { FILTER_CAT_GFF }      from '../modules/filter_cat_gff/main'
include { RUN_COMPARISON }      from '../modules/run_comparison/main'
include { QC_METRICS }          from '../subworkflows/qc_metrics/main'
include { AGGREGATE }           from '../subworkflows/aggregate/main'

workflow ENSEMBL_CAT_COMPARISON {
    take:
    assemblies_ch  // Channel of [assembly_accession, sample_name]
    ensg_lookup    // Channel with path to ENSG lookup file

    main:
    // Fetch GENCODE GTF only when GRCh38 divergence is requested.
    if (params.run_grch38_divergence) {
        gencode_ch = Channel.value(params.gencode_version)
        if (params.gencode_gtf) {
            gencode_gtf_ch = Channel.fromPath(params.gencode_gtf, checkIfExists: true)
        } else {
            FETCH_GENCODE_GTF(gencode_ch)
            gencode_gtf_ch = FETCH_GENCODE_GTF.out.gtf
        }
    } else {
        gencode_gtf_ch = Channel.empty()
    }

    // Fetch GFF files and assembly reports
    FETCH_FILES(assemblies_ch)
    // Rename Ensembl GFF chromosomes to match CAT (GenBank accessions)
    // This is needed for GFFcompare, track hubs, and any tool that compares
    // the two GFFs by chromosome name. The assembly report maps
    // Assigned-Molecule (1, X) -> GenBank-Accn (CM089370.1, CM089395.1).
    FETCH_FILES.out.paired
        .map { accession, sample, ensembl_gff, cat_gff, assembly_report ->
            tuple(accession, sample, ensembl_gff, assembly_report)
        }
        .set { rename_input }
    RENAME_ENSEMBL_GFF(rename_input)

    // Build channel with renamed Ensembl GFF replacing original
    renamed_paired = RENAME_ENSEMBL_GFF.out.gff
        .join(
            FETCH_FILES.out.paired.map { acc, sample, ens, cat, report ->
                tuple(acc, sample, cat)
            },
            by: [0, 1]
        )
        .map { acc, sample, renamed_ens, cat ->
            tuple(acc, sample, renamed_ens, cat)
        }

    // Optionally filter CAT GFF before comparison, e.g. for no-Kinnex/no-ULC
    // sensitivity runs. The filtered CAT GFF then goes through the normal RBH
    // comparison and QC path.
    if (params.filter_cat_gff) {
        FILTER_CAT_GFF(
            renamed_paired.map { acc, sample, renamed_ens, cat ->
                tuple(acc, sample, cat)
            }
        )
        analysis_paired = renamed_paired
            .map { acc, sample, renamed_ens, cat -> tuple(acc, sample, renamed_ens) }
            .join(FILTER_CAT_GFF.out.gff, by: [0, 1])
            .map { acc, sample, renamed_ens, filtered_cat ->
                tuple(acc, sample, renamed_ens, filtered_cat)
            }
        filtered_cat_gffs_ch = FILTER_CAT_GFF.out.gff
        filtered_cat_audits_ch = FILTER_CAT_GFF.out.audit
    } else {
        analysis_paired = renamed_paired
        filtered_cat_gffs_ch = Channel.empty()
        filtered_cat_audits_ch = Channel.empty()
    }

    // RUN_COMPARISON expects [accession, sample, ensembl_gff, cat_gff].
    RUN_COMPARISON(analysis_paired)

    // Prepare data for QC metrics using renamed Ensembl GFF
    // Combine GFF files with comparison results
    qc_input = analysis_paired
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

    // Run QC metrics analyses (now includes GRCh38 divergence)
    QC_METRICS(qc_input, ensg_lookup, gencode_gtf_ch)

    gene_presence_summaries_ch = Channel.empty()
    sankey_summaries_ch = Channel.empty()
    coding_integrity_summaries_ch = Channel.empty()
    transcript_count_summaries_ch = Channel.empty()
    divergence_summaries_ch = Channel.empty()
    intron_chain_biotype_summaries_ch = Channel.empty()

    if (params.run_aggregate) {
        // Aggregate selected per-assembly results into intermediate spreadsheets.
        collected_gene_presence = QC_METRICS.out.gene_presence
            .map { accession, sample, tsv -> tsv }
            .collect()

        collected_rbh = RUN_COMPARISON.out.rbh
            .map { accession, tsv -> tsv }
            .collect()

        collected_transcript = QC_METRICS.out.transcript_concordance
            .map { accession, sample, tsv -> tsv }
            .collect()

        collected_coding = QC_METRICS.out.coding_integrity
            .map { accession, sample, tsv -> tsv }
            .collect()

        collected_divergence = QC_METRICS.out.grch38_divergence
            .map { accession, sample, tsv -> tsv }
            .collect()

        collected_multi_mapping = QC_METRICS.out.multi_mapping
            .map { accession, sample, tsv -> tsv }
            .collect()

        collected_gene_counts = QC_METRICS.out.gene_transcript_counts
            .map { accession, sample, tsv -> tsv }
            .collect()

        collected_cat_gene_counts = QC_METRICS.out.cat_gene_transcript_counts
            .map { accession, sample, tsv -> tsv }
            .collect()

        collected_gff_feature_counts = QC_METRICS.out.gff_feature_counts
            .map { accession, sample, tsv -> tsv }
            .collect()
        collected_gff_gene_metrics = QC_METRICS.out.gff_gene_metrics
            .map { accession, sample, tsv -> tsv }
            .collect()
        collected_gff_tx_metrics = QC_METRICS.out.gff_tx_metrics
            .map { accession, sample, tsv -> tsv }
            .collect()

        AGGREGATE(
            collected_gene_presence,
            collected_rbh,
            collected_transcript,
            collected_coding,
            collected_divergence,
            collected_multi_mapping,
            collected_gene_counts,
            collected_cat_gene_counts,
            collected_gff_feature_counts,
            collected_gff_gene_metrics,
            collected_gff_tx_metrics
        )

        gene_presence_summaries_ch = AGGREGATE.out.gene_presence_summaries
        sankey_summaries_ch = AGGREGATE.out.sankey_summaries
        coding_integrity_summaries_ch = AGGREGATE.out.coding_integrity_summaries
        transcript_count_summaries_ch = AGGREGATE.out.transcript_count_summaries
        divergence_summaries_ch = AGGREGATE.out.divergence_summaries
        intron_chain_biotype_summaries_ch = AGGREGATE.out.intron_chain_biotype_summaries
    }

    emit:
    renamed_ensembl_gffs = RENAME_ENSEMBL_GFF.out.gff
    filtered_cat_gffs = filtered_cat_gffs_ch
    filtered_cat_audits = filtered_cat_audits_ch
    rbh = RUN_COMPARISON.out.rbh
    all_pairs = RUN_COMPARISON.out.all_pairs
    logs = RUN_COMPARISON.out.log
    transcript_concordance = QC_METRICS.out.transcript_concordance
    coding_integrity = QC_METRICS.out.coding_integrity
    gene_presence = QC_METRICS.out.gene_presence
    multi_mapping = QC_METRICS.out.multi_mapping
    grch38_divergence = QC_METRICS.out.grch38_divergence
    gene_transcript_counts = QC_METRICS.out.gene_transcript_counts
    gff_feature_counts = QC_METRICS.out.gff_feature_counts
    gff_gene_metrics   = QC_METRICS.out.gff_gene_metrics
    gff_tx_metrics     = QC_METRICS.out.gff_tx_metrics
    // Aggregated intermediate spreadsheets
    gene_presence_summaries = gene_presence_summaries_ch
    sankey_summaries = sankey_summaries_ch
    coding_integrity_summaries = coding_integrity_summaries_ch
    transcript_count_summaries = transcript_count_summaries_ch
    divergence_summaries = divergence_summaries_ch
    intron_chain_biotype_summaries = intron_chain_biotype_summaries_ch
}
