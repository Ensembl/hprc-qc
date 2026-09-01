include { TRANSCRIPT_CONCORDANCE } from '../../modules/transcript_concordance/main'
include { CODING_INTEGRITY }      from '../../modules/coding_integrity/main'
include { GENE_PRESENCE }         from '../../modules/gene_presence/main'
include { MULTI_MAPPING }         from '../../modules/multi_mapping/main'
include { GRCH38_DIVERGENCE }     from '../../modules/grch38_divergence/main'
include { COUNT_GFF_TRANSCRIPTS } from '../../modules/count_gff_transcripts/main'
include { GFF_COMPARE }           from '../../modules/gff_compare/main'
include { COUNT_CAT_GFF_TRANSCRIPTS } from '../../modules/count_cat_gff_transcripts/main'
include { GFF_TO_GTF as ENS_GFF_TO_GTF }       from '../../modules/gff_to_gtf/main'
include { GFF_TO_GTF as CAT_GFF_TO_GTF }       from '../../modules/gff_to_gtf/main'
include { GFFCOMPARE_RUN }   from '../../modules/gffcompare_run2/main'
include { GFFCOMPARE_PARSE } from '../../modules/gffcompare_parse/main'

workflow QC_METRICS {
    take:
    comparison_results  // Channel of [assembly_accession, sample_name, ensembl_gff, cat_gff, rbh_tsv, all_pairs_tsv]
    ensg_lookup         // Channel with path to ENSG lookup file
    gencode_gtf         // Channel with path to GENCODE GTF (for GRCh38 divergence)

    main:
    // Prepare inputs for each QC process
    comparison_results
        .multiMap { accession, sample, ensembl_gff, cat_gff, rbh_tsv, all_pairs_tsv ->
            transcript: tuple(accession, sample, ensembl_gff, cat_gff, rbh_tsv)
            coding: tuple(accession, sample, ensembl_gff, cat_gff, rbh_tsv)
            presence: tuple(accession, sample, ensembl_gff, cat_gff)
            multi: tuple(accession, sample, all_pairs_tsv)
            divergence: tuple(accession, sample, ensembl_gff, cat_gff, rbh_tsv)
            gff_only: tuple(accession, sample, ensembl_gff)
            cat_gff_only: tuple(accession, sample, cat_gff)
            gff_pair: tuple(accession, sample, ensembl_gff, cat_gff)
        }
        .set { qc_inputs }

    // Run selected QC analyses in parallel. Defaults run everything; sensitivity
    // runs can disable heavy endpoints and keep only the outputs needed.
    if (params.run_transcript_concordance) {
        TRANSCRIPT_CONCORDANCE(qc_inputs.transcript)
        transcript_concordance_ch = TRANSCRIPT_CONCORDANCE.out.metrics
    } else {
        transcript_concordance_ch = Channel.empty()
    }

    if (params.run_coding_integrity) {
        CODING_INTEGRITY(qc_inputs.coding)
        coding_integrity_ch = CODING_INTEGRITY.out.metrics
    } else {
        coding_integrity_ch = Channel.empty()
    }

    if (params.run_gene_presence.toString().toBoolean()) {
        GENE_PRESENCE(qc_inputs.presence.combine(ensg_lookup))
        gene_presence_ch = GENE_PRESENCE.out.metrics
    } else {
        gene_presence_ch = Channel.empty()
    }

    if (params.run_multi_mapping) {
        MULTI_MAPPING(qc_inputs.multi)
        multi_mapping_ch = MULTI_MAPPING.out.metrics
    } else {
        multi_mapping_ch = Channel.empty()
    }

    if (params.run_grch38_divergence.toString().toBoolean()) {
        GRCH38_DIVERGENCE(qc_inputs.divergence, gencode_gtf.first())
        grch38_divergence_ch = GRCH38_DIVERGENCE.out.metrics
    } else {
        grch38_divergence_ch = Channel.empty()
    }

    if (params.run_gene_transcript_counts) {
        COUNT_GFF_TRANSCRIPTS(qc_inputs.gff_only)
        gene_transcript_counts_ch = COUNT_GFF_TRANSCRIPTS.out.counts
    } else {
        gene_transcript_counts_ch = Channel.empty()
    }

    if (params.run_cat_gene_transcript_counts) {
        COUNT_CAT_GFF_TRANSCRIPTS(qc_inputs.cat_gff_only)
        cat_gene_transcript_counts_ch = COUNT_CAT_GFF_TRANSCRIPTS.out.cat_counts
    } else {
        cat_gene_transcript_counts_ch = Channel.empty()
    }

    if (params.run_gff_feature_metrics) {
        GFF_COMPARE(qc_inputs.gff_pair)
        gff_feature_counts_ch = GFF_COMPARE.out.features
        gff_gene_metrics_ch = GFF_COMPARE.out.gene_metrics
        gff_tx_metrics_ch = GFF_COMPARE.out.tx_metrics
    } else {
        gff_feature_counts_ch = Channel.empty()
        gff_gene_metrics_ch = Channel.empty()
        gff_tx_metrics_ch = Channel.empty()
    }

    // Optionally run gffcompare-based transcript overlap analysis (--run_gffcompare)
    if (params.run_gffcompare) {
        ENS_GFF_TO_GTF(qc_inputs.gff_only)
        CAT_GFF_TO_GTF(qc_inputs.cat_gff_only)

        ENS_GFF_TO_GTF.out.gtf
            .map { acc, sample, gtf -> tuple([acc, sample], gtf) }
            .join(CAT_GFF_TO_GTF.out.gtf.map { acc, sample, gtf -> tuple([acc, sample], gtf) })
            .map { key, ens_path, cat_path -> tuple(key[0], key[1], ens_path, cat_path) }
            .flatMap { acc, sample, ens_path, cat_path ->
                [
                    tuple(acc, sample, 'Ensembl_to_CAT', ens_path, cat_path),
                    tuple(acc, sample, 'CAT_to_Ensembl', cat_path, ens_path)
                ]
            }
            .set { gffcmp_jobs }

        GFFCOMPARE_RUN(gffcmp_jobs)
        GFFCOMPARE_PARSE(GFFCOMPARE_RUN.out.tmap)

        gffcompare_tmap_ch         = GFFCOMPARE_RUN.out.tmap
        gffcompare_stats_ch        = GFFCOMPARE_RUN.out.stats
        gffcompare_class_counts_ch = GFFCOMPARE_PARSE.out.counts
    } else {
        gffcompare_tmap_ch         = Channel.empty()
        gffcompare_stats_ch        = Channel.empty()
        gffcompare_class_counts_ch = Channel.empty()
    }

    emit:
    transcript_concordance     = transcript_concordance_ch
    coding_integrity           = coding_integrity_ch
    gene_presence              = gene_presence_ch
    multi_mapping              = multi_mapping_ch
    grch38_divergence          = grch38_divergence_ch
    gene_transcript_counts     = gene_transcript_counts_ch
    cat_gene_transcript_counts = cat_gene_transcript_counts_ch
    gff_feature_counts         = gff_feature_counts_ch
    gff_gene_metrics           = gff_gene_metrics_ch
    gff_tx_metrics             = gff_tx_metrics_ch
    gffcompare_tmap            = gffcompare_tmap_ch
    gffcompare_stats           = gffcompare_stats_ch
    gffcompare_class_counts    = gffcompare_class_counts_ch
}
