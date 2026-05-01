include { AGGREGATE_GENE_PRESENCE }      from '../../modules/aggregate_results/main'
include { AGGREGATE_SANKEY }              from '../../modules/aggregate_results/main'
include { AGGREGATE_CODING_INTEGRITY }    from '../../modules/aggregate_results/main'
include { AGGREGATE_TRANSCRIPT_COUNTS }   from '../../modules/aggregate_results/main'
include { AGGREGATE_DIVERGENCE }          from '../../modules/aggregate_results/main'
include { AGGREGATE_CONCORDANCE_VS_REF }  from '../../modules/aggregate_results/main'
include { AGGREGATE_INTRON_CHAIN_BY_BIOTYPE } from '../../modules/aggregate_results/main'

workflow AGGREGATE {
    take:
    gene_presence_files          // Channel of gene_presence TSV files (collected)
    rbh_files                    // Channel of RBH TSV files (collected)
    transcript_concordance_files // Channel of transcript_concordance TSV files (collected)
    coding_integrity_files       // Channel of coding_integrity TSV files (collected)
    divergence_files             // Channel of grch38_divergence TSV files (collected)
    multi_mapping_files          // Channel of multi_mapping TSV files (collected)
    gene_transcript_count_files  // Channel of Ensembl gene_transcript_counts TSV files (collected)
    cat_gene_transcript_count_files // Channel of CAT gene_transcript_counts TSV files (collected)
    gff_feature_counts // Channel of per-assembly feature_counts.tsv (collected)
    gff_gene_metrics   // Channel of per-assembly gene_metrics.tsv (collected)
    gff_tx_metrics     // Channel of per-assembly tx_metrics.tsv (collected)

    main:
    // Run selected aggregations in parallel. Defaults run everything.
    if (params.aggregate_gene_presence) {
        AGGREGATE_GENE_PRESENCE(gene_presence_files)
        gene_presence_summaries_ch = AGGREGATE_GENE_PRESENCE.out.summaries
    } else {
        gene_presence_summaries_ch = Channel.empty()
    }

    if (params.aggregate_sankey) {
        AGGREGATE_SANKEY(
            gene_presence_files,
            rbh_files,
            transcript_concordance_files,
            coding_integrity_files
        )
        sankey_summaries_ch = AGGREGATE_SANKEY.out.summaries
    } else {
        sankey_summaries_ch = Channel.empty()
    }

    if (params.aggregate_coding_integrity) {
        AGGREGATE_CODING_INTEGRITY(coding_integrity_files)
        coding_integrity_summaries_ch = AGGREGATE_CODING_INTEGRITY.out.summaries
    } else {
        coding_integrity_summaries_ch = Channel.empty()
    }

    if (params.aggregate_transcript_counts) {
        AGGREGATE_TRANSCRIPT_COUNTS(transcript_concordance_files)
        transcript_count_summaries_ch = AGGREGATE_TRANSCRIPT_COUNTS.out.summaries
    } else {
        transcript_count_summaries_ch = Channel.empty()
    }

    if (params.aggregate_intron_chain_by_biotype) {
        AGGREGATE_INTRON_CHAIN_BY_BIOTYPE(
            transcript_concordance_files,
            gene_transcript_count_files,
            cat_gene_transcript_count_files
        )
        intron_chain_biotype_summaries_ch = AGGREGATE_INTRON_CHAIN_BY_BIOTYPE.out.summaries
    } else {
        intron_chain_biotype_summaries_ch = Channel.empty()
    }

    // No aggregator yet for gff_compare: the notebook will read them directly from qc_metrics

    if (params.aggregate_divergence) {
        AGGREGATE_DIVERGENCE(divergence_files)
        divergence_summaries_ch = AGGREGATE_DIVERGENCE.out.summaries
    } else {
        divergence_summaries_ch = Channel.empty()
    }

    if (params.aggregate_concordance_vs_ref) {
        AGGREGATE_CONCORDANCE_VS_REF(
            rbh_files,
            transcript_concordance_files,
            coding_integrity_files,
            divergence_files,
            multi_mapping_files,
        )
        concordance_vs_ref_summaries_ch = AGGREGATE_CONCORDANCE_VS_REF.out.summaries
    } else {
        concordance_vs_ref_summaries_ch = Channel.empty()
    }

    emit:
    gene_presence_summaries      = gene_presence_summaries_ch
    sankey_summaries             = sankey_summaries_ch
    coding_integrity_summaries   = coding_integrity_summaries_ch
    transcript_count_summaries   = transcript_count_summaries_ch
    divergence_summaries         = divergence_summaries_ch
    concordance_vs_ref_summaries = concordance_vs_ref_summaries_ch
    intron_chain_biotype_summaries = intron_chain_biotype_summaries_ch
}
