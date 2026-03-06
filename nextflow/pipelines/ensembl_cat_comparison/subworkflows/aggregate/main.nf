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
    gene_transcript_count_files  // Channel of gene_transcript_counts TSV files (collected)

    main:
    // Run all aggregations in parallel
    AGGREGATE_GENE_PRESENCE(gene_presence_files)

    AGGREGATE_SANKEY(
        gene_presence_files,
        rbh_files,
        transcript_concordance_files,
        coding_integrity_files
    )

    AGGREGATE_CODING_INTEGRITY(coding_integrity_files)

    AGGREGATE_TRANSCRIPT_COUNTS(transcript_concordance_files)

    AGGREGATE_INTRON_CHAIN_BY_BIOTYPE(
        transcript_concordance_files,
        gene_transcript_count_files
    )

    AGGREGATE_DIVERGENCE(divergence_files)

    AGGREGATE_CONCORDANCE_VS_REF(
        rbh_files,
        transcript_concordance_files,
        coding_integrity_files,
        divergence_files,
        multi_mapping_files,
    )

    emit:
    gene_presence_summaries      = AGGREGATE_GENE_PRESENCE.out.summaries
    sankey_summaries             = AGGREGATE_SANKEY.out.summaries
    coding_integrity_summaries   = AGGREGATE_CODING_INTEGRITY.out.summaries
    transcript_count_summaries   = AGGREGATE_TRANSCRIPT_COUNTS.out.summaries
    divergence_summaries         = AGGREGATE_DIVERGENCE.out.summaries
    concordance_vs_ref_summaries = AGGREGATE_CONCORDANCE_VS_REF.out.summaries
    intron_chain_biotype_summaries = AGGREGATE_INTRON_CHAIN_BY_BIOTYPE.out.summaries
}
