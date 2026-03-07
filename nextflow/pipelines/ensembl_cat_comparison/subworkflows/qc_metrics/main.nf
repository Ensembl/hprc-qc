include { TRANSCRIPT_CONCORDANCE } from '../../modules/transcript_concordance/main'
include { CODING_INTEGRITY }      from '../../modules/coding_integrity/main'
include { GENE_PRESENCE }         from '../../modules/gene_presence/main'
include { MULTI_MAPPING }         from '../../modules/multi_mapping/main'
include { GRCH38_DIVERGENCE }     from '../../modules/grch38_divergence/main'
include { COUNT_GFF_TRANSCRIPTS } from '../../modules/count_gff_transcripts/main'
include { COUNT_CAT_GFF_TRANSCRIPTS } from '../../modules/count_cat_gff_transcripts/main'

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
        }
        .set { qc_inputs }

    // Run QC analyses in parallel
    TRANSCRIPT_CONCORDANCE(qc_inputs.transcript)
    CODING_INTEGRITY(qc_inputs.coding)
    GENE_PRESENCE(qc_inputs.presence.combine(ensg_lookup))
    MULTI_MAPPING(qc_inputs.multi)
    GRCH38_DIVERGENCE(qc_inputs.divergence, gencode_gtf.first())
    COUNT_GFF_TRANSCRIPTS(qc_inputs.gff_only)
    COUNT_CAT_GFF_TRANSCRIPTS(qc_inputs.cat_gff_only)

    emit:
    transcript_concordance = TRANSCRIPT_CONCORDANCE.out.metrics
    coding_integrity = CODING_INTEGRITY.out.metrics
    gene_presence = GENE_PRESENCE.out.metrics
    multi_mapping = MULTI_MAPPING.out.metrics
    grch38_divergence = GRCH38_DIVERGENCE.out.metrics
    gene_transcript_counts = COUNT_GFF_TRANSCRIPTS.out.counts
    cat_gene_transcript_counts = COUNT_CAT_GFF_TRANSCRIPTS.out.cat_counts
}
