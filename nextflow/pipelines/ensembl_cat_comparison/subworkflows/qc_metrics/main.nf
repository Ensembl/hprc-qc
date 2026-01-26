include { TRANSCRIPT_CONCORDANCE } from '../../modules/transcript_concordance/main'
include { CODING_INTEGRITY }      from '../../modules/coding_integrity/main'
include { GENE_PRESENCE }         from '../../modules/gene_presence/main'
include { MULTI_MAPPING }         from '../../modules/multi_mapping/main'

workflow QC_METRICS {
    take:
    comparison_results  // Channel of [assembly_accession, sample_name, ensembl_gff, cat_gff, rbh_tsv, all_pairs_tsv]

    main:
    // Prepare inputs for each QC process
    comparison_results
        .multiMap { accession, sample, ensembl_gff, cat_gff, rbh_tsv, all_pairs_tsv ->
            transcript: tuple(accession, sample, ensembl_gff, cat_gff, rbh_tsv)
            coding: tuple(accession, sample, ensembl_gff, cat_gff, rbh_tsv)
            presence: tuple(accession, sample, ensembl_gff, cat_gff)
            multi: tuple(accession, sample, all_pairs_tsv)
        }
        .set { qc_inputs }

    // Run QC analyses in parallel
    TRANSCRIPT_CONCORDANCE(qc_inputs.transcript)
    CODING_INTEGRITY(qc_inputs.coding)
    GENE_PRESENCE(qc_inputs.presence)
    MULTI_MAPPING(qc_inputs.multi)

    emit:
    transcript_concordance = TRANSCRIPT_CONCORDANCE.out.metrics
    coding_integrity = CODING_INTEGRITY.out.metrics
    gene_presence = GENE_PRESENCE.out.metrics
    multi_mapping = MULTI_MAPPING.out.metrics
}
