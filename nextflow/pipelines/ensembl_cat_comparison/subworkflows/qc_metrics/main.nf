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
include { GFFCOMPARE_RUN }   from '../../modules/gffcompare_run/main'
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

    // Run QC analyses in parallel
    TRANSCRIPT_CONCORDANCE(qc_inputs.transcript)
    CODING_INTEGRITY(qc_inputs.coding)
    GENE_PRESENCE(qc_inputs.presence.combine(ensg_lookup))
    MULTI_MAPPING(qc_inputs.multi)
    GRCH38_DIVERGENCE(qc_inputs.divergence, gencode_gtf.first())
    COUNT_GFF_TRANSCRIPTS(qc_inputs.gff_only)
    COUNT_CAT_GFF_TRANSCRIPTS(qc_inputs.cat_gff_only)
    GFF_COMPARE(qc_inputs.gff_pair)

    // Convert both GFFs to GTF for gffcompare compatibility
    ENS_GFF_TO_GTF(qc_inputs.gff_only)
    ENS_GFF_TO_GTF.out.gtf.set { ens_gtf }

    CAT_GFF_TO_GTF(qc_inputs.cat_gff_only)
    CAT_GFF_TO_GTF.out.gtf.set { cat_gtf }

    // Join Ensembl and CAT GTFs by (assembly_accession, sample_name)
    ens_gtf
        .map { acc, sample, gtf -> tuple([acc, sample], gtf) }
        .set { ens_gtf_kv }

    cat_gtf
        .map { acc, sample, gtf -> tuple([acc, sample], gtf) }
        .set { cat_gtf_kv }

    ens_gtf_kv
        .join(cat_gtf_kv)
        .map { key, ens_path, cat_path -> tuple(key[0], key[1], ens_path, cat_path) }
        .set { gtf_pairs }

    // Create two direction jobs per pair
    gtf_pairs
        .flatMap { acc, sample, ens_path, cat_path ->
            [
                tuple(acc, sample, 'Ensembl_to_CAT', ens_path, cat_path),
                tuple(acc, sample, 'CAT_to_Ensembl', cat_path, ens_path)
            ]
        }
        .set { gffcmp_jobs }

    // Run gffcompare per direction and parse class_code distributions
    GFFCOMPARE_RUN(gffcmp_jobs)
    GFFCOMPARE_PARSE(GFFCOMPARE_RUN.out.tmap)

    emit:
    transcript_concordance = TRANSCRIPT_CONCORDANCE.out.metrics
    coding_integrity = CODING_INTEGRITY.out.metrics
    gene_presence = GENE_PRESENCE.out.metrics
    multi_mapping = MULTI_MAPPING.out.metrics
    grch38_divergence = GRCH38_DIVERGENCE.out.metrics
    gene_transcript_counts = COUNT_GFF_TRANSCRIPTS.out.counts
    cat_gene_transcript_counts = COUNT_CAT_GFF_TRANSCRIPTS.out.cat_counts
    gff_feature_counts = GFF_COMPARE.out.features
    gff_gene_metrics   = GFF_COMPARE.out.gene_metrics
    gff_tx_metrics     = GFF_COMPARE.out.tx_metrics
    gffcompare_tmap         = GFFCOMPARE_RUN.out.tmap
    gffcompare_stats        = GFFCOMPARE_RUN.out.stats
    gffcompare_class_counts = GFFCOMPARE_PARSE.out.counts
}
