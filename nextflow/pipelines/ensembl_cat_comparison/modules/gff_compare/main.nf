process GFF_COMPARE {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'conda-forge::python=3.11 conda-forge::pandas bioconda::gffcompare bioconda::gffread'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(ensembl_gff), path(cat_gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("gff_compare/feature_counts.tsv"), emit: features
    tuple val(assembly_accession), val(sample_name), path("gff_compare/gene_metrics.tsv"), emit: gene_metrics
    tuple val(assembly_accession), val(sample_name), path("gff_compare/tx_metrics.tsv"), emit: tx_metrics

    script:
    """
    mkdir -p gff_compare
    compare_gff_features.py \
        --ensembl-gff ${ensembl_gff} \
        --cat-gff ${cat_gff} \
        --output-dir gff_compare \
        --assembly-accession ${assembly_accession} \
        --sample-name ${sample_name}
    """

    stub:
    """
    mkdir -p gff_compare
    touch gff_compare/feature_counts.tsv
    touch gff_compare/gene_metrics.tsv
    touch gff_compare/tx_metrics.tsv
    """
}
