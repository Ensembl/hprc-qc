process GENE_PRESENCE {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'conda-forge::python=3.11'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(ensembl_gff), path(cat_gff), path(ensg_lookup)

    output:
    tuple val(assembly_accession), val(sample_name), path("*_gene_presence.tsv"), emit: metrics

    script:
    """
    compare_gene_presence.py \\
        --ensembl-gff ${ensembl_gff} \\
        --cat-gff ${cat_gff} \\
        --output ${assembly_accession}_gene_presence.tsv \\
        --assembly-accession ${assembly_accession} \\
        --sample-name ${sample_name} \\
        --ensg-lookup ${ensg_lookup}
    """

    stub:
    """
    touch ${assembly_accession}_gene_presence.tsv
    """
}
