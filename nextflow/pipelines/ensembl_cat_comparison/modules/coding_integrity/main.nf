process CODING_INTEGRITY {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'conda-forge::python=3.11'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(ensembl_gff), path(cat_gff), path(rbh_pairs)

    output:
    tuple val(assembly_accession), val(sample_name), path("*_coding_integrity.tsv"), emit: metrics

    script:
    """
    assess_coding_integrity.py \\
        --ensembl-gff ${ensembl_gff} \\
        --cat-gff ${cat_gff} \\
        --rbh-pairs ${rbh_pairs} \\
        --output ${assembly_accession}_coding_integrity.tsv \\
        --assembly-accession ${assembly_accession} \\
        --sample-name ${sample_name}
    """

    stub:
    """
    touch ${assembly_accession}_coding_integrity.tsv
    """
}
