process GRCH38_DIVERGENCE {
    tag "${assembly_accession}_${sample_name}"
    label 'process_medium'
    conda 'conda-forge::python=3.11'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(ensembl_gff), path(cat_gff), path(rbh_pairs)
    path(gencode_gtf)

    output:
    tuple val(assembly_accession), val(sample_name), path("*_grch38_divergence.tsv"), emit: metrics

    script:
    """
    compute_grch38_divergence.py \
        --ensembl-gff ${ensembl_gff} \
        --cat-gff ${cat_gff} \
        --rbh-pairs ${rbh_pairs} \
        --gencode-gtf ${gencode_gtf} \
        --assembly-accession ${assembly_accession} \
        --sample-name ${sample_name} \
        --output ${assembly_accession}_grch38_divergence.tsv
    """

    stub:
    """
    touch ${assembly_accession}_grch38_divergence.tsv
    """
}
