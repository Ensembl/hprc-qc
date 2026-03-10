process MULTI_MAPPING {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'conda-forge::python=3.11'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(all_pairs)

    output:
    tuple val(assembly_accession), val(sample_name), path("*_multi_mapping.tsv"), emit: metrics

    script:
    """
    analyze_multi_mapping.py \\
        --all-pairs ${all_pairs} \\
        --output ${assembly_accession}_multi_mapping.tsv \\
        --assembly-accession ${assembly_accession} \\
        --sample-name ${sample_name}
    """

    stub:
    """
    touch ${assembly_accession}_multi_mapping.tsv
    """
}
