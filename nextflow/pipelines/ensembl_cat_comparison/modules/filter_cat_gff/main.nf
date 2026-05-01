process FILTER_CAT_GFF {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'conda-forge::python=3.11'

    publishDir "${params.outdir}/filtered_cat_gff/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(cat_gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("*filtered_cat.gff3.gz"), emit: gff

    script:
    def excludeArgs = params.filter_cat_exclude_biotypes
        .toString()
        .split(',')
        .collect { it.trim() }
        .findAll { it }
        .collect { "--exclude-biotype ${it}" }
        .join(' ')
    """
    filter_gff_by_biotype.py \\
        --input-gff ${cat_gff} \\
        --output-gff ${assembly_accession}.filtered_cat.gff3.gz \\
        ${excludeArgs}
    """

    stub:
    """
    touch ${assembly_accession}.filtered_cat.gff3.gz
    """
}
