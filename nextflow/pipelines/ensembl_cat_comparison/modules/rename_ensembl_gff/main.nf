process RENAME_ENSEMBL_GFF {
    tag "${assembly_accession}"
    label 'process_low'
    conda 'conda-forge::python=3.11'

    publishDir "${params.outdir}/renamed_gffs/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(ensembl_gff), path(assembly_report)

    output:
    tuple val(assembly_accession), val(sample_name), path("${assembly_accession}.ensembl.renamed.gff3.gz"), emit: gff

    script:
    """
    rename_gff_chroms.py \\
        --gff ${ensembl_gff} \\
        --assembly-report ${assembly_report} \\
        --output ${assembly_accession}.ensembl.renamed.gff3.gz
    """

    stub:
    """
    cp ${ensembl_gff} ${assembly_accession}.ensembl.renamed.gff3.gz
    """
}
