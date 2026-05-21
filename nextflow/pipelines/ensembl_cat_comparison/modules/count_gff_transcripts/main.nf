process COUNT_GFF_TRANSCRIPTS {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'conda-forge::python=3.11'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(ensembl_gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("*_gene_transcript_counts.tsv"), emit: counts

    script:
    """
    count_gff_transcripts.py \\
        --gff ${ensembl_gff} \\
        --output ${assembly_accession}_gene_transcript_counts.tsv \\
        --assembly-accession ${assembly_accession}
    """

    stub:
    """
    touch ${assembly_accession}_gene_transcript_counts.tsv
    """
}
