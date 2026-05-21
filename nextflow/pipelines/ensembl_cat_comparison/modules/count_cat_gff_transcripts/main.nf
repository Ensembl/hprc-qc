process COUNT_CAT_GFF_TRANSCRIPTS {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'conda-forge::python=3.11'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(cat_gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("*_cat_gene_transcript_counts.tsv"), emit: cat_counts

    script:
    """
    count_gff_transcripts.py \
        --gff ${cat_gff} \
        --output ${assembly_accession}_cat_gene_transcript_counts.tsv \
        --assembly-accession ${assembly_accession}
    """

    stub:
    """
    touch ${assembly_accession}_cat_gene_transcript_counts.tsv
    """
}
