process GFF_TO_GTF {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'bioconda::gffread'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("gff_to_gtf/*.gtf"), emit: gtf

    script:
    """
    mkdir -p gff_to_gtf
    # Convert GFF3 to GTF using gffread for gffcompare compatibility
    gffread -T -o "gff_to_gtf/${assembly_accession}.${sample_name}.gtf" "${gff}"
    """

    stub:
    """
    mkdir -p gff_to_gtf
    touch "gff_to_gtf/${assembly_accession}.${sample_name}.gtf"
    """
}
