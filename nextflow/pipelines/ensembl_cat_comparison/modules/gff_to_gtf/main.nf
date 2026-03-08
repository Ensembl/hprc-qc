process GFF_TO_GTF {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'bioconda::gffread'
    // Containerized gffread; override with --gffread_container if needed
    container "${params.gffread_container ?: 'docker://quay.io/biocontainers/gffread:0.12.7--hd03093a_1'}"

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy', overwrite: true

    input:
    tuple val(assembly_accession), val(sample_name), path(gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("gff_to_gtf/${assembly_accession}.${sample_name}.gtf"), emit: gtf

    script:
    """
    gff_to_gtf.sh "${gff}" "gff_to_gtf/${assembly_accession}.${sample_name}.gtf"
    """

    stub:
    """
    mkdir -p gff_to_gtf
    touch "gff_to_gtf/${assembly_accession}.${sample_name}.gtf"
    """
}
