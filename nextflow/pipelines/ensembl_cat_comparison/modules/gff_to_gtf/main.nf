process GFF_TO_GTF {
    tag "${assembly_accession}_${sample_name}"
    // GFF3 -> GTF can be memory-heavy; use medium resources
    label 'process_medium'
    conda 'bioconda::gffread=0.12.7'
    container 'quay.io/biocontainers/gffread:0.12.7--hd03093a_1'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy', overwrite: true

    input:
    tuple val(assembly_accession), val(sample_name), path(gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("gff_to_gtf/${assembly_accession}.${sample_name}.gtf"), emit: gtf

    script:
    """
    mkdir -p gff_to_gtf
    # Use a temporary uncompressed file to reduce peak memory vs pipe on some systems
    gzip -dc "${gff}" > gff_to_gtf/in.gff3 || cp -f "${gff}" gff_to_gtf/in.gff3
    gffread -E -F -T -o "gff_to_gtf/${assembly_accession}.${sample_name}.gtf" gff_to_gtf/in.gff3
    rm -f gff_to_gtf/in.gff3
    """

    stub:
    """
    mkdir -p gff_to_gtf
    touch "gff_to_gtf/${assembly_accession}.${sample_name}.gtf"
    """
}
