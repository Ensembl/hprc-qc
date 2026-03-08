process GFF_TO_GTF {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'bioconda::gffread'
    container "quay.io/biocontainers/gffread:0.9.12--0"

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("gff_to_gtf/${assembly_accession}.${sample_name}.gtf"), emit: gtf

    script:
    """
    set -euo pipefail
    mkdir -p gff_to_gtf

    gzip -dc "${gff}" > input.gff3

    gffread -T -o "gff_to_gtf/${assembly_accession}.${sample_name}.gtf" input.gff3
    """

    stub:
    """
    mkdir -p gff_to_gtf
    touch "gff_to_gtf/${assembly_accession}.${sample_name}.gtf"
    """
}