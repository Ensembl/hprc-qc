process GFFCOMPARE_PARSE {
    tag "${assembly_accession}_${sample_name}_${direction}"
    label 'process_low'
    conda 'conda-forge::python=3.11'
    container 'python:3.11-slim'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), val(direction), path(tmap)

    output:
    tuple val(assembly_accession), val(sample_name), val(direction), path("gffcompare/class_counts_${direction}.tsv"), emit: counts

    script:
    """
    set -euo pipefail
    mkdir -p gffcompare
    parse_gffcompare_tmap.py \
      --tmap "${tmap}" \
      --assembly-accession "${assembly_accession}" \
      --sample-name "${sample_name}" \
      --direction "${direction}" \
      --output "gffcompare/class_counts_${direction}.tsv"
    """

    stub:
    """
    mkdir -p gffcompare
    printf "assembly_accession\tsample_name\tdirection\tclass_code\tn_transcripts\tdenominator\tpct\n" > "gffcompare/class_counts_${direction}.tsv"
    """
}
