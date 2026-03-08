process GFFCOMPARE_RUN {
    tag "${assembly_accession}_${sample_name}_${direction}"
    label 'process_low'
    conda 'bioconda::gffcompare'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), val(direction), path(ref_gtf, stageAs: 'ref.gtf'), path(qry_gtf, stageAs: 'qry.gtf')

    output:
    tuple val(assembly_accession), val(sample_name), val(direction), path("gffcompare/tmap_${direction}.tsv"), emit: tmap
    tuple val(assembly_accession), val(sample_name), val(direction), path("gffcompare/stats_${direction}.txt"), emit: stats

    script:
    """
    set -euo pipefail
    mkdir -p gffcompare
    # Run gffcompare with the reference and query GTFs, using direction as prefix
    gffcompare -r ref.gtf -o "${direction}" qry.gtf || true

    # Normalize outputs and copy to standardized paths
    if [ -f "${direction}.tmap" ]; then
        cp "${direction}.tmap" "gffcompare/tmap_${direction}.tsv"
    else
        : > "gffcompare/tmap_${direction}.tsv"
    fi

    if [ -f "${direction}.stats" ]; then
        cp "${direction}.stats" "gffcompare/stats_${direction}.txt"
    else
        : > "gffcompare/stats_${direction}.txt"
    fi
    """

    stub:
    """
    mkdir -p gffcompare
    : > "gffcompare/tmap_${direction}.tsv"
    : > "gffcompare/stats_${direction}.txt"
    """
}
