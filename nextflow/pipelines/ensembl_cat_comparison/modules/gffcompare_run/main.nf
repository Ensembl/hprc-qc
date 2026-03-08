process GFFCOMPARE_RUN {
    tag "${assembly_accession}_${sample_name}_${direction}"
    label 'process_low'
    conda 'bioconda::gffcompare'
    // Containerized gffcompare; override with --gffcompare_container if needed
    container "quay.io/biocontainers/gffcompare:0.10.6--h2d50403_0"

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

    # Locate outputs robustly across gffcompare versions without relying on ls/head
    tmap_src=""
    for f in "${direction}"*.tmap *.tmap; do
      if [ -f "$f" ]; then tmap_src="$f"; break; fi
    done
    stats_src=""
    if [ -f "${direction}.stats" ]; then
      stats_src="${direction}.stats"
    else
      for f in *.stats; do
        if [ -f "$f" ]; then stats_src="$f"; break; fi
      done
    fi

    # Normalize outputs and copy to standardized paths
    if [ -n "${tmap_src:-}" ] && [ -f "$tmap_src" ]; then
        cp "$tmap_src" "gffcompare/tmap_${direction}.tsv"
    else
        : > "gffcompare/tmap_${direction}.tsv"
    fi

    if [ -n "${stats_src:-}" ] && [ -f "$stats_src" ]; then
        cp "$stats_src" "gffcompare/stats_${direction}.txt"
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
