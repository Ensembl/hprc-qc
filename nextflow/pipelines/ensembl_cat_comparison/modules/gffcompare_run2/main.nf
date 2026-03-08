process GFFCOMPARE_RUN {
    tag "${assembly_accession}_${sample_name}_${direction}"
    label 'process_low'
    conda 'bioconda::gffcompare'
    container 'quay.io/biocontainers/gffcompare:0.12.6--h4ac6f70_2'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), val(direction), path(ref_gtf, stageAs: 'ref.gtf'), path(qry_gtf, stageAs: 'qry.gtf')

    output:
    tuple val(assembly_accession), val(sample_name), val(direction), path("gffcompare/tmap_${direction}.tsv"), emit: tmap
    tuple val(assembly_accession), val(sample_name), val(direction), path("gffcompare/stats_${direction}.txt"), emit: stats

    script:
    """
    mkdir -p gffcompare
    gffcompare -r ref.gtf -o ${direction} qry.gtf || true
    cp ${direction}.qry.gtf.tmap gffcompare/tmap_${direction}.tsv 2>/dev/null || : > gffcompare/tmap_${direction}.tsv
    cp ${direction}.stats gffcompare/stats_${direction}.txt 2>/dev/null || : > gffcompare/stats_${direction}.txt
    """

    stub:
    """
    mkdir -p gffcompare
    : > gffcompare/tmap_${direction}.tsv
    : > gffcompare/stats_${direction}.txt
    """
}
