process GFFCOMPARE_RUN {
    tag "${assembly_accession}_${sample_name}_${direction}"
    label 'process_low'
    conda 'bioconda::gffcompare'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), val(direction), path(ref_gtf), path(qry_gtf)

    output:
    tuple val(assembly_accession), val(sample_name), val(direction), path("gffcompare/tmap_${direction}.tsv"), emit: tmap
    tuple val(assembly_accession), val(sample_name), val(direction), path("gffcompare/stats_${direction}.txt"), emit: stats

    script:
    """
    mkdir -p gffcompare
    gffcompare -r ${ref_gtf} -o ${direction} ${qry_gtf} || true
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
