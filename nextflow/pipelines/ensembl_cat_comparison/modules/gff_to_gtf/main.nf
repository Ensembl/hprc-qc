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
    set -euo pipefail
    mkdir -p gff_to_gtf

    in_gff="input.gff3"
    case "${gff}" in
        *.gz) gzip -dc "${gff}" > "$in_gff" ;;
        *)    cp -f "${gff}" "$in_gff" ;;
    esac

    # Convert GFF3 to GTF using gffread for gffcompare compatibility
    # -E: tolerate minor format issues; -F: keep attributes
    gffread -E -F -T -o "gff_to_gtf/${assembly_accession}.${sample_name}.gtf" "$in_gff"

    # Sanity report to help diagnose empty tmap cases
    exons=$(grep -c $'\texon\t' "gff_to_gtf/${assembly_accession}.${sample_name}.gtf" || true)
    tids=$(grep -c 'transcript_id "' "gff_to_gtf/${assembly_accession}.${sample_name}.gtf" || true)
    {
      echo "input: ${gff}"
      echo "gtf: gff_to_gtf/${assembly_accession}.${sample_name}.gtf"
      echo "exon_lines=${exons}"
      echo "transcript_id_tags=${tids}"
    } > "gff_to_gtf/${assembly_accession}.${sample_name}.convert.log"
    """

    stub:
    """
    mkdir -p gff_to_gtf
    touch "gff_to_gtf/${assembly_accession}.${sample_name}.gtf"
    """
}
