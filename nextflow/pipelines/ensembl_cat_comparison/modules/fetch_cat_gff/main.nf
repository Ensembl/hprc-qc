process FETCH_CAT_GFF {
    tag "$sample_name"
    label 'process_low'

    input:
    tuple val(assembly_accession), val(sample_name)

    output:
    tuple val(assembly_accession), val(sample_name), path("*.gff3.gz"), emit: gff

    script:
    """
    BASE_SAMPLE=\$(echo "${sample_name}" | sed 's/_pat\$//' | sed 's/_mat\$//')
    S3_PATH="s3://human-pangenomics/working/HPRC/\${BASE_SAMPLE}/assemblies/release2/annotation/cat"

    PATTERNS=(
        "${sample_name}_hprc_r2_v1.0.1_cat_v1.1.gff3.gz"
        "${sample_name}_cat_v1.1.gff3.gz"
        "${sample_name}_cat.gff3.gz"
    )

    for pattern in "\${PATTERNS[@]}"; do
        if aws s3 cp \${S3_PATH}/\${pattern} ${sample_name}_cat.gff3.gz --no-sign-request 2>/dev/null; then
            exit 0
        fi
    done
    """

    stub:
    """
    touch ${sample_name}_cat.gff3.gz
    """
}
