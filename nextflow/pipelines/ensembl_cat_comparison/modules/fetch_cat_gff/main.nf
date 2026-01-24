process FETCH_CAT_GFF {
    tag "$sample_name"
    label 'process_low'
    conda 'conda-forge::awscli'

    input:
    tuple val(assembly_accession), val(sample_name)

    output:
    tuple val(assembly_accession), val(sample_name), path("*.gff3.gz"), emit: gff

    script:
    """
    BASE_SAMPLE=\$(echo "${sample_name}" | sed 's/_pat\$//' | sed 's/_mat\$//')
    S3_PATH="s3://human-pangenomics/working/HPRC/\${BASE_SAMPLE}/assemblies/release2/annotation/cat"

    CAT_FILE="${sample_name}_hprc_r2_v1.0.1_cat_v1.1.gff3.gz"

    aws s3 cp \${S3_PATH}/\${CAT_FILE} ${sample_name}_cat.gff3.gz --no-sign-request

    if [ ! -f "${sample_name}_cat.gff3.gz" ]; then
        echo "ERROR: Failed to download \${CAT_FILE}" >&2
        echo "Listing available files:" >&2
        aws s3 ls \${S3_PATH}/ --no-sign-request >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${sample_name}_cat.gff3.gz
    """
}
