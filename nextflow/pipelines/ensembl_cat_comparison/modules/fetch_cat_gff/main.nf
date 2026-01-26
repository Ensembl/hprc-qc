process FETCH_CAT_GFF {
    tag "$sample_name"
    label 'process_download'
    conda 'conda-forge::awscli'
    container 'oras://community.wave.seqera.io/library/awscli:9678a829f8ce915d'
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.cat_cache_dir}", mode: 'copy', enabled: params.cat_cache_dir != null

    input:
    tuple val(assembly_accession), val(sample_name)

    output:
    tuple val(assembly_accession), val(sample_name), path("*.gff3.gz"), emit: gff

    script:
    def cached_file = params.cat_cache_dir ? "${params.cat_cache_dir}/${sample_name}_cat.gff3.gz" : null
    if (cached_file && file(cached_file).exists()) {
        """
        echo "Using cached CAT GFF for ${sample_name}" >&2
        ln -s ${cached_file} ${sample_name}_cat.gff3.gz
        """
    } else {
    """
    # Extract base sample name by removing haplotype suffix
    BASE_SAMPLE=\$(echo "${sample_name}" | sed 's/_pat\$//' | sed 's/_mat\$//' | sed 's/_hap1\$//' | sed 's/_hap2\$//')
    S3_PATH="s3://human-pangenomics/working/HPRC/\${BASE_SAMPLE}/assemblies/release2/annotation/cat"

    # Try the original sample name first
    CAT_FILE="${sample_name}_hprc_r2_v1.0.1_cat_v1.1.gff3.gz"

    echo "Attempting to download: \${CAT_FILE}" >&2
    aws s3 cp \${S3_PATH}/\${CAT_FILE} ${sample_name}_cat.gff3.gz --no-sign-request 2>/dev/null

    # If that fails, try alternative haplotype naming (pat->hap1, mat->hap2)
    if [ ! -f "${sample_name}_cat.gff3.gz" ]; then
        echo "Primary file not found, trying alternative naming..." >&2

        # Convert naming convention
        ALT_SAMPLE="${sample_name}"
        if echo "${sample_name}" | grep -q "_pat\$"; then
            ALT_SAMPLE=\$(echo "${sample_name}" | sed 's/_pat\$/_hap1/')
        elif echo "${sample_name}" | grep -q "_mat\$"; then
            ALT_SAMPLE=\$(echo "${sample_name}" | sed 's/_mat\$/_hap2/')
        elif echo "${sample_name}" | grep -q "_hap1\$"; then
            ALT_SAMPLE=\$(echo "${sample_name}" | sed 's/_hap1\$/_pat/')
        elif echo "${sample_name}" | grep -q "_hap2\$"; then
            ALT_SAMPLE=\$(echo "${sample_name}" | sed 's/_hap2\$/_mat/')
        fi

        CAT_FILE="\${ALT_SAMPLE}_hprc_r2_v1.0.1_cat_v1.1.gff3.gz"
        echo "Attempting to download: \${CAT_FILE}" >&2
        aws s3 cp \${S3_PATH}/\${CAT_FILE} ${sample_name}_cat.gff3.gz --no-sign-request 2>/dev/null
    fi

    # Final check - if still not found, list available files and exit
    if [ ! -f "${sample_name}_cat.gff3.gz" ]; then
        echo "ERROR: Failed to download CAT GFF for ${sample_name}" >&2
        echo "Tried: ${sample_name}_hprc_r2_v1.0.1_cat_v1.1.gff3.gz and \${CAT_FILE}" >&2
        echo "Listing available files in \${S3_PATH}:" >&2
        aws s3 ls \${S3_PATH}/ --no-sign-request >&2
        exit 1
    fi

    echo "Successfully downloaded CAT GFF for ${sample_name}" >&2
        """
    }

    stub:
    """
    touch ${sample_name}_cat.gff3.gz
    """
}
