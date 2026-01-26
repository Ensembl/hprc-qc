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
    BASE_SAMPLE=\$(echo "${sample_name}" | sed 's/_pat\$//' | sed 's/_mat\$//')
    S3_PATH="s3://human-pangenomics/working/HPRC/\${BASE_SAMPLE}/assemblies/release2/annotation/cat"

    # Try original sample name first (pat/mat naming)
    CAT_FILE="${sample_name}_hprc_r2_v1.0.1_cat_v1.1.gff3.gz"
    aws s3 cp \${S3_PATH}/\${CAT_FILE} ${sample_name}_cat.gff3.gz --no-sign-request 2>/dev/null

    # If that fails, try hap1/hap2 naming
    if [ ! -f "${sample_name}_cat.gff3.gz" ]; then
        echo "pat/mat naming failed, trying hap1/hap2 naming..." >&2

        # Determine haplotype number based on suffix
        if echo "${sample_name}" | grep -q "_pat\$"; then
            HAP_FILE="\${BASE_SAMPLE}_hap1_hprc_r2_v1.0.1_cat_v1.1.gff3.gz"
        elif echo "${sample_name}" | grep -q "_mat\$"; then
            HAP_FILE="\${BASE_SAMPLE}_hap2_hprc_r2_v1.0.1_cat_v1.1.gff3.gz"
        else
            echo "ERROR: Sample name doesn't match expected pattern: ${sample_name}" >&2
            exit 1
        fi

        aws s3 cp \${S3_PATH}/\${HAP_FILE} ${sample_name}_cat.gff3.gz --no-sign-request
    fi

    # Final check
    if [ ! -f "${sample_name}_cat.gff3.gz" ]; then
        echo "ERROR: Failed to download CAT GFF for ${sample_name}" >&2
        echo "Tried: \${CAT_FILE} and \${HAP_FILE}" >&2
        echo "Listing available files:" >&2
        aws s3 ls \${S3_PATH}/ --no-sign-request >&2
        exit 1
    fi
        """
    }

    stub:
    """
    touch ${sample_name}_cat.gff3.gz
    """
}
