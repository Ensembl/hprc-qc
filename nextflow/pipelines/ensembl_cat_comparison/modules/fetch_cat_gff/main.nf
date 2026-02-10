process FETCH_CAT_GFF {
    tag "$sample_name"
    label 'process_download'
    conda 'conda-forge::awscli'
    container 'oras://community.wave.seqera.io/library/awscli:9678a829f8ce915d'
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.cat_cache_dir}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name)

    output:
    tuple val(assembly_accession), val(sample_name), path("*.gff3.gz"), emit: gff

    script:
    """
    set +e

    # Extract base sample name by removing haplotype suffix
    BASE_SAMPLE=\$(echo "${sample_name}" | sed 's/_pat\$//' | sed 's/_mat\$//' | sed 's/_hap1\$//' | sed 's/_hap2\$//')

    # Try both HPRC and HPRC_PLUS paths
    S3_PATHS="s3://human-pangenomics/working/HPRC/\${BASE_SAMPLE}/assemblies/release2/annotation/cat s3://human-pangenomics/working/HPRC_PLUS/\${BASE_SAMPLE}/assemblies/release2/annotation/cat"

    # Get list of alternative sample names (both pat/mat and hap1/hap2)
    SAMPLE_NAMES="${sample_name}"
    if echo "${sample_name}" | grep -q "_pat\$"; then
        ALT_NAME=\$(echo "${sample_name}" | sed 's/_pat\$/_hap1/')
        SAMPLE_NAMES="${sample_name} \${ALT_NAME}"
    elif echo "${sample_name}" | grep -q "_mat\$"; then
        ALT_NAME=\$(echo "${sample_name}" | sed 's/_mat\$/_hap2/')
        SAMPLE_NAMES="${sample_name} \${ALT_NAME}"
    elif echo "${sample_name}" | grep -q "_hap1\$"; then
        ALT_NAME=\$(echo "${sample_name}" | sed 's/_hap1\$/_pat/')
        SAMPLE_NAMES="${sample_name} \${ALT_NAME}"
    elif echo "${sample_name}" | grep -q "_hap2\$"; then
        ALT_NAME=\$(echo "${sample_name}" | sed 's/_hap2\$/_mat/')
        SAMPLE_NAMES="${sample_name} \${ALT_NAME}"
    fi

    # Try different version patterns (most recent first)
    VERSIONS="v1.1.0 v1.0.1 v1.0.0"

    # Try each combination of S3 path, sample name, and version
    for S3_PATH in \${S3_PATHS}; do
        echo "Searching for CAT GFF in \${S3_PATH}" >&2
        for SAMPLE_VAR in \${SAMPLE_NAMES}; do
            for VERSION in \${VERSIONS}; do
                CAT_FILE="\${SAMPLE_VAR}_hprc_r2_\${VERSION}_cat_v1.1.gff3.gz"
                echo "Trying: \${CAT_FILE}" >&2

                if aws s3 cp \${S3_PATH}/\${CAT_FILE} ${sample_name}_cat.gff3.gz --no-sign-request 2>/dev/null; then
                    if [ -f "${sample_name}_cat.gff3.gz" ] && [ -s "${sample_name}_cat.gff3.gz" ]; then
                        echo "Successfully downloaded: \${CAT_FILE} from \${S3_PATH}" >&2
                        exit 0
                    fi
                fi
            done
        done
    done

    # If nothing worked, list available files and exit
    echo "ERROR: Failed to download CAT GFF for ${sample_name}" >&2
    echo "Tried all combinations of paths (HPRC, HPRC_PLUS), sample names (\${SAMPLE_NAMES}), and versions (\${VERSIONS})" >&2
    for S3_PATH in \${S3_PATHS}; do
        echo "Listing available files in \${S3_PATH}:" >&2
        aws s3 ls \${S3_PATH}/ --no-sign-request >&2 || echo "  (path does not exist)" >&2
    done
    exit 1
    """

    stub:
    """
    touch ${sample_name}_cat.gff3.gz
    """
}
