process FETCH_CAT_GFF {
    tag "$sample_name"
    label 'process_download'
    conda 'conda-forge::curl'
    container 'oras://community.wave.seqera.io/library/curl:8.11.1--d090f1b1e3ed7548'
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

    BASE_URL="https://public.gi.ucsc.edu/~pnhebbar/share/hprc/catv2.0_genes"

    # Get list of sample name variants (both pat/mat and hap1/hap2)
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

    echo "Downloading CAT v2.0 GFF for ${sample_name}" >&2

    # Try each sample name variant
    for SAMPLE_VAR in \${SAMPLE_NAMES}; do
        CAT_FILE="\${SAMPLE_VAR}_CAT_v2.0.gff3.gz"
        URL="\${BASE_URL}/\${CAT_FILE}"
        echo "Trying: \${URL}" >&2

        if curl -fsSL -o ${sample_name}_cat.gff3.gz "\${URL}" 2>/dev/null; then
            if [ -f "${sample_name}_cat.gff3.gz" ] && [ -s "${sample_name}_cat.gff3.gz" ]; then
                echo "Successfully downloaded: \${CAT_FILE}" >&2
                exit 0
            fi
        fi
    done

    # If nothing worked, exit with error
    echo "ERROR: Failed to download CAT GFF for ${sample_name}" >&2
    echo "Tried sample name variants: \${SAMPLE_NAMES}" >&2
    echo "Expected URL pattern: \${BASE_URL}/{sample}_CAT_v2.0.gff3.gz" >&2
    exit 1
    """

    stub:
    """
    touch ${sample_name}_cat.gff3.gz
    """
}
