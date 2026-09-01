process FETCH_CAT_GFF {
    tag "$sample_name"
    label 'process_download'
    // Use a base image that includes bash and curl; avoids apt-get under Singularity
    container 'docker://rockylinux:9'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    publishDir "${params.cat_cache_dir}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), val(cat_gff_url), val(cat_gff_path)

    output:
    tuple val(assembly_accession), val(sample_name), path("*.gff3.gz", optional: true), emit: gff

    script:
    """
    set -e

    if [ -n "${cat_gff_path}" ]; then
        echo "Using local CAT GFF: ${cat_gff_path}" >&2
        case "${cat_gff_path}" in
            *.gz)
                cp "${cat_gff_path}" ${sample_name}_cat.gff3.gz
                ;;
            *)
                gzip -c "${cat_gff_path}" > ${sample_name}_cat.gff3.gz
                ;;
        esac
        test -s ${sample_name}_cat.gff3.gz
        exit 0
    fi

    if [ -n "${cat_gff_url}" ]; then
        echo "Downloading CAT GFF from explicit URL: ${cat_gff_url}" >&2
        case "${cat_gff_url}" in
            *.gz|*.gz\\?*|*.gz#*)
                curl -fsSL -o ${sample_name}_cat.gff3.gz "${cat_gff_url}"
                ;;
            *)
                curl -fsSL "${cat_gff_url}" | awk '/^##FASTA/ { exit } { print }' | gzip -c > ${sample_name}_cat.gff3.gz
                ;;
        esac
        test -s ${sample_name}_cat.gff3.gz
        exit 0
    fi

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
        CAT_FILE="\${SAMPLE_VAR}_consensus.gff3.gz"
        URL="\${BASE_URL}/\${CAT_FILE}"
        echo "Trying: \${URL}" >&2

        if curl -fsSL -o ${sample_name}_cat.gff3.gz "\${URL}" 2>/dev/null; then
            if [ -f "${sample_name}_cat.gff3.gz" ] && [ -s "${sample_name}_cat.gff3.gz" ]; then
                echo "Successfully downloaded: \${CAT_FILE}" >&2
                exit 0
            fi
        fi
    done

    # If nothing worked, log warning and exit gracefully (no output file produced)
    echo "WARNING: CAT GFF not found for ${sample_name}, skipping" >&2
    echo "Tried sample name variants: \${SAMPLE_NAMES}" >&2
    exit 0
    """

    stub:
    """
    touch ${sample_name}_cat.gff3.gz
    """
}
