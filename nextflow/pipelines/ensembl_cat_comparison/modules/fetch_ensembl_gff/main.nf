process FETCH_ENSEMBL_GFF {
    tag "$assembly_accession"
    label 'process_download'
    // Use a base image that includes bash and curl; avoids apt-get under Singularity
    container 'docker://rockylinux:9'

    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.ensembl_cache_dir}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), val(ensembl_gff_url), val(ensembl_gff_path)

    output:
    tuple val(assembly_accession), val(sample_name), path("${assembly_accession}.gff3.gz"), emit: gff

    script:
    """
    if [ -n "${ensembl_gff_path}" ]; then
        echo "Using local Ensembl GFF: ${ensembl_gff_path}" >&2
        case "${ensembl_gff_path}" in
            *.gz)
                cp "${ensembl_gff_path}" "${assembly_accession}.gff3.gz"
                ;;
            *)
                gzip -c "${ensembl_gff_path}" > "${assembly_accession}.gff3.gz"
                ;;
        esac
        test -s "${assembly_accession}.gff3.gz"
        exit 0
    fi

    if [ -n "${ensembl_gff_url}" ]; then
        echo "Downloading Ensembl GFF from explicit URL: ${ensembl_gff_url}" >&2
        curl -f -L -o "${assembly_accession}.gff3.gz" "${ensembl_gff_url}"
        test -s "${assembly_accession}.gff3.gz"
        exit 0
    fi

    # 1. Define the base location for this assembly
    BASE_URL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/${assembly_accession}/ensembl/geneset"

    echo "Fetching available versions from \${BASE_URL}..." >&2

    # 2. Dynamically list directories.
    # Ensembl's HTTPS server returns an HTML index. We extract the links ending in /
    VERSIONS=\$(curl -sL \${BASE_URL}/ | \\
        grep -oP 'href="\\K[0-9]{4}_[0-9]{2}(?=/")' | \\
        sort -rV)

    if [ -z "\${VERSIONS}" ]; then
        echo "ERROR: Could not find any version directories at \${BASE_URL}" >&2
        exit 1
    fi

    # 3. Iterate through discovered versions until one works
    for VERSION in \${VERSIONS}; do
        URL="\${BASE_URL}/\${VERSION}/genes.gff3.gz"
        echo "Attempting to download version \${VERSION}: \${URL}" >&2

        if curl -f -L -o "${assembly_accession}.gff3.gz" "\${URL}"; then
            # Verify file is not empty
            if [ -s "${assembly_accession}.gff3.gz" ]; then
                echo "Successfully downloaded Ensembl GFF version \${VERSION}" >&2
                exit 0
            fi
        fi
    done

    echo "ERROR: All discovery attempts failed for ${assembly_accession}" >&2
    exit 1
    """

    stub:
    """
    touch ${assembly_accession}.gff3.gz
    """
}
