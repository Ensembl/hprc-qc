process FETCH_ENSEMBL_GFF {
    tag "$assembly_accession"
    label 'process_low'

    input:
    tuple val(assembly_accession), val(sample_name)

    output:
    tuple val(assembly_accession), val(sample_name), path("*.gff3.gz"), emit: gff

    script:
    """
    BASE_URL="https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens"
    ASSEMBLY_URL="\${BASE_URL}/${assembly_accession}/ensembl/geneset"

    VERSIONS=\$(curl -s \${ASSEMBLY_URL}/ | grep -oP '\\d{4}_\\d{2}' | sort -r | head -5)

    if [ -z "\${VERSIONS}" ]; then
        VERSIONS="2025_08 2024_10 2024_08 2024_06 2024_04"
    fi

    for VERSION in \${VERSIONS}; do
        URL="\${ASSEMBLY_URL}/\${VERSION}/genes.gff3.gz"
        if curl -f -L -o ${assembly_accession}.gff3.gz \${URL} 2>&1; then
            if [ -f "${assembly_accession}.gff3.gz" ] && [ -s "${assembly_accession}.gff3.gz" ]; then
                exit 0
            fi
        fi
    done
    """

    stub:
    """
    touch ${assembly_accession}.gff3.gz
    """
}
