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

    for VERSION in 2025_08 2024_10 2024_08 2024_06; do
        if curl -f -L -o ${assembly_accession}.gff3.gz \${BASE_URL}/${assembly_accession}/ensembl/geneset/\${VERSION}/genes.gff3.gz 2>/dev/null; then
            exit 0
        fi
    done
    """

    stub:
    """
    touch ${assembly_accession}.gff3.gz
    """
}
