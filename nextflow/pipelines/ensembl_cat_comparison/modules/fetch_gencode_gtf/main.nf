process FETCH_GENCODE_GTF {
    tag "gencode_v${gencode_version}"
    label 'process_download'
    container 'quay.io/biocontainers/curl:8.7.1--he654da7_0'

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    val(gencode_version)

    output:
    path("gencode.v${gencode_version}.annotation.gtf.gz"), emit: gtf

    script:
    """
    curl -fsSL -o gencode.v${gencode_version}.annotation.gtf.gz \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${gencode_version}/gencode.v${gencode_version}.annotation.gtf.gz"
    """

    stub:
    """
    touch gencode.v${gencode_version}.annotation.gtf.gz
    """
}
