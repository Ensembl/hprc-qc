process FETCH_GENCODE_GTF {
    tag "gencode_v${gencode_version}"
    label 'process_download'
    // Use a base image that includes bash and curl; avoids apt-get under Singularity
    container 'docker://rockylinux:9'

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
