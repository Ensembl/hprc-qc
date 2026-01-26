process FETCH_ASSEMBLY_REPORT {
    tag "$assembly_accession"
    label 'process_download'
    conda 'bioconda::biopython conda-forge::requests'
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.assembly_reports_dir}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name)

    output:
    tuple val(assembly_accession), val(sample_name), path("${assembly_accession}_assembly_report.txt"), emit: report

    script:
    """
    fetch_ncbi_assembly_report.py ${assembly_accession} ${assembly_accession}_assembly_report.txt
    """

    stub:
    """
    touch ${assembly_accession}_assembly_report.txt
    """
}
