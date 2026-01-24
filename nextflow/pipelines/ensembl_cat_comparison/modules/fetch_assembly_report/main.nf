process FETCH_ASSEMBLY_REPORT {
    tag "$assembly_accession"
    label 'process_download'
    conda 'bioconda::biopython conda-forge::requests'
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.assembly_reports_dir}", mode: 'copy', enabled: params.assembly_reports_dir != null

    input:
    tuple val(assembly_accession), val(sample_name)

    output:
    tuple val(assembly_accession), val(sample_name), path("${assembly_accession}_assembly_report.txt"), emit: report

    script:
    def cached_file = params.assembly_reports_dir ? "${params.assembly_reports_dir}/${assembly_accession}_assembly_report.txt" : null
    if (cached_file && file(cached_file).exists()) {
        """
        echo "Using cached assembly report for ${assembly_accession}" >&2
        ln -s ${cached_file} ${assembly_accession}_assembly_report.txt
        """
    } else {
        """
        fetch_ncbi_assembly_report.py ${assembly_accession} ${assembly_accession}_assembly_report.txt
        """
    }

    stub:
    """
    touch ${assembly_accession}_assembly_report.txt
    """
}
