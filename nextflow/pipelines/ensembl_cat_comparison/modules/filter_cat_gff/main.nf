process FILTER_CAT_GFF {
    tag "${assembly_accession}_${sample_name}"
    label 'process_low'
    conda 'conda-forge::python=3.11'

    publishDir "${params.outdir}/filtered_cat_gff/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), path(cat_gff)

    output:
    tuple val(assembly_accession), val(sample_name), path("*filtered_cat.gff3.gz"), emit: gff
    tuple val(assembly_accession), val(sample_name), path("*filter_audit.tsv"), emit: audit

    script:
    def shellQuote = { value ->
        def s = value.toString()
        return "'" + s.replace("'", "'\"'\"'") + "'"
    }
    def csvArgs = { value, flag ->
        value
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .collect { "${flag} ${shellQuote(it)}" }
            .join(' ')
    }
    def boolArg = { value, flag ->
        value.toString().toBoolean() ? flag : ''
    }
    def filterName = params.filter_cat_filter_name ?: 'filtered_cat'
    def excludeArgs = [
        csvArgs(params.filter_cat_exclude_biotypes ?: '', '--exclude-biotype'),
        csvArgs(params.filter_cat_exclude_sources ?: '', '--exclude-source'),
        csvArgs(params.filter_cat_exclude_transcript_modes ?: '', '--exclude-transcript-mode'),
        csvArgs(params.filter_cat_exclude_gene_name_regex ?: '', '--exclude-gene-name-regex'),
        boolArg(params.filter_cat_exclude_hyphenated_gene_names ?: false, '--exclude-hyphenated-gene-names'),
        csvArgs(params.filter_cat_exclude_attr_nonempty ?: '', '--exclude-attr-nonempty'),
        csvArgs(params.filter_cat_exclude_attr_equals ?: '', '--exclude-attr-equals'),
        csvArgs(params.filter_cat_exclude_attr_in ?: '', '--exclude-attr-in'),
    ].findAll { it }.join(' ')
    """
    filter_gff_by_biotype.py \\
        --input-gff ${cat_gff} \\
        --output-gff ${assembly_accession}.${filterName}.filtered_cat.gff3.gz \\
        --audit-tsv ${assembly_accession}.${filterName}.filter_audit.tsv \\
        --filter-name ${shellQuote(filterName)} \\
        ${excludeArgs}
    """

    stub:
    """
    touch ${assembly_accession}.${params.filter_cat_filter_name ?: 'filtered_cat'}.filtered_cat.gff3.gz
    touch ${assembly_accession}.${params.filter_cat_filter_name ?: 'filtered_cat'}.filter_audit.tsv
    """
}
