process AGGREGATE_GENE_PRESENCE {
    label 'process_very_high'
    conda 'conda-forge::python=3.11 conda-forge::pandas'

    publishDir "${params.outdir}/intermediate_spreadsheets", mode: 'copy'

    input:
    path(gene_presence_files)

    output:
    path("gene_presence/*.tsv"), emit: summaries

    script:
    """
    mkdir -p gene_presence_input gene_presence
    for f in ${gene_presence_files}; do
        cp \$f gene_presence_input/
    done
    aggregate_gene_presence.py \
        --input-dir gene_presence_input \
        --output-dir gene_presence
    """

    stub:
    """
    mkdir -p gene_presence
    touch gene_presence/sankey_level1_gene_presence.tsv
    touch gene_presence/gene_presence_per_assembly_summary.tsv
    touch gene_presence/biotype_enrichment_exclusive_genes.tsv
    touch gene_presence/biotype_enrichment_summary.tsv
    """
}


process AGGREGATE_SANKEY {
    label 'process_very_high'
    conda 'conda-forge::python=3.11 conda-forge::pandas'

    publishDir "${params.outdir}/intermediate_spreadsheets", mode: 'copy'

    input:
    path(gene_presence_files)
    path(rbh_files)
    path(transcript_concordance_files)
    path(coding_integrity_files)

    output:
    path("sankey/*.tsv"), emit: summaries

    script:
    """
    mkdir -p gp_dir rbh_dir tc_dir ci_dir sankey
    for f in ${gene_presence_files}; do cp \$f gp_dir/; done
    for f in ${rbh_files}; do cp \$f rbh_dir/; done
    for f in ${transcript_concordance_files}; do cp \$f tc_dir/; done
    for f in ${coding_integrity_files}; do cp \$f ci_dir/; done

    aggregate_sankey_flows.py \
        --gene-presence-dir gp_dir \
        --rbh-dir rbh_dir \
        --transcript-concordance-dir tc_dir \
        --coding-integrity-dir ci_dir \
        --output-dir sankey
    """

    stub:
    """
    mkdir -p sankey
    touch sankey/sankey_flow_counts.tsv
    touch sankey/sankey_per_assembly_flows.tsv
    """
}


process AGGREGATE_CODING_INTEGRITY {
    label 'process_very_high'
    conda 'conda-forge::python=3.11 conda-forge::pandas'

    publishDir "${params.outdir}/intermediate_spreadsheets", mode: 'copy'

    input:
    path(coding_integrity_files)

    output:
    path("coding_integrity/*.tsv"), emit: summaries

    script:
    """
    mkdir -p ci_input coding_integrity
    for f in ${coding_integrity_files}; do cp \$f ci_input/; done
    aggregate_coding_integrity.py \
        --input-dir ci_input \
        --output-dir coding_integrity
    """

    stub:
    """
    mkdir -p coding_integrity
    touch coding_integrity/cds_concordance_matrices.tsv
    touch coding_integrity/coding_integrity_per_assembly.tsv
    touch coding_integrity/cds_classification_distribution.tsv
    """
}


process AGGREGATE_TRANSCRIPT_COUNTS {
    label 'process_very_high'
    conda 'conda-forge::python=3.11 conda-forge::pandas'

    publishDir "${params.outdir}/intermediate_spreadsheets", mode: 'copy'

    input:
    path(transcript_concordance_files)

    output:
    path("transcript_counts/*.tsv"), emit: summaries

    script:
    """
    mkdir -p tc_input transcript_counts
    for f in ${transcript_concordance_files}; do cp \$f tc_input/; done
    aggregate_transcript_counts.py \
        --input-dir tc_input \
        --output-dir transcript_counts
    """

    stub:
    """
    mkdir -p transcript_counts
    touch transcript_counts/transcript_count_scatter_data.tsv
    touch transcript_counts/transcript_count_per_assembly.tsv
    """
}


process AGGREGATE_DIVERGENCE {
    label 'process_very_high'
    conda 'conda-forge::python=3.11 conda-forge::pandas'

    publishDir "${params.outdir}/intermediate_spreadsheets", mode: 'copy'

    input:
    path(divergence_files)

    output:
    path("divergence/*.tsv"), emit: summaries

    script:
    """
    mkdir -p div_input divergence
    for f in ${divergence_files}; do cp \$f div_input/; done
    aggregate_grch38_divergence.py \
        --input-dir div_input \
        --output-dir divergence
    """

    stub:
    """
    mkdir -p divergence
    touch divergence/grch38_divergence_cross_tab.tsv
    touch divergence/grch38_divergence_per_assembly.tsv
    touch divergence/grch38_divergence_by_biotype.tsv
    """
}
