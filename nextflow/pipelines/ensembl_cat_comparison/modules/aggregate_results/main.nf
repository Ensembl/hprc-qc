process AGGREGATE_GENE_PRESENCE {
    label 'process_very_high'
    conda 'conda-forge::python=3.11 conda-forge::pandas'

    publishDir "${params.outdir}/intermediate_spreadsheets", mode: 'copy'

    input:
    path(gene_presence_files, stageAs: 'gp_input/*')

    output:
    path("gene_presence/*.tsv"), emit: summaries

    script:
    """
    mkdir -p gene_presence
    echo "DEBUG: Files in gp_input:" >&2
    ls gp_input/ >&2
    aggregate_gene_presence.py \
        --input-dir gp_input \
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
    path(gene_presence_files, stageAs: 'gp_staged/*')
    path(rbh_files, stageAs: 'rbh_staged/*')
    path(transcript_concordance_files, stageAs: 'tc_staged/*')
    path(coding_integrity_files, stageAs: 'ci_staged/*')

    output:
    path("sankey/*.tsv"), emit: summaries

    script:
    """
    mkdir -p sankey
    echo "DEBUG: gene presence files in gp_staged:" >&2
    ls gp_staged/ >&2 || echo "gp_staged dir empty/missing" >&2
    echo "DEBUG: RBH files in rbh_staged:" >&2
    ls rbh_staged/ >&2 || echo "rbh_staged dir empty/missing" >&2
    echo "DEBUG: transcript concordance files in tc_staged:" >&2
    ls tc_staged/ >&2 || echo "tc_staged dir empty/missing" >&2
    echo "DEBUG: coding integrity files in ci_staged:" >&2
    ls ci_staged/ >&2 || echo "ci_staged dir empty/missing" >&2

    aggregate_sankey_flows.py \
        --gene-presence-dir gp_staged \
        --rbh-dir rbh_staged \
        --transcript-concordance-dir tc_staged \
        --coding-integrity-dir ci_staged \
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
    path(coding_integrity_files, stageAs: 'ci_input/*')

    output:
    path("coding_integrity/*.tsv"), emit: summaries

    script:
    """
    mkdir -p coding_integrity
    echo "DEBUG: Files in ci_input:" >&2
    ls ci_input/ >&2
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
    path(transcript_concordance_files, stageAs: 'tc_input/*')

    output:
    path("transcript_counts/*.tsv"), emit: summaries

    script:
    """
    mkdir -p transcript_counts
    echo "DEBUG: Files in tc_input:" >&2
    ls tc_input/ >&2
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
    path(divergence_files, stageAs: 'div_input/*')

    output:
    path("divergence/*.tsv"), emit: summaries

    script:
    """
    mkdir -p divergence
    echo "DEBUG: Files in div_input:" >&2
    ls div_input/ >&2
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
    touch divergence/grch38_divergence_delta_distributions.tsv
    touch divergence/grch38_divergence_per_assembly_delta_summary.tsv
    touch divergence/grch38_cds_change_type_cross_tab.tsv
    touch divergence/grch38_cds_change_type_by_biotype.tsv
    """
}


process AGGREGATE_CONCORDANCE_VS_REF {
    label 'process_very_high'
    conda 'conda-forge::python=3.11 conda-forge::pandas'

    publishDir "${params.outdir}/intermediate_spreadsheets", mode: 'copy'

    input:
    path(rbh_files, stageAs: 'rbh_staged/*')
    path(transcript_concordance_files, stageAs: 'tc_staged/*')
    path(coding_integrity_files, stageAs: 'ci_staged/*')
    path(divergence_files, stageAs: 'div_staged/*')
    path(multi_mapping_files, stageAs: 'mm_staged/*')

    output:
    path("sankey_plus_divergence/*.tsv"), emit: summaries

    script:
    """
    mkdir -p sankey_plus_divergence
    aggregate_concordance_vs_grch38.py \
        --rbh-dir rbh_staged \
        --transcript-concordance-dir tc_staged \
        --coding-integrity-dir ci_staged \
        --divergence-dir div_staged \
        --multi-mapping-dir mm_staged \
        --rbh-coverage-threshold 0.95 \
        --output-dir sankey_plus_divergence
    """

    stub:
    """
    mkdir -p sankey_plus_divergence
    touch sankey_plus_divergence/links_per_assembly.tsv
    touch sankey_plus_divergence/flows_per_assembly.tsv
    touch sankey_plus_divergence/flows_medians.tsv
    """
}


process AGGREGATE_INTRON_CHAIN_BY_BIOTYPE {
    label 'process_very_high'
    conda 'conda-forge::python=3.11 conda-forge::pandas conda-forge::numpy'

    publishDir "${params.outdir}/intermediate_spreadsheets", mode: 'copy'

    input:
    path(transcript_concordance_files, stageAs: 'tc_input/*')
    path(gene_transcript_count_files, stageAs: 'gc_input/*')

    output:
    path("intron_chain/*.tsv"), emit: summaries

    script:
    """
    mkdir -p intron_chain
    aggregate_intron_chain_by_biotype.py \\
        --transcript-concordance-dir tc_input \\
        --output-dir intron_chain \\
        --all-ensembl-genes-dir gc_input
    """

    stub:
    """
    mkdir -p intron_chain
    touch intron_chain/intron_chain_by_biotype_per_assembly.tsv
    touch intron_chain/intron_chain_full_denom_per_assembly.tsv
    touch intron_chain/jaccard_by_biotype_per_assembly.tsv
    """
}
