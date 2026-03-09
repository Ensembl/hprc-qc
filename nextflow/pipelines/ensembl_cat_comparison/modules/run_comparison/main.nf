process RUN_COMPARISON {
    tag "$assembly_accession"
    label 'process_medium'
    conda 'conda-forge::pandas bioconda::pyranges conda-forge::numpy'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'

    publishDir "${params.outdir}/results/${assembly_accession}", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(assembly_accession), val(sample_name), path(ensembl_gff), path(cat_gff)

    output:
    tuple val(assembly_accession), path("*.gene_pairs_rbh.tsv"), emit: rbh
    tuple val(assembly_accession), path("*.gene_pairs_all.tsv"), emit: all_pairs
    path("*.log"), emit: log

    script:
    """
    # Create a minimal working directory structure for the script
    mkdir -p dummy_dirs/{ensembl,cat}

    # Symlink GFF files to expected locations
    ln -s \$(readlink -f ${ensembl_gff}) dummy_dirs/ensembl/
    ln -s \$(readlink -f ${cat_gff}) dummy_dirs/cat/

    # Create minimal index files for the script
    echo "Assembly Accession,Sample Name" > assemblies.csv
    echo "${assembly_accession},${sample_name}" >> assemblies.csv

    echo "Sample Name,GFF File" > cat_index.csv
    echo "${sample_name},${cat_gff.name}" >> cat_index.csv

    # Run the comparison script (bin/ is automatically on PATH)
    # Chromosome normalization is 'none' because the Ensembl GFF has already
    # been renamed upstream (RENAME_ENSEMBL_GFF) to use GenBank accessions
    # matching the CAT GFF.
    hprc_ensembl_cat_overlap.py \\
        --ensembl-gff ${ensembl_gff} \\
        --cat-gff ${cat_gff} \\
        --ensembl-gff-dir dummy_dirs/ensembl \\
        --cat-anno-dir dummy_dirs/cat \\
        --cat-index cat_index.csv \\
        --assemblies-index assemblies.csv \\
        --output-prefix ${assembly_accession} \\
        --contig-normalization none \\
        2>&1 | tee ${assembly_accession}.log

    # Check if output files were created
    if [ ! -f "${assembly_accession}.gene_pairs_rbh.tsv" ]; then
        echo "ERROR: RBH file not created for ${assembly_accession}" >> ${assembly_accession}.log
        # Create empty placeholder to prevent pipeline failure
        echo -e "ensembl_gene_id\tcat_gene_id" > ${assembly_accession}.gene_pairs_rbh.tsv
        echo -e "ensembl_gene_id\tcat_gene_id" > ${assembly_accession}.gene_pairs_all.tsv
        exit 1
    fi
    """

    stub:
    """
    touch ${assembly_accession}.gene_pairs_rbh.tsv
    touch ${assembly_accession}.gene_pairs_all.tsv
    touch ${assembly_accession}.log
    """
}
