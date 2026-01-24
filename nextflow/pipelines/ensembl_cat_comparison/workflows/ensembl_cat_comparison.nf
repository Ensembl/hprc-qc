include { FETCH_FILES }      from '../subworkflows/fetch_files/main'
include { RUN_COMPARISON }   from '../modules/run_comparison/main'

workflow ENSEMBL_CAT_COMPARISON {
    take:
    assemblies_ch  // Channel of [assembly_accession, sample_name]

    main:
    // Fetch GFF files and assembly reports
    FETCH_FILES(assemblies_ch)

    // FETCH_FILES returns [accession, sample, ensembl_gff, cat_gff, assembly_report]
    // Pass directly to RUN_COMPARISON which expects the same tuple structure
    RUN_COMPARISON(FETCH_FILES.out.paired)

    emit:
    rbh = RUN_COMPARISON.out.rbh
    all_pairs = RUN_COMPARISON.out.all_pairs
    logs = RUN_COMPARISON.out.log
}
