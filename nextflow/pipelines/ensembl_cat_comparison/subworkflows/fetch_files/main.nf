include { FETCH_ENSEMBL_GFF }     from '../../modules/fetch_ensembl_gff/main'
include { FETCH_CAT_GFF }         from '../../modules/fetch_cat_gff/main'
include { FETCH_ASSEMBLY_REPORT } from '../../modules/fetch_assembly_report/main'

workflow FETCH_FILES {
    take:
    assemblies_ch  // Channel of [assembly_accession, sample_name]

    main:
    // Fetch all files in parallel
    FETCH_ENSEMBL_GFF(assemblies_ch)
    FETCH_CAT_GFF(assemblies_ch)
    FETCH_ASSEMBLY_REPORT(assemblies_ch)

    // Join all files by assembly_accession and sample_name
    paired_files = FETCH_ENSEMBL_GFF.out.gff
        .join(FETCH_CAT_GFF.out.gff, by: [0, 1])
        .join(FETCH_ASSEMBLY_REPORT.out.report, by: [0, 1])
        .map { accession, sample, ensembl_gff, cat_gff, assembly_report ->
            tuple(accession, sample, ensembl_gff, cat_gff, assembly_report)
        }

    emit:
    paired = paired_files  // [assembly_accession, sample_name, ensembl_gff, cat_gff, assembly_report]
}
