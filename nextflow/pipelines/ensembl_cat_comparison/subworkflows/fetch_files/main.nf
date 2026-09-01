include { FETCH_ENSEMBL_GFF }     from '../../modules/fetch_ensembl_gff/main'
include { FETCH_CAT_GFF }         from '../../modules/fetch_cat_gff/main'
include { FETCH_ASSEMBLY_REPORT } from '../../modules/fetch_assembly_report/main'

workflow FETCH_FILES {
    take:
    assemblies_ch  // Channel of [assembly_accession, sample_name, ensembl_gff_url, cat_gff_url, ensembl_gff_path, cat_gff_path]

    main:
    // Split input channel for parallel processing
    assemblies_ch
        .multiMap { accession, sample, ensembl_gff_url, cat_gff_url, ensembl_gff_path, cat_gff_path ->
            ensembl: tuple(accession, sample, ensembl_gff_url ?: '', ensembl_gff_path ?: '')
            cat: tuple(accession, sample, cat_gff_url ?: '', cat_gff_path ?: '')
            report: tuple(accession, sample)
        }
        .set { split_ch }

    // Load cached Ensembl GFFs and separate cached from missing
    split_ch.ensembl
        .map { accession, sample, ensembl_gff_url, ensembl_gff_path ->
            def cached = file("${params.ensembl_cache_dir}/${accession}.gff3.gz")
            tuple(accession, sample, ensembl_gff_url, ensembl_gff_path, cached.exists() ? cached : null)
        }
        .branch { accession, sample, ensembl_gff_url, ensembl_gff_path, cached_file ->
            cached: cached_file != null
                return tuple(accession, sample, cached_file)
            missing: true
                return tuple(accession, sample, ensembl_gff_url, ensembl_gff_path)
        }
        .set { ensembl_split }

    // Fetch missing Ensembl GFFs
    FETCH_ENSEMBL_GFF(ensembl_split.missing)

    // Combine cached and fetched Ensembl GFFs
    ensembl_gffs = ensembl_split.cached.mix(FETCH_ENSEMBL_GFF.out.gff)

    // Load cached CAT GFFs and separate cached from missing
    split_ch.cat
        .map { accession, sample, cat_gff_url, cat_gff_path ->
            def cached = file("${params.cat_cache_dir}/${sample}_cat.gff3.gz")
            tuple(accession, sample, cat_gff_url, cat_gff_path, cached.exists() ? cached : null)
        }
        .branch { accession, sample, cat_gff_url, cat_gff_path, cached_file ->
            cached: cached_file != null
                return tuple(accession, sample, cached_file)
            missing: true
                return tuple(accession, sample, cat_gff_url, cat_gff_path)
        }
        .set { cat_split }

    // Fetch missing CAT GFFs
    FETCH_CAT_GFF(cat_split.missing)

    // Combine cached and fetched CAT GFFs
    cat_gffs = cat_split.cached.mix(FETCH_CAT_GFF.out.gff)

    // Load cached assembly reports and separate cached from missing
    split_ch.report
        .map { accession, sample ->
            def cached = file("${params.assembly_reports_dir}/${accession}_assembly_report.txt")
            tuple(accession, sample, cached.exists() ? cached : null)
        }
        .branch { accession, sample, cached_file ->
            cached: cached_file != null
                return tuple(accession, sample, cached_file)
            missing: true
                return tuple(accession, sample)
        }
        .set { report_split }

    // Fetch missing assembly reports
    FETCH_ASSEMBLY_REPORT(report_split.missing)

    // Combine cached and fetched assembly reports
    assembly_reports = report_split.cached.mix(FETCH_ASSEMBLY_REPORT.out.report)

    // Join all files by assembly_accession and sample_name
    paired_files = ensembl_gffs
        .join(cat_gffs, by: [0, 1])
        .join(assembly_reports, by: [0, 1])
        .map { accession, sample, ensembl_gff, cat_gff, assembly_report ->
            tuple(accession, sample, ensembl_gff, cat_gff, assembly_report)
        }

    emit:
    paired = paired_files  // [assembly_accession, sample_name, ensembl_gff, cat_gff, assembly_report]
}
