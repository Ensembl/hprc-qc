nextflow.enable.dsl=2

include { GFF_TO_GTF as ENS_GFF_TO_GTF } from './nextflow/pipelines/ensembl_cat_comparison/modules/gff_to_gtf/main'
include { GFF_TO_GTF as CAT_GFF_TO_GTF } from './nextflow/pipelines/ensembl_cat_comparison/modules/gff_to_gtf/main'
include { GFFCOMPARE_RUN }   from './nextflow/pipelines/ensembl_cat_comparison/modules/gffcompare_run/main'
include { GFFCOMPARE_PARSE } from './nextflow/pipelines/ensembl_cat_comparison/modules/gffcompare_parse/main'

workflow {
    Channel.of(
        tuple('GCA_TEST', 'CHM13', file('CHM13.r113.v1.gff3'))
    ).set { ens_gff }

    Channel.of(
        tuple('GCA_TEST', 'CHM13', file('CHM13v2.0_CAT_Liftoff_v1.numeric.gff3.gz'))
    ).set { cat_gff }

    ENS_GFF_TO_GTF(ens_gff)
    CAT_GFF_TO_GTF(cat_gff)

    ens = ENS_GFF_TO_GTF.out.gtf.map { acc, sample, gtf -> tuple([acc,sample], gtf) }
    cat = CAT_GFF_TO_GTF.out.gtf.map { acc, sample, gtf -> tuple([acc,sample], gtf) }

    ens.join(cat).map { key, ens_gtf, cat_gtf -> tuple(key[0], key[1], ens_gtf, cat_gtf) }.set { pairs }

    jobs = pairs.flatMap { acc, sample, ens_gtf, cat_gtf -> [
        tuple(acc, sample, 'Ensembl_to_CAT', ens_gtf, cat_gtf),
        tuple(acc, sample, 'CAT_to_Ensembl', cat_gtf, ens_gtf)
    ] }

    GFFCOMPARE_RUN(jobs)
    GFFCOMPARE_PARSE(GFFCOMPARE_RUN.out.tmap)
}
