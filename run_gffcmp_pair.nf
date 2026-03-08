nextflow.enable.dsl = 2

include { GFF_TO_GTF as ENS_GFF_TO_GTF } from './nextflow/pipelines/ensembl_cat_comparison/modules/gff_to_gtf/main'
include { GFF_TO_GTF as CAT_GFF_TO_GTF } from './nextflow/pipelines/ensembl_cat_comparison/modules/gff_to_gtf/main'
include { GFFCOMPARE_RUN }   from './nextflow/pipelines/ensembl_cat_comparison/modules/gffcompare_run2/main'
include { GFFCOMPARE_PARSE } from './nextflow/pipelines/ensembl_cat_comparison/modules/gffcompare_parse/main'

params.outdir = params.outdir ?: 'results_gffcmp2'

workflow {
  // Two pseudo-assemblies using local CHM13 files
  ens = Channel.of(
    tuple('TEST_A','CHM13', file('CHM13.r113.v1.ucsc.gff3.gz')),
    tuple('TEST_B','CHM13', file('CHM13.r113.v1.ucsc.gff3.gz'))
  )
  cat = Channel.of(
    tuple('TEST_A','CHM13', file('CHM13v2.0_CAT_Liftoff_v1.numeric.sorted.gff3.gz')),
    tuple('TEST_B','CHM13', file('CHM13v2.0_CAT_Liftoff_v1.numeric.sorted.gff3.gz'))
  )

  ENS_GFF_TO_GTF(ens)
  CAT_GFF_TO_GTF(cat)

  ens_gtf = ENS_GFF_TO_GTF.out.gtf.map{ acc,samp,gtf -> tuple([acc,samp], gtf) }
  cat_gtf = CAT_GFF_TO_GTF.out.gtf.map{ acc,samp,gtf -> tuple([acc,samp], gtf) }

  ens_gtf.join(cat_gtf)
    .map{ key, ens_gtf_path, cat_gtf_path -> tuple(key[0], key[1], ens_gtf_path, cat_gtf_path) }
    .flatMap{ acc, samp, ens_gtf_path, cat_gtf_path ->
      [
        tuple(acc, samp, 'Ensembl_to_CAT', ens_gtf_path, cat_gtf_path),
        tuple(acc, samp, 'CAT_to_Ensembl', cat_gtf_path, ens_gtf_path)
      ]
    }
    .set{ jobs }

  GFFCOMPARE_RUN(jobs)
  GFFCOMPARE_PARSE(GFFCOMPARE_RUN.out.tmap)
}
