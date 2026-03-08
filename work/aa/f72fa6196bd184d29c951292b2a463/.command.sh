#!/bin/bash -ue
mkdir -p gffcompare
    printf "assembly_accession	sample_name	direction	class_code	n_transcripts	denominator	pct
" > "gffcompare/class_counts_Ensembl_to_CAT.tsv"
