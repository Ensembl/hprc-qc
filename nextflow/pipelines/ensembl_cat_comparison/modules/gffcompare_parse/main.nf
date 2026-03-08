process GFFCOMPARE_PARSE {
    tag "${assembly_accession}_${sample_name}_${direction}"
    label 'process_low'
    conda 'conda-forge::coreutils'

    publishDir "${params.outdir}/qc_metrics/${assembly_accession}", mode: 'copy'

    input:
    tuple val(assembly_accession), val(sample_name), val(direction), path(tmap)

    output:
    tuple val(assembly_accession), val(sample_name), val(direction), path("gffcompare/class_counts_${direction}.tsv"), emit: counts

    script:
    """
    set -euo pipefail
    mkdir -p gffcompare

    out="gffcompare/class_counts_${direction}.tsv"

    # If tmap missing or empty -> header-only file
    if [ ! -s "${tmap}" ]; then
        printf "assembly_accession\tsample_name\tdirection\tclass_code\tn_transcripts\tdenominator\tpct\n" > "${out}"
        exit 0
    fi

    # Determine class_code column index from header, then count occurrences
    # Compute denominator as number of data rows (excluding header)
    awk -v acc="${assembly_accession}" -v samp="${sample_name}" -v dirn="${direction}" '
        BEGIN { FS=OFS="\t" }
        NR==1 {
            for (i=1; i<=NF; i++) {
                if (\$i=="class_code") { c=i }
            }
            next
        }
        NR>1 {
            denom++
            code = (c>0 ? \$c : "")
            counts[code]++
        }
        END {
            # Print header
            print "assembly_accession","sample_name","direction","class_code","n_transcripts","denominator","pct" > "tmp_counts.tsv"
            if (denom==0) exit 0
            # Dump unsorted; will sort with coreutils sort for determinism
            for (k in counts) {
                n = counts[k]
                pct = (denom>0 ? (100.0*n/denom) : 0)
                printf("%s\t%s\t%s\t%s\t%d\t%d\t%.4f\n", acc, samp, dirn, k, n, denom, pct) >> "tmp_counts.tsv"
            }
        }
    ' "${tmap}"

    # Ensure deterministic class_code ordering
    { head -n 1 tmp_counts.tsv && tail -n +2 tmp_counts.tsv | LC_ALL=C sort -t \$'\t' -k4,4; } > "${out}"
    rm -f tmp_counts.tsv || true
    """

    stub:
    """
    mkdir -p gffcompare
    printf "assembly_accession\tsample_name\tdirection\tclass_code\tn_transcripts\tdenominator\tpct\n" > "gffcompare/class_counts_${direction}.tsv"
    """
}
