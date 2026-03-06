#!/usr/bin/env python3
"""
Count transcripts per gene from an Ensembl GFF file.

Outputs a TSV with gene_id, gene_name, biotype, n_transcripts for every gene
in the GFF.  This provides the "total Ensembl transcripts" denominator needed
to compute full-denominator transcript concordance (including genes that have
no CAT counterpart).

Reuses the same TRANSCRIPT_TYPES and GFF parsing pattern as
calculate_transcript_concordance.py.

Usage:
    count_gff_transcripts.py \
        --gff <ensembl.gff3[.gz]> \
        --output <gene_transcript_counts.tsv> \
        --assembly-accession <GCA_...>
"""

import argparse
import gzip
import sys
from collections import defaultdict


# Same transcript types as calculate_transcript_concordance.py
TRANSCRIPT_TYPES = {
    'transcript', 'mRNA', 'lnc_RNA', 'ncRNA',
    'miRNA', 'snoRNA', 'snRNA', 'tRNA', 'rRNA',
    'pseudogenic_transcript', 'antisense_RNA',
    'guide_RNA', 'scRNA', 'vault_RNA', 'Y_RNA',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Count transcripts per gene from Ensembl GFF"
    )
    parser.add_argument("--gff", required=True,
                        help="Ensembl GFF3 file (plain or gzipped)")
    parser.add_argument("--output", required=True,
                        help="Output TSV")
    parser.add_argument("--assembly-accession", required=True,
                        help="Assembly accession for output column")
    return parser.parse_args()


def open_maybe_gzip(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_attributes(attr_string):
    attrs = {}
    for item in attr_string.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key] = val
    return attrs


def main():
    args = parse_args()

    # Pass 1: collect gene metadata (biotype, name)
    genes = {}  # gene_id -> {biotype, name}
    tx_counts = defaultdict(int)  # gene_id -> n_transcripts

    print(f"Parsing {args.gff}", file=sys.stderr)
    with open_maybe_gzip(args.gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            attrs = parse_attributes(parts[8])

            # Collect gene records
            if feature_type == 'gene' or feature_type.endswith('gene'):
                gene_id = attrs.get('ID', attrs.get('gene_id', ''))
                biotype = attrs.get('biotype', attrs.get('gene_biotype', 'unknown'))
                name = attrs.get('Name', attrs.get('gene_name', gene_id))
                if gene_id:
                    genes[gene_id] = {'biotype': biotype, 'name': name}

            # Count transcripts per parent gene
            elif feature_type in TRANSCRIPT_TYPES:
                parent = attrs.get('Parent', attrs.get('gene_id', ''))
                if parent:
                    tx_counts[parent] += 1

    # Write output
    with open(args.output, 'w') as out:
        out.write('\t'.join([
            'assembly_accession', 'gene_id', 'gene_name',
            'biotype', 'n_transcripts',
        ]) + '\n')

        for gene_id, meta in sorted(genes.items()):
            n_tx = tx_counts.get(gene_id, 0)
            out.write('\t'.join([
                args.assembly_accession,
                gene_id,
                meta['name'],
                meta['biotype'],
                str(n_tx),
            ]) + '\n')

    n_genes = len(genes)
    n_with_tx = sum(1 for g in genes if tx_counts.get(g, 0) > 0)
    total_tx = sum(tx_counts.values())
    print(f"  {n_genes} genes, {n_with_tx} with transcripts, "
          f"{total_tx} total transcripts", file=sys.stderr)


if __name__ == "__main__":
    main()
