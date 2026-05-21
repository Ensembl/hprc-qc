#!/usr/bin/env python3
"""
Compare gene presence/absence by gene name between Ensembl and CAT.

Creates a presence/absence matrix for all named genes, useful for identifying
genes that are missing in one annotation set but present in the other.

Usage:
    compare_gene_presence.py \\
        --ensembl-gff <path> \\
        --cat-gff <path> \\
        --output <path> \\
        --assembly-accession <accession> \\
        --sample-name <sample>
"""

import argparse
import gzip
import sys
from typing import Dict, Set
from urllib.parse import unquote



def parse_args():
    parser = argparse.ArgumentParser(description="Compare gene presence/absence by name")
    parser.add_argument("--ensembl-gff", required=True, help="Path to Ensembl GFF3 file")
    parser.add_argument("--cat-gff", required=True, help="Path to CAT GFF3 file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--assembly-accession", required=True, help="Assembly accession")
    parser.add_argument("--sample-name", required=True, help="Sample name")
    parser.add_argument("--ensg-lookup", required=False, help="Path to TSV mapping transcript ID to ENSG ID")
    return parser.parse_args()


def extract_gene_name(attrs: Dict[str, str]) -> str:
    """Extract gene name, handling Ensembl's encoded description format."""
    for key in ['Name', 'gene_name', 'gene']:
        if key in attrs:
            return attrs[key]

    desc = attrs.get('description', '')
    if desc:
        decoded = unquote(desc)
        for part in decoded.split(';'):
            if part.startswith('parent_gene_display_xref='):
                return part.split('=', 1)[1]

    return None


def open_maybe_gzip(path: str):
    """Open file, handling gzip compression."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF3 attributes into dictionary."""
    attrs = {}
    for item in attr_string.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key] = val
    return attrs


def load_ensg_lookup(path: str) -> Dict[str, str]:
    """Load transcript ID -> ENSG ID mapping from TSV."""
    lookup = {}
    if not path:
        return lookup

    print(f"Loading ENSG lookup from {path}", file=sys.stderr)
    with open_maybe_gzip(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # trans_id -> ensg_id
                lookup[parts[0]] = parts[1]
    return lookup


def load_named_genes(gff_path: str, ensg_lookup: Dict[str, str] = None) -> Dict[str, Dict]:
    """
    Load genes with names from GFF.

    Returns:
        {gene_name: {id, biotype, chrom, start, end}}
    """
    named_genes = {}
    ensg_lookup = ensg_lookup or {}

    with open_maybe_gzip(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]

            # Only process gene features
            if feature_type != 'gene' and not feature_type.endswith('gene'):
                continue

            attrs = parse_attributes(parts[8])

            gene_name = extract_gene_name(attrs)

            # If no gene name, try to use lookup table with transcript ID
            if not gene_name:
                # Try getting transcript ID from ID or Parent or other attributes
                # In GFF3 'gene' feature, ID is usually the gene ID, but sometimes features act as transcripts
                # For CAT/Ensembl, if it's a 'gene' feature, ID *should* be gene ID.
                # But user request implies we might have transcript IDs here or need to look up by *some* ID.
                # "transcript stable id in col 1 and ensg in col 2"
                # If this feature is a gene, we might need to check if its ID is in the lookup?
                # Or maybe the feature is actually a transcript masked as gene?
                # Let's check ID against lookup
                feature_id = attrs.get('ID', '').replace('gene:', '').replace('transcript:', '')
                if feature_id in ensg_lookup:
                    gene_name = ensg_lookup[feature_id]

            if not gene_name:
                # Fall back to stable ID so nothing is invisible
                gene_name = attrs.get('gene_id', attrs.get('ID', '').replace('gene:', ''))

                # Check if this fallback ID is in lookup (e.g. if it was a transcript ID)
                if gene_name in ensg_lookup:
                    gene_name = ensg_lookup[gene_name]

                if not gene_name:
                    continue

            if not gene_name:
                continue  # Skip genes without names

            # Get gene ID
            gene_id = attrs.get('ID', attrs.get('gene_id', ''))

            # Get biotype
            biotype = attrs.get('biotype', attrs.get('gene_biotype', 'unknown'))

            # Store gene info (if duplicate names, keep first occurrence)
            if gene_name not in named_genes:
                named_genes[gene_name] = {
                    'id': gene_id,
                    'biotype': biotype,
                    'chrom': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4])
                }

    return named_genes


def main():
    args = parse_args()

    ensg_lookup = load_ensg_lookup(args.ensg_lookup)

    print(f"Loading Ensembl named genes from {args.ensembl_gff}", file=sys.stderr)
    ensembl_genes = load_named_genes(args.ensembl_gff, ensg_lookup)
    print(f"Found {len(ensembl_genes)} named genes in Ensembl", file=sys.stderr)

    print(f"Loading CAT named genes from {args.cat_gff}", file=sys.stderr)
    cat_genes = load_named_genes(args.cat_gff, ensg_lookup)
    print(f"Found {len(cat_genes)} named genes in CAT", file=sys.stderr)

    # Get all unique gene names
    all_names = set(ensembl_genes.keys()) | set(cat_genes.keys())
    print(f"Total unique gene names: {len(all_names)}", file=sys.stderr)

    # Classify genes
    both = set(ensembl_genes.keys()) & set(cat_genes.keys())
    ensembl_only = set(ensembl_genes.keys()) - set(cat_genes.keys())
    cat_only = set(cat_genes.keys()) - set(ensembl_genes.keys())

    print(f"Genes in both: {len(both)}", file=sys.stderr)
    print(f"Ensembl only: {len(ensembl_only)}", file=sys.stderr)
    print(f"CAT only: {len(cat_only)}", file=sys.stderr)

    # Write output
    with open(args.output, 'w') as out:
        # Header
        out.write('\t'.join([
            'assembly_accession',
            'sample_name',
            'gene_name',
            'present_in_ensembl',
            'present_in_cat',
            'ensembl_gene_id',
            'cat_gene_id',
            'ensembl_biotype',
            'cat_biotype',
            'ensembl_chrom',
            'cat_chrom'
        ]) + '\n')

        # Write all genes
        for gene_name in sorted(all_names):
            in_ensembl = gene_name in ensembl_genes
            in_cat = gene_name in cat_genes

            e_info = ensembl_genes.get(gene_name, {})
            c_info = cat_genes.get(gene_name, {})

            out.write('\t'.join([
                args.assembly_accession,
                args.sample_name,
                gene_name,
                str(in_ensembl),
                str(in_cat),
                e_info.get('id', ''),
                c_info.get('id', ''),
                e_info.get('biotype', ''),
                c_info.get('biotype', ''),
                e_info.get('chrom', ''),
                c_info.get('chrom', '')
            ]) + '\n')

    print(f"Wrote gene presence/absence to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
