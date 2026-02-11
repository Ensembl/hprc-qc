#!/usr/bin/env python3
"""
Assess coding sequence integrity for protein-coding RBH gene pairs.

Compares start codons, stop codons, and CDS lengths to detect potential
annotation errors or frameshifts.

Usage:
    assess_coding_integrity.py \\
        --ensembl-gff <path> \\
        --cat-gff <path> \\
        --rbh-pairs <path> \\
        --output <path> \\
        --assembly-accession <accession> \\
        --sample-name <sample>
"""

import argparse
import gzip
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


def parse_args():
    parser = argparse.ArgumentParser(description="Assess coding integrity for RBH gene pairs")
    parser.add_argument("--ensembl-gff", required=True, help="Path to Ensembl GFF3 file")
    parser.add_argument("--cat-gff", required=True, help="Path to CAT GFF3 file")
    parser.add_argument("--rbh-pairs", required=True, help="Path to RBH pairs TSV (with biotype column)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--assembly-accession", required=True, help="Assembly accession")
    parser.add_argument("--sample-name", required=True, help="Sample name")
    return parser.parse_args()


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


def load_cds_features(gff_path: str) -> Dict[str, List[Dict]]:
    """
    Load CDS features grouped by gene.

    Returns:
        {gene_id: [{chrom, start, end, strand, phase, parent}, ...]}
    """
    cds_by_gene = defaultdict(list)
    gene_to_transcript = {}  # Map transcript_id to gene_id

    # First pass: map transcripts to genes
    with open_maybe_gzip(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            attrs = parse_attributes(parts[8])

            if feature_type in ['transcript', 'mRNA', 'lnc_RNA', 'ncRNA']:
                tid = attrs.get('ID', attrs.get('transcript_id'))
                parent = attrs.get('Parent', attrs.get('gene_id'))
                if tid and parent:
                    gene_to_transcript[tid] = parent

    # Second pass: collect CDS features
    with open_maybe_gzip(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]

            if feature_type == 'CDS':
                attrs = parse_attributes(parts[8])
                parent = attrs.get('Parent', attrs.get('transcript_id'))

                if parent:
                    gene_id = gene_to_transcript.get(parent, parent)

                    cds_by_gene[gene_id].append({
                        'chrom': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6],
                        'phase': parts[7],
                        'parent': parent
                    })

    return dict(cds_by_gene)


def calculate_cds_length(cds_list: List[Dict]) -> int:
    """Calculate total CDS length from list of CDS features."""
    total = 0
    for cds in cds_list:
        total += (cds['end'] - cds['start'] + 1)
    return total


def get_start_stop_positions(cds_list: List[Dict]) -> Tuple[Optional[int], Optional[int]]:
    """
    Get start and stop codon positions from CDS features.

    Returns:
        (start_pos, stop_pos) - coordinates of start and stop codons
    """
    if not cds_list:
        return None, None

    # Sort by position
    sorted_cds = sorted(cds_list, key=lambda x: x['start'])

    strand = sorted_cds[0]['strand']

    if strand == '+':
        # Start codon at beginning of first CDS
        start_pos = sorted_cds[0]['start']
        # Stop codon at end of last CDS
        stop_pos = sorted_cds[-1]['end']
    else:
        # Reverse strand: start at end, stop at beginning
        start_pos = sorted_cds[-1]['end']
        stop_pos = sorted_cds[0]['start']

    return start_pos, stop_pos


def load_rbh_pairs(rbh_path: str) -> List[Tuple[str, str, str, str]]:
    """
    Load RBH gene pairs from TSV with biotype information.

    Returns:
        [(ensembl_id, cat_id, ensembl_biotype, cat_biotype), ...]
    """
    pairs = []
    with open(rbh_path, 'r') as f:
        header = f.readline().strip().split('\t')

        # Find column indices
        try:
            ens_id_idx = header.index('ensembl_id')
            cat_id_idx = header.index('cat_id')
            ens_bio_idx = header.index('ensembl_biotype') if 'ensembl_biotype' in header else None
            cat_bio_idx = header.index('cat_biotype') if 'cat_biotype' in header else None
        except ValueError as e:
            print(f"Error: Required column not found in RBH file: {e}", file=sys.stderr)
            sys.exit(1)

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(ens_id_idx, cat_id_idx):
                ensembl_id = parts[ens_id_idx]
                cat_id = parts[cat_id_idx]
                ensembl_bio = parts[ens_bio_idx] if ens_bio_idx is not None else 'unknown'
                cat_bio = parts[cat_bio_idx] if cat_bio_idx is not None else 'unknown'
                pairs.append((ensembl_id, cat_id, ensembl_bio, cat_bio))

    return pairs


def assess_pair(ensembl_cds: List[Dict], cat_cds: List[Dict],
                tolerance: int = 3) -> Dict:
    """
    Assess coding integrity for a gene pair.

    Args:
        tolerance: Allow start/stop positions to differ by this many bp
    """
    result = {
        'has_ensembl_cds': len(ensembl_cds) > 0,
        'has_cat_cds': len(cat_cds) > 0,
        'cds_length_ensembl': 0,
        'cds_length_cat': 0,
        'start_codon_match': False,
        'stop_codon_match': False,
        'frameshift_detected': False,
        'length_difference': 0
    }

    if not ensembl_cds or not cat_cds:
        return result

    # Calculate lengths
    e_len = calculate_cds_length(ensembl_cds)
    c_len = calculate_cds_length(cat_cds)

    result['cds_length_ensembl'] = e_len
    result['cds_length_cat'] = c_len
    result['length_difference'] = abs(e_len - c_len)

    # Check for frameshift (length difference not divisible by 3)
    if result['length_difference'] > 0 and result['length_difference'] % 3 != 0:
        result['frameshift_detected'] = True

    # Get start/stop positions
    e_start, e_stop = get_start_stop_positions(ensembl_cds)
    c_start, c_stop = get_start_stop_positions(cat_cds)

    if e_start is not None and c_start is not None:
        if abs(e_start - c_start) <= tolerance:
            result['start_codon_match'] = True

    if e_stop is not None and c_stop is not None:
        if abs(e_stop - c_stop) <= tolerance:
            result['stop_codon_match'] = True

    return result


def main():
    args = parse_args()

    print(f"Loading Ensembl CDS features from {args.ensembl_gff}", file=sys.stderr)
    ensembl_cds = load_cds_features(args.ensembl_gff)

    print(f"Loading CAT CDS features from {args.cat_gff}", file=sys.stderr)
    cat_cds = load_cds_features(args.cat_gff)

    print(f"Loading RBH pairs from {args.rbh_pairs}", file=sys.stderr)
    rbh_pairs = load_rbh_pairs(args.rbh_pairs)

    print(f"Analyzing {len(rbh_pairs)} RBH gene pairs", file=sys.stderr)

    # Filter to protein-coding genes
    protein_coding_pairs = [
        (eid, cid, ebio, cbio) for eid, cid, ebio, cbio in rbh_pairs
        if 'protein_coding' in ebio.lower() or 'protein_coding' in cbio.lower()
    ]

    print(f"Found {len(protein_coding_pairs)} protein-coding gene pairs", file=sys.stderr)

    # Write output
    with open(args.output, 'w') as out:
        # Header
        out.write('\t'.join([
            'assembly_accession',
            'sample_name',
            'ensembl_gene_id',
            'cat_gene_id',
            'ensembl_biotype',
            'cat_biotype',
            'has_ensembl_cds',
            'has_cat_cds',
            'cds_length_ensembl',
            'cds_length_cat',
            'length_difference',
            'start_codon_match',
            'stop_codon_match',
            'frameshift_detected'
        ]) + '\n')

        # Process each protein-coding RBH pair
        for ensembl_id, cat_id, ensembl_bio, cat_bio in protein_coding_pairs:
            e_cds = ensembl_cds.get(ensembl_id, [])
            c_cds = cat_cds.get(cat_id, [])

            metrics = assess_pair(e_cds, c_cds, tolerance=3)

            out.write('\t'.join([
                args.assembly_accession,
                args.sample_name,
                ensembl_id,
                cat_id,
                ensembl_bio,
                cat_bio,
                str(metrics['has_ensembl_cds']),
                str(metrics['has_cat_cds']),
                str(metrics['cds_length_ensembl']),
                str(metrics['cds_length_cat']),
                str(metrics['length_difference']),
                str(metrics['start_codon_match']),
                str(metrics['stop_codon_match']),
                str(metrics['frameshift_detected'])
            ]) + '\n')

    print(f"Wrote coding integrity assessment to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
