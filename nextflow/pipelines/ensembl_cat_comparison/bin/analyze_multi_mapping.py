#!/usr/bin/env python3
"""
Analyze multi-mapping patterns in gene pair overlaps.

Identifies 1-to-many and many-to-1 relationships between Ensembl and CAT genes,
which can indicate paralogs, segmental duplications, or annotation issues.

Usage:
    analyze_multi_mapping.py \\
        --all-pairs <path> \\
        --output <path> \\
        --assembly-accession <accession> \\
        --sample-name <sample>
"""

import argparse
import sys
from collections import defaultdict
from typing import Dict, List, Tuple


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze multi-mapping gene pairs")
    parser.add_argument("--all-pairs", required=True, help="Path to all gene pairs TSV")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--assembly-accession", required=True, help="Assembly accession")
    parser.add_argument("--sample-name", required=True, help="Sample name")
    parser.add_argument("--min-overlap", type=float, default=0.1,
                       help="Minimum overlap fraction to consider (default: 0.1)")
    return parser.parse_args()


def load_gene_pairs(pairs_path: str, min_overlap: float) -> List[Dict]:
    """
    Load gene pair overlaps from TSV.

    Returns:
        List of dicts with overlap information
    """
    pairs = []

    with open(pairs_path, 'r') as f:
        header = f.readline().strip().split('\t')

        # Find column indices
        col_idx = {}
        for i, col in enumerate(header):
            col_idx[col] = i

        # Required columns
        required = ['ensembl_gene_id', 'cat_gene_id', 'frac_ensembl_covered', 'frac_cat_covered']
        for col in required:
            if col not in col_idx:
                print(f"Error: Required column '{col}' not found in input", file=sys.stderr)
                sys.exit(1)

        # Optional columns
        optional = ['ensembl_biotype', 'cat_biotype', 'overlap_bp', 'classification', 'is_rbh']

        for line in f:
            parts = line.strip().split('\t')

            e_frac = float(parts[col_idx['frac_ensembl_covered']])
            c_frac = float(parts[col_idx['frac_cat_covered']])

            # Filter by minimum overlap
            if e_frac < min_overlap and c_frac < min_overlap:
                continue

            pair_info = {
                'ensembl_id': parts[col_idx['ensembl_gene_id']],
                'cat_id': parts[col_idx['cat_gene_id']],
                'frac_ensembl': e_frac,
                'frac_cat': c_frac,
                'max_frac': max(e_frac, c_frac)
            }

            # Add optional fields
            for col in optional:
                if col in col_idx and col_idx[col] < len(parts):
                    pair_info[col] = parts[col_idx[col]]

            pairs.append(pair_info)

    return pairs


def analyze_multi_mapping(pairs: List[Dict]) -> Tuple[Dict, Dict]:
    """
    Analyze multi-mapping patterns.

    Returns:
        (ensembl_multi, cat_multi) - dicts mapping gene IDs to their matches
    """
    # Group by gene IDs
    ensembl_matches = defaultdict(list)
    cat_matches = defaultdict(list)

    for pair in pairs:
        ensembl_matches[pair['ensembl_id']].append(pair)
        cat_matches[pair['cat_id']].append(pair)

    # Filter to multi-mapping only
    ensembl_multi = {eid: matches for eid, matches in ensembl_matches.items() if len(matches) > 1}
    cat_multi = {cid: matches for cid, matches in cat_matches.items() if len(matches) > 1}

    return ensembl_multi, cat_multi


def classify_match_type(matches: List[Dict], source_id_key: str) -> str:
    """
    Classify the type of multi-mapping.

    Returns:
        'tandem' - all on same chromosome
        'dispersed' - on different chromosomes
        'mixed' - combination
    """
    # Get chromosomes (if available)
    chroms = set()
    for match in matches:
        # This would require chromosome info in the pairs file
        # For now, just return 'unknown'
        pass

    return 'unknown'


def main():
    args = parse_args()

    print(f"Loading gene pairs from {args.all_pairs}", file=sys.stderr)
    pairs = load_gene_pairs(args.all_pairs, args.min_overlap)
    print(f"Loaded {len(pairs)} gene pairs with overlap >= {args.min_overlap}", file=sys.stderr)

    print("Analyzing multi-mapping patterns...", file=sys.stderr)
    ensembl_multi, cat_multi = analyze_multi_mapping(pairs)

    print(f"Ensembl genes with multiple CAT matches: {len(ensembl_multi)}", file=sys.stderr)
    print(f"CAT genes with multiple Ensembl matches: {len(cat_multi)}", file=sys.stderr)

    # Write output
    with open(args.output, 'w') as out:
        # Header
        out.write('\t'.join([
            'assembly_accession',
            'sample_name',
            'gene_id',
            'source',
            'n_matches',
            'match_type',
            'matched_gene_ids',
            'matched_biotypes',
            'overlap_fractions',
            'max_overlap_fraction',
            'has_rbh_match'
        ]) + '\n')

        # Write Ensembl multi-mapping genes
        for ens_id, matches in sorted(ensembl_multi.items()):
            cat_ids = [m['cat_id'] for m in matches]
            cat_biotypes = [m.get('cat_biotype', '') for m in matches]
            fracs = [f"{m['frac_ensembl']:.3f}" for m in matches]
            max_frac = max(m['frac_ensembl'] for m in matches)

            # Check if any match is RBH
            has_rbh = any(m.get('is_rbh', 'False') == 'True' for m in matches)

            match_type = '1-to-many'

            out.write('\t'.join([
                args.assembly_accession,
                args.sample_name,
                ens_id,
                'ensembl',
                str(len(matches)),
                match_type,
                ';'.join(cat_ids),
                ';'.join(cat_biotypes),
                ';'.join(fracs),
                f"{max_frac:.3f}",
                str(has_rbh)
            ]) + '\n')

        # Write CAT multi-mapping genes
        for cat_id, matches in sorted(cat_multi.items()):
            ens_ids = [m['ensembl_id'] for m in matches]
            ens_biotypes = [m.get('ensembl_biotype', '') for m in matches]
            fracs = [f"{m['frac_cat']:.3f}" for m in matches]
            max_frac = max(m['frac_cat'] for m in matches)

            # Check if any match is RBH
            has_rbh = any(m.get('is_rbh', 'False') == 'True' for m in matches)

            match_type = 'many-to-1'

            out.write('\t'.join([
                args.assembly_accession,
                args.sample_name,
                cat_id,
                'cat',
                str(len(matches)),
                match_type,
                ';'.join(ens_ids),
                ';'.join(ens_biotypes),
                ';'.join(fracs),
                f"{max_frac:.3f}",
                str(has_rbh)
            ]) + '\n')

    print(f"Wrote multi-mapping analysis to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
