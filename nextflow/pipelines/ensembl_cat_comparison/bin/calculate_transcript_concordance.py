#!/usr/bin/env python3
"""
Calculate transcript-level concordance for RBH gene pairs.

For each RBH gene pair, compares transcript structures (exon coordinates)
and reports exact matches, partial matches, and unmatched transcripts.

Usage:
    calculate_transcript_concordance.py \\
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
from typing import Dict, List, Tuple, Set


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate transcript concordance for RBH gene pairs")
    parser.add_argument("--ensembl-gff", required=True, help="Path to Ensembl GFF3 file")
    parser.add_argument("--cat-gff", required=True, help="Path to CAT GFF3 file")
    parser.add_argument("--rbh-pairs", required=True, help="Path to RBH pairs TSV")
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


def load_transcripts_and_exons(gff_path: str, source_label: str) -> Tuple[Dict, Dict]:
    """
    Load transcripts and their exon structures from GFF.

    Returns:
        transcripts: {gene_id: {transcript_id: metadata}}
        exons: {transcript_id: [(start, end), ...]} - sorted exon coordinates
    """
    transcripts_by_gene = defaultdict(dict)
    exons_by_transcript = defaultdict(list)

    with open_maybe_gzip(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            attrs = parse_attributes(parts[8])

            # Handle transcripts
            if feature_type in ['transcript', 'mRNA', 'lnc_RNA', 'ncRNA']:
                tid = attrs.get('ID', attrs.get('transcript_id'))
                parent = attrs.get('Parent', attrs.get('gene_id'))

                if tid and parent:
                    transcripts_by_gene[parent][tid] = {
                        'chrom': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6]
                    }

            # Handle exons
            elif feature_type == 'exon':
                parent = attrs.get('Parent', attrs.get('transcript_id'))
                if parent:
                    start = int(parts[3])
                    end = int(parts[4])
                    exons_by_transcript[parent].append((start, end))

    # Sort exons for each transcript
    for tid in exons_by_transcript:
        exons_by_transcript[tid].sort()

    return dict(transcripts_by_gene), dict(exons_by_transcript)


def load_rbh_pairs(rbh_path: str) -> List[Tuple[str, str]]:
    """Load RBH gene pairs from TSV."""
    pairs = []
    with open(rbh_path, 'r') as f:
        header = f.readline().strip().split('\t')

        # Find column indices
        try:
            ens_id_idx = header.index('ensembl_gene_id')
            cat_id_idx = header.index('cat_gene_id')
        except ValueError as e:
            print(f"Error: Required column not found in RBH file: {e}", file=sys.stderr)
            print(f"Available columns: {header}", file=sys.stderr)
            sys.exit(1)

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(ens_id_idx, cat_id_idx):
                ensembl_id = parts[ens_id_idx]
                cat_id = parts[cat_id_idx]
                pairs.append((ensembl_id, cat_id))
    return pairs


def compare_exon_structures(exons1: List[Tuple[int, int]],
                            exons2: List[Tuple[int, int]],
                            tolerance: int = 0) -> str:
    """
    Compare two exon structures.

    Returns:
        'exact' - all exons match exactly (within tolerance)
        'partial' - some exons match
        'none' - no matching structure
    """
    if not exons1 or not exons2:
        return 'none'

    # Convert to sets for comparison
    set1 = set(exons1)
    set2 = set(exons2)

    if tolerance == 0:
        if set1 == set2:
            return 'exact'
    else:
        # Check with tolerance
        all_match = True
        for e1 in exons1:
            found = False
            for e2 in exons2:
                if abs(e1[0] - e2[0]) <= tolerance and abs(e1[1] - e2[1]) <= tolerance:
                    found = True
                    break
            if not found:
                all_match = False
                break

        if all_match and len(exons1) == len(exons2):
            return 'exact'

    # Check for partial overlap
    overlap = len(set1 & set2)
    if overlap > 0:
        return 'partial'

    return 'none'


def calculate_concordance(ensembl_txs: Dict, ensembl_exons: Dict,
                          cat_txs: Dict, cat_exons: Dict,
                          ensembl_gene: str, cat_gene: str) -> Dict:
    """Calculate transcript concordance metrics for a gene pair."""

    e_transcripts = ensembl_txs.get(ensembl_gene, {})
    c_transcripts = cat_txs.get(cat_gene, {})

    n_ensembl = len(e_transcripts)
    n_cat = len(c_transcripts)

    # Find best matches for each Ensembl transcript
    exact_matches = 0
    partial_matches = 0
    matched_ensembl_ids = set()
    matched_cat_ids = set()

    for e_tid, e_meta in e_transcripts.items():
        e_exon_list = ensembl_exons.get(e_tid, [])

        best_match = 'none'
        best_cat_tid = None

        for c_tid, c_meta in c_transcripts.items():
            c_exon_list = cat_exons.get(c_tid, [])

            match_type = compare_exon_structures(e_exon_list, c_exon_list, tolerance=3)

            if match_type == 'exact':
                best_match = 'exact'
                best_cat_tid = c_tid
                break
            elif match_type == 'partial' and best_match != 'exact':
                best_match = 'partial'
                best_cat_tid = c_tid

        if best_match == 'exact':
            exact_matches += 1
            matched_ensembl_ids.add(e_tid)
            if best_cat_tid:
                matched_cat_ids.add(best_cat_tid)
        elif best_match == 'partial':
            partial_matches += 1
            matched_ensembl_ids.add(e_tid)
            if best_cat_tid:
                matched_cat_ids.add(best_cat_tid)

    # Calculate concordance rate
    concordance_rate = exact_matches / n_ensembl if n_ensembl > 0 else 0.0

    return {
        'n_ensembl_transcripts': n_ensembl,
        'n_cat_transcripts': n_cat,
        'n_exact_matches': exact_matches,
        'n_partial_matches': partial_matches,
        'n_unmatched_ensembl': n_ensembl - len(matched_ensembl_ids),
        'n_unmatched_cat': n_cat - len(matched_cat_ids),
        'transcript_concordance_rate': round(concordance_rate, 4)
    }


def main():
    args = parse_args()

    print(f"Loading Ensembl transcripts from {args.ensembl_gff}", file=sys.stderr)
    ensembl_txs, ensembl_exons = load_transcripts_and_exons(args.ensembl_gff, 'ensembl')

    print(f"Loading CAT transcripts from {args.cat_gff}", file=sys.stderr)
    cat_txs, cat_exons = load_transcripts_and_exons(args.cat_gff, 'cat')

    print(f"Loading RBH pairs from {args.rbh_pairs}", file=sys.stderr)
    rbh_pairs = load_rbh_pairs(args.rbh_pairs)

    print(f"Analyzing {len(rbh_pairs)} RBH gene pairs", file=sys.stderr)

    # Write output
    with open(args.output, 'w') as out:
        # Header
        out.write('\t'.join([
            'assembly_accession',
            'sample_name',
            'ensembl_gene_id',
            'cat_gene_id',
            'n_ensembl_transcripts',
            'n_cat_transcripts',
            'n_exact_matches',
            'n_partial_matches',
            'n_unmatched_ensembl',
            'n_unmatched_cat',
            'transcript_concordance_rate'
        ]) + '\n')

        # Process each RBH pair
        for ensembl_gene, cat_gene in rbh_pairs:
            metrics = calculate_concordance(
                ensembl_txs, ensembl_exons,
                cat_txs, cat_exons,
                ensembl_gene, cat_gene
            )

            out.write('\t'.join([
                args.assembly_accession,
                args.sample_name,
                ensembl_gene,
                cat_gene,
                str(metrics['n_ensembl_transcripts']),
                str(metrics['n_cat_transcripts']),
                str(metrics['n_exact_matches']),
                str(metrics['n_partial_matches']),
                str(metrics['n_unmatched_ensembl']),
                str(metrics['n_unmatched_cat']),
                str(metrics['transcript_concordance_rate'])
            ]) + '\n')

    print(f"Wrote transcript concordance metrics to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
