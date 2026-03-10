#!/usr/bin/env python3
"""
Calculate transcript-level concordance for gene pairs.

For each gene pair, compares transcript structures (exon coordinates)
and reports exact matches, partial matches, and unmatched transcripts in
BOTH directions (Ensembl→CAT and CAT→Ensembl).

Accepts both RBH-only pairs (ensembl_id/cat_id columns) and all-pairs
files (ensembl_gene_id/cat_gene_id columns, with ensembl_biotype).

Usage:
    calculate_transcript_concordance.py \\
        --ensembl-gff <path> \\
        --cat-gff <path> \\
        --pairs <path> \\
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
    parser = argparse.ArgumentParser(description="Calculate transcript concordance for gene pairs")
    parser.add_argument("--ensembl-gff", required=True, help="Path to Ensembl GFF3 file")
    parser.add_argument("--cat-gff", required=True, help="Path to CAT GFF3 file")
    parser.add_argument("--pairs", dest="pairs",
                        help="Path to gene pairs TSV (ensembl_gene_id/cat_gene_id or ensembl_id/cat_id columns)")
    parser.add_argument("--rbh-pairs", dest="pairs",
                        help="Alias for --pairs (backward compatibility)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--assembly-accession", required=True, help="Assembly accession")
    parser.add_argument("--sample-name", required=True, help="Sample name")
    args = parser.parse_args()
    if args.pairs is None:
        parser.error("One of --pairs or --rbh-pairs is required")
    return args


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


# Feature types recognised as transcript-level records
TRANSCRIPT_TYPES = {
    'transcript', 'mRNA', 'lnc_RNA', 'ncRNA',
    'miRNA', 'snoRNA', 'snRNA', 'tRNA', 'rRNA',
    'pseudogenic_transcript', 'antisense_RNA',
    'guide_RNA', 'scRNA', 'vault_RNA', 'Y_RNA',
}


def load_transcripts_and_exons(gff_path: str) -> Tuple[Dict, Dict]:
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

            if feature_type in TRANSCRIPT_TYPES:
                tid = attrs.get('ID', attrs.get('transcript_id'))
                parent = attrs.get('Parent', attrs.get('gene_id'))

                if tid and parent:
                    transcripts_by_gene[parent][tid] = {
                        'chrom': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6]
                    }

            elif feature_type == 'exon':
                parent = attrs.get('Parent', attrs.get('transcript_id'))
                if parent:
                    exons_by_transcript[parent].append((int(parts[3]), int(parts[4])))

    for tid in exons_by_transcript:
        exons_by_transcript[tid].sort()

    return dict(transcripts_by_gene), dict(exons_by_transcript)


def load_pairs(pairs_path: str) -> List[Tuple[str, str, str]]:
    """
    Load gene pairs from TSV, handling two column-name conventions:
      - all-pairs format:  ensembl_gene_id / cat_gene_id / ensembl_biotype
      - rbh-pairs format:  ensembl_id      / cat_id      / (ensembl_biotype optional)

    Returns list of (ensembl_gene_id, cat_gene_id, ensembl_biotype).
    """
    pairs = []
    with open(pairs_path, 'r') as f:
        header = f.readline().strip().split('\t')

        ens_col = next((c for c in ['ensembl_gene_id', 'ensembl_id'] if c in header), None)
        cat_col = next((c for c in ['cat_gene_id', 'cat_id'] if c in header), None)

        if ens_col is None or cat_col is None:
            print(f"Error: could not find gene ID columns in {pairs_path}", file=sys.stderr)
            print(f"Available columns: {header}", file=sys.stderr)
            sys.exit(1)

        ens_idx = header.index(ens_col)
        cat_idx = header.index(cat_col)
        bio_idx = header.index('ensembl_biotype') if 'ensembl_biotype' in header else None

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(ens_idx, cat_idx):
                biotype = parts[bio_idx] if bio_idx is not None and bio_idx < len(parts) else ''
                pairs.append((parts[ens_idx], parts[cat_idx], biotype))
    return pairs


def get_introns(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Derive introns from sorted exon list."""
    introns = []
    for i in range(len(exons) - 1):
        introns.append((exons[i][1] + 1, exons[i + 1][0] - 1))
    return introns


def calculate_jaccard(exons1: List[Tuple[int, int]],
                      exons2: List[Tuple[int, int]]) -> float:
    """
    Calculate Jaccard index of two sorted exon interval lists.
    Uses interval arithmetic — O(n+m), not O(total_bp).
    """
    if not exons1 or not exons2:
        return 0.0

    def total_bp(exons):
        return sum(e - s + 1 for s, e in exons)

    def intersection_bp(a: List[Tuple[int, int]],
                        b: List[Tuple[int, int]]) -> int:
        total = 0
        j = 0
        for s1, e1 in a:
            while j < len(b) and b[j][1] < s1:
                j += 1
            k = j
            while k < len(b) and b[k][0] <= e1:
                s2, e2 = b[k]
                total += max(0, min(e1, e2) - max(s1, s2) + 1)
                k += 1
        return total

    inter = intersection_bp(exons1, exons2)
    union = total_bp(exons1) + total_bp(exons2) - inter
    return inter / union if union > 0 else 0.0


def compare_structures(exons1: List[Tuple[int, int]],
                       exons2: List[Tuple[int, int]]) -> Dict:
    """Compare two transcript structures with detailed classification."""
    result = {'classification': 'No_Match', 'jaccard': 0.0}

    if not exons1 or not exons2:
        return result

    if exons1 == exons2:
        result['classification'] = 'Exact_Match'
        result['jaccard'] = 1.0
        return result

    result['jaccard'] = calculate_jaccard(exons1, exons2)

    introns1 = get_introns(exons1)
    introns2 = get_introns(exons2)

    if not introns1 and not introns2:
        if result['jaccard'] > 0:
            result['classification'] = 'Single_Exon_Overlap'
        return result

    if introns1 == introns2:
        result['classification'] = 'Intron_Match'
        return result

    set_i1 = set(introns1)
    set_i2 = set(introns2)

    if set_i1.issubset(set_i2):
        result['classification'] = 'Intron_Subset'
        return result
    if set_i2.issubset(set_i1):
        result['classification'] = 'Intron_Superset'
        return result

    common_introns = set_i1 & set_i2
    if common_introns:
        match_start = (introns1[0] == introns2[0])
        match_end = (introns1[-1] == introns2[-1])
        if match_start and match_end:
            result['classification'] = 'Internal_Mismatch'
        elif match_start:
            result['classification'] = 'Partial_5_Prime_Match'
        elif match_end:
            result['classification'] = 'Partial_3_Prime_Match'
        else:
            result['classification'] = 'Internal_Match'
    else:
        if result['jaccard'] > 0:
            result['classification'] = 'Exon_Overlap'

    return result


_PRIORITY = [
    'Exact_Match', 'Intron_Match', 'Intron_Superset', 'Intron_Subset',
    'Partial_5_Prime_Match', 'Partial_3_Prime_Match', 'Internal_Match',
    'Internal_Mismatch', 'Exon_Overlap', 'Single_Exon_Overlap', 'No_Match',
]

def _score(cls: str) -> int:
    try:
        return _PRIORITY.index(cls)
    except ValueError:
        return 99


def _best_match_for_query(query_txs: Dict, query_exons: Dict,
                           target_txs: Dict, target_exons: Dict,
                           query_gene: str, target_gene: str) -> Dict:
    """
    For each transcript in query_gene, find the best-matching transcript in
    target_gene and accumulate counts.  Direction-agnostic.
    """
    q_transcripts = query_txs.get(query_gene, {})
    t_transcripts = target_txs.get(target_gene, {})

    n_query = len(q_transcripts)
    n_target = len(t_transcripts)

    counts = {
        'Exact_Match': 0,
        'Intron_Match': 0,
        'Intron_Subset': 0,
        'Intron_Superset': 0,
        'Partial_5_Prime': 0,
        'Partial_3_Prime': 0,
        'Other_Partial': 0,
        'No_Match': 0,
    }
    jaccard_sum = 0.0

    for q_tid, q_meta in q_transcripts.items():
        q_exon_list = query_exons.get(q_tid, [])
        q_strand = q_meta['strand']

        best_class = 'No_Match'
        best_jaccard = 0.0

        for t_tid in t_transcripts:
            t_exon_list = target_exons.get(t_tid, [])
            res = compare_structures(q_exon_list, t_exon_list)
            cls = res['classification']
            jac = res['jaccard']

            if _score(cls) < _score(best_class):
                best_class = cls
                best_jaccard = jac
            elif _score(cls) == _score(best_class) and jac > best_jaccard:
                best_jaccard = jac

        # Strand-aware partial relabelling
        final_class = best_class
        if best_class == 'Partial_5_Prime_Match':
            final_class = 'Partial_3_Prime' if q_strand == '-' else 'Partial_5_Prime'
        elif best_class == 'Partial_3_Prime_Match':
            final_class = 'Partial_5_Prime' if q_strand == '-' else 'Partial_3_Prime'
        elif best_class in {'Internal_Match', 'Internal_Mismatch',
                            'Exon_Overlap', 'Single_Exon_Overlap'}:
            final_class = 'Other_Partial'

        counts[final_class if final_class in counts else 'Other_Partial'] += 1
        jaccard_sum += best_jaccard

    concordance_rate = counts['Exact_Match'] / n_query if n_query > 0 else 0.0
    avg_jaccard = jaccard_sum / n_query if n_query > 0 else 0.0

    return {
        'n_query': n_query,
        'n_target': n_target,
        'n_exact': counts['Exact_Match'],
        'n_intron_match': counts['Intron_Match'],
        'n_subset': counts['Intron_Subset'],
        'n_superset': counts['Intron_Superset'],
        'n_partial_5': counts['Partial_5_Prime'],
        'n_partial_3': counts['Partial_3_Prime'],
        'n_other_partial': counts['Other_Partial'],
        'n_unmatched': counts['No_Match'],
        'concordance_rate': round(concordance_rate, 4),
        'avg_jaccard': round(avg_jaccard, 4),
    }


def main():
    args = parse_args()

    print(f"Loading Ensembl transcripts from {args.ensembl_gff}", file=sys.stderr)
    ensembl_txs, ensembl_exons = load_transcripts_and_exons(args.ensembl_gff)

    print(f"Loading CAT transcripts from {args.cat_gff}", file=sys.stderr)
    cat_txs, cat_exons = load_transcripts_and_exons(args.cat_gff)

    print(f"Loading gene pairs from {args.pairs}", file=sys.stderr)
    pairs = load_pairs(args.pairs)

    print(f"Analysing {len(pairs)} gene pairs (bidirectional)...", file=sys.stderr)

    with open(args.output, 'w') as out:
        out.write('\t'.join([
            'assembly_accession', 'sample_name',
            'ensembl_gene_id', 'cat_gene_id', 'ensembl_biotype',
            # Ensembl → CAT
            'n_ensembl_transcripts', 'n_cat_transcripts',
            'n_ens_exact', 'n_ens_intron_match',
            'n_ens_subset', 'n_ens_superset',
            'n_ens_partial_5', 'n_ens_partial_3',
            'n_ens_other_partial', 'n_ens_unmatched',
            'ens_to_cat_concordance_rate', 'ens_to_cat_avg_jaccard',
            # CAT → Ensembl
            'n_cat_exact', 'n_cat_intron_match',
            'n_cat_subset', 'n_cat_superset',
            'n_cat_partial_5', 'n_cat_partial_3',
            'n_cat_other_partial', 'n_cat_unmatched',
            'cat_to_ens_concordance_rate', 'cat_to_ens_avg_jaccard',
            # Symmetric summary
            'transcript_concordance_rate',  # mean of both directions
            'avg_jaccard_index',            # mean of both directions
        ]) + '\n')

        for ensembl_gene, cat_gene, biotype in pairs:
            fwd = _best_match_for_query(
                ensembl_txs, ensembl_exons,
                cat_txs, cat_exons,
                ensembl_gene, cat_gene
            )
            rev = _best_match_for_query(
                cat_txs, cat_exons,
                ensembl_txs, ensembl_exons,
                cat_gene, ensembl_gene
            )

            sym_concordance = round(
                (fwd['concordance_rate'] + rev['concordance_rate']) / 2, 4
            )
            sym_jaccard = round(
                (fwd['avg_jaccard'] + rev['avg_jaccard']) / 2, 4
            )

            out.write('\t'.join([
                args.assembly_accession, args.sample_name,
                ensembl_gene, cat_gene, biotype,
                str(fwd['n_query']), str(fwd['n_target']),
                str(fwd['n_exact']), str(fwd['n_intron_match']),
                str(fwd['n_subset']), str(fwd['n_superset']),
                str(fwd['n_partial_5']), str(fwd['n_partial_3']),
                str(fwd['n_other_partial']), str(fwd['n_unmatched']),
                str(fwd['concordance_rate']), str(fwd['avg_jaccard']),
                str(rev['n_exact']), str(rev['n_intron_match']),
                str(rev['n_subset']), str(rev['n_superset']),
                str(rev['n_partial_5']), str(rev['n_partial_3']),
                str(rev['n_other_partial']), str(rev['n_unmatched']),
                str(rev['concordance_rate']), str(rev['avg_jaccard']),
                str(sym_concordance), str(sym_jaccard),
            ]) + '\n')

    print(f"Wrote bidirectional transcript concordance to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
