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
            ens_id_idx = header.index('ensembl_id')
            cat_id_idx = header.index('cat_id')
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


def get_introns(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Derive introns from sorted exon list.
    Returns list of (start, end) tuples for introns.
    Intron start = exon1 end + 1
    Intron end = exon2 start - 1
    """
    introns = []
    if len(exons) < 2:
        return introns
        
    for i in range(len(exons) - 1):
        intron_start = exons[i][1] + 1
        intron_end = exons[i+1][0] - 1
        introns.append((intron_start, intron_end))
    return introns


def calculate_jaccard(exons1: List[Tuple[int, int]], exons2: List[Tuple[int, int]]) -> float:
    """Calculate Jaccard Index of exon intervals at base-pair level."""
    if not exons1 or not exons2:
        return 0.0
        
    # Simple set-based approximation for now, or precise interval calculation?
    # Precise interval calculation is better.
    
    def get_coverage(exons):
        covered = set()
        for start, end in exons:
            covered.update(range(start, end + 1))
        return covered

    set1 = get_coverage(exons1)
    set2 = get_coverage(exons2)
    
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    
    return intersection / union if union > 0 else 0.0


def compare_structures(exons1: List[Tuple[int, int]],
                       exons2: List[Tuple[int, int]],
                       tolerance: int = 0) -> Dict:
    """
    Compare two transcript structures with detailed classification.
    """
    result = {
        'classification': 'No_Match',
        'jaccard': 0.0
    }
    
    if not exons1 or not exons2:
        return result

    # 1. Exact Exon Match
    # Check sets first for speed (ignoring order/dupes which shouldn't happen)
    if exons1 == exons2:
        result['classification'] = 'Exact_Match'
        result['jaccard'] = 1.0
        return result

    # Calculate Jaccard
    result['jaccard'] = calculate_jaccard(exons1, exons2)
    
    # 2. Intron Chain Analysis
    introns1 = get_introns(exons1)
    introns2 = get_introns(exons2)
    
    if not introns1 and not introns2:
        # Both single-exon. If not exact match (checked above), determine overlap
        if result['jaccard'] > 0:
            result['classification'] = 'Single_Exon_Overlap'
        return result

    if introns1 == introns2:
        # Introns match exactly, but exons differed (checked above) -> UTR diffs
        result['classification'] = 'Intron_Match'
        return result
        
    # Check for Subset/Superset
    set_i1 = set(introns1)
    set_i2 = set(introns2)
    
    if set_i1.issubset(set_i2):
        result['classification'] = 'Intron_Subset' # 1 is inside 2
        return result
    if set_i2.issubset(set_i1):
        result['classification'] = 'Intron_Superset' # 1 contains 2
        return result
        
    # Check for 5' or 3' Partial Matches (if multi-intron)
    # We require at least one intron to match to call it a structural match
    common_introns = set_i1 & set_i2
    
    if common_introns:
        # Determine location of matches
        # 5' Match: First introns match
        # 3' Match: Last introns match
        
        # Sort to check ends
        # Note: introns1/2 are already sorted by coordinate. 
        # But "5' end" depends on Strand!
        # This function doesn't know strand, assuming logic handles input order consistent with genomic coords.
        # Arguments passed to this function are genomic sorted. 
        # So low coord = 5' for forward, 3' for reverse.
        # We will classify as "Partial_5_Prime" or "Partial_3_Prime" relative to genomic coords first,
        # and caller can swap label if strand is minus? 
        # Actually, let's just stick to "Genomic_5_Prime_Match" etc or check simple list ends.
        
        match_start = (introns1[0] == introns2[0])
        match_end = (introns1[-1] == introns2[-1])
        
        if match_start and match_end:
            # Middle differs?
            result['classification'] = 'Internal_Mismatch'
        elif match_start:
            result['classification'] = 'Partial_5_Prime_Match' # Genomic 5'
        elif match_end:
            result['classification'] = 'Partial_3_Prime_Match' # Genomic 3'
        else:
            result['classification'] = 'Internal_Match' # Some internal introns match
    
    else:
        # No shared introns, but maybe overlap?
        if result['jaccard'] > 0:
            result['classification'] = 'Exon_Overlap'
            
    return result


def calculate_concordance(ensembl_txs: Dict, ensembl_exons: Dict,
                          cat_txs: Dict, cat_exons: Dict,
                          ensembl_gene: str, cat_gene: str) -> Dict:
    """Calculate transcript concordance metrics for a gene pair."""

    e_transcripts = ensembl_txs.get(ensembl_gene, {})
    c_transcripts = cat_txs.get(cat_gene, {})

    n_ensembl = len(e_transcripts)
    n_cat = len(c_transcripts)

    # Categories
    counts = {
        'Exact_Match': 0,
        'Intron_Match': 0,
        'Intron_Subset': 0,
        'Intron_Superset': 0,
        'Partial_5_Prime': 0,  # Combined strand-aware
        'Partial_3_Prime': 0,  # Combined strand-aware
        'Other_Partial': 0,
        'No_Match': 0
    }
    
    best_jaccard_sum = 0.0
    matched_ensembl = set()

    for e_tid, e_meta in e_transcripts.items():
        e_exon_list = ensembl_exons.get(e_tid, [])
        e_strand = e_meta['strand']
        
        best_class = 'No_Match'
        best_jaccard = 0.0
        
        # Priority for "Best Match": Exact > Intron > Superset > Subset > Partials
        priority = ['Exact_Match', 'Intron_Match', 'Intron_Superset', 'Intron_Subset', 
                   'Partial_5_Prime_Match', 'Partial_3_Prime_Match', 'Internal_Match', 
                   'Internal_Mismatch', 'Exon_Overlap', 'Single_Exon_Overlap', 'No_Match']
        
        def get_score(cls):
            try: return priority.index(cls) 
            except: return 99

        for c_tid, c_meta in c_transcripts.items():
            c_exon_list = cat_exons.get(c_tid, [])
            
            res = compare_structures(e_exon_list, c_exon_list)
            cls = res['classification']
            jac = res['jaccard']
            
            if get_score(cls) < get_score(best_class):
                best_class = cls
                best_jaccard = jac
            elif get_score(cls) == get_score(best_class) and jac > best_jaccard:
                best_jaccard = jac

        # Post-process strand for partials
        # compare_structures returns "Partial_5_Prime_Match" referring to Genomic 5' (low coord)
        # If strand is '-', low coord is 3' end of transcript.
        final_class = best_class
        if e_strand == '-':
            if best_class == 'Partial_5_Prime_Match':
                final_class = 'Partial_3_Prime' # Genomic 5' is Transcript 3' on minus strand
            elif best_class == 'Partial_3_Prime_Match':
                final_class = 'Partial_5_Prime' # Genomic 3' is Transcript 5' on minus strand
        else:
            if best_class == 'Partial_5_Prime_Match':
                final_class = 'Partial_5_Prime'
            elif best_class == 'Partial_3_Prime_Match':
                final_class = 'Partial_3_Prime'
        
        # Simplify other partials
        if final_class in ['Internal_Match', 'Internal_Mismatch', 'Exon_Overlap', 'Single_Exon_Overlap']:
            final_class = 'Other_Partial'
            
        if final_class in counts:
            counts[final_class] += 1
            matched_ensembl.add(e_tid)
        elif final_class == 'No_Match':
            counts['No_Match'] += 1
        else:
            # Fallback for anything else (should be mapped above)
            counts['Other_Partial'] += 1
            matched_ensembl.add(e_tid)
            
        best_jaccard_sum += best_jaccard

    # Calculate concordance rate (Exact matches only)
    concordance_rate = counts['Exact_Match'] / n_ensembl if n_ensembl > 0 else 0.0
    avg_jaccard = best_jaccard_sum / n_ensembl if n_ensembl > 0 else 0.0

    return {
        'n_ensembl_transcripts': n_ensembl,
        'n_cat_transcripts': n_cat,
        'n_exact_matches': counts['Exact_Match'],
        'n_intron_matches': counts['Intron_Match'],
        'n_subset': counts['Intron_Subset'],
        'n_superset': counts['Intron_Superset'],
        'n_partial_5': counts['Partial_5_Prime'],
        'n_partial_3': counts['Partial_3_Prime'],
        'n_other_partial': counts['Other_Partial'],
        'n_unmatched': counts['No_Match'],
        'transcript_concordance_rate': round(concordance_rate, 4),
        'avg_jaccard_index': round(avg_jaccard, 4)
    }


def main():
    args = parse_args()

    print(f"Loading Ensembl transcripts from {args.ensembl_gff}", file=sys.stderr)
    ensembl_txs, ensembl_exons = load_transcripts_and_exons(args.ensembl_gff, 'ensembl')

    print(f"Loading CAT transcripts from {args.cat_gff}", file=sys.stderr)
    cat_txs, cat_exons = load_transcripts_and_exons(args.cat_gff, 'cat')

    print(f"Loading RBH pairs from {args.rbh_pairs}", file=sys.stderr)
    rbh_pairs = load_rbh_pairs(args.rbh_pairs)

    print(f"Analyzing {len(rbh_pairs)} RBH gene pairs with enhanced logic...", file=sys.stderr)

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
            'n_intron_matches',
            'n_subset',
            'n_superset',
            'n_partial_5',
            'n_partial_3',
            'n_other_partial',
            'n_unmatched',
            'transcript_concordance_rate',
            'avg_jaccard_index'
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
                str(metrics['n_intron_matches']),
                str(metrics['n_subset']),
                str(metrics['n_superset']),
                str(metrics['n_partial_5']),
                str(metrics['n_partial_3']),
                str(metrics['n_other_partial']),
                str(metrics['n_unmatched']),
                str(metrics['transcript_concordance_rate']),
                str(metrics['avg_jaccard_index'])
            ]) + '\n')

    print(f"Wrote enhanced transcript concordance metrics to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
