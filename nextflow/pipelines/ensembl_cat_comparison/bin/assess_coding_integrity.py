#!/usr/bin/env python3
"""
Assess coding sequence integrity for protein-coding RBH gene pairs.

For each RBH pair, identifies the best-matching Ensembl/CAT transcript pair
(maximising CDS exon overlap) then compares start codons, stop codons, and
CDS lengths on that pair alone — avoiding inflation from multi-isoform genes.

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
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_attributes(attr_string: str) -> Dict[str, str]:
    attrs = {}
    for item in attr_string.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key] = val
    return attrs


TRANSCRIPT_TYPES = {
    'transcript', 'mRNA', 'lnc_RNA', 'ncRNA',
    'miRNA', 'snoRNA', 'snRNA', 'tRNA', 'rRNA',
    'pseudogenic_transcript',
}


def load_cds_by_transcript(gff_path: str) -> Tuple[Dict[str, List[str]], Dict[str, List[Dict]]]:
    """
    Load CDS features grouped by transcript, and gene → transcript ID mapping.

    Returns:
        gene_to_transcripts: {gene_id: [transcript_id, ...]}
        cds_by_transcript:   {transcript_id: [{chrom, start, end, strand, phase}, ...]}
    """
    gene_to_transcripts: Dict[str, List[str]] = defaultdict(list)
    cds_by_transcript: Dict[str, List[Dict]] = defaultdict(list)
    transcript_to_gene: Dict[str, str] = {}

    # Pass 1: build transcript → gene map
    with open_maybe_gzip(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] not in TRANSCRIPT_TYPES:
                continue
            attrs = parse_attributes(parts[8])
            tid = attrs.get('ID', attrs.get('transcript_id'))
            parent = attrs.get('Parent', attrs.get('gene_id'))
            if tid and parent:
                transcript_to_gene[tid] = parent
                gene_to_transcripts[parent].append(tid)

    # Pass 2: collect CDS per transcript
    with open_maybe_gzip(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attrs = parse_attributes(parts[8])
            parent = attrs.get('Parent', attrs.get('transcript_id'))
            if parent and parent in transcript_to_gene:
                cds_by_transcript[parent].append({
                    'chrom': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'phase': parts[7],
                })

    for tid in cds_by_transcript:
        cds_by_transcript[tid].sort(key=lambda x: x['start'])

    return dict(gene_to_transcripts), dict(cds_by_transcript)


def cds_intersection_bp(cds_a: List[Dict], cds_b: List[Dict]) -> int:
    """Total overlapping base-pairs between two sorted CDS interval lists."""
    total = 0
    j = 0
    for seg_a in cds_a:
        s1, e1 = seg_a['start'], seg_a['end']
        while j < len(cds_b) and cds_b[j]['end'] < s1:
            j += 1
        k = j
        while k < len(cds_b) and cds_b[k]['start'] <= e1:
            s2, e2 = cds_b[k]['start'], cds_b[k]['end']
            total += max(0, min(e1, e2) - max(s1, s2) + 1)
            k += 1
    return total


def find_best_transcript_pair(
    ens_transcripts: List[str],
    cat_transcripts: List[str],
    ens_cds: Dict[str, List[Dict]],
    cat_cds: Dict[str, List[Dict]],
) -> Tuple[Optional[str], Optional[str]]:
    """
    Return the (ensembl_tx, cat_tx) pair with the greatest CDS exon overlap.
    Falls back to the longest-CDS transcript in each gene when no spatial
    overlap is found.
    """
    ens_with_cds = [t for t in ens_transcripts if ens_cds.get(t)]
    cat_with_cds = [t for t in cat_transcripts if cat_cds.get(t)]

    if not ens_with_cds or not cat_with_cds:
        return None, None

    best_ens, best_cat, best_overlap = None, None, -1
    for e_tid in ens_with_cds:
        e_cds = ens_cds[e_tid]
        for c_tid in cat_with_cds:
            overlap = cds_intersection_bp(e_cds, cat_cds[c_tid])
            if overlap > best_overlap:
                best_overlap = overlap
                best_ens, best_cat = e_tid, c_tid

    if best_overlap == 0:
        # No spatial overlap — fall back to longest CDS in each gene
        def cds_len(cds_list):
            return sum(s['end'] - s['start'] + 1 for s in cds_list)
        best_ens = max(ens_with_cds, key=lambda t: cds_len(ens_cds[t]))
        best_cat = max(cat_with_cds, key=lambda t: cds_len(cat_cds[t]))

    return best_ens, best_cat


def calculate_cds_length(cds_list: List[Dict]) -> int:
    return sum(s['end'] - s['start'] + 1 for s in cds_list)


def get_start_stop_positions(cds_list: List[Dict]) -> Tuple[Optional[int], Optional[int]]:
    """
    Return (start_codon_pos, stop_codon_pos) from a sorted single-transcript CDS list.
    On '+': start = left edge of first exon, stop = right edge of last.
    On '-': start = right edge of last exon, stop = left edge of first.
    """
    if not cds_list:
        return None, None
    strand = cds_list[0]['strand']
    if strand == '+':
        return cds_list[0]['start'], cds_list[-1]['end']
    else:
        return cds_list[-1]['end'], cds_list[0]['start']


def load_rbh_pairs(rbh_path: str):
    pairs = []
    with open(rbh_path, 'r') as f:
        header = f.readline().strip().split('\t')
        try:
            ens_id_idx = header.index('ensembl_id')
            cat_id_idx = header.index('cat_id')
            ens_bio_idx = header.index('ensembl_biotype') if 'ensembl_biotype' in header else None
            cat_bio_idx = header.index('cat_biotype') if 'cat_biotype' in header else None
        except ValueError as e:
            print(f"Error: Required column not found: {e}", file=sys.stderr)
            sys.exit(1)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(ens_id_idx, cat_id_idx):
                pairs.append((
                    parts[ens_id_idx],
                    parts[cat_id_idx],
                    parts[ens_bio_idx] if ens_bio_idx is not None else 'unknown',
                    parts[cat_bio_idx] if cat_bio_idx is not None else 'unknown',
                ))
    return pairs


def assess_pair(e_cds: List[Dict], c_cds: List[Dict],
                tolerance: int = 3) -> Dict:
    """Assess coding integrity for a single matched transcript pair."""
    result = {
        'has_ensembl_cds': len(e_cds) > 0,
        'has_cat_cds': len(c_cds) > 0,
        'cds_length_ensembl': 0,
        'cds_length_cat': 0,
        'percent_coverage': 0.0,
        'start_codon_match': False,
        'stop_codon_match': False,
        'frameshift_detected': False,
        'length_difference': 0,
        'classification': 'No_CDS',
    }

    if not e_cds or not c_cds:
        return result

    e_len = calculate_cds_length(e_cds)
    c_len = calculate_cds_length(c_cds)
    result['cds_length_ensembl'] = e_len
    result['cds_length_cat'] = c_len
    result['length_difference'] = abs(e_len - c_len)
    if e_len > 0:
        result['percent_coverage'] = (c_len / e_len) * 100.0

    e_start, e_stop = get_start_stop_positions(e_cds)
    c_start, c_stop = get_start_stop_positions(c_cds)

    start_match = (e_start is not None and c_start is not None
                   and abs(e_start - c_start) <= tolerance)
    stop_match = (e_stop is not None and c_stop is not None
                  and abs(e_stop - c_stop) <= tolerance)
    result['start_codon_match'] = start_match
    result['stop_codon_match'] = stop_match

    is_frameshift = (result['length_difference'] % 3 != 0)
    result['frameshift_detected'] = is_frameshift

    if start_match and stop_match:
        result['classification'] = ('Full_Match' if e_len == c_len
                                    else 'Internal_Frameshift' if is_frameshift
                                    else 'Internal_Indel')
    elif not start_match and not stop_match:
        result['classification'] = 'Complex_Mismatch'
    elif not start_match and stop_match:
        result['classification'] = 'N_term_Truncation' if c_len < e_len else 'N_term_Extension'
    elif start_match and not stop_match:
        result['classification'] = 'C_term_Truncation' if c_len < e_len else 'C_term_Extension'

    return result


def main():
    args = parse_args()

    print(f"Loading Ensembl CDS by transcript from {args.ensembl_gff}", file=sys.stderr)
    ens_gene_to_tx, ens_cds = load_cds_by_transcript(args.ensembl_gff)

    print(f"Loading CAT CDS by transcript from {args.cat_gff}", file=sys.stderr)
    cat_gene_to_tx, cat_cds = load_cds_by_transcript(args.cat_gff)

    print(f"Loading RBH pairs from {args.rbh_pairs}", file=sys.stderr)
    rbh_pairs = load_rbh_pairs(args.rbh_pairs)

    protein_coding_pairs = [
        (eid, cid, ebio, cbio) for eid, cid, ebio, cbio in rbh_pairs
        if 'protein_coding' in ebio.lower() or 'protein_coding' in cbio.lower()
    ]
    print(f"Found {len(protein_coding_pairs)} protein-coding pairs; "
          "comparing best-matched transcript pairs", file=sys.stderr)

    with open(args.output, 'w') as out:
        out.write('\t'.join([
            'assembly_accession', 'sample_name',
            'ensembl_gene_id', 'cat_gene_id',
            'ensembl_biotype', 'cat_biotype',
            'ensembl_transcript_id', 'cat_transcript_id',
            'has_ensembl_cds', 'has_cat_cds',
            'classification',
            'percent_coverage',
            'cds_length_ensembl', 'cds_length_cat',
            'length_difference',
            'start_codon_match', 'stop_codon_match',
            'frameshift_detected',
        ]) + '\n')

        for ensembl_id, cat_id, ensembl_bio, cat_bio in protein_coding_pairs:
            ens_txs = ens_gene_to_tx.get(ensembl_id, [])
            cat_txs = cat_gene_to_tx.get(cat_id, [])

            best_ens_tx, best_cat_tx = find_best_transcript_pair(
                ens_txs, cat_txs, ens_cds, cat_cds
            )

            e_cds_list = ens_cds.get(best_ens_tx, []) if best_ens_tx else []
            c_cds_list = cat_cds.get(best_cat_tx, []) if best_cat_tx else []

            metrics = assess_pair(e_cds_list, c_cds_list, tolerance=3)

            out.write('\t'.join([
                args.assembly_accession,
                args.sample_name,
                ensembl_id,
                cat_id,
                ensembl_bio,
                cat_bio,
                best_ens_tx or '',
                best_cat_tx or '',
                str(metrics['has_ensembl_cds']),
                str(metrics['has_cat_cds']),
                metrics['classification'],
                f"{metrics['percent_coverage']:.2f}",
                str(metrics['cds_length_ensembl']),
                str(metrics['cds_length_cat']),
                str(metrics['length_difference']),
                str(metrics['start_codon_match']),
                str(metrics['stop_codon_match']),
                str(metrics['frameshift_detected']),
            ]) + '\n')

    print(f"Wrote coding integrity assessment to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
