#!/usr/bin/env python3
"""
Compare Ensembl and CAT projected annotations against GRCh38 GENCODE reference
to identify concordant and discordant divergences (MAIN-2).

For each RBH gene pair, computes valid cross-assembly metrics (lengths and counts —
never raw coordinates, which are incomparable across assembly spaces) and classifies
divergence relative to the GRCh38 reference.

Metrics computed per annotation source (ref/ens/cat):
  - Exon count, transcript length, exon length stats (mean/median/min/max)
  - Intron count and intron length stats
  - CDS length (coding genes only)

Divergence classification per source vs reference:
  - structure_match: exon count within tolerance
  - length_match: transcript length within tolerance
  - cds_match: CDS length within tolerance (coding only)
  - coding_status_match: both coding or both non-coding
  - cds_frameshift: CDS length delta not divisible by 3

Summary divergence_category (4-way, consistent with prior analysis):
  both_agree_reference | both_agree_diverged |
  ensembl_specific_divergence | cat_specific_divergence

Usage:
    compute_grch38_divergence.py \\
        --ensembl-gff <path> \\
        --cat-gff <path> \\
        --rbh-pairs <path> \\
        --gencode-gtf <path> \\
        --assembly-accession <accession> \\
        --sample-name <sample> \\
        --output <path>
"""

import argparse
import gzip
import statistics
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


def parse_args():
    parser = argparse.ArgumentParser(description="Compute GRCh38 divergence for RBH gene pairs")
    parser.add_argument("--ensembl-gff", required=True)
    parser.add_argument("--cat-gff", required=True)
    parser.add_argument("--rbh-pairs", required=True)
    parser.add_argument("--gencode-gtf", required=True, help="GENCODE v47 GTF for GRCh38 reference")
    parser.add_argument("--output", required=True)
    parser.add_argument("--assembly-accession", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--exon-count-tolerance", type=int, default=0,
                        help="Allowed exon count difference (default: 0)")
    parser.add_argument("--length-tolerance", type=int, default=0,
                        help="Allowed transcript length difference in bp (default: 0)")
    parser.add_argument("--cds-length-tolerance", type=int, default=0,
                        help="Allowed CDS length difference in bp (default: 0)")
    return parser.parse_args()


def open_maybe_gzip(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    attrs = {}
    for item in attr_string.strip().split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key.strip()] = val.strip()
    return attrs


def parse_gtf_attributes(attr_string: str) -> Dict[str, str]:
    attrs = {}
    for item in attr_string.strip().rstrip(';').split(';'):
        item = item.strip()
        if ' ' in item:
            key, val = item.split(' ', 1)
            attrs[key] = val.strip('"')
    return attrs


TRANSCRIPT_TYPES = {
    'transcript', 'mRNA', 'lnc_RNA', 'ncRNA',
    'miRNA', 'snoRNA', 'snRNA', 'tRNA', 'rRNA',
    'pseudogenic_transcript', 'antisense_RNA',
    'guide_RNA', 'scRNA', 'vault_RNA', 'Y_RNA',
}


def _exon_stats(exon_list: List[Tuple[int, int]]) -> Dict:
    """Compute length stats from a sorted list of (start, end) exon tuples (1-based inclusive)."""
    if not exon_list:
        return {
            'n_exons': 0,
            'transcript_length': 0,
            'mean_exon_length': 0.0,
            'median_exon_length': 0.0,
            'min_exon_length': 0,
            'max_exon_length': 0,
            'n_introns': 0,
            'mean_intron_length': None,
            'median_intron_length': None,
            'min_intron_length': None,
            'max_intron_length': None,
        }

    lengths = [e - s + 1 for s, e in exon_list]
    transcript_length = sum(lengths)
    n_exons = len(lengths)

    intron_lengths = []
    for i in range(len(exon_list) - 1):
        intron_len = exon_list[i + 1][0] - exon_list[i][1] - 1
        if intron_len >= 0:
            intron_lengths.append(intron_len)

    return {
        'n_exons': n_exons,
        'transcript_length': transcript_length,
        'mean_exon_length': statistics.mean(lengths),
        'median_exon_length': statistics.median(lengths),
        'min_exon_length': min(lengths),
        'max_exon_length': max(lengths),
        'n_introns': len(intron_lengths),
        'mean_intron_length': statistics.mean(intron_lengths) if intron_lengths else None,
        'median_intron_length': statistics.median(intron_lengths) if intron_lengths else None,
        'min_intron_length': min(intron_lengths) if intron_lengths else None,
        'max_intron_length': max(intron_lengths) if intron_lengths else None,
    }


def load_gene_models_gff(gff_path: str) -> Dict[str, Dict]:
    """Load gene models from GFF3. Returns gene_id -> model dict with length/count metrics."""
    gene_to_tx = defaultdict(list)
    tx_to_exons = defaultdict(list)
    tx_to_cds = defaultdict(list)
    tx_to_gene = {}

    with open_maybe_gzip(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            ftype = parts[2]
            attrs = parse_gff_attributes(parts[8])

            if ftype in TRANSCRIPT_TYPES:
                tid = attrs.get('ID', attrs.get('transcript_id'))
                parent = attrs.get('Parent', attrs.get('gene_id'))
                if tid and parent:
                    gene_to_tx[parent].append(tid)
                    tx_to_gene[tid] = parent

            elif ftype == 'exon':
                parent = attrs.get('Parent', attrs.get('transcript_id'))
                if parent:
                    tx_to_exons[parent].append((int(parts[3]), int(parts[4])))

            elif ftype == 'CDS':
                parent = attrs.get('Parent', attrs.get('transcript_id'))
                if parent:
                    tx_to_cds[parent].append((int(parts[3]), int(parts[4])))

    models = {}
    for gene_id, txs in gene_to_tx.items():
        best_tx = None
        best_cds_len = 0
        best_exon_count = 0

        for tx in txs:
            exons = tx_to_exons.get(tx, [])
            cds = tx_to_cds.get(tx, [])
            cds_len = sum(e - s + 1 for s, e in cds)

            if cds_len > best_cds_len or (cds_len == best_cds_len and len(exons) > best_exon_count):
                best_tx = tx
                best_cds_len = cds_len
                best_exon_count = len(exons)

        if best_tx is None and txs:
            best_tx = txs[0]

        if best_tx:
            exons = sorted(tx_to_exons.get(best_tx, []))
            cds = tx_to_cds.get(best_tx, [])
            cds_len = sum(e - s + 1 for s, e in cds)

            model = _exon_stats(exons)
            model['cds_length'] = cds_len
            model['is_coding'] = cds_len > 0
            models[gene_id] = model

    return models


def load_gencode_models(gtf_path: str) -> Dict[str, Dict]:
    """Load reference gene models from GENCODE GTF keyed by gene_name.

    GENCODE GTF CDS features include the stop codon in their coordinate range,
    whereas Ensembl/CAT GFF3 CDS features exclude it.  To normalise to the
    GFF3 convention (CDS = coding sequence without stop codon) we subtract
    any annotated stop_codon intervals from the CDS length.
    """
    gene_name_map = {}
    gene_to_tx = defaultdict(list)
    tx_to_exons = defaultdict(list)
    tx_to_cds = defaultdict(list)
    tx_to_stop = defaultdict(list)   # stop_codon intervals to subtract
    gene_id_to_biotype = {}

    with open_maybe_gzip(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            ftype = parts[2]
            attrs = parse_gtf_attributes(parts[8])

            if ftype == 'gene':
                gname = attrs.get('gene_name', '')
                gid = attrs.get('gene_id', '')
                biotype = attrs.get('gene_type', '')
                if gname:
                    gene_name_map[gname] = gid
                    gene_id_to_biotype[gid] = biotype

            elif ftype == 'transcript':
                tid = attrs.get('transcript_id', '')
                gid = attrs.get('gene_id', '')
                if tid and gid:
                    gene_to_tx[gid].append(tid)

            elif ftype == 'exon':
                tid = attrs.get('transcript_id', '')
                if tid:
                    tx_to_exons[tid].append((int(parts[3]), int(parts[4])))

            elif ftype == 'CDS':
                tid = attrs.get('transcript_id', '')
                if tid:
                    tx_to_cds[tid].append((int(parts[3]), int(parts[4])))

            elif ftype == 'stop_codon':
                tid = attrs.get('transcript_id', '')
                if tid:
                    tx_to_stop[tid].append((int(parts[3]), int(parts[4])))

    models = {}
    for gname, gid in gene_name_map.items():
        txs = gene_to_tx.get(gid, [])

        best_tx = None
        best_cds_len = 0
        best_exon_count = 0

        for tx in txs:
            exons = tx_to_exons.get(tx, [])
            cds = tx_to_cds.get(tx, [])
            cds_len = sum(e - s + 1 for s, e in cds)

            if cds_len > best_cds_len or (cds_len == best_cds_len and len(exons) > best_exon_count):
                best_tx = tx
                best_cds_len = cds_len
                best_exon_count = len(exons)

        if best_tx is None and txs:
            best_tx = txs[0]

        if best_tx:
            exons = sorted(tx_to_exons.get(best_tx, []))
            cds = tx_to_cds.get(best_tx, [])
            stop_len = sum(e - s + 1 for s, e in tx_to_stop.get(best_tx, []))
            cds_len = max(0, sum(e - s + 1 for s, e in cds) - stop_len)

            model = _exon_stats(exons)
            model['cds_length'] = cds_len
            model['is_coding'] = cds_len > 0
            model['biotype'] = gene_id_to_biotype.get(gid, '')
            models[gname] = model

    return models


def _fmt(val) -> str:
    """Format a value for TSV output: None -> '', float -> 2 decimal places."""
    if val is None:
        return ''
    if isinstance(val, float):
        return f'{val:.4f}'
    return str(val)


def compare_to_ref(projected: Optional[Dict], reference: Dict,
                   exon_tol: int, len_tol: int, cds_tol: int) -> Dict:
    """
    Compare a projected gene model against the GRCh38 reference using valid
    cross-assembly metrics (counts and lengths only — no coordinate comparison).

    Returns a dict of signed deltas, match flags, direction labels, and an overall
    diverged flag.
    """
    if projected is None:
        return {
            'delta_n_exons': None,
            'delta_transcript_length': None,
            'delta_mean_exon_length': None,
            'delta_cds_length': None,
            'cds_frameshift': None,
            'cds_change_type': None,
            'structure_match': None,
            'length_match': None,
            'cds_match': None,
            'coding_status_match': None,
            'length_direction': 'missing',
            'cds_direction': 'missing',
            'diverged': None,
        }

    delta_n_exons = projected['n_exons'] - reference['n_exons']
    delta_tx_len = projected['transcript_length'] - reference['transcript_length']
    delta_mean_exon = (
        (projected['mean_exon_length'] - reference['mean_exon_length'])
        if reference['n_exons'] > 0 and projected['n_exons'] > 0
        else None
    )
    delta_cds = projected['cds_length'] - reference['cds_length']

    structure_match = abs(delta_n_exons) <= exon_tol
    length_match = abs(delta_tx_len) <= len_tol
    coding_status_match = projected['is_coding'] == reference['is_coding']

    # CDS match only meaningful if both are coding
    if reference['is_coding'] and projected['is_coding']:
        cds_match = abs(delta_cds) <= cds_tol
        cds_frameshift = (delta_cds != 0) and (abs(delta_cds) % 3 != 0)
        if delta_cds == 0:
            cds_direction = 'match'
            cds_change_type = 'exact_match'
        elif cds_frameshift:
            cds_direction = 'frameshift'
            cds_change_type = 'frameshift_shorter' if delta_cds < 0 else 'frameshift_longer'
        elif delta_cds > 0:
            cds_direction = 'longer'
            cds_change_type = 'in_frame_longer'
        else:
            cds_direction = 'shorter'
            cds_change_type = 'in_frame_shorter'
    elif not reference['is_coding'] and not projected['is_coding']:
        cds_match = True
        cds_frameshift = False
        cds_direction = 'non_coding'
        cds_change_type = 'non_coding'
    else:
        # Coding status changed
        cds_match = False
        cds_frameshift = False
        cds_direction = 'coding_status_change'
        cds_change_type = 'coding_lost' if reference['is_coding'] else 'coding_gained'

    if delta_tx_len == 0:
        length_direction = 'same'
    elif delta_tx_len > 0:
        length_direction = 'longer'
    else:
        length_direction = 'shorter'

    diverged = (
        not structure_match
        or not length_match
        or not cds_match
        or not coding_status_match
    )

    return {
        'delta_n_exons': delta_n_exons,
        'delta_transcript_length': delta_tx_len,
        'delta_mean_exon_length': delta_mean_exon,
        'delta_cds_length': delta_cds,
        'cds_frameshift': cds_frameshift,
        'cds_change_type': cds_change_type,
        'structure_match': structure_match,
        'length_match': length_match,
        'cds_match': cds_match,
        'coding_status_match': coding_status_match,
        'length_direction': length_direction,
        'cds_direction': cds_direction,
        'diverged': diverged,
    }


def load_rbh_pairs(rbh_path: str) -> list:
    pairs = []
    with open(rbh_path, 'r') as f:
        header = f.readline().strip().split('\t')
        col_idx = {col: i for i, col in enumerate(header)}

        ens_idx = col_idx.get('ensembl_id', col_idx.get('ensembl_gene_id'))
        cat_idx = col_idx.get('cat_id', col_idx.get('cat_gene_id'))
        name_idx = col_idx.get('gene_name', col_idx.get('ensembl_name'))
        ens_bio_idx = col_idx.get('ensembl_biotype')
        cat_bio_idx = col_idx.get('cat_biotype')

        if ens_idx is None or cat_idx is None:
            print(f"ERROR: Cannot find ensembl_id/cat_id columns in {rbh_path}", file=sys.stderr)
            sys.exit(1)

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(ens_idx, cat_idx):
                continue
            pairs.append({
                'ensembl_id': parts[ens_idx],
                'cat_id': parts[cat_idx],
                'gene_name': parts[name_idx] if name_idx is not None and name_idx < len(parts) else '',
                'ensembl_biotype': parts[ens_bio_idx] if ens_bio_idx is not None and ens_bio_idx < len(parts) else '',
                'cat_biotype': parts[cat_bio_idx] if cat_bio_idx is not None and cat_bio_idx < len(parts) else '',
            })
    return pairs


# Column definitions for the output TSV
METRIC_COLS = [
    'n_exons', 'transcript_length',
    'mean_exon_length', 'median_exon_length', 'min_exon_length', 'max_exon_length',
    'n_introns', 'mean_intron_length', 'median_intron_length',
    'min_intron_length', 'max_intron_length',
    'cds_length',
]

DELTA_COLS = [
    'delta_n_exons', 'delta_transcript_length', 'delta_mean_exon_length',
    'delta_cds_length', 'cds_frameshift', 'cds_change_type',
    'structure_match', 'length_match', 'cds_match', 'coding_status_match',
    'length_direction', 'cds_direction', 'diverged',
]

HEADER = (
    ['assembly_accession', 'sample_name',
     'ensembl_gene_id', 'cat_gene_id', 'gene_name',
     'ensembl_biotype', 'ref_biotype', 'is_coding']
    + [f'ref_{c}' for c in METRIC_COLS]
    + [f'ens_{c}' for c in METRIC_COLS]
    + [f'cat_{c}' for c in METRIC_COLS]
    + [f'ens_{c}' for c in DELTA_COLS]
    + [f'cat_{c}' for c in DELTA_COLS]
    + ['ens_cat_delta_n_exons', 'ens_cat_delta_transcript_length', 'ens_cat_delta_cds_length']
    + ['divergence_category']
)


def model_metric_values(model: Optional[Dict]) -> List[str]:
    """Return formatted metric values for one annotation source."""
    if model is None:
        return [''] * len(METRIC_COLS)
    return [_fmt(model.get(c)) for c in METRIC_COLS]


def delta_values(d: Dict) -> List[str]:
    """Return formatted delta/flag values for one comparison."""
    return [_fmt(d.get(c)) for c in DELTA_COLS]


def main():
    args = parse_args()

    print(f"Loading GENCODE v47 reference from {args.gencode_gtf}", file=sys.stderr)
    ref_models = load_gencode_models(args.gencode_gtf)
    print(f"Loaded {len(ref_models)} reference gene models", file=sys.stderr)

    print(f"Loading Ensembl projections from {args.ensembl_gff}", file=sys.stderr)
    ens_models = load_gene_models_gff(args.ensembl_gff)
    print(f"Loaded {len(ens_models)} Ensembl gene models", file=sys.stderr)

    print(f"Loading CAT projections from {args.cat_gff}", file=sys.stderr)
    cat_models = load_gene_models_gff(args.cat_gff)
    print(f"Loaded {len(cat_models)} CAT gene models", file=sys.stderr)

    print(f"Loading RBH pairs from {args.rbh_pairs}", file=sys.stderr)
    rbh_pairs = load_rbh_pairs(args.rbh_pairs)

    exon_tol = args.exon_count_tolerance
    len_tol = args.length_tolerance
    cds_tol = args.cds_length_tolerance

    n_written = 0
    n_no_ref = 0

    with open(args.output, 'w') as out:
        out.write('\t'.join(HEADER) + '\n')

        for pair in rbh_pairs:
            ens_id = pair['ensembl_id']
            cat_id = pair['cat_id']
            gene_name = pair['gene_name']

            ref_model = ref_models.get(gene_name) if gene_name else None
            if ref_model is None:
                n_no_ref += 1
                continue

            ens_model = ens_models.get(ens_id)
            cat_model = cat_models.get(cat_id)

            ens_cmp = compare_to_ref(ens_model, ref_model, exon_tol, len_tol, cds_tol)
            cat_cmp = compare_to_ref(cat_model, ref_model, exon_tol, len_tol, cds_tol)

            # 4-way divergence category
            ens_div = ens_cmp['diverged']
            cat_div = cat_cmp['diverged']

            if ens_div is None or cat_div is None:
                category = 'insufficient_data'
            elif not ens_div and not cat_div:
                category = 'both_agree_reference'
            elif ens_div and cat_div:
                category = 'both_agree_diverged'
            elif ens_div:
                category = 'ensembl_specific_divergence'
            else:
                category = 'cat_specific_divergence'

            # Ensembl vs CAT direct deltas
            if ens_model and cat_model:
                ens_cat_delta_n = ens_model['n_exons'] - cat_model['n_exons']
                ens_cat_delta_tx = ens_model['transcript_length'] - cat_model['transcript_length']
                ens_cat_delta_cds = ens_model['cds_length'] - cat_model['cds_length']
            else:
                ens_cat_delta_n = ens_cat_delta_tx = ens_cat_delta_cds = None

            row = (
                [args.assembly_accession, args.sample_name,
                 ens_id, cat_id, gene_name,
                 pair.get('ensembl_biotype', ''),
                 ref_model.get('biotype', ''),
                 str(ref_model.get('is_coding', ''))]
                + model_metric_values(ref_model)
                + model_metric_values(ens_model)
                + model_metric_values(cat_model)
                + delta_values(ens_cmp)
                + delta_values(cat_cmp)
                + [_fmt(ens_cat_delta_n), _fmt(ens_cat_delta_tx), _fmt(ens_cat_delta_cds)]
                + [category]
            )
            out.write('\t'.join(row) + '\n')
            n_written += 1

    print(
        f"Wrote {n_written} divergence records ({n_no_ref} pairs had no reference match)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
