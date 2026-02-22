#!/usr/bin/env python3
"""
Compare Ensembl and CAT projected annotations against GRCh38 GENCODE reference
to identify concordant and discordant divergences (MAIN-2).

For each RBH gene pair, classifies:
  1. Both methods agree - same as reference
  2. Both methods agree - diverged from reference (high-confidence genuine variation)
  3. Ensembl diverged, CAT agrees with reference (Ensembl-specific)
  4. CAT diverged, Ensembl agrees with reference (CAT-specific)

Divergence is assessed on transcript structure (exon count) and CDS properties
(start codon position, stop codon position, CDS length).

Usage:
    compute_grch38_divergence.py \
        --ensembl-gff <path> \
        --cat-gff <path> \
        --rbh-pairs <path> \
        --gencode-gtf <path> \
        --assembly-accession <accession> \
        --sample-name <sample> \
        --output <path>
"""

import argparse
import gzip
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


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
    parser.add_argument("--cds-boundary-tolerance", type=int, default=9,
                        help="Allowed CDS boundary shift in bp (default: 9, codon-preserving)")
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


def load_gene_models_gff(gff_path: str) -> Dict[str, Dict]:
    """Load gene models from GFF3: gene_id -> {n_transcripts, n_exons_per_tx, cds_length, ...}"""
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
                    tx_to_cds[parent].append((int(parts[3]), int(parts[4]), parts[6]))

    # Build gene models
    models = {}
    for gene_id, txs in gene_to_tx.items():
        # Use longest transcript as representative
        best_tx = None
        best_exon_count = 0
        best_cds_len = 0

        for tx in txs:
            exons = sorted(tx_to_exons.get(tx, []))
            cds = sorted(tx_to_cds.get(tx, []))
            cds_len = sum(e - s + 1 for s, e, _ in cds)

            if cds_len > best_cds_len or (cds_len == best_cds_len and len(exons) > best_exon_count):
                best_tx = tx
                best_exon_count = len(exons)
                best_cds_len = cds_len

        if best_tx:
            exons = sorted(tx_to_exons.get(best_tx, []))
            cds = sorted(tx_to_cds.get(best_tx, []))
            strand = cds[0][2] if cds else '+'

            cds_start = None
            cds_stop = None
            if cds:
                if strand == '+':
                    cds_start = cds[0][0]
                    cds_stop = cds[-1][1]
                else:
                    cds_start = cds[-1][1]
                    cds_stop = cds[0][0]

            models[gene_id] = {
                'n_transcripts': len(txs),
                'n_exons': len(exons),
                'cds_length': best_cds_len,
                'cds_start': cds_start,
                'cds_stop': cds_stop,
                'transcript_id': best_tx,
                'is_coding': best_cds_len > 0,
            }

    return models


def load_gencode_models(gtf_path: str) -> Dict[str, Dict]:
    """Load reference gene models from GENCODE GTF by gene_name."""
    gene_name_map = {}  # gene_name -> gene_id
    gene_to_tx = defaultdict(list)
    tx_to_exons = defaultdict(list)
    tx_to_cds = defaultdict(list)

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
                    gene_name_map[gname] = {
                        'gene_id': gid,
                        'biotype': biotype,
                    }

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
                    tx_to_cds[tid].append((int(parts[3]), int(parts[4]), parts[6]))

    # Build models keyed by gene_name
    models = {}
    for gname, ginfo in gene_name_map.items():
        gid = ginfo['gene_id']
        txs = gene_to_tx.get(gid, [])

        best_tx = None
        best_exon_count = 0
        best_cds_len = 0

        for tx in txs:
            exons = sorted(tx_to_exons.get(tx, []))
            cds = sorted(tx_to_cds.get(tx, []))
            cds_len = sum(e - s + 1 for s, e, _ in cds)

            if cds_len > best_cds_len or (cds_len == best_cds_len and len(exons) > best_exon_count):
                best_tx = tx
                best_exon_count = len(exons)
                best_cds_len = cds_len

        if best_tx:
            exons = sorted(tx_to_exons.get(best_tx, []))
            cds = sorted(tx_to_cds.get(best_tx, []))
            strand = cds[0][2] if cds else '+'

            cds_start = None
            cds_stop = None
            if cds:
                if strand == '+':
                    cds_start = cds[0][0]
                    cds_stop = cds[-1][1]
                else:
                    cds_start = cds[-1][1]
                    cds_stop = cds[0][0]

            models[gname] = {
                'n_transcripts': len(txs),
                'n_exons': best_exon_count,
                'cds_length': best_cds_len,
                'cds_start': cds_start,
                'cds_stop': cds_stop,
                'biotype': ginfo['biotype'],
                'is_coding': best_cds_len > 0,
            }

    return models


def classify_divergence(projected: Dict, reference: Dict,
                        exon_tol: int, cds_tol: int) -> Tuple[bool, bool, str]:
    """
    Compare a projected gene model against the GRCh38 reference.

    Returns:
        (structure_diverged, cds_diverged, detail_string)
    """
    if not projected or not reference:
        return False, False, 'missing_data'

    # Transcript structure divergence: exon count difference
    exon_diff = abs(projected['n_exons'] - reference['n_exons'])
    structure_diverged = exon_diff > exon_tol

    # CDS divergence
    cds_diverged = False
    cds_detail = 'non_coding'

    if projected['is_coding'] and reference['is_coding']:
        len_diff = abs(projected['cds_length'] - reference['cds_length'])

        # CDS boundary shift
        boundary_diverged = False
        if projected['cds_start'] is not None and reference['cds_start'] is not None:
            if abs(projected['cds_start'] - reference['cds_start']) > cds_tol:
                boundary_diverged = True
        if projected['cds_stop'] is not None and reference['cds_stop'] is not None:
            if abs(projected['cds_stop'] - reference['cds_stop']) > cds_tol:
                boundary_diverged = True

        cds_diverged = boundary_diverged or len_diff > cds_tol
        cds_detail = 'cds_match' if not cds_diverged else 'cds_diverged'
    elif projected['is_coding'] != reference['is_coding']:
        cds_diverged = True
        cds_detail = 'coding_status_change'

    detail = []
    if structure_diverged:
        detail.append(f'exon_diff={exon_diff}')
    if cds_diverged:
        detail.append(cds_detail)
    if not detail:
        detail.append('identical')

    return structure_diverged, cds_diverged, ';'.join(detail)


def load_rbh_pairs_with_names(rbh_path: str) -> list:
    """Load RBH pairs and extract gene names if available."""
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


def main():
    args = parse_args()

    print(f"Loading GENCODE v47 reference from {args.gencode_gtf}", file=sys.stderr)
    ref_models = load_gencode_models(args.gencode_gtf)
    print(f"Loaded {len(ref_models)} reference gene models", file=sys.stderr)

    print(f"Loading Ensembl projections from {args.ensembl_gff}", file=sys.stderr)
    ens_models = load_gene_models_gff(args.ensembl_gff)

    print(f"Loading CAT projections from {args.cat_gff}", file=sys.stderr)
    cat_models = load_gene_models_gff(args.cat_gff)

    print(f"Loading RBH pairs from {args.rbh_pairs}", file=sys.stderr)
    rbh_pairs = load_rbh_pairs_with_names(args.rbh_pairs)

    exon_tol = args.exon_count_tolerance
    cds_tol = args.cds_boundary_tolerance

    with open(args.output, 'w') as out:
        out.write('\t'.join([
            'assembly_accession', 'sample_name',
            'ensembl_gene_id', 'cat_gene_id', 'gene_name',
            'ensembl_biotype', 'ref_biotype',
            'ens_structure_diverged', 'ens_cds_diverged', 'ens_detail',
            'cat_structure_diverged', 'cat_cds_diverged', 'cat_detail',
            'divergence_category',
            'ref_n_exons', 'ens_n_exons', 'cat_n_exons',
            'ref_cds_length', 'ens_cds_length', 'cat_cds_length',
        ]) + '\n')

        n_written = 0
        n_no_ref = 0

        for pair in rbh_pairs:
            ens_id = pair['ensembl_id']
            cat_id = pair['cat_id']
            gene_name = pair['gene_name']

            ens_model = ens_models.get(ens_id)
            cat_model = cat_models.get(cat_id)

            # Look up reference by gene name
            ref_model = ref_models.get(gene_name) if gene_name else None

            if ref_model is None:
                n_no_ref += 1
                continue

            # Classify Ensembl vs reference
            ens_struct_div, ens_cds_div, ens_detail = classify_divergence(
                ens_model, ref_model, exon_tol, cds_tol
            )

            # Classify CAT vs reference
            cat_struct_div, cat_cds_div, cat_detail = classify_divergence(
                cat_model, ref_model, exon_tol, cds_tol
            )

            # Cross-tabulate into 4 categories
            ens_diverged = ens_struct_div or ens_cds_div
            cat_diverged = cat_struct_div or cat_cds_div

            if not ens_diverged and not cat_diverged:
                category = 'both_agree_reference'
            elif ens_diverged and cat_diverged:
                category = 'both_agree_diverged'
            elif ens_diverged and not cat_diverged:
                category = 'ensembl_specific_divergence'
            else:
                category = 'cat_specific_divergence'

            out.write('\t'.join([
                args.assembly_accession, args.sample_name,
                ens_id, cat_id, gene_name,
                pair.get('ensembl_biotype', ''),
                ref_model.get('biotype', ''),
                str(ens_struct_div), str(ens_cds_div), ens_detail,
                str(cat_struct_div), str(cat_cds_div), cat_detail,
                category,
                str(ref_model.get('n_exons', 0)),
                str(ens_model.get('n_exons', 0) if ens_model else 0),
                str(cat_model.get('n_exons', 0) if cat_model else 0),
                str(ref_model.get('cds_length', 0)),
                str(ens_model.get('cds_length', 0) if ens_model else 0),
                str(cat_model.get('cds_length', 0) if cat_model else 0),
            ]) + '\n')
            n_written += 1

        print(f"Wrote {n_written} divergence records ({n_no_ref} pairs had no reference match)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
