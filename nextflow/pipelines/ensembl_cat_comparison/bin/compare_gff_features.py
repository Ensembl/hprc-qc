#!/usr/bin/env python3
"""
Compare Ensembl and CAT GFF feature content per assembly.

Outputs three TSVs in the given output directory:
  - feature_counts.tsv  : per source (Ensembl/CAT) global feature counts
  - gene_metrics.tsv    : per gene metrics (n_transcripts, gene_span_bp)
  - tx_metrics.tsv      : per transcript metrics (n_exons, tx_len_bp, cds_len_bp)

Usage:
    compare_gff_features.py \
        --ensembl-gff <path> \
        --cat-gff <path> \
        --output-dir <dir> \
        --assembly-accession <acc> \
        --sample-name <sample>
"""

import argparse
import gzip
import os
from collections import defaultdict


TRANSCRIPT_TYPES = {
    'transcript', 'mRNA', 'lnc_RNA', 'ncRNA',
    'miRNA', 'snoRNA', 'snRNA', 'tRNA', 'rRNA',
    'pseudogenic_transcript', 'antisense_RNA',
    'guide_RNA', 'scRNA', 'vault_RNA', 'Y_RNA',
}


def parse_args():
    ap = argparse.ArgumentParser(description="Compare GFF feature content")
    ap.add_argument("--ensembl-gff", required=True)
    ap.add_argument("--cat-gff", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--assembly-accession", required=True)
    ap.add_argument("--sample-name", required=True)
    return ap.parse_args()


def open_maybe_gzip(path: str):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')


def parse_attributes(attr: str):
    d = {}
    for item in attr.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            d[k] = v
    return d


def load_gff(path: str, source: str):
    # Containers
    gene_info = {}  # gene_id -> {biotype, start, end}
    gene_to_txs = defaultdict(set)  # gene_id -> {tx_id}
    tx_to_gene = {}  # tx_id -> gene_id
    tx_exons = defaultdict(list)  # tx_id -> [(s,e)]
    tx_cds = defaultdict(list)    # tx_id -> [(s,e)]

    n_features = defaultdict(int)

    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqid, source_field, ftype, start, end, score, strand, phase, attrs = parts
            start = int(start); end = int(end)
            a = parse_attributes(attrs)

            if ftype == 'gene' or ftype.endswith('gene'):
                gid = a.get('ID', a.get('gene_id'))
                if not gid:
                    continue
                bio = a.get('biotype', a.get('gene_biotype', 'unknown'))
                prev = gene_info.get(gid)
                if prev is None:
                    gene_info[gid] = {'biotype': bio, 'start': start, 'end': end}
                else:
                    # Expand span if repeated entries appear
                    prev['start'] = min(prev['start'], start)
                    prev['end'] = max(prev['end'], end)
                n_features['gene'] += 1

            elif ftype in TRANSCRIPT_TYPES:
                tid = a.get('ID', a.get('transcript_id'))
                gid = a.get('Parent', a.get('gene_id'))
                if not tid or not gid:
                    continue
                gene_to_txs[gid].add(tid)
                tx_to_gene[tid] = gid
                n_features['transcript'] += 1

            elif ftype == 'exon':
                tid = a.get('Parent', a.get('transcript_id'))
                if tid:
                    tx_exons[tid].append((start, end))
                    n_features['exon'] += 1

            elif ftype == 'CDS':
                tid = a.get('Parent', a.get('transcript_id'))
                if tid:
                    tx_cds[tid].append((start, end))
                    n_features['cds'] += 1

            elif ftype == 'five_prime_UTR':
                n_features['utr5'] += 1
            elif ftype == 'three_prime_UTR':
                n_features['utr3'] += 1

    # Normalize exon/CDS sorting
    for tid in list(tx_exons.keys()):
        tx_exons[tid].sort()
    for tid in list(tx_cds.keys()):
        tx_cds[tid].sort()

    return {
        'gene_info': gene_info,
        'gene_to_txs': gene_to_txs,
        'tx_to_gene': tx_to_gene,
        'tx_exons': tx_exons,
        'tx_cds': tx_cds,
        'feature_counts': n_features,
    }


def length_sum(intervals):
    return sum(e - s + 1 for s, e in intervals) if intervals else 0


def group_biotype(b: str) -> str:
    b = str(b or '').lower()
    if 'protein_coding' in b:
        return 'protein_coding'
    if 'lncrna' in b or 'lnc_rna' in b:
        return 'lncRNA'
    if 'pseudogene' in b or 'pseudogenic' in b:
        return 'pseudogene'
    if any(x in b for x in ['snrna','snorna','mirna','trna','rrna','ncrna','antisense','tec','guide_rna','scrna','vault_rna','y_rna']):
        return 'other_ncRNA'
    return 'other'


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    results = []
    for label, gff in [('Ensembl', args.ensembl_gff), ('CAT', args.cat_gff)]:
        r = load_gff(gff, label)
        results.append((label, r))

    # Write feature_counts.tsv
    with open(os.path.join(args.output_dir, 'feature_counts.tsv'), 'w') as out:
        out.write('\t'.join(['assembly_accession','sample_name','source','n_genes','n_transcripts','n_exon','n_cds','n_utr5','n_utr3']) + '\n')
        for label, r in results:
            fc = r['feature_counts']
            out.write('\t'.join([
                args.assembly_accession,
                args.sample_name,
                label,
                str(fc.get('gene', 0)),
                str(fc.get('transcript', 0)),
                str(fc.get('exon', 0)),
                str(fc.get('cds', 0)),
                str(fc.get('utr5', 0)),
                str(fc.get('utr3', 0)),
            ]) + '\n')

    # Write gene_metrics.tsv
    with open(os.path.join(args.output_dir, 'gene_metrics.tsv'), 'w') as out:
        out.write('\t'.join(['assembly_accession','sample_name','source','gene_id','biotype','biotype_group','n_transcripts','gene_span_bp']) + '\n')
        for label, r in results:
            gi = r['gene_info']; gt = r['gene_to_txs']
            for gid, meta in gi.items():
                bio = meta.get('biotype', 'unknown')
                span = (meta.get('end', 0) - meta.get('start', 1) + 1) if meta else 0
                n_tx = len(gt.get(gid, ()))
                out.write('\t'.join([
                    args.assembly_accession, args.sample_name, label,
                    gid, bio, group_biotype(bio), str(n_tx), str(span)
                ]) + '\n')

    # Write tx_metrics.tsv
    with open(os.path.join(args.output_dir, 'tx_metrics.tsv'), 'w') as out:
        out.write('\t'.join(['assembly_accession','sample_name','source','gene_id','transcript_id','biotype','biotype_group','n_exons','tx_len_bp','cds_len_bp']) + '\n')
        for label, r in results:
            gi = r['gene_info']; tx_ex = r['tx_exons']; tx_c = r['tx_cds']; t2g = r['tx_to_gene']
            for tid, exons in tx_ex.items():
                gid = t2g.get(tid)
                if not gid:
                    continue
                bio = gi.get(gid, {}).get('biotype', 'unknown')
                out.write('\t'.join([
                    args.assembly_accession, args.sample_name, label,
                    gid, tid, bio, group_biotype(bio),
                    str(len(exons)), str(length_sum(exons)), str(length_sum(tx_c.get(tid, [])))
                ]) + '\n')

if __name__ == '__main__':
    main()
