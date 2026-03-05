#!/usr/bin/env python3
"""
Aggregate per-assembly transcript concordance into intron-chain classification
distributions by biotype, plus Jaccard index quantiles.

For each gene pair the script determines the *best* intron-chain classification
across its Ensembl transcripts (Ensembl→CAT direction). Classifications are
ordered by priority:

    Exact_Match > Intron_Match > Intron_Subset > Intron_Superset
    > Partial_5 > Partial_3 > Other_Partial > No_Match

The gene-level best class is then aggregated per assembly × biotype.

Produces two output files:
  - intron_chain_by_biotype_per_assembly.tsv
        Columns: assembly_accession, biotype, classification, n_genes, pct
  - jaccard_by_biotype_per_assembly.tsv
        Columns: assembly_accession, biotype, n_genes, mean, median,
                 p5, p25, p75, p95

Usage:
    aggregate_intron_chain_by_biotype.py \
        --transcript-concordance-dir <dir> \
        --output-dir <output_dir>
"""

import argparse
import csv
import gc
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ACCESSION_RE = re.compile(r'(GC[AF]_\d+\.\d+)')

# Priority ordering (lower index = better match)
CLASSIFICATION_PRIORITY = [
    'Exact_Match',
    'Intron_Match',
    'Intron_Subset',
    'Intron_Superset',
    'Partial_5',
    'Partial_3',
    'Other_Partial',
    'No_Match',
]
_PRIORITY_MAP = {c: i for i, c in enumerate(CLASSIFICATION_PRIORITY)}

# Biotype grouping
BIOTYPE_ORDER = ['protein_coding', 'lncRNA', 'pseudogene', 'other_ncRNA', 'other']


def group_biotype(b: str) -> str:
    b = str(b or '').lower()
    if 'protein_coding' in b:
        return 'protein_coding'
    if 'lncrna' in b or 'lnc_rna' in b:
        return 'lncRNA'
    if 'pseudogene' in b or 'pseudogenic' in b:
        return 'pseudogene'
    if any(x in b for x in ['snrna', 'snorna', 'mirna', 'trna', 'rrna',
                             'ncrna', 'antisense', 'tec', 'guide_rna',
                             'scrna', 'vault_rna', 'y_rna']):
        return 'other_ncRNA'
    return 'other'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate intron-chain classification by biotype"
    )
    parser.add_argument("--transcript-concordance-dir", required=True,
                        help="Directory containing *_transcript_concordance.tsv files")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for aggregated TSVs")
    return parser.parse_args()


def extract_accession(filepath):
    """Extract assembly accession from filename."""
    m = ACCESSION_RE.search(os.path.basename(filepath))
    if m:
        return m.group(1)
    # Fallback: try first row
    try:
        with open(filepath) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            row = next(reader, None)
            if row and 'assembly_accession' in row:
                return row['assembly_accession']
    except Exception:
        pass
    return None


def find_files(directory):
    """Find transcript concordance files and index by accession."""
    accession_to_file = {}
    for root, _, files in os.walk(directory):
        for fn in sorted(files):
            if fn.endswith('_transcript_concordance.tsv'):
                path = os.path.join(root, fn)
                acc = extract_accession(path)
                if acc:
                    accession_to_file[acc] = path
    return accession_to_file


def best_gene_classification(row):
    """
    Determine the best intron-chain classification for a gene pair
    based on the Ensembl→CAT direction counts.

    The columns are: n_ens_exact, n_ens_intron_match, n_ens_subset,
    n_ens_superset, n_ens_partial_5, n_ens_partial_3, n_ens_other_partial,
    n_ens_unmatched.

    Returns the highest-priority classification that has count >= 1.
    """
    col_map = [
        ('n_ens_exact',         'Exact_Match'),
        ('n_ens_intron_match',  'Intron_Match'),
        ('n_ens_subset',        'Intron_Subset'),
        ('n_ens_superset',      'Intron_Superset'),
        ('n_ens_partial_5',     'Partial_5'),
        ('n_ens_partial_3',     'Partial_3'),
        ('n_ens_other_partial', 'Other_Partial'),
    ]
    for col, cls in col_map:
        if int(row.get(col, 0)) >= 1:
            return cls
    return 'No_Match'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_assembly(filepath, accession):
    """
    Process a single assembly's transcript concordance file.

    Returns:
        class_rows: list of dicts with assembly_accession, biotype,
                    classification, n_genes
        jaccard_rows: list of dicts with assembly_accession, biotype,
                      n_genes, mean, median, p5, p25, p75, p95
    """
    df = pd.read_csv(filepath, sep='\t',
                     usecols=['ensembl_biotype',
                              'n_ens_exact', 'n_ens_intron_match',
                              'n_ens_subset', 'n_ens_superset',
                              'n_ens_partial_5', 'n_ens_partial_3',
                              'n_ens_other_partial', 'n_ens_unmatched',
                              'avg_jaccard_index'],
                     dtype={'ensembl_biotype': str})

    if df.empty:
        return [], []

    # Group biotype
    df['biotype'] = df['ensembl_biotype'].map(group_biotype)

    # Assign best gene-level classification
    df['best_class'] = df.apply(best_gene_classification, axis=1)

    # --- Classification counts per biotype ---
    class_rows = []
    for biotype in BIOTYPE_ORDER:
        bio_df = df[df['biotype'] == biotype]
        n_total = len(bio_df)
        if n_total == 0:
            continue
        for cls in CLASSIFICATION_PRIORITY:
            n = int((bio_df['best_class'] == cls).sum())
            class_rows.append({
                'assembly_accession': accession,
                'biotype': biotype,
                'classification': cls,
                'n_genes': n,
                'pct': round(100.0 * n / n_total, 4) if n_total > 0 else 0.0,
            })

    # --- Jaccard quantiles per biotype ---
    jaccard_rows = []
    for biotype in BIOTYPE_ORDER:
        bio_df = df[df['biotype'] == biotype]
        jac = pd.to_numeric(bio_df['avg_jaccard_index'], errors='coerce').dropna()
        if len(jac) == 0:
            continue
        jaccard_rows.append({
            'assembly_accession': accession,
            'biotype': biotype,
            'n_genes': len(jac),
            'mean': round(float(jac.mean()), 6),
            'median': round(float(jac.median()), 6),
            'p5': round(float(jac.quantile(0.05)), 6),
            'p25': round(float(jac.quantile(0.25)), 6),
            'p75': round(float(jac.quantile(0.75)), 6),
            'p95': round(float(jac.quantile(0.95)), 6),
        })

    return class_rows, jaccard_rows


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Index files
    acc_files = find_files(args.transcript_concordance_dir)
    print(f"Found {len(acc_files)} transcript concordance files", file=sys.stderr)

    if not acc_files:
        print("ERROR: No *_transcript_concordance.tsv files found.", file=sys.stderr)
        sys.exit(1)

    all_class_rows = []
    all_jaccard_rows = []

    for i, (acc, fpath) in enumerate(sorted(acc_files.items()), 1):
        if i % 50 == 0 or i == 1:
            print(f"  Processing {i}/{len(acc_files)}: {acc}", file=sys.stderr)
        try:
            c_rows, j_rows = process_assembly(fpath, acc)
            all_class_rows.extend(c_rows)
            all_jaccard_rows.extend(j_rows)
        except Exception as e:
            print(f"  WARNING: failed on {acc}: {e}", file=sys.stderr)

        if i % 100 == 0:
            gc.collect()

    # Write classification output
    class_df = pd.DataFrame(all_class_rows)
    class_out = os.path.join(args.output_dir, 'intron_chain_by_biotype_per_assembly.tsv')
    class_df.to_csv(class_out, sep='\t', index=False)
    print(f"Wrote {len(class_df)} rows to {class_out}", file=sys.stderr)

    # Write Jaccard output
    jac_df = pd.DataFrame(all_jaccard_rows)
    jac_out = os.path.join(args.output_dir, 'jaccard_by_biotype_per_assembly.tsv')
    jac_df.to_csv(jac_out, sep='\t', index=False)
    print(f"Wrote {len(jac_df)} rows to {jac_out}", file=sys.stderr)

    # Print summary
    if not class_df.empty:
        n_asm = class_df['assembly_accession'].nunique()
        print(f"\nSummary: {n_asm} assemblies", file=sys.stderr)
        # Median % per biotype × classification
        medians = (
            class_df.groupby(['biotype', 'classification'])['pct']
            .median()
            .unstack(fill_value=0)
            .reindex(index=BIOTYPE_ORDER, columns=CLASSIFICATION_PRIORITY, fill_value=0)
        )
        print("\nMedian % across assemblies:", file=sys.stderr)
        print(medians.round(1).to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
