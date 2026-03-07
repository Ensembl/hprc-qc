#!/usr/bin/env python3
"""
Aggregate per-assembly transcript concordance into intron-chain classification
distributions by biotype, plus Jaccard index quantiles.

Operates at the **transcript level**: the per-gene concordance TSV already
contains counts of how many transcripts fall into each classification
(Exact_Match, Intron_Match, …, No_Match) in BOTH directions:
  - Ensembl → CAT (n_ens_*):  "can each Ensembl transcript find a CAT match?"
  - CAT → Ensembl (n_cat_*):  "can each CAT transcript find an Ensembl match?"

Both directions are aggregated.  A `direction` column distinguishes them.

This avoids the optimistic bias of the alternative gene-level approach, where
a gene with 1/10 exact-matching transcripts would be counted the same as a
gene with 10/10 exact matches.

Produces three output files (always):
  - intron_chain_by_biotype_per_assembly.tsv
        Columns: assembly_accession, direction, biotype, classification,
                 n_transcripts, pct
        Denominator: transcripts from paired genes only.
  - jaccard_by_biotype_per_assembly.tsv
        Columns: assembly_accession, biotype, n_genes, mean, median,
                 p5, p25, p75, p95
  - transcript_count_ratio_per_assembly.tsv
        Columns: assembly_accession, biotype, n_ensembl_tx, n_cat_tx,
                 ratio_cat_to_ens

And one optional output (when --all-ensembl-genes-dir is provided):
  - intron_chain_full_denom_per_assembly.tsv
        Same schema as above but denominator includes ALL Ensembl transcripts.
        Transcripts from genes with no CAT pair appear as "Gene_Not_Shared".
        (Ensembl→CAT direction only, since we don't have a full CAT gene list.)

Usage:
    aggregate_intron_chain_by_biotype.py \
        --transcript-concordance-dir <dir> \
        --output-dir <output_dir> \
        [--all-ensembl-genes-dir <dir>]   # from count_gff_transcripts.py
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
    parser.add_argument("--all-ensembl-genes-dir", default=None,
                        help="Optional: directory of *_gene_transcript_counts.tsv "
                             "files (from count_gff_transcripts.py). When provided, "
                             "produces a full-denominator output that includes "
                             "transcripts from unmatched genes as 'Gene_Not_Shared'.")
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


# Mapping from concordance-TSV column names to classification labels
_COL_TO_CLASS_ENS = [
    ('n_ens_exact',         'Exact_Match'),
    ('n_ens_intron_match',  'Intron_Match'),
    ('n_ens_subset',        'Intron_Subset'),
    ('n_ens_superset',      'Intron_Superset'),
    ('n_ens_partial_5',     'Partial_5'),
    ('n_ens_partial_3',     'Partial_3'),
    ('n_ens_other_partial', 'Other_Partial'),
    ('n_ens_unmatched',     'No_Match'),
]

# CAT → Ensembl direction (reverse)
# Note: Subset/Superset swap in reverse direction
_COL_TO_CLASS_CAT = [
    ('n_cat_exact',         'Exact_Match'),
    ('n_cat_intron_match',  'Intron_Match'),
    ('n_cat_subset',        'Intron_Subset'),
    ('n_cat_superset',      'Intron_Superset'),
    ('n_cat_partial_5',     'Partial_5'),
    ('n_cat_partial_3',     'Partial_3'),
    ('n_cat_other_partial', 'Other_Partial'),
    ('n_cat_unmatched',     'No_Match'),
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def find_gene_count_files(directory):
    """Find gene transcript count files and index by accession."""
    accession_to_file = {}
    for root, _, files in os.walk(directory):
        for fn in sorted(files):
            if fn.endswith('_gene_transcript_counts.tsv'):
                path = os.path.join(root, fn)
                m = ACCESSION_RE.search(fn)
                if m:
                    accession_to_file[m.group(1)] = path
    return accession_to_file


def load_all_ensembl_genes(filepath):
    """
    Load the output of count_gff_transcripts.py.

    Returns:
        DataFrame with columns: gene_id, biotype (grouped), n_transcripts
    """
    df = pd.read_csv(filepath, sep='\t',
                     usecols=['gene_id', 'biotype', 'n_transcripts'],
                     dtype={'gene_id': str, 'biotype': str})
    df['biotype_group'] = df['biotype'].map(group_biotype)
    return df


def process_assembly(filepath, accession):
    """
    Process a single assembly's transcript concordance file.

    Operates at the **transcript level**: sums the per-gene transcript counts
    across all genes in each biotype, in BOTH directions:
      - Ensembl → CAT (n_ens_*): "can each Ensembl transcript find a CAT match?"
      - CAT → Ensembl (n_cat_*): "can each CAT transcript find an Ensembl match?"

    A gene with 1 exact + 9 unmatched contributes 1 Exact_Match transcript and
    9 No_Match transcripts, rather than being counted as one "Exact_Match gene."

    Returns:
        class_rows: list of dicts with assembly_accession, direction, biotype,
                    classification, n_transcripts, pct
        jaccard_rows: list of dicts with assembly_accession, biotype,
                      n_genes, mean, median, p5, p25, p75, p95
        ratio_rows: list of dicts with assembly_accession, biotype,
                    n_ensembl_tx, n_cat_tx, ratio_cat_to_ens
        paired_gene_ids: set of ensembl_gene_id values present in the file
    """
    ens_count_cols = [col for col, _ in _COL_TO_CLASS_ENS]
    cat_count_cols = [col for col, _ in _COL_TO_CLASS_CAT]
    all_count_cols = ens_count_cols + cat_count_cols

    # Also read n_ensembl_transcripts and n_cat_transcripts for ratio output
    usecols = (['ensembl_gene_id', 'ensembl_biotype',
                'n_ensembl_transcripts', 'n_cat_transcripts'] +
               all_count_cols + ['avg_jaccard_index'])

    df = pd.read_csv(filepath, sep='\t', usecols=usecols,
                     dtype={'ensembl_gene_id': str, 'ensembl_biotype': str})

    if df.empty:
        return [], [], [], set()

    # Track which Ensembl gene IDs are paired (for full-denom computation)
    paired_gene_ids = set(df['ensembl_gene_id'].dropna().unique())

    # Ensure count columns are numeric
    for col in all_count_cols + ['n_ensembl_transcripts', 'n_cat_transcripts']:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

    # Group biotype
    df['biotype'] = df['ensembl_biotype'].map(group_biotype)

    # --- Transcript-level classification counts per biotype, BOTH directions ---
    directions = [
        ('Ensembl_to_CAT', _COL_TO_CLASS_ENS),
        ('CAT_to_Ensembl', _COL_TO_CLASS_CAT),
    ]

    class_rows = []
    for direction_label, col_to_class in directions:
        for biotype in BIOTYPE_ORDER:
            bio_df = df[df['biotype'] == biotype]
            if bio_df.empty:
                continue

            # Sum transcript counts across all genes in this biotype
            totals = {cls: int(bio_df[col].sum()) for col, cls in col_to_class}
            n_total_tx = sum(totals.values())

            if n_total_tx == 0:
                continue

            for cls in CLASSIFICATION_PRIORITY:
                n = totals[cls]
                class_rows.append({
                    'assembly_accession': accession,
                    'direction': direction_label,
                    'biotype': biotype,
                    'classification': cls,
                    'n_transcripts': n,
                    'pct': round(100.0 * n / n_total_tx, 4),
                })

    # --- Transcript count ratios per biotype ---
    ratio_rows = []
    for biotype in BIOTYPE_ORDER:
        bio_df = df[df['biotype'] == biotype]
        if bio_df.empty:
            continue
        n_ens = int(bio_df['n_ensembl_transcripts'].sum())
        n_cat = int(bio_df['n_cat_transcripts'].sum())
        if n_ens == 0:
            continue
        ratio_rows.append({
            'assembly_accession': accession,
            'biotype': biotype,
            'n_ensembl_tx': n_ens,
            'n_cat_tx': n_cat,
            'ratio_cat_to_ens': round(n_cat / n_ens, 4),
        })

    # --- Jaccard quantiles per biotype (stays gene-level) ---
    # avg_jaccard_index is the mean Jaccard across a gene's transcripts,
    # so the distribution here shows per-gene average overlap.
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

    return class_rows, jaccard_rows, ratio_rows, paired_gene_ids


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Index concordance files
    acc_files = find_files(args.transcript_concordance_dir)
    print(f"Found {len(acc_files)} transcript concordance files", file=sys.stderr)

    if not acc_files:
        print("ERROR: No *_transcript_concordance.tsv files found.", file=sys.stderr)
        sys.exit(1)

    # Optionally index all-ensembl-genes count files
    gene_count_files = {}
    if args.all_ensembl_genes_dir:
        gene_count_files = find_gene_count_files(args.all_ensembl_genes_dir)
        print(f"Found {len(gene_count_files)} gene transcript count files",
              file=sys.stderr)

    all_class_rows = []
    all_jaccard_rows = []
    all_ratio_rows = []
    all_full_denom_rows = []

    for i, (acc, fpath) in enumerate(sorted(acc_files.items()), 1):
        if i % 50 == 0 or i == 1:
            print(f"  Processing {i}/{len(acc_files)}: {acc}", file=sys.stderr)
        try:
            c_rows, j_rows, r_rows, paired_ids = process_assembly(fpath, acc)
            all_class_rows.extend(c_rows)
            all_jaccard_rows.extend(j_rows)
            all_ratio_rows.extend(r_rows)

            # Full-denominator: add unmatched-gene transcripts
            # (Ensembl→CAT direction only, since we don't have a full CAT gene list)
            if acc in gene_count_files:
                all_genes = load_all_ensembl_genes(gene_count_files[acc])
                unmatched = all_genes[~all_genes['gene_id'].isin(paired_ids)]

                # Sum unmatched transcripts per biotype
                unmatched_tx = (
                    unmatched.groupby('biotype_group')['n_transcripts']
                    .sum()
                    .to_dict()
                )

                # Build full-denom rows: start from Ensembl→CAT paired totals
                # then add Gene_Not_Shared and recompute pct
                paired_totals = defaultdict(lambda: defaultdict(int))
                ens_rows = [r for r in c_rows
                            if r['direction'] == 'Ensembl_to_CAT']
                for row in ens_rows:
                    bio = row['biotype']
                    paired_totals[bio][row['classification']] = row['n_transcripts']

                for biotype in BIOTYPE_ORDER:
                    paired_by_cls = paired_totals.get(biotype, {})
                    n_paired_tx = sum(paired_by_cls.values())
                    n_unmatched_tx = unmatched_tx.get(biotype, 0)
                    n_total = n_paired_tx + n_unmatched_tx

                    if n_total == 0:
                        continue

                    # Existing classifications
                    for cls in CLASSIFICATION_PRIORITY:
                        n = paired_by_cls.get(cls, 0)
                        all_full_denom_rows.append({
                            'assembly_accession': acc,
                            'direction': 'Ensembl_to_CAT',
                            'biotype': biotype,
                            'classification': cls,
                            'n_transcripts': n,
                            'pct': round(100.0 * n / n_total, 4),
                        })

                    # New category: Gene_Not_Shared
                    all_full_denom_rows.append({
                        'assembly_accession': acc,
                        'direction': 'Ensembl_to_CAT',
                        'biotype': biotype,
                        'classification': 'Gene_Not_Shared',
                        'n_transcripts': n_unmatched_tx,
                        'pct': round(100.0 * n_unmatched_tx / n_total, 4),
                    })

        except Exception as e:
            print(f"  WARNING: failed on {acc}: {e}", file=sys.stderr)

        if i % 100 == 0:
            gc.collect()

    # Write classification output (paired genes only, both directions)
    class_df = pd.DataFrame(all_class_rows)
    class_out = os.path.join(args.output_dir, 'intron_chain_by_biotype_per_assembly.tsv')
    class_df.to_csv(class_out, sep='\t', index=False)
    print(f"Wrote {len(class_df)} rows to {class_out}", file=sys.stderr)

    # Write Jaccard output
    jac_df = pd.DataFrame(all_jaccard_rows)
    jac_out = os.path.join(args.output_dir, 'jaccard_by_biotype_per_assembly.tsv')
    jac_df.to_csv(jac_out, sep='\t', index=False)
    print(f"Wrote {len(jac_df)} rows to {jac_out}", file=sys.stderr)

    # Write transcript count ratio output
    ratio_df = pd.DataFrame(all_ratio_rows)
    ratio_out = os.path.join(args.output_dir,
                             'transcript_count_ratio_per_assembly.tsv')
    ratio_df.to_csv(ratio_out, sep='\t', index=False)
    print(f"Wrote {len(ratio_df)} rows to {ratio_out}", file=sys.stderr)

    # Write full-denominator output (if gene count files were provided)
    if all_full_denom_rows:
        full_df = pd.DataFrame(all_full_denom_rows)
        full_out = os.path.join(
            args.output_dir, 'intron_chain_full_denom_per_assembly.tsv')
        full_df.to_csv(full_out, sep='\t', index=False)
        print(f"Wrote {len(full_df)} rows to {full_out}", file=sys.stderr)

    # Print summary
    if not class_df.empty:
        n_asm = class_df['assembly_accession'].nunique()
        print(f"\nSummary: {n_asm} assemblies (transcript-level aggregation)",
              file=sys.stderr)

        for direction in ['Ensembl_to_CAT', 'CAT_to_Ensembl']:
            dir_df = class_df[class_df['direction'] == direction]
            if dir_df.empty:
                continue
            medians = (
                dir_df.groupby(['biotype', 'classification'])['pct']
                .median()
                .unstack(fill_value=0)
                .reindex(index=BIOTYPE_ORDER,
                         columns=CLASSIFICATION_PRIORITY, fill_value=0)
            )
            print(f"\nMedian % of transcripts — {direction} (paired genes):",
                  file=sys.stderr)
            print(medians.round(1).to_string(), file=sys.stderr)

        # Transcript count ratios
        if not ratio_df.empty:
            ratio_med = (
                ratio_df.groupby('biotype')['ratio_cat_to_ens']
                .median()
                .reindex(BIOTYPE_ORDER)
            )
            print("\nMedian CAT/Ensembl transcript count ratio by biotype:",
                  file=sys.stderr)
            print(ratio_med.round(2).to_string(), file=sys.stderr)

        if all_full_denom_rows:
            full_cls = CLASSIFICATION_PRIORITY + ['Gene_Not_Shared']
            full_medians = (
                full_df.groupby(['biotype', 'classification'])['pct']
                .median()
                .unstack(fill_value=0)
                .reindex(index=BIOTYPE_ORDER, columns=full_cls, fill_value=0)
            )
            print("\nMedian % of transcripts (full denominator, Ensembl→CAT):",
                  file=sys.stderr)
            print(full_medians.round(1).to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
