#!/usr/bin/env python3
"""
Aggregate per-assembly GRCh38 divergence results into cohort-level summaries for MAIN-2.

Produces:
  - grch38_divergence_cross_tab.tsv        - 2x2 category counts (agreement vs method-specific)
  - grch38_divergence_per_assembly.tsv     - per-assembly divergence category counts
  - grch38_divergence_by_biotype.tsv       - breakdown by biotype

Usage:
    aggregate_grch38_divergence.py \
        --input-dir <dir> \
        --output-dir <output_dir>
"""

import argparse
import glob
import os
import sys
from collections import Counter

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate GRCh38 divergence data")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*_grch38_divergence.tsv")))
    if not tsv_files:
        print("ERROR: No grch38_divergence TSV files found", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(tsv_files)} divergence files", file=sys.stderr)

    assembly_rows = []
    category_by_biotype = []
    all_category_counts = Counter()

    for f in tsv_files:
        try:
            df = pd.read_csv(f, sep='\t', dtype=str)
        except Exception as e:
            print(f"WARNING: Could not read {f}: {e}", file=sys.stderr)
            continue

        if df.empty:
            continue

        accession = df['assembly_accession'].iloc[0]
        sample = df['sample_name'].iloc[0]
        n_total = len(df)

        cats = df['divergence_category'].value_counts()
        n_both_ref = int(cats.get('both_agree_reference', 0))
        n_both_div = int(cats.get('both_agree_diverged', 0))
        n_ens_spec = int(cats.get('ensembl_specific_divergence', 0))
        n_cat_spec = int(cats.get('cat_specific_divergence', 0))

        for cat, cnt in cats.items():
            all_category_counts[cat] += int(cnt)

        assembly_rows.append({
            'assembly_accession': accession,
            'sample_name': sample,
            'n_genes_assessed': n_total,
            'n_both_agree_reference': n_both_ref,
            'n_both_agree_diverged': n_both_div,
            'n_ensembl_specific': n_ens_spec,
            'n_cat_specific': n_cat_spec,
            'pct_both_agree_reference': round(n_both_ref / n_total * 100, 2) if n_total > 0 else 0,
            'pct_both_agree_diverged': round(n_both_div / n_total * 100, 2) if n_total > 0 else 0,
            'pct_ensembl_specific': round(n_ens_spec / n_total * 100, 2) if n_total > 0 else 0,
            'pct_cat_specific': round(n_cat_spec / n_total * 100, 2) if n_total > 0 else 0,
        })

        # Biotype stratification
        if 'ref_biotype' in df.columns:
            for (biotype, category), count in df.groupby(['ref_biotype', 'divergence_category']).size().items():
                category_by_biotype.append({
                    'assembly_accession': accession,
                    'biotype': biotype,
                    'divergence_category': category,
                    'count': int(count),
                })

    # Write per-assembly summary
    pd.DataFrame(assembly_rows).to_csv(
        os.path.join(args.output_dir, 'grch38_divergence_per_assembly.tsv'),
        sep='\t', index=False
    )

    # Write aggregate cross-tabulation
    total = sum(all_category_counts.values())
    cross_tab = []
    for cat in ['both_agree_reference', 'both_agree_diverged',
                'ensembl_specific_divergence', 'cat_specific_divergence']:
        cnt = all_category_counts.get(cat, 0)
        cross_tab.append({
            'divergence_category': cat,
            'total_gene_assembly_instances': cnt,
            'pct_of_total': round(cnt / total * 100, 2) if total > 0 else 0,
        })

    pd.DataFrame(cross_tab).to_csv(
        os.path.join(args.output_dir, 'grch38_divergence_cross_tab.tsv'),
        sep='\t', index=False
    )

    # Write biotype breakdown
    if category_by_biotype:
        bio_df = pd.DataFrame(category_by_biotype)
        # Aggregate across assemblies
        bio_agg = bio_df.groupby(['biotype', 'divergence_category'])['count'].sum().reset_index()
        bio_agg.to_csv(
            os.path.join(args.output_dir, 'grch38_divergence_by_biotype.tsv'),
            sep='\t', index=False
        )

    print(f"Wrote divergence aggregation for {len(assembly_rows)} assemblies", file=sys.stderr)


if __name__ == "__main__":
    main()
