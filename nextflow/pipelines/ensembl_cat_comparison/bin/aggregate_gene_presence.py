#!/usr/bin/env python3
"""
Aggregate per-assembly gene presence/absence files into cohort-level summaries.

Produces:
  1. sankey_level1_gene_presence.tsv  - per-gene counts across assemblies
  2. biotype_enrichment_exclusive_genes.tsv - biotype breakdown of method-exclusive genes
  3. gene_presence_per_assembly_summary.tsv - per-assembly presence counts

Usage:
    aggregate_gene_presence.py \
        --input-dir <dir_with_per_assembly_tsvs> \
        --output-dir <output_dir>
"""

import argparse
import glob
import os
import sys
from collections import Counter, defaultdict

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate gene presence/absence across assemblies")
    parser.add_argument("--input-dir", required=True, help="Directory containing per-assembly gene_presence TSV files")
    parser.add_argument("--output-dir", required=True, help="Output directory for aggregated results")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*_gene_presence.tsv")))
    if not tsv_files:
        print("ERROR: No gene_presence TSV files found", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(tsv_files)} gene presence files", file=sys.stderr)

    # --- Per-assembly summary ---
    assembly_summaries = []

    # --- Per-gene aggregation ---
    gene_both_count = Counter()          # gene_name -> n_assemblies where both
    gene_ensembl_only_count = Counter()  # gene_name -> n_assemblies where ensembl only
    gene_cat_only_count = Counter()      # gene_name -> n_assemblies where cat only
    gene_biotype_ensembl = {}            # gene_name -> biotype from ensembl
    gene_biotype_cat = {}                # gene_name -> biotype from cat

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

        # Per-assembly counts
        both = ((df['present_in_ensembl'] == 'True') & (df['present_in_cat'] == 'True')).sum()
        ens_only = ((df['present_in_ensembl'] == 'True') & (df['present_in_cat'] == 'False')).sum()
        cat_only = ((df['present_in_ensembl'] == 'False') & (df['present_in_cat'] == 'True')).sum()

        assembly_summaries.append({
            'assembly_accession': accession,
            'sample_name': sample,
            'n_both': int(both),
            'n_ensembl_only': int(ens_only),
            'n_cat_only': int(cat_only),
            'n_total': int(both + ens_only + cat_only),
            'pct_both': round(both / (both + ens_only + cat_only) * 100, 2) if (both + ens_only + cat_only) > 0 else 0,
        })

        # Per-gene aggregation
        for _, row in df.iterrows():
            gn = row['gene_name']
            in_ens = row['present_in_ensembl'] == 'True'
            in_cat = row['present_in_cat'] == 'True'

            if in_ens and in_cat:
                gene_both_count[gn] += 1
            elif in_ens:
                gene_ensembl_only_count[gn] += 1
            elif in_cat:
                gene_cat_only_count[gn] += 1

            # Store biotypes (first seen wins)
            if in_ens and gn not in gene_biotype_ensembl:
                gene_biotype_ensembl[gn] = row.get('ensembl_biotype', 'unknown')
            if in_cat and gn not in gene_biotype_cat:
                gene_biotype_cat[gn] = row.get('cat_biotype', 'unknown')

    n_assemblies = len(tsv_files)

    # --- Write per-assembly summary ---
    pd.DataFrame(assembly_summaries).to_csv(
        os.path.join(args.output_dir, 'gene_presence_per_assembly_summary.tsv'),
        sep='\t', index=False
    )
    print(f"Wrote per-assembly summary ({len(assembly_summaries)} assemblies)", file=sys.stderr)

    # --- Write per-gene Sankey Level 1 data ---
    all_genes = set(gene_both_count) | set(gene_ensembl_only_count) | set(gene_cat_only_count)
    gene_rows = []
    for gn in sorted(all_genes):
        n_both = gene_both_count.get(gn, 0)
        n_ens = gene_ensembl_only_count.get(gn, 0)
        n_cat = gene_cat_only_count.get(gn, 0)
        total_seen = n_both + n_ens + n_cat

        # Classify the gene across the cohort
        if n_both == total_seen:
            cohort_status = 'always_both'
        elif n_both > 0:
            cohort_status = 'mixed'
        elif n_ens == total_seen:
            cohort_status = 'always_ensembl_only'
        elif n_cat == total_seen:
            cohort_status = 'always_cat_only'
        else:
            cohort_status = 'mixed_exclusive'

        gene_rows.append({
            'gene_name': gn,
            'n_assemblies_both': n_both,
            'n_assemblies_ensembl_only': n_ens,
            'n_assemblies_cat_only': n_cat,
            'n_assemblies_total': total_seen,
            'pct_assemblies_both': round(n_both / total_seen * 100, 2) if total_seen > 0 else 0,
            'cohort_status': cohort_status,
            'ensembl_biotype': gene_biotype_ensembl.get(gn, ''),
            'cat_biotype': gene_biotype_cat.get(gn, ''),
        })

    pd.DataFrame(gene_rows).to_csv(
        os.path.join(args.output_dir, 'sankey_level1_gene_presence.tsv'),
        sep='\t', index=False
    )
    print(f"Wrote Sankey Level 1 gene presence ({len(gene_rows)} genes)", file=sys.stderr)

    # --- Biotype enrichment for method-exclusive genes (SUPP-B) ---
    # Compare biotype distribution of exclusive vs concordant genes
    biotype_rows = []

    # Genes that are ensembl-only in >50% of assemblies where they appear
    for gn in all_genes:
        n_both = gene_both_count.get(gn, 0)
        n_ens = gene_ensembl_only_count.get(gn, 0)
        n_cat = gene_cat_only_count.get(gn, 0)
        total = n_both + n_ens + n_cat
        if total == 0:
            continue

        if n_ens / total > 0.5 and n_both == 0:
            category = 'ensembl_exclusive'
        elif n_cat / total > 0.5 and n_both == 0:
            category = 'cat_exclusive'
        elif n_both / total > 0.5:
            category = 'concordant'
        else:
            category = 'variable'

        biotype_rows.append({
            'gene_name': gn,
            'category': category,
            'ensembl_biotype': gene_biotype_ensembl.get(gn, ''),
            'cat_biotype': gene_biotype_cat.get(gn, ''),
            'n_assemblies_both': n_both,
            'n_assemblies_ensembl_only': n_ens,
            'n_assemblies_cat_only': n_cat,
        })

    biotype_df = pd.DataFrame(biotype_rows)
    biotype_df.to_csv(
        os.path.join(args.output_dir, 'biotype_enrichment_exclusive_genes.tsv'),
        sep='\t', index=False
    )

    # Summary table of biotype counts per category
    if not biotype_df.empty:
        # Use ensembl_biotype for ensembl-exclusive, cat_biotype for cat-exclusive, ensembl_biotype for concordant
        summary_rows = []
        for category in ['ensembl_exclusive', 'cat_exclusive', 'concordant', 'variable']:
            subset = biotype_df[biotype_df['category'] == category]
            bio_col = 'cat_biotype' if category == 'cat_exclusive' else 'ensembl_biotype'
            for biotype, count in subset[bio_col].value_counts().items():
                summary_rows.append({
                    'category': category,
                    'biotype': biotype,
                    'n_genes': count,
                    'pct_of_category': round(count / len(subset) * 100, 2) if len(subset) > 0 else 0,
                })

        pd.DataFrame(summary_rows).to_csv(
            os.path.join(args.output_dir, 'biotype_enrichment_summary.tsv'),
            sep='\t', index=False
        )

    print("Aggregation complete", file=sys.stderr)


if __name__ == "__main__":
    main()
