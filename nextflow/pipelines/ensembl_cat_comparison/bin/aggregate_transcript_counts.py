#!/usr/bin/env python3
"""
Aggregate transcript count data for SUPP-D scatter/Bland-Altman plots.

For each RBH gene pair across all assemblies, extracts the number of
transcripts annotated by Ensembl vs CAT. This data is already available
in the transcript concordance output (n_ensembl_transcripts, n_cat_transcripts).

Produces:
  - transcript_count_scatter_data.tsv  - per-gene-pair transcript counts
  - transcript_count_per_assembly.tsv  - per-assembly summary statistics

Usage:
    aggregate_transcript_counts.py \
        --input-dir <dir> \
        --output-dir <output_dir>
"""

import argparse
import glob
import os
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate transcript count concordance data")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--max-scatter-points", type=int, default=100000,
                        help="Max data points for scatter plot (sampled if exceeded)")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*_transcript_concordance.tsv")))
    if not tsv_files:
        print("ERROR: No transcript_concordance TSV files found", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(tsv_files)} transcript concordance files", file=sys.stderr)

    all_scatter_data = []
    assembly_summaries = []

    for f in tsv_files:
        try:
            df = pd.read_csv(f, sep='\t')
        except Exception as e:
            print(f"WARNING: Could not read {f}: {e}", file=sys.stderr)
            continue

        if df.empty:
            continue

        accession = df['assembly_accession'].iloc[0]
        sample = df['sample_name'].iloc[0]

        if 'n_ensembl_transcripts' not in df.columns or 'n_cat_transcripts' not in df.columns:
            print(f"WARNING: {f} missing transcript count columns", file=sys.stderr)
            continue

        ens_counts = df['n_ensembl_transcripts'].astype(int)
        cat_counts = df['n_cat_transcripts'].astype(int)

        # Per-assembly summary
        assembly_summaries.append({
            'assembly_accession': accession,
            'sample_name': sample,
            'n_rbh_pairs': len(df),
            'ens_tx_mean': round(ens_counts.mean(), 2),
            'ens_tx_median': ens_counts.median(),
            'cat_tx_mean': round(cat_counts.mean(), 2),
            'cat_tx_median': cat_counts.median(),
            'mean_difference': round((ens_counts - cat_counts).mean(), 2),
            'correlation': round(ens_counts.corr(cat_counts), 4) if len(df) > 1 else 0,
            'n_exact_count_match': int((ens_counts == cat_counts).sum()),
            'pct_exact_count_match': round((ens_counts == cat_counts).sum() / len(df) * 100, 2) if len(df) > 0 else 0,
        })

        # Collect scatter data
        for _, row in df.iterrows():
            all_scatter_data.append({
                'assembly_accession': accession,
                'ensembl_gene_id': row['ensembl_gene_id'],
                'cat_gene_id': row['cat_gene_id'],
                'n_ensembl_transcripts': int(row['n_ensembl_transcripts']),
                'n_cat_transcripts': int(row['n_cat_transcripts']),
                'transcript_concordance_rate': float(row.get('transcript_concordance_rate', 0)),
            })

    # Write per-assembly summary
    pd.DataFrame(assembly_summaries).to_csv(
        os.path.join(args.output_dir, 'transcript_count_per_assembly.tsv'),
        sep='\t', index=False
    )

    # Write scatter data (sample if too large)
    scatter_df = pd.DataFrame(all_scatter_data)
    if len(scatter_df) > args.max_scatter_points:
        print(f"Sampling {args.max_scatter_points} from {len(scatter_df)} data points", file=sys.stderr)
        scatter_df = scatter_df.sample(n=args.max_scatter_points, random_state=42)

    scatter_df.to_csv(
        os.path.join(args.output_dir, 'transcript_count_scatter_data.tsv'),
        sep='\t', index=False
    )

    print(f"Wrote transcript count data for {len(assembly_summaries)} assemblies "
          f"({len(scatter_df)} scatter points)", file=sys.stderr)


if __name__ == "__main__":
    main()
