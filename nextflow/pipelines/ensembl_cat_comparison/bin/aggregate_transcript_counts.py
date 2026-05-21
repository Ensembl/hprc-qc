#!/usr/bin/env python3
"""
Aggregate transcript count data for SUPP-D scatter/Bland-Altman plots.

For each RBH gene pair across all assemblies, extracts the number of
transcripts annotated by Ensembl vs CAT. This data is already available
in the transcript concordance output (n_ensembl_transcripts, n_cat_transcripts).

Memory-efficient: streams scatter data to disk one assembly at a time,
then reservoir-samples if the file exceeds --max-scatter-points.

Produces:
  - transcript_count_scatter_data.tsv  - per-gene-pair transcript counts
  - transcript_count_per_assembly.tsv  - per-assembly summary statistics

Usage:
    aggregate_transcript_counts.py \
        --input-dir <dir> \
        --output-dir <output_dir>
"""

import argparse
import csv
import gc
import glob
import os
import random
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate transcript count concordance data")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--max-scatter-points", type=int, default=100000,
                        help="Max data points for scatter plot (sampled if exceeded)")
    return parser.parse_args()


SCATTER_FIELDS = [
    'assembly_accession', 'ensembl_gene_id', 'cat_gene_id',
    'n_ensembl_transcripts', 'n_cat_transcripts', 'transcript_concordance_rate',
]


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*_transcript_concordance.tsv")))
    if not tsv_files:
        tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*.tsv")))
    if not tsv_files:
        print("ERROR: No transcript_concordance TSV files found in " + args.input_dir, file=sys.stderr)
        try:
            print(f"  Directory contents: {os.listdir(args.input_dir)[:20]}", file=sys.stderr)
        except Exception:
            pass
        sys.exit(1)

    print(f"Found {len(tsv_files)} transcript concordance files", file=sys.stderr)

    # --- Pass 1: per-assembly summaries + stream scatter data to temp file ---
    assembly_summaries = []
    scatter_tmp = os.path.join(args.output_dir, '_scatter_all.tmp.tsv')
    total_scatter_rows = 0

    with open(scatter_tmp, 'w', newline='') as scatter_fh:
        writer = csv.DictWriter(scatter_fh, fieldnames=SCATTER_FIELDS, delimiter='\t')
        writer.writeheader()

        for i, f in enumerate(tsv_files):
            if (i + 1) % 50 == 0:
                print(f"  Processing file {i + 1}/{len(tsv_files)}", file=sys.stderr)

            try:
                df = pd.read_csv(f, sep='\t', engine='python')
            except Exception as e:
                print(f"WARNING: Could not read {f}: {e}", file=sys.stderr)
                continue

            if df.empty:
                continue

            accession = df['assembly_accession'].iloc[0] if 'assembly_accession' in df.columns else os.path.basename(f).split('_')[0]
            sample = df['sample_name'].iloc[0] if 'sample_name' in df.columns else ''

            if 'n_ensembl_transcripts' not in df.columns or 'n_cat_transcripts' not in df.columns:
                print(f"WARNING: {f} missing transcript count columns", file=sys.stderr)
                del df
                gc.collect()
                continue

            ens_counts = df['n_ensembl_transcripts'].astype(int)
            cat_counts = df['n_cat_transcripts'].astype(int)

            # Per-assembly summary (lightweight dict, not a DataFrame)
            assembly_summaries.append({
                'assembly_accession': accession,
                'sample_name': sample,
                'n_rbh_pairs': len(df),
                'ens_tx_mean': round(ens_counts.mean(), 2),
                'ens_tx_median': float(ens_counts.median()),
                'cat_tx_mean': round(cat_counts.mean(), 2),
                'cat_tx_median': float(cat_counts.median()),
                'mean_difference': round((ens_counts - cat_counts).mean(), 2),
                'correlation': round(ens_counts.corr(cat_counts), 4) if len(df) > 1 else 0,
                'n_exact_count_match': int((ens_counts == cat_counts).sum()),
                'pct_exact_count_match': round((ens_counts == cat_counts).sum() / len(df) * 100, 2) if len(df) > 0 else 0,
            })

            # Stream scatter rows to disk
            for _, row in df[['ensembl_gene_id', 'cat_gene_id',
                              'n_ensembl_transcripts', 'n_cat_transcripts',
                              'transcript_concordance_rate']].iterrows():
                writer.writerow({
                    'assembly_accession': accession,
                    'ensembl_gene_id': row['ensembl_gene_id'],
                    'cat_gene_id': row['cat_gene_id'],
                    'n_ensembl_transcripts': int(row['n_ensembl_transcripts']),
                    'n_cat_transcripts': int(row['n_cat_transcripts']),
                    'transcript_concordance_rate': float(row.get('transcript_concordance_rate', 0)),
                })
                total_scatter_rows += 1

            del df  # free memory immediately
            gc.collect()

    # Write per-assembly summary
    pd.DataFrame(assembly_summaries).to_csv(
        os.path.join(args.output_dir, 'transcript_count_per_assembly.tsv'),
        sep='\t', index=False
    )

    # --- Pass 2: reservoir sample scatter data if too large ---
    scatter_out = os.path.join(args.output_dir, 'transcript_count_scatter_data.tsv')

    if total_scatter_rows <= args.max_scatter_points:
        # Small enough — just rename
        os.rename(scatter_tmp, scatter_out)
        n_final = total_scatter_rows
    else:
        # Reservoir sample: read line-by-line, keep at most max_scatter_points
        print(f"Reservoir sampling {args.max_scatter_points} from {total_scatter_rows} scatter points",
              file=sys.stderr)
        random.seed(42)
        reservoir = []
        with open(scatter_tmp) as fh:
            header = fh.readline()
            for i, line in enumerate(fh):
                if i < args.max_scatter_points:
                    reservoir.append(line)
                else:
                    j = random.randint(0, i)
                    if j < args.max_scatter_points:
                        reservoir[j] = line

        with open(scatter_out, 'w') as out:
            out.write(header)
            for line in reservoir:
                out.write(line)

        n_final = len(reservoir)
        os.remove(scatter_tmp)

    print(f"Wrote transcript count data for {len(assembly_summaries)} assemblies "
          f"({n_final} scatter points)", file=sys.stderr)


if __name__ == "__main__":
    main()
