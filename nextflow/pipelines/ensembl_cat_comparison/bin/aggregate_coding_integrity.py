#!/usr/bin/env python3
"""
Aggregate per-assembly coding integrity files into concordance matrices for SUPP-C.

Produces:
  - cds_concordance_matrices.tsv     - 3x3 heatmap data (start/stop/frame: Ensembl vs CAT)
  - coding_integrity_per_assembly.tsv - per-assembly CDS summary stats

Usage:
    aggregate_coding_integrity.py \
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
    parser = argparse.ArgumentParser(description="Aggregate coding integrity into concordance matrices")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*_coding_integrity.tsv")))
    if not tsv_files:
        tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*.tsv")))
    if not tsv_files:
        print("ERROR: No coding_integrity TSV files found in " + args.input_dir, file=sys.stderr)
        try:
            print(f"  Directory contents: {os.listdir(args.input_dir)[:20]}", file=sys.stderr)
        except Exception:
            pass
        sys.exit(1)

    print(f"Found {len(tsv_files)} coding integrity files", file=sys.stderr)

    # Concordance matrix counters
    # For each CDS property, we count (ensembl_call, cat_call) combinations
    start_matrix = Counter()  # (ens_start_match_status, cat_start_match_status) -> count
    stop_matrix = Counter()
    frame_matrix = Counter()

    # We track per-assembly summaries
    assembly_rows = []

    # Classification counts across all assemblies
    classification_counts = Counter()

    for i, f in enumerate(tsv_files):
        if (i + 1) % 50 == 0:
            print(f"  Processing file {i + 1}/{len(tsv_files)}", file=sys.stderr)

        try:
            df = pd.read_csv(f, sep='\t', dtype=str, engine='python')
        except Exception as e:
            print(f"WARNING: Could not read {f}: {e}", file=sys.stderr)
            continue

        if df.empty:
            continue

        accession = df['assembly_accession'].iloc[0]
        sample = df['sample_name'].iloc[0]

        n_total = len(df)
        has_ens_cds = df['has_ensembl_cds'].astype(str) == 'True'
        has_cat_cds = df['has_cat_cds'].astype(str) == 'True'
        has_both = has_ens_cds & has_cat_cds
        n_has_both_cds = int(has_both.sum())
        n_no_cds = int((~has_both).sum())

        # Classification breakdown (vectorised)
        for cls, count in df['classification'].value_counts().items():
            classification_counts[cls] += count

        # Concordance matrices — vectorised instead of iterrows
        start_matrix[('no_cds', 'no_cds')] += n_no_cds
        stop_matrix[('no_cds', 'no_cds')] += n_no_cds
        frame_matrix[('no_cds', 'no_cds')] += n_no_cds

        if n_has_both_cds > 0:
            cds_df = df.loc[has_both]

            start_match_mask = cds_df['start_codon_match'].astype(str) == 'True'
            start_matrix[('match', 'start')] += int(start_match_mask.sum())
            start_matrix[('mismatch', 'start')] += int((~start_match_mask).sum())

            stop_match_mask = cds_df['stop_codon_match'].astype(str) == 'True'
            stop_matrix[('match', 'stop')] += int(stop_match_mask.sum())
            stop_matrix[('mismatch', 'stop')] += int((~stop_match_mask).sum())

            no_fs_mask = cds_df['frameshift_detected'].astype(str) == 'False'
            frame_matrix[('intact', 'frame')] += int(no_fs_mask.sum())
            frame_matrix[('frameshift', 'frame')] += int((~no_fs_mask).sum())

        # Per-assembly summary (vectorised)
        n_start_match = int((df['start_codon_match'].astype(str) == 'True').sum()) if 'start_codon_match' in df.columns else 0
        n_stop_match = int((df['stop_codon_match'].astype(str) == 'True').sum()) if 'stop_codon_match' in df.columns else 0
        n_no_frameshift = int((df['frameshift_detected'].astype(str) == 'False').sum()) if 'frameshift_detected' in df.columns else 0
        n_full_match = int((df['classification'] == 'Full_Match').sum()) if 'classification' in df.columns else 0

        assembly_rows.append({
            'assembly_accession': accession,
            'sample_name': sample,
            'n_coding_pairs': int(n_total),
            'n_has_both_cds': n_has_both_cds,
            'n_start_match': n_start_match,
            'n_stop_match': n_stop_match,
            'n_frame_intact': n_no_frameshift,
            'n_full_match': n_full_match,
            'pct_start_match': round(n_start_match / n_has_both_cds * 100, 2) if n_has_both_cds > 0 else 0,
            'pct_stop_match': round(n_stop_match / n_has_both_cds * 100, 2) if n_has_both_cds > 0 else 0,
            'pct_frame_intact': round(n_no_frameshift / n_has_both_cds * 100, 2) if n_has_both_cds > 0 else 0,
            'pct_full_match': round(n_full_match / n_total * 100, 2) if n_total > 0 else 0,
        })

        del df  # free memory

    # Write per-assembly summary
    pd.DataFrame(assembly_rows).to_csv(
        os.path.join(args.output_dir, 'coding_integrity_per_assembly.tsv'),
        sep='\t', index=False
    )

    # Write concordance matrix data
    # Reshape into a format suitable for heatmap plotting
    matrix_rows = []

    # Start codon concordance
    for (status, _), count in start_matrix.items():
        matrix_rows.append({
            'cds_property': 'start_codon',
            'status': status,
            'count': count,
        })

    # Stop codon concordance
    for (status, _), count in stop_matrix.items():
        matrix_rows.append({
            'cds_property': 'stop_codon',
            'status': status,
            'count': count,
        })

    # Frame concordance
    for (status, _), count in frame_matrix.items():
        matrix_rows.append({
            'cds_property': 'reading_frame',
            'status': status,
            'count': count,
        })

    pd.DataFrame(matrix_rows).to_csv(
        os.path.join(args.output_dir, 'cds_concordance_matrices.tsv'),
        sep='\t', index=False
    )

    # Write classification distribution
    cls_rows = [{'classification': cls, 'count': cnt}
                for cls, cnt in sorted(classification_counts.items(), key=lambda x: -x[1])]
    pd.DataFrame(cls_rows).to_csv(
        os.path.join(args.output_dir, 'cds_classification_distribution.tsv'),
        sep='\t', index=False
    )

    print(f"Wrote coding integrity aggregation for {len(assembly_rows)} assemblies", file=sys.stderr)


if __name__ == "__main__":
    main()
