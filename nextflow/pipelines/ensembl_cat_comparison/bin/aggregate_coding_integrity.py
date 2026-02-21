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
        print("ERROR: No coding_integrity TSV files found", file=sys.stderr)
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
        n_has_both_cds = ((df['has_ensembl_cds'].astype(str) == 'True') &
                          (df['has_cat_cds'].astype(str) == 'True')).sum()

        # Classification breakdown
        for cls, count in df['classification'].value_counts().items():
            classification_counts[cls] += count

        # Per-pair concordance data for matrices
        for _, row in df.iterrows():
            has_ens = str(row.get('has_ensembl_cds', 'False')) == 'True'
            has_cat = str(row.get('has_cat_cds', 'False')) == 'True'

            if not has_ens or not has_cat:
                # Both need CDS for meaningful comparison
                start_matrix[('no_cds', 'no_cds')] += 1
                stop_matrix[('no_cds', 'no_cds')] += 1
                frame_matrix[('no_cds', 'no_cds')] += 1
                continue

            ens_start = 'present' if str(row.get('start_codon_match', 'False')) == 'True' else 'absent'
            cat_start = 'present' if str(row.get('start_codon_match', 'False')) == 'True' else 'absent'

            ens_stop = 'present' if str(row.get('stop_codon_match', 'False')) == 'True' else 'absent'
            cat_stop = 'present' if str(row.get('stop_codon_match', 'False')) == 'True' else 'absent'

            ens_frame = 'intact' if str(row.get('frameshift_detected', 'False')) == 'False' else 'disrupted'
            cat_frame = 'intact' if str(row.get('frameshift_detected', 'False')) == 'False' else 'disrupted'

            # For the concordance matrix, we classify each property as match/mismatch
            # Start codon: match means both agree on position
            start_match = str(row.get('start_codon_match', 'False')) == 'True'
            stop_match = str(row.get('stop_codon_match', 'False')) == 'True'
            no_frameshift = str(row.get('frameshift_detected', 'False')) == 'False'

            start_matrix[('match' if start_match else 'mismatch', 'start')] += 1
            stop_matrix[('match' if stop_match else 'mismatch', 'stop')] += 1
            frame_matrix[('intact' if no_frameshift else 'frameshift', 'frame')] += 1

        # Per-assembly summary
        n_start_match = (df['start_codon_match'].astype(str) == 'True').sum() if 'start_codon_match' in df.columns else 0
        n_stop_match = (df['stop_codon_match'].astype(str) == 'True').sum() if 'stop_codon_match' in df.columns else 0
        n_no_frameshift = (df['frameshift_detected'].astype(str) == 'False').sum() if 'frameshift_detected' in df.columns else 0
        n_full_match = (df['classification'] == 'Full_Match').sum() if 'classification' in df.columns else 0

        assembly_rows.append({
            'assembly_accession': accession,
            'sample_name': sample,
            'n_coding_pairs': int(n_total),
            'n_has_both_cds': int(n_has_both_cds),
            'n_start_match': int(n_start_match),
            'n_stop_match': int(n_stop_match),
            'n_frame_intact': int(n_no_frameshift),
            'n_full_match': int(n_full_match),
            'pct_start_match': round(n_start_match / n_has_both_cds * 100, 2) if n_has_both_cds > 0 else 0,
            'pct_stop_match': round(n_stop_match / n_has_both_cds * 100, 2) if n_has_both_cds > 0 else 0,
            'pct_frame_intact': round(n_no_frameshift / n_has_both_cds * 100, 2) if n_has_both_cds > 0 else 0,
            'pct_full_match': round(n_full_match / n_total * 100, 2) if n_total > 0 else 0,
        })

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
