#!/usr/bin/env python3
"""
Aggregate per-gene per-assembly transcript concordance "problem" frequencies.

Inputs: one or more *_transcript_concordance.tsv files produced by
nextflow/pipelines/ensembl_cat_comparison/modules/transcript_concordance.

For each gene pair row, compute per-direction frequencies for the
non-ideal classifications ("problems"):
  - Partial_5, Partial_3, Other_Partial, Intron_Subset, Intron_Superset, No_Match

Optionally also include the non-problem categories for context
  - Exact_Match, Intron_Match

Output schema (wide, tidy by gene×assembly):
  assembly_accession, sample_name, ensembl_gene_id, cat_gene_id, ensembl_biotype,
  n_ensembl_transcripts, n_cat_transcripts,
  ens_freq_exact, ens_freq_intron_match, ens_freq_subset, ens_freq_superset,
  ens_freq_partial_5, ens_freq_partial_3, ens_freq_other_partial, ens_freq_unmatched,
  ens_problem_rate,
  cat_freq_exact, cat_freq_intron_match, cat_freq_subset, cat_freq_superset,
  cat_freq_partial_5, cat_freq_partial_3, cat_freq_other_partial, cat_freq_unmatched,
  cat_problem_rate

Usage:
    aggregate_gene_assembly_problem_freq.py \
        --input-dir <dir_with_concordance_tsvs> \
        --output <output_tsv> \
        [--include-context] [--long <long_output_tsv>]
"""

import argparse
import os
import sys
import glob
import pandas as pd


PROBLEM_COLSETS = {
    'ens': [
        ('n_ens_partial_5',     'ens_freq_partial_5'),
        ('n_ens_partial_3',     'ens_freq_partial_3'),
        ('n_ens_other_partial', 'ens_freq_other_partial'),
        ('n_ens_subset',        'ens_freq_subset'),
        ('n_ens_superset',      'ens_freq_superset'),
        ('n_ens_unmatched',     'ens_freq_unmatched'),
    ],
    'cat': [
        ('n_cat_partial_5',     'cat_freq_partial_5'),
        ('n_cat_partial_3',     'cat_freq_partial_3'),
        ('n_cat_other_partial', 'cat_freq_other_partial'),
        ('n_cat_subset',        'cat_freq_subset'),
        ('n_cat_superset',      'cat_freq_superset'),
        ('n_cat_unmatched',     'cat_freq_unmatched'),
    ],
}

CONTEXT_COLSETS = {
    'ens': [
        ('n_ens_exact',        'ens_freq_exact'),
        ('n_ens_intron_match', 'ens_freq_intron_match'),
    ],
    'cat': [
        ('n_cat_exact',        'cat_freq_exact'),
        ('n_cat_intron_match', 'cat_freq_intron_match'),
    ],
}


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate per-gene problem frequencies across assemblies")
    p.add_argument("--input-dir", required=True, help="Directory to scan for *_transcript_concordance.tsv (recursive)")
    p.add_argument("--output", required=True, help="Output TSV path")
    p.add_argument("--include-context", action="store_true",
                   help="Also include exact/intron_match frequency columns")
    p.add_argument("--long", dest="long_out", default=None,
                   help="Optional: write a tidy long-form table here with columns: "
                        "assembly_accession, sample_name, ensembl_gene_id, cat_gene_id, ensembl_biotype, "
                        "direction, classification, frequency")
    return p.parse_args()


def find_concordance_files(root: str):
    files = []
    for path in glob.glob(os.path.join(root, "**", "*_transcript_concordance.tsv"), recursive=True):
        if os.path.isfile(path):
            files.append(path)
    return sorted(files)


def ensure_numeric(df: pd.DataFrame, cols):
    for c in cols:
        if c not in df.columns:
            df[c] = 0
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)
    return df


def compute_freqs(row, count_col, denom):
    n = float(row.get(count_col, 0))
    d = float(row.get(denom, 0))
    if d <= 0:
        return 0.0
    return round(n / d, 6)


def process_file(path: str, include_context: bool) -> pd.DataFrame:
    # Columns we need; tolerate missing by filling with 0
    base_cols = [
        'assembly_accession', 'sample_name',
        'ensembl_gene_id', 'cat_gene_id', 'ensembl_biotype',
        'n_ensembl_transcripts', 'n_cat_transcripts',
    ]

    count_cols = [c for c, _ in PROBLEM_COLSETS['ens'] + PROBLEM_COLSETS['cat']]
    if include_context:
        count_cols += [c for c, _ in CONTEXT_COLSETS['ens'] + CONTEXT_COLSETS['cat']]

    usecols = base_cols + count_cols

    try:
        df = pd.read_csv(path, sep='\t', usecols=lambda c: c in usecols)
    except Exception as e:
        print(f"WARNING: Skipping {path}: {e}", file=sys.stderr)
        return pd.DataFrame()

    if df.empty:
        return df

    # Fill missing essential columns if absent in partial files
    for c in base_cols:
        if c not in df.columns:
            df[c] = '' if c in {'assembly_accession', 'sample_name', 'ensembl_biotype'} else 0

    ensure_numeric(df, ['n_ensembl_transcripts', 'n_cat_transcripts'] + count_cols)

    out_rows = []
    for _, r in df.iterrows():
        rec = {
            'assembly_accession': r['assembly_accession'],
            'sample_name': r['sample_name'],
            'ensembl_gene_id': r['ensembl_gene_id'],
            'cat_gene_id': r.get('cat_gene_id', ''),
            'ensembl_biotype': r.get('ensembl_biotype', ''),
            'n_ensembl_transcripts': int(r.get('n_ensembl_transcripts', 0)),
            'n_cat_transcripts': int(r.get('n_cat_transcripts', 0)),
        }

        # Ensembl→CAT problem frequencies
        for c, out_name in PROBLEM_COLSETS['ens']:
            rec[out_name] = compute_freqs(r, c, 'n_ensembl_transcripts')

        # CAT→Ensembl problem frequencies
        for c, out_name in PROBLEM_COLSETS['cat']:
            rec[out_name] = compute_freqs(r, c, 'n_cat_transcripts')

        # Optional context (non-problem) frequencies
        if include_context:
            for c, out_name in CONTEXT_COLSETS['ens']:
                rec[out_name] = compute_freqs(r, c, 'n_ensembl_transcripts')
            for c, out_name in CONTEXT_COLSETS['cat']:
                rec[out_name] = compute_freqs(r, c, 'n_cat_transcripts')

        # Summed problem rate per direction (bounds-safe)
        rec['ens_problem_rate'] = round(
            sum(rec[k] for _, k in PROBLEM_COLSETS['ens']), 6
        )
        rec['cat_problem_rate'] = round(
            sum(rec[k] for _, k in PROBLEM_COLSETS['cat']), 6
        )

        out_rows.append(rec)

    return pd.DataFrame(out_rows)


def main():
    args = parse_args()
    files = find_concordance_files(args.input_dir)
    if not files:
        print(f"ERROR: No *_transcript_concordance.tsv found under {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    all_rows = []
    for i, f in enumerate(files, 1):
        if i == 1 or i % 50 == 0:
            print(f"  Processing {i}/{len(files)}: {os.path.basename(f)}", file=sys.stderr)
        df = process_file(f, args.include_context)
        if not df.empty:
            all_rows.append(df)

    if not all_rows:
        print("ERROR: No data rows found (files may be empty).", file=sys.stderr)
        sys.exit(1)

    out_df = pd.concat(all_rows, ignore_index=True)

    # Wide table (per gene×assembly)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(out_df)} rows to {args.output}", file=sys.stderr)

    # Optional long-form table
    if args.long_out:
        long_rows = []
        problem_cols = [c for _, c in PROBLEM_COLSETS['ens']] + [c for _, c in PROBLEM_COLSETS['cat']]
        context_cols = []
        if args.include_context:
            context_cols = [c for _, c in CONTEXT_COLSETS['ens']] + [c for _, c in CONTEXT_COLSETS['cat']]

        def add_rows(row, prefix, labelmap):
            for col, pretty in labelmap:
                long_rows.append({
                    'assembly_accession': row['assembly_accession'],
                    'sample_name': row['sample_name'],
                    'ensembl_gene_id': row['ensembl_gene_id'],
                    'cat_gene_id': row.get('cat_gene_id', ''),
                    'ensembl_biotype': row.get('ensembl_biotype', ''),
                    'direction': 'Ensembl_to_CAT' if prefix == 'ens' else 'CAT_to_Ensembl',
                    'classification': pretty.replace('ens_freq_', '').replace('cat_freq_', ''),
                    'frequency': row[pretty]
                })

        # Build mapping of raw count label → output frequency column name used above
        ens_map = PROBLEM_COLSETS['ens'] + (CONTEXT_COLSETS['ens'] if args.include_context else [])
        cat_map = PROBLEM_COLSETS['cat'] + (CONTEXT_COLSETS['cat'] if args.include_context else [])

        for _, row in out_df.iterrows():
            add_rows(row, 'ens', ens_map)
            add_rows(row, 'cat', cat_map)

        long_df = pd.DataFrame(long_rows)
        os.makedirs(os.path.dirname(os.path.abspath(args.long_out)), exist_ok=True)
        long_df.to_csv(args.long_out, sep='\t', index=False)
        print(f"Wrote {len(long_df)} long-form rows to {args.long_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
