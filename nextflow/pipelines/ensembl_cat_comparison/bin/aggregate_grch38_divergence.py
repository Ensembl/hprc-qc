#!/usr/bin/env python3
"""
Aggregate per-assembly GRCh38 divergence results into cohort-level summaries for MAIN-2.

Produces:
  - grch38_divergence_per_assembly.tsv          - per-assembly category counts + median deltas
  - grch38_divergence_cross_tab.tsv             - cohort-level category counts
  - grch38_divergence_by_biotype.tsv            - category counts + median deltas by biotype
  - grch38_divergence_delta_distributions.tsv   - per-biotype delta distribution stats (for R plots)
  - grch38_divergence_per_assembly_delta_summary.tsv - per-assembly x biotype median deltas (for heatmaps)

Usage:
    aggregate_grch38_divergence.py \\
        --input-dir <dir> \\
        --output-dir <output_dir>
"""

import argparse
import gc
import glob
import math
import os
import sys
from collections import Counter, defaultdict

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate GRCh38 divergence data")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


NUMERIC_DELTA_COLS = [
    'ens_delta_n_exons',
    'ens_delta_transcript_length',
    'ens_delta_cds_length',
    'cat_delta_n_exons',
    'cat_delta_transcript_length',
    'cat_delta_cds_length',
]

# Columns to summarise per assembly in the per-assembly table
ASSEMBLY_MEDIAN_COLS = [
    ('ens_delta_n_exons', 'median_ens_delta_n_exons'),
    ('cat_delta_n_exons', 'median_cat_delta_n_exons'),
    ('ens_delta_transcript_length', 'median_ens_delta_transcript_length'),
    ('cat_delta_transcript_length', 'median_cat_delta_transcript_length'),
    ('ens_delta_cds_length', 'median_ens_delta_cds_length'),
    ('cat_delta_cds_length', 'median_cat_delta_cds_length'),
]


def _median(vals):
    """Median of a list of numeric values (ignoring NaN/None)."""
    clean = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not clean:
        return None
    clean.sort()
    n = len(clean)
    mid = n // 2
    if n % 2 == 0:
        return (clean[mid - 1] + clean[mid]) / 2
    return clean[mid]


def _mean(vals):
    clean = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not clean:
        return None
    return sum(clean) / len(clean)


def _sd(vals):
    clean = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if len(clean) < 2:
        return None
    m = sum(clean) / len(clean)
    return math.sqrt(sum((x - m) ** 2 for x in clean) / (len(clean) - 1))


def _percentile(vals, p):
    clean = sorted(v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v)))
    if not clean:
        return None
    idx = (len(clean) - 1) * p / 100
    lo = int(idx)
    hi = lo + 1
    if hi >= len(clean):
        return clean[-1]
    frac = idx - lo
    return clean[lo] + frac * (clean[hi] - clean[lo])


def _fmt(val, decimals=4):
    if val is None:
        return ''
    if isinstance(val, float):
        return f'{val:.{decimals}f}'
    return str(val)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*_grch38_divergence.tsv")))
    if not tsv_files:
        tsv_files = sorted(glob.glob(os.path.join(args.input_dir, "*.tsv")))
    if not tsv_files:
        print("ERROR: No grch38_divergence TSV files found in " + args.input_dir, file=sys.stderr)
        try:
            print(f"  Directory contents: {os.listdir(args.input_dir)[:20]}", file=sys.stderr)
        except Exception:
            pass
        sys.exit(1)

    print(f"Found {len(tsv_files)} divergence files", file=sys.stderr)

    # -----------------------------------------------------------------------
    # Accumulators
    # -----------------------------------------------------------------------

    assembly_rows = []
    all_category_counts = Counter()

    # biotype -> divergence_category -> count
    biotype_category_counts = defaultdict(Counter)

    # biotype -> delta_col -> [values across all assemblies]
    biotype_delta_vals = defaultdict(lambda: defaultdict(list))

    # biotype -> source ('ensembl'|'cat') -> cds_change_type -> count
    biotype_cds_change_counts = defaultdict(lambda: defaultdict(Counter))

    # (assembly_accession, sample_name, biotype) -> delta_col -> [values]
    asm_biotype_delta_vals = defaultdict(lambda: defaultdict(list))

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

        # Convert numeric delta columns
        for col in NUMERIC_DELTA_COLS:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        cats = df['divergence_category'].value_counts()
        n_both_ref = int(cats.get('both_agree_reference', 0))
        n_both_div = int(cats.get('both_agree_diverged', 0))
        n_ens_spec = int(cats.get('ensembl_specific_divergence', 0))
        n_cat_spec = int(cats.get('cat_specific_divergence', 0))

        for cat, cnt in cats.items():
            all_category_counts[cat] += int(cnt)

        # Coding/non-coding split
        if 'is_coding' in df.columns:
            is_coding_mask = df['is_coding'].str.lower() == 'true'
            n_coding = int(is_coding_mask.sum())
            n_noncoding = int((~is_coding_mask).sum())
        else:
            n_coding = n_noncoding = None

        # Frameshift counts
        n_ens_fs = n_cat_fs = None
        if 'ens_cds_frameshift' in df.columns:
            n_ens_fs = int((df['ens_cds_frameshift'].str.lower() == 'true').sum())
        if 'cat_cds_frameshift' in df.columns:
            n_cat_fs = int((df['cat_cds_frameshift'].str.lower() == 'true').sum())

        # Per-assembly medians
        asm_row = {
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
            'n_coding_genes_assessed': n_coding,
            'n_noncoding_genes_assessed': n_noncoding,
            'n_ens_frameshift': n_ens_fs,
            'n_cat_frameshift': n_cat_fs,
        }

        for src_col, out_col in ASSEMBLY_MEDIAN_COLS:
            if src_col in df.columns:
                asm_row[out_col] = round(df[src_col].median(), 4) if not df[src_col].isna().all() else None
            else:
                asm_row[out_col] = None

        assembly_rows.append(asm_row)

        # Biotype stratification
        if 'ref_biotype' in df.columns:
            for biotype, grp in df.groupby('ref_biotype'):
                biotype = str(biotype)
                for cat, cnt in grp['divergence_category'].value_counts().items():
                    biotype_category_counts[biotype][cat] += int(cnt)

                for col in NUMERIC_DELTA_COLS:
                    if col in grp.columns:
                        vals = grp[col].dropna().tolist()
                        biotype_delta_vals[biotype][col].extend(vals)
                        asm_biotype_delta_vals[(accession, sample, biotype)][col].extend(vals)

                # CDS change type tracking
                for src, col in [('ensembl', 'ens_cds_change_type'), ('cat', 'cat_cds_change_type')]:
                    if col in grp.columns:
                        for chg_type, cnt in grp[col].value_counts().items():
                            if pd.notna(chg_type):
                                biotype_cds_change_counts[biotype][src][str(chg_type)] += int(cnt)

        del df
        gc.collect()

    # -----------------------------------------------------------------------
    # Write per-assembly summary
    # -----------------------------------------------------------------------
    asm_cols = [
        'assembly_accession', 'sample_name', 'n_genes_assessed',
        'n_both_agree_reference', 'n_both_agree_diverged',
        'n_ensembl_specific', 'n_cat_specific',
        'pct_both_agree_reference', 'pct_both_agree_diverged',
        'pct_ensembl_specific', 'pct_cat_specific',
        'n_coding_genes_assessed', 'n_noncoding_genes_assessed',
        'n_ens_frameshift', 'n_cat_frameshift',
        'median_ens_delta_n_exons', 'median_cat_delta_n_exons',
        'median_ens_delta_transcript_length', 'median_cat_delta_transcript_length',
        'median_ens_delta_cds_length', 'median_cat_delta_cds_length',
    ]
    pd.DataFrame(assembly_rows)[asm_cols].to_csv(
        os.path.join(args.output_dir, 'grch38_divergence_per_assembly.tsv'),
        sep='\t', index=False,
    )

    # -----------------------------------------------------------------------
    # Write cohort-level cross-tabulation
    # -----------------------------------------------------------------------
    total = sum(all_category_counts.values())
    cross_tab = []
    for cat in ['both_agree_reference', 'both_agree_diverged',
                'ensembl_specific_divergence', 'cat_specific_divergence',
                'insufficient_data']:
        cnt = all_category_counts.get(cat, 0)
        cross_tab.append({
            'divergence_category': cat,
            'total_gene_assembly_instances': cnt,
            'pct_of_total': round(cnt / total * 100, 2) if total > 0 else 0,
        })
    pd.DataFrame(cross_tab).to_csv(
        os.path.join(args.output_dir, 'grch38_divergence_cross_tab.tsv'),
        sep='\t', index=False,
    )

    # -----------------------------------------------------------------------
    # Write biotype breakdown with median deltas
    # -----------------------------------------------------------------------
    bio_rows = []
    for biotype in sorted(biotype_category_counts):
        for cat in ['both_agree_reference', 'both_agree_diverged',
                    'ensembl_specific_divergence', 'cat_specific_divergence',
                    'insufficient_data']:
            cnt = biotype_category_counts[biotype].get(cat, 0)
            if cnt == 0:
                continue
            row = {
                'biotype': biotype,
                'divergence_category': cat,
                'count': cnt,
            }
            for col in NUMERIC_DELTA_COLS:
                vals = biotype_delta_vals[biotype].get(col, [])
                row[f'median_{col}'] = round(_median(vals), 4) if vals else None
            bio_rows.append(row)

    if bio_rows:
        bio_cols = (['biotype', 'divergence_category', 'count']
                    + [f'median_{c}' for c in NUMERIC_DELTA_COLS])
        pd.DataFrame(bio_rows)[bio_cols].to_csv(
            os.path.join(args.output_dir, 'grch38_divergence_by_biotype.tsv'),
            sep='\t', index=False,
        )

    # -----------------------------------------------------------------------
    # Write delta distributions (for violin/density plots in R)
    # -----------------------------------------------------------------------
    dist_rows = []
    for biotype in sorted(biotype_delta_vals):
        for col in NUMERIC_DELTA_COLS:
            vals = biotype_delta_vals[biotype].get(col, [])
            if not vals:
                continue
            # Parse source from column name (ens_ or cat_)
            source = 'ensembl' if col.startswith('ens_') else 'cat'
            metric = col[len('ens_'):] if col.startswith('ens_') else col[len('cat_'):]
            dist_rows.append({
                'biotype': biotype,
                'metric': metric,
                'source': source,
                'n': len(vals),
                'mean': _fmt(_mean(vals)),
                'median': _fmt(_median(vals)),
                'sd': _fmt(_sd(vals)),
                'p5': _fmt(_percentile(vals, 5)),
                'p25': _fmt(_percentile(vals, 25)),
                'p75': _fmt(_percentile(vals, 75)),
                'p95': _fmt(_percentile(vals, 95)),
            })

    if dist_rows:
        dist_cols = ['biotype', 'metric', 'source', 'n',
                     'mean', 'median', 'sd', 'p5', 'p25', 'p75', 'p95']
        pd.DataFrame(dist_rows)[dist_cols].to_csv(
            os.path.join(args.output_dir, 'grch38_divergence_delta_distributions.tsv'),
            sep='\t', index=False,
        )

    # -----------------------------------------------------------------------
    # Write CDS change type breakdown by biotype (for protein-coding analysis)
    # -----------------------------------------------------------------------
    cds_change_rows = []
    for biotype in sorted(biotype_cds_change_counts):
        for source in ['ensembl', 'cat']:
            for cds_type, cnt in sorted(biotype_cds_change_counts[biotype][source].items()):
                cds_change_rows.append({
                    'biotype': biotype,
                    'source': source,
                    'cds_change_type': cds_type,
                    'count': cnt,
                })

    if cds_change_rows:
        pd.DataFrame(cds_change_rows).to_csv(
            os.path.join(args.output_dir, 'grch38_cds_change_type_by_biotype.tsv'),
            sep='\t', index=False,
        )

    # -----------------------------------------------------------------------
    # Write per-assembly x biotype delta summary (for heatmaps)
    # -----------------------------------------------------------------------
    asm_bio_rows = []
    for (accession, sample, biotype) in sorted(asm_biotype_delta_vals):
        col_data = asm_biotype_delta_vals[(accession, sample, biotype)]
        n_genes = max((len(v) for v in col_data.values()), default=0)
        row = {
            'assembly_accession': accession,
            'sample_name': sample,
            'biotype': biotype,
            'n_genes': n_genes,
        }
        for col in NUMERIC_DELTA_COLS:
            vals = col_data.get(col, [])
            row[f'median_{col}'] = round(_median(vals), 4) if vals else None
        asm_bio_rows.append(row)

    if asm_bio_rows:
        asm_bio_cols = (['assembly_accession', 'sample_name', 'biotype', 'n_genes']
                        + [f'median_{c}' for c in NUMERIC_DELTA_COLS])
        pd.DataFrame(asm_bio_rows)[asm_bio_cols].to_csv(
            os.path.join(args.output_dir, 'grch38_divergence_per_assembly_delta_summary.tsv'),
            sep='\t', index=False,
        )

    print(
        f"Wrote divergence aggregation for {len(assembly_rows)} assemblies",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
