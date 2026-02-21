#!/usr/bin/env python3
"""
Aggregate per-assembly QC metrics into Sankey/funnel flow counts for MAIN-1.

Levels:
  1. Gene presence (both | ensembl_only | cat_only)
  2. RBH status (rbh_found | no_rbh) - subset of 'both'
  3. Transcript concordance (full_match | partial | none) - subset of 'rbh_found'
  4. CDS integrity (intact | partial_disruption | disrupted) - subset of coding RBH

Produces:
  - sankey_flow_counts.tsv           - aggregate flow counts per level
  - sankey_per_assembly_flows.tsv    - per-assembly flow counts for distributions

Usage:
    aggregate_sankey_flows.py \
        --gene-presence-dir <dir> \
        --rbh-dir <dir> \
        --transcript-concordance-dir <dir> \
        --coding-integrity-dir <dir> \
        --output-dir <output_dir> \
        [--rbh-coverage-threshold 0.95]
"""

import argparse
import glob
import os
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate Sankey flow counts across assemblies")
    parser.add_argument("--gene-presence-dir", required=True)
    parser.add_argument("--rbh-dir", required=True)
    parser.add_argument("--transcript-concordance-dir", required=True)
    parser.add_argument("--coding-integrity-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--rbh-coverage-threshold", type=float, default=0.95,
                        help="Min reciprocal coverage fraction for RBH to pass Level 2 (default: 0.95)")
    return parser.parse_args()


def load_tsvs(directory, pattern, label):
    """Load all matching TSVs from a directory."""
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        print(f"WARNING: No {label} files found in {directory}", file=sys.stderr)
        return {}

    dfs = {}
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty:
                continue
            accession = df['assembly_accession'].iloc[0] if 'assembly_accession' in df.columns else os.path.basename(f).split('_')[0]
            dfs[accession] = df
        except Exception as e:
            print(f"WARNING: Could not read {f}: {e}", file=sys.stderr)
    return dfs


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    threshold = args.rbh_coverage_threshold

    # Load all per-assembly data
    presence_dfs = load_tsvs(args.gene_presence_dir, "*_gene_presence.tsv", "gene presence")
    rbh_dfs = load_tsvs(args.rbh_dir, "*.gene_pairs_rbh.tsv", "RBH")
    tc_dfs = load_tsvs(args.transcript_concordance_dir, "*_transcript_concordance.tsv", "transcript concordance")
    ci_dfs = load_tsvs(args.coding_integrity_dir, "*_coding_integrity.tsv", "coding integrity")

    # Get common assemblies
    all_accessions = set(presence_dfs.keys()) & set(rbh_dfs.keys())
    print(f"Processing {len(all_accessions)} assemblies with both presence and RBH data", file=sys.stderr)

    per_assembly_rows = []

    for accession in sorted(all_accessions):
        pres = presence_dfs[accession]
        rbh = rbh_dfs[accession]

        # Level 1: Gene presence
        n_both = ((pres['present_in_ensembl'] == 'True') & (pres['present_in_cat'] == 'True')).sum() if 'present_in_ensembl' in pres.columns else 0
        # Handle string/bool
        if n_both == 0 and 'present_in_ensembl' in pres.columns:
            n_both = ((pres['present_in_ensembl'].astype(str) == 'True') & (pres['present_in_cat'].astype(str) == 'True')).sum()
        n_ens_only = ((pres['present_in_ensembl'].astype(str) == 'True') & (pres['present_in_cat'].astype(str) == 'False')).sum()
        n_cat_only = ((pres['present_in_ensembl'].astype(str) == 'False') & (pres['present_in_cat'].astype(str) == 'True')).sum()
        n_total = n_both + n_ens_only + n_cat_only

        # Level 2: RBH status
        n_rbh_total = len(rbh)
        if 'frac_ensembl_covered' in rbh.columns and 'frac_cat_covered' in rbh.columns:
            rbh_pass = rbh[
                (rbh['frac_ensembl_covered'].astype(float) >= threshold) &
                (rbh['frac_cat_covered'].astype(float) >= threshold)
            ]
            n_rbh_pass = len(rbh_pass)
        else:
            n_rbh_pass = n_rbh_total

        n_rbh_fail = n_rbh_total - n_rbh_pass

        # Level 3: Transcript concordance (subset of RBH)
        tc = tc_dfs.get(accession)
        n_tx_full = 0
        n_tx_partial = 0
        n_tx_none = 0
        if tc is not None and not tc.empty:
            if 'transcript_concordance_rate' in tc.columns:
                rates = tc['transcript_concordance_rate'].astype(float)
                n_tx_full = (rates >= 1.0).sum()
                n_tx_partial = ((rates > 0) & (rates < 1.0)).sum()
                n_tx_none = (rates == 0).sum()
            n_tx_total = n_tx_full + n_tx_partial + n_tx_none
        else:
            n_tx_total = 0

        # Level 4: CDS integrity (subset of coding RBH)
        ci = ci_dfs.get(accession)
        n_cds_intact = 0
        n_cds_partial = 0
        n_cds_disrupted = 0
        n_cds_total = 0
        if ci is not None and not ci.empty:
            n_cds_total = len(ci)
            if 'classification' in ci.columns:
                n_cds_intact = (ci['classification'] == 'Full_Match').sum()
                disrupted_classes = {'Internal_Frameshift', 'Complex_Mismatch', 'No_CDS'}
                n_cds_disrupted = ci['classification'].isin(disrupted_classes).sum()
                n_cds_partial = n_cds_total - n_cds_intact - n_cds_disrupted

        sample = pres['sample_name'].iloc[0] if 'sample_name' in pres.columns else ''
        per_assembly_rows.append({
            'assembly_accession': accession,
            'sample_name': sample,
            # Level 1
            'l1_total': int(n_total),
            'l1_both': int(n_both),
            'l1_ensembl_only': int(n_ens_only),
            'l1_cat_only': int(n_cat_only),
            'l1_pct_both': round(n_both / n_total * 100, 2) if n_total > 0 else 0,
            # Level 2
            'l2_rbh_total': int(n_rbh_total),
            'l2_rbh_pass': int(n_rbh_pass),
            'l2_rbh_fail': int(n_rbh_fail),
            'l2_pct_pass': round(n_rbh_pass / n_rbh_total * 100, 2) if n_rbh_total > 0 else 0,
            # Level 3
            'l3_tx_total': int(n_tx_total),
            'l3_tx_full': int(n_tx_full),
            'l3_tx_partial': int(n_tx_partial),
            'l3_tx_none': int(n_tx_none),
            'l3_pct_full': round(n_tx_full / n_tx_total * 100, 2) if n_tx_total > 0 else 0,
            # Level 4
            'l4_cds_total': int(n_cds_total),
            'l4_cds_intact': int(n_cds_intact),
            'l4_cds_partial': int(n_cds_partial),
            'l4_cds_disrupted': int(n_cds_disrupted),
            'l4_pct_intact': round(n_cds_intact / n_cds_total * 100, 2) if n_cds_total > 0 else 0,
        })

    per_asm_df = pd.DataFrame(per_assembly_rows)
    per_asm_df.to_csv(
        os.path.join(args.output_dir, 'sankey_per_assembly_flows.tsv'),
        sep='\t', index=False
    )

    # Aggregate across assemblies (medians and totals)
    if not per_asm_df.empty:
        agg = {
            'n_assemblies': len(per_asm_df),
            'rbh_coverage_threshold': threshold,
        }
        for col in per_asm_df.columns:
            if col.startswith('l') and col[1].isdigit():
                agg[f'{col}_median'] = round(per_asm_df[col].median(), 2)
                agg[f'{col}_mean'] = round(per_asm_df[col].mean(), 2)
                agg[f'{col}_q25'] = round(per_asm_df[col].quantile(0.25), 2)
                agg[f'{col}_q75'] = round(per_asm_df[col].quantile(0.75), 2)
                agg[f'{col}_sum'] = int(per_asm_df[col].sum())

        pd.DataFrame([agg]).to_csv(
            os.path.join(args.output_dir, 'sankey_flow_counts.tsv'),
            sep='\t', index=False
        )

    print(f"Wrote Sankey flow data for {len(per_assembly_rows)} assemblies", file=sys.stderr)


if __name__ == "__main__":
    main()
