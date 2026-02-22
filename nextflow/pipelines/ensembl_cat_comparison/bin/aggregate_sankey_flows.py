#!/usr/bin/env python3
"""
Aggregate per-assembly QC metrics into Sankey/funnel flow counts for MAIN-1.

Levels:
  1. Gene presence (both | ensembl_only | cat_only)
  2. RBH status (rbh_found | no_rbh) - subset of 'both'
  3. Transcript concordance (full_match | partial | none) - subset of 'rbh_found'
  4. CDS integrity (intact | partial_disruption | disrupted) - subset of coding RBH

Memory-efficient: indexes files first, then loads one assembly at a time.

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
import csv
import gc
import glob
import os
import re
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


# Regex for GCA/GCF accession with version: GCA_018466835.2
ACCESSION_RE = re.compile(r'(GC[AF]_\d+\.\d+)')


def extract_accession_from_filename(filepath):
    """Extract assembly accession from filename."""
    basename = os.path.basename(filepath)
    m = ACCESSION_RE.search(basename)
    if m:
        return m.group(1)
    return None


def extract_accession_from_header(filepath):
    """Read just the first data row to get assembly_accession, without loading full file."""
    try:
        with open(filepath) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            row = next(reader, None)
            if row and 'assembly_accession' in row:
                return row['assembly_accession']
    except Exception:
        pass
    return None


def index_files(directory, patterns, label):
    """Build accession -> filepath mapping without loading data.

    Tries glob patterns in order; falls back to *.tsv.
    Extracts accession from filename first, then from file header if needed.
    """
    files = []
    for pattern in patterns:
        files = sorted(glob.glob(os.path.join(directory, pattern)))
        if files:
            print(f"Found {len(files)} {label} files with pattern '{pattern}'", file=sys.stderr)
            break

    if not files:
        files = sorted(glob.glob(os.path.join(directory, "*.tsv")))
        if files:
            print(f"WARNING: No {label} files with specific patterns, "
                  f"using all {len(files)} TSV files in {directory}", file=sys.stderr)
        else:
            print(f"WARNING: No {label} TSV files found in {directory}", file=sys.stderr)
            try:
                contents = os.listdir(directory)
                print(f"  Directory contents: {contents[:20]}", file=sys.stderr)
            except Exception as e:
                print(f"  Cannot list directory: {e}", file=sys.stderr)
            return {}

    index = {}
    for f in files:
        accession = extract_accession_from_filename(f)
        if not accession:
            accession = extract_accession_from_header(f)
        if accession:
            index[accession] = f
        else:
            print(f"WARNING: Could not determine accession for {f}", file=sys.stderr)

    print(f"Indexed {len(index)} {label} files: {sorted(list(index.keys()))[:5]}...", file=sys.stderr)
    return index


def group_biotype(biotype: str) -> str:
    """Collapse detailed biotype strings into broad display categories."""
    b = biotype.lower()
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


def load_tsv(filepath):
    """Load a single TSV file with robust engine."""
    try:
        return pd.read_csv(filepath, sep='\t', engine='python')
    except Exception as e:
        print(f"WARNING: Could not read {filepath}: {e}", file=sys.stderr)
        return None


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    threshold = args.rbh_coverage_threshold

    # Build file indexes (no data loaded yet)
    presence_idx = index_files(args.gene_presence_dir,
                               ["*_gene_presence.tsv", "*gene_presence*"],
                               "gene presence")
    rbh_idx = index_files(args.rbh_dir,
                          ["*.gene_pairs_rbh.tsv", "*rbh*", "*gene_pairs*"],
                          "RBH")
    tc_idx = index_files(args.transcript_concordance_dir,
                         ["*_transcript_concordance.tsv", "*transcript*"],
                         "transcript concordance")
    ci_idx = index_files(args.coding_integrity_dir,
                         ["*_coding_integrity.tsv", "*coding_integrity*"],
                         "coding integrity")

    # Assemblies that have at least gene presence + RBH
    all_accessions = sorted(set(presence_idx.keys()) & set(rbh_idx.keys()))
    print(f"Presence accessions: {sorted(list(presence_idx.keys()))[:5]}", file=sys.stderr)
    print(f"RBH accessions: {sorted(list(rbh_idx.keys()))[:5]}", file=sys.stderr)
    print(f"Processing {len(all_accessions)} assemblies with both presence and RBH data", file=sys.stderr)

    # Stream results: write per-assembly rows as we go
    out_path = os.path.join(args.output_dir, 'sankey_per_assembly_flows.tsv')
    header_written = False
    fieldnames = None
    n_written = 0

    with open(out_path, 'w', newline='') as out_fh:
        writer = None

        for i, accession in enumerate(all_accessions):
            if (i + 1) % 50 == 0:
                print(f"  Processing assembly {i + 1}/{len(all_accessions)}: {accession}", file=sys.stderr)

            # Load only this assembly's files
            pres = load_tsv(presence_idx[accession])
            rbh = load_tsv(rbh_idx[accession])

            if pres is None or pres.empty or rbh is None or rbh.empty:
                continue

            # Level 1: Gene presence (all genes, including ENSG-fallback names)
            n_both = ((pres['present_in_ensembl'].astype(str) == 'True') &
                      (pres['present_in_cat'].astype(str) == 'True')).sum()
            n_ens_only = ((pres['present_in_ensembl'].astype(str) == 'True') &
                          (pres['present_in_cat'].astype(str) == 'False')).sum()
            n_cat_only = ((pres['present_in_ensembl'].astype(str) == 'False') &
                          (pres['present_in_cat'].astype(str) == 'True')).sum()
            n_total = n_both + n_ens_only + n_cat_only

            # Level 1 (named only): exclude ENSG-fallback gene names
            if 'gene_name' in pres.columns:
                named = pres[~pres['gene_name'].astype(str).str.startswith('ENSG', na=False)]
            else:
                named = pres
            n_named_both = ((named['present_in_ensembl'].astype(str) == 'True') &
                            (named['present_in_cat'].astype(str) == 'True')).sum()
            n_named_ens_only = ((named['present_in_ensembl'].astype(str) == 'True') &
                                (named['present_in_cat'].astype(str) == 'False')).sum()
            n_named_cat_only = ((named['present_in_ensembl'].astype(str) == 'False') &
                                (named['present_in_cat'].astype(str) == 'True')).sum()
            n_named_total = n_named_both + n_named_ens_only + n_named_cat_only

            sample = pres['sample_name'].iloc[0] if 'sample_name' in pres.columns else ''
            del pres

            # Level 2: RBH status — build per-pair key index for downstream filtering
            n_rbh_total = len(rbh)
            rbh_keys = rbh['ensembl_id'].astype(str) + '|' + rbh['cat_id'].astype(str)
            bio_col = 'ensembl_biotype' if 'ensembl_biotype' in rbh.columns else None
            raw_bios = rbh[bio_col].astype(str) if bio_col else pd.Series([''] * len(rbh), index=rbh.index)
            biotype_by_key = dict(zip(rbh_keys, raw_bios.map(group_biotype)))

            if 'frac_ensembl_covered' in rbh.columns and 'frac_cat_covered' in rbh.columns:
                pass_mask = ((rbh['frac_ensembl_covered'].astype(float) >= threshold) &
                             (rbh['frac_cat_covered'].astype(float) >= threshold))
                pass_pair_keys = set(rbh_keys[pass_mask])
            else:
                pass_pair_keys = set(rbh_keys)

            n_rbh_pass = len(pass_pair_keys)
            n_rbh_fail = n_rbh_total - n_rbh_pass
            del rbh

            # Level 3: Transcript concordance — restricted to RBH pass genes only
            n_tx_full = n_tx_partial = n_tx_none = 0
            full_concordant_keys = set()  # (ens|cat) pairs with concordance_rate == 1.0

            if accession in tc_idx:
                tc = load_tsv(tc_idx[accession])
                if tc is not None and not tc.empty and 'transcript_concordance_rate' in tc.columns:
                    tc_keys = tc['ensembl_gene_id'].astype(str) + '|' + tc['cat_gene_id'].astype(str)
                    tc_pass_mask = tc_keys.isin(pass_pair_keys)
                    rates = tc['transcript_concordance_rate'].astype(float)

                    rates_pass = rates[tc_pass_mask]
                    n_tx_full    = int((rates_pass >= 1.0).sum())
                    n_tx_partial = int(((rates_pass > 0) & (rates_pass < 1.0)).sum())
                    n_tx_none    = int((rates_pass == 0).sum())

                    full_concordant_keys = set(tc_keys[tc_pass_mask & (rates >= 1.0)])
                del tc

            n_tx_total = n_tx_full + n_tx_partial + n_tx_none

            # Level 3b: biotype breakdown of fully concordant genes
            n_l3b_coding = n_l3b_lncrna = n_l3b_pseudogene = n_l3b_other_ncrna = n_l3b_other = 0
            for key in full_concordant_keys:
                bio = biotype_by_key.get(key, 'other')
                if bio == 'protein_coding':
                    n_l3b_coding += 1
                elif bio == 'lncRNA':
                    n_l3b_lncrna += 1
                elif bio == 'pseudogene':
                    n_l3b_pseudogene += 1
                elif bio == 'other_ncRNA':
                    n_l3b_other_ncrna += 1
                else:
                    n_l3b_other += 1
            n_l3b_total = n_l3b_coding + n_l3b_lncrna + n_l3b_pseudogene + n_l3b_other_ncrna + n_l3b_other

            # Level 4: CDS integrity — restricted to protein-coding full-concordant genes only
            n_cds_intact = n_cds_partial = n_cds_disrupted = n_cds_total = 0
            if accession in ci_idx:
                ci = load_tsv(ci_idx[accession])
                if ci is not None and not ci.empty:
                    ci_keys = ci['ensembl_gene_id'].astype(str) + '|' + ci['cat_gene_id'].astype(str)
                    ci_filtered = ci[ci_keys.isin(full_concordant_keys)]
                    n_cds_total = len(ci_filtered)
                    if 'classification' in ci_filtered.columns and n_cds_total > 0:
                        n_cds_intact = int((ci_filtered['classification'] == 'Full_Match').sum())
                        disrupted_classes = {'Internal_Frameshift', 'Complex_Mismatch', 'No_CDS'}
                        n_cds_disrupted = int(ci_filtered['classification'].isin(disrupted_classes).sum())
                        n_cds_partial = n_cds_total - n_cds_intact - n_cds_disrupted
                del ci

            gc.collect()

            row = {
                'assembly_accession': accession,
                'sample_name': sample,
                # Level 1 (all genes, including ENSG-fallback)
                'l1_total': int(n_total),
                'l1_both': int(n_both),
                'l1_ensembl_only': int(n_ens_only),
                'l1_cat_only': int(n_cat_only),
                'l1_pct_both': round(n_both / n_total * 100, 2) if n_total > 0 else 0,
                # Level 1 named only (ENSG-fallback names excluded)
                'l1_named_total': int(n_named_total),
                'l1_named_both': int(n_named_both),
                'l1_named_ensembl_only': int(n_named_ens_only),
                'l1_named_cat_only': int(n_named_cat_only),
                'l1_named_pct_both': round(n_named_both / n_named_total * 100, 2) if n_named_total > 0 else 0,
                # Level 2: RBH pass/fail (denominator = all RBH pairs)
                'l2_rbh_total': int(n_rbh_total),
                'l2_rbh_pass': int(n_rbh_pass),
                'l2_rbh_fail': int(n_rbh_fail),
                'l2_pct_pass': round(n_rbh_pass / n_rbh_total * 100, 2) if n_rbh_total > 0 else 0,
                # Level 3: transcript concordance (denominator = RBH pass)
                'l3_tx_total': int(n_tx_total),
                'l3_tx_full': int(n_tx_full),
                'l3_tx_partial': int(n_tx_partial),
                'l3_tx_none': int(n_tx_none),
                'l3_pct_full': round(n_tx_full / n_tx_total * 100, 2) if n_tx_total > 0 else 0,
                # Level 3b: biotype breakdown (denominator = L3 full concordant)
                'l3b_total': int(n_l3b_total),
                'l3b_protein_coding': int(n_l3b_coding),
                'l3b_lncrna': int(n_l3b_lncrna),
                'l3b_pseudogene': int(n_l3b_pseudogene),
                'l3b_other_ncrna': int(n_l3b_other_ncrna),
                'l3b_other': int(n_l3b_other),
                'l3b_pct_protein_coding': round(n_l3b_coding / n_l3b_total * 100, 2) if n_l3b_total > 0 else 0,
                # Level 4: CDS integrity (denominator = coding full-concordant genes)
                'l4_cds_total': int(n_cds_total),
                'l4_cds_intact': int(n_cds_intact),
                'l4_cds_partial': int(n_cds_partial),
                'l4_cds_disrupted': int(n_cds_disrupted),
                'l4_pct_intact': round(n_cds_intact / n_cds_total * 100, 2) if n_cds_total > 0 else 0,
            }

            # Write row
            if not header_written:
                fieldnames = list(row.keys())
                writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                header_written = True
            writer.writerow(row)
            n_written += 1

    # If nothing was written, create file with just a header
    if not header_written:
        print("WARNING: No assembly data found - writing empty output files", file=sys.stderr)
        with open(out_path, 'w') as f:
            f.write('\t'.join([
                'assembly_accession', 'sample_name',
                'l1_total', 'l1_both', 'l1_ensembl_only', 'l1_cat_only', 'l1_pct_both',
                'l1_named_total', 'l1_named_both', 'l1_named_ensembl_only', 'l1_named_cat_only', 'l1_named_pct_both',
                'l2_rbh_total', 'l2_rbh_pass', 'l2_rbh_fail', 'l2_pct_pass',
                'l3_tx_total', 'l3_tx_full', 'l3_tx_partial', 'l3_tx_none', 'l3_pct_full',
                'l3b_total', 'l3b_protein_coding', 'l3b_lncrna', 'l3b_pseudogene',
                'l3b_other_ncrna', 'l3b_other', 'l3b_pct_protein_coding',
                'l4_cds_total', 'l4_cds_intact', 'l4_cds_partial', 'l4_cds_disrupted', 'l4_pct_intact',
            ]) + '\n')

    # Now read back the per-assembly file for aggregate summary
    per_asm_df = pd.read_csv(out_path, sep='\t', engine='python')

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
    else:
        pd.DataFrame(columns=['n_assemblies', 'rbh_coverage_threshold']).to_csv(
            os.path.join(args.output_dir, 'sankey_flow_counts.tsv'),
            sep='\t', index=False
        )

    print(f"Wrote Sankey flow data for {n_written} assemblies", file=sys.stderr)


if __name__ == "__main__":
    main()
