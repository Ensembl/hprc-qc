#!/usr/bin/env python3
"""
Aggregate adjacent-rung flows linking RBH → transcript concordance → biotype → (CDS) → GRCh38 divergence.

Inputs (directories of per-assembly TSVs):
  --rbh-dir, --transcript-concordance-dir, --coding-integrity-dir, --divergence-dir

Outputs (to --output-dir):
  - links_per_assembly.tsv   one row per RBH pair with rung assignments
  - flows_per_assembly.tsv   per-assembly counts for each adjacent rung pair (A,B,C figures)
  - flows_medians.tsv        cohort medians for each adjacent rung pair

Notes:
  - RBH pass threshold defaults to 0.95 on both coverage fractions.
  - Biotype bucketing matches aggregate_sankey_flows.py.
  - CDS rung is protein_coding only.
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict, Counter

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate flows for integrated GRCh38 Sankey")
    p.add_argument('--rbh-dir', required=True)
    p.add_argument('--transcript-concordance-dir', required=True)
    p.add_argument('--coding-integrity-dir', required=True)
    p.add_argument('--divergence-dir', required=True)
    p.add_argument('--multi-mapping-dir', default=None,
                   help='Optional directory of *_multi_mapping.tsv files. When provided, '
                        'ens_multi_mapping and cat_multi_mapping columns are populated.')
    p.add_argument('--rbh-coverage-threshold', type=float, default=0.95)
    p.add_argument('--output-dir', required=True)
    return p.parse_args()


def group_biotype(biotype: str) -> str:
    b = str(biotype or '').lower()
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


_ACC_RE = re.compile(r'^(GCA_\d+\.\d+)', re.IGNORECASE)

def load_first(tsv_dir: str, pattern: str) -> dict:
    """Index per-assembly files by accession.

    Handles both underscore-separated names (e.g. GCA_018466835.2_transcript_concordance.tsv)
    and dot-separated names produced by hprc_ensembl_cat_overlap.py
    (e.g. GCA_018466835.2.gene_pairs_rbh.tsv).
    """
    idx = {}
    for f in glob.glob(os.path.join(tsv_dir, pattern)):
        m = _ACC_RE.match(os.path.basename(f))
        if m:
            idx[m.group(1)] = f
        else:
            print(f"WARNING: could not parse accession from {os.path.basename(f)}", file=sys.stderr)
    return idx


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Index files by assembly accession (prefix of filename)
    rbh_idx = load_first(args.rbh_dir, "*.gene_pairs_rbh.tsv") or load_first(args.rbh_dir, "*rbh*.tsv")
    tc_idx = load_first(args.transcript_concordance_dir, "*_transcript_concordance.tsv")
    ci_idx = load_first(args.coding_integrity_dir, "*_coding_integrity.tsv")
    div_idx = load_first(args.divergence_dir, "*_grch38_divergence.tsv")
    mm_idx = load_first(args.multi_mapping_dir, "*_multi_mapping.tsv") if args.multi_mapping_dir else {}

    accessions = sorted(set(rbh_idx) & set(tc_idx) & set(div_idx))
    if not accessions:
        print("ERROR: No overlapping assemblies across required inputs", file=sys.stderr)
        sys.exit(1)

    # Outputs
    links_rows = []
    flows_rows = []  # per assembly, per adjacency key counts

    thr = float(args.rbh_coverage_threshold)

    for acc in accessions:
        # RBH — used to classify pass/fail status and supply cat_biotype for each pair
        rbh = pd.read_csv(rbh_idx[acc], sep='\t', dtype=str)
        rbh['pair_key'] = rbh['ensembl_id'].astype(str) + '|' + rbh['cat_id'].astype(str)
        if {'frac_ensembl_covered','frac_cat_covered'}.issubset(rbh.columns):
            pass_mask = (rbh['frac_ensembl_covered'].astype(float) >= thr) & (rbh['frac_cat_covered'].astype(float) >= thr)
        else:
            pass_mask = pd.Series(True, index=rbh.index)
        rbh_pass_keys = set(rbh.loc[pass_mask, 'pair_key'])
        rbh_all_keys  = set(rbh['pair_key'])
        sample = rbh['sample_name'].iloc[0] if 'sample_name' in rbh.columns else ''
        # cat_biotype from RBH file (indexed by pair_key for O(1) lookup)
        rbh_cat_bio = (rbh.drop_duplicates('pair_key').set_index('pair_key')['cat_biotype']
                       if 'cat_biotype' in rbh.columns else pd.Series(dtype=str))

        # Multi-mapping classifications (optional — null if file not provided or gene not multi-mapped)
        mm_ens_class: dict = {}  # ensembl_gene_id → classification
        mm_cat_class: dict = {}  # cat_gene_id       → classification
        if acc in mm_idx:
            try:
                mm = pd.read_csv(mm_idx[acc], sep='\t', dtype=str)
                if 'source' in mm.columns and 'gene_id' in mm.columns and 'classification' in mm.columns:
                    mm_ens_class = dict(zip(
                        mm.loc[mm['source'] == 'ensembl', 'gene_id'],
                        mm.loc[mm['source'] == 'ensembl', 'classification']
                    ))
                    mm_cat_class = dict(zip(
                        mm.loc[mm['source'] == 'cat', 'gene_id'],
                        mm.loc[mm['source'] == 'cat', 'classification']
                    ))
            except Exception as e:
                print(f"WARNING [{acc}]: could not load multi-mapping file: {e}", file=sys.stderr)

        # Transcript concordance — all RBH gene pairs (pass + fail).
        # tc carries ensembl_biotype directly (added in calculate_transcript_concordance.py).
        tc = pd.read_csv(tc_idx[acc], sep='\t', dtype=str)
        tc['pair_key'] = tc['ensembl_gene_id'].astype(str) + '|' + tc['cat_gene_id'].astype(str)
        tc_indexed = tc.set_index('pair_key')
        tc_rate = tc_indexed['transcript_concordance_rate'].astype(float)
        tc_biotype = tc_indexed['ensembl_biotype'].map(group_biotype) if 'ensembl_biotype' in tc_indexed.columns else pd.Series(dtype=str)

        # Validate pair_key equivalence: TC is driven by RBH, so every TC pair should be in rbh_all_keys.
        # A mismatch here indicates a column-naming or file-mismatch bug rather than biology.
        orphan_tc_keys = set(tc_indexed.index) - rbh_all_keys
        if orphan_tc_keys:
            print(f"WARNING [{acc}]: {len(orphan_tc_keys)} TC pair_keys not found in RBH file "
                  f"(will be assigned rbh_status='no_rbh'). "
                  f"Verify that transcript_concordance was run with gene_pairs_rbh.tsv.",
                  file=sys.stderr)

        # Coding integrity (optional per assembly)
        ci = None
        if acc in ci_idx:
            try:
                ci = pd.read_csv(ci_idx[acc], sep='\t', dtype=str)
                ci['pair_key'] = ci['ensembl_gene_id'].astype(str) + '|' + ci['cat_gene_id'].astype(str)
                ci = ci.set_index('pair_key')
            except Exception:
                ci = None

        # Divergence vs GRCh38 (RBH pairs only — non-RBH pairs get null divergence)
        dv = pd.read_csv(div_idx[acc], sep='\t', dtype=str)
        dv['pair_key'] = dv['ensembl_gene_id'].astype(str) + '|' + dv['cat_gene_id'].astype(str)
        dv = dv.set_index('pair_key')

        # Build link rows for ALL pairs present in the transcript concordance file
        for k in tc_indexed.index:
            # rbh_status: 'pass' | 'fail' | 'no_rbh'
            if k in rbh_pass_keys:
                rbh_status = 'pass'
            elif k in rbh_all_keys:
                rbh_status = 'fail'
            else:
                rbh_status = 'no_rbh'

            v = float(tc_rate[k])
            txc = 'full' if v >= 1.0 else ('partial' if v > 0 else 'none')

            bio = tc_biotype.get(k, 'other') if k in tc_biotype.index else 'other'

            cds_class = None
            if bio == 'protein_coding' and ci is not None and k in ci.index:
                cls = ci.loc[k]['classification'] if 'classification' in ci.columns else None
                if cls == 'Full_Match':
                    cds_class = 'intact'
                elif cls in {'Internal_Frameshift', 'Complex_Mismatch', 'No_CDS'}:
                    cds_class = 'disrupted'
                else:
                    cds_class = 'partial'

            # Divergence — only available for RBH pairs
            if k in dv.index:
                div_cat    = dv.loc[k]['divergence_category'] if 'divergence_category' in dv.columns else None
                gene_name  = dv.loc[k]['gene_name'] if 'gene_name' in dv.columns else ''
                ens_cds_ct = dv.loc[k]['ens_cds_change_type'] if 'ens_cds_change_type' in dv.columns else ''
                cat_cds_ct = dv.loc[k]['cat_cds_change_type'] if 'cat_cds_change_type' in dv.columns else ''
            else:
                div_cat = gene_name = ens_cds_ct = cat_cds_ct = None

            # cat_biotype from RBH file (same pair_key construction, values are identical)
            cat_bio_raw = rbh_cat_bio.get(k) if k in rbh_cat_bio.index else None

            # Multi-mapping classification — null means 1-to-1 (not flagged as multi-mapped)
            ens_id, cat_id = k.split('|', 1)
            ens_mm = mm_ens_class.get(ens_id) or None
            cat_mm = mm_cat_class.get(cat_id) or None

            links_rows.append({
                'assembly_accession': acc,
                'sample_name': sample,
                'pair_key': k,
                'ensembl_gene_id': ens_id,
                'cat_gene_id': cat_id,
                'gene_name': gene_name,
                'rbh_status': rbh_status,
                'tx_concordance': txc,
                'biotype': bio,
                'cat_biotype': cat_bio_raw,
                'cds_class': cds_class,
                'ens_cds_change_type': ens_cds_ct,
                'cat_cds_change_type': cat_cds_ct,
                'divergence_category': div_cat,
                'ens_multi_mapping': ens_mm,
                'cat_multi_mapping': cat_mm,
            })

        # Per-assembly flows for A/B/C Sankey diagrams.
        # Flows are restricted to RBH pairs only (pass + fail) to keep the
        # funnel semantics: gene presence → RBH → concordance → biotype → GRCh38.
        df_links = pd.DataFrame([r for r in links_rows
                                  if r['assembly_accession'] == acc
                                  and r['rbh_status'] in {'pass', 'fail'}])
        if df_links.empty:
            continue

        def add_counts(level_from, level_to, name):
            sub = df_links.dropna(subset=[level_from, level_to])
            counts = sub.groupby([level_from, level_to]).size().reset_index(name='n')
            for _, row in counts.iterrows():
                flows_rows.append({
                    'assembly_accession': acc,
                    'sample_name': sample,
                    'adjacency': name,
                    'from': row[level_from],
                    'to': row[level_to],
                    'n': int(row['n']),
                })

        # Figure A adjacencies
        add_counts('rbh_status', 'tx_concordance', 'A_RBH_to_Concordance')
        add_counts('tx_concordance', 'biotype', 'A_Concordance_to_Biotype')
        add_counts('biotype', 'divergence_category', 'A_Biotype_to_GRCh38')

        # Figure B adjacencies (protein_coding path only for the last two)
        pc = df_links[df_links['biotype'] == 'protein_coding'].copy()
        if not pc.empty:
            def add_counts_pc(level_from, level_to, name):
                sub = pc.dropna(subset=[level_from, level_to])
                counts = sub.groupby([level_from, level_to]).size().reset_index(name='n')
                for _, row in counts.iterrows():
                    flows_rows.append({
                        'assembly_accession': acc,
                        'sample_name': sample,
                        'adjacency': name,
                        'from': row[level_from],
                        'to': row[level_to],
                        'n': int(row['n']),
                    })

            add_counts_pc('biotype', 'cds_class', 'B_Biotype_to_CDS')
            add_counts_pc('cds_class', 'divergence_category', 'B_CDS_to_GRCh38')

        # Figure C adjacencies (overall + pc-only)
        add_counts('tx_concordance', 'divergence_category', 'C_Concordance_to_GRCh38')
        if not pc.empty:
            dfc = pc
            counts = dfc.dropna(subset=['tx_concordance', 'divergence_category']).\
                groupby(['tx_concordance', 'divergence_category']).size().reset_index(name='n')
            for _, row in counts.iterrows():
                flows_rows.append({
                    'assembly_accession': acc,
                    'sample_name': sample,
                    'adjacency': 'C_Concordance_to_GRCh38_protein_coding',
                    'from': row['tx_concordance'],
                    'to': row['divergence_category'],
                    'n': int(row['n']),
                })

    # Write links_per_assembly
    links_df = pd.DataFrame(links_rows)
    links_df.to_csv(os.path.join(args.output_dir, 'links_per_assembly.tsv'), sep='\t', index=False)

    # Write flows_per_assembly
    flows_df = pd.DataFrame(flows_rows)
    flows_df.to_csv(os.path.join(args.output_dir, 'flows_per_assembly.tsv'), sep='\t', index=False)

    # Compute medians per adjacency/from/to across assemblies
    if not flows_df.empty:
        med_rows = []
        for (adj, frm, to), sub in flows_df.groupby(['adjacency', 'from', 'to']):
            med = sub['n'].median()
            mean = sub['n'].mean()
            q25 = sub['n'].quantile(0.25)
            q75 = sub['n'].quantile(0.75)
            med_rows.append({
                'adjacency': adj,
                'from': frm,
                'to': to,
                'median': round(float(med), 2),
                'mean': round(float(mean), 2),
                'q25': round(float(q25), 2),
                'q75': round(float(q75), 2),
            })
        pd.DataFrame(med_rows).to_csv(os.path.join(args.output_dir, 'flows_medians.tsv'), sep='\t', index=False)
    else:
        pd.DataFrame(columns=['adjacency', 'from', 'to', 'median', 'mean', 'q25', 'q75']).to_csv(
            os.path.join(args.output_dir, 'flows_medians.tsv'), sep='\t', index=False
        )

    print(f"Wrote integrated GRCh38 Sankey tables for {len(accessions)} assemblies", file=sys.stderr)


if __name__ == '__main__':
    main()

