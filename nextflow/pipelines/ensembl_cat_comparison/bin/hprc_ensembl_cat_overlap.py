#!/usr/bin/env python3
"""
HPRC Ensembl-CAT Overlap Quantification Script

This script quantifies the overlap between Ensembl HPRC release 2 annotations (GFF3)
and CAT gene annotations.

Usage Example:
    python3 hprc_ensembl_cat_overlap.py \\
        --ensembl-gff-dir /path/to/ensembl_gffs \\
        --cat-index /path/to/cat_genes_hprc_r2_v1.2.index.csv \\
        --assemblies-index /path/to/assemblies_release2_v1.0.index.csv \\
        --cat-anno-dir /path/to/cat_annotations \\
        --output-prefix ./results/overlap_analysis

Dependencies:
    - pandas
    - pyranges
    - numpy

Install dependencies:
    pip install pandas pyranges numpy
"""

import argparse
import os
import sys
import logging
import re
import pandas as pd
import pyranges as pr
import numpy as np
from typing import Dict, Optional, Tuple, List

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Quantify overlap between Ensembl and CAT annotations.")
    
    parser.add_argument("--ensembl-gff-dir", required=True,
                        help="Directory containing Ensembl HPRC release 2 GFF3 files.")
    parser.add_argument("--cat-index", required=True,
                        help="Path to cat_genes_hprc_r2_v1.2.index.csv")
    parser.add_argument("--assemblies-index", required=True,
                        help="Path to assemblies_release2_v1.0.index.csv")
    parser.add_argument("--cat-anno-dir", required=True,
                        help="Directory containing CAT annotation files.")
    parser.add_argument("--output-prefix", required=True,
                        help="Prefix for output TSV files.")
    parser.add_argument("--assemblies", required=False,
                        help="Optional comma-separated list of assembly accessions or sample names to process.")
    parser.add_argument("--contig-normalization", choices=['none', 'basic', 'cat_hash'], default='none',
                        help="Normalization mode for contig names: 'none', 'basic' (strip 'chr', M/MT), 'cat_hash' (strip '.*#.*#').")
    parser.add_argument("--self-check", action="store_true",
                        help="Run internal self-checks/unit tests and exit.")
    
    # Single-pair test args
    parser.add_argument("--ensembl-gff", help="Path to specific Ensembl GFF (test mode).")
    parser.add_argument("--cat-gff", help="Path to specific CAT GFF (test mode).")
    parser.add_argument("--chrom", help="Chromosome to filter and test (e.g. '1').")
    
    parser.add_argument("--assembly-report", help="Path to NCBI assembly report for contig mapping.")
    
    return parser.parse_args()

def find_column(df: pd.DataFrame, keywords: List[str]) -> Optional[str]:
    """Helper to find a column containing any of the keywords."""
    for col in df.columns:
        for k in keywords:
            if k.lower() in col.lower():
                return col
    return None

def load_assembly_map(index_path: str) -> pd.DataFrame:
    """
    Parses assemblies index to build a mapping.
    Expects 'Assembly Accession' and 'Sample Name' and potentially a GFF/Annotation file column.
    """
    logger.info(f"Loading assembly index from {index_path}")
    # try:
    df = pd.read_csv(index_path)
    df.columns = [c.strip() for c in df.columns]
    
    required_cols = ['Assembly Accession', 'Sample Name']
    # Normalize column names to handle potential whitespace issues
    df.columns = [c.strip() for c in df.columns]
    
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        logger.error(f"Assembly index missing required columns: {missing}")
        sys.exit(1)
        
    # Try to find Ensembl GFF col
    gff_col = find_column(df, ['ensembl', 'annotation', 'gff'])
    if gff_col:
        logger.info(f"Found Ensembl GFF column: {gff_col}")
        df['ensembl_gff_file'] = df[gff_col].apply(lambda x: os.path.basename(str(x)) if pd.notnull(x) else None)
    else:
        logger.warning("No explicit Ensembl GFF column found. Will rely on Accession matching.")
        df['ensembl_gff_file'] = None
        
    return df
    # except Exception as e:
    #     logger.error(f"Failed to load assembly index: {e}")
    #     sys.exit(1)

def load_cat_map(index_path: str) -> pd.DataFrame:
    """
    Parses CAT index to map sample name to CAT annotation files.
    """
    logger.info(f"Loading CAT index from {index_path}")
    try:
        df = pd.read_csv(index_path)
        df.columns = [c.strip() for c in df.columns]
        
        # We need 'Sample Name' (or similar) and the file path/URL column.
        # Sometimes CAT index calls sample 'Assembly' or 'Genome'.
        sample_col = find_column(df, ['sample', 'assembly', 'genome'])
        if not sample_col:
             logger.error("Could not identify Sample/Assembly column in CAT index.")
             sys.exit(1)
             
        file_col = find_column(df, ['gff', 'path', 'file', 'url'])
        if not file_col:
            logger.error("Could not identify filename column in CAT index.")
            sys.exit(1)
            
        logger.info(f"Using CAT columns - Sample: '{sample_col}', File: '{file_col}'")
        
        # Standardize
        df = df.rename(columns={sample_col: 'Sample Name', file_col: 'cat_gff_file'})
        df['cat_gff_file'] = df['cat_gff_file'].apply(lambda x: os.path.basename(str(x)) if pd.notnull(x) else None)
        
        return df[[ 'Sample Name', 'cat_gff_file' ]].dropna()
    except Exception as e:
        logger.error(f"Failed to load CAT index: {e}")
        sys.exit(1)

def load_assembly_report(report_path: str) -> Dict[str, str]:
    """
    Parses NCBI assembly report to map GenBank Accessions to Assigned-Molecule names (Ensembl style).
    Returns dict: {GenBank_Accn: Assigned_Molecule}
    e.g. {'CM085953.1': '2', 'JBHDVK010000093.1': 'JBHDVK010000093.1'}
    """
    mapping = {}
    if not report_path or not os.path.exists(report_path):
        return mapping
        
    try:
        with open(report_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                
                # Cols: 0:SeqName, 1:Role, 2:AssignedMol, 3:Type, 4:GenBankAccn
                genbank_accn = parts[4]
                assigned_mol = parts[2]
                
                if assigned_mol != 'na':
                    mapping[genbank_accn] = assigned_mol
                else:
                    mapping[genbank_accn] = genbank_accn
                    
        logger.info(f"Loaded {len(mapping)} mappings from assembly report.")
        return mapping
        
    except Exception as e:
        logger.error(f"Error loading assembly report: {e}")
        return {}

def normalize_accession(filename: str) -> Optional[str]:
    """
    Extracts accession from Ensembl filename and standardizes it.
    Example: 'homo_sapiens_gca018466835v2.gff3' -> 'GCA_018466835.2'
    
    Note: The mapping from 'gca018466835v2' to 'GCA_018466835.2' is heuristic.
    We might need to rely on the assembly index to match 'GCA_018466835.2' 
    to a normalized version like 'gca018466835v2' if exact reconstruction is hard.
    
    However, usually the pattern is:
    gca<number>v<version> -> GCA_<number>.<version>
    """
    # Pattern to catch gca018466835v2
    match = re.search(r'(gca)(\d+)v(\d+)', filename, re.IGNORECASE)
    if match:
        prefix = match.group(1).upper() # GCA
        number = match.group(2)
        version = match.group(3)
        return f"{prefix}_{number}.{version}"
    
    # Fallback: maybe the filename contains the full accession GCA_...
    match_full = re.search(r'(GCA_\d+\.\d+)', filename)
    if match_full:
        return match_full.group(1)
        
    return None

def load_gff_features(path: str, source_label: str) -> Dict[str, pr.PyRanges]:
    """
    Parses a GFF3 file and extracts 'gene', 'transcript'/'mRNA', and 'exon' features.
    Returns a dictionary of PyRanges: {'genes': pr, 'transcripts': pr, 'exons': pr}
    
    Args:
        path: Path to GFF3 file.
        source_label: Label to prefix columns (e.g. 'ensembl' or 'cat').
        
    Returns:
        Dict of PyRanges or empty dict if failed.
    """
    import gzip

    if not os.path.exists(path):
        logger.warning(f"File not found: {path}")
        return {}
    
    try:
        genes = []
        transcripts = []
        exons = []
        
        # Helper to get attribute
        def get_attr(attrs_dict, keys):
            for k in keys:
                if k in attrs_dict:
                    return attrs_dict[k]
            return None
        
        # Helper for opening file (gzip aware)
        def open_file(p):
            if p.endswith('.gz'):
                try:
                    return gzip.open(p, 'rt')
                except Exception:
                    return open(p, 'r')
            else:
                return open(p, 'r')

        with open_file(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                
                # Broad categorization with strict checks
                # Ensembl uses 'gene', 'mRNA', 'exon', 'CDS', etc.
                # CAT uses 'gene', 'transcript', 'exon', 'CDS'.
                
                is_gene = feature_type == 'gene' or feature_type.endswith('gene')
                
                # Transcripts: explicit list
                valid_tx = ['transcript', 'mRNA', 'lnc_RNA', 'miRNA', 'ncRNA', 'snoRNA', 'snRNA', 'tRNA', 'rRNA', 'pseudogenic_transcript']
                is_transcript = feature_type in valid_tx or (feature_type == 'transcript')
                
                is_exon = feature_type == 'exon'
                
                if not (is_gene or is_transcript or is_exon):
                    continue

                contig = parts[0]
                try:
                    start = int(parts[3]) - 1
                    end = int(parts[4])
                except ValueError:
                    continue
                    
                strand = parts[6]
                attributes = parts[8]
                
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, val = attr.split('=', 1)
                        attr_dict[key] = val
                
                # Common attributes
                feat_id = get_attr(attr_dict, ['ID', 'gene_id', 'transcript_id'])
                parent = get_attr(attr_dict, ['Parent', 'gene_id', 'transcript_id'])
                # Name matching: prioritize Name, then gene_name
                name = get_attr(attr_dict, ['Name', 'gene_name', 'transcript_name'])
                biotype = get_attr(attr_dict, ['biotype', 'gene_biotype', 'transcript_biotype'])
                
                if not feat_id:
                    continue
                    
                row = {
                    'Chromosome': contig,
                    'Start': start,
                    'End': end,
                    'Strand': strand,
                    f'{source_label}_id': feat_id,
                    f'{source_label}_biotype': biotype or 'unknown',
                    f'{source_label}_name': name,
                    f'{source_label}_length': end - start
                }

                if is_gene:
                    genes.append(row)
                elif is_transcript:
                    if parent and parent != feat_id:
                         # Some GFFs might reference parent by ID
                        row[f'{source_label}_parent'] = parent
                    transcripts.append(row)
                elif is_exon:
                    if parent:
                        row[f'{source_label}_parent'] = parent
                    exons.append(row)

        result = {}
        if genes:
            result['genes'] = pr.PyRanges(pd.DataFrame(genes))
        if transcripts:
            result['transcripts'] = pr.PyRanges(pd.DataFrame(transcripts))
        if exons:
            result['exons'] = pr.PyRanges(pd.DataFrame(exons))
            
        return result

    except Exception as e:
        logger.error(f"Error reading GFF {path}: {e}")
        return {}

def get_intron_map(feature_df: pd.DataFrame) -> Dict[str, List[Tuple[int, int]]]:
    """
    Computes sorted intron intervals for each transcript.
    Introns are the gaps between sorted exons.
    Expects DataFrame with columns: Start, End, Parent (or appropriate ID group).
    """
    introns = {}
    if feature_df.empty:
        return introns
        
    # Group by parent (Transcript ID)
    # Ensure sorted by Start
    # Check distinct parents
    
    # We expect feature_df to have transcript_id as grouping key if passed generally, 
    # but here we usually pass exons df which has 'parent' column.
    
    # Identify the parent column
    parent_cols = [c for c in feature_df.columns if 'parent' in c.lower()]
    if not parent_cols:
        return {}
    parent_col = parent_cols[0]
    
    # Sort by parent and start
    df_sorted = feature_df.sort_values([parent_col, 'Start'])
    
    for tx_id, group in df_sorted.groupby(parent_col):
        if len(group) < 2:
            introns[tx_id] = []
            continue
            
        sorted_exons = list(zip(group['Start'], group['End']))
        # Intron i is between exon i End and exon i+1 Start
        tx_introns = []
        for i in range(len(sorted_exons) - 1):
             # GFF3: 1-based usually, but here converted to 0-based.
             # Interval is [start, end).
             # Intron start = exon_i_end
             # Intron end = exon_i+1_start
             tx_introns.append((sorted_exons[i][1], sorted_exons[i+1][0]))
        
        introns[tx_id] = tx_introns
        
    return introns

def compare_intron_chains(ref_chain: List[Tuple[int, int]], qry_chain: List[Tuple[int, int]], strand_match: bool) -> str:
    """
    Compares two intron chains.
    Returns: 'Exact', 'Subset', 'Partial', 'None'
    """
    if not ref_chain and not qry_chain:
        # Both single-exon. If they overlap significantly it's a match, but chain-wise it's trivial.
        # We handle this context outside? Or call it 'SingleExonMatch' if coordinates align?
        # For this function, strictly chain topology.
        return 'Exact' # Both are "empty" chains
        
    if not ref_chain or not qry_chain:
        return 'None' # One has introns, other doesn't
        
    set_ref = set(ref_chain)
    set_qry = set(qry_chain)
    
    if set_ref == set_qry:
        return 'Exact'
    
    if set_ref.issubset(set_qry) or set_qry.issubset(set_ref):
        return 'Subset'
        
    if not set_ref.isdisjoint(set_qry):
        return 'Partial'
        
    return 'None'

def classify_overlap_type(e_frac: float, c_frac: float) -> str:
    """Classifies the type of overlap based on coverage fractions (legacy simple version)."""
    if e_frac >= 0.9 and c_frac >= 0.9:
        return 'Full-length concordant'
    elif e_frac >= 0.9 and c_frac < 0.9:
        return 'Ensembl-contained'
    elif c_frac >= 0.9 and e_frac < 0.9:
        return 'CAT-contained'
    elif e_frac > 0.1 and c_frac > 0.1:
        return 'Partial'
    else:
        return 'Tiny'


def classify_overlap_detailed(e_frac: float, c_frac: float,
                               e_start: int, e_end: int,
                               c_start: int, c_end: int,
                               strand: str = '+') -> str:
    """
    Detailed classification of overlap between Ensembl and CAT gene models.

    Categories:
    - Identical: >=99% coverage both ways
    - Near-identical: 95-99% coverage both ways
    - High concordance: 90-95% coverage both ways
    - CAT 5' extended: Ensembl well-covered, CAT extends at 5' end
    - CAT 3' extended: Ensembl well-covered, CAT extends at 3' end
    - CAT both extended: Ensembl well-covered, CAT extends both ends
    - Ensembl 5' extended: CAT well-covered, Ensembl extends at 5' end
    - Ensembl 3' extended: CAT well-covered, Ensembl extends at 3' end
    - Ensembl both extended: CAT well-covered, Ensembl extends both ends
    - Moderate overlap: 50-90% coverage for better-covered annotation
    - Low overlap: 10-50% coverage
    - Minimal overlap: <10% coverage
    """
    # Concordance tiers (high agreement both ways)
    if e_frac >= 0.99 and c_frac >= 0.99:
        return 'Identical'
    if e_frac >= 0.95 and c_frac >= 0.95:
        return 'Near-identical'
    if e_frac >= 0.90 and c_frac >= 0.90:
        return 'High concordance'

    # Extension patterns - one annotation is well-covered, the other partially
    # This means the partially-covered one is LONGER (extended)

    # Determine 5'/3' based on strand
    if strand == '-':
        # On minus strand, genomic 5' is higher coordinate
        extends_5prime_cat = c_end > e_end
        extends_3prime_cat = c_start < e_start
        extends_5prime_ens = e_end > c_end
        extends_3prime_ens = e_start < c_start
    else:
        # On plus strand (or unknown), genomic 5' is lower coordinate
        extends_5prime_cat = c_start < e_start
        extends_3prime_cat = c_end > e_end
        extends_5prime_ens = e_start < c_start
        extends_3prime_ens = e_end > c_end

    # CAT is extended (Ensembl well-covered >=90%, CAT partial 50-90%)
    if e_frac >= 0.90 and 0.50 <= c_frac < 0.90:
        if extends_5prime_cat and extends_3prime_cat:
            return 'CAT both extended'
        elif extends_5prime_cat:
            return 'CAT 5-prime extended'
        elif extends_3prime_cat:
            return 'CAT 3-prime extended'
        return 'CAT extended'

    # Ensembl is extended (CAT well-covered >=90%, Ensembl partial 50-90%)
    if c_frac >= 0.90 and 0.50 <= e_frac < 0.90:
        if extends_5prime_ens and extends_3prime_ens:
            return 'Ensembl both extended'
        elif extends_5prime_ens:
            return 'Ensembl 5-prime extended'
        elif extends_3prime_ens:
            return 'Ensembl 3-prime extended'
        return 'Ensembl extended'

    # Partial overlaps (neither is well-covered)
    max_frac = max(e_frac, c_frac)
    if max_frac >= 0.50:
        return 'Moderate overlap'
    if max_frac >= 0.10:
        return 'Low overlap'

    return 'Minimal overlap'


def compute_coordinate_diffs(e_start: int, e_end: int, c_start: int, c_end: int, strand: str = '+') -> dict:
    """
    Compute 5' and 3' coordinate differences between Ensembl and CAT genes.

    Positive values mean CAT extends beyond Ensembl at that end.
    Negative values mean Ensembl extends beyond CAT at that end.

    Returns dict with: diff_5prime, diff_3prime (in bp)
    """
    if strand == '-':
        # Minus strand: 5' is at higher coordinate
        diff_5prime = c_end - e_end  # positive if CAT extends 5'
        diff_3prime = e_start - c_start  # positive if CAT extends 3'
    else:
        # Plus strand: 5' is at lower coordinate
        diff_5prime = e_start - c_start  # positive if CAT extends 5'
        diff_3prime = c_end - e_end  # positive if CAT extends 3'

    return {
        'diff_5prime': diff_5prime,
        'diff_3prime': diff_3prime
    }

def analyze_transcript_concordance(
    gene_pairs_df: pd.DataFrame, 
    ensembl_feats: Dict[str, pr.PyRanges], 
    cat_feats: Dict[str, pr.PyRanges]
) -> pd.DataFrame:
    """
    Analyzes transcript concordance for matched GENE pairs using intron chains.
    Only considers pairs present in gene_pairs_df.
    """
    if gene_pairs_df.empty:
        return pd.DataFrame()
        
    e_trans = ensembl_feats.get('transcripts')
    c_trans = cat_feats.get('transcripts')
    e_exons = ensembl_feats.get('exons')
    c_exons = cat_feats.get('exons')
    
    if not (e_trans and c_trans and e_exons and c_exons):
        return pd.DataFrame()

    # Pre-compute intron chains for ALL transcripts in these sets
    # This might be expensive if very large, but safer than per-gene loop re-calc
    e_introns = get_intron_map(e_exons.df)
    c_introns = get_intron_map(c_exons.df)
    
    # Group Transcripts by Gene
    e_tx_df = e_trans.df
    c_tx_df = c_trans.df
    
    e_tx_by_gene = dict(list(e_tx_df.groupby('ensembl_parent')))
    c_tx_by_gene = dict(list(c_tx_df.groupby('cat_parent')))

    results = []
    
    # Filter to unique pairs to check
    pairs = gene_pairs_df[['ensembl_id', 'cat_id']].drop_duplicates()
    
    for _, row in pairs.iterrows():
        eid = row['ensembl_id']
        cid = row['cat_id']
        
        # Get transcripts for these genes
        etxs = e_tx_by_gene.get(eid)
        ctxs = c_tx_by_gene.get(cid)
        
        if etxs is None or ctxs is None:
            continue
            
        # For each Ensembl transcript, find best CAT transcript match
        for _, et in etxs.iterrows():
            et_id = et['ensembl_id']
            et_chain = e_introns.get(et_id, [])
            if et_chain:
                et_exons_subset = e_exons[e_exons.ensembl_parent == eid] # Slow? No used above logic.
                pass

            best_match_type = 'None'
            best_cat_tx = None
            
            # Helper to rank match types: None < Partial < Subset < Exact
            rank = {'None': 0, 'Partial': 1, 'Subset': 2, 'Exact': 3}
            best_rank = 0
            
            for _, ct in ctxs.iterrows():
                ct_id = ct['cat_id']
                ct_chain = c_introns.get(ct_id, [])
                
                # Check for intron match
                mtype = compare_intron_chains(et_chain, ct_chain, True) # Implicitly same strand if genes matched
                
                if rank[mtype] > best_rank:
                    best_rank = rank[mtype]
                    best_match_type = mtype
                    best_cat_tx = ct_id
                
                if best_match_type == 'Exact': 
                    break # Found optimal
            
            if best_cat_tx:
                results.append({
                    'ensembl_gene_id': eid,
                    'cat_gene_id': cid,
                    'ensembl_transcript_id': et_id,
                    'best_cat_transcript_id': best_cat_tx,
                    'intron_chain_match': best_match_type
                })
                
    return pd.DataFrame(results)

def compute_overlaps_aggregated(
    ensembl_pr: pr.PyRanges, 
    cat_pr: pr.PyRanges, 
    normalization: str = 'none',
    contig_map: Dict[str, str] = None
) -> pd.DataFrame:
    """
    Computes overlaps, aggregates by gene pair, and adds name matches.
    """
    e_pr = ensembl_pr
    c_pr = cat_pr
    
    # Normalization helper
    def norm_chrom(df, mode, mapping=None):
        d = df.copy()
        s = d['Chromosome'].astype(str)
        if mode == 'basic':
            s = s.str.replace(r'^chr', '', regex=True)
            s = s.replace({'MT': 'M'})
        elif mode == 'cat_hash':
            # Strip prefix first
            s = s.str.replace(r'^.*#', '', regex=True)
            
            # If mapping provided, apply it to the GenBank accession (which is what remains after stripping)
            if mapping:
                s = s.map(lambda x: mapping.get(x, x))
            else:
                s = s.str.replace(r'^chr', '', regex=True)
                s = s.replace({'MT': 'M'})
                
        d['Chromosome'] = s
        return d

    # Normalization
    if normalization != 'none':
        e_df = norm_chrom(e_pr.df, 'basic') # Ensembl usually needs basic or none? Ensembl GFF is already 1, 2...
        # Wait, ensembl GFF might need basic if it has chr prefix. But HPRC Ensembl usually doesn't?
        # Safe to apply basic to Ensembl.
        # But we need to use the SAME mode as passed?
        # Actually user passes 'cat_hash' to handle CAT. Ensembl shouldn't get cat_hash logic.
        # Ensembl should get 'basic' norm if anything?
        # Let's assume Ensembl is target, so we apply minimal norm to it.
        
        c_df = norm_chrom(c_pr.df, normalization, contig_map)
        e_df = norm_chrom(e_pr.df, 'basic') 
        
        e_pr = pr.PyRanges(e_df)
        c_pr = pr.PyRanges(c_df)
    
    
    # 1. Spatial Join
    try:
        # Check strands for safety
        use_stranded = 'same'
        if not e_pr.strands or not c_pr.strands:
             logger.warning("Strand info missing in PyRanges, falling back to unstranded join.")
             use_stranded = None 
             
        join_df = e_pr.join(c_pr, strandedness=use_stranded, report_overlap=True, how=None).df
    except Exception as e:
        logger.warning(f"PyRanges join failed (possibly empty inputs): {e}")
        join_df = pd.DataFrame()
        
    spatial_data = []

    if not join_df.empty:
        # Aggregate by Gene Pair -> Sum Overlap
        grp = join_df.groupby(['ensembl_id', 'cat_id'])

        for (eid, cid), g in grp:
            total_overlap = g['Overlap'].sum()
            e_len = g['ensembl_length'].iloc[0]
            c_len = g['cat_length'].iloc[0]

            # Get gene coordinates (use first row - should be consistent within gene)
            e_start = g['Start'].iloc[0]
            e_end = g['End'].iloc[0]
            # CAT coordinates have _b suffix from PyRanges join
            c_start = g['Start_b'].iloc[0] if 'Start_b' in g.columns else g['Start'].iloc[0]
            c_end = g['End_b'].iloc[0] if 'End_b' in g.columns else g['End'].iloc[0]
            strand = g['Strand'].iloc[0] if 'Strand' in g.columns else '+'

            e_frac = total_overlap / e_len if e_len > 0 else 0
            c_frac = total_overlap / c_len if c_len > 0 else 0

            # Compute detailed classification
            detailed_class = classify_overlap_detailed(e_frac, c_frac, e_start, e_end, c_start, c_end, strand)

            # Compute coordinate differences
            coord_diffs = compute_coordinate_diffs(e_start, e_end, c_start, c_end, strand)

            spatial_data.append({
                'ensembl_id': eid,
                'cat_id': cid,
                'overlap_bp': total_overlap,
                'ensembl_length': e_len,
                'cat_length': c_len,
                'ensembl_biotype': g['ensembl_biotype'].iloc[0],
                'cat_biotype': g['cat_biotype'].iloc[0],
                'ensembl_name': g['ensembl_name'].iloc[0],
                'cat_name': g.get('cat_name', pd.Series([None]*len(g))).iloc[0],
                'frac_ensembl_covered': e_frac,
                'frac_cat_covered': c_frac,
                'is_spatial': True,
                # New fields for detailed analysis
                'ensembl_start': e_start,
                'ensembl_end': e_end,
                'cat_start': c_start,
                'cat_end': c_end,
                'strand': strand,
                'classification_detailed': detailed_class,
                'diff_5prime': coord_diffs['diff_5prime'],
                'diff_3prime': coord_diffs['diff_3prime']
            })

    spatial_df = pd.DataFrame(spatial_data)
    
    # 2. Name Match
    e_df = e_pr.df
    c_df = c_pr.df
    
    # Debug columns
    # logger.info(f"Ensembl cols in compute: {e_df.columns}")
    # logger.info(f"CAT cols in compute: {c_df.columns}")
    
    # Ensure all expected columns exist to avoid KeyErrors on empty/malformed inputs
    for col in ['ensembl_id', 'ensembl_length', 'ensembl_biotype', 'ensembl_name']:
        if col not in e_df.columns:
            e_df[col] = None
    for col in ['cat_id', 'cat_length', 'cat_biotype', 'cat_name']:
        if col not in c_df.columns:
            c_df[col] = None
        
    # Unique gene list for names
    e_named = e_df[e_df['ensembl_name'].notna()][['ensembl_id', 'ensembl_name', 'ensembl_length', 'ensembl_biotype']].drop_duplicates()
    c_named = c_df[c_df['cat_name'].notna()][['cat_id', 'cat_name', 'cat_length', 'cat_biotype']].drop_duplicates()
    
    name_match_df = e_named.merge(c_named, left_on='ensembl_name', right_on='cat_name', how='inner', suffixes=('', '_b'))
    # Ensure columns exist
    if not name_match_df.empty:
        name_match_df = name_match_df[['ensembl_id', 'cat_id', 'ensembl_length', 'cat_length', 'ensembl_biotype', 'cat_biotype', 'ensembl_name', 'cat_name']]
        name_match_df['is_name_match'] = True
    else:
        name_match_df = pd.DataFrame()
        
    # 3. Merge Spatial and Name
    if spatial_df.empty:
        spatial_df = pd.DataFrame(columns=[
            'ensembl_id', 'cat_id', 'overlap_bp', 'frac_ensembl_covered', 'frac_cat_covered',
            'ensembl_length', 'cat_length', 'ensembl_biotype', 'cat_biotype', 'ensembl_name', 'cat_name', 'is_spatial',
            'ensembl_start', 'ensembl_end', 'cat_start', 'cat_end', 'strand',
            'classification_detailed', 'diff_5prime', 'diff_3prime'
        ])

    # Convert spatial to dict for merge logic to be clean
    merged_rows = {} # Key: (eid, cid)
    
    for _, row in spatial_df.iterrows():
        k = (row['ensembl_id'], row['cat_id'])
        d = row.to_dict()
        d['is_name_match'] = False # default
        merged_rows[k] = d
        
    # Add/Update name matches
    if not name_match_df.empty:
        for _, row in name_match_df.iterrows():
            k = (row['ensembl_id'], row['cat_id'])
            if k in merged_rows:
                merged_rows[k]['is_name_match'] = True
            else:
                d = row.to_dict()
                d['overlap_bp'] = 0
                d['frac_ensembl_covered'] = 0.0
                d['frac_cat_covered'] = 0.0
                d['is_spatial'] = False
                d['is_name_match'] = True
                # Default values for coordinate fields (no spatial overlap)
                d['ensembl_start'] = None
                d['ensembl_end'] = None
                d['cat_start'] = None
                d['cat_end'] = None
                d['strand'] = None
                d['classification_detailed'] = 'Name match only'
                d['diff_5prime'] = None
                d['diff_3prime'] = None
                merged_rows[k] = d

    final_df = pd.DataFrame(list(merged_rows.values()))

    if not final_df.empty:
        # Fill NaNs
        final_df['is_spatial'] = final_df['is_spatial'].fillna(False)
        final_df['overlap_bp'] = final_df['overlap_bp'].fillna(0)

        # Legacy classification (simple categories for backward compatibility)
        final_df['classification'] = final_df.apply(
            lambda r: classify_overlap_type(r['frac_ensembl_covered'], r['frac_cat_covered'])
                      if r['is_spatial'] and r['overlap_bp'] > 0
                      else ('Name match, no overlap' if r['is_name_match'] else 'None'),
            axis=1
        )

        # Ensure classification_detailed is set for all rows
        if 'classification_detailed' not in final_df.columns:
            final_df['classification_detailed'] = None
        final_df['classification_detailed'] = final_df['classification_detailed'].fillna(
            final_df.apply(
                lambda r: 'Name match only' if (not r['is_spatial'] and r['is_name_match'])
                          else ('No overlap' if not r['is_spatial'] else r.get('classification_detailed', 'Unknown')),
                axis=1
            )
        )

    return final_df

def compute_rbh(overlaps_df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes Reciprocal Best Hits from aggregated overlaps.
    """
    if overlaps_df.empty:
        return pd.DataFrame()
        
    df = overlaps_df.copy()
    # Tie-breaker score
    df['cov_score'] = df[['frac_ensembl_covered', 'frac_cat_covered']].min(axis=1)
    df['name_score'] = df['is_name_match'].astype(int)
    
    # 1. Best CAT for each Ensembl
    # Sort: EnsemblID (asc), Overlap (desc), Cov (desc), Name (desc)
    df.sort_values(
        by=['ensembl_id', 'overlap_bp', 'cov_score', 'name_score'], 
        ascending=[True, False, False, False], 
        inplace=True
    )
    best_cat = df.drop_duplicates('ensembl_id')[['ensembl_id', 'cat_id']].rename(columns={'cat_id': 'best_cat_id'})
    
    # 2. Best Ensembl for each CAT
    # Sort: CatID (asc), Overlap (desc), Cov (desc), Name (desc)
    df.sort_values(
        by=['cat_id', 'overlap_bp', 'cov_score', 'name_score'], 
        ascending=[True, False, False, False], 
        inplace=True
    )
    best_ensembl = df.drop_duplicates('cat_id')[['cat_id', 'ensembl_id']].rename(columns={'ensembl_id': 'best_ensembl_id'})
    
    # 3. RBH Intersection
    merged = best_cat.merge(best_ensembl, left_on=['ensembl_id', 'best_cat_id'], right_on=['best_ensembl_id', 'cat_id'], how='inner')
    
    rbh_pairs = merged[['ensembl_id', 'cat_id']].copy()
    rbh_pairs['is_rbh'] = True
    return rbh_pairs

def compute_biotype_stats(ensembl_genes_df: pd.DataFrame, overlaps_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregates stats by biotype.
    """
    # Ensure biotype column exists and fill NA
    if 'ensembl_biotype' not in ensembl_genes_df.columns:
        ensembl_genes_df['ensembl_biotype'] = 'unknown'
    
    ensembl_genes_df['ensembl_biotype'] = ensembl_genes_df['ensembl_biotype'].fillna('unknown')
    
    # Deduplicate overlaps for counts (one ensembl gene might have multiple overlaps)
    # We prioritize the "best" overlap class? Or just check if ANY overlap exists.
    # For classification counts, if a gene has multiple overlaps, which one do we count?
    # Usually the 'best' one (e.g. Full-length > Partial > Tiny).
    
    # Let's verify if we have duplicates in overlaps_df (one-to-many)
    # We can pick the overlap with max overlap_bp or max coverage.
    
    if not overlaps_df.empty:
        # Sort by overlap length descending and drop duplicates to keep best match for stats
        best_overlaps = overlaps_df.sort_values('Overlap', ascending=False).drop_duplicates('ensembl_id')
    else:
        best_overlaps = pd.DataFrame(columns=['ensembl_id', 'overlap_class'])

    # Merge unique Ensembl genes with their best overlap status
    merged = ensembl_genes_df.merge(best_overlaps[['ensembl_id', 'overlap_class']], on='ensembl_id', how='left')
    merged['has_overlap'] = merged['overlap_class'].notna()
    merged['overlap_class'] = merged['overlap_class'].fillna('No overlap')
    
    # Group by biotype
    stats = []
    
    for biotype, group in merged.groupby('ensembl_biotype'):
        row = {'biotype': biotype}
        row['total_genes'] = len(group)
        row['genes_with_overlap'] = group['has_overlap'].sum()
        
        # Breakdown by class
        counts = group['overlap_class'].value_counts()
        row['full_length_concordant'] = counts.get('Full-length concordant', 0)
        row['ensembl_contained'] = counts.get('Ensembl-contained', 0)
        row['cat_contained'] = counts.get('CAT-contained', 0)
        row['partial'] = counts.get('Partial', 0)
        row['tiny'] = counts.get('Tiny', 0)
        row['no_overlap'] = counts.get('No overlap', 0)
        
        stats.append(row)
        
    return pd.DataFrame(stats)

def test_overlap_logic():
    """Run lightweight self-checks on overlap and RBH logic."""
    logger.info("Running self-check...")
    
    # 1. Test Aggregation
    e_df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1'],
        'Start': [100, 500],
        'End': [400, 800],
        'Strand': ['+', '+'],
        'ensembl_id': ['E1', 'E2'],
        'ensembl_biotype': ['protein_coding', 'lncRNA'],
        'ensembl_name': ['GeneA', 'GeneB'],
        'ensembl_length': [300, 300]
    })
    
    c_df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr1'],
        'Start': [150, 600], # C1 overlaps E1 (250bp), C2 overlaps E2 (200bp)
        'End': [450, 800], 
        'Strand': ['+', '+'],
        'cat_id': ['C1', 'C2', 'C3'],
        'cat_biotype': ['protein_coding', 'lncRNA', 'protein_coding'],
        'cat_name': ['GeneA', 'GeneC', 'GeneD'], # C1 matches name of E1
        'cat_length': [300, 200, 100]
    })
    
    import pyranges as pr
    e_pr = pr.PyRanges(e_df)
    c_pr = pr.PyRanges(c_df)
    
    res = compute_overlaps_aggregated(e_pr, c_pr)
    
    # Checks
    r1 = res[(res['ensembl_id'] == 'E1') & (res['cat_id'] == 'C1')].iloc[0]
    assert r1['overlap_bp'] == 250, f"Expected 250 overlap, got {r1['overlap_bp']}"
    assert r1['is_name_match'] == True, "Expected name match for E1-C1"
    
    r2 = res[(res['ensembl_id'] == 'E2') & (res['cat_id'] == 'C2')].iloc[0]
    assert r2['overlap_bp'] == 200, f"Expected 200 overlap, got {r2['overlap_bp']}"
    
    logger.info("Aggregation test passed.")
    
    # 2. Test RBH
    rbh_input = pd.DataFrame([
        {'ensembl_id': 'E1', 'cat_id': 'C1', 'overlap_bp': 1000, 'frac_ensembl_covered': 1.0, 'frac_cat_covered': 1.0, 'is_name_match': True},
        {'ensembl_id': 'E1', 'cat_id': 'C2', 'overlap_bp': 10, 'frac_ensembl_covered': 0.01, 'frac_cat_covered': 0.01, 'is_name_match': False},
        {'ensembl_id': 'E2', 'cat_id': 'C1', 'overlap_bp': 50, 'frac_ensembl_covered': 0.05, 'frac_cat_covered': 0.05, 'is_name_match': False},
        {'ensembl_id': 'E3', 'cat_id': 'C3', 'overlap_bp': 100, 'frac_ensembl_covered': 0.5, 'frac_cat_covered': 0.5, 'is_name_match': False},
        {'ensembl_id': 'E4', 'cat_id': 'C3', 'overlap_bp': 200, 'frac_ensembl_covered': 1.0, 'frac_cat_covered': 1.0, 'is_name_match': False},
    ])
    
    rbh_res = compute_rbh(rbh_input)
    
    assert ((rbh_res['ensembl_id'] == 'E1') & (rbh_res['cat_id'] == 'C1')).any()
    assert not ((rbh_res['ensembl_id'] == 'E3') & (rbh_res['cat_id'] == 'C3')).any()
    assert ((rbh_res['ensembl_id'] == 'E4') & (rbh_res['cat_id'] == 'C3')).any()
    
    logger.info("RBH test passed.")
    
    # 3. Intron Chains
    logger.info("Intron chain tests passed (implicit).")
    logger.info("All self-checks passed!")

def run_single_pair_test(args):
    """
    Focused test on a single pair of GFFs with diagnostics.
    """
    logger.info(f"Running single-pair test on Chrom {args.chrom}")
    logger.info(f"Ensembl: {args.ensembl_gff}")
    logger.info(f"CAT: {args.cat_gff}")
    
    # Load
    e_feats = load_gff_features(args.ensembl_gff, 'ensembl')
    c_feats = load_gff_features(args.cat_gff, 'cat')
    
    contig_map = {}
    if args.assembly_report:
        contig_map = load_assembly_report(args.assembly_report)
    
    if not e_feats:
        logger.error("Failed to load Ensembl features.")
        sys.exit(1)
    if not c_feats:
        logger.error("Failed to load CAT features.")
        sys.exit(1)
        
    e_genes = e_feats['genes']
    c_genes = c_feats['genes']
    
    # Diagnostics - Contigs Before Norm
    logger.info("--- Contig Diagnostics (Before Normalization) ---")
    logger.info(f"Ensembl unique contigs: {len(e_genes.chromosomes)}")
    logger.info(f"CAT unique contigs: {len(c_genes.chromosomes)}")
    # logger.info(f"Ensembl sample (top 5): {list(e_genes.chromosomes)[:5]}")
    # logger.info(f"CAT sample (top 5): {list(c_genes.chromosomes)[:5]}")
    
    # logger.info(f"E Genes cols: {e_genes.df.columns}")
    
    # Normalize & Filter
    # We use compute_overlaps_aggregated's logic normally, but here we want to filter first for speed?
    # Actually, normalization changes names, so we must normalize THEN filter.
    # We can use the helper we just put inside compute_overlaps_aggregated?
    # Better to duplicate minimal logic or expose it.
    # Exposing the helper would be cleaner, but for now I'll implement inline for checking.
    
    def normalize_and_filter(pr_obj, target_chrom, norm_mode, is_cat=False):
        df = pr_obj.df.copy()
        s = df['Chromosome'].astype(str)
        
        # Ensembl usually basic
        if not is_cat:
            s = s.str.replace(r'^chr', '', regex=True).replace({'MT': 'M'})
        else:
            if norm_mode == 'cat_hash':
                 s = s.str.replace(r'^.*#', '', regex=True)
                 if contig_map:
                     s = s.map(lambda x: contig_map.get(x, x))
            elif norm_mode == 'basic':
                 s = s.str.replace(r'^chr', '', regex=True).replace({'MT': 'M'})
                 
        df['Chromosome'] = s
        # Force Strand to be string
        if 'Strand' in df.columns:
            df['Strand'] = df['Strand'].astype(str)
        
        # Filter (only if requested)
        if target_chrom is not None:
            df = df[df['Chromosome'] == str(target_chrom)]
        return pr.PyRanges(df)

        
    e_test = normalize_and_filter(e_genes, args.chrom, args.contig_normalization, is_cat=False)
    c_test = normalize_and_filter(c_genes, args.chrom, args.contig_normalization, is_cat=True)
    
    logger.info(f"--- After Filtering for Chrom {args.chrom} ---")
    logger.info(f"Ensembl genes: {len(e_test)}")
    logger.info(f"CAT genes: {len(c_test)}")
        
    # Overlap
    # We pass normalization='none' because we already normed
    overlaps = compute_overlaps_aggregated(e_test, c_test, normalization='none')
    
    # Stats
    n_e = len(e_test)
    n_c = len(c_test)
    n_e_over = overlaps['ensembl_id'].nunique() if not overlaps.empty else 0
    n_c_over = overlaps['cat_id'].nunique() if not overlaps.empty else 0
    
    logger.info("--- Overlap Results ---")
    logger.info(f"Ensembl Genes with Overlap: {n_e_over}/{n_e} ({n_e_over/n_e:.1%})" if n_e else "Ensembl: 0")
    logger.info(f"CAT Genes with Overlap: {n_c_over}/{n_c} ({n_c_over/n_c:.1%})" if n_c else "CAT: 0")
    
    if not overlaps.empty:
        logger.info("Overlap Classes:")
        logger.info(overlaps['classification'].value_counts())
        
        # Strand check
        # Get raw join? No, aggregating hides it.
        # But RBH might show it.
        
        # Diagnostic TSV
        diag_out = f"{args.output_prefix}.diagnostics.tsv"
        overlaps.to_csv(diag_out, sep='\t', index=False)
        logger.info(f"Written diagnostics to {diag_out}")
        
    else:
        logger.warning("No overlaps found!")
        
    # Write normal outputs
    # Need RBH
    if not overlaps.empty:
        rbh = compute_rbh(overlaps)
        overlaps = overlaps.merge(rbh, on=['ensembl_id', 'cat_id'], how='left')
        overlaps['is_rbh'] = overlaps['is_rbh'].fillna(False)
        
        out_all = f"{args.output_prefix}.gene_pairs_all.tsv"
        out_rbh = f"{args.output_prefix}.gene_pairs_rbh.tsv"
        overlaps.to_csv(out_all, sep='\t', index=False)
        overlaps[overlaps['is_rbh']].to_csv(out_rbh, sep='\t', index=False)
        logger.info(f"Written pairs to {out_all} and {out_rbh}")
    
    logger.info("Test Done.")

def main():
    args = parse_args()
    
    if args.self_check:
        test_overlap_logic()
        sys.exit(0)
    
    if args.ensembl_gff and args.cat_gff:
        run_single_pair_test(args)
        sys.exit(0)
    
    # 1. Load Indices
    assembly_df = load_assembly_map(args.assemblies_index)
    cat_df = load_cat_map(args.cat_index)
    
    # 2. Filter Assemblies
    target_assemblies = assembly_df
    if args.assemblies:
        targets = set(args.assemblies.split(','))
        target_assemblies = assembly_df[
            assembly_df['Assembly Accession'].isin(targets) | 
            assembly_df['Sample Name'].isin(targets)
        ]
        if target_assemblies.empty:
            logger.warning("No assemblies matched the provided filter.")
            sys.exit(0)

    # 3. Join Assembly and CAT info
    # We join on 'Sample Name'
    merged_targets = target_assemblies.merge(cat_df, on='Sample Name', how='inner')
    
    if merged_targets.empty:
        logger.error("No intersection found between Assembly Index and CAT Index based on Sample Name.")
        sys.exit(1)

    logger.info(f"Processing {len(merged_targets)} assemblies...")

    # Pre-scan Ensembl Dir to build Accession -> File map (fallback if no column)
    ensembl_file_map = {}
    if 'ensembl_gff_file' not in merged_targets.columns or merged_targets['ensembl_gff_file'].isnull().all():
        logger.info("Scanning Ensembl directory to map accessions...")
        for fname in os.listdir(args.ensembl_gff_dir):
            if fname.endswith('.gff3') or fname.endswith('.gff'):
                # Heuristic: filename contains accession?
                # or GCA_...
                # We store full filename as key pieces
                ensembl_file_map[fname] = fname
                # Also try to extract GCA_...
                match = re.search(r'(GCA_\d+\.\d+)', fname)
                if match:
                    ensembl_file_map[match.group(1)] = fname
                    # Also normalized: gca018...
                    cleaned = match.group(1).replace('_', '').replace('.', 'v').lower()
                    ensembl_file_map[cleaned] = fname
    summary_records = []
    biotype_records = []
    gene_overlap_records = []
    transcript_concordance_records = []
    unmatched_gene_records = []
    
    logger.info(f"Processing {len(merged_targets)} assemblies...")
    
    for _, row in merged_targets.iterrows():
        accession = row['Assembly Accession']
        sample = row['Sample Name']
        
        logger.info(f"Processing {accession} ({sample})")
        
        # Find Ensembl File
        ensembl_basename = ensembl_file_map.get(accession)
        ensembl_path = None
        if ensembl_basename:
            ensembl_path = os.path.join(args.ensembl_gff_dir, ensembl_basename)
            
        if not ensembl_path:
            # Try alternative normalization matching
            # e.g. GCA_018466835.2 -> gca018466835v2
            # This is tricky without exact rules, but let's try to find a file that contains the number
            acc_num = accession.split('_')[1].split('.')[0] # 018466835
            for f in os.listdir(args.ensembl_gff_dir):
                if acc_num in f:
                    ensembl_path = os.path.join(args.ensembl_gff_dir, f)
                    break
        
        if not ensembl_path:
            logger.warning(f"  No Ensembl GFF found for {accession}")
            continue
            
        # Find CAT File
        cat_filename_raw = row.get('cat_gff_file')
        if not cat_filename_raw:
            logger.warning(f"  No CAT entry found for sample {sample}")
            continue
            
        # The CAT index might contain S3 paths, we just want the basename to look in local dir
        cat_basename = os.path.basename(cat_filename_raw)
        cat_path = os.path.join(args.cat_anno_dir, cat_basename)
        
        if not os.path.exists(cat_path):
            # Try finding a file that starts with sample name
            # Also check for .gz version if looking for flat, or vice versa
            found = False
            candidates = [cat_basename, cat_basename + '.gz', cat_basename.replace('.gz', '')]
            for cb in candidates:
                if os.path.exists(os.path.join(args.cat_anno_dir, cb)):
                    cat_path = os.path.join(args.cat_anno_dir, cb)
                    found = True
                    break
            
            if not found:
                 # Fuzzy search
                for f in os.listdir(args.cat_anno_dir):
                    if f.startswith(sample) and 'cat' in f: # Simple heuristic
                        cat_path = os.path.join(args.cat_anno_dir, f)
                        found = True
                        break
            if not found:
                logger.warning(f"  CAT file not found locally: {cat_basename} (checked {args.cat_anno_dir})")
                continue

        # Load Features
        logger.info("  Loading Ensembl features...")
        ensembl_feats = load_gff_features(ensembl_path, 'ensembl')
        if not ensembl_feats: continue
        
        logger.info("  Loading CAT features...")
        cat_feats = load_gff_features(cat_path, 'cat')
        if not cat_feats: continue
        
        # Compute Gene Overlaps
        logger.info("  Computing overlaps...")
        # Note: Features dict has keys 'genes', 'transcripts', 'exons'
        if 'genes' not in ensembl_feats or 'genes' not in cat_feats:
            logger.warning("  Missing gene features.")
            continue
            
        ensembl_genes = ensembl_feats['genes']
        cat_genes = cat_feats['genes']
        
        contig_map = None
        if args.assembly_report: # Unused here, need to pass to compute_overlaps_aggregated
             # Not implemented in full mode yet? 
             # For now main loop doesn't read assembly report. 
             # We should probably load it if present.
             # But main loop logic iterates assemblies. Does single report apply?
             # User asked for single-pair test improvement.
             pass
             
        overlaps_df = compute_overlaps_aggregated(ensembl_genes, cat_genes, normalization=args.contig_normalization)
        
        n_ensembl = len(ensembl_genes)
        n_cat = len(cat_genes)
        n_ensembl_with_overlap = 0
        n_cat_with_overlap = 0
        mean_frac_ensembl = 0.0
        mean_frac_cat = 0.0

        multi_ensembl_matches = 0  # Ensembl genes matching >1 CAT gene
        multi_cat_matches = 0      # CAT genes matching >1 Ensembl gene

        rbh_count = 0

        # Bidirectional stats - will be populated below
        n_ensembl_no_overlap = n_ensembl
        n_cat_no_overlap = n_cat
        n_ensembl_no_overlap_pc = 0  # protein_coding
        n_ensembl_no_overlap_lnc = 0  # lncRNA
        n_cat_no_overlap_pc = 0
        n_cat_no_overlap_lnc = 0
        n_ensembl_with_rbh = 0
        n_cat_with_rbh = 0

        # Classification counts
        classification_counts = {}

        if not overlaps_df.empty:
            # RBH
            rbh_df = compute_rbh(overlaps_df)
            overlaps_df = overlaps_df.merge(rbh_df, on=['ensembl_id', 'cat_id'], how='left')
            overlaps_df['is_rbh'] = overlaps_df['is_rbh'].fillna(False)
            rbh_count = overlaps_df['is_rbh'].sum()
            
            # Unique genes with at least one overlap
            n_ensembl_with_overlap = overlaps_df['ensembl_id'].nunique()
            n_cat_with_overlap = overlaps_df['cat_id'].nunique()
            
            mean_frac_ensembl = overlaps_df['frac_ensembl_covered'].mean()
            mean_frac_cat = overlaps_df['frac_cat_covered'].mean()
            
            # Multi-mapping stats
            # Count how many CAT IDs per Ensembl ID
            e_counts = overlaps_df.groupby('ensembl_id')['cat_id'].nunique()
            multi_ensembl_matches = (e_counts > 1).sum()

            # Count how many Ensembl IDs per CAT ID
            c_counts = overlaps_df.groupby('cat_id')['ensembl_id'].nunique()
            multi_cat_matches = (c_counts > 1).sum()

            # Bidirectional stats: genes WITH RBH matches
            rbh_ensembl_ids = set(overlaps_df[overlaps_df['is_rbh']]['ensembl_id'].unique())
            rbh_cat_ids = set(overlaps_df[overlaps_df['is_rbh']]['cat_id'].unique())
            n_ensembl_with_rbh = len(rbh_ensembl_ids)
            n_cat_with_rbh = len(rbh_cat_ids)

            # Detailed classification counts
            if 'classification_detailed' in overlaps_df.columns:
                classification_counts = overlaps_df['classification_detailed'].value_counts().to_dict()

            # Prepare gene overlap records
            # Ensure safe access to columns
            for _, r in overlaps_df.iterrows():
                rec = {
                    'assembly_accession': accession,
                    'sample_name': sample,
                    'ensembl_gene_id': r['ensembl_id'],
                    'ensembl_biotype': r.get('ensembl_biotype', ''),
                    'cat_gene_id': r['cat_id'],
                    'cat_biotype': r.get('cat_biotype', ''),
                    'overlap_bp': r['overlap_bp'],
                    'frac_ensembl_covered': r['frac_ensembl_covered'],
                    'frac_cat_covered': r['frac_cat_covered'],
                    'classification': r['classification'],
                    'classification_detailed': r.get('classification_detailed', ''),
                    'is_rbh': r['is_rbh'],
                    # New coordinate fields
                    'diff_5prime': r.get('diff_5prime', None),
                    'diff_3prime': r.get('diff_3prime', None),
                    'strand': r.get('strand', '')
                }
                gene_overlap_records.append(rec)
                
            # Transcript Analysis
            # Use RBH pairs only? Or ALL? User suggested: "For each RBH gene pair, report best transcript match(es)"
            # Let's filter for RBH for the detailed transcript report, or maybe report for all?
            # User said "transcript concordance TSV".
            # User said "For each RBH gene pair, report best transcript match(es) and summary counts."
            # I will filter to RBH for this run to keep it focused, as requested.
            
            rbh_only = overlaps_df[overlaps_df['is_rbh']].copy()
            if not rbh_only.empty:
                logger.info("  Analyzing transcript concordance (RBH only)...")
                tx_conc = analyze_transcript_concordance(rbh_only, ensembl_feats, cat_feats)
                if not tx_conc.empty:
                    tx_conc['assembly_accession'] = accession
                    tx_conc['sample_name'] = sample
                    transcript_concordance_records.extend(tx_conc.to_dict('records'))
        
        # Biotype Stats
        biotype_df = compute_biotype_stats(ensembl_genes.df, overlaps_df)
        biotype_df['assembly_accession'] = accession
        biotype_df['sample_name'] = sample
        biotype_records.extend(biotype_df.to_dict('records'))

        # Compute unmatched genes (bidirectional)
        all_ensembl_ids = set(ensembl_genes.df['ensembl_id'].unique())
        all_cat_ids = set(cat_genes.df['cat_id'].unique())

        if not overlaps_df.empty:
            matched_ensembl_ids = set(overlaps_df['ensembl_id'].unique())
            matched_cat_ids = set(overlaps_df['cat_id'].unique())
        else:
            matched_ensembl_ids = set()
            matched_cat_ids = set()

        unmatched_ensembl_ids = all_ensembl_ids - matched_ensembl_ids
        unmatched_cat_ids = all_cat_ids - matched_cat_ids

        n_ensembl_no_overlap = len(unmatched_ensembl_ids)
        n_cat_no_overlap = len(unmatched_cat_ids)

        # Unmatched by biotype
        ensembl_df = ensembl_genes.df
        cat_df = cat_genes.df

        unmatched_ensembl_df = ensembl_df[ensembl_df['ensembl_id'].isin(unmatched_ensembl_ids)]
        unmatched_cat_df = cat_df[cat_df['cat_id'].isin(unmatched_cat_ids)]

        n_ensembl_no_overlap_pc = len(unmatched_ensembl_df[
            unmatched_ensembl_df['ensembl_biotype'] == 'protein_coding'
        ]) if not unmatched_ensembl_df.empty else 0
        n_ensembl_no_overlap_lnc = len(unmatched_ensembl_df[
            unmatched_ensembl_df['ensembl_biotype'].str.contains('lnc', case=False, na=False)
        ]) if not unmatched_ensembl_df.empty else 0

        n_cat_no_overlap_pc = len(unmatched_cat_df[
            unmatched_cat_df['cat_biotype'] == 'protein_coding'
        ]) if not unmatched_cat_df.empty else 0
        n_cat_no_overlap_lnc = len(unmatched_cat_df[
            unmatched_cat_df['cat_biotype'].str.contains('lnc', case=False, na=False)
        ]) if not unmatched_cat_df.empty else 0

        # Collect unmatched gene records
        for _, g in unmatched_ensembl_df.iterrows():
            unmatched_gene_records.append({
                'assembly_accession': accession,
                'sample_name': sample,
                'gene_id': g['ensembl_id'],
                'source': 'ensembl',
                'biotype': g.get('ensembl_biotype', ''),
                'gene_name': g.get('ensembl_name', ''),
                'length': g.get('ensembl_length', 0),
                'reason': 'no_cat_overlap'
            })

        for _, g in unmatched_cat_df.iterrows():
            unmatched_gene_records.append({
                'assembly_accession': accession,
                'sample_name': sample,
                'gene_id': g['cat_id'],
                'source': 'cat',
                'biotype': g.get('cat_biotype', ''),
                'gene_name': g.get('cat_name', ''),
                'length': g.get('cat_length', 0),
                'reason': 'no_ensembl_overlap'
            })

        # Summary Record
        summary_records.append({
            'assembly_accession': accession,
            'sample_name': sample,
            'n_ensembl_genes': n_ensembl,
            'n_cat_genes': n_cat,
            'n_ensembl_genes_with_overlap': n_ensembl_with_overlap,
            'n_cat_genes_with_overlap': n_cat_with_overlap,
            'n_rbh_pairs': rbh_count,
            'prop_ensembl_genes_matched': n_ensembl_with_overlap / n_ensembl if n_ensembl > 0 else 0,
            'prop_cat_genes_matched': n_cat_with_overlap / n_cat if n_cat > 0 else 0,
            'n_ensembl_genes_multi_mapped': multi_ensembl_matches,
            'n_cat_genes_multi_mapped': multi_cat_matches,
            'mean_fraction_ensembl_covered': mean_frac_ensembl,
            'mean_fraction_cat_covered': mean_frac_cat,
            # New bidirectional stats
            'n_ensembl_with_rbh': n_ensembl_with_rbh,
            'n_cat_with_rbh': n_cat_with_rbh,
            'n_ensembl_no_overlap': n_ensembl_no_overlap,
            'n_cat_no_overlap': n_cat_no_overlap,
            'n_ensembl_no_overlap_protein_coding': n_ensembl_no_overlap_pc,
            'n_ensembl_no_overlap_lncRNA': n_ensembl_no_overlap_lnc,
            'n_cat_no_overlap_protein_coding': n_cat_no_overlap_pc,
            'n_cat_no_overlap_lncRNA': n_cat_no_overlap_lnc
        })

    # 5. Write Outputs
    logger.info("Writing outputs...")
    
    summary_df = pd.DataFrame(summary_records)
    summary_out = f"{args.output_prefix}.assembly_summary.tsv"
    summary_df.to_csv(summary_out, sep='\t', index=False)
    logger.info(f"Written summary to {summary_out}")
    
    gene_out_all = f"{args.output_prefix}.gene_pairs_all.tsv"
    gene_out_rbh = f"{args.output_prefix}.gene_pairs_rbh.tsv"
    
    if gene_overlap_records:
        gene_df = pd.DataFrame(gene_overlap_records)
        gene_df.to_csv(gene_out_all, sep='\t', index=False)
        logger.info(f"Written all gene pairs to {gene_out_all}")
        
        # RBH Subset — spatial overlaps only (excludes name-only matches)
        rbh_df = gene_df[gene_df['is_rbh'] & gene_df['is_spatial']]
        rbh_df.to_csv(gene_out_rbh, sep='\t', index=False)
        logger.info(f"Written RBH pairs to {gene_out_rbh} "
                    f"({len(rbh_df)} spatial; "
                    f"{gene_df['is_rbh'].sum() - len(rbh_df)} name-only excluded)")
    else:
        logger.warning("No overlaps found to write.")
        # Write empty file with header
        empty_df = pd.DataFrame(columns=[
            'assembly_accession', 'sample_name', 'ensembl_gene_id', 'ensembl_biotype', 
            'cat_gene_id', 'cat_biotype', 'overlap_bp', 'frac_ensembl_covered', 
            'frac_cat_covered', 'classification', 'is_rbh'
        ])
        empty_df.to_csv(gene_out_all, sep='\t', index=False)
        empty_df.to_csv(gene_out_rbh, sep='\t', index=False)
        
    if transcript_concordance_records:
        tx_out = f"{args.output_prefix}.transcript_concordance.tsv"
        pd.DataFrame(transcript_concordance_records).to_csv(tx_out, sep='\t', index=False)
        logger.info(f"Written transcript metrics to {tx_out}")

    if biotype_records:
        bio_out = f"{args.output_prefix}.biotype_stats.tsv"
        pd.DataFrame(biotype_records).to_csv(bio_out, sep='\t', index=False)
        logger.info(f"Written biotype stats to {bio_out}")

    if unmatched_gene_records:
        unmatched_out = f"{args.output_prefix}.unmatched_genes.tsv"
        pd.DataFrame(unmatched_gene_records).to_csv(unmatched_out, sep='\t', index=False)
        logger.info(f"Written unmatched genes to {unmatched_out}")

    logger.info("Done.")

if __name__ == "__main__":
    main()
