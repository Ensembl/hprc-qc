#!/usr/bin/env python3
"""
Convert HPRC metadata CSV to pipeline sample sheet format.

Supports multiple input formats:
1. HPRC_2 - all assemblies.csv format
2. HPRC_2 - original.csv format

Output format:
    assembly_accession,sample_name
    GCA_041900255.1,HG00408_pat
    GCA_041900245.1,HG00408_mat
"""
import csv
import sys

if len(sys.argv) < 2:
    print("Usage: convert_hprc_metadata_to_samplesheet.py <input_csv>", file=sys.stderr)
    sys.exit(1)

print("assembly_accession,sample_name")

with open(sys.argv[1], 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Try different column name formats
        # Format 1: "HPRC_2 - all assemblies.csv" (first column is sample ID, may have leading space)
        sample_id = row.get(' ', '').strip() or row.get('', '').strip() or row.get('sample_id', '').strip()
        haplotype_num = row.get('haplotype', '').strip()
        accession = row.get('genbank_accession', '').strip() or row.get('accession', '').strip()
        hap_text = row.get('paternal/ maternal', '').strip()

        # Skip if no sample ID or accession
        if not sample_id or not accession:
            continue

        # Skip #N/A accessions or invalid accessions
        if accession == '#N/A' or not accession.startswith('GCA_'):
            continue

        # Determine suffix from various possible sources
        suffix = None

        # Try haplotype number first (1 = pat, 2 = mat)
        if haplotype_num == '1':
            suffix = 'pat'
        elif haplotype_num == '2':
            suffix = 'mat'
        # Try paternal/maternal text field
        elif hap_text and hap_text.lower().startswith('pat'):
            suffix = 'pat'
        elif hap_text and hap_text.lower().startswith('mat'):
            suffix = 'mat'
        # Try assembly_name field (e.g., "NA18879_hap1_...")
        else:
            assembly_name = row.get('assembly_name', '').strip()
            if 'hap1' in assembly_name.lower() or '_1_' in assembly_name or assembly_name.endswith('_1'):
                suffix = 'pat'
            elif 'hap2' in assembly_name.lower() or '_2_' in assembly_name or assembly_name.endswith('_2'):
                suffix = 'mat'

        if not suffix:
            continue

        sample_name = f"{sample_id}_{suffix}"
        print(f"{accession},{sample_name}")
