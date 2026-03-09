#!/usr/bin/env python3
"""
Rename chromosome/contig names in an Ensembl GFF3 file to match GenBank accessions
using an NCBI assembly report.

Ensembl GFFs use Assigned-Molecule names (1, 2, X, MT) while CAT GFFs use
GenBank accessions (CM089370.1, CM089371.1). This script builds a reverse
mapping from the assembly report and rewrites the GFF3 seqid column.

Usage:
    rename_gff_chroms.py --gff input.gff3.gz --assembly-report report.txt --output renamed.gff3.gz
"""
import argparse
import gzip
import logging
import shutil
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_assembly_report(report_path: str) -> dict:
    """
    Parse NCBI assembly report to build reverse mapping:
        Assigned-Molecule -> GenBank-Accn

    Assembly report columns (tab-separated, # comment lines):
        0: Sequence-Name
        1: Sequence-Role
        2: Assigned-Molecule
        3: Assigned-Molecule-loc/type
        4: GenBank-Accn
        5: Relationship
        6: RefSeq-Accn
        7: Assembly-Unit
        8: Sequence-Length
        9: UCSC-style-name

    For chromosomes: Assigned-Molecule = '1', '2', 'X', etc.
    For scaffolds:   Assigned-Molecule = 'na' -> use GenBank-Accn as-is

    Returns dict: {assigned_molecule: genbank_accn}
    Only includes entries where assigned_molecule != genbank_accn (actual renames needed).
    """
    mapping = {}
    n_total = 0
    n_rename = 0

    with open(report_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            n_total += 1
            assigned_mol = parts[2]
            genbank_accn = parts[4]

            if genbank_accn == "na":
                continue

            # For scaffolds where assigned_mol is 'na', the seqid in the
            # Ensembl GFF is already the GenBank accession — no rename needed.
            if assigned_mol == "na":
                continue

            # Only add mapping if it's actually a rename
            if assigned_mol != genbank_accn:
                mapping[assigned_mol] = genbank_accn
                n_rename += 1

    logger.info(
        f"Assembly report: {n_total} sequences, {n_rename} need renaming "
        f"(e.g., chromosome names -> GenBank accessions)"
    )
    if mapping:
        examples = list(mapping.items())[:3]
        logger.info(f"  Sample mappings: {examples}")

    return mapping


def rename_gff(gff_path: str, mapping: dict, output_path: str) -> None:
    """
    Read a GFF3 file (possibly gzipped), rename seqids in column 1
    using the mapping, and write gzipped output.
    """
    n_lines = 0
    n_renamed = 0
    n_data_lines = 0
    renamed_seqids = set()

    # Determine if input is gzipped
    open_func = gzip.open if gff_path.endswith(".gz") else open
    open_mode = "rt" if gff_path.endswith(".gz") else "r"

    with open_func(gff_path, open_mode) as fin, gzip.open(output_path, "wt") as fout:
        for line in fin:
            n_lines += 1

            # Handle ##sequence-region pragma: ##sequence-region seqid start end
            if line.startswith("##sequence-region"):
                parts = line.strip().split()
                if len(parts) >= 2 and parts[1] in mapping:
                    parts[1] = mapping[parts[1]]
                    fout.write(" ".join(parts) + "\n")
                    n_renamed += 1
                    continue
                fout.write(line)
                continue

            # Pass through other comment/pragma lines unchanged
            if line.startswith("#"):
                fout.write(line)
                continue

            # Data line: tab-separated, column 0 is seqid
            parts = line.split("\t", 1)
            if len(parts) < 2:
                fout.write(line)
                continue

            n_data_lines += 1
            seqid = parts[0]
            if seqid in mapping:
                renamed_seqids.add(seqid)
                fout.write(mapping[seqid] + "\t" + parts[1])
                n_renamed += 1
            else:
                fout.write(line)

    logger.info(
        f"Processed {n_lines} lines ({n_data_lines} data lines), "
        f"renamed {n_renamed} lines"
    )
    if renamed_seqids:
        logger.info(f"  Renamed seqids: {sorted(renamed_seqids)}")


def main():
    parser = argparse.ArgumentParser(
        description="Rename Ensembl GFF3 chromosome names to GenBank accessions"
    )
    parser.add_argument(
        "--gff", required=True, help="Input GFF3 file (plain or .gz)"
    )
    parser.add_argument(
        "--assembly-report",
        required=True,
        help="NCBI assembly report file",
    )
    parser.add_argument(
        "--output", required=True, help="Output GFF3 file (.gz)"
    )
    args = parser.parse_args()

    # Check inputs exist
    if not Path(args.gff).exists():
        logger.error(f"GFF file not found: {args.gff}")
        sys.exit(1)

    report_path = Path(args.assembly_report)
    if not report_path.exists() or report_path.stat().st_size == 0:
        logger.warning(
            f"Assembly report missing or empty: {args.assembly_report}. "
            f"Copying GFF unchanged."
        )
        shutil.copy2(args.gff, args.output)
        return

    # Check if report is just a placeholder comment
    with open(args.assembly_report, "r") as f:
        first_line = f.readline().strip()
    if first_line.startswith("# Assembly report not found"):
        logger.warning(
            f"Assembly report is a placeholder: {first_line}. "
            f"Copying GFF unchanged."
        )
        shutil.copy2(args.gff, args.output)
        return

    # Parse assembly report
    mapping = parse_assembly_report(args.assembly_report)

    if not mapping:
        logger.info("No renames needed (all seqids already match). Copying GFF unchanged.")
        shutil.copy2(args.gff, args.output)
        return

    # Rename and write
    rename_gff(args.gff, mapping, args.output)
    logger.info(f"Written renamed GFF to {args.output}")


if __name__ == "__main__":
    main()
