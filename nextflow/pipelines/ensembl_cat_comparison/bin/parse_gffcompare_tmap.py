#!/usr/bin/env python3
import argparse
import csv
import os
from collections import Counter

HEADER = [
    "assembly_accession",
    "sample_name",
    "direction",
    "class_code",
    "n_transcripts",
    "denominator",
    "pct",
]

def main():
    ap = argparse.ArgumentParser(description="Parse gffcompare .tmap to class_code counts with per-direction denominator.")
    ap.add_argument("--tmap", required=True, help="Path to .tmap file")
    ap.add_argument("--assembly-accession", required=True)
    ap.add_argument("--sample-name", required=True)
    ap.add_argument("--direction", required=True, choices=["Ensembl_to_CAT", "CAT_to_Ensembl"])
    ap.add_argument("--output", required=True, help="Output TSV path")
    args = ap.parse_args()

    # Handle missing/empty input: header-only output
    if not os.path.exists(args.tmap) or os.path.getsize(args.tmap) == 0:
        with open(args.output, "w", newline="") as outfh:
            outfh.write("\t".join(HEADER) + "\n")
        return

    # Read header and rows
    with open(args.tmap, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            header = []
            rows = []
        else:
            rows = list(reader)

    # Identify class_code column; default to empty code if not found
    try:
        cidx = header.index("class_code")
    except ValueError:
        cidx = None

    # Count class codes and compute denominator
    denom = 0
    counts = Counter()
    for row in rows:
        if not row:
            continue
        denom += 1
        code = row[cidx] if cidx is not None and cidx < len(row) else ""
        counts[code] += 1

    # Write output
    with open(args.output, "w", newline="") as outfh:
        outfh.write("\t".join(HEADER) + "\n")
        if denom == 0:
            return
        for code in sorted(counts.keys()):
            n = counts[code]
            pct = 100.0 * n / denom if denom else 0.0
            outfh.write("\t".join([
                args.assembly_accession,
                args.sample_name,
                args.direction,
                str(code),
                str(n),
                str(denom),
                f"{pct:.4f}",
            ]) + "\n")

if __name__ == "__main__":
    main()
