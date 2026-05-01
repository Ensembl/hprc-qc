#!/usr/bin/env python3
"""Filter GFF3 records by gene/transcript biotype.

The main use case is sensitivity analysis for CAT annotations: remove
`unknown_likely_coding` gene/transcript models, plus their child features,
before re-running the normal Ensembl-vs-CAT comparison workflow.
"""

import argparse
import gzip
from pathlib import Path


TRANSCRIPT_TYPES = {
    "transcript",
    "mRNA",
    "lnc_RNA",
    "ncRNA",
    "miRNA",
    "snoRNA",
    "snRNA",
    "tRNA",
    "rRNA",
    "pseudogenic_transcript",
    "antisense_RNA",
    "guide_RNA",
    "scRNA",
    "vault_RNA",
    "Y_RNA",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove GFF3 gene/transcript models by biotype."
    )
    parser.add_argument("--input-gff", required=True, help="Input GFF3(.gz)")
    parser.add_argument("--output-gff", required=True, help="Filtered output GFF3(.gz)")
    parser.add_argument(
        "--exclude-biotype",
        action="append",
        required=True,
        help="Biotype to remove. May be supplied multiple times.",
    )
    return parser.parse_args()


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_attrs(attr_string):
    attrs = {}
    for item in attr_string.strip().split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def attr_biotype(attrs):
    return (
        attrs.get("biotype")
        or attrs.get("gene_biotype")
        or attrs.get("transcript_biotype")
        or ""
    ).strip()


def feature_id(attrs):
    return attrs.get("ID") or attrs.get("gene_id") or attrs.get("transcript_id")


def parent_ids(attrs):
    raw = attrs.get("Parent") or attrs.get("gene_id") or attrs.get("transcript_id") or ""
    return [x for x in str(raw).split(",") if x]


def is_gene_type(feature_type):
    return feature_type == "gene" or feature_type.endswith("gene")


def main():
    args = parse_args()
    input_gff = Path(args.input_gff)
    output_gff = Path(args.output_gff)
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    excluded_biotypes = set(args.exclude_biotype)

    skip_genes = set()
    skip_transcripts = set()

    # Pass 1: identify explicitly excluded genes/transcripts.
    with open_text(input_gff, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            attrs = parse_attrs(parts[8])
            if attr_biotype(attrs) not in excluded_biotypes:
                continue

            fid = feature_id(attrs)
            if is_gene_type(feature_type):
                if fid:
                    skip_genes.add(fid)
            elif feature_type in TRANSCRIPT_TYPES:
                if fid:
                    skip_transcripts.add(fid)
                skip_genes.update(parent_ids(attrs))

    # Pass 2: any transcript under a skipped gene should be skipped too.
    with open_text(input_gff, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in TRANSCRIPT_TYPES:
                continue

            attrs = parse_attrs(parts[8])
            if any(parent in skip_genes for parent in parent_ids(attrs)):
                fid = feature_id(attrs)
                if fid:
                    skip_transcripts.add(fid)

    kept = 0
    skipped = 0
    with open_text(input_gff, "rt") as in_handle, open_text(output_gff, "wt") as out_handle:
        for line in in_handle:
            if not line or line.startswith("#"):
                out_handle.write(line)
                kept += 1
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                out_handle.write(line)
                kept += 1
                continue

            attrs = parse_attrs(parts[8])
            fid = feature_id(attrs)
            parents = parent_ids(attrs)
            remove = (
                fid in skip_genes
                or fid in skip_transcripts
                or any(parent in skip_genes or parent in skip_transcripts for parent in parents)
            )
            if remove:
                skipped += 1
                continue

            out_handle.write(line)
            kept += 1

    print(
        f"Filtered {input_gff}: excluded biotypes={sorted(excluded_biotypes)} "
        f"skip_genes={len(skip_genes)} skip_transcripts={len(skip_transcripts)} "
        f"kept_records={kept} skipped_records={skipped}"
    )


if __name__ == "__main__":
    main()
