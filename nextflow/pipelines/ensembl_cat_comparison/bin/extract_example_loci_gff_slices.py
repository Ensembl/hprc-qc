#!/usr/bin/env python3
"""
Extract CAT and Ensembl GFF feature slices for selected example loci.

This is intended for the cluster, where the full CAT and Ensembl GFF collections
are available. The manifest should contain one row per example locus with at
least `assembly_accession`, `sample_name`, `ensembl_gene_id`, `cat_gene_id`,
`gene_name`, and `example_type` columns.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--ensembl-gff-dir", required=True, type=Path)
    parser.add_argument("--cat-gff-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional maximum manifest rows to process; 0 means all rows.",
    )
    return parser.parse_args()


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def parse_attrs(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in attr.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
    return out


def find_one(root: Path, patterns: list[str]) -> Path | None:
    for pattern in patterns:
        hits = sorted(root.rglob(pattern))
        if hits:
            return hits[0]
    return None


def find_ensembl_gff(root: Path, accession: str) -> Path | None:
    compact = accession.replace("_", "").replace(".", "").lower()
    versioned = accession.lower().replace("gca_", "gca").replace(".", "v")
    return find_one(
        root,
        [
            f"*{accession}*.gff3*",
            f"*{versioned}*.gff3*",
            f"*{compact}*.gff3*",
        ],
    )


def find_cat_gff(root: Path, sample_name: str) -> Path | None:
    return find_one(
        root,
        [
            f"*{sample_name}*cat*.gff3*",
            f"*{sample_name}*.gff3*",
        ],
    )


def normalise_id(value: str) -> str:
    return value.replace("gene:", "").replace("transcript:", "")


def extract_gene_model(gff_path: Path, gene_id: str) -> list[str]:
    """Return gene, transcript, exon, CDS, and UTR lines for one GFF gene ID."""
    target = normalise_id(str(gene_id))
    if not target:
        return []

    gene_lines: list[str] = []
    transcript_ids: set[str] = set()
    child_lines: list[str] = []

    with open_maybe_gzip(gff_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attrs(fields[8])
            feature_id = normalise_id(attrs.get("ID", ""))
            parent = normalise_id(attrs.get("Parent", ""))

            if feature_id == target:
                gene_lines.append(line.rstrip("\n"))
            elif parent == target:
                transcript_ids.add(feature_id)
                child_lines.append(line.rstrip("\n"))

    with open_maybe_gzip(gff_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attrs(fields[8])
            parent_ids = {normalise_id(x) for x in attrs.get("Parent", "").split(",")}
            if transcript_ids & parent_ids:
                child_lines.append(line.rstrip("\n"))

    return gene_lines + child_lines


def safe_name(text: str) -> str:
    keep = []
    for char in str(text):
        keep.append(char if char.isalnum() or char in "._-" else "_")
    return "".join(keep)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    if args.limit:
        manifest = manifest.head(args.limit)

    status_rows: list[dict[str, object]] = []
    for idx, row in manifest.iterrows():
        prefix = f"{idx:03d}_{safe_name(row['example_type'])}_{safe_name(row['gene_name'])}_{safe_name(row['assembly_accession'])}"
        ensembl_gff = find_ensembl_gff(args.ensembl_gff_dir, row["assembly_accession"])
        cat_gff = find_cat_gff(args.cat_gff_dir, row["sample_name"])

        cat_lines = extract_gene_model(cat_gff, row["cat_gene_id"]) if cat_gff else []
        ens_lines = extract_gene_model(ensembl_gff, row["ensembl_gene_id"]) if ensembl_gff else []

        cat_out = args.output_dir / f"{prefix}.CAT.gff3"
        ens_out = args.output_dir / f"{prefix}.Ensembl.gff3"
        cat_out.write_text("\n".join(cat_lines) + ("\n" if cat_lines else ""))
        ens_out.write_text("\n".join(ens_lines) + ("\n" if ens_lines else ""))

        status_rows.append(
            {
                **row.to_dict(),
                "ensembl_gff_found": bool(ensembl_gff),
                "cat_gff_found": bool(cat_gff),
                "ensembl_feature_lines": len(ens_lines),
                "cat_feature_lines": len(cat_lines),
                "ensembl_slice": str(ens_out),
                "cat_slice": str(cat_out),
            }
        )

    pd.DataFrame(status_rows).to_csv(args.output_dir / "extract_status.tsv", sep="\t", index=False)
    print(f"Wrote GFF slices and status to: {args.output_dir}")


if __name__ == "__main__":
    main()
