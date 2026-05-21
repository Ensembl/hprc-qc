#!/usr/bin/env python3
"""Filter GFF3 gene/transcript models by configurable predicates.

The original use case was removing CAT `unknown_likely_coding` models before
rerunning the normal Ensembl-vs-CAT comparison workflow.  The script now keeps
that behaviour while supporting additional sensitivity filters, such as
hyphenated/readthrough-like gene names, collapsed source genes, extra paralogs,
specific transcript modes, sources, or arbitrary attribute predicates.

Filtering is model-aware:
  * a matching gene removes the gene and all child transcripts/features;
  * a matching transcript removes that transcript and all child features;
  * a gene whose transcripts are all removed is also removed.

An audit TSV is written so filtered runs can be interpreted without guessing
what was removed.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import re
from pathlib import Path
from typing import Iterable


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

EMPTY_VALUES = {"", "N/A", "NA", "nan", "None", "."}
NAME_ATTRS = ("gene_name", "Name", "source_gene_common_name")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove GFF3 gene/transcript models by configurable predicates."
    )
    parser.add_argument("--input-gff", required=True, help="Input GFF3(.gz)")
    parser.add_argument("--output-gff", required=True, help="Filtered output GFF3(.gz)")
    parser.add_argument(
        "--audit-tsv",
        required=True,
        help="Output TSV summarising removed models and records.",
    )
    parser.add_argument(
        "--filter-name",
        default="filtered_cat",
        help="Human-readable filter name to include in audit output.",
    )
    parser.add_argument(
        "--exclude-biotype",
        action="append",
        default=[],
        help="Gene/transcript biotype to remove. May be supplied multiple times.",
    )
    parser.add_argument(
        "--exclude-source",
        action="append",
        default=[],
        help="GFF source column value to remove, e.g. Liftoff. May be supplied multiple times.",
    )
    parser.add_argument(
        "--exclude-transcript-mode",
        action="append",
        default=[],
        help="Value in transcript_modes to remove. May be supplied multiple times.",
    )
    parser.add_argument(
        "--exclude-gene-name-regex",
        action="append",
        default=[],
        help=(
            "Regex matched against gene_name, Name and source_gene_common_name. "
            "May be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--exclude-gene-name-file",
        action="append",
        default=[],
        help=(
            "Remove models where gene_name, Name or source_gene_common_name matches "
            "one value from a newline-delimited file. Blank lines and # comments are ignored."
        ),
    )
    parser.add_argument(
        "--exclude-hyphenated-gene-names",
        action="store_true",
        help=(
            "Remove models whose gene_name/Name/source_gene_common_name contains a hyphen. "
            "This matches the broad readthrough/hyphenated source category; it is not a "
            "pure readthrough-only biological definition."
        ),
    )
    parser.add_argument(
        "--exclude-attr-nonempty",
        action="append",
        default=[],
        help=(
            "Remove models where this attribute is present and not N/A/nan/empty. "
            "May be supplied multiple times, e.g. collapsed_gene_ids."
        ),
    )
    parser.add_argument(
        "--exclude-attr-equals",
        action="append",
        default=[],
        help=(
            "Remove models where KEY=VALUE, formatted as KEY=VALUE. "
            "May be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--exclude-attr-in",
        action="append",
        default=[],
        help=(
            "Remove models where KEY has one of VALUEs, formatted as KEY=VALUE1|VALUE2. "
            "For comma-delimited GFF attributes, any listed token can match."
        ),
    )
    parser.add_argument(
        "--exclude-attr-in-file",
        action="append",
        default=[],
        help=(
            "Remove models where KEY has one value from FILE, formatted as KEY=FILE. "
            "For comma-delimited GFF attributes, any listed token can match. "
            "May be supplied multiple times."
        ),
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


def attr_tokens(value: str | None) -> set[str]:
    if value is None:
        return set()
    return {x.strip() for x in str(value).split(",") if x.strip()}


def attr_values(value: str | None) -> set[str]:
    if value is None:
        return set()
    raw = str(value).strip()
    if not raw:
        return set()
    return {raw} | attr_tokens(raw)


def attr_nonempty(attrs, key):
    value = attrs.get(key)
    return value is not None and str(value).strip() not in EMPTY_VALUES


def feature_id(attrs):
    return attrs.get("ID") or attrs.get("gene_id") or attrs.get("transcript_id")


def gene_id(attrs):
    return attrs.get("gene_id")


def parent_ids(attrs):
    raw = attrs.get("Parent") or attrs.get("gene_id") or attrs.get("transcript_id") or ""
    return [x for x in str(raw).split(",") if x]


def is_gene_type(feature_type):
    return feature_type == "gene" or feature_type.endswith("gene")


def model_names(attrs) -> list[str]:
    return [attrs[k] for k in NAME_ATTRS if attrs.get(k)]


def parse_attr_equals(specs: Iterable[str]) -> list[tuple[str, str]]:
    parsed = []
    for spec in specs:
        if "=" not in spec:
            raise SystemExit(f"--exclude-attr-equals expects KEY=VALUE, got: {spec}")
        key, value = spec.split("=", 1)
        parsed.append((key.strip(), value.strip()))
    return parsed


def parse_attr_in(specs: Iterable[str]) -> list[tuple[str, set[str]]]:
    parsed = []
    for spec in specs:
        if "=" in spec:
            key, values = spec.split("=", 1)
        elif ":" in spec:
            key, values = spec.split(":", 1)
        else:
            raise SystemExit(f"--exclude-attr-in expects KEY=V1|V2, got: {spec}")
        parsed.append((key.strip(), {v.strip() for v in values.split("|") if v.strip()}))
    return parsed


def load_value_file(path: str) -> set[str]:
    values = set()
    with open_text(path, "rt") as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            values.add(value)
    return values


def load_value_files(paths: Iterable[str]) -> set[str]:
    values = set()
    for path in paths:
        values.update(load_value_file(path))
    return values


def parse_attr_in_files(specs: Iterable[str]) -> list[tuple[str, str, set[str]]]:
    parsed = []
    for spec in specs:
        if "=" in spec:
            key, path = spec.split("=", 1)
        elif ":" in spec:
            key, path = spec.split(":", 1)
        else:
            raise SystemExit(f"--exclude-attr-in-file expects KEY=FILE, got: {spec}")
        key = key.strip()
        path = path.strip()
        if not key or not path:
            raise SystemExit(f"--exclude-attr-in-file expects KEY=FILE, got: {spec}")
        parsed.append((key, Path(path).name, load_value_file(path)))
    return parsed


class PredicateSet:
    def __init__(self, args):
        self.excluded_biotypes = set(args.exclude_biotype)
        self.excluded_sources = set(args.exclude_source)
        self.excluded_transcript_modes = set(args.exclude_transcript_mode)
        self.gene_name_regexes = [re.compile(p) for p in args.exclude_gene_name_regex]
        self.gene_name_file_values = load_value_files(args.exclude_gene_name_file)
        self.exclude_hyphenated_gene_names = bool(args.exclude_hyphenated_gene_names)
        self.attr_nonempty_keys = list(args.exclude_attr_nonempty)
        self.attr_equals = parse_attr_equals(args.exclude_attr_equals)
        self.attr_in = parse_attr_in(args.exclude_attr_in)
        self.attr_in_files = parse_attr_in_files(args.exclude_attr_in_file)

    def reasons(self, source, feature_type, attrs) -> set[str]:
        reasons = set()
        biotype = attr_biotype(attrs)
        if biotype in self.excluded_biotypes:
            reasons.add(f"biotype:{biotype}")

        if source in self.excluded_sources:
            reasons.add(f"source:{source}")

        modes = attr_tokens(attrs.get("transcript_modes"))
        for mode in sorted(modes & self.excluded_transcript_modes):
            reasons.add(f"transcript_mode:{mode}")

        names = model_names(attrs)
        for rx in self.gene_name_regexes:
            if any(rx.search(name) for name in names):
                reasons.add(f"gene_name_regex:{rx.pattern}")

        if self.gene_name_file_values and any(name in self.gene_name_file_values for name in names):
            reasons.add("gene_name_file")

        if self.exclude_hyphenated_gene_names and any("-" in name for name in names):
            reasons.add("hyphenated_gene_name")

        for key in self.attr_nonempty_keys:
            if attr_nonempty(attrs, key):
                reasons.add(f"attr_nonempty:{key}")

        for key, value in self.attr_equals:
            if attrs.get(key) == value:
                reasons.add(f"attr_equals:{key}={value}")

        for key, values in self.attr_in:
            if attr_values(attrs.get(key)) & values:
                reasons.add(f"attr_in:{key}={'|'.join(sorted(values))}")

        for key, label, values in self.attr_in_files:
            if attr_values(attrs.get(key)) & values:
                reasons.add(f"attr_in_file:{key}:{label}")

        return reasons


def model_summary(feature_type, source, attrs):
    return {
        "feature_type": feature_type,
        "source": source,
        "gene_biotype": attrs.get("gene_biotype") or attrs.get("biotype") or "",
        "transcript_biotype": attrs.get("transcript_biotype") or attrs.get("biotype") or "",
        "transcript_modes": attrs.get("transcript_modes") or "",
        "transcript_class": attrs.get("transcript_class") or "",
        "gene_name": attrs.get("gene_name") or attrs.get("Name") or attrs.get("source_gene_common_name") or "",
    }


def audit_key(section, reasons, summary):
    return (
        section,
        ",".join(sorted(reasons)),
        summary["feature_type"],
        summary["source"],
        summary["gene_biotype"],
        summary["transcript_biotype"],
        summary["transcript_modes"],
        summary["transcript_class"],
    )


def write_audit(path, filter_name, audit_counts):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "filter_name",
                "section",
                "reasons",
                "feature_type",
                "source",
                "gene_biotype",
                "transcript_biotype",
                "transcript_modes",
                "transcript_class",
                "n_records_or_models",
            ]
        )
        for key, count in sorted(audit_counts.items()):
            section, reasons, feature_type, source, gene_bt, tx_bt, modes, tx_class = key
            writer.writerow(
                [
                    filter_name,
                    section,
                    reasons,
                    feature_type,
                    source,
                    gene_bt,
                    tx_bt,
                    modes,
                    tx_class,
                    count,
                ]
            )


def main():
    args = parse_args()
    input_gff = Path(args.input_gff)
    output_gff = Path(args.output_gff)
    audit_tsv = Path(args.audit_tsv)
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    audit_tsv.parent.mkdir(parents=True, exist_ok=True)

    predicates = PredicateSet(args)
    skip_genes: dict[str, set[str]] = collections.defaultdict(set)
    skip_transcripts: dict[str, set[str]] = collections.defaultdict(set)
    gene_to_transcripts: dict[str, set[str]] = collections.defaultdict(set)
    gene_summaries = {}
    transcript_summaries = {}

    # Pass 1: identify genes/transcripts directly matching predicates, and
    # collect gene->transcript membership.
    with open_text(input_gff, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            source = parts[1]
            feature_type = parts[2]
            attrs = parse_attrs(parts[8])
            fid = feature_id(attrs)
            reasons = predicates.reasons(source, feature_type, attrs)

            if is_gene_type(feature_type):
                gid = fid or gene_id(attrs)
                if gid:
                    gene_summaries[gid] = model_summary(feature_type, source, attrs)
                    if reasons:
                        skip_genes[gid].update(reasons)

            elif feature_type in TRANSCRIPT_TYPES:
                tid = fid or attrs.get("transcript_id")
                parents = parent_ids(attrs)
                for parent in parents:
                    gene_to_transcripts[parent].add(tid)
                if tid:
                    transcript_summaries[tid] = model_summary(feature_type, source, attrs)
                    if reasons:
                        skip_transcripts[tid].update(reasons)

    # Pass 2: transcripts under skipped genes should be skipped, and genes whose
    # complete transcript set is skipped should also be removed.
    for gene, transcripts in gene_to_transcripts.items():
        if gene in skip_genes:
            for tx in transcripts:
                skip_transcripts[tx].update({f"parent_gene:{r}" for r in skip_genes[gene]})

    changed = True
    while changed:
        changed = False
        for gene, transcripts in gene_to_transcripts.items():
            if gene in skip_genes or not transcripts:
                continue
            if transcripts and all(tx in skip_transcripts for tx in transcripts):
                skip_genes[gene].add("all_child_transcripts_removed")
                changed = True

    # Audit removed models before record-level output.
    audit_counts = collections.Counter()
    for gid, reasons in skip_genes.items():
        summary = gene_summaries.get(
            gid,
            {
                "feature_type": "gene",
                "source": "",
                "gene_biotype": "",
                "transcript_biotype": "",
                "transcript_modes": "",
                "transcript_class": "",
                "gene_name": "",
            },
        )
        audit_counts[audit_key("removed_model", reasons, summary)] += 1
    for tid, reasons in skip_transcripts.items():
        summary = transcript_summaries.get(
            tid,
            {
                "feature_type": "transcript",
                "source": "",
                "gene_biotype": "",
                "transcript_biotype": "",
                "transcript_modes": "",
                "transcript_class": "",
                "gene_name": "",
            },
        )
        audit_counts[audit_key("removed_model", reasons, summary)] += 1

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

            source = parts[1]
            feature_type = parts[2]
            attrs = parse_attrs(parts[8])
            fid = feature_id(attrs)
            parents = parent_ids(attrs)
            reasons = set()

            if fid in skip_genes:
                reasons.update(skip_genes[fid])
            if fid in skip_transcripts:
                reasons.update(skip_transcripts[fid])
            for parent in parents:
                if parent in skip_genes:
                    reasons.update({f"parent_gene:{r}" for r in skip_genes[parent]})
                if parent in skip_transcripts:
                    reasons.update({f"parent_transcript:{r}" for r in skip_transcripts[parent]})

            if reasons:
                skipped += 1
                audit_counts[audit_key("removed_record", reasons, model_summary(feature_type, source, attrs))] += 1
                continue

            out_handle.write(line)
            kept += 1

    write_audit(audit_tsv, args.filter_name, audit_counts)
    print(
        f"Filtered {input_gff}: filter_name={args.filter_name} "
        f"skip_genes={len(skip_genes)} skip_transcripts={len(skip_transcripts)} "
        f"kept_records={kept} skipped_records={skipped} audit={audit_tsv}"
    )


if __name__ == "__main__":
    main()
