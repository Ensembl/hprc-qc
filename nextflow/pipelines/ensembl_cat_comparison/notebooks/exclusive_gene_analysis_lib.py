from __future__ import annotations

import os
import re
import shlex
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_OUTPUT_DIR = Path(
    "/hps/nobackup/flicek/ensembl/genebuild/jackt/hprc/hprc-qc/results"
)
CATALOG_NAME_CLASSES = {"standard_symbol", "orf_symbol", "hyphenated_symbol"}

_CLONE_LOCUS_RE = re.compile(r"^(AC|AL|AP|BX|CT|FP|FO|CR|AF|AD)\d")
_ORF_RE = re.compile(r"^C\d+orf")


def classify_gene_name(name: object) -> str:
    if pd.isna(name):
        return "missing"
    text = str(name)
    if text.startswith("gene:ENSG") or text.startswith("ENSG"):
        return "ENSG_fallback"
    if text.startswith("MSTRG."):
        return "MSTRG"
    if text.startswith("CHM13_G"):
        return "CHM13_gene_id"
    if _CLONE_LOCUS_RE.match(text):
        return "clone_locus"
    if _ORF_RE.match(text):
        return "orf_symbol"
    if "-" in text:
        return "hyphenated_symbol"
    return "standard_symbol"


def catalog_like_mask(values: pd.Series) -> pd.Series:
    return values.map(classify_gene_name).isin(CATALOG_NAME_CLASSES)


def first_nonempty(series: pd.Series) -> str:
    cleaned = series.dropna().astype(str)
    cleaned = cleaned[cleaned != ""]
    return cleaned.iloc[0] if len(cleaned) else ""


def parse_gff_attributes(attr_string: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def _parse_gff3_gene_lines(gff_path: Path) -> list[str]:
    """Return only gene-feature lines from a gzipped GFF3 using zcat+awk (far faster than Python gzip)."""
    result = subprocess.run(
        f"zcat {shlex.quote(str(gff_path))} | awk -F'\\t' 'NF==9 && $3==\"gene\"'",
        shell=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.splitlines()


def resolve_paths(
    output_dir: str | Path | None = None,
    cat_cache_dir: str | Path | None = None,
    ensembl_cache_dir: str | Path | None = None,
) -> dict[str, Path]:
    output_dir = Path(os.getenv("HPRC_QC_OUTPUT_DIR", output_dir or DEFAULT_OUTPUT_DIR))
    cat_cache_dir = Path(
        os.getenv("HPRC_QC_CAT_CACHE_DIR", cat_cache_dir or output_dir.parent / "cache" / "cat")
    )
    ensembl_cache_dir = Path(
        os.getenv(
            "HPRC_QC_ENSEMBL_CACHE_DIR",
            ensembl_cache_dir or output_dir.parent / "cache" / "ensembl",
        )
    )

    intermediate_dir = output_dir / "intermediate_spreadsheets"
    analysis_dir = intermediate_dir / "exclusive_gene_analysis"
    cache_dir = analysis_dir / "cache"
    figure_dir = analysis_dir / "figures"

    for directory in (analysis_dir, cache_dir, figure_dir):
        directory.mkdir(parents=True, exist_ok=True)

    return {
        "output_dir": output_dir,
        "qc_dir": output_dir / "qc_metrics",
        # gene_pairs_all files are published under results/, not qc_metrics/
        "results_dir": output_dir / "results",
        "intermediate_dir": intermediate_dir,
        "analysis_dir": analysis_dir,
        "cache_dir": cache_dir,
        "figure_dir": figure_dir,
        "cat_cache_dir": cat_cache_dir,
        "ensembl_cache_dir": ensembl_cache_dir,
    }


def discover_files(paths: dict[str, Path]) -> dict[str, list[Path]]:
    qc_dir = paths["qc_dir"]
    results_dir = paths["results_dir"]
    gene_presence_files = sorted(qc_dir.glob("*/*_gene_presence.tsv"))
    ensembl_gene_count_files = sorted(
        path
        for path in qc_dir.glob("*/*_gene_transcript_counts.tsv")
        if not path.name.endswith("_cat_gene_transcript_counts.tsv")
    )
    cat_gene_count_files = sorted(qc_dir.glob("*/*_cat_gene_transcript_counts.tsv"))
    # gene_pairs_all: published to results/{accession}/{accession}.gene_pairs_all.tsv
    gene_pairs_all_files = sorted(results_dir.glob("*/*.gene_pairs_all.tsv"))
    # multi_mapping: published to qc_metrics/{accession}/{accession}_multi_mapping.tsv
    multi_mapping_files = sorted(qc_dir.glob("*/*_multi_mapping.tsv"))
    return {
        "gene_presence_files": gene_presence_files,
        "ensembl_gene_count_files": ensembl_gene_count_files,
        "cat_gene_count_files": cat_gene_count_files,
        "gene_pairs_all_files": gene_pairs_all_files,
        "multi_mapping_files": multi_mapping_files,
    }


def load_gene_presence(
    paths: dict[str, Path], files: dict[str, list[Path]], force_reload: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame]:
    assembly_cache = paths["cache_dir"] / "assembly_presence_summary.parquet"
    exclusive_cache = paths["cache_dir"] / "exclusive_observations.parquet"

    if assembly_cache.exists() and exclusive_cache.exists() and not force_reload:
        return pd.read_parquet(assembly_cache), pd.read_parquet(exclusive_cache)

    # Per-file chunk caching: if a run crashes mid-way, the next run resumes from the last saved chunk
    chunk_dir = paths["cache_dir"] / "presence_chunks"
    chunk_dir.mkdir(exist_ok=True)

    assembly_rows: list[dict[str, object]] = []

    for i, path in enumerate(files["gene_presence_files"]):
        excl_chunk = chunk_dir / f"excl_{i:05d}.parquet"
        asm_chunk = chunk_dir / f"asm_{i:05d}.parquet"

        if excl_chunk.exists() and asm_chunk.exists() and not force_reload:
            asm_row = pd.read_parquet(asm_chunk)
            if not asm_row.empty:
                assembly_rows.append(asm_row.iloc[0].to_dict())
            continue

        df = pd.read_csv(path, sep="\t", dtype=str)
        if df.empty:
            pd.DataFrame().to_parquet(excl_chunk)
            pd.DataFrame().to_parquet(asm_chunk)
            continue

        df["present_in_ensembl"] = df["present_in_ensembl"].eq("True")
        df["present_in_cat"] = df["present_in_cat"].eq("True")
        df["name_class"] = df["gene_name"].map(classify_gene_name)
        df["is_named"] = df["name_class"] != "ENSG_fallback"
        df["catalog_like"] = df["name_class"].isin(CATALOG_NAME_CLASSES)

        both = df["present_in_ensembl"] & df["present_in_cat"]
        ens_only = df["present_in_ensembl"] & ~df["present_in_cat"]
        cat_only = ~df["present_in_ensembl"] & df["present_in_cat"]

        row: dict[str, object] = {
            "assembly_accession": df["assembly_accession"].iloc[0],
            "sample_name": df["sample_name"].iloc[0],
            "n_total": int(len(df)),
            "n_both": int(both.sum()),
            "n_ensembl_only": int(ens_only.sum()),
            "n_cat_only": int(cat_only.sum()),
            "n_named_total": int(df["is_named"].sum()),
            "n_named_both": int((both & df["is_named"]).sum()),
            "n_named_ensembl_only": int((ens_only & df["is_named"]).sum()),
            "n_named_cat_only": int((cat_only & df["is_named"]).sum()),
            "n_catalog_like_total": int(df["catalog_like"].sum()),
            "n_catalog_like_both": int((both & df["catalog_like"]).sum()),
            "n_catalog_like_ensembl_only": int((ens_only & df["catalog_like"]).sum()),
            "n_catalog_like_cat_only": int((cat_only & df["catalog_like"]).sum()),
        }
        pd.DataFrame([row]).to_parquet(asm_chunk, index=False)
        assembly_rows.append(row)

        exclusive = df.loc[ens_only | cat_only, [
            "assembly_accession",
            "sample_name",
            "gene_name",
            "name_class",
            "catalog_like",
            "ensembl_gene_id",
            "cat_gene_id",
            "ensembl_biotype",
            "cat_biotype",
            "ensembl_chrom",
            "cat_chrom",
        ]].copy()
        exclusive["exclusive_side"] = np.where(cat_only.loc[exclusive.index], "cat_only", "ensembl_only")
        exclusive.to_parquet(excl_chunk, index=False)

    assembly_summary = pd.DataFrame(assembly_rows).sort_values("assembly_accession").reset_index(drop=True)
    excl_chunks = [pd.read_parquet(p) for p in sorted(chunk_dir.glob("excl_*.parquet"))]
    excl_chunks = [c for c in excl_chunks if not c.empty]
    exclusive_obs = pd.concat(excl_chunks, ignore_index=True) if excl_chunks else pd.DataFrame()

    assembly_summary.to_parquet(assembly_cache, index=False)
    exclusive_obs.to_parquet(exclusive_cache, index=False)
    assembly_summary.to_csv(paths["analysis_dir"] / "assembly_presence_summary.tsv", sep="\t", index=False)
    exclusive_obs.to_csv(paths["analysis_dir"] / "exclusive_observations.tsv.gz", sep="\t", index=False)
    return assembly_summary, exclusive_obs


def _summarise_gene_count_file(path: Path, annotation: str) -> dict[str, object]:
    df = pd.read_csv(path, sep="\t", dtype={"gene_name": str, "biotype": str})
    df["name_class"] = df["gene_name"].map(classify_gene_name)
    df["non_fallback"] = df["name_class"] != "ENSG_fallback"
    df["catalog_like"] = df["name_class"].isin(CATALOG_NAME_CLASSES)
    protein_coding = df["biotype"].eq("protein_coding")

    return {
        "assembly_accession": df["assembly_accession"].iloc[0],
        f"n_genes_{annotation}_total": int(len(df)),
        f"n_genes_{annotation}_non_fallback": int(df["non_fallback"].sum()),
        f"n_genes_{annotation}_catalog_like": int(df["catalog_like"].sum()),
        f"n_genes_{annotation}_protein_coding": int(protein_coding.sum()),
        f"n_genes_{annotation}_protein_coding_catalog_like": int((protein_coding & df["catalog_like"]).sum()),
    }


def load_gene_counts(
    paths: dict[str, Path], files: dict[str, list[Path]], force_reload: bool = False
) -> pd.DataFrame:
    cache_path = paths["cache_dir"] / "assembly_gene_count_summary.parquet"
    if cache_path.exists() and not force_reload:
        return pd.read_parquet(cache_path)

    ens_rows = [_summarise_gene_count_file(path, "ens") for path in files["ensembl_gene_count_files"]]
    cat_rows = [_summarise_gene_count_file(path, "cat") for path in files["cat_gene_count_files"]]

    ens_df = pd.DataFrame(ens_rows)
    cat_df = pd.DataFrame(cat_rows)
    merged = ens_df.merge(cat_df, on="assembly_accession", how="outer").sort_values("assembly_accession")

    merged["delta_total"] = merged["n_genes_cat_total"] - merged["n_genes_ens_total"]
    merged["delta_non_fallback"] = merged["n_genes_cat_non_fallback"] - merged["n_genes_ens_non_fallback"]
    merged["delta_catalog_like"] = merged["n_genes_cat_catalog_like"] - merged["n_genes_ens_catalog_like"]
    merged["delta_protein_coding"] = (
        merged["n_genes_cat_protein_coding"] - merged["n_genes_ens_protein_coding"]
    )
    merged["delta_protein_coding_catalog_like"] = (
        merged["n_genes_cat_protein_coding_catalog_like"]
        - merged["n_genes_ens_protein_coding_catalog_like"]
    )

    merged.to_parquet(cache_path, index=False)
    merged.to_csv(paths["analysis_dir"] / "assembly_gene_count_summary.tsv", sep="\t", index=False)
    return merged


def build_exclusive_rollups(
    paths: dict[str, Path], exclusive_obs: pd.DataFrame, force_reload: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rollup_cache = paths["cache_dir"] / "exclusive_gene_rollup.parquet"
    name_cache = paths["cache_dir"] / "exclusive_name_class_summary.parquet"
    biotype_cache = paths["cache_dir"] / "exclusive_biotype_summary.parquet"

    if rollup_cache.exists() and name_cache.exists() and biotype_cache.exists() and not force_reload:
        return (
            pd.read_parquet(rollup_cache),
            pd.read_parquet(name_cache),
            pd.read_parquet(biotype_cache),
        )

    rollup = (
        exclusive_obs.groupby("gene_name", dropna=False)
        .agg(
            name_class=("name_class", "first"),
            catalog_like=("catalog_like", "max"),
            n_assemblies_cat_only=("exclusive_side", lambda s: int((s == "cat_only").sum())),
            n_assemblies_ensembl_only=("exclusive_side", lambda s: int((s == "ensembl_only").sum())),
            n_assemblies_exclusive=("assembly_accession", "nunique"),
            ensembl_biotype=("ensembl_biotype", first_nonempty),
            cat_biotype=("cat_biotype", first_nonempty),
            ensembl_gene_id=("ensembl_gene_id", first_nonempty),
            cat_gene_id=("cat_gene_id", first_nonempty),
            ensembl_chrom=("ensembl_chrom", first_nonempty),
            cat_chrom=("cat_chrom", first_nonempty),
        )
        .reset_index()
    )
    rollup["dominant_side"] = np.where(
        rollup["n_assemblies_cat_only"] > rollup["n_assemblies_ensembl_only"],
        "cat_only",
        np.where(
            rollup["n_assemblies_cat_only"] < rollup["n_assemblies_ensembl_only"],
            "ensembl_only",
            "balanced",
        ),
    )

    name_summary = (
        exclusive_obs.groupby(["exclusive_side", "name_class"], dropna=False)
        .size()
        .reset_index(name="n_observations")
        .sort_values(["exclusive_side", "n_observations"], ascending=[True, False])
    )

    biotype_frame = exclusive_obs.copy()
    biotype_frame["exclusive_biotype"] = np.where(
        biotype_frame["exclusive_side"].eq("cat_only"),
        biotype_frame["cat_biotype"],
        biotype_frame["ensembl_biotype"],
    )
    biotype_summary = (
        biotype_frame.groupby(["exclusive_side", "exclusive_biotype"], dropna=False)
        .size()
        .reset_index(name="n_observations")
        .sort_values(["exclusive_side", "n_observations"], ascending=[True, False])
    )

    rollup.to_parquet(rollup_cache, index=False)
    name_summary.to_parquet(name_cache, index=False)
    biotype_summary.to_parquet(biotype_cache, index=False)

    rollup.to_csv(paths["analysis_dir"] / "exclusive_gene_rollup.tsv", sep="\t", index=False)
    name_summary.to_csv(paths["analysis_dir"] / "exclusive_name_class_summary.tsv", sep="\t", index=False)
    biotype_summary.to_csv(paths["analysis_dir"] / "exclusive_biotype_summary.tsv", sep="\t", index=False)
    return rollup, name_summary, biotype_summary


_KEEP_PAIRS_COLS = [
    "assembly_accession",
    "ensembl_gene_id",
    "cat_gene_id",
    "ensembl_name",
    "cat_name",
    "ensembl_biotype",
    "cat_biotype",
    "frac_ensembl_covered",
    "frac_cat_covered",
    "is_rbh",
    "is_spatial",
]


def _normalise_pairs_columns(df: pd.DataFrame, assembly_accession: str) -> pd.DataFrame:
    """Normalise gene_pairs_all column names and add assembly_accession if missing.

    The overlap script uses 'ensembl_id'/'cat_id' in its compute path but
    'ensembl_gene_id'/'cat_gene_id' in the empty-file fallback.  Both produce
    identical string values; we standardise to ensembl_gene_id/cat_gene_id to
    match the rest of the QC system.
    """
    renames: dict[str, str] = {}
    if "ensembl_id" in df.columns and "ensembl_gene_id" not in df.columns:
        renames["ensembl_id"] = "ensembl_gene_id"
    if "cat_id" in df.columns and "cat_gene_id" not in df.columns:
        renames["cat_id"] = "cat_gene_id"
    if renames:
        df = df.rename(columns=renames)
    if "assembly_accession" not in df.columns:
        df = df.assign(assembly_accession=assembly_accession)
    # is_spatial may be absent in name-match-only output versions
    if "is_spatial" not in df.columns:
        df = df.assign(is_spatial=True)
    # Retain only the columns we need (ignore extra columns silently)
    present = [c for c in _KEEP_PAIRS_COLS if c in df.columns]
    return df[present].copy()


def _stream_gene_pairs_all(
    paths: dict[str, Path],
    files: dict[str, list[Path]],
    target_cat_ids: set | None = None,
    target_ens_ids: set | None = None,
    cache_name: str = "pairs_query",
    force_reload: bool = False,
) -> pd.DataFrame:
    """Stream gene_pairs_all files one at a time, keeping only rows matching target IDs.

    Unlike a full load, this never materialises the entire dataset in memory.
    Rows are kept when cat_gene_id is in target_cat_ids OR ensembl_gene_id is in
    target_ens_ids.  If neither set is supplied all rows are kept (avoid on large
    cohorts — tens of millions of rows across 462 assemblies).

    Results are cached to ``{cache_name}.parquet`` so a second call with the same
    ``cache_name`` and ``force_reload=False`` returns instantly.
    """
    cache_path = paths["cache_dir"] / f"{cache_name}.parquet"
    if cache_path.exists() and not force_reload:
        return pd.read_parquet(cache_path)

    chunks: list[pd.DataFrame] = []
    for path in files.get("gene_pairs_all_files", []):
        accession = path.name.replace(".gene_pairs_all.tsv", "")
        try:
            df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        except Exception:
            continue
        if df.empty:
            continue

        df = _normalise_pairs_columns(df, accession)

        # Filter to target IDs
        if target_cat_ids is not None or target_ens_ids is not None:
            mask = pd.Series(False, index=df.index)
            if target_cat_ids is not None and "cat_gene_id" in df.columns:
                mask |= df["cat_gene_id"].isin(target_cat_ids)
            if target_ens_ids is not None and "ensembl_gene_id" in df.columns:
                mask |= df["ensembl_gene_id"].isin(target_ens_ids)
            df = df[mask]

        if df.empty:
            continue

        for col in ("frac_ensembl_covered", "frac_cat_covered"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        for col in ("is_rbh", "is_spatial"):
            if col in df.columns:
                df[col] = df[col].map({"True": True, "False": False}).fillna(False)

        chunks.append(df)

    result = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=_KEEP_PAIRS_COLS)
    result.to_parquet(cache_path, index=False)
    return result


_MM_COLS = [
    "assembly_accession", "sample_name", "gene_id", "source", "n_matches",
    "relationship", "classification", "matched_gene_ids", "matched_biotypes",
    "overlap_fractions", "total_coverage_fraction", "has_rbh_match",
]


def _stream_multi_mapping(
    paths: dict[str, Path],
    files: dict[str, list[Path]],
    target_gene_ids: set | None = None,
    cache_name: str = "mm_query",
    force_reload: bool = False,
) -> pd.DataFrame:
    """Stream multi_mapping files one at a time, keeping only rows matching target gene_ids.

    If target_gene_ids is None all rows are kept (avoid on large cohorts).
    Results are cached to ``{cache_name}.parquet``.
    """
    cache_path = paths["cache_dir"] / f"{cache_name}.parquet"
    if cache_path.exists() and not force_reload:
        return pd.read_parquet(cache_path)

    chunks: list[pd.DataFrame] = []
    for path in files.get("multi_mapping_files", []):
        try:
            df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        except Exception:
            continue
        if df.empty:
            continue

        if target_gene_ids is not None and "gene_id" in df.columns:
            df = df[df["gene_id"].isin(target_gene_ids)]
        if df.empty:
            continue

        for col in ("n_matches", "total_coverage_fraction"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        if "has_rbh_match" in df.columns:
            df["has_rbh_match"] = df["has_rbh_match"].map({"True": True, "False": False}).fillna(False)

        chunks.append(df)

    result = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=_MM_COLS)
    result.to_parquet(cache_path, index=False)
    return result


def lookup_ensembl_grch38_locations(gene_names: list[str]) -> pd.DataFrame:
    """Query the Ensembl REST API for GRCh38 genomic locations of gene symbols.

    Returns a DataFrame with columns:
        gene_name, grch38_seq_region, grch38_start, grch38_end, grch38_strand,
        grch38_assembly_name, is_primary_assembly.

    seq_region will be a chromosome number ('2', 'X') for primary-assembly genes
    or an alt-locus/patch name (e.g. 'CHR_HSCHR2_1_CTG1') for non-primary loci.
    is_primary_assembly is True when seq_region is a bare chromosome number.
    """
    import json
    import urllib.request

    rows: list[dict[str, object]] = []
    batch_size = 200  # Ensembl REST API limit per POST
    for batch_start in range(0, len(gene_names), batch_size):
        batch = gene_names[batch_start : batch_start + batch_size]
        url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
        body = json.dumps({"symbols": batch}).encode()
        req = urllib.request.Request(
            url,
            data=body,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                data: dict = json.loads(resp.read())
        except Exception as exc:
            print(f"Warning: Ensembl REST lookup failed for batch starting {batch_start}: {exc}")
            continue

        for name, info in data.items():
            if not isinstance(info, dict):
                rows.append({"gene_name": name, "grch38_seq_region": None})
                continue
            seq_region = info.get("seq_region_name", "")
            rows.append(
                {
                    "gene_name": name,
                    "grch38_seq_region": seq_region,
                    "grch38_start": info.get("start"),
                    "grch38_end": info.get("end"),
                    "grch38_strand": info.get("strand"),
                    "grch38_assembly_name": info.get("assembly_name", ""),
                    # Primary assembly chromosomes are bare digits or single letters
                    "is_primary_assembly": bool(
                        seq_region and re.match(r"^(\d{1,2}|X|Y|MT)$", str(seq_region))
                    ),
                }
            )

    return pd.DataFrame(rows)


def run_followup_investigations(
    results: dict[str, object],
    force_reload: bool = False,
) -> dict[str, object]:
    """Cross-reference exclusive genes against gene_pairs_all and multi_mapping data.

    Streams gene_pairs_all and multi_mapping one file at a time, keeping only rows
    relevant to the four investigations.  This avoids loading the full multi-GB
    dataset into memory — the target ID sets are computed from already-loaded
    DataFrames before any file I/O begins, then a single filtered pass is made
    through each file collection.

    Cached query results are written to the cache directory so a second call with
    the same ``force_reload=False`` returns immediately.

    Returns a dict with keys:
        readthrough_in_pairs – pairs_all subset for readthrough CAT gene IDs
        readthrough_multi    – multi_mapping rows for readthrough CAT gene IDs
        ulc_in_pairs         – pairs_all subset for unknown_likely_coding CAT IDs
        ens_only_grch38      – GRCh38 location lookup for top Ensembl-only exclusives
        residual_cat_only    – catalog-like CAT-only exclusives classified by pairs status
        residual_ens_only    – catalog-like Ensembl-only exclusives classified by pairs status
    """
    paths: dict[str, Path] = results["paths"]
    files: dict[str, list[Path]] = results["files"]
    cat_only: pd.DataFrame = results["cat_only_enriched"]
    exclusive_obs: pd.DataFrame = results["exclusive_obs"]

    # ── Compute target ID sets from already-loaded DataFrames (no file I/O) ──

    # 1. Readthrough: hyphenated CAT-source genes
    readthrough_ids: set[str] = set(
        cat_only.loc[
            cat_only["name_class"].eq("hyphenated_symbol") & cat_only["cat_source"].eq("CAT"),
            "cat_gene_id",
        ]
        .dropna()
        .astype(str)
    )

    # 2. unknown_likely_coding
    ulc_ids: set[str] = set(
        cat_only.loc[
            cat_only["cat_biotype"].eq("unknown_likely_coding"),
            "cat_gene_id",
        ]
        .dropna()
        .astype(str)
    )

    # 3. Ensembl-only catalog-like observations
    ens_only_catalog_obs = exclusive_obs[
        exclusive_obs["exclusive_side"].eq("ensembl_only") & exclusive_obs["catalog_like"]
    ]
    ens_only_catalog_ids: set[str] = set(
        ens_only_catalog_obs["ensembl_gene_id"].dropna().astype(str)
    ) if "ensembl_gene_id" in ens_only_catalog_obs.columns else set()

    # 4. Residual catalog-like CAT-only obs (not readthrough, not ULC)
    cat_only_catalog_obs = exclusive_obs[
        exclusive_obs["exclusive_side"].eq("cat_only") & exclusive_obs["catalog_like"]
    ]
    cat_only_catalog_ids: set[str] = (
        set(cat_only_catalog_obs["cat_gene_id"].dropna().astype(str))
        if "cat_gene_id" in cat_only_catalog_obs.columns
        else set()
    ) - readthrough_ids - ulc_ids

    # ── ONE streaming pass through gene_pairs_all ─────────────────────────────
    # Union of all CAT IDs we need across investigations 1, 2, and 4
    all_target_cat = readthrough_ids | ulc_ids | cat_only_catalog_ids
    pairs_subset = _stream_gene_pairs_all(
        paths, files,
        target_cat_ids=all_target_cat if all_target_cat else None,
        target_ens_ids=ens_only_catalog_ids if ens_only_catalog_ids else None,
        cache_name="followup_pairs_subset",
        force_reload=force_reload,
    )

    # ── ONE streaming pass through multi_mapping (readthrough IDs only) ───────
    mm_readthrough = _stream_multi_mapping(
        paths, files,
        target_gene_ids=readthrough_ids if readthrough_ids else None,
        cache_name="followup_mm_readthrough",
        force_reload=force_reload,
    )

    # ── Split pairs_subset per investigation ──────────────────────────────────

    # Investigation 1 — readthrough
    if not pairs_subset.empty and "cat_gene_id" in pairs_subset.columns:
        rt_in_pairs = pairs_subset[pairs_subset["cat_gene_id"].isin(readthrough_ids)].copy()
    else:
        rt_in_pairs = pd.DataFrame()

    rt_multi = (
        mm_readthrough[mm_readthrough["source"].eq("cat")].copy()
        if not mm_readthrough.empty and "source" in mm_readthrough.columns
        else pd.DataFrame()
    )

    # Investigation 2 — unknown_likely_coding
    if not pairs_subset.empty and "cat_gene_id" in pairs_subset.columns:
        ulc_in_pairs = pairs_subset[pairs_subset["cat_gene_id"].isin(ulc_ids)].copy()
    else:
        ulc_in_pairs = pd.DataFrame()

    # ── 3. Ensembl-only exclusive gene GRCh38 lookup ─────────────────────────
    ens_only_counts = (
        ens_only_catalog_obs.groupby("gene_name", dropna=False)
        .agg(n_assemblies=("assembly_accession", "nunique"), ensembl_chrom=("ensembl_chrom", first_nonempty))
        .reset_index()
        .sort_values("n_assemblies", ascending=False)
    )
    top_ens_only = ens_only_counts.head(100)["gene_name"].tolist()
    if top_ens_only:
        grch38_locs = lookup_ensembl_grch38_locations(top_ens_only)
        ens_only_grch38 = ens_only_counts.merge(grch38_locs, on="gene_name", how="left")
    else:
        ens_only_grch38 = ens_only_counts.copy()

    # ── 4. Residual catalog-like exclusives ───────────────────────────────────
    def _classify_against_pairs(
        obs: pd.DataFrame, id_col: str, pairs: pd.DataFrame, pairs_id_col: str
    ) -> pd.DataFrame:
        """Label each exclusive gene as 'in_pairs_all_non_rbh', 'rbh_only', or 'absent'."""
        if obs.empty:
            return obs.assign(pairs_status="absent")
        if pairs.empty or pairs_id_col not in pairs.columns:
            return obs.assign(pairs_status="pairs_data_unavailable")
        ids_non_rbh = set(
            pairs.loc[~pairs["is_rbh"].eq(True), pairs_id_col].dropna().astype(str)
        ) if "is_rbh" in pairs.columns else set()
        ids_rbh_only = set(
            pairs.loc[pairs["is_rbh"].eq(True), pairs_id_col].dropna().astype(str)
        ) if "is_rbh" in pairs.columns else set()
        ids_in_pairs = set(pairs[pairs_id_col].dropna().astype(str))

        def classify(gid: object) -> str:
            gid = str(gid) if pd.notna(gid) else ""
            if gid in ids_non_rbh:
                return "in_pairs_all_non_rbh"
            if gid in ids_rbh_only:
                return "rbh_only"
            if gid in ids_in_pairs:
                return "in_pairs_all"
            return "absent"

        return obs.assign(pairs_status=obs[id_col].map(classify))

    residual_cat = cat_only_catalog_obs[
        cat_only_catalog_obs.get("cat_gene_id", pd.Series(dtype=str)).isin(cat_only_catalog_ids)
    ] if "cat_gene_id" in cat_only_catalog_obs.columns else cat_only_catalog_obs

    pairs_for_cat = (
        pairs_subset[pairs_subset["cat_gene_id"].isin(cat_only_catalog_ids)]
        if not pairs_subset.empty and "cat_gene_id" in pairs_subset.columns
        else pd.DataFrame()
    )
    pairs_for_ens = (
        pairs_subset[pairs_subset["ensembl_gene_id"].isin(ens_only_catalog_ids)]
        if not pairs_subset.empty and "ensembl_gene_id" in pairs_subset.columns
        else pd.DataFrame()
    )

    residual_cat_classified = _classify_against_pairs(
        residual_cat, "cat_gene_id", pairs_for_cat, "cat_gene_id"
    )
    residual_ens_classified = _classify_against_pairs(
        ens_only_catalog_obs, "ensembl_gene_id", pairs_for_ens, "ensembl_gene_id"
    )

    return {
        "readthrough_in_pairs": rt_in_pairs,
        "readthrough_multi": rt_multi,
        "ulc_in_pairs": ulc_in_pairs,
        "ens_only_grch38": ens_only_grch38,
        "residual_cat_only": residual_cat_classified,
        "residual_ens_only": residual_ens_classified,
    }


def parse_cat_only_metadata(
    paths: dict[str, Path], exclusive_obs: pd.DataFrame, force_reload: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cat_meta_cache = paths["cache_dir"] / "cat_only_gene_metadata.parquet"
    cat_only_cache = paths["cache_dir"] / "cat_only_enriched.parquet"
    source_summary_cache = paths["cache_dir"] / "cat_only_source_summary.parquet"

    if (
        cat_meta_cache.exists()
        and cat_only_cache.exists()
        and source_summary_cache.exists()
        and not force_reload
    ):
        return (
            pd.read_parquet(cat_meta_cache),
            pd.read_parquet(cat_only_cache),
            pd.read_parquet(source_summary_cache),
        )

    cat_only = exclusive_obs[exclusive_obs["exclusive_side"].eq("cat_only")].copy()
    targets = cat_only[["sample_name", "cat_gene_id", "gene_name"]].drop_duplicates()

    # Per-sample caching: saves each sample's extracted rows so a crash can resume mid-loop
    sample_cache_dir = paths["cache_dir"] / "cat_gff_samples"
    sample_cache_dir.mkdir(exist_ok=True)

    rows: list[dict[str, object]] = []
    for sample_name, target in targets.groupby("sample_name", dropna=False):
        sample_cache = sample_cache_dir / f"{sample_name}.parquet"
        if sample_cache.exists() and not force_reload:
            cached = pd.read_parquet(sample_cache)
            if not cached.empty:
                rows.extend(cached.to_dict("records"))
            continue

        gff_path = paths["cat_cache_dir"] / f"{sample_name}_cat.gff3.gz"
        if not gff_path.exists():
            pd.DataFrame().to_parquet(sample_cache)
            continue

        wanted_gene_ids = set(target["cat_gene_id"].dropna().astype(str))
        wanted_gene_names = set(target["gene_name"].dropna().astype(str))

        # _parse_gff3_gene_lines uses zcat+awk to return only gene-feature lines,
        # which is orders of magnitude faster than Python-level gzip iteration
        sample_rows: list[dict[str, object]] = []
        for line in _parse_gff3_gene_lines(gff_path):
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            attrs = parse_gff_attributes(parts[8])
            gene_id = attrs.get("ID", attrs.get("gene_id", ""))
            gene_name = attrs.get("gene_name", attrs.get("Name", attrs.get("gene", gene_id)))

            if gene_id not in wanted_gene_ids and gene_name not in wanted_gene_names:
                continue

            sample_rows.append(
                {
                    "sample_name": sample_name,
                    "cat_gene_id": gene_id,
                    "gene_name": gene_name,
                    "cat_source": parts[1],
                    "cat_source_biotype": attrs.get("gene_biotype", attrs.get("biotype", "unknown")),
                    "cat_source_chrom": parts[0],
                }
            )

        sample_df = pd.DataFrame(sample_rows) if sample_rows else pd.DataFrame()
        sample_df.to_parquet(sample_cache)
        rows.extend(sample_rows)

    cat_meta = pd.DataFrame(rows).drop_duplicates(["sample_name", "cat_gene_id", "gene_name"])

    cat_only = cat_only.merge(
        cat_meta[["sample_name", "cat_gene_id", "cat_source", "cat_source_biotype", "cat_source_chrom"]],
        on=["sample_name", "cat_gene_id"],
        how="left",
    )
    cat_only["cat_source"] = cat_only["cat_source"].fillna("UNRESOLVED")
    cat_only["cat_source_biotype"] = cat_only["cat_source_biotype"].fillna(cat_only["cat_biotype"])
    cat_only["cat_source_chrom"] = cat_only["cat_source_chrom"].fillna(cat_only["cat_chrom"])

    source_summary = (
        cat_only.groupby(["cat_source", "name_class", "cat_biotype"], dropna=False)
        .size()
        .reset_index(name="n_observations")
        .sort_values("n_observations", ascending=False)
    )

    cat_meta.to_parquet(cat_meta_cache, index=False)
    cat_only.to_parquet(cat_only_cache, index=False)
    source_summary.to_parquet(source_summary_cache, index=False)

    cat_only.to_csv(paths["analysis_dir"] / "cat_only_enriched.tsv.gz", sep="\t", index=False)
    source_summary.to_csv(paths["analysis_dir"] / "cat_only_source_summary.tsv", sep="\t", index=False)
    return cat_meta, cat_only, source_summary


def run_analysis(
    output_dir: str | Path | None = None,
    cat_cache_dir: str | Path | None = None,
    ensembl_cache_dir: str | Path | None = None,
    force_reload: bool | None = None,
) -> dict[str, object]:
    if force_reload is None:
        force_reload = os.getenv("HPRC_QC_FORCE_RELOAD", "0") == "1"

    paths = resolve_paths(output_dir, cat_cache_dir, ensembl_cache_dir)
    files = discover_files(paths)
    assembly_presence, exclusive_obs = load_gene_presence(paths, files, force_reload=force_reload)
    assembly_gene_counts = load_gene_counts(paths, files, force_reload=force_reload)
    assembly_overview = assembly_presence.merge(assembly_gene_counts, on="assembly_accession", how="left")
    assembly_overview.to_csv(paths["analysis_dir"] / "assembly_overview.tsv", sep="\t", index=False)

    exclusive_rollup, name_summary, biotype_summary = build_exclusive_rollups(
        paths, exclusive_obs, force_reload=force_reload
    )
    cat_meta, cat_only_enriched, source_summary = parse_cat_only_metadata(
        paths, exclusive_obs, force_reload=force_reload
    )

    return {
        "paths": paths,
        "files": files,
        "assembly_presence": assembly_presence,
        "assembly_gene_counts": assembly_gene_counts,
        "assembly_overview": assembly_overview,
        "exclusive_obs": exclusive_obs,
        "exclusive_rollup": exclusive_rollup,
        "exclusive_name_summary": name_summary,
        "exclusive_biotype_summary": biotype_summary,
        "cat_only_metadata": cat_meta,
        "cat_only_enriched": cat_only_enriched,
        "cat_only_source_summary": source_summary,
    }


def make_overview_figure(results: dict[str, object], top_n: int = 20) -> tuple[plt.Figure, np.ndarray]:
    assembly = results["assembly_overview"].copy()
    cat_only = results["cat_only_enriched"].copy()
    exclusive_rollup = results["exclusive_rollup"].copy()
    figure_dir = results["paths"]["figure_dir"]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    ax = axes[0]
    x = assembly["n_genes_ens_protein_coding_catalog_like"]
    y = assembly["n_genes_cat_protein_coding_catalog_like"]
    ax.scatter(x, y, s=18, alpha=0.7, color="#1f77b4")
    low = min(x.min(), y.min())
    high = max(x.max(), y.max())
    ax.plot([low, high], [low, high], color="black", lw=1, ls="--")
    ax.set_xlabel("Ensembl protein-coding catalog-like genes")
    ax.set_ylabel("CAT protein-coding catalog-like genes")
    ax.set_title("Per-assembly canonical protein-coding counts")

    ax = axes[1]
    top_cat = (
        cat_only[cat_only["cat_biotype"].eq("protein_coding") & cat_only["catalog_like"]]
        .groupby(["gene_name", "cat_source"], dropna=False)
        .size()
        .reset_index(name="n_assemblies")
        .sort_values("n_assemblies", ascending=False)
        .head(top_n)
        .sort_values("n_assemblies")
    )
    colors = top_cat["cat_source"].map({"CAT": "#d62728", "Liftoff": "#2ca02c"}).fillna("#7f7f7f")
    ax.barh(range(len(top_cat)), top_cat["n_assemblies"], color=colors)
    ax.set_yticks(range(len(top_cat)))
    ax.set_yticklabels(top_cat["gene_name"], fontsize=8)
    ax.set_xlabel("CAT-only assemblies")
    ax.set_title("Top CAT-only catalog-like protein-coding genes")

    fig.tight_layout()
    fig.savefig(figure_dir / "exclusive_gene_overview.png", dpi=200, bbox_inches="tight")
    return fig, axes
