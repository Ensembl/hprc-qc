from __future__ import annotations

import gzip
import os
import re
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
        "intermediate_dir": intermediate_dir,
        "analysis_dir": analysis_dir,
        "cache_dir": cache_dir,
        "figure_dir": figure_dir,
        "cat_cache_dir": cat_cache_dir,
        "ensembl_cache_dir": ensembl_cache_dir,
    }


def discover_files(paths: dict[str, Path]) -> dict[str, list[Path]]:
    qc_dir = paths["qc_dir"]
    gene_presence_files = sorted(qc_dir.glob("*/*_gene_presence.tsv"))
    ensembl_gene_count_files = sorted(
        path
        for path in qc_dir.glob("*/*_gene_transcript_counts.tsv")
        if not path.name.endswith("_cat_gene_transcript_counts.tsv")
    )
    cat_gene_count_files = sorted(qc_dir.glob("*/*_cat_gene_transcript_counts.tsv"))
    return {
        "gene_presence_files": gene_presence_files,
        "ensembl_gene_count_files": ensembl_gene_count_files,
        "cat_gene_count_files": cat_gene_count_files,
    }


def load_gene_presence(
    paths: dict[str, Path], files: dict[str, list[Path]], force_reload: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame]:
    assembly_cache = paths["cache_dir"] / "assembly_presence_summary.parquet"
    exclusive_cache = paths["cache_dir"] / "exclusive_observations.parquet"

    if assembly_cache.exists() and exclusive_cache.exists() and not force_reload:
        return pd.read_parquet(assembly_cache), pd.read_parquet(exclusive_cache)

    assembly_rows: list[dict[str, object]] = []
    exclusive_chunks: list[pd.DataFrame] = []

    for path in files["gene_presence_files"]:
        df = pd.read_csv(path, sep="\t", dtype=str)
        if df.empty:
            continue

        df["present_in_ensembl"] = df["present_in_ensembl"].eq("True")
        df["present_in_cat"] = df["present_in_cat"].eq("True")
        df["name_class"] = df["gene_name"].map(classify_gene_name)
        df["is_named"] = df["name_class"] != "ENSG_fallback"
        df["catalog_like"] = df["name_class"].isin(CATALOG_NAME_CLASSES)

        both = df["present_in_ensembl"] & df["present_in_cat"]
        ens_only = df["present_in_ensembl"] & ~df["present_in_cat"]
        cat_only = ~df["present_in_ensembl"] & df["present_in_cat"]

        assembly_rows.append(
            {
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
        )

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
        exclusive_chunks.append(exclusive)

    assembly_summary = pd.DataFrame(assembly_rows).sort_values("assembly_accession").reset_index(drop=True)
    exclusive_obs = pd.concat(exclusive_chunks, ignore_index=True) if exclusive_chunks else pd.DataFrame()

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

    rows: list[dict[str, object]] = []
    for sample_name, target in targets.groupby("sample_name", dropna=False):
        gff_path = paths["cat_cache_dir"] / f"{sample_name}_cat.gff3.gz"
        if not gff_path.exists():
            continue

        wanted_gene_ids = set(target["cat_gene_id"].dropna().astype(str))
        wanted_gene_names = set(target["gene_name"].dropna().astype(str))

        with gzip.open(gff_path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue

                attrs = parse_gff_attributes(parts[8])
                gene_id = attrs.get("ID", attrs.get("gene_id", ""))
                gene_name = attrs.get("gene_name", attrs.get("Name", attrs.get("gene", gene_id)))

                if gene_id not in wanted_gene_ids and gene_name not in wanted_gene_names:
                    continue

                rows.append(
                    {
                        "sample_name": sample_name,
                        "cat_gene_id": gene_id,
                        "gene_name": gene_name,
                        "cat_source": parts[1],
                        "cat_source_biotype": attrs.get("gene_biotype", attrs.get("biotype", "unknown")),
                        "cat_source_chrom": parts[0],
                    }
                )

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
