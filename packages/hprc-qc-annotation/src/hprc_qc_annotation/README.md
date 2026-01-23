# hprc_ensembl_cat_overlap.py

Quantify overlap between **Ensembl HPRC Release 2** gene annotations and **CAT** gene annotations across HPRC assemblies.

The script compares gene models using coordinate-based overlap, derives reciprocal best-hit gene pairs, and summarizes concordance at the gene, transcript, biotype, and assembly levels.

---

## Overview

For each HPRC Release 2 assembly where both annotations are available, the script produces:

- Per-assembly summary metrics
- Per gene-pair overlap records (Ensembl ↔ CAT)
- Reciprocal Best Hit (RBH) gene pairs
- Transcript concordance within RBH gene pairs (intron-chain based)
- Biotype-stratified overlap statistics

The core philosophy is **coordinates first, names second**:
genomic overlap on the same contig and strand is the primary signal; gene names are used only as supporting evidence.

---

## Dependencies

### Python

- Python 3.x (3.9+ recommended)

### Python packages

- `pandas`
- `pyranges`
- `numpy` (indirect dependency)

Install dependencies with:

```bash
pip install pandas pyranges numpy
```

---

## Inputs

### Required index files (bulk mode)

#### Assemblies index CSV

Example: `assemblies_release2_v1.0.index.csv`

Must contain:
- `Assembly Accession`
- `Sample Name`

May optionally include a column pointing directly to the Ensembl GFF filename.

---

#### CAT index CSV

Example: `cat_genes_hprc_r2_v1.2.index.csv`

Must map each `Sample Name` to a CAT GFF filename.

---

### Annotation directories

- `--ensembl-gff-dir`  
  Directory containing Ensembl GFFs (typically named by assembly accession)

- `--cat-anno-dir`  
  Directory containing CAT GFFs (typically named by sample)

---

## Usage

### Bulk mode (many assemblies)

```bash
python3 hprc_ensembl_cat_overlap.py \
  --ensembl-gff-dir /path/to/ensembl_gffs \
  --cat-index /path/to/cat_index.csv \
  --assemblies-index /path/to/assemblies_index.csv \
  --cat-anno-dir /path/to/cat_annos \
  --output-prefix ./results/overlap \
  --contig-normalization basic
```

Restrict to a subset (for parallelization):

```bash
--assemblies GCA_018466835.2,HG002
```

---

### Single-pair / region test mode

```bash
python3 hprc_ensembl_cat_overlap.py \
  --ensembl-gff ensembl.gff3 \
  --cat-gff cat.gff3 \
  --chrom 2 \
  --output-prefix ./results/chr2_test \
  --contig-normalization basic
```

`--chrom` filters on the GFF seqid and works for any contig or scaffold name, not just chromosomes.

---

## Contig normalization

- `none`: no changes
- `basic`: strip `chr`, convert `MT` → `M`
- `cat_hash`: keep only the final field after `#` in CAT contig names

Assembly-report-based contig mapping is supported only in single-pair mode.

---

## Outputs

- `PREFIX.assembly_summary.tsv`
- `PREFIX.gene_pairs_all.tsv`
- `PREFIX.gene_pairs_rbh.tsv`
- `PREFIX.transcript_concordance.tsv`
- `PREFIX.biotype_stats.tsv`

---

## Parallelization model

Parallelization is external. Use `--assemblies` to split work across jobs (Snakemake / Slurm array jobs).

---

## Outstanding implementation gaps

- Assembly report mapping not applied in bulk mode
- Gene-level comparison is span-based
- Multi-overlap interpretation is simplistic
- Gene name matching is exact-string only
- Transcript concordance is structure-only (no CDS awareness)
- No per-contig or regional summaries

---

## Diagnostic workflow example

```bash
python3 hprc_ensembl_cat_overlap.py \
  --ensembl-gff ensembl.gff3 \
  --cat-gff cat.gff3 \
  --chrom CM000663.2 \
  --output-prefix ./results/diag_region \
  --contig-normalization cat_hash
```
