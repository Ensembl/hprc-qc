#!/usr/bin/env python3
"""Execute all poster figure notebooks and collect vector/raster exports."""

from pathlib import Path
import os
import shutil
import subprocess
import sys

ROOT = Path(os.environ.get('HPRC_QC_ROOT', Path(__file__).resolve().parent)).resolve()
OUT = Path(os.environ.get('HPRC_POSTER_EXPORT_DIR', ROOT / 'results/poster_vector_exports')).resolve()
NOTEBOOK_DIR = ROOT / 'nextflow/pipelines/ensembl_cat_comparison/notebooks'
PATCHED_DIR = ROOT / 'poster_exports_v1/notebooks_vector'
EXECUTED_DIR = OUT / 'executed_notebooks'

def run(args):
    print('+', ' '.join(str(x) for x in args), flush=True)
    subprocess.run(args, cwd=ROOT, check=True)

OUT.mkdir(parents=True, exist_ok=True)
EXECUTED_DIR.mkdir(parents=True, exist_ok=True)
run([sys.executable, str(ROOT / 'make_poster_vector_notebooks.py')])

notebooks = [
    'annotation_exploration.ipynb',
    'figure_main1_sankey.ipynb',
    'figure_main_concordance_v2.ipynb',
    'figure_supp_a_projection_rates.ipynb',
    'figure_supp_b_biotype_exclusives.ipynb',
    'figure_supp_c_cds_heatmaps.ipynb',
]

for name in notebooks:
    run([
        'jupyter', 'nbconvert', '--to', 'notebook', '--execute',
        '--ExecutePreprocessor.timeout=-1',
        f'--ExecutePreprocessor.cwd={NOTEBOOK_DIR}',
        '--output-dir', str(EXECUTED_DIR), '--output', name,
        str(PATCHED_DIR / name),
    ])

# The export cells use the notebook working directory so all original relative
# data paths remain valid. Consolidate their products into the requested folder.
generated = NOTEBOOK_DIR / 'poster_vector_exports'
if generated.exists():
    for item in generated.iterdir():
        target = OUT / item.name
        if item.is_file():
            shutil.copy2(item, target)

print(f'Poster exports collected in {OUT}')
