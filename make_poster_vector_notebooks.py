#!/usr/bin/env python3
"""Create reproducible notebook copies with vector/publication exports added."""

from pathlib import Path
import json
import os
import subprocess

ROOT = Path(os.environ.get('HPRC_QC_ROOT', Path(__file__).resolve().parent)).resolve()
OUT = ROOT / 'poster_exports_v1' / 'notebooks_vector'
OUT.mkdir(parents=True, exist_ok=True)

TARGETS = [
    ('origin/marker2_comparison', 'nextflow/pipelines/ensembl_cat_comparison/notebooks/annotation_exploration.ipynb',
     {9: 'concordance_by_chromosome', 11: 'discordant_gene_density', 16: 'gene_level_cds_agreement',
      27: 'coding_loss_asymmetry', 33: 'ensembl_cds_rescue_candidates', 36: 'joint_exon_count_differences',
      38: 'POLR2A_example', 49: 'pangenome_core_shell_cloud'}),
    ('origin/marker2_comparison', 'nextflow/pipelines/ensembl_cat_comparison/notebooks/figure_main_concordance_v2.ipynb',
     {28: 'transcript_count_ratio'}),
    ('origin/marker2_comparison', 'nextflow/pipelines/ensembl_cat_comparison/notebooks/figure_supp_a_projection_rates.ipynb',
     {4: 'projection_rates'}),
    ('origin/marker2_comparison', 'nextflow/pipelines/ensembl_cat_comparison/notebooks/figure_supp_b_biotype_exclusives.ipynb',
     {3: 'method_exclusive_gene_composition'}),
    ('origin/marker2_comparison', 'nextflow/pipelines/ensembl_cat_comparison/notebooks/figure_supp_c_cds_heatmaps.ipynb',
     {3: 'cds_integrity_heatmaps'}),
    ('origin/feature/ensembl_cat_compare_workflow', 'nextflow/pipelines/ensembl_cat_comparison/notebooks/figure_main1_sankey.ipynb',
     {17: 'concordance_ladder'}),
]

def source_from_git(branch, path):
    return subprocess.check_output(['git', 'show', f'{branch}:{path}'], cwd=ROOT)

def export_block(name):
    return f'''\n\n# --- Poster publication export: {name} ---\n# SVG/PDF preserve vector lines and editable text; PNG is a 600 dpi fallback.\n_poster_fig = globals().get("fig")\nif _poster_fig is not None:\n    _poster_dir = Path("poster_vector_exports")\n    _poster_dir.mkdir(parents=True, exist_ok=True)\n    _poster_fig.savefig(_poster_dir / "{name}.svg", format="svg", bbox_inches="tight", pad_inches=0.05)\n    _poster_fig.savefig(_poster_dir / "{name}.pdf", format="pdf", bbox_inches="tight", pad_inches=0.05)\n    _poster_fig.savefig(_poster_dir / "{name}.png", format="png", dpi=600, bbox_inches="tight", pad_inches=0.05)\n    print("Saved poster exports:", "{name}.svg", "{name}.pdf", "{name}.png")\n'''

for branch, path, cells in TARGETS:
    data = json.loads(source_from_git(branch, path))
    # Keep the source notebook's default relative path, but allow HPC runs to
    # point at the existing large supplementary-data handoff directory.
    if path.endswith('/annotation_exploration.ipynb'):
        for cell in data['cells']:
            if cell.get('cell_type') == 'code':
                source = ''.join(cell.get('source', []))
                if "DATA  = Path('../../../../data/supplementary_tables')" in source:
                    source = source.replace(
                        "DATA  = Path('../../../../data/supplementary_tables')",
                        "DATA  = Path(os.environ.get('HPRC_POSTER_DATA_DIR', '../../../../data/supplementary_tables'))",
                    )
                    source = 'import os\n' + source if 'import os\n' not in source else source
                    cell['source'] = source
                    break
    for idx, name in cells.items():
        if idx >= len(data['cells']) or data['cells'][idx].get('cell_type') != 'code':
            raise RuntimeError(f'Expected code cell {idx} in {path}')
        source = ''.join(data['cells'][idx].get('source', []))
        if 'Poster publication export:' not in source:
            data['cells'][idx]['source'] = source + export_block(name)
    out_name = Path(path).name
    data.setdefault('metadata', {}).setdefault('poster_export', {})['source_branch'] = branch
    data['metadata']['poster_export']['source_path'] = path
    data['metadata']['poster_export']['export_format'] = 'SVG + PDF + PNG@600dpi'
    (OUT / out_name).write_text(json.dumps(data, indent=1) + '\n', encoding='utf-8')

readme = OUT / 'README.md'
readme.write_text(
    '# Publication-export notebook copies\n\n'
    'These notebooks are generated from the named Git branches without modifying those branches.\n'
    'Each poster figure cell now saves SVG, PDF and 600-dpi PNG into a `poster_vector_exports/`\n'
    'directory when the notebook is run in an environment with its normal data paths.\n\n'
    'The SVG/PDF files are the preferred Inkscape sources. The PNG files are only fallbacks.\n'
    'The exact source branch and notebook path are stored in each notebook metadata.\n',
    encoding='utf-8'
)
print(f'Created {len(TARGETS)} publication-export notebook copies in {OUT}')
