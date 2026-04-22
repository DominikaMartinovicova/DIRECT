#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_kontextual_fig2.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Reproduce Figure 2A, 2B, and 2D style plots from Ameen et al. 2025
#   (Cell Reports Methods) using local DIRECT data.
#
#   Panel A  — full sample spatial scatter, all cell types coloured using the
#              project palette.
#   Panel B  — L-function vs Kontextual curves over radii with simulation
#              ribbons (mean ± std of permutations).
#   Panel D  — spatial scatter for a chosen cell pair: source (type_i),
#              target (type_j), and remaining parent cells only.
#              Cells outside the parent pool are omitted entirely.
#
#   Inputs:
#       --sample        sample ID matching the per_sample results folder name
#       --cell_pairs    one or more "from__to" strings (space-separated)
#       --adata         path to combined h5ad (for spatial coords + cell types)
#       --input_dir     per-sample kontextual results directory
#       --celltype_col  obs column with cell type labels
#       --o_plots       output directory for plots
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#   python3 plot_kontextual_fig2.py \
#       --sample T23_004535_110005_1 \
#       --cell_pairs Tumor_cells__T_cell_CD8_functional B_cell__T_cell_regulatory \
#       --adata data/combined/Neutro_Epi_extImm_pooled_A_EM_N_combined_adatas_for_analysis_w_v1.7.h5ad \
#       --input_dir results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_sample \
#       --o_plots plots/...


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os
import pickle
import argparse
import anndata as ad


def parse_args():
    parser = argparse.ArgumentParser(
        prog='python3 plot_kontextual_fig2.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='Recreate Kontextual paper Fig 2A/B/D style plots.')
    parser.add_argument('--sample', dest='sample', type=str, required=True,
                        help='Sample ID (folder name in input_dir).')
    parser.add_argument('--cell_pairs', dest='cell_pairs', type=str, nargs='+', required=True,
                        help='One or more cell pairs as "from__to".')
    parser.add_argument('--adata', dest='adata_path', type=str, required=True,
                        help='Path to combined h5ad file.')
    parser.add_argument('--input_dir', dest='input_dir', type=str, required=True,
                        help='Per-sample kontextual results directory.')
    parser.add_argument('--celltype_col', dest='celltype_col', type=str,
                        default='Neutro_Epi_extImm_pooled_A_EM_N',
                        help='obs column with cell type labels.')
    parser.add_argument('--o_plots', dest='output_dir_plots', type=str, required=True,
                        help='Output directory for plots.')
    return parser.parse_args()


args       = parse_args()
sample     = args.sample
cell_pairs = args.cell_pairs
adata_path = args.adata_path
input_dir  = args.input_dir
ct_col     = args.celltype_col
out_dir    = args.output_dir_plots

os.makedirs(out_dir, exist_ok=True)

print(f'Sample: {sample}')
print(f'Cell pairs: {cell_pairs}')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Colour palette
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
celltype_palette = {
    "B_cell": "#1f77b4",
    "DC_mature": "#d62728",
    "Endothelial_cell": "#ff9896",
    "Epithelial_cell": "#c7c7c7",
    "Macrophage": "#ffbb78",
    "Macrophage_alveolar": "#ff7f0e",
    "Mast_cell": "#8b8b12",
    "Monocyte_classical": "#bcbd22",
    "Monocyte_non-classical": "#dbdb8d",
    "NAN": "#9edae5",
    "NK_cell": "#8c564b",
    "Plasma_cell": "#aec7e8",
    "Stromal": "#c49c94",
    "TAN": "#17becf",
    "T_cell_CD4": "#f7b6d2",
    "T_cell_CD8_functional": "#2ca02c",
    "T_cell_CD8_terminally_exhausted": "#98df8a",
    "T_cell_NK-like": "#c5b0d5",
    "T_cell_regulatory": "#e377c2",
    "Tumor_cells": "#7f7f7f",
    "cDC1": "#800020",
    "cDC2": "#CF4265",
    "pDC": "#9467bd",
}

# Article-style colours for Panel D roles (source=dark, target=red, parent=light gray)
COLOR_SOURCE_DEFAULT = '#1a1a1a'   # near-black — source cells
COLOR_TARGET_DEFAULT = '#E63946'   # red        — target cells
COLOR_PARENT         = '#BDBDBD'   # light gray — remaining parent cells


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Load data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2a. Spatial coordinates and cell type labels for this sample
#------------------------------------------------------------------------------
print('Loading adata...')
adata = ad.read_h5ad(adata_path, backed='r')
sample_mask = adata.obs['sample'] == sample
if sample_mask.sum() == 0:
    raise ValueError(f'Sample {sample} not found in adata.')

sample_adata = adata[sample_mask]
coords     = sample_adata.obsm['spatial']
cell_types = sample_adata.obs[ct_col].values
adata.file.close()

df = pd.DataFrame({'x': coords[:, 0], 'y': coords[:, 1], 'cell_type': cell_types})
print(f'  {len(df)} cells loaded. Cell types: {df["cell_type"].nunique()}')

# 2b. Kontextual pkl for this sample
#------------------------------------------------------------------------------
pkl_path = os.path.join(input_dir, sample, 'kontextual_ripleys', 'kontextual_ripleys_L.pkl')
if not os.path.exists(pkl_path):
    raise FileNotFoundError(f'pkl not found: {pkl_path}')

with open(pkl_path, 'rb') as f:
    kont_data = pickle.load(f)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Panel A — full sample spatial scatter, all cell types (Figure 2A equivalent)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plot_full_scatter(df, sample, celltype_palette, out_dir):
    all_celltypes = sorted(df['cell_type'].unique())

    fig, ax = plt.subplots(figsize=(8, 8))

    for ct in all_celltypes:
        mask  = df['cell_type'] == ct
        color = celltype_palette.get(ct, '#aaaaaa')
        ax.scatter(df.loc[mask, 'x'], df.loc[mask, 'y'],
                   c=color, s=3, alpha=0.8, linewidths=0,
                   rasterized=True)

    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_xlabel('x (µm)', fontsize=10)
    ax.set_ylabel('y (µm)', fontsize=10)
    ax.set_title(sample, fontsize=11)
    ax.tick_params(labelsize=8)

    legend_handles = [
        mpatches.Patch(color=celltype_palette.get(ct, '#aaaaaa'), label=ct)
        for ct in all_celltypes
    ]
    ax.legend(handles=legend_handles, fontsize=6, loc='upper right',
              framealpha=0.8, bbox_to_anchor=(1.38, 1.0),
              title='Cell type', title_fontsize=7)

    plt.tight_layout()
    out = os.path.join(out_dir, f'{sample}_fig2A_full_scatter.svg')
    plt.savefig(out, format='svg', bbox_inches='tight', dpi=300)
    plt.close()
    print(f'  Saved: {out}')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Panel B — L-function vs Kontextual curves (Figure 2B equivalent)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plot_curves(pair_data, sample, type_i, type_j, cell_pair, out_dir):
    r           = pair_data['r']
    parent_name = pair_data['parent_name']

    # Centre L values at 0 (subtract r, since L = sqrt(K/pi))
    L_orig_c      = pair_data['L_original']           - r
    L_orig_mean_c = pair_data['L_original_sim_mean']  - r
    L_orig_std    = pair_data['L_original_sim_std']
    L_kont_c      = pair_data['L_kontextual']          - r
    L_kont_mean_c = pair_data['L_kontextual_sim_mean'] - r
    L_kont_std    = pair_data['L_kontextual_sim_std']

    COLOR_LFUNC = '#00B4D8'   # cyan-blue — L-function (matches article)
    COLOR_KONT  = '#E63946'   # red       — Kontextual (matches article)

    fig, ax = plt.subplots(figsize=(6, 4.5))

    ax.plot(r, L_orig_c, color=COLOR_LFUNC, lw=2, label='L-function')
    ax.fill_between(r,
                    L_orig_mean_c - L_orig_std,
                    L_orig_mean_c + L_orig_std,
                    color=COLOR_LFUNC, alpha=0.2)

    ax.plot(r, L_kont_c, color=COLOR_KONT, lw=2, label='Kontextual')
    ax.fill_between(r,
                    L_kont_mean_c - L_kont_std,
                    L_kont_mean_c + L_kont_std,
                    color=COLOR_KONT, alpha=0.2)

    ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.6)

    ax.set_xlabel('Radius r (µm)', fontsize=11)
    ax.set_ylabel('Relationship value (L − r)', fontsize=11)
    ax.set_title(f'{sample}\n{type_i} → {type_j}  |  context: {parent_name}', fontsize=9)
    ax.legend(fontsize=10, framealpha=0.7)

    plt.tight_layout()
    out = os.path.join(out_dir, f'{sample}_{cell_pair}_fig2B_curves.svg')
    plt.savefig(out, format='svg', bbox_inches='tight', dpi=300)
    plt.close()
    print(f'  Saved: {out}')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5 Panel D — highlighted pair scatter, parent context only (Figure 2D equiv.)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plot_pair_scatter(df, pair_data, sample, type_i, type_j, cell_pair,
                      celltype_palette, out_dir):
    parent_name     = pair_data['parent_name']
    perm_pool_types = pair_data['perm_pool_types']   # all types in parent context

    remaining_parent = [ct for ct in perm_pool_types if ct not in (type_i, type_j)]

    def classify_role(ct):
        if ct == type_i:
            return 'source'
        elif ct == type_j:
            return 'target'
        elif ct in remaining_parent:
            return 'parent'
        else:
            return None   # excluded

    df = df.copy()
    df['role'] = df['cell_type'].map(classify_role)

    # Keep only cells within the parent context + the source/target cells
    df_plot = df[df['role'].notna()].copy()

    n_source = (df_plot['role'] == 'source').sum()
    n_target = (df_plot['role'] == 'target').sum()
    n_parent = (df_plot['role'] == 'parent').sum()

    # Use palette colours for source and target; article-style remaining parent
    #color_source = celltype_palette.get(type_i, COLOR_SOURCE_DEFAULT)
    #color_target = celltype_palette.get(type_j, COLOR_TARGET_DEFAULT)
    color_source = COLOR_SOURCE_DEFAULT
    color_target = COLOR_TARGET_DEFAULT

    role_order  = ['parent', 'source', 'target']
    role_colors = {'parent': COLOR_PARENT, 'source': color_source, 'target': color_target}
    role_sizes  = {'parent': 2,            'source': 3,            'target': 3}
    role_alphas = {'parent': 0.8,          'source': 0.85,         'target': 0.85}
    role_labels = {
        'parent': f'Remaining {parent_name} (n={n_parent:,})',
        'source': f'{type_i} (n={n_source:,})',
        'target': f'{type_j} (n={n_target:,})',
    }

    fig, ax = plt.subplots(figsize=(8, 8))

    for role in role_order:
        mask = df_plot['role'] == role
        ax.scatter(df_plot.loc[mask, 'x'], df_plot.loc[mask, 'y'],
                   c=role_colors[role],
                   s=role_sizes[role],
                   alpha=role_alphas[role],
                   linewidths=0,
                   rasterized=True)

    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_xlabel('x (µm)', fontsize=10)
    ax.set_ylabel('y (µm)', fontsize=10)
    ax.set_title(f'{sample}\n{type_i}  →  {type_j}  |  context: {parent_name}',
                 fontsize=10)
    ax.tick_params(labelsize=8)

    legend_handles = [
        mpatches.Patch(color=role_colors[r], label=role_labels[r])
        for r in role_order
    ]
    ax.legend(handles=legend_handles, fontsize=8, loc='upper right', framealpha=0.8)

    # Stats annotation
    kont_integral = pair_data['integral_kontextual']
    z_kont        = pair_data['z_kontextual']
    p_kont        = pair_data['p_kontextual']
    direction     = 'localised' if kont_integral > 0 else 'dispersed'
    ax.text(0.02, 0.02,
            f'integral_kontextual = {kont_integral:.1f}\n'
            f'z = {z_kont:.2f},  p = {p_kont:.4f}  ({direction})',
            transform=ax.transAxes, fontsize=8, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.75))

    plt.tight_layout()
    out = os.path.join(out_dir, f'{sample}_{cell_pair}_fig2D_pair_scatter.svg')
    plt.savefig(out, format='svg', bbox_inches='tight', dpi=300)
    plt.close()
    print(f'  Saved: {out}')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 6 Run
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Panel A — once per sample
print('Plotting Panel A (full scatter, all cell types)...')
plot_full_scatter(df, sample, celltype_palette, out_dir)

# Panels B and D — once per cell pair
for cell_pair in cell_pairs:
    type_i, type_j = cell_pair.split('__')
    print(f'\nCell pair: {type_i} → {type_j}')

    pair_key = (type_i, type_j)
    if pair_key not in kont_data:
        available = [f'{k[0]}__{k[1]}' for k in kont_data.keys()]
        print(f'  WARNING: pair {pair_key} not found in pkl. Skipping.')
        print('  Available pairs (first 20):\n  ' + '\n  '.join(available[:20]))
        continue

    pair_data = kont_data[pair_key]
    print(f'  Parent context: {pair_data["parent_name"]}')
    print(f'  n_source: {pair_data["n_source_cells"]}, n_target: {pair_data["n_target_cells"]}')

    print('  Plotting Panel B (curves)...')
    plot_curves(pair_data, sample, type_i, type_j, cell_pair, out_dir)

    print('  Plotting Panel D (pair scatter)...')
    plot_pair_scatter(df, pair_data, sample, type_i, type_j, cell_pair,
                      celltype_palette, out_dir)

print('\nDone.')
