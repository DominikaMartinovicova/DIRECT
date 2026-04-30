#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# stat_analysis_kontextual.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Compare kontextual Ripley's L statistics between groups.
#   Supports flexible comparison variables and optional stratification.
#   Samples are merged by cell-count-weighted average of integral_kontextual
#   at the finest grain needed: (T_number x compare_col x split_col x cell_pair).
#   This preserves sample_type distinctions within the same patient when
#   comparing or splitting by sample_type.
#
#   0 Import libraries and set parameters
#   1 Load data (once, shared across all combinations)
#       a. kontextual_ripleys_summary.tsv per sample
#       b. cell counts from h5ad
#       c. samples metadata
#   2–4 For each (compare_col, split_col) combination:
#       2. Weighted aggregation
#       3. Statistical comparison (Mann-Whitney U, FDR correction)
#          - run once per split level if split_col is given
#       4. Visualisation
#          a. Volcano plot (one per split level)
#          b. Heatmap of median integral_kontextual per group (one per split level)
#          c. Boxplots for significant cell pairs
#             (side-by-side subplots per split level if split_col is given)
#
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage (CLI mode — uncomment argparse block and comment out hardcoded block):
#   python3 stat_analysis_kontextual.py \
#       --input_dir  <per_sample_results_dir> \
#       --metadata   <samples_metadata.csv> \
#       --adata      <combined.h5ad> \
#       --compare    <MPR|sample_type|treatment> \
#       [--groups    <group1> <group2>] \
#       [--split_by  <MPR|sample_type|treatment>] \
#       [--sample_type Resection|Biopsy|both] \
#       [--exclude_v17] \
#       --o_results  <output_results_dir> \
#       --o_plots    <output_plots_dir>


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and set parameters
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import glob
import argparse
import anndata as ad
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statannotations.Annotator import Annotator

VALID_COLS = ['MPR', 'sample_type', 'treatment']


# def parse_args():
#     "Parse inputs from commandline and return as Namespace object."
#     parser = argparse.ArgumentParser(
#         prog='python3 stat_analysis_kontextual.py',
#         formatter_class=argparse.RawTextHelpFormatter,
#         description='  Compare kontextual Ripley statistics between groups.')
#     parser.add_argument('--input_dir', dest='input_dir', type=str, required=True,
#                         help='Directory containing per-sample kontextual result folders.')
#     parser.add_argument('--metadata', dest='metadata', type=str, required=True,
#                         help='Path to samples_metadata.csv.')
#     parser.add_argument('--adata', dest='adata', type=str, required=True,
#                         help='Path to analyzed h5ad file (for cell counts).')
#     parser.add_argument('--compare', dest='compare_col', type=str, required=True,
#                         choices=VALID_COLS,
#                         help='Metadata column defining the two groups to compare.')
#     parser.add_argument('--groups', dest='groups', type=str, nargs=2, default=None,
#                         metavar=('GROUP1', 'GROUP2'),
#                         help='Two values of --compare to test (default: auto-detect from data).')
#     parser.add_argument('--split_by', dest='split_col', type=str, default=None,
#                         choices=VALID_COLS,
#                         help='Optional: stratify by this column; produces one comparison\n'
#                              'per level with side-by-side boxplots.')
#     parser.add_argument('--sample_type', dest='sample_type', type=str,
#                         default=None, choices=['Resection', 'Biopsy', 'both'],
#                         help='Pre-filter to specific sample type (default: both).')
#     parser.add_argument('--exclude_v17', action='store_true',
#                         help='Exclude treatment scheme v1.7 samples.')
#     parser.add_argument('--o_results', dest='output_dir_results', type=str, required=True,
#                         help='Output directory for results (TSV files).')
#     parser.add_argument('--o_plots', dest='output_dir_plots', type=str, required=True,
#                         help='Output directory for plots.')
#     return parser.parse_args()
#
#
# args = parse_args()
# input_dir          = args.input_dir
# metadata_path      = args.metadata
# adata_path         = args.adata
# groups_arg         = args.groups
# sample_type_filter = args.sample_type
# exclude_v17        = args.exclude_v17
# output_dir_results = args.output_dir_results
# output_dir_plots   = args.output_dir_plots
# combinations_to_run = [(args.compare_col, args.split_col)]

# --- Hardcoded parameters (comment out and uncomment argparse block above to use CLI) ---

_BASE = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT'
input_dir          = f'{_BASE}/results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_sample'
metadata_path      = f'{_BASE}/data/adata_per_sample/Neutro_Epi_extImm_pooled_A_EM_N/samples_metadata.csv'
adata_path         = (f'{_BASE}/data/combined/Neutro_Epi_extImm_pooled_A_EM_N_combined_adatas_for_analysis_w_v1.7.h5ad')
groups_arg         = None           # e.g. ['<90', '>=90'] or None to auto-detect per combination
exclude_v17        = True
sample_type_filter = None    # 'Resection' | 'Biopsy' | 'both' | None
output_dir_results = f'{_BASE}/results/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/kontextual/Tnumber'
output_dir_plots   = f'{_BASE}/plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/kontextual/Tnumber'
combinations_to_run = [
    # (compare_col,   split_col)     comment out any lines to skip
    ('MPR',          None),
    ('MPR',          'sample_type'),
    ('MPR',          'treatment'),
    ('sample_type',  None),
    ('sample_type',  'MPR'),
    ('sample_type',  'treatment'),
    ('treatment',    None),
    ('treatment',    'MPR'),
    ('treatment',    'sample_type'),
]
# -----------------------------------------------------------------------------------------

os.makedirs(output_dir_results, exist_ok=True)
os.makedirs(output_dir_plots, exist_ok=True)

sns.set_style("whitegrid")

suffix       = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
stype_suffix = sample_type_filter if sample_type_filter else 'both'

print(f'Sample type filter: {stype_suffix}')
print(f'Excluding v1.7:     {exclude_v17}')
print(f'Combinations to run: {len(combinations_to_run)}')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Load data (shared across all combinations)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1a. Load kontextual_ripleys_summary.tsv for every sample
#------------------------------------------------------------------------------
summary_files = glob.glob(
    os.path.join(input_dir, '*/kontextual_ripleys/kontextual_ripleys_summary.tsv'))
print(f'Found {len(summary_files)} kontextual summary files.')

dfs = []
for f in summary_files:
    sample_id = f.split(os.sep)[-3]
    df = pd.read_csv(f, sep='\t')
    df['sample'] = sample_id
    dfs.append(df)

kontextual_df = pd.concat(dfs, ignore_index=True)
kontextual_df['cell_pair'] = kontextual_df['from'] + '__' + kontextual_df['to']
print(f'Cell pairs found: {kontextual_df["cell_pair"].nunique()}')

# 1b. Cell counts per sample from h5ad
#------------------------------------------------------------------------------
print('Loading h5ad for cell counts (obs only)...')
adata = ad.read_h5ad(adata_path, backed='r')
cell_counts = adata.obs['sample'].value_counts().rename('n_cells').reset_index()
cell_counts.columns = ['sample', 'n_cells']
adata.file.close()

# 1c. Samples metadata
#------------------------------------------------------------------------------
meta = pd.read_csv(metadata_path)

if sample_type_filter and sample_type_filter != 'both':
    meta = meta[meta['sample_type'] == sample_type_filter].copy()
    print(f'Samples after {sample_type_filter} filter: {len(meta)}')

if exclude_v17:
    meta = meta[meta['treatment_scheme'] != 'v1.7'].copy()
    print(f'Samples after v1.7 exclusion: {len(meta)}')

meta_cols = [c for c in ['sample', 'T_number', 'MPR', 'sample_type', 'treatment',
                          'treatment_scheme', 'pt_id'] if c in meta.columns]

merged = (kontextual_df
          .merge(meta[meta_cols], on='sample', how='inner')
          .merge(cell_counts, on='sample', how='left'))

missing_counts = merged['n_cells'].isna().sum()
if missing_counts > 0:
    print(f'WARNING: {missing_counts} rows have missing cell counts; using equal weights.')
    merged['n_cells'] = merged['n_cells'].fillna(1)

print(f'Samples retained after metadata merge: {merged["sample"].nunique()}')
print(f'T_numbers retained: {merged["T_number"].nunique()}')


# Helper functions
#------------------------------------------------------------------------------
def weighted_mean(values, weights):
    w = np.array(weights, dtype=float)
    v = np.array(values, dtype=float)
    mask = ~np.isnan(v)
    if mask.sum() == 0:
        return np.nan
    return np.average(v[mask], weights=w[mask])


def run_stats(df, compare_col, groups):
    "Mann-Whitney U per cell pair between groups[0] and groups[1]; FDR-BH correction."
    results = []
    for cell_pair, grp in df.groupby('cell_pair'):
        g1 = grp[grp[compare_col] == groups[0]]['integral_kontextual'].dropna()
        g2 = grp[grp[compare_col] == groups[1]]['integral_kontextual'].dropna()
        if len(g1) < 2 or len(g2) < 2:
            continue
        stat, pval = mannwhitneyu(g1, g2, alternative='two-sided')
        results.append({
            'cell_pair':     cell_pair,
            'group1':        groups[0],
            'group2':        groups[1],
            'n_group1':      len(g1),
            'n_group2':      len(g2),
            'median_group1': g1.median(),
            'median_group2': g2.median(),
            'median_diff':   g2.median() - g1.median(),
            'statistic':     stat,
            'p_value':       pval,
        })
    stat_df = pd.DataFrame(results)
    if stat_df.empty:
        return stat_df
    reject, pvals_adj, _, _ = multipletests(stat_df['p_value'], method='fdr_bh')
    stat_df['p_adj']       = pvals_adj
    stat_df['significant'] = reject
    return stat_df.sort_values('p_adj')


def plot_volcano(stat_df, groups, compare_col, title, out_path):
    if stat_df.empty:
        return
    plot_df = stat_df.dropna(subset=['median_diff'])
    if plot_df.empty:
        return

    y_vals      = -np.log10(plot_df['p_adj'].clip(lower=1e-10))
    sig_threshold = -np.log10(0.05)
    near_thresh = (~plot_df['significant']) & y_vals.between(1, sig_threshold)
    colors = [
        'red'    if sig else
        'orange' if near else
        'grey'
        for sig, near in zip(plot_df['significant'], near_thresh)
    ]
    sig_df  = plot_df[plot_df['significant']]
    near_df = plot_df[near_thresh]
    xlabel  = f'Median diff ({compare_col}): {groups[1]} − {groups[0]}'
    ylabel  = '-log10(FDR-adjusted p-value)'

    fig, (ax_lin, ax_log) = plt.subplots(1, 2, figsize=(18, 7), sharey=True)

    for ax, log_scale in ((ax_lin, False), (ax_log, True)):
        ax.scatter(plot_df['median_diff'], y_vals,
                   c=colors, alpha=0.7, s=40, linewidths=0)
        ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.6)
        ax.axvline(0,               color='black', linestyle='--', linewidth=0.8, alpha=0.6)
        for _, row in sig_df.iterrows():
            ax.annotate(
                row['cell_pair'].replace('__', '\n→ '),
                xy=(row['median_diff'], -np.log10(max(row['p_adj'], 1e-10))),
                fontsize=9, ha='center', va='bottom',
                xytext=(0, 4), textcoords='offset points')
        for _, row in near_df.iterrows():
            ax.annotate(
                row['cell_pair'].replace('__', '\n→ '),
                xy=(row['median_diff'], -np.log10(max(row['p_adj'], 1e-10))),
                fontsize=8, ha='center', va='bottom',
                xytext=(0, 4), textcoords='offset points')
        if log_scale:
            # symlog handles values that cross zero; linthresh sets the linear region width
            linthresh = max(1.0, np.abs(plot_df['median_diff']).quantile(0.05))
            ax.set_xscale('symlog', linthresh=linthresh)
            ax.set_xlabel(f'{xlabel} [symlog scale]', fontsize=13)
        else:
            ax.set_xlabel(xlabel, fontsize=13)
        ax.set_ylabel(ylabel, fontsize=13)


    ax_lin.set_title(f'{title}\n(linear x)', fontsize=12)
    ax_log.set_title(f'{title}\n(symlog x)', fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, format='svg', bbox_inches='tight')
    plt.close()


def plot_heatmap(patient_sub, groups, compare_col, title, out_path, row_order=None):
    if patient_sub.empty:
        return
    heatmap_data = (patient_sub
                    .groupby([compare_col, 'cell_pair'])['integral_kontextual']
                    .median()
                    .unstack(compare_col))
    heatmap_data = heatmap_data[[c for c in groups if c in heatmap_data.columns]]
    if row_order is not None:
        heatmap_data = heatmap_data.reindex(row_order)
    else:
        heatmap_data = heatmap_data.sort_index()
    vmax = heatmap_data.abs().quantile(0.95).max()
    fig, ax = plt.subplots(figsize=(4, max(6, len(heatmap_data) * 0.28)))
    sns.heatmap(
        heatmap_data, cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
        ax=ax, linewidths=0.3,
        cbar_kws={'label': 'Median integral_kontextual'})
    ax.set_title(title, fontsize=10)
    ax.set_ylabel('Cell pair')
    ax.set_xlabel(compare_col)
    plt.tight_layout()
    plt.savefig(out_path, format='svg', bbox_inches='tight')
    plt.close()


def plot_combined_heatmap(patient_df, split_levels, split_col, groups, compare_col,
                          title, out_path, row_order=None):
    "Single heatmap with columns grouped by split level: [sv1_g1, sv1_g2, sv2_g1, sv2_g2, ...]."
    parts = []
    for sv in split_levels:
        sub = (patient_df[patient_df[split_col] == sv]
               if sv is not None else patient_df)
        sub = sub[sub[compare_col].isin(groups)]
        if sub.empty:
            continue
        part = (sub.groupby([compare_col, 'cell_pair'])['integral_kontextual']
                   .median()
                   .unstack(compare_col))
        part = part[[c for c in groups if c in part.columns]]
        part.columns = [f'{sv}\n{g}' for g in part.columns]
        parts.append(part)

    if not parts:
        return

    combined = pd.concat(parts, axis=1)
    combined = combined.reindex(row_order) if row_order is not None else combined.sort_index()

    vmax = combined.abs().quantile(0.95).max()
    n_cols = len(combined.columns)
    fig, ax = plt.subplots(figsize=(max(4, n_cols * 1.2), max(6, len(combined) * 0.28)))
    sns.heatmap(
        combined, cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
        ax=ax, linewidths=0.3,
        cbar_kws={'label': 'Median integral_kontextual'})

    # Vertical separators between split-level groups
    n_groups = len(groups)
    for i in range(1, len(parts)):
        ax.axvline(i * n_groups, color='black', linewidth=2)

    ax.set_title(title, fontsize=10)
    ax.set_ylabel('Cell pair')
    ax.set_xlabel('')
    plt.tight_layout()
    plt.savefig(out_path, format='svg', bbox_inches='tight')
    plt.close()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2–4 Loop over all (compare_col, split_col) combinations
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot_dir = output_dir_plots

for compare_col, split_col in combinations_to_run:
    print(f'\n{"="*70}')
    print(f'Compare: {compare_col}   Split by: {split_col}')
    print(f'{"="*70}')

    if compare_col == split_col:
        print('SKIP: compare_col and split_col are the same.')
        continue
    if exclude_v17 and 'treatment' in (compare_col, split_col):
        print('SKIP: treatment comparison not possible without v1.7 samples (only milder present).')
        continue
    if compare_col == 'sample_type' and sample_type_filter and sample_type_filter != 'both':
        print(f'WARNING: compare=sample_type but sample_type filter={sample_type_filter}; '
              f'only one group will be present.')

    split_tag = f'_by_{split_col}' if split_col else ''
    file_tag  = f'{compare_col}{split_tag}_{stype_suffix}_{suffix}'

    # 2. Weighted aggregation
    # -------------------------------------------------------------------------
    # Aggregation keys: T_number + compare_col + split_col (if any) + cell_pair.
    # Keeping compare_col and split_col as separate dimensions preserves e.g.
    # the biopsy/resection distinction within the same patient when sample_type
    # is used as either the comparison or the split variable.
    agg_key_cols = list(dict.fromkeys(
        ['T_number', compare_col] + ([split_col] if split_col else []) + ['cell_pair']
    ))
    carry_cols = [c for c in ['MPR', 'sample_type', 'treatment', 'pt_id']
                  if c in merged.columns and c not in agg_key_cols]

    agg_rows = []
    for keys, grp in merged.groupby(agg_key_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(agg_key_cols, keys))
        for c in carry_cols:
            row[c] = grp[c].iloc[0]
        row['integral_kontextual']     = weighted_mean(grp['integral_kontextual'],     grp['n_cells'])
        row['integral_kontextual_abs'] = weighted_mean(grp['integral_kontextual_abs'], grp['n_cells'])
        row['z_kontextual']            = weighted_mean(grp['z_kontextual'],            grp['n_cells'])
        row['integral_original']       = weighted_mean(grp['integral_original'],       grp['n_cells'])
        row['z_original']              = weighted_mean(grp['z_original'],              grp['n_cells'])
        row['n_samples']   = len(grp)
        row['total_cells'] = grp['n_cells'].sum()
        agg_rows.append(row)

    patient_df = pd.DataFrame(agg_rows)
    patient_df.to_csv(
        os.path.join(output_dir_results, f'kontextual_patient_level_{file_tag}.tsv'),
        sep='\t', index=False)
    print(f'Patient-level data saved. Shape: {patient_df.shape}')

    # 3. Statistical comparison
    # -------------------------------------------------------------------------
    all_groups = sorted(patient_df[compare_col].dropna().unique())
    if groups_arg:
        groups = list(groups_arg)
        missing = [g for g in groups if g not in all_groups]
        if missing:
            print(f'SKIP: specified groups {missing} not found in {compare_col}: {all_groups}')
            continue
    else:
        if len(all_groups) < 2:
            print(f'SKIP: need ≥2 values in {compare_col}, found: {all_groups}')
            continue
        groups = all_groups[:2]
        print(f'Auto-detected groups for {compare_col}: {groups}')

    split_levels = (sorted(patient_df[split_col].dropna().unique())
                    if split_col else [None])

    stats_by_split = {}
    for sv in split_levels:
        sub = (patient_df[patient_df[split_col] == sv].copy()
               if sv is not None else patient_df.copy())
        sub = sub[sub[compare_col].isin(groups)]

        stat_df = run_stats(sub, compare_col, groups)
        stats_by_split[sv] = stat_df

        label = f'{split_col}={sv}' if sv is not None else 'all'
        if stat_df.empty:
            print(f'  [{label}] No testable cell pairs.')
            continue

        n_sig = stat_df['significant'].sum()
        print(f'  [{label}] Significant cell pairs (FDR < 0.05): {n_sig} / {len(stat_df)}')
        if n_sig > 0:
            print(stat_df[stat_df['significant']][
                ['cell_pair', 'median_diff', 'p_value', 'p_adj']].to_string())

        split_label = f'_{sv}' if sv is not None else ''
        stat_df.to_csv(
            os.path.join(output_dir_results, f'kontextual_stats_{file_tag}{split_label}.tsv'),
            sep='\t', index=False)

    # 4. Visualisation
    # -------------------------------------------------------------------------

    # 4a & 4b: Volcano per split level; single combined heatmap (or one if no split)
    heatmap_row_order = sorted(patient_df['cell_pair'].unique())
    base_title = f'({stype_suffix}, {suffix})'

    for sv in split_levels:
        stat_df   = stats_by_split[sv]
        split_label   = f'_{sv}' if sv is not None else ''
        title_context = f' [{split_col}={sv}]' if sv is not None else ''

        plot_volcano(
            stat_df, groups, compare_col,
            title=f'Kontextual Ripley: {groups[0]} vs {groups[1]}{title_context}\n{base_title}',
            out_path=os.path.join(plot_dir, f'volcano_{file_tag}{split_label}.svg'))

    if split_col is not None:
        # Combined heatmap: columns grouped as [sv1_g1, sv1_g2, sv2_g1, sv2_g2, ...]
        plot_combined_heatmap(
            patient_df, split_levels, split_col, groups, compare_col,
            title=f'Median integral_kontextual — by {split_col}\n{base_title}',
            out_path=os.path.join(plot_dir, f'heatmap_{file_tag}.svg'),
            row_order=heatmap_row_order)
    else:
        sub = patient_df[patient_df[compare_col].isin(groups)].copy()
        plot_heatmap(
            sub, groups, compare_col,
            title=f'Median integral_kontextual\n{base_title}',
            out_path=os.path.join(plot_dir, f'heatmap_{file_tag}.svg'),
            row_order=heatmap_row_order)

    # 4c: Boxplots for significant cell pairs
    #     No split  → single figure
    #     With split → side-by-side subplots, one per split level
    all_sig_pairs = []
    for sv in split_levels:
        stat_df = stats_by_split[sv]
        if not stat_df.empty:
            all_sig_pairs.extend(stat_df[stat_df['significant']]['cell_pair'].tolist())
    sig_pairs_ordered = list(dict.fromkeys(all_sig_pairs))

    if not sig_pairs_ordered:
        print('  No significant cell pairs; skipping boxplots.')
        continue

    box_cols = (['T_number', compare_col] + ([split_col] if split_col else []) +
                ['cell_pair', 'integral_kontextual'])
    box_df = (patient_df[box_cols].copy()
              .dropna(subset=['integral_kontextual', compare_col])
              .query(f'{compare_col} in @groups')
              .query('cell_pair in @sig_pairs_ordered'))

    n_splits = len(split_levels)
    n_pairs  = len(sig_pairs_ordered)
    panel_w  = max(6, n_pairs * 0.9)
    fig, axes = plt.subplots(1, n_splits,
                             figsize=(panel_w * n_splits, 6),
                             squeeze=False)

    for col_idx, sv in enumerate(split_levels):
        ax  = axes[0, col_idx]
        sub = (box_df[box_df[split_col] == sv].copy()
               if sv is not None else box_df.copy())

        if sub.empty:
            ax.set_visible(False)
            continue

        sns.boxplot(
            data=sub, x='cell_pair', y='integral_kontextual',
            hue=compare_col, hue_order=groups, palette='Paired',
            ax=ax, order=sig_pairs_ordered)

        stat_df = stats_by_split[sv]
        if not stat_df.empty:
            sig_here = stat_df[stat_df['significant'] &
                               stat_df['cell_pair'].isin(sig_pairs_ordered)]
            if not sig_here.empty:
                annot_pairs = [((cp, groups[0]), (cp, groups[1]))
                               for cp in sig_here['cell_pair']]
                annot_pvals = sig_here.set_index('cell_pair').loc[
                    sig_here['cell_pair'], 'p_adj'].tolist()
                try:
                    annot = Annotator(ax, annot_pairs, data=sub,
                                      x='cell_pair', y='integral_kontextual',
                                      hue=compare_col, order=sig_pairs_ordered)
                    annot.configure(text_format='star')
                    annot.set_pvalues_and_annotate(annot_pvals)
                except Exception as e:
                    print(f'  Annotation failed for {split_col}={sv}: {e}')

        subplot_title = (f'{split_col}={sv}' if sv is not None
                         else f'{compare_col}: {groups[0]} vs {groups[1]}')
        ax.set_title(subplot_title, fontsize=10)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        ax.set_xlabel('Cell pair')
        ax.set_ylabel('integral_kontextual (weighted avg)')
        ax.legend(title=compare_col, loc='upper right')

    fig.suptitle(
        f'Kontextual Ripley — significant pairs\n({stype_suffix}, {suffix})',
        fontsize=11)
    plt.tight_layout()
    plt.savefig(
        os.path.join(plot_dir, f'boxplot_{file_tag}_significant.svg'),
        format='svg', bbox_inches='tight')
    plt.close()

print('\nDone. All plots and results saved.')
