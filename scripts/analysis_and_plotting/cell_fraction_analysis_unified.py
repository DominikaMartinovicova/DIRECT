#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# cell_fraction_analysis_unified.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Unified cell type fraction analysis script.
#   Merges functionality from:
#       - cell_fraction_analysis copy.py   (Biopsy vs Resection, full plot suite)
#       - cell_fraction_analysis.py        (T_number-level, fold change, lineplots)
#       - cell_fraction_analysis_core_wise.py (structure-wise comparison)
#
#   Two analysis modes:
#       bvr       - Compare Biopsy vs Resection cell type fractions
#       structure - Compare fractions between tissue structures within resection
#                   samples (e.g. tumor core vs tumor bed vs margin)
#       both      - Run both modes
#
#   Plot types (mix and match with --plot_types):
#       box         - Boxplot comparing fractions between two groups
#       line        - Stripplot with lines connecting paired samples
#       foldchange  - Log2 fold change per patient (bvr mode only)
#       shift       - Absolute difference per patient: group2 - group1 (bvr only)
#       composition - Composition boxplot across all samples (bvr mode only)
#
#   Each plot type is run for:
#       (a) All cell types
#       (b) Immune cell types only (re-normalised to immune-only sum)
#   and optionally stratified by each --categories column (e.g. MPR, treatment).
#
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage examples:
#
#   # Biopsy vs Resection — boxplot and fold change:
#   python3 cell_fraction_analysis_unified.py \
#       -i /path/to/adata.h5ad \
#       --phen_level cell_type \
#       --celltype_list /path/to/celltypes.txt \
#       --analysis_mode bvr \
#       --plot_types box foldchange line \
#       --categories MPR treatment \
#       --exclude_v17 \
#       --output_dir_plots /path/to/plots/ \
#       --output_dir_results /path/to/results/
#
#   # Structure-wise (core-wise) comparison within resection samples:
#   python3 cell_fraction_analysis_unified.py \
#       -i /path/to/adata.h5ad \
#       --phen_level cell_type \
#       --celltype_list /path/to/celltypes.txt \
#       --analysis_mode structure \
#       --structure_cols structure structure_core \
#       --plot_types box \
#       --categories MPR \
#       --output_dir_plots /path/to/plots/ \
#       --output_dir_results /path/to/results/
#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0  Imports and argument parsing
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import os
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import wilcoxon, mannwhitneyu
from statannotations.Annotator import Annotator

warnings.filterwarnings("ignore")

# Metadata columns expected in adata.obs (missing ones are silently skipped)
_META_COLS = [
    'sample', 'pt_id', 'sample_type', 'disease_stage', 'T_number',
    'regression', 'treatment_scheme', 'MPR', 'treatment',
    'structure', 'structure_core',
]

# Cell types considered non-immune (excluded from immune-only analysis)
_NON_IMMUNE = {
    'Epithelial_cell', 'Fibroblast', 'Endothelial_cell',
    'Pericyte', 'Stromal', 'Tumor_cells',
}


# ------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        prog='python3 cell_fraction_analysis_unified.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description='Unified cell type fraction analysis script.')

    parser.add_argument(
        '-i', dest='input', type=str, required=True,
        help='Path to combined .h5ad file')
    parser.add_argument(
        '--phen_level', dest='phen_level', type=str, required=True,
        help='Key for cell type annotation in adata.obs')
    parser.add_argument(
        '--celltype_list', dest='celltype_list', type=str, required=True,
        help='Path to TSV file listing cell types to include (one per line, no header)')
    parser.add_argument(
        '--exclude_v17', action='store_true',
        help='Exclude samples with v1.7 treatment scheme')
    parser.add_argument(
        '--analysis_mode', dest='analysis_mode',
        choices=['bvr', 'structure', 'both'], default='bvr',
        help=(
            'Analysis mode (default: bvr):\n'
            '  bvr       - Biopsy vs Resection\n'
            '  structure - Structure-wise comparison within resection\n'
            '  both      - Run both modes'))
    parser.add_argument(
        '--plot_types', dest='plot_types', nargs='+',
        choices=['box', 'line', 'foldchange', 'shift', 'composition', 'within_sampletype'],
        default=['box', 'foldchange'],
        help=(
            'Plot types to generate, space-separated (default: box foldchange):\n'
            '  box               - Boxplot: fractions per group\n'
            '  line              - Stripplot with paired connecting lines\n'
            '  foldchange        - Log2 fold change per patient [bvr only]\n'
            '  shift             - Absolute difference per patient [bvr only]\n'
            '  composition       - Composition boxplot across all samples [bvr only]\n'
            '  within_sampletype - Fractions within Biopsy/Resection split by\n'
            '                      category (requires --categories) [bvr only]'))
    parser.add_argument(
        '--stat_test', dest='stat_test', default='mannwhitneyu',
        choices=['mannwhitneyu', 'wilcoxon'],
        help=(
            'Statistical test for independent-samples comparisons (default: mannwhitneyu).\n'
            'Note: the line plot always uses wilcoxon (paired test) regardless of this setting.'))
    parser.add_argument(
        '--groupby_key', dest='groupby_key', default='T_number',
        help='obs column used to compute per-unit fractions (default: T_number)')
    parser.add_argument(
        '--categories', dest='categories', nargs='*', default=None,
        help=(
            'Optional metadata columns to stratify analyses by\n'
            '(e.g. --categories MPR treatment). Each is run in addition\n'
            'to the unstratified analysis.'))
    parser.add_argument(
        '--structure_cols', dest='structure_cols', nargs='+',
        default=['structure', 'structure_core'],
        help='obs columns defining tissue structures [structure mode]\n(default: structure structure_core)')
    parser.add_argument(
        '--exclude_core3', action='store_true',
        help='Exclude core_3 observations when computing fractions [bvr mode]')
    parser.add_argument(
        '-o', '--output_dir_results', dest='output_dir_results',
        type=str, required=True,
        help='Output directory for statistical results (CSV files)')
    parser.add_argument(
        '--output_dir_plots', dest='output_dir_plots',
        type=str, required=True,
        help='Output directory for plots (SVG files)')

    return parser.parse_args()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1  Data helpers
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def compute_fractions(adata, groupby_key, phen_level):
    """
    Compute per-unit cell type fractions and attach available metadata.

    Parameters
    ----------
    adata       : AnnData
    groupby_key : str   obs column to iterate over (e.g. 'T_number', 'sample')
    phen_level  : str   obs column containing cell type labels

    Returns
    -------
    DataFrame with one row per unit, cell type columns (0–1 fractions),
    and any metadata columns that exist in adata.obs.
    """
    meta_present = [m for m in _META_COLS if m in adata.obs.columns]
    acc = {}

    for unit in adata.obs[groupby_key].dropna().unique():
        sub = adata[adata.obs[groupby_key] == unit]
        fracs = sub.obs[phen_level].value_counts() / sub.shape[0]
        row = fracs.to_dict()
        for meta in meta_present:
            row[meta] = sub.obs[meta].iloc[0]
        acc[unit] = row

    df = pd.DataFrame.from_dict(acc, orient='index').fillna(0)
    df.index.name = groupby_key
    return df.reset_index()


def filter_v17(df, exclude_v17):
    """Remove v1.7 treatment scheme rows if requested."""
    if exclude_v17 and 'treatment_scheme' in df.columns:
        before = len(df)
        df = df[~df['treatment_scheme'].str.contains('v1.7', na=False)].copy()
        print(f'  v1.7 exclusion: {before} -> {len(df)} rows')
    return df


def get_paired(df, group_col='sample_type'):
    """Keep only patients that appear in both groups (paired analysis)."""
    return df.groupby('pt_id').filter(lambda x: x[group_col].nunique() == 2).copy()


def aggregate_per_patient(df, group_col, cell_type_list):
    """
    Average cell type fractions per patient per group.
    Used before computing per-patient fold change / shift.

    Returns a DataFrame with one row per (pt_id, group_col).
    """
    meta_cols = [c for c in ['pt_id', group_col, 'MPR', 'treatment',
                              'structure', 'structure_core', 'regression',
                              'treatment_scheme']
                 if c in df.columns]
    # Keep only one metadata row per (pt_id, group_col)
    meta = df[meta_cols].drop_duplicates(subset=['pt_id', group_col])
    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    ct_mean = df.groupby(['pt_id', group_col])[avail_ct].mean()
    result = ct_mean.reset_index().merge(meta, on=['pt_id', group_col], how='left')
    return result


def get_immune_fractions(df, cell_type_list, preserve_cols):
    """
    Subset df to immune cell types and re-normalise each row to sum to 1.

    Parameters
    ----------
    df             : DataFrame with cell type fraction columns
    cell_type_list : list of cell type column names currently active
    preserve_cols  : list of metadata column names to carry over

    Returns
    -------
    (df_immune, immune_ct_list)
    """
    to_exclude = _NON_IMMUNE.intersection(set(cell_type_list))
    immune_ct_list = [ct for ct in cell_type_list if ct not in to_exclude and ct in df.columns]
    if not immune_ct_list:
        print('  Warning: no immune cell types remain after exclusion.')
        return None, []

    df_cells = df[immune_ct_list].copy()
    row_sums = df_cells.sum(axis=1).replace(0, np.nan)
    df_cells = df_cells.div(row_sums, axis=0)

    avail_meta = [c for c in preserve_cols if c in df.columns]
    df_cells[avail_meta] = df[avail_meta].values
    return df_cells, immune_ct_list


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2  Statistical testing helpers
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def stat_testing_two_groups(df, cell_cols, stat_test, group, groups):
    """
    Run stat_test between two groups for each cell type.

    Returns
    -------
    stat_df      : raw results DataFrame
    stat_df_annot: formatted for statannotations.Annotator
    """
    results = []
    for ct in cell_cols:
        g1 = df[df[group] == groups[0]][ct].dropna()
        g2 = df[df[group] == groups[1]][ct].dropna()
        if g1.empty or g2.empty:
            continue
        stat, p = stat_test(g1, g2)
        results.append({'cell_type': ct, 'statistic': stat, 'p_value': p})

    if not results:
        return pd.DataFrame(), pd.DataFrame()

    stat_df = pd.DataFrame(results)
    stat_df_annot = (
        stat_df
        .rename(columns={'cell_type': 'variable', 'p_value': 'pval'})
        .assign(group1=groups[0], group2=groups[1])
        [['variable', 'group1', 'group2', 'pval']]
    )
    return stat_df, stat_df_annot


def annotate_significant(ax, stat_df_annot, data, x, y, hue, alpha=0.05):
    """Add significance stars to ax for p < alpha."""
    if stat_df_annot.empty:
        return
    sig = stat_df_annot[stat_df_annot['pval'] < alpha].reset_index(drop=True)
    if sig.empty:
        return
    pairs = [((r.variable, r.group1), (r.variable, r.group2)) for _, r in sig.iterrows()]
    try:
        annot = Annotator(ax, pairs, data=data, x=x, y=y, hue=hue)
        annot.configure(text_format='star')
        annot.set_pvalues_and_annotate(sig['pval'])
    except Exception as e:
        print(f'  Warning: annotation failed ({e})')


def save_stat_results(stat_df, output_dir, prefix, group, category, stat_test, immune, exclude_v17):
    """Save statistical results DataFrame to CSV."""
    if stat_df.empty:
        return
    suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    imm = '_immune' if immune else ''
    cat = f'_{category}' if category else ''
    fname = f'{prefix}{imm}_{group}{cat}_{stat_test.__name__}_{suffix}.csv'
    fpath = os.path.join(output_dir, fname)
    stat_df.to_csv(fpath, index=False)
    print(f'  Stats: {fpath}')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3  Filename builder and shared plot utilities
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def build_plot_path(output_dir, prefix, group, category, stat_test, immune, exclude_v17):
    suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    imm = '_immune' if immune else ''
    cat = f'_{category}' if category else ''
    return os.path.join(output_dir, f'{prefix}{imm}_{group}{cat}_{stat_test.__name__}_{suffix}.svg')


def _save_and_close(fpath):
    plt.tight_layout()
    plt.savefig(fpath, format='svg', bbox_inches='tight')
    plt.close()
    print(f'  Plot:  {fpath}')


def _set_labels_and_title(g_or_ax, cell_type_list, title, key, legend_title='', use_catplot=False):
    """Apply rotated x-tick labels, axis labels and title."""
    if use_catplot:
        g_or_ax.set_xticklabels(rotation=45, ha='right')
        g_or_ax.set_xlabels('Cell type')
        g_or_ax.set_ylabels(key)
        plt.suptitle(title, y=1.03)
        if g_or_ax.legend is not None:
            g_or_ax.legend.set_title(legend_title)
            g_or_ax.legend.set_loc('upper right')
    else:
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(key)


def _build_title(base, groups, category, stat_test, immune, exclude_v17):
    imm = 'immune ' if immune else ''
    title = f'{imm}{base}: {groups[0]} vs {groups[1]}'
    if category:
        title += f' | split by {category}'
    if exclude_v17:
        title += ' (excl. v1.7)'
    title += f' ({stat_test.__name__})'
    return title


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4  Plot functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------------------------------
def plot_box(df, cell_type_list, group, category, id_col,
             output_dir_plots, output_dir_results,
             stat_test, immune, exclude_v17, key='Cell fraction'):
    """
    Boxplot of cell type fractions, comparing two groups.
    Optionally faceted by category.
    """
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip box] fewer than 2 groups in "{group}": {groups}')
        return

    id_vars = [id_col, group] + ([category] if category else [])
    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    melt = df[id_vars + avail_ct].melt(id_vars=id_vars, var_name='cell_type', value_name=key)

    col_order = sorted([g for g in df[category].dropna().unique()
                        if not str(g).isdigit()]) if category else None

    g = sns.catplot(melt, x='cell_type', y=key, hue=group, hue_order=groups,
                    col=category, col_order=col_order,
                    kind='box', palette='tab20', height=6, aspect=1.5)

    axes = g.axes.flat if category else [g.ax]
    facet_data = list(g.facet_data()) if category else [(None, df)]

    for ax, (_, subdata) in zip(axes, facet_data):
        if category:
            subset_df = subdata.pivot(index=id_vars, columns='cell_type', values=key).reset_index()
            subset_melt = subdata
        else:
            subset_df = subdata[id_vars + avail_ct]
            subset_melt = subset_df.melt(id_vars=id_vars, value_vars=avail_ct,
                                          var_name='cell_type', value_name=key)
        stat_df, stat_df_annot = stat_testing_two_groups(subset_df, avail_ct, stat_test, group, groups)
        annotate_significant(ax, stat_df_annot, subset_melt, 'cell_type', key, group)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results, 'stats_box', group,
                              category, stat_test, immune, exclude_v17)

    title = _build_title('Cell fractions', groups, category, stat_test, immune, exclude_v17)
    _set_labels_and_title(g, avail_ct, title, key, legend_title=group, use_catplot=True)
    fpath = build_plot_path(output_dir_plots, 'cf_box', group, category, stat_test, immune, exclude_v17)
    _save_and_close(fpath)


# ------------------------------------------------------------------------------
def plot_line(df, cell_type_list, group, category, id_col,
              output_dir_plots, output_dir_results,
              stat_test, immune, exclude_v17, key='Cell fraction'):
    """
    Stripplot with lines connecting paired samples per patient (id_col groups).
    Blue lines = increase, red lines = decrease between groups.
    """
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip line] fewer than 2 groups in "{group}": {groups}')
        return

    id_vars = [id_col, group] + ([category] if category else [])
    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    melt = df[id_vars + avail_ct].melt(id_vars=id_vars, var_name='cell_type', value_name=key)
    col_order = sorted([g for g in df[category].dropna().unique()
                        if not str(g).isdigit()]) if category else None

    g = sns.catplot(melt, x='cell_type', y=key, hue=group, hue_order=groups,
                    col=category, col_order=col_order,
                    kind='strip', palette={groups[0]: 'gray', groups[1]: 'black'},
                    jitter=False, dodge=True, height=6, aspect=1.5, size=4)

    axes = g.axes.flat if category else [g.ax]
    facet_data = list(g.facet_data()) if category else [(None, df)]
    offsets = np.linspace(-0.2, 0.2, len(groups))

    for ax, (_, subdata) in zip(axes, facet_data):
        if category:
            subset_df = subdata.pivot(index=id_vars, columns='cell_type', values=key).reset_index()
            subset_melt = subdata
        else:
            subset_df = subdata[id_vars + avail_ct]
            subset_melt = subset_df.melt(id_vars=id_vars, value_vars=avail_ct,
                                          var_name='cell_type', value_name=key)

        # Draw connecting lines between paired observations
        for i, cell in enumerate(avail_ct):
            cell_data = subset_melt[subset_melt['cell_type'] == cell]
            for pt, pt_df in cell_data.groupby(id_col):
                pt_df = pt_df.set_index(group)
                if groups[0] not in pt_df.index or groups[1] not in pt_df.index:
                    continue
                y1 = pt_df.loc[groups[0], key]
                y2 = pt_df.loc[groups[1], key]
                color = 'blue' if y2 > y1 else 'red'
                ax.plot([i + offsets[0], i + offsets[1]], [y1, y2],
                        color=color, alpha=0.6, linewidth=1)

        stat_df, stat_df_annot = stat_testing_two_groups(subset_df, avail_ct, stat_test, group, groups)
        annotate_significant(ax, stat_df_annot, subset_melt, 'cell_type', key, group)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results, 'stats_line', group,
                              category, stat_test, immune, exclude_v17)

    title = _build_title('Cell fractions (paired)', groups, category, stat_test, immune, exclude_v17)
    _set_labels_and_title(g, avail_ct, title, key, legend_title=group, use_catplot=True)
    fpath = build_plot_path(output_dir_plots, 'cf_line', group, category, stat_test, immune, exclude_v17)
    _save_and_close(fpath)


# ------------------------------------------------------------------------------
def plot_foldchange(df, cell_type_list, group, category,
                    output_dir_plots, output_dir_results,
                    stat_test, immune, exclude_v17, key='Cell fraction'):
    """
    Log2 fold change boxplot: log2(group[1] / group[0]) per patient.
    df must already be aggregated to one row per (pt_id, group).
    """
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip foldchange] fewer than 2 groups in "{group}": {groups}')
        return

    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    PSEUDO = 1e-4

    df_ref = df[df[group] == groups[0]].set_index('pt_id')[avail_ct] + PSEUDO
    df_tgt = (df[df[group] == groups[1]]
              .set_index('pt_id')
              .reindex(df_ref.index)[avail_ct] + PSEUDO)

    fc_df = np.log2(df_tgt.div(df_ref)).replace([np.inf, -np.inf], np.nan)

    # Attach metadata
    meta_cols = ['pt_id'] + ([category] if category else [])
    avail_meta = [c for c in meta_cols if c in df.columns]
    meta_map = df[avail_meta].drop_duplicates(subset=['pt_id'])
    fc_df = fc_df.reset_index().merge(meta_map, on='pt_id', how='left')

    plt.figure(figsize=(12, 6))
    if category is None:
        fc_melt = fc_df[avail_ct].melt(var_name='cell_type', value_name=key)
        ax = sns.boxplot(data=fc_melt, x='cell_type', y=key, palette='tab20')
    else:
        fc_melt = fc_df.melt(id_vars=['pt_id', category], value_vars=avail_ct,
                              var_name='cell_type', value_name=key)
        cat_order = sorted(fc_melt[category].dropna().unique())
        ax = sns.boxplot(data=fc_melt, x='cell_type', y=key, hue=category,
                         hue_order=cat_order, palette='tab20')
        # Compare between categories (independent samples)
        stat_df, stat_df_annot = stat_testing_two_groups(
            fc_df, avail_ct, stat_test, category, cat_order[:2])
        annotate_significant(ax, stat_df_annot, fc_melt, 'cell_type', key, category)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results, 'stats_foldchange', group,
                              category, stat_test, immune, exclude_v17)

    plt.axhline(y=0, color='red', linestyle='--', linewidth=1)

    imm = 'immune ' if immune else ''
    title = f'Log2 FC {imm}cell fractions: {groups[1]} / {groups[0]}'
    if category:
        title += f' | split by {category}'
    if exclude_v17:
        title += ' (excl. v1.7)'
    title += f' ({stat_test.__name__})'
    plt.title(title)
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Cell type')
    plt.ylabel(f'Log2 FC ({groups[1]} / {groups[0]})')
    fpath = build_plot_path(output_dir_plots, 'cf_fc_box', group, category, stat_test, immune, exclude_v17)
    _save_and_close(fpath)


# ------------------------------------------------------------------------------
def plot_shift(df, cell_type_list, group, category,
               output_dir_plots, output_dir_results,
               stat_test, immune, exclude_v17, key='Cell fraction'):
    """
    Absolute difference boxplot: group[1] - group[0] per patient.
    df must already be aggregated to one row per (pt_id, group).
    """
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip shift] fewer than 2 groups in "{group}": {groups}')
        return

    avail_ct = [ct for ct in cell_type_list if ct in df.columns]

    df_ref = df[df[group] == groups[0]].set_index('pt_id')[avail_ct]
    df_tgt = (df[df[group] == groups[1]]
              .set_index('pt_id')
              .reindex(df_ref.index)[avail_ct])
    diff_df = (df_tgt - df_ref).reset_index()

    # Attach metadata
    meta_cols = ['pt_id'] + ([category] if category else [])
    avail_meta = [c for c in meta_cols if c in df.columns]
    meta_map = df[avail_meta].drop_duplicates(subset=['pt_id'])
    diff_df = diff_df.merge(meta_map, on='pt_id', how='left')

    plt.figure(figsize=(12, 6))
    if category is None:
        diff_melt = diff_df[avail_ct].melt(var_name='cell_type', value_name=key)
        ax = sns.boxplot(data=diff_melt, x='cell_type', y=key, palette='tab20')
    else:
        diff_melt = diff_df.melt(id_vars=['pt_id', category], value_vars=avail_ct,
                                  var_name='cell_type', value_name=key)
        cat_order = sorted(diff_melt[category].dropna().unique())
        diff_melt[category] = pd.Categorical(diff_melt[category],
                                              categories=cat_order, ordered=True)
        ax = sns.boxplot(data=diff_melt, x='cell_type', y=key, hue=category, palette='tab20')
        stat_df, stat_df_annot = stat_testing_two_groups(
            diff_df, avail_ct, stat_test, category, cat_order[:2])
        annotate_significant(ax, stat_df_annot, diff_melt, 'cell_type', key, category)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results, 'stats_shift', group,
                              category, stat_test, immune, exclude_v17)

    plt.axhline(y=0, color='red', linestyle='--', linewidth=1)

    imm = 'immune ' if immune else ''
    title = f'Fraction shift ({groups[1]} - {groups[0]}) {imm}'
    if category:
        title += f' | split by {category}'
    if exclude_v17:
        title += ' (excl. v1.7)'
    title += f' ({stat_test.__name__})'
    plt.title(title)
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Cell type')
    plt.ylabel(f'Fraction shift ({groups[1]} - {groups[0]})')
    fpath = build_plot_path(output_dir_plots, 'cf_shift_box', group, category, stat_test, immune, exclude_v17)
    _save_and_close(fpath)


# ------------------------------------------------------------------------------
def plot_composition(df, cell_type_list, group, category,
                     output_dir_plots, output_dir_results,
                     stat_test, immune, exclude_v17, key='Cell fraction'):
    """
    Composition boxplot: raw fractions across all samples per group.
    No pairing required.
    """
    groups = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])
    if len(groups) < 2:
        print(f'  [skip composition] fewer than 2 groups in "{group}": {groups}')
        return

    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    id_vars = ['pt_id', group] + ([category] if category else [])
    melt = df[id_vars + avail_ct].melt(id_vars=id_vars, var_name='cell_type', value_name=key)

    if category is None:
        plt.figure(figsize=(12, 6))
        ax = sns.boxplot(data=melt, x='cell_type', y=key, hue=group,
                         hue_order=groups, palette='tab20', fill=True)
        stat_df, stat_df_annot = stat_testing_two_groups(df, avail_ct, stat_test, group, groups)
        annotate_significant(ax, stat_df_annot, melt, 'cell_type', key, group)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results, 'stats_composition', group,
                              category, stat_test, immune, exclude_v17)
        title = _build_title('Cell fraction composition', groups, category, stat_test, immune, exclude_v17)
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(key)
        plt.legend(title=group, loc='upper right')
        fpath = build_plot_path(output_dir_plots, 'cf_composition_box', group, category,
                                stat_test, immune, exclude_v17)
        _save_and_close(fpath)

        # Also save a zoomed-in version
        plt.figure(figsize=(12, 6))
        ax2 = sns.boxplot(data=melt, x='cell_type', y=key, hue=group,
                          hue_order=groups, palette='tab20', fill=True)
        ax2.set_ylim(-0.01, 0.1)
        plt.title(title + ' (zoomed)')
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(key)
        plt.legend(title=group, loc='upper right')
        zoomed = fpath.replace('.svg', '_zoom_0_0.1.svg')
        _save_and_close(zoomed)

    else:
        col_order = sorted([g for g in df[category].dropna().unique()
                            if not str(g).isdigit()])
        g = sns.catplot(melt, x='cell_type', y=key, hue=group, hue_order=groups,
                        col=category, col_order=col_order,
                        kind='box', palette='tab20', height=6, aspect=1.5)
        for ax, (_, subdata) in zip(g.axes.flat, g.facet_data()):
            subset_df = subdata.pivot(index=id_vars, columns='cell_type', values=key).reset_index()
            stat_df, stat_df_annot = stat_testing_two_groups(subset_df, avail_ct, stat_test, group, groups)
            annotate_significant(ax, stat_df_annot, subdata, 'cell_type', key, group)

        title = _build_title('Cell fraction composition', groups, category, stat_test, immune, exclude_v17)
        _set_labels_and_title(g, avail_ct, title, key, legend_title=group, use_catplot=True)
        fpath = build_plot_path(output_dir_plots, 'cf_composition_box', group, category,
                                stat_test, immune, exclude_v17)
        _save_and_close(fpath)


# ------------------------------------------------------------------------------
def plot_within_sampletype(df, cell_type_list, group, category,
                           output_dir_plots, output_dir_results,
                           stat_test, immune, exclude_v17, key='Cell fraction'):
    """
    For each value of `group` (e.g. Biopsy, Resection), plot cell type fractions
    split by `category` (e.g. MPR). This answers: within Resection samples, do
    MPR and non-MPR patients differ in their cell type fractions?

    Requires category to be set; silently skips if category is None.
    """
    if category is None:
        print('  [skip within_sampletype] no category specified')
        return

    avail_ct = [ct for ct in cell_type_list if ct in df.columns]
    group_vals = sorted([g for g in df[group].dropna().unique() if not str(g).isdigit()])

    for sample_val in group_vals:
        sub = df[df[group] == sample_val].copy()
        cat_vals = sorted([c for c in sub[category].dropna().unique() if not str(c).isdigit()])
        if len(cat_vals) < 2:
            print(f'  [skip within_sampletype] fewer than 2 values of {category} in {sample_val}')
            continue

        id_vars = ['pt_id', category]
        melt = sub[id_vars + avail_ct].melt(id_vars=id_vars, var_name='cell_type', value_name=key)

        plt.figure(figsize=(12, 6))
        ax = sns.boxplot(data=melt, x='cell_type', y=key, hue=category,
                         hue_order=cat_vals, palette='tab20')

        stat_df, stat_df_annot = stat_testing_two_groups(
            sub, avail_ct, stat_test, category, cat_vals[:2])
        annotate_significant(ax, stat_df_annot, melt, 'cell_type', key, category)
        if not stat_df.empty:
            save_stat_results(stat_df, output_dir_results, f'stats_within_{sample_val}',
                              category, None, stat_test, immune, exclude_v17)

        suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
        imm = '_immune' if immune else ''
        title = f'{"Immune " if immune else ""}Cell fractions in {sample_val} | split by {category}'
        if exclude_v17:
            title += ' (excl. v1.7)'
        title += f' ({stat_test.__name__})'
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(key)
        plt.legend(title=category, loc='upper right')

        fname = f'cf_within_{sample_val}{imm}_{category}_{stat_test.__name__}_{suffix}.svg'
        fpath = os.path.join(output_dir_plots, fname)
        _save_and_close(fpath)

        # Zoomed version
        plt.figure(figsize=(12, 6))
        ax2 = sns.boxplot(data=melt, x='cell_type', y=key, hue=category,
                          hue_order=cat_vals, palette='tab20')
        ax2.set_ylim(-0.01, 0.1)
        plt.title(title + ' (zoomed)')
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(key)
        plt.legend(title=category, loc='upper right')
        _save_and_close(fpath.replace('.svg', '_zoom_0_0.1.svg'))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5  Analysis runners
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Maps plot type name -> function
_PLOT_FUNCS = {
    'box':               plot_box,
    'line':              plot_line,
    'foldchange':        plot_foldchange,
    'shift':             plot_shift,
    'composition':       plot_composition,
    'within_sampletype': plot_within_sampletype,
}

# These require per-patient aggregated data (not per T_number/sample)
_PATIENT_LEVEL_PLOTS = {'foldchange', 'shift'}

# These are only meaningful for BVR mode (not structure-wise)
_BVR_ONLY_PLOTS = {'foldchange', 'shift', 'composition', 'within_sampletype'}

# Line plot is always paired -> always use wilcoxon regardless of --stat_test
_PAIRED_PLOTS = {'line'}


def _run_plots(df, df_pt, cell_type_list, group, categories, id_col,
               plot_types, output_dir_plots, output_dir_results,
               stat_test, immune, exclude_v17, mode_label=''):
    """
    Dispatch selected plot types for all category combinations.

    df     : row-per-unit data (T_number or sample) for box/line/composition
    df_pt  : row-per-patient data (averaged) for foldchange/shift
    stat_test : default statistical test (used for all plots except line,
                which always uses wilcoxon because it shows paired data)
    """
    print(f'\n  --- {mode_label} | immune={immune} ---')
    for ptype in plot_types:
        func = _PLOT_FUNCS[ptype]
        use_pt_level = ptype in _PATIENT_LEVEL_PLOTS
        # Line plot visualises paired data -> always use the paired wilcoxon test
        effective_stat_test = wilcoxon if ptype in _PAIRED_PLOTS else stat_test

        for category in categories:
            print(f'  [{ptype}] group={group}, category={category}')
            data = df_pt if use_pt_level else df
            kwargs = dict(
                cell_type_list=cell_type_list,
                group=group,
                category=category,
                output_dir_plots=output_dir_plots,
                output_dir_results=output_dir_results,
                stat_test=effective_stat_test,
                immune=immune,
                exclude_v17=exclude_v17,
            )
            if ptype in {'box', 'line'}:
                kwargs['id_col'] = id_col
            func(data, **kwargs)


def run_bvr_analysis(fractions_df, cell_type_list, args):
    """
    Biopsy vs Resection analysis.
    Runs all selected plot types for all cells and immune-only.
    """
    print('\n======  Biopsy vs Resection  ======')
    group = 'sample_type'
    id_col = args.groupby_key
    categories = [None] + (args.categories or [])
    plot_types = args.plot_types
    stat_test = wilcoxon if args.stat_test == 'wilcoxon' else mannwhitneyu

    df = filter_v17(fractions_df, args.exclude_v17)

    # For box/line/composition: keep only samples from patients with both sample types
    paired_df = get_paired(df, group_col=group)
    print(f'  Paired patients (both Biopsy + Resection): {paired_df["pt_id"].nunique()}')

    # For foldchange/shift: aggregate to one row per patient per sample_type
    paired_pt_df = aggregate_per_patient(paired_df, group, cell_type_list)

    # -- All cell types --
    _run_plots(
        df=paired_df,
        df_pt=paired_pt_df,
        cell_type_list=cell_type_list,
        group=group,
        categories=categories,
        id_col=id_col,
        plot_types=plot_types,
        output_dir_plots=args.output_dir_plots,
        output_dir_results=args.output_dir_results,
        stat_test=stat_test,
        immune=False,
        exclude_v17=args.exclude_v17,
        mode_label='All cell types',
    )

    # -- Immune cell types only --
    preserve_meta = ['pt_id', group] + (args.categories or [])
    df_immune, immune_ct_list = get_immune_fractions(paired_df, cell_type_list, preserve_meta)
    if df_immune is not None and immune_ct_list:
        df_immune_pt = aggregate_per_patient(df_immune, group, immune_ct_list)
        _run_plots(
            df=df_immune,
            df_pt=df_immune_pt,
            cell_type_list=immune_ct_list,
            group=group,
            categories=categories,
            id_col=id_col,
            plot_types=plot_types,
            output_dir_plots=args.output_dir_plots,
            output_dir_results=args.output_dir_results,
            stat_test=stat_test,
            immune=True,
            exclude_v17=args.exclude_v17,
            mode_label='Immune cell types only',
        )


def run_structure_analysis(fractions_df, cell_type_list, args):
    """
    Structure-wise analysis within resection samples.
    Pools (averages) multiple TMA cores per patient per structure type,
    then compares between structures.
    """
    print('\n======  Structure-wise analysis  ======')
    categories = [None] + (args.categories or [])
    # Exclude plot types that only make sense for BVR (paired patient comparisons)
    plot_types = [pt for pt in args.plot_types if pt not in _BVR_ONLY_PLOTS]
    stat_test = wilcoxon if args.stat_test == 'wilcoxon' else mannwhitneyu

    df = filter_v17(fractions_df, args.exclude_v17)
    res_df = df[df['sample_type'] == 'Resection'].copy()
    print(f'  Resection samples: {res_df.shape[0]}')

    # Metadata needed for merging back after groupby
    meta_cols = [c for c in ['pt_id', 'MPR', 'treatment', 'regression', 'treatment_scheme',
                              'disease_stage', 'sample_type']
                 if c in res_df.columns and c not in ['structure', 'structure_core']]
    meta_df = res_df[meta_cols].drop_duplicates(subset=['pt_id'])

    for structure_col in args.structure_cols:
        if structure_col not in res_df.columns:
            print(f'  Warning: "{structure_col}" not found in data, skipping.')
            continue

        print(f'\n  Structure column: {structure_col}')
        # Pool samples: average fractions per patient per structure type
        pooled = res_df.groupby(['pt_id', structure_col], as_index=False).mean(numeric_only=True)
        pooled = pooled.merge(meta_df, on='pt_id', how='left')

        # -- All cell types --
        _run_plots(
            df=pooled,
            df_pt=pooled,  # not used (no patient-level plots in structure mode)
            cell_type_list=cell_type_list,
            group=structure_col,
            categories=categories,
            id_col='pt_id',
            plot_types=plot_types,
            output_dir_plots=args.output_dir_plots,
            output_dir_results=args.output_dir_results,
            stat_test=stat_test,
            immune=False,
            exclude_v17=args.exclude_v17,
            mode_label=f'All cell types | {structure_col}',
        )

        # Also run with groups swapped (category as group, structure as facet)
        for category in [c for c in categories if c is not None]:
            print(f'  [{structure_col}] swapped: group={category}, category={structure_col}')
            for ptype in [pt for pt in plot_types if pt in {'box'}]:
                plot_box(pooled, cell_type_list, group=category, category=structure_col,
                         id_col='pt_id',
                         output_dir_plots=args.output_dir_plots,
                         output_dir_results=args.output_dir_results,
                         stat_test=stat_test, immune=False,
                         exclude_v17=args.exclude_v17)

        # -- Immune cell types only --
        preserve_meta = ['pt_id', structure_col] + (args.categories or [])
        df_immune, immune_ct_list = get_immune_fractions(pooled, cell_type_list, preserve_meta)
        if df_immune is not None and immune_ct_list:
            _run_plots(
                df=df_immune,
                df_pt=df_immune,
                cell_type_list=immune_ct_list,
                group=structure_col,
                categories=categories,
                id_col='pt_id',
                plot_types=plot_types,
                output_dir_plots=args.output_dir_plots,
                output_dir_results=args.output_dir_results,
                stat_test=stat_test,
                immune=True,
                exclude_v17=args.exclude_v17,
                mode_label=f'Immune cell types | {structure_col}',
            )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 6  Main
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def main():
    args = parse_args()

    os.makedirs(args.output_dir_plots, exist_ok=True)
    os.makedirs(args.output_dir_results, exist_ok=True)
    sns.set_style('whitegrid')

    # Print run summary
    print('=' * 60)
    print('Cell fraction analysis — unified script')
    print('=' * 60)
    print(f'  Input:           {args.input}')
    print(f'  Phenotype level: {args.phen_level}')
    print(f'  Groupby key:     {args.groupby_key}')
    print(f'  Analysis mode:   {args.analysis_mode}')
    print(f'  Plot types:      {args.plot_types}')
    print(f'  Categories:      {args.categories}')
    print(f'  Exclude v1.7:    {args.exclude_v17}')
    print(f'  Stat test:       {args.stat_test} (line plot always uses wilcoxon)')
    if args.analysis_mode in ('structure', 'both'):
        print(f'  Structure cols:  {args.structure_cols}')
    print(f'  Output plots:    {args.output_dir_plots}')
    print(f'  Output results:  {args.output_dir_results}')
    print('=' * 60)

    # Load data
    print('\nReading data...')
    adata = sc.read_h5ad(args.input)
    print(adata)

    cell_type_list = pd.read_csv(args.celltype_list, header=None, sep='\t')[0].tolist()
    print(f'\nCell types loaded: {len(cell_type_list)}')

    # Optionally exclude core_3 for BVR mode
    if args.analysis_mode in ('bvr', 'both') and args.exclude_core3:
        before = adata.n_obs
        adata = adata[adata.obs['structure_core'] != 'core_3', :]
        print(f'  Excluded core_3: {before} -> {adata.n_obs} cells')

    # Compute fractions
    print(f'\nComputing fractions per {args.groupby_key}...')
    fractions_df = compute_fractions(adata, args.groupby_key, args.phen_level)
    print(f'  Fractions dataframe: {fractions_df.shape}')

    # Run selected analyses
    if args.analysis_mode in ('bvr', 'both'):
        run_bvr_analysis(fractions_df, cell_type_list, args)

    if args.analysis_mode in ('structure', 'both'):
        if args.groupby_key != 'sample':
            print('\n  Note: structure mode is typically run with --groupby_key sample.')
        run_structure_analysis(fractions_df, cell_type_list, args)

    print('\nDone.')


if __name__ == '__main__':
    main()
